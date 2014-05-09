## Functions that simplify data parsing and model creation and fitting using the package

library("igraph")
library("pheatmap")
library("parallel")

source("R/randomLHS.r"); # Latin Hypercube Sampling

# Global variable to have more outputs
verbose = FALSE;
debug = TRUE;

create_model <- function(model.links="links", data.stimulation="data", basal_activity = "basal.dat", data.variation="", multi_thread = FALSE)
{
# Creates a parametrized model from an experiment file and the network structure
# It requires the file network_reverse_engineering-X.X/r_binding/fitmodel/R/generate_model.R of the fitmodel package
# model.links should be an adjacency list file representing the network
# Experiment file data.stimulation syntax should be as follows, with one line per replicate
#          stimulator                |          inhibitor                |                         type                       | [one column per measured nodes]
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
# stimulator name or solvant if none | inhibitor name or solvant if none | c for control, b for blank and t for experiment |    measure for the condition

    ### READ DATA
    print("Reading data")
    # Creation of the model structure object
    links = read.delim(model.links, header=FALSE)
    model.structure=getModelStructure(links)
    model_graph = graph.edgelist(as.matrix(links))
    # Plot the network in a file
    pdf(gsub(".tab$", ".pdf", model.links))
    plot.igraph(model_graph, edge.arrow.size=0.5, layout=layout.fruchterman.reingold.grid)
    dev.off()

    # Read the experiment design and extract the values
    use_midas = FALSE
    if (grepl(".data$", data.stimulation)) {
        data.file = read.delim(data.stimulation)
        data.file[data.file=="Medium"] = "DMSO"
        begin_measure = 4
        # Indicate where the conditions are
        conditions = c(1, 2)
        data.values = data.file[, colnames(data.file) %in% model.structure$names]
    } else if (grepl(".csv$", data.stimulation)) {
        use_midas = TRUE
        data.file = read.delim(data.stimulation, sep=",")
        begin_measure = which(grepl("^DA.", colnames(data.file)))
        data.file = data.file[-begin_measure] # Delete the DA field which is not used
        begin_measure = begin_measure[1] # If there were several DA fields
        # Indicate where the conditions are
        conditions = 2:(begin_measure-1)

        # Extract the measurements of nodes in the network
        data.values = data.file[,grepl("^DV", colnames(data.file))]
        colnames(data.file) = gsub("^[A-Z]{2}.", "", colnames(data.file))
        colnames(data.values) = gsub("^[A-Z]{2}.", "", colnames(data.values))
        data.values = data.values[, colnames(data.values) %in% model.structure$names]
    }

    # Means of basal activity of the network and of the blank fixation of the antibodies
    unstim.values = colMeans(data.values[data.file$type=="c",])
    blank.values = colMeans(data.values[data.file$type=="blank",])
    blank.values[is.nan(blank.values)] = 0; # For sample without blank values

    # Calculates the mean and standard deviation for each condition 
    mean.values = aggregate(as.list(data.values),by=data.file[,1:(begin_measure-1)],mean);
    sd.values = aggregate(as.list(data.values),by=data.file[,1:(begin_measure-1)],sd);
    print("Data used :")
    print(mean.values)

    # Separate values and perturbation
    data.stim = mean.values[mean.values$type=="t",begin_measure:dim(mean.values)[2]];
    data.perturb = mean.values[mean.values$type=="t", conditions]

    ### CALCULATE ERROR MODEL
    if (grepl("\\.cv$", data.variation) || grepl("\\.var$", data.variation)) {
        # We use the CV file if there is one
        # The format and the order of the conditions are assumed to be the same as the data file
        print("Using var file")
        if (use_midas) {
            variation.file = read.delim(data.variation, sep=",")
            pre_cv = variation.file[, grepl("^DV", colnames(variation.file))]
            colnames(pre_cv) = gsub("^[A-Z]{2}.", "", colnames(pre_cv))
            cv.values = aggregate(as.list( pre_cv[colnames(pre_cv) %in% model.structure$names] ), by=data.file[,1:(begin_measure-1)], mean);
        } else {
            variation.file = read.delim(data.variation)
            cv.values = aggregate(as.list(variation.file[, colnames(data.values) %in% model.structure$names]),by=data.file[,1:(begin_measure-1)],mean);
        }
        cv.stim = cv.values[cv.values$type=="t", begin_measure:dim(cv.values)[2]];
        error = matrix(rep(blank.values,each=dim(data.stim)[1]),nrow=dim(data.stim)[1]) + cv.stim * data.stim
    } else {
    # Define the lower and default error threshold
        min.cv=0.1;     # parameters: minimal cv=0.1
        default.cv=0.3; # parameters: default cv if there are only 2 replicates

        # Calculate error percentage
        cv.values = sd.values[begin_measure:dim(sd.values)[2]] / mean.values[begin_measure:dim(sd.values)[2]];
        # Values to close to the blank are removed because the error is not due to antibody specific binding
        cv.values[!mean.values[,begin_measure:dim(mean.values)[2]] > 2 * matrix(rep(blank.values,each=dim(mean.values)[1]), nrow=dim(mean.values)[1])] = NA;
            
        # Generation of error percentage, one cv per antibody calculated using all the replicates available, default.cv if there is only two replicate to calculate the cv
        cv = colMeans(cv.values,na.rm=TRUE)
        cv[cv<min.cv] = min.cv;
        cv[is.nan(cv)|is.na(cv)]=default.cv;

        if (FALSE) { #"Multiline comment"
        for (i in 1:dim(cv.values)[2]) {
            count = 0;
            for (j in 1:dim(cv.values)[1]) {
                if (!is.na(cv[i][j]) | !is.nan(cv[i][j])) {
                    count = count + 1;
                }
            }
            if (count <= 2) {
                cv[i] = default.cv
                print("Defaulted");
            }
        }
        }

        if (verbose) {
            print("Error model :");
            for (i in 1:length(cv)) {
                print(paste(colnames(data.values)[i], " : ", cv[i]));
            }
        }
        error = matrix(rep(blank.values,each=dim(data.stim)[1]),nrow=dim(data.stim)[1])+matrix(rep(cv,each=dim(data.stim)[1]),nrow=dim(data.stim)[1])*data.stim
    }


### SET UP DATA OBJECT

    data=new(fitmodel::Data);
    data$set_unstim_data (matrix(rep(unstim.values,each=dim(data.stim)[1]),nrow=dim(data.stim)[1]));
    data$set_scale( data$unstim_data );
    data$set_stim_data( as.matrix(data.stim) );
    data$set_error( as.matrix(error ));


### EXTRACT EXPERIMENTAL DESIGN

# Extraction of stimulated, inhibited and measured nodes
    if (use_midas) {
        names = colnames(mean.values)[conditions]
        stim_names = names[grepl("[^i]$", names)];
        stim.nodes = as.character( stim_names[ stim_names %in% model.structure$names] )
        names = gsub("i$", "", names[grepl("i$", names)])
        inhib.nodes = as.character( names[names %in% model.structure$names] )
    } else {
        stim.nodes = as.character(unique(mean.values$stimulator[mean.values$stimulator %in% model.structure$names]));
        inhib.nodes = as.character(unique(mean.values$inhibitor[mean.values$inhibitor %in% model.structure$names]));
    }
    measured.nodes=colnames(data.stim);

## Identification of nodes with basal activity
    basal.activity=as.character(read.delim(basal_activity,header=FALSE)[,1]);

# Inhibition and stimulation vectors for each experiment
    if (use_midas) {
        stimuli = as.matrix(data.perturb[stim.nodes])
    } else {
        stimuli=matrix(0,ncol=length(stim.nodes),nrow=dim(data.perturb)[1])
        for (i in 1:length(stim.nodes)) {
            stimuli[grepl(stim.nodes[i], data.perturb$stimulator),i]=1;
        }
    }
    if (verbose) {
        print("Stimulated nodes");
        print(stim.nodes);
        print(stimuli);
    }
    if (use_midas) {
        inhibitor = as.matrix(data.perturb[paste(inhib.nodes, "i", sep="")])
    } else {
        inhibitor=matrix(0,ncol=length(inhib.nodes),nrow=dim(data.perturb)[1])
        if (length(inhib.nodes) > 0) { # Usefull for artificial networks
            for (i in 1:length(inhib.nodes)) {
                inhibitor[grepl(inhib.nodes[i], data.perturb$inhibitor),i]=1;
            }
        }
    }
    if (verbose) {
        print("Inhibited nodes");
        print(inhib.nodes);
        print(inhibitor);
    }

# Experimental design
    expdes=getExperimentalDesign(model.structure,stim.nodes,inhib.nodes,measured.nodes,stimuli,inhibitor,basal.activity);

### MODEL SETUP
    model = new(fitmodel::Model);
    model$setModel(expdes, model.structure);
## INITIAL FIT
    nb_samples = 100000;
    print (paste("Initializing the model parameters… (", nb_samples, " random samplings)", sep=""))
    samples = qnorm(randomLHS(nb_samples, model$nr_of_parameters()));
    # Parallelized version uses all cores but one to keep control
    results = parallel_initialisation(model, expdes, model.structure, data, samples, detectCores()-1)
    
# Choice of the best fit
    params = results$params;
    residuals = results$residuals;
    if (debug) {
        # Print the 20 smallest residuals to check if the optimum has been found several times
        print(sort(residuals)[1:20])
    }
    init_params = params[,order(residuals)[1]]
    init_residual = residuals[order(residuals)[1]]

    print("Model simulation :")
    print(model$simulate(data, init_params)$prediction);


# Information required to run the model (including the model itself)
    model_description = list()
    model_description$model = model;
    model_description$design = expdes;
    model_description$structure = model.structure;
    model_description$data = data;
    model_description$parameters = init_params;
    model_description$bestfit = init_residual;
    model_description$fileName = data.stimulation;
    model_description$basal = basal.activity;

    return(model_description);
}

# Initialise the parameters with a one core processing
classic_initialisation <- function(model, data, samples) {
    for (i in 1:nb_samples) {
        result = model$fitmodel( data, samples[i,] )
        residuals = c(residuals,result$residuals);
        params = cbind(params,result$parameter)
        if (i %% (nb_samples/20) == 0) {
            print(paste(i %/% (nb_samples/20), "/ 20 initialization done."))
        }
    }
    if (debug) {
        print(sort(residuals))
    }
    results = list()
    results$residuals = resisuals
    results$params = params
}

# Parallel initialisation of the parameters
parallel_initialisation <- function(model, expdes, structure, data, samples, NB_THREADS) {
    # Put it under list format, as randomLHS only provides a matrix
    parallel_sampling = list()
    for (i in 1:dim(samples)[1]) {
        parallel_sampling[[i]] = samples[i,]
    }
    # Parallel initialisations
    # The number of cores used depends on the ability of the detectCores function to detect them
    parallel_results = mclapply(parallel_sampling, function(params, data, model) { model$fitmodel(data, params) }, data, model, mc.cores=NB_THREADS);

    # Reorder the results to get the same output as the linear function
    results = list()
    results$residuals = c()
    results$params = c()
    for (entry in parallel_results) {
        results$residuals = c(results$residuals, entry$residuals)
        results$params = cbind(results$params, entry$parameter)
    }

    return(results)
}

# Computes the profile likelihood and the parameters relationships of each parameters in the model
# Returns a list with the profiles for each parameters
profile_likelihood <- function(model_description, nb_points=10000 , in_file=FALSE) {
### Get the information from the model description
    model = model_description$model;
    model.structure = model_description$structure;
    data = model_description$data;

    init_params = model_description$parameters;
    initial_response = model$getLocalResponseFromParameter( init_params )
    init_residual = model_description$bestfit;
    print(paste(length(init_params), " paths to evaluate"));

    profiles_list = list();
    #if (MULTITHREAD) {
##
    #}
    for (path in 1:length(init_params)) {
        lprofile = model$profileLikelihood(data, init_params, path, nb_points);
        # Residuals bigger than the simultaneous threshold are useless and would extend y axis, hidding the information
        lprofile$residuals[path, lprofile$residuals[path,] > 1.1 * lprofile$thresholds[2]] = 1.1 * lprofile$thresholds[2];
        # Simplify the name of the path
        lprofile$path = simplify_path_name(lprofile$path);
        if (FALSE) {
        # Collapse the abscisse of the big residuals to the closest valid abscisse
        last_correct = 0; to_correct = c();
        for (entry in 1:length(lprofile$residuals[path,])) {
            if (lprofile$residuals[path, entry] >= 1.1 * lprofile$thresholds[2]) {
        # Collapsing flat zones does not work for some reason, the flat zones are identified but not collapsed
        #for (entry in 1:(length(lprofile$residuals[path,])-1) ) {
          #  if (abs(lprofile$residuals[path, entry] - lprofile$residuals[path, entry+1]) < 0.001) {
                if (last_correct != 0) {
                    lprofile$explored[entry] = lprofile$explored[last_correct];
                } else {
                    to_correct = c(to_correct, entry); # the abscisse can't be corrected as long as no valid abscisse has been encountered
                }
            } else if ( (entry+1) < length(lprofile$residuals[path,]) && lprofile$residuals[path, entry+1] < 1.1 * lprofile$thresholds[2] ) {
                if (last_correct == 0) {
                    for (i in to_correct) {
                        lprofile$explored[i] = lprofile$explored[entry];
                    }
                }
                last_correct = entry;
            }
        }
        }

        print(paste0("Parameter ", lprofile$path, " decided (", path, "/", length(init_params), ")"));

        profiles_list[[path]] = lprofile;
    }

    ### Print results
    print(paste("Residual =", init_residual))
    print_error_intervals(profiles_list)

    # Print the adjacency matrix and the inhibitors values
    print("Response matrix : ");
    print(model.structure$names)
    local_response = model$getLocalResponseFromParameter(init_params);
    print("Local response : ");
    print(local_response$local_response)
    print("Inhibitors :");
    print(local_response$inhibitors);

    return(profiles_list);
}

# Plots heatmaps of the model prediction against the data weighted by the error
plot_model_accuracy <- function(model_description, data_name = "default") {
    # Calculate the mismatch
    model = model_description$model;
    data = model_description$data;
    error = data$error;
    stim_data = data$stim_data
    init_params = model_description$parameter;

    mismatch = (stim_data - model$simulate(data, init_params)$prediction) / error;

    # Rebuild the conditions from the design
    nodes = model_description$structure$names
    design = model_description$design
    treatments = c()
    for (row in 1:nrow(mismatch)) {
        stim_names = nodes[design$stim_nodes[which(design$stimuli[row,]==1)]+1];
        inhib_names = nodes[design$inhib_nodes[which(design$inhibitor[row,]==1)]+1];
        if (length(inhib_names) > 0) {
            paste(inhib_names, "i", sep="")
        }
        treatments = c(treatments, paste(c(stim_names, inhib_names), collapse="+", sep="") );
    }
    print(treatments)
    colnames(mismatch) = nodes[design$measured_nodes + 1];
    rownames(mismatch) = treatments;

    pdf(paste0("accuracy_heatmap_", data_name, ".pdf"))
    bk = unique( c(seq(min(mismatch), 0, length=50), seq(0, max(mismatch), length=50)) )
    pheatmap(mismatch, color=colorRampPalette(c("blue", "black", "red"))(length(bk-1)), breaks = bk, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", cluster_rows=F, cluster_col=F, clustering_method="ward", display_numbers = T)
    dev.off()
}

# Print each path value with its error
print_error_intervals <- function(profiles_list) {
    print("Parameters :");
        print(profiles_list[1]$lower_pointwise)
    for (i in 1:length(profiles_list)) {
        if (profiles_list[[i]]$lower_pointwise && profiles_list[[i]]$upper_pointwise) {
            print(paste( profiles_list[[i]]$path, "=", profiles_list[[i]]$value, "(", profiles_list[[i]]$lower_error[1], "-", profiles_list[[i]]$upper_error[1], ")"));
        } else if (profiles_list[[i]]$lower_pointwise) {
            print(paste( profiles_list[[i]]$path, "=", profiles_list[[i]]$value, "(", profiles_list[[i]]$lower_error[1], "- ni )"));
        } else if (profiles_list[[i]]$upper_pointwise) {
            print(paste( profiles_list[[i]]$path, "=", profiles_list[[i]]$value, "( ni -", profiles_list[[i]]$upper_error[1], ")"));
        } else {
            print(paste( profiles_list[[i]]$path, "=", profiles_list[[i]]$value, "(non identifiable)"));
        }
    }
}

# Print each path value
print_parameters <- function(model_description) {
    model = model_description$model;
    parameters = model_description$parameters;

    paths = model$getParametersLinks()
    for (i in 1:length(paths)) {
        print (paste(simplify_path_name(paths[i]), ":", parameters[i]))
    }
}

# Simplify the path name from the r_t_f* form to a f->t-> form
simplify_path_name <- function (path_name) {
    # Inhibitors do not need modification
    if (grepl("i$", path_name) || grepl("^i", path_name)) {
        return(path_name)
    }

    # Prepare a matrix 2*nb_nodes
    elements = c();
    nodes = unique(unlist(strsplit(path_name, "\\*|_|r|\\^")))
    for (node in nodes) {
        if (node != "" && !grepl("\\(-1\\)", node)) {
            elements = rbind(elements, c("", ""))
            rownames(elements)[dim(elements)[1]] = node
        }
    }
    # Indicate for each node it previous and following node in the path
    for (link in unlist(strsplit(path_name, "\\*")) ) {
        if (grepl("\\(-1\\)", link)) { # If the power is -1, we reverse the link
            nodes = unlist(strsplit(link, "_|\\^"))[3:2]
        } else {
            nodes = unlist(strsplit(link, "_"))[2:3]
        }
        elements[nodes[1], 1] = nodes[2]
        elements[nodes[2], 2] = nodes[1]
    }
    # Look for the node without predecessor
    selected = 1
    while (elements[selected, 1] != "") {
        selected = selected + 1;
    }
    # Build the simple path
    node = elements[selected, 2];
    simple_path = paste0(elements[node, 1], "->", elements[selected, 2]);
    while(elements[node, 2] != "") {
        node = elements[node, 2];
        simple_path = paste0(simple_path, "->", node);
    }

    return(simple_path);
}

# Separates the profiles whether they are identifiables or not
classify_profiles <- function (profiles_list) {
    ni_profiles = list()
    i_profiles = list();

    for (lprofile in profiles_list) {
        # Parameters are non identifiable if their profile likelihood does not reach the low threshold on both sides of the minimum 
        if (lprofile$lower_pointwise && lprofile$upper_pointwise) {
            i_profiles[[length(i_profiles)+1]] <- lprofile;
        }
        else {
            ni_profiles[[length(ni_profiles)+1]] <- lprofile;
        }
    }
    sorted_profiles = list(i_profiles, ni_profiles);
    print("Profiles sorted")
    return(sorted_profiles);
}

# Plots the functionnal relation between each non identifiable parameter and the profile likelihood of all parameters
ni_pf_plot <- function(profiles_list, init_residual=0, data_name="default") {
    # Sort the profiles to print differently whether they are identifiable or not
    sorted_profiles = classify_profiles(profiles_list)
    i_profiles = sorted_profiles[[1]];
    ni_profiles = sorted_profiles[[2]];

# Non identifiables
    nbni = length(ni_profiles);
    print(paste(nbni, "non identifiable paths"));
    # Compute the dimension, minimal size if there are not enough non identifiables
    if (nbni > 3) {
        dimension = 2 * nbni + 1;
    }
    else {
        dimension = 7;
    }

    pdf(paste0("NIplot_", data_name, ".pdf"), height=dimension, width=dimension);
    #limx = i_profiles[[1]]$
    if (nbni > 0) {
        margin = c(2, 2, 0.5, 0.5);
        par(mfcol=c(nbni, nbni), mar=margin, lab=c(3, 3, 4));
        for (ni in 1:nbni) {
            # Set the y margins
            if (ni == 1) {margin[2]=4;}
            else {margin[2]=2;}

            for (j in 1:nbni) {
                # Set the x margins
                if (j == nbni) {margin[1]=4;}
                else {margin[1]=2;}

                par(mar=margin);
                if (j == ni) {
                    limy = c( range(ni_profiles[[ni]]$residuals[ni_profiles[[j]]$pathid,])[1], ni_profiles[[ni]]$thresholds[2] * 1.1);
                }
                else {
                    limy = range(ni_profiles[[ni]]$residuals[ni_profiles[[j]]$pathid,])
                }
                #print (limy)
                ### Modification of limy if it reaches Inf, should not have to be done
                if (abs(limy[2]) == Inf) { limy[2] = 10000; print("...")}
                if (abs(limy[1]) == Inf) { limy[1] = 0}

                plot(ni_profiles[[ni]]$explored, ni_profiles[[ni]]$residuals[ni_profiles[[j]]$pathid,], xlab="", ylab="", type="l", col=abs(ni-j)+1, ylim = limy);
                # Print labels on the right and bottom
                if (ni == 1) { title(ylab=ni_profiles[[j]]$path, las=2); }
                if (j == nbni) { title(xlab=ni_profiles[[ni]]$path); }
                if (j == ni) {
                    # Could be accelerated with two points instead of hundreds
                    lines( ni_profiles[[ni]]$explored, rep(ni_profiles[[ni]]$thresholds[1], length(ni_profiles[[ni]]$explored)), lty=2, col="grey" );
                    lines( ni_profiles[[ni]]$explored, rep(ni_profiles[[ni]]$thresholds[2], length(ni_profiles[[ni]]$explored)), lty=2, col="grey" );
                    if (init_residual != 0) { lines( rep(ni_profiles[[ni]]$value, length(-5:100)), (1 + -5:100/100) * init_residual, col="red"); }
                }
            }
            print(paste("Non identifiable path", ni_profiles[[ni]]$path, "plotted"));
        }
    }

# Identifiables
    nbid = length(i_profiles);
    print(paste(nbid, "identifiable paths"));
    par( mfcol=c(1, 2), mar=c(3, 2, 0, 1), oma=c(0, 0, 2, 0) );
    if (nbid > 0) {
        for (id in 1:nbid) {
            # Plot the profile likelihood
            plot(i_profiles[[id]]$explored, i_profiles[[id]]$residuals[i_profiles[[id]]$pathid, ], type="l", sub=paste(i_profiles[[id]]$path, "profile"));
            lines( i_profiles[[id]]$explored, rep(i_profiles[[id]]$thresholds[1], length(i_profiles[[id]]$explored)), lty=2, col="grey" );
            lines( i_profiles[[id]]$explored, rep(i_profiles[[id]]$thresholds[2], length(i_profiles[[id]]$explored)), lty=2, col="grey" );
            if (init_residual != 0) { lines( rep(i_profiles[[id]]$value, length(-5:100)), (1 + -5:100/100) * init_residual, col="red"); }

            plot(1, type="n", xlim=range(i_profiles[[id]]$explored), ylim=range( i_profiles[[id]]$residuals[-i_profiles[[id]]$pathid], na.rm=T) );
            # Plot the functionnal relation
            for (i in 1:dim(i_profiles[[id]]$residuals)[1]) {
                if (i != i_profiles[[id]]$pathid) {
                    lines(i_profiles[[id]]$explored, i_profiles[[id]]$residuals[i, ], sub="Functionnal relation", col=i);
                }
            }
            title (main=i_profiles[[id]]$path, outer=T);
            print(paste("Identifiable path", i_profiles[[id]]$path, "plotted") );
        }
    }

    dev.off();

}

# Custom paste functions
pastecoma <- function(...) {
    return(paste(..., sep=","))
}
pastetab <- function(...) {
    return(paste(..., sep="\t"))
}

# TODO
# Exports the model in a file 
export_model <- function(model_description, file_name="model") {
    print("Not ready yet")
    # Add an extension
    if (!grepl(".mra$", file_name)) {
        file_name = paste0(file_name, ".mra")
    }
    handle = file(file_name, open="w")
    for (row in 1:dim(model_description$data)[1]) {
        text = pastecoma("", "")
        for (col in 1:dim(model_description$data)[2]) {

        }
    }
}

# TODO
# Import model from a file
import_model <- function(file_name) {
    print("Not implemented yet")

#    model_description = list()
#    model_description$model = model;
#    model_description$design = expdes;
#    model_description$structure = model.structure;
#    model_description$data = data;
#    model_description$parameters = init_params;
#    model_description$bestfit = init_residual;
#    model_description$fileName = data.stimulation;
#    model_description$basal = basal.activity;
}



