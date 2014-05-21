## Functions that simplify data parsing and model creation and fitting using the package

library("igraph")
library("pheatmap")
library("parallel")

source("R/randomLHS.r"); # Latin Hypercube Sampling

# Global variable to have more outputs
verbose = FALSE;
debug = TRUE;

create_model <- function(model.links="links", data.stimulation="data", basal_activity = "basal.dat", data.variation="", cores=0, nb_init=1000)
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
    nb_samples = nb_init;
    print (paste("Initializing the model parameters… (", nb_samples, " random samplings) with ", cores, " cores", sep=""))
    samples = qnorm(randomLHS(nb_samples, model$nr_of_parameters()));
    # Parallelized version uses all cores but one to keep control
    if (cores == 0) {
        results = parallel_initialisation(model, expdes, model.structure, data, samples, detectCores()-1)
    } else {
        results = parallel_initialisation(model, expdes, model.structure, data, samples, cores)
    }
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
    #model_description$name = ;

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
        # Simplify the name of the paths, keep the old one
        lprofile$old_path = lprofile$path;
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
    pheatmap(mismatch, color=colorRampPalette(c("blue", "black", "red"))(length(bk-1)), breaks = bk, cluster_rows=F, cluster_col=F, display_numbers = T)
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

# Print the value of each path
print_parameters <- function(model_description) {
    model = model_description$model;
    parameters = model_description$parameters;

    print("Parameters :")
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

# Selection of a minimal model
select_minimal_model <- function(model_description, accuracy=0.95) {
    model = model_description$model;
    init_params = model_description$parameters;
    initial_response = model$getLocalResponseFromParameter( init_params );
    expdes = model_description$design;
    model_structure = model_description$structure;
    data = model_description$data;

    print("Performing model reduction…");
    reduce=TRUE;
    while (reduce) {
        links.to.test=which(adj==1)
        params=c();
        residuals=c();
        ranks=c();
        newadj=adj;

        # Each link is removed and the best of those networks is compared to the previous model
        newadj=adj;
        for (i in links.to.test) {
            newadj[i]=0;
            model_structure$setAdjacencymatrix( newadj );
            model$setModel ( expdes, model_structure );
            paramstmp = model$getParameterFromLocalResponse(initial.response$local_response, initial.response$inhibitors);
            result = model$fitmodel( data,paramstmp)
            response.matrix = model$getLocalResponseFromParameter( result$parameter )
            residuals = c(residuals,result$residuals);   
            params = cbind(params,c(response.matrix));
            new_rank = model$modelRank();
            ranks = c(ranks, new_rank);

            if (verbose) {
                dr = rank - new_rank;
                print(paste("old :", rank, ", new : ", new_rank));
                deltares = residuals[length(residuals)]-initresidual;
                print(paste(model_structure$names[(i-1) %/% dim(adj)[1]+1], "->", model_structure$names[(i-1) %% dim(adj)[1]+1], ": Delta residual = ", deltares, "; Delta rank = ", dr, ", p-value = ", pchisq(deltares, df=dr) ));
            }

            newadj[i]=1; ## Slightly accelerate the computation
        }
        
        order.res=order(residuals);
        # The loss of degree of freedom is equal to the difference in the ranks of the matrices
        new_rank = ranks[order.res[1]];
        dr = rank - new_rank;
        deltares = residuals[order.res[1]]-initresidual;
        # Some boundary cases might give low improvement of the fit
        if (deltares < 0) { deltares = -deltares; warning(paste("Negative delta residual :", deltares)) }
        if (deltares < qchisq(accuracy, df=dr)) {
            adj[links.to.test[order.res[1]]]=0;
            rank = new_rank;
            initial.response=params[,order.res[1]];
            print(paste0("remove ",
                      model_structure$names[((links.to.test[order.res[1]]-1) %/% (dim(adj)[1])) +1], "->", # Line
                      model_structure$names[((links.to.test[order.res[1]]-1) %% (dim(adj)[1])) +1])); # Column (+1 because of the modulo and the R matrices starting by 1 instead of 0)

            print(paste( "residual = ", residuals[order.res[1]], ", Delta residual = ", residuals[order.res[1]]-initresidual, ",  p-value = ", pchisq(deltares, df=dr) ));
            print("------------------------------------------------------------------------------------------------------");

        } else {
          reduce=FALSE;
        }
    }
    #print(adj)
    # We recover the final model
    model_description$model_structure$setAdjacencymatrix(adj);
    model_description$model$setModel(expdes, model_description$model_structure);

    return(model_description);
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

    for (row in 1:dim(model_description$data$stim_data)[1]) {
        text = pastecoma("", "")
        for (col in 1:dim(model_description$data)[2]) {

        }
    }
}

# TODO
# Import model from a file
import_model <- function(file_name) {
    print("Not implemented yet")

    model_description = list()
    model = new(fitmodel::model);

    file = readLines(file_name)
    lnb = 1;
    if (!grepl("^[NF]", file[lnb])) {
        stop("This is not a mra model file.")
    }
    # Model name (cell line, network, ...)
    if (grepl("^F", file[lnb])) {
        model_description$name = gsub("^F( |\t)", "", file[1]);
        lnb = lnb + 1;
    }

    # Names of the nodes of the model
    nodes = c();
    while (grepl("^N", file[lnb])) {
        nodes = c(nodes, gsub("^N( |\t)", "", file[lnb]))
        lnb = lnb + 1;
    }

    # Get the format of the network and put it in a modelStructure object
    if (grepl("[mM](atrix|ATRIX)", file[lnb])) { # Adjacency matrix
        lnb=lnb+1;
        links = c();
        headers = unlist(strsplit(",| |\t", file[lnb]))
        headers = headers[3:length(headers)]
        lnb=lnb+1;
        while (grepl("^L"), file[lnb]) {
            line = unlist(strsplit(",| |\t", file[lnb]))

            lnb=lnb+1;
        }
    } else if (grepl("[Aa](djacency|DJACENCY)", file[lnb])) { # Adjacency list
        lnb = lnb + 1;
        print("Adjacency list not implemented yet")
    } else { # List of links
        if (grepl("list|LIST", file[lnb])) { lnb = lnb + 1; }
        print("Link list not implemented yet")
        links = c()
        while (grepl("^L", file[lnb])) {
            line = gsub("^L( |\t)", "", file[lnb]);
            line = unlist(strsplit(line, " |\t"));
            links = rbind(links, c(line[1], line[2]));

            lnb = lnb + 1;
        }
        model_structure = getModelStructure(links);
        model_description$structure = model_structure;
    }

    if (grepl("^I", file[lnb])) {
        while (grepl("^I", file[lnb])) {
            line = unlist(strsplit(" |\t", file[lnb]))


            lnb = lnb + 1;
        }
    } else {
        print("No inhibitors found !")
    }

    # Get the unstimulated data
#    model_description$data = data;

#    model_description$model = model;
#    model_description$design = expdes;
#    model_description$parameters = init_params;
#    model_description$bestfit = init_residual;
#    model_description$fileName = data.stimulation;
#    model_description$basal = basal.activity;
}

simulateModelWithTargets <- function(model_description, ) {
}

# TODO
# Get the target simulations in a MIDAS design-like or list format
# Returns a matrix in a MIDAS measure-like format
simulateModel <- function(model_description, targets, readouts = "all") {
    design = model_description$design;
    nodes = model_description$structure$names;

    stim_nodes = design$stim_nodes;
    inhib_nodes = design$inhib_nodes;
    measured_nodes = design$measured_nodes;
    # Get the experimental design constraints on the prediction capacity (index+1 because the arrays start at 0)
    perturbables = nodes[ 1 + unique(c(stim_nodes, inhib_nodes)) ];
    perID = 1 + unique(c(stim_nodes, inhib_nodes));
    measurables = nodes[ 1 + unique(c(measured_nodes, stim_nodes, inhib_nodes)) ];
    measID = 1 + unique(c(measured_nodes, stim_nodes, inhib_nodes));

    # Get the names of the nodes in the network to stimulate and inhibit, and the matrix of the perturbation is MIDAS-like format
    if (is.matrix(targets)) {
        # Already in matrix form
        inhibitions = gsub("i$", "", targets[grepl("i$", colnames(targets))] )
        stimulations = targets[!grepl("i$", colnames(targets))]
        target_matrix = targets;
    } else if (is.list(targets)) { # TODO distinguish between numeric and character
        # List of perturbation giving nodes names
        target_names = unlist(targets)
        inhibitions = unique(c( inhibitions, gsub("i$", "", target_names[grepl("i$", target_names)]) ))
        stimulations = unique(c( stimulations, target_names[!grepl("i$", target_names)] ))
        target_matrix = rep(0, length(c(stimulations, inhibitions))
        colnames(target_matrix) = c(stimulations, inhibitions);

        for (combination in targets) {
            line = 
        }
        colnames(target_matrix) = c(stimulations, inhibitions);
    }

    inhib_nodes = c();
    inhibitor = c();
    # Set the inhibition matrices
    for (node in inhibitions) {
        if (!(node %in% perturbables)) {
            plot(paste0("Node ", node, " is not in the network and won't be used"))
        } else {
            inhib_nodes = cbind(inhib_nodes, which(nodes == node)-1)
            inhibitor = cbind(inhibitors, target_matrix[, which(nodes == node)])
        }
    }

    stim_nodes = c();
    # Set the stimuli matrices
    for (node in stimulations) {
        if (!(node %in% perturbables)) {
            plot(paste0("Node ", node, " is not in the network and won't be used"))
        } else {

        }
    }

    # Set up the model for the simulation
    model = new(fitmodel::Model)
    model$setModel( design, model_description$structure )
    params = model$getParameterFromLocalResponse(model_description$model$getLocalResponse(model_description$parameters))
    prediction = model$simulate(model_description$data, params)$prediction;
    # Only unstim_data is required for the simulation

    rm(model) # Free the memory

}

# Create the perturbation matrix for a set of perturbations, building all n-combinations of stimulators with all m-combinations of inhibitors. Add the cases with only stimulations and only inhibitions
get_matrix_combination <- function (perturbations, inhib_combo = 2, stim_combo = 1) {
    stimulators = perturbations[!grepl("i$", perturbations)]
    if (stim_combo > length(stimulators) ) {
        stop ("Not enough stimulations to build the combinations")
    }
    inhibitors = perturbations[grepl("i$", perturbations)]
    if (inhib_combo > length(inhibitors) ) {
        stop ("Not enough inhibitions to build the combinations")
    }

    # Create the inhibition matrix
    inhib_combos = build_combo(seq(length(inhibitors)), inhib_combo, c())
    tmp = rep(0, length(inhibitors))
    inhib_matrix = c();
    for (i in 1:nrow(inhib_combos)) {
        inhib_matrix = rbind(inhib_matrix, tmp);
        inhib_matrix[i, inhib_combos[i,]] = 1;
    }
    colnames(inhib_matrix) = inhibitors;

    # Create the stimulation matrix
    stim_combos = build_combo(seq(length(stimulators)), stim_combo, c())
    tmp = rep(0, length(stimulators))
    stim_matrix = c();
    for (i in 1:nrow(stim_combos)) {
        stim_matrix = rbind(stim_matrix, tmp);
        stim_matrix[i, stim_combos[i,]] = 1;
    }
    colnames(stim_matrix) = stimulators;

    # Merge the two matrices to get all the combinations
    perturbation_matrix = inhib_matrix;
    ## Start merging inhibitions only and stimulations only
    for (i in 1:ncol(stim_matrix)) {
        perturbation_matrix = cbind(perturbation_matrix, 0)
    }
    colnames(perturbation_matrix) = c(inhibitors, stimulators);
    tmp = rep(0, nrow(stim_matrix))
    tmp2 = c()
    for (i in 1:ncol(inhib_matrix)) {
        tmp2 = cbind(tmp2, tmp)
    }
    perturbation_matrix = rbind(perturbation_matrix, cbind(tmp2, stim_matrix))
    ## Combination of perturbations and inhibitions
    for (i in 1:nrow(inhib_matrix)) {
        perturbation_matrix = rbind(perturbation_matrix, cbind( matrix(rep(inhib_matrix[i,], nrow(stim_matrix)), nrow(stim_matrix), ncol(inhib_matrix) , byrow=T), stim_matrix ))
    }
    rownames(perturbation_matrix) = NULL;

    return(perturbation_matrix)
}

# Recursively build the n choose k combinations for a set
build_combo <- function (symbols, remaining_steps, to_extend) {
    if (remaining_steps <= 0) {
        return(to_extend);
    } else if (length(symbols) < remaining_steps) {
        return(c());
    }
    final = c()
    # Add each symbol and deepen the recursion with remaining symbols beyond the selected one
    for ( index in seq(length(symbols)) ) {
        extension = cbind(to_extend, symbols[index])
        extension = build_combo(symbols[ -(1:index) ], remaining_steps-1, extension)
        final = rbind(final, extension)
    }
    return(final);
}





