################################ create_model.R #########################################
## Functions that simplify data parsing and model creation and fitting using the package

#library("igraph")
#library("pheatmap")
#library("parallel")

#source("R/randomLHS.r"); # Latin Hypercube Sampling

# Global variable to have more outputs
verbose = FALSE
debug = TRUE

# Creates a parameterised model from experiment files and the network structure
createModel <- function(model_links="links", data.stimulation="data", basal_file = "basal.dat", data.variation="", cores=1, inits=1000, init_distribution=F) {

    # Creation of the model structure object
    links = read.delim(model_links, header=FALSE)
    model_structure=getModelStructure(links)
    # Plot the network in a file # TODO find a better visualisation tool
    model_graph = graph.edgelist(as.matrix(links))
    pdf(paste0( "graph_", gsub("/(\\w+/)+(\\w+).tab$", "\\2.pdf", model_links) ))
    plot.igraph(model_graph, edge.arrow.size=0.5, layout=layout.fruchterman.reingold.grid)
    dev.off()

    core = extractModelCore(model_structure, as.character(read.delim(basal_file,header=FALSE)[,1]), data.stimulation, data.variation)
    expdes = core$design
    data = core$data
    model_structure = core$structure
    basal_activity = core$basal

### MODEL SETUP
    model = new(fitmodel::Model)
    model$setModel(expdes, model_structure)
## INITIAL FIT
    print (paste("Initializing the model parameters… (", inits, " random samplings) with ", cores, " cores", sep=""))
    samples = qnorm(randomLHS(inits, model$nr_of_parameters()))
    # Parallelized version uses all cores but one to keep control
    if (cores == 0) {
        cores = detectCores()-1
    }
    results = parallel_initialisation(model, expdes, model_structure, data, samples, cores)
    # Choice of the best fit
    params = results$params
    residuals = results$residuals
    if (init_distribution) { hist(log(residuals), breaks="fd") }
    if (debug) {
        # Print the 20 smallest residuals to check if the optimum has been found several times
        print(sort(residuals)[1:20])
    }

    init_params = params[,order(residuals)[1]]
    init_residual = residuals[order(residuals)[1]]

    print("Model simulation with optimal parameters :")
    print(model$simulate(data, init_params)$prediction)


# Information required to run the model (including the model itself)
    model_description = list()
    # Objects to build the model
    model_description$model = model
    model_description$design = expdes
    model_description$structure = model_structure
    model_description$basal = basal_activity
    model_description$data = data
    # Remember the optimal parameters
    model_description$parameters = init_params
    model_description$bestfit = init_residual
    # Name and infos of the model
    model_description$name = data.stimulation
    model_description$infos = c(paste0(inits, " samplings"), paste0( sort("Best residuals : "), paste0(sort(residuals)[1:5], collapse=" ") ))
    # Values that can't be defined without the profile likelihood
    model_description$param_range = list()
    model_description$lower_values = c()
    model_description$upper_values = c()

    return(model_description)
}

# Extracts the data, the experimental design and the structure from the input files
extractModelCore <- function(model_structure, basal_activity, data_filename, var_file="") {
# links must be a matrix of links [node1, node2]
# Experiment file dta_file syntax should be as follows, with one line per replicate
#          stimulator                |          inhibitor                |                         type                       | [one column per measured nodes]
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
# stimulator name or solvant if none | inhibitor name or solvant if none | c for control, b for blank and t for experiment |    measure for the condition
    ### READ DATA
    print("Reading data")
    # Read the experiment design and extract the values
    use_midas = FALSE
    if (grepl(".data$", data_filename)) {
        data_file = read.delim(data_filename)
        data_file[data_file=="Medium"] = "DMSO"
        begin_measure = 4
        # Indicate where the conditions are
        conditions = c(1, 2)
        data.values = data_file[, colnames(data_file) %in% model_structure$names]
    } else if (grepl(".csv$", data_filename)) {
        use_midas = TRUE
        data_file = read.delim(data_filename, sep=",")
        begin_measure = which(grepl("^DA.", colnames(data_file)))
        data_file = data_file[-begin_measure] # Delete the DA field which is not used
        begin_measure = begin_measure[1] # If there were several DA fields
        # Indicate where the conditions are
        conditions = 2:(begin_measure-1)

        # Extract the measurements of nodes in the network
        data.values = data_file[,grepl("^DV", colnames(data_file))]
        colnames(data_file) = gsub("^[A-Z]{2}.", "", colnames(data_file))
        colnames(data.values) = gsub("^[A-Z]{2}.", "", colnames(data.values))
        data.values = data.values[, colnames(data.values) %in% model_structure$names]
    }

    # Means of basal activity of the network and of the blank fixation of the antibodies
    unstim.values = colMeans(data.values[data_file$type=="c",])
    lapply(unstim.values, function(X) { if (is.nan(X)|is.na(X)) {stop("Unstimulated data are recquired to simulate the network")}})
    blank.values = colMeans(data.values[data_file$type=="blank",])
    blank.values[is.nan(blank.values)] = 0; # For conditions without blank values

    # Calculates the mean and standard deviation for each condition 
    mean.values = aggregate(as.list(data.values),by=data_file[,1:(begin_measure-1)],mean)
    sd.values = aggregate(as.list(data.values),by=data_file[,1:(begin_measure-1)],sd)
    print("Data used :")
    print(mean.values)

    # Separate values and perturbation
    data.stim = mean.values[mean.values$type=="t",begin_measure:dim(mean.values)[2]]
    data.perturb = mean.values[mean.values$type=="t", conditions]

    ### CALCULATE ERROR MODEL
    if (grepl("\\.cv$", var_file) || grepl("\\.var$", var_file)) {
        # We use the CV file if there is one
        # The format and the order of the conditions are assumed to be the same as the data file
        print("Using var file")
        if (use_midas) {
            variation.file = read.delim(var_file, sep=",")
            pre_cv = variation.file[, grepl("^DV", colnames(variation.file))]
            colnames(pre_cv) = gsub("^[A-Z]{2}.", "", colnames(pre_cv))
            cv.values = aggregate(as.list( pre_cv[colnames(pre_cv) %in% model_structure$names] ), by=data_file[,1:(begin_measure-1)], mean)
        } else {
            variation.file = read.delim(var_file)
            cv.values = aggregate(as.list(variation.file[, colnames(data.values) %in% model_structure$names]),by=data_file[,1:(begin_measure-1)],mean)
        }
        cv.stim = cv.values[cv.values$type=="t", begin_measure:dim(cv.values)[2]]
        error = matrix(rep(blank.values,each=dim(data.stim)[1]),nrow=dim(data.stim)[1]) + cv.stim * data.stim
    } else {
    # Define the lower and default error threshold
        min.cv=0.1;     # parameters: minimal cv=0.1
        default.cv=0.3; # parameters: default cv if there are only 2 replicates

        # Calculate error percentage
        cv.values = sd.values[begin_measure:dim(sd.values)[2]] / mean.values[begin_measure:dim(sd.values)[2]]
        # Values to close to the blank are removed because the error is not due to antibody specific binding
        cv.values[!mean.values[,begin_measure:dim(mean.values)[2]] > 2 * matrix(rep(blank.values,each=dim(mean.values)[1]), nrow=dim(mean.values)[1])] = NA
            
        # Generation of error percentage, one cv per antibody calculated using all the replicates available, default.cv if there is only two replicate to calculate the cv
        cv = colMeans(cv.values,na.rm=TRUE)
        cv[cv<min.cv] = min.cv
        cv[is.nan(cv)|is.na(cv)]=default.cv

        if (FALSE) { #"Multiline comment"
        for (i in 1:dim(cv.values)[2]) {
            count = 0
            for (j in 1:dim(cv.values)[1]) {
                if (!is.na(cv[i][j]) | !is.nan(cv[i][j])) {
                    count = count + 1
                }
            }
            if (count <= 2) {
                cv[i] = default.cv
                print("Defaulted")
            }
        }
        }

        if (verbose) {
            print("Error model :")
            for (i in 1:length(cv)) {
                print(paste(colnames(data.values)[i], " : ", cv[i]))
            }
        }
        error = matrix(rep(blank.values,each=dim(data.stim)[1]),nrow=dim(data.stim)[1])+matrix(rep(cv,each=dim(data.stim)[1]),nrow=dim(data.stim)[1])*data.stim
    }

### SET UP DATA OBJECT

    data=new(fitmodel::Data)
    data$set_unstim_data (matrix(rep(unstim.values,each=dim(data.stim)[1]),nrow=dim(data.stim)[1]))
    data$set_scale( data$unstim_data )
    data$set_stim_data( as.matrix(data.stim) )
    data$set_error( as.matrix(error ))

### EXTRACT EXPERIMENTAL DESIGN

# Extraction of stimulated, inhibited and measured nodes
    if (use_midas) {
        names = colnames(mean.values)[conditions]
        stim_names = names[grepl("[^i]$", names)]
        stim.nodes = as.character( stim_names[ stim_names %in% model_structure$names] )
        names = gsub("i$", "", names[grepl("i$", names)])
        inhib.nodes = as.character( names[names %in% model_structure$names] )
    } else {
        stim.nodes = as.character(unique(mean.values$stimulator[mean.values$stimulator %in% model_structure$names]))
        inhib.nodes = as.character(unique(mean.values$inhibitor[mean.values$inhibitor %in% model_structure$names]))
    }
    measured.nodes=colnames(data.stim)

# Inhibition and stimulation vectors for each experiment
    if (use_midas) {
        stimuli = as.matrix(data.perturb[stim.nodes])
    } else {
        stimuli=matrix(0,ncol=length(stim.nodes),nrow=dim(data.perturb)[1])
        for (i in 1:length(stim.nodes)) {
            stimuli[grepl(stim.nodes[i], data.perturb$stimulator),i]=1
        }
    }
    if (verbose) {
        print("Stimulated nodes")
        print(stim.nodes)
        print(stimuli)
    }
    if (use_midas) {
        inhibitor = as.matrix(data.perturb[paste(inhib.nodes, "i", sep="")])
    } else {
        inhibitor=matrix(0,ncol=length(inhib.nodes),nrow=dim(data.perturb)[1])
        if (length(inhib.nodes) > 0) { # Usefull for artificial networks
            for (i in 1:length(inhib.nodes)) {
                inhibitor[grepl(inhib.nodes[i], data.perturb$inhibitor),i]=1
            }
        }
    }
    if (verbose) {
        print("Inhibited nodes")
        print(inhib.nodes)
        print(inhibitor)
    }

# Experimental design
    expdes=getExperimentalDesign(model_structure,stim.nodes,inhib.nodes,measured.nodes,stimuli,inhibitor,basal_activity)

    core = list()
    core$design = expdes
    core$data = data
    core$structure = model_structure
    core$basal = basal_activity

    return(core)
}

# Build a model from a file, and import data for this model
rebuildModel <- function(model_file, data_file, var_file="") {
    if (!grep(".mra$", model_file)) {
        stop("The model file does not have the mra extension")
    }
    model = importModel(model_file)
    #links = matrix(rep(model$structure$names, 2), ncol=2)
    core = extractModelCore(model$structure, model$basal, data_file, var_file)

    model$data = core$data
    model$bestfit = model$model$fitmodel(model$data, model$parameters)

    return(model)
}

# Initialise the parameters with a one core processing
# needs corrections, not used anyway
classic_initialisation <- function(model, data, samples) {
    for (i in 1:nb_samples) {
        result = model$fitmodel( data, samples[i,] )
        residuals = c(residuals,result$residuals)
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
    return(results)
}

# Parallel initialisation of the parameters
parallel_initialisation <- function(model, expdes, structure, data, samples, NB_CORES) {
    # Put it under list format, as randomLHS only provides a matrix
    parallel_sampling = list()
    for (i in 1:dim(samples)[1]) {
        parallel_sampling[[i]] = samples[i,]
    }
    # Parallel initialisations
    # The number of cores used depends on the ability of the detectCores function to detect them
    parallel_results = mclapply(parallel_sampling, function(params, data, model) { model$fitmodel(data, params) }, data, model, mc.cores=NB_CORES)

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
