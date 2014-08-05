################################ create_model.R #########################################
## Functions that simplify data parsing and model creation and fitting using the package

#' @import igraph
#' @import pheatmap
#' @import parallel

#source("R/randomLHS.r"); # Latin Hypercube Sampling

# Global variable to have more outputs
verbose = FALSE
debug = TRUE

#' Creates a parameterised model from experiment files and the network structure
#' @param model_links Path to the file containing the network structure, either in matrix form or in list of links form. Extension .tab expected
#' @param data.stimulation Path to the file containing the data in MRA_MIDAS format. Extension .csv expected.
#' @param basal_file Path to the file indicating the nodes without basal activity. Extension .dat expected.
## CHECK THE IMPLEMENTATION FOR THE NO BASAL
#' @param data.variation Path to the file containing the coefficient of variation for each measurement in MRA_MIDAS format. If it is not provided, the function uses the replicates of the data.stimulation file to determine a variance per probe (i.e antibody/DNA fragment/...). Extension .var expected.
#' @param cores Number of cores that should be used for the computation
#' @param inits Number of initialisation steps which should be performed (see method for the exact meaning of this value)
#' @param init_distribution Whether the distribution of the residuals and the parameters deduced by correlation should be plotted or not
#' @param method Method to be used for the initialisation, available methods are :
#'      random : Perform a Latin Hypercube Sampling to choose \emph{inits} starting points then perform a gradient descent to find the local optimum for each of those points.
#'      correlation : Deduce some parameters from the correlation between the measurements for the target node and all of its input nodes, then perform random to find the other parameters. Recommended, very efficient for small datasets.
#'      genetic : Genetic algorithm with mutation only. \emph{inits} is the total number of points sampled.
#'      annealing : Simulated annealing and gradient descent on the best result. \emph{inits} is the maximum number of iteration without any change before the algorithm decides it reached the best value. Use not recommended.
#' @return A list describing the model and its best fit
#' @export
#' @seealso importModel, exportModel, rebuildModel
#' @author Mathurin Dorel \mail{dorel@horus.ens.fr}
#' @examples
#' model = createModel("links.tab", "data_MIDAS.csv", "basal.dat") # Produces a model for the network described in links.tab using the data in data_MIDES.csv
#' model = createModel("links.tab", "data_MIDAS.csv", "basal.dat", "variation.var") # Uses the variation from a variation file
#' model = createModel("links.tab", "data_MIDAS.csv", "basal.dat", cores = detectCores()) # Uses all cores available (with the package parallel)
#' model = createModel("links.tab", "data_MIDAS.csv", "basal.dat", inits = 1000000) # Uses more initialisations for a complex network
createModel <- function(model_links, data.stimulation, basal_file, data.variation="", cores=1, inits=1000, init_distribution=F, method="default") {

    # Creation of the model structure object
    links = read.delim(model_links, header=FALSE)
    model_structure=getModelStructure(links)
    # Plot the network in a file # TODO find a better visualisation tool
    model_graph = graph.edgelist(as.matrix(links))
    name = unlist(strsplit(model_links, "/"))
    name = name[length(name)]
    pdf(paste0( "graph_", gsub(" ", "_", gsub(".tab$", ".pdf", name)) ))
    plot.igraph(model_graph, edge.arrow.size=0.5, layout=layout.fruchterman.reingold.grid)
    dev.off()

# TODO extraction of the basal activity with different format
    basal_activity = as.character(read.delim(basal_file,header=FALSE)[,1])
    # Create a valid basal activity list from a list of the nodes that do not have basal activity
# TODO Find a day to switch
# IDEA Consider it is basal if more than 2/3 of the nodes and not basal if less than 1/3, ask otherwise
#   basal = c()
#   for (node in model_structure$names) {
#       if (!(node %in% no_basal)) {
#           basal = c(basal, node)
#       }
#   }
    core = extractModelCore(model_structure, basal_activity , data.stimulation, data.variation)
    expdes = core$design
    data = core$data
    model_structure = core$structure

    # MODEL SETUP
    model = new(fitmodel::Model)
    model$setModel(expdes, model_structure)

    # INITIAL FIT
    # Parallelized version uses all cores but one to keep control
    if (method == "default") {
        method = "correlation"
    }
    if (cores == 0) { cores = detectCores()-1 }
    print (paste("Initializing the model parameters… (", inits, " random samplings) with ", cores, " cores, method : ", method, sep=""))
    # Different sampling methods
    if (method == "correlation") {
        samples = sampleWithCorrelation(model, core, inits, plot=init_distribution, sd=2)$samples
        results = parallel_initialisation(model, expdes, data, samples, cores)
    } else if (method == "random" || method == "sample") {
        samples = qnorm(randomLHS(inits, model$nr_of_parameters()), sd=2)
        results = parallel_initialisation(model, expdes, data, samples, cores)
    } else if (method == "genetic" || method == "explore" || method == "deep") {
        results = deep_initialisation(model, core, cores, 3, inits, init_distribution)
    } else if (method == "annealing" || method == "SA") {
        results = parallelAnnealing(model, core, inits, cores, init_distribution)
    } else {
        stop("The selected initialisation method does not exist (valid methods are correlation, random, explore, and annealing)")
    }
    print("Initial fits completed")

    # Choice of the best fit
    params = results$params
    residuals = results$residuals
    #write("All residuals : ", stderr())
    #write(residuals, stderr())
    print(paste( sum(is.na(residuals)), "NA and ", sum(is.infinite(residuals)), "infinite residuals" ))
    # Refit the best residuals to make sure it converged # TODO WITH A MORE PRECISE DESCENT (more max steps and better delta-p)
    order_id = order(residuals)[1:20]
    for ( i in 1:min(20, length(residuals)) ) {
        p_res = residuals[order_id[i]]
        if (!is.na(p_res) && !is.infinite(p_res)) {
            repeat {
                result = model$fitmodel(core$data, params[order_id[i],])
                if (p_res == result$residuals) {
                    break
                }
                p_res = result$residuals
            }
            params[order_id[i],] = result$parameter
            residuals[order_id[i]] = result$residuals
        }
    }
    if (init_distribution) { hist(log(residuals, base=10), breaks="fd", main="Distribution of the residuals") }
    if (debug) {
        # Print the 20 smallest residuals to check if the optimum has been found several times
        print("Best residuals :")
        print(sort(residuals)[1:20])
    }

    best = order(residuals)[1]
    init_params = params[best,]
    init_residual = residuals[best]
    if (verbose) {
        print("Model simulation with optimal parameters :")
        print(model$simulate(data, init_params)$prediction)
    }

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
    model_description$infos = c(paste0(inits, " samplings"), paste0( sort("Best residuals : "), paste0(sort(residuals)[1:5], collapse=" ") ), paste0("Method : ", method))
    # Values that can't be defined without the profile likelihood
    model_description$param_range = list()
    model_description$lower_values = c()
    model_description$upper_values = c()

    return(model_description)
}

#' Perform an initialisation of the model 
initModel <- function(model, core, nb_inits, nb_cores=1) {
    # Parallelized version uses all cores but one to keep control
    if (cores == 0) { cores = detectCores()-1 }
    print (paste("Initializing the model parameters… (", inits, " random samplings) with ", cores, " cores", sep=""))
    # Different sampling methods
    if (method == "default") {
        method = "correlation"
    }
    if (method == "correlation") {
        samples = sampleWithCorrelation(model, core, inits, plot=init_distribution, sd=2)$samples
        results = parallel_initialisation(model, expdes, data, samples, cores)
    } else if (method == "random" || method == "sample") {
        samples = qnorm(randomLHS(inits, model$nr_of_parameters()), sd=2)
        results = parallel_initialisation(model, expdes, data, samples, cores)
    } else if (method == "genetic" || method == "explore" || method == "deep") {
        results = deep_initialisation(model, core, cores, 3, inits, init_distribution)
    } else if (method == "annealing" || method == "SA") {
        results = parallelAnnealing(model, core, inits, cores, init_distribution)
    } else {
        stop("The selected initialisation method does not exist (valid methods are correlation, random, explore, and annealing)")
    }

    return(results)
}

# Perform several simulated annealing, one per core
parallelAnnealing <- function(model, core, max_it, nb_cores, do_plot=F) {
    correlation = parametersFromCorrelation(model, core, plot=do_plot)
    correlation_list = list()
    for (i in 1:nb_cores) {
        correlation_list[[i]] = correlation
    }
    # Perform several annealing to select the best
    annealings = mclapply(correlation_list, simulated_annealing_wrapper, model, core$data, max_it, 0, mc.cores=nb_cores)
    results = list()
    results$residuals = c()
    results$params = c()
    for (annealing in annealings) {
        results$residuals = c(results$residuals, annealing$residuals)
        results$params = rbind(results$params, annealing$parameter)
    }

    return(results)
}

# Wrapper for the C++ function because Rcpp does not take into account optionnal arguments (and R crashes if there is not exactly the correct number of arguments)
# and for the parallelisation
simulated_annealing_wrapper <- function(correlated_parameters, model, data, max_it=0, max_depth=0) {
    return(model$annealingFit(data, correlated_parameters, max_it, max_depth))
}

# TODO put it in deep_initialisation
# Produces shifted random samples with some parameters inferred from correlations in the measurements
#' @param shift A vector giving the shift to be applied to the normal distribution. Must have the same size as the number of parameters
sampleWithCorrelation <- function(model, core, nb_samples, shift=0, correlated="", sd=2, plot=F) {
    if (!is.list(correlated)) {
        correlated = correlate_parameters(model, core, plot)
    }
    random_samples = qnorm(randomLHS(nb_samples, model$nr_of_parameters()-length(correlated$list)), sd=sd)
    samples = c()
    j=1
    for (i in 1:model$nr_of_parameters()) {
        if (i %in% correlated$list) {
            samples = cbind(samples, rep(correlated$values[which(correlated$list==i)], times=nb_samples))
        } else {
            if ( length(shift) == model$nr_of_parameters() ) {
                samples = cbind(samples, random_samples[,j] + shift[i])
            } else {
                samples = cbind(samples, random_samples[,j])
            }
            j = j+1
        }
    }

    cor_samples = list()
    cor_samples$samples = samples
    cor_samples$cor = correlated
    print("Sampling terminated.")
    return(cor_samples)
}

# Gives a parameter vector with a value for correlated parameters and 0 otherwise
parametersFromCorrelation <- function(model, core, plot=F) {
    correlated = correlate_parameters(model, core, plot)
    param_vector = c()
    j = 1
    for (i in 1:model$nr_of_parameters()) {
        if (i %in% correlated$list) {
            param_vector = c(param_vector, correlated$values[j])
            j = j + 1
        } else {
            param_vector = c(param_vector, 0)
        }
    }
    return(param_vector)
}

# Indentify the links that can be deduced by a simple correlation, calculate them, and return their index and the corresponding parameter vector with the other values set to 0
correlate_parameters <- function(model, core, plot=F) {
    print("Looking for identifiable correlations...")
    expdes = core$design
    model_structure = core$structure
    data = core$data
    # Collect the nodes that can be identified by correlation
    measured_nodes = expdes$measured_nodes+1;
    valid_nodes = c()
    upstreams = list()
    for (node in measured_nodes) {
        valid = TRUE
        senders = c()
        for (sender in 1:length(model_structure$names)) {
            if (model_structure$adjacencyMatrix[node, sender] == 1) {
                # All sender nodes must be measured for the input links to be indentifiable
                if (sender %in% measured_nodes) {
                    senders = c(senders, sender)
                } else {
                    valid = FALSE
                    if (verbose) { print(paste(model_structure$names[node], "excluded because", model_structure$names[sender], "is not measured")) }
                    break
                }
            }
        }
        if (valid) {
            upstreams[[node]] = senders
            valid_nodes = c(valid_nodes, node)
        }
    }
    print(paste(length(valid_nodes), "target nodes found correlatable with their inputs"))

    # Collect the values and perform the correlation
    params_matrix = matrix(0, ncol=length(model_structure$names), nrow=length(model_structure$names)) # To store the values
    for (node in valid_nodes) {
        use = rep(TRUE, times=nrow(expdes$inhibitor))
        # We do not use the conditions where one of the sender nodes is inhibited (since it decorelates the measurements)
        for ( sender in which((expdes$inhib_nodes+1) %in% upstreams[[node]]) ) {
            for ( i in 1:nrow(expdes$inhibitor) ) {
                if (expdes$inhibitor[i, sender] == 1) {
                    use[i] = FALSE
                }
            }
        }

        # Collect the measured data for the regression
        node_mes = which( measured_nodes == node)
        measurements = log(data$stim_data[use, node_mes] / data$unstim_data[1, node_mes] )
        regression = paste0("lm(measurements[,", 1, '] ~ ')
        first = TRUE
        condition = paste(model_structure$names[node], "and")
        for (sender in upstreams[[node]] ) {
            mes_index = which(measured_nodes==sender)
            measurements = cbind(measurements, log(data$stim_data[use, mes_index] / data$unstim_data[1,mes_index]) )
            if (first) {
                first = FALSE
            } else {
                regression = paste0(regression, "+")
                condition = paste(condition, "+")
            }
            regression = paste0(regression, "measurements[,", 1+which(upstreams[[node]] == sender), "]")
            condition = paste(condition, model_structure$names[sender])
        }
        # Perform the regression and put the values in the adjacency matrix
        regression = paste0(regression, ")")
        result = eval(parse(text=regression))
        condition = paste(condition, ".\n R^2 =", signif(summary(result)$r.squared, 2) )
        for (sender in 1:length(upstreams[[node]])) {
            params_matrix[ node, upstreams[[node]][sender] ] = result$coefficients[sender+1]
            # TODO add the information on the quality of the fit
        }
        if (plot) {
            # Plot the data and the fit
            if (length(upstreams[[node]]) == 1) {
                plot(measurements[,2], measurements[,1], main=condition)
                lines(measurements[,2], result$coefficients[1] + result$coefficients[2] * measurements[,2])
            } else {
                plot(1:nrow(measurements), measurements[,1], pch=4, main=condition)
                for (measure in 1:nrow(measurements)) {
                    fitted = result$coefficients[1]
                    for (sender in 2:ncol(measurements)) {
                        fitted = fitted + result$coefficients[sender] * measurements[measure, sender]
                    }
                    points(measure, fitted, col="blue", pch=20)
                    lines(rep(measure, 2), c(fitted, measurements[measure,1]), col="red")
                }
            }
        }
    }
    if (verbose) {
        print(model_structure$names)
        print(params_matrix)
    }

    # Each link identified by correlation will correspond to one parameter
    # We identify which parameter it is and return it with its value
    params_vector = model$getParameterFromLocalResponse( params_matrix, rep(0, times=ncol(expdes$inhibitor)) )
    correlated = list()
    correlated$list = c()
    correlated$values = c()
    links = model$getParametersLinks()
    for ( param in 1:length(params_vector) ) {
        if (!is.nan(params_vector[param]) && params_vector[param] != 0) {
            correlated$list = c(correlated$list, param)
            correlated$values = c(correlated$values, params_vector[param])
            print(paste0("Correlated link ", simplify_path_name(links[param]), " = ", params_vector[param]))
        }
    }
    return(correlated)
}

# Generates a random number between 0 and 1
rand <- function(decimals=4) {
    return(sample(1:10^decimals, size=1) / 10^decimals)
}

# TODO integrate it into the creation function
# Explore the parameter space with series of random initialisations and selection after each step of the best fits (mutation only genetic algorithm)
MIN_SAMPLE_SIZE = 100 # TODO Find a heuristic value depending on the number of parameters to estimate
deep_initialisation <- function (model, core, NB_CORES, depth=3, totalSamples=100, plot=F) { # TODO Chose to give the total number or the per sampling

    if (depth < 1) { depth = 1 } # We do at least one initialisation
    #if (nSamples < MIN_SAMPLE_SIZE) { nSamples = MIN_SAMPLE_SIZE }
    nSamples = ceiling(sqrt(totalSamples / depth))
    print(paste("Mutations per point =", nSamples))
    
    # The first iteration does not use any shift
    kept = matrix(0, ncol=model$nr_of_parameters())
    correlation = correlate_parameters(model, core, plot)
    for (i in 1:depth) {
        print(paste("Depth :", i))
        # Create the new samples, with a random initialisation shifted by the previous steps local minimum and recombinations of those previous best steps
        # We reduce the exploration range at each step
        samples = c()
        for (j in 1:nrow(kept)) {
            # Generate half of the new samples by simple mutation
            new_sample = sampleWithCorrelation(model, core, ceiling(nSamples * 0.5), kept[j,], correlation, sd=3/i)
            samples = rbind( samples, new_sample$samples) 
            # Recombination of the individual j with some of the others
            for (k in 1:ceiling(nSamples * 0.5)) {
                new_sample$samples[k,] = (kept[j,] + kept[sample(1:nrow(kept), size=1),]) / 2
            }
            samples = rbind( samples, new_sample$samples) 
        }
        results = parallel_initialisation(model, core$design, core$data, samples, NB_CORES)
        if (plot) { hist(log(results$residuals, base=10), breaks="fd") }
        # Keep the number of samples, and select the best parameters sets for the next iteration
        kept = results$params[order(results$residuals)[1:nSamples],]
    }

    # Finish by checking that sampling around the optimum with a high variation still finds several times the optimum
    new_sample = sampleWithCorrelation(model, core, nSamples, kept[1,], correlation, sd=3)$samples
    if (plot) { hist(log(parallel_initialisation(model, core$design, core$data, new_sample, NB_CORES)$residuals, base=10), breaks="fd") }

    return(results) # Return the last set of fit
}

# Parallel initialisation of the parameters
parallel_initialisation <- function(model, expdes, data, samples, NB_CORES) {
    # Put the samples under list format, as mclapply only take "one dimension" objects
    parallel_sampling = list()
    nb_samples = dim(samples)[1]
    length(parallel_sampling) = nb_samples # Avoid to resize the list multiple times
    for (i in 1:nb_samples) {
        parallel_sampling[[i]] = samples[i,]
    }
    # Parallel initialisations
    # The number of cores used depends on the ability of the detectCores function to detect them
    fitmodel_wrapper <- function(params, data, model) {
        init = proc.time()[3]
        result = model$fitmodel(data, params)
        if (verbose) {
            write(paste(signif(proc.time()[3]-init, 3), "s for the descent, residual =", result$residuals), stderr())
        }
        return(result)
    }

    # Since the aggregation in the end of the parallel calculations takes a lot of time for a big list, we calculate by block and take the best result of each block
    if (nb_samples > 10000) {
        print("Using the block version for the parallel initialisation")
        best_results = list()
        best_results$residuals = c()
        best_results$params = c()
        for (i in 1:(nb_samples %/% 10000)) {
            parallel_results = mclapply(parallel_sampling[ (10000*(i-1)+1):(10000*i) ], fitmodel_wrapper, data, model, mc.cores=NB_CORES)

            # Reorder the results to get the same output as the linear function
            results = list()
            results$residuals = c()
            results$params = c()
            for (entry in parallel_results) {
                results$residuals = c(results$residuals, entry$residuals)
                results$params = rbind(results$params, entry$parameter)
            }
            # Only keep the best fits
            best = order(results$residuals)[1:20]
            write(paste(sum(is.na(results$residuals)), "NAs"), stderr())
            best_results$residuals = c(best_results$residuals, results$residuals[best])
            best_results$params = rbind(best_results$params, results$params[best,])
            # We make it so that the size of best_results never execeeds 10000
            if (i %% 500 == 0) {
                best = order(best_results$residuals)[1:20]
                best_results$residuals = best_results$residuals[best]
                best_results$params = best_results$params[best,]
            }
        }
        # The last block is smaller
        if (nb_samples %% 10000 != 0) {
            parallel_results = mclapply(parallel_sampling[(nb_samples %/% 10000):nb_samples], fitmodel_wrapper, data, model, mc.cores=NB_CORES)

            results = list()
            results$residuals = c()
            results$params = c()
            for (entry in parallel_results) {
                results$residuals = c(results$residuals, entry$residuals)
                results$params = rbind(results$params, entry$parameter)
            }

            best = order(results$residuals)[1:20]
            best_results$residuals = c(best_results$residuals, results$residuals[best])
            best_results$params = rbind(best_results$params, results$params[best,])
        }
    } else {
        parallel_results = mclapply(parallel_sampling, fitmodel_wrapper, data, model, mc.cores=NB_CORES)

        # Reorder the best_results to get the same output as the linear function
        best_results = list()
        best_results$residuals = c()
        best_results$params = c()
        for (entry in parallel_results) {
            best_results$residuals = c(best_results$residuals, entry$residuals)
            best_results$params = rbind(best_results$params, entry$parameter)
        }
    }


    return(best_results)
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

# Detect the format of the structure file and extract the structure of the network
extractStructure = function(model_links, names="") {
    # Read the file, and extract the matrix corresponding to it
    structure_file = readLines(model_links)
    splitlist = strsplit(structure_file, ",|->|;|\\ |\t")
    if (length(splitlist) < 3) {
        stop("This is a trivial network, you don't need a simulation !")
    }
    struct_matrix = c()
    for (i in 1:length(splitlist)) {
        struct_matrix = suppressWarnings( rbind(struct_matrix, unlist(splitlist[[i]])) )
    }

    # Detect if it is an list of links or an adjacency matrix
    if (ncol(struct_matrix) == 2) {
        # If it is a list of links and the number of links is indicated at the beginning, remove it 
        if (suppressWarnings( is.na(as.numeric(struct_matrix[1, 1])) )) {
            struct_matrix = struct_matrix[,-1]
        }
        links_list = struct_matrix
    } else {
        # Remove the number of nodes
        if ( length(unique(struct_matrix[1,])) == 1 && struct_matrix[1, 1] != 1) {
            struct_matrix = struct_matrix[,-1]
        }
        # Set the column names if they are present, otherwise use the names from another file, and in last resort use numbers
        if (suppressWarnings( is.na(as.numeric(struct_matrix[1,1])) )) {
            colnames(struct_matrix) = struct_matrix[1,]
            struct_matrix = struct_matrix[-1,]
            if (suppressWarnings( is.na(as.numeric(struct_matrix[1,1])) )) {
                struct_matrix = struct_matrix[,-1]
            }
            rownames(struct_matrix) = colnames(struct_matrix)
        } else if (length(name) == ncol(struct_matrix)) {
            colnames(struct_matrix) = names
            rownames(struct_matrix) = names
        } else {
            print("No names were provided for the nodes, using numbers instead")
            colnames(struct_matrix) = 1:ncol(struct_matrix)
            rownames(struct_matrix) = 1:ncol(struct_matrix)
        }
        # Check that the dimensions are correct
        if ( ncol(struct_matrix) != nrow(struct_matrix) ) {
            stop("The adjacency matrix is not square. Abort...")
        }

        links_list = c();
        for (i in 1:nrow(struct_matrix)) {
            for (j in 1:ncol(struct_matrix)) {
                if (as.numeric(struct_matrix[i,j]) == 1) {
                    links_list = rbind(links_list, c(colnames(struct_matrix)[j], rownames(struct_matrix)[i]))
                }
            }
        }

    }

    # Plot the graph of the network
    model_graph = graph.edgelist(as.matrix(links_list))
    name = unlist(strsplit(model_links, "/"))
    name = name[length(name)]
    pdf(paste0("graph_", gsub(".tab$", ".pdf", name)))
    plot.igraph(model_graph, edge.arrow.size=0.5, layout=layout.fruchterman.reingold.grid)
    dev.off()
    
    model_structure=getModelStructure(links_list)

    return(model_structure)
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
        not_included = (colnames(data_file)[!(colnames(data_file) %in% model_structure$names)])[-(1:(begin_measure-1))]
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
        not_included = colnames(data.values)[!(colnames(data.values) %in% model_structure$names)]
        data.values = data.values[, colnames(data.values) %in% model_structure$names]
    }
    # Warn for the measured not that have not been found in the network
    if (length(not_included) > 0) {
        print(paste(not_included , "measurement is not in the network structure (could be a mispelling)" ))
    }

    # Means of basal activity of the network and of the blank fixation of the antibodies
    unstim.values = colMeans(data.values[data_file$type=="c",])
    lapply(unstim.values, function(X) { if (is.nan(X)|is.na(X)) {stop("Unstimulated data are required to simulate the network")}})
    blank.values = colMeans(data.values[data_file$type=="blank",])
    blank.values[is.nan(blank.values)] = 0; # For conditions without blank values

    # Calculates the mean and standard deviation for each condition 
    mean.values = aggregate(as.list(data.values),by=data_file[,1:(begin_measure-1)],mean)
    sd.values = aggregate(as.list(data.values),by=data_file[,1:(begin_measure-1)],sd)
    if (verbose) {
        print("Data used :")
        print(mean.values)
    }

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
            # Check that the number of samples is the same for the measurements and the variation
            if (nrow(pre_cv) != nrow(data_file)) { stop("Different number of experiments for the variation and the measurement files") }
            # Check that the names in the measurements file and in the variation file match
            matchError = FALSE
            for (name in colnames(pre_cv)) {
                if (!(name %in% colnames(data_file))) {
                    matchError = TRUE
                    write(paste0(name, " from the variation file is not in the measurement file"), stderr())
                }
            }
            if (matchError) { stop("Names of the variation and measurement files do not match") }
            cv.values = aggregate(as.list( pre_cv[colnames(pre_cv) %in% model_structure$names] ), by=data_file[,1:(begin_measure-1)], mean)
        } else {
            variation.file = read.delim(var_file)
            cv.values = aggregate(as.list(variation.file[, colnames(data.values) %in% model_structure$names]),by=data_file[,1:(begin_measure-1)], mean)
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
    model$bestfit = model$model$fitmodel(model$data, model$parameters)$residuals

    return(model)
}
