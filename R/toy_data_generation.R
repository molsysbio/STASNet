############################ fake_data_generation.R ################################
# Generate fake data for fake models for testing

#' Generate a fake network tab file
#'
#' @param ninputs An integer, the number of inputs in the model
#' @param p_link Number or vector of number between 0 and 1. Probability of linking of each node in a layer with each node in the upper layer.
#' @param layers A vector of integers, the size of the vector is the number of layers besides the receptor layer, and the integers are the number of nodes in each layer
#' @param nfeedbacks An integer, the number of feedback (link from a bottom layer to an upper layer)
#' @param p_layer_link Number or vector of number. Probability of two nodes in a layer being directly linked
#' @param fname A filename where the adjacency list should be saved, "" for no saving.
#' @return A 2-columns matrix representing the network as an adjacency list
#' @export
#' @family toy data generation
generateToyNetwork <- function(ninputs, p_link=0.2, layers=c(3, 3), nfeedbacks=0, p_layer_link=0, fname="") {
    nlayers = length(layers)
    if (is.numeric(p_link)) { p_link = rep(p_link, nlayers) }
    if (is.numeric(p_layer_link)) { p_layer_link = rep(p_layer_link, nlayers) }
    if (ninputs > 26 || any(layers > 26)) { stop("The number of nodes per layer is limited to 26") }

    inputs = paste0("0", LETTERS[1:ninputs])
    links_list = cbind(paste0("I", LETTERS[1:ninputs]), inputs)
    network = list()
    network[[nlayers]] = c()
    for (ii in 1:nlayers) {
        network[[ii]] = paste0(ii, LETTERS[1:layers[ii]])
    }
    links_list = addLayerLinks( links_list, inputs, network[[1]], p_link[1])
    links_list = addLayerLinks( links_list, network[[1]], network[[1]], p_layer_link[1], FALSE )
    for (ll in 2:nlayers) {
        links_list = addLayerLinks( links_list, network[[ll-1]], network[[ll]], p_link[ll] )
        links_list = addLayerLinks( links_list, network[[ll]], network[[ll]], p_layer_link[ll], FALSE )
    }
    # Add feedbacks
    if (nfeedbacks > 0) {
        for (ii in 1:nfeedbacks) {
            sl = sample(0:nlayers, 2)
            sl = sort(sl)
            l1 = network[[ sl[2] ]]
            if (sl[1] == 0) {
                l2 = inputs
            } else {
                l2 = network[[ sl[1] ]]
            }
            links_list = rbind( links_list, c( l1[sample(1:length(l1),1)], l2[sample(1:length(l2),1)] ) )
        }
    }


    # Save the network in a file
    if (fname != "") {
        if (!grepl(".tab$", fname)) {
            fname = paste0(fname, ".tab")
        }
        write.table(links_list, fname, FALSE, FALSE, sep=" ", row.names=FALSE, col.names=FALSE)
    }

    return(links_list)
}

#' Generate a list of links between two layers with a probability of link between two nodes of p_link, with possibly at least one link to each node of the bottom layer
#'
#' @param links Initial list of links
#' @param lup List of nodes of the upper layer
#' @param ldown List of nodes of the lower layer
#' @param p_link Probability of linking two nodes between the layer
#' @param connect_bottom Whether each node of the lower layer must have at least on link
#' @family toy data generation
addLayerLinks <- function(links, lup, ldown, p_link=0.2, connect_bottom = TRUE) {
    if (p_link > 0) {
        for (ndown in ldown) {
            linkup = FALSE
            for (nup in lup) {
                if (runif(1, 0, 1) < p_link) {
                    links = rbind( links, c(nup, ndown) )
                    linkup = TRUE
                }
            }
            if (!linkup && connect_bottom) {
                links = rbind( links, c( lup[sample(1:length(lup), 1)], ndown ) )
            }
        }
    }
    return(links)
}

#' Generate an experimental design for a toy network
#'
#' Generate an experimental design for a toy network, with all stimulus applied, some nodes inhibited and some nodes measured
#' @param network A 2-columns matrix. The toy network as an adjacency list
#' @param nmes Number of measured nodes
#' @param ninh Number of inhibited nodes
#' @param stim_combo Number of simultaneous stimulations
#' @param inhib_combo Number of simultaneous inhibitions
#' @return An experimental design object
#' @seealso \code{\link{getCombinationMatrix}}, \code{\link{generateToyNetwork}}
#' @export
#' @family toy data generation
# TODO @param clever Optimize measurements and inhibitions in order to maximise the number of identifiable nodes
generateToyDesign <- function(network, nmes=4, ninh=2, stim_combo=1, inhib_combo=1) {
    clever = FALSE

    all_nodes = unique(as.vector(network))
    inputs = all_nodes[grepl("^I.", all_nodes)]
    nodes = all_nodes[!grepl("^I", all_nodes)]

    structure = getModelStructure(network)

    if (!clever) {
        measured = sample(nodes, nmes)
        inhibited = sample(nodes, ninh)
        stimulated = inputs
    }

    perturbations = getCombinationMatrix( c(stimulated, paste0(inhibited, "i")), inhib_combo, stim_combo )
    inh_nodes = grepl("i$", colnames(perturbations))
    inhibitions = perturbations[,inh_nodes]
    stimulations = perturbations[,!inh_nodes]

    design = getExperimentalDesign(structure, stimulated, inhibited, measured, stimulations, inhibitions, nodes)
    return(design)
}

#' Get an adjacency matrix from a file
#'
#' Get the adjacency matrix from a various formats of network files
#' @param network_file The name of the file containing a network structure as an adjacency matrix, or an adjacency list.
#' @return The weighted adjacency matrix corresponding to the network
readNetworkAdj <- function(network_file) {
    if (!is.matrix(network_file)) {
        #network_file = as.matrix(read.csv(network_file, header=F))
        network_split = strsplit(readLines(network_file), ",|->|;|\\ |\t")
        network_file = t(sapply(network_split, function(X){X}))
        # Adjacency list
        if (ncol(network_file) <= 3) {
            values = rep(1, nrow(network_file))
            if (ncol(network_file) == 3) {
                values = as.numeric(network_file[,3])
                network_file = network_file[,1:2]
            }
            nodes = unique(as.character(network_file))
            nnodes = length(nodes)
            adm = matrix(0, ncol=nnodes, nrow=nnodes, dimnames=list(nodes, nodes))
            for (ii in 1:nnodes) { adm[ii,ii]=-1 }
            for (rr in 1:nrow(network_file)) {
                adm[network_file[rr,1],network_file[rr,2]] = values[rr]
            }
        } else { # Adjacency matrix
            if ( any(is.na(suppressWarnings(as.numeric(network_file)))) ) {
                colnames(network_file) = network_file[1,]
                network_file = matrix( as.numeric(network_file[-1,]), ncol=ncol(network_file), dimnames=list(NULL, colnames(network_file)) )
            } else if (all( colnames(network_file) == paste0("V", 1:ncol(network_file)) )) {
                colnames(network_file) = NULL
            }
            if (ncol(network_file) != nrow(network_file)) {
                stop("The adjacency matrix has incorrect dimensions, number of lines and columns do not match")
            }
            rownames(network_file) = colnames(network_file)
            adm = matrix( as.numeric(network_file), ncol=ncol(network_file), dimnames=list(NULL, colnames(network_file)) )
        }
    } else {
        adm = network_file
    }
    invisible(adm)
}

#' Create simulated data
#'
#' Simulate data with noise from a network and a perturbation scheme
#' @param input_network The input network. A 2 or 3 columns matrix representing an adjacency list, or an adjacency matrix. Alternatively, the name of a csv file where such matrix is writen (without headers for an adjacency list, with the nodes names in the first line for an adjacency matrix). In case of a 2 columns adjacency list, the values of all links is assumed to be 1, if 3 columns, the 3rd columns is used as coefficient.
#' @param perturbations The experimental design. A matrix where the column names are the names of the perturbation: NODE for a stimulation and NODEi for an inhibition. Alternatively, the name of csv file where the matrix is writen, the first line is used as column names.
#' @param measured A list of nodes to be the measured nodes
#' @param inhibitions A single value, a list of values or NA. Values in ]0, -inf] to use for the inhibition, representing the log2-fold change in activity of the node (alternatively, a value between 0 and 1 representing the fraction of activity remaining after inhibition compared to basal). If NA, the values fitted for the inhibition will be used, or -1 if an inhibition is requested for a node that was not inhibited in the experiment.
#' @param noise A numeric between 0 and 1 giving the noise level, interpreted as a coefficient of variation. 0 means no noise at all.
#' @param replicates Number of replicates to simulate
#' @export
#' @seealso \code{\link{simulateModel}}
#' @family toy data generation
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
createSimulation <- function(input_network, perturbations="", measured="", inhibitions=0.5, noise=0, replicates=3) {
    adm = readNetworkAdj(input_network)
    structure = extractStructure(input_network)
    input_nodes = which(apply(structure$adjacencyMatrix, 1, sum)==0)
    if (length(input_nodes)==0) { input_nodes = 1 }
    output_nodes = which(apply(structure$adjacencyMatrix, 2, sum)==0)

    # If perturbations or measured are not set, stimulate all input nodes,  measure all the others, and inhibit all the others but the output nodes
    if (!is.matrix(perturbations)) {
        if (perturbations == "") {
            perturbations = getCombinationMatrix( c(structure$names[input_nodes], paste0(structure$names[-c(input_nodes, output_nodes)], "i")), 1, 1 )
        } else if (length(perturbations) == 1) {
            pfile = perturbations
            perturbations = as.matrix(read.csv(perturbations))
            if (nrow(perturbations) == 1 || ncol(perturbations) == 1) {
                perturbations = getCombinationMatrix(unlist(strsplit(readLines(pfile), "\t| |,|;")), 1, 1)
            }
        } else {
            perturbations = getCombinationMatrix(perturbations, 1, 1)
        }
    }

    if (length(measured) == 1) {
        if (measured == "") {
            measured = structure$names[-input_nodes]
        } else {
            measured = unlist(strsplit( readLines(measured), ",|->|;|\\ |\t" ))
        }
    }
    measured = unlist(measured)

    # Extract stimuli and inhibitors for the experimental design
    stim_nodes = colnames(perturbations)[grepl("[^i]$", colnames(perturbations))]
    stimuli = matrix(perturbations[,stim_nodes], ncol=length(stim_nodes), dimnames=list(NULL, stim_nodes))
    inhib_nodes = sub("i$", "", colnames(perturbations)[grepl("i$", colnames(perturbations))])
    inhibitors = matrix(perturbations[,paste0(inhib_nodes, "i")], ncol=length(inhib_nodes), dimnames=list(NULL, paste0(inhib_nodes, "i")))
    expdes = getExperimentalDesign(structure, stim_nodes, inhib_nodes, measured, stimuli, inhibitors, structure$names)
    if (length(inhibitions) != length(inhib_nodes)) {
        if (length(inhibitions) != 1) { warning("Incorrect number of inhibitions, the first one will be repeated") }
        inhibitions = rep(inhibitions[1], length(inhib_nodes))
    }

    data = new(STASNet:::Data)
    data$set_unstim_data( matrix( rep(10, nrow(stimuli)*length(measured)), nrow=nrow(stimuli) ))
    data$set_scale( matrix( rep(0, nrow(stimuli)*length(measured)), nrow=nrow(stimuli) ))

    model = new(STASNet:::Model)
    model$setModel(expdes, structure)
    mra_model = MRAmodel(model, expdes, structure, basal=structure$names, parameters=model$getParameterFromLocalResponse(adm, inhibitions), data=data)

    simulated_data = simulateModel( mra_model, inhibition_effect=inhibitions )$bestfit
    colnames(simulated_data) = paste0("DV:", colnames(simulated_data))
    colnames(stimuli) = paste0("TR:", colnames(stimuli))
    colnames(inhibitors) = paste0("TR:", colnames(inhibitors))
    results = list(model=mra_model, simulation=simulated_data, noise_free_simulation=simulated_data)
    results$noise_free_simulation = cbind(stimuli, inhibitors, results$noise_free_simulation)
    if (noise > 0) {
        replicates = floor(replicates)
        if (replicates < 1) { replicates=1 }
        if (replicates > 1) { # Simulate replicate experiments
            stimuli = apply(stimuli, 2, rep, replicates)
            simulated_data = apply(simulated_data, 2, rep, replicates)
            inhibitors = apply(inhibitors, 2, rep, replicates)
        }
        noise_coef = rnorm(length(simulated_data), sd=noise)
        noise_coef[noise_coef > 0.99] = 0.99
        noise_coef[noise_coef < -0.99] = -0.99
        noisy_simulation = simulated_data * (1 + noise_coef)
        results$simulation = cbind(stimuli, inhibitors, noisy_simulation)
        results$simulation = results$simulation[order(matrix( 1:(nrow(results$simulation)), ncol=replicates, byrow=T )),] # Put replicates next to each other

    } else {
        results$simulation = results$noise_free_simulation
    }

    return(results)
}

