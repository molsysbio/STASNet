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
    links_list = addLinks( links_list, inputs, network[[1]], p_link[1])
    links_list = addLinks( links_list, network[[1]], network[[1]], p_layer_link[1], FALSE )
    for (ll in 2:nlayers) {
        links_list = addLinks( links_list, network[[ll-1]], network[[ll]], p_link[ll] )
        links_list = addLinks( links_list, network[[ll]], network[[ll]], p_layer_link[ll], FALSE )
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
addLinks <- function(links, lup, ldown, p_link=0.2, connect_bottom = TRUE) {
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

