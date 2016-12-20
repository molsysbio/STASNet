########################### import_export.R ##########################
# Functions use to import and export models

# Custom paste functions
pastecoma <- function(...) {
    return(paste(..., sep=",", collapse=","))
}
pastetab <- function(...) {
    return(paste(..., sep="\t", collapse="\t"))
}

# Helper function to determine if an object is a string
is.string <- function(x) {
    if (is.character(x) && length(x) == 1) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#' Exports the model in a file 
#'
#' Export an MRAmodel object in a .mra file
#' @param model_description An MRAmodel object
#' @param file_name Name of the output file
#' @param export_data Whether the CVs should be incorporated to the .mra file or not
#' @return Nothing
#' @export
#' @seealso importModel, rebuildModel
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
exportModel <- function(model_description, file_name="mra_model", export_data=FALSE) {
    # Add an extension
    if (!grepl(".mra$", file_name)) {
        file_name = paste0(file_name, ".mra")
    }
    # For the controls
    nb_params = length(model_description$parameters)

    handle = file(file_name, open="w")
    # Header of the model, with its name and extra infos
    writeLines(paste0("H ", model_description$name), handle)
    for (info in model_description$infos) {
        writeLines(paste0("H ", info), handle)
    }
    
    # Bestfit and bestfitscore
    writeLines(paste0("BF ", model_description$bestfit), handle)
    writeLines(paste0("BFS ", model_description$bestfitscore), handle)
    writeLines(paste0("RS ", paste(model_description$Rscores, collapse=" ")), handle)
    # Unused readouts and perturbations
    if (length(model_description$unused_perturbations) > 0) {
        writeLines(paste0("UP ", paste(model_description$unused_perturbations, collapse=" ")), handle)
    }
    if (length(model_description$unused_readouts) > 0) {
        writeLines(paste0("UR ", paste(model_description$unused_readouts, collapse=" ")), handle)
    }
    writeLines(paste0("MCV ", model_description$min_cv), handle)
    writeLines(paste0("DCV ", model_description$default_cv), handle)

    # Names of the nodes, with info on basal activity
    for (name in model_description$structure$names) {
        line = paste0("N ", name ," ")
        if (name %in% model_description$basal) {
            line = paste(line, 1)
        } else {
            line = paste(line, 0)
        }
        writeLines(line, handle)
    }

    # Copy adjacency matrix for the structure of the network
    adj = model_description$structure$adjacencyMatrix
    for (r in 1:nrow(adj)) {
        line = paste0(adj[r,], collapse=" ")
        line = paste0("M ", line)
        writeLines(line, handle)
    }

    # Write the values of the parameters and the extreme sets
    for (i in 1:length(model_description$parameters)) {
        line = paste0("P ", model_description$parameters[i])
        # Write the range provided by the profile likelihood if both limits are there
        if (length(model_description$lower_values) == nb_params && length(model_description$upper_values) == nb_params) {
            line = paste(line, model_description$lower_values[i], model_description$upper_values[i], sep=" ")
        }
        writeLines(line, handle)
        # Write the parameters sets provided by the profile likelihood if it is there
        if (length(model_description$param_range) == length(model_description$parameters) && length(model_description$param_range[[i]]) == 2) {
            if (!is.na(model_description$param_range[[i]]$low_set[1])) {
                line = paste0(model_description$param_range[[i]]$low_set, collapse=" ")
            } else {
                line = "NA"
            }
            line = paste0("PL ", line)
            writeLines(line, handle)
            if (!is.na(model_description$param_range[[i]]$high_set[1])) {
                line = paste0(model_description$param_range[[i]]$high_set, collapse=" ")
            } else {
                line = "NA"
            }
            line = paste0("PU ", line)
            writeLines(line, handle)
        }
    }
    
    # Write the experimental design
    design = model_description$design
    ## Inhibitions
    line = paste0(model_description$structure$names[1 + design$inhib_nodes], collapse = " ")
    writeLines(paste0("IN ", line) , handle)
    for (r in 1:nrow(design$inhibitor)) {
        line = paste0(design$inhibitor[r,], collapse = " ")
        writeLines(paste0("I ", line), handle)
    }
    ## Stimulations
    line = paste0(model_description$structure$names[1 + design$stim_nodes], collapse = " ")
    writeLines(paste0("SN ", line) , handle)
    for (r in 1:nrow(design$stimuli)) {
        line = paste0(design$stimuli[r,], collapse = " ")
        writeLines(paste0("S ", line), handle)
    }
    # Measured nodes with the unstimulated value
    for (i in 1:length(design$measured_nodes)) {
        writeLines(paste0("MN ", model_description$structure$names[1 + design$measured_nodes[i]], " ", model_description$data$unstim_data[1, i]), handle)
    }

    # Write the raw data used for the model in the file
    if (export_data) {
        # CV matrix
        if (length(model_description$cv) > 1) {
            for (r in 1:nrow(model_description$cv)) {
                writeLines(paste0("CV ", paste0(model_description$cv[r,], collapse=" ")), handle)
            }
        }
    }

    close(handle)
}

#' Import model from a file
#'
#' Import an MRAmodel object from a .mra file
#' @param file_name Name of the .mra file
#' @return An MRAmodel object
#' @export
#' @seealso exportModel, rebuildModel
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
importModel <- function(file_name) {

    if (!grepl(".mra", file_name)) {
        warning("This file does not have the expected .mra extension. Trying to extract a model anyway...")
    }

    file = readLines(file_name)
    lnb = 1
    if (!grepl("^[NH]", file[lnb])) {
        stop("This is not a valid mra model file.")
    }

    # Model name (cell line, network, ...)
    infos = ""
    if (grepl("^H", file[lnb])) {
        name = gsub("^H[A-Z]?( |\t)", "", file[lnb])
        lnb = lnb + 1
    } else {
        name = ""
    }
    while (grepl("^H", file[lnb])) {
        infos = c(infos, gsub("^H[A-Z]?( |\t)", "", file[lnb]))
        lnb = lnb + 1
    }
    # Model fitting residual and score
    bestfit = NA
    if (grepl("^BF", file[lnb])) {
        bestfit = as.numeric(gsub("^BF( |\t)", "", file[lnb]))
        lnb = lnb + 1
    }
    bestfitscore = NA
    if (grepl("^BFS", file[lnb])) {
        bestfitscore = as.numeric(gsub("^BFS( |\t)", "", file[lnb]))
        lnb = lnb + 1
    }
    Rscores = NA
    if (grepl("^RS", file[lnb])) {
        Rscores = gsub("^RS( |\t)", "", file[lnb])
        Rscores = as.numeric( unlist(strsplit(Rscores, " ")) )
        lnb = lnb + 1
    }
    # Model fitting modification performed on the data matrix
    unused_perturbations = c()
    if (grepl("^UP", file[lnb])) {
        unused_perturbations = gsub("^UP( |\t)", "", file[lnb])
        unused_perturbations = unlist(strsplit(unused_perturbations, " |\t"))
        lnb = lnb + 1
    }
    unused_readouts = c()
    if (grepl("^UR", file[lnb])) {
        unused_readouts = gsub("^UR( |\t)", "", file[lnb])
        unused_readouts = unlist(strsplit(unused_perturbations, " |\t"))
        lnb = lnb + 1
    }
    # CV settings
    min_cv = 0.1
    if (grepl("^MCV", file[lnb])) {
        min_cv = as.numeric(gsub("^MCV ", "", file[lnb]))
        lnb = lnb + 1
    }
    default_cv = 0.3
    if (grepl("^DCV", file[lnb])) {
        default_cv = as.numeric(gsub("^DCV ", "", file[lnb]))
        lnb = lnb + 1
    }

    # Names of the nodes of the model and nodes with basal activity
    nodes = c()
    basal_activity = c()
    while (grepl("^N", file[lnb])) {
        line = unlist(strsplit(file[lnb], "( |\t)+"))
        nodes = c(nodes, line[2])
        if ( as.numeric(line[3]) == 1 ) {
            basal_activity = c(basal_activity, line[2])
        }
        lnb = lnb + 1
    }
    basal = basal_activity

    # Get the format of the network and put it in a modelStructure object
    links_matrix = c()
    links_list = c()
    # Adjacency matrix
    if (grepl("^M", file[lnb])) {
        while (grepl("^M", file[lnb])) {
            line = unlist(strsplit(file[lnb], ",| |\t|;"))
            links_matrix = rbind(links_matrix, as.numeric( line[2:length(line)] ))
            lnb=lnb+1
        }
    }
#    else if (grepl("[Aa](djacency|DJACENCY)", file[lnb])) { # Adjacency list
#        lnb = lnb + 1
#        message("Adjacency list not implemented yet")
#    } else { # List of links
#        if (grepl("list|LIST", file[lnb])) { lnb = lnb + 1; }
#        message("Link list not implemented yet")
#        while (grepl("^L", file[lnb])) {
#            line = gsub("^L( |\t)", "", file[lnb])
#            line = unlist(strsplit(line, " |\t"))
#            links_matrix = rbind(links_matrix, c(line[1], line[2]))

#            lnb = lnb + 1
#        }
#    }
    # Convert from the matrix form to the link list form to get the model structure
    for (r in 1:nrow(links_matrix)) {
        for (c in 1:ncol(links_matrix)) {
            if (links_matrix[r, c] != 0) {
                links_list = rbind(links_list, c(nodes[c], nodes[r]))
            }
        }
    }
    structure = getModelStructure(links_list)

    if (!grepl("^P", file[lnb])) {
        stop("This mra file is not valid, the parameters lines should start with a P")
    }
    # Collect the values of the parameters with the limit cases provided by the profile likelihood if available
    parameters = c()
    lower_values = c()
    upper_values = c()
    param_range = list()
    id = 1
    while (grepl("^P", file[lnb])) {
        line = unlist(strsplit(file[lnb], " +|\t|;")) # PV fitted_value lower_value upper_value
        parameters = c( parameters, as.numeric(line[2]) )
        if (length(line) > 2) {
            # NA will be introduced if there is no limit
            lower_values = c(lower_values, suppressWarnings(as.numeric(line[3])) )
            upper_values = c(upper_values, suppressWarnings(as.numeric(line[4])) )
        }
        lnb = lnb + 1
        # Parameters sets for the extreme values of the confidence interval for the parameter
        if (grepl("^PL|^PU", file[lnb])) {
            param_range[[id]] = list()
        }
        while (grepl("^PL|^PU", file[lnb])) {
            line = unlist(strsplit(file[lnb], " +|\t|;")) # One parameter set
            line = suppressWarnings(as.numeric( line[2:length(line)] ))
            if (grepl("^PL", file[lnb])) {
                param_range[[id]]$low_set = line
            }
            if (grepl("^PU", file[lnb])) {
                param_range[[id]]$high_set = line
            }
            lnb = lnb + 1
        }
        id = id + 1; # Parameter index
    }

    # Collect the experimental design to build the equations
    ## List of inhibited nodes by C++ index
    if (grepl("^IN", file[lnb])) {
        line = unlist(strsplit(file[lnb], " +|\t|;"))
        inhib_nodes = line[2:length(line)]
        lnb = lnb + 1
    } else {
        stop("This mra file is not valid, the experimental design lines should be Inhibited Nodes, Inhibition matrix, Stimulated nodes, Stimulation matrix, Measured nodes (with unstimulated value)")
    }
    ## Inhibitions for each measurement
    inhibitions = c()
    while (grepl("^I", file[lnb])) {
        line = unlist(strsplit(file[lnb], " +|\t|;"))
        line = as.numeric(line[2:length(line)])
        lnb = lnb + 1

        inhibitions = rbind(inhibitions, line)
    }
    ## Index of the inhibited nodes
    if (grepl("^SN", file[lnb])) {
        line = unlist(strsplit(file[lnb], " +|\t|;"))
        stim_nodes = line[2:length(line)]
        lnb = lnb + 1
    }
    ## Stimuli for each measurement
    stimuli = c()
    while (grepl("^S", file[lnb])) {
        line = unlist(strsplit(file[lnb], " +|\t|;"))
        line = as.numeric(line[2:length(line)])
        lnb = lnb + 1

        stimuli = rbind(stimuli, line)
    }
    ## Measured nodes with their unstimulated value
    unstim_data = c()
    measured_nodes = c()
    while (grepl("^MN", file[lnb])) {
        line = unlist(strsplit(file[lnb], " +|\t|;"))
        lnb = lnb + 1

        measured_nodes = c(measured_nodes, line[2])
        unstim_data = c(unstim_data, as.numeric(line[3]))
    }

    # Set up the experimental design and the model
    expDes = getExperimentalDesign(structure, stim_nodes, inhib_nodes, measured_nodes, stimuli, inhibitions, basal_activity)
    design = expDes
    model = new(STASNet:::Model)
    model$setModel( expDes, structure )

    # Get the unstimulated data
    data = new(STASNet:::Data)
    data$set_unstim_data( matrix(rep(unstim_data, each = nrow(stimuli)), nrow = nrow(stimuli)) )

    # Get the cv values if they are present
    cv_values = c()
    while(grepl("^CV", file[lnb])) {
        line = unlist(strsplit(file[lnb], " +|\t|;"))
        lnb = lnb + 1

        cv_values = rbind(cv_values, line[2:length(line)])
    }
    cv = cv_values
# TODO import the data, and calculate the base fit

    model_description = MRAmodel(model, design, structure, basal, data, cv, parameters, bestfit, name, infos, param_range, lower_values, upper_values, unused_perturbations, unused_readouts, min_cv, default_cv)
    model_description$bestfitscore = bestfitscore
    meas_nodes = getMeasuredNodesNames(model_description)
    if ( length(Rscores) == length(meas_nodes) ) { names(Rscores) = meas_nodes }
    model_description$Rscores = Rscores
    return(model_description)
}

#' Read a MIDAS file
#' 
#' Read a MIDAS file and returns it as a matrix containing the experiments as rows and the readouts as colums
#' @param fname Name of the MIDAS file
#' @export
readMIDAS <- function(fname) {
    data_file = read.delim(fname, sep=",")
    checkMIDAS(data_file)

    measures = data.matrix(data_file[grepl("^DV", colnames(data_file))])
    treatments = data.matrix(data_file[grepl("^TR", colnames(data_file))])
    rownames(measures) = sapply(1:nrow(measures), function(rr) { paste0(paste0("", colnames(treatments)[as.logical(treatments[rr,])]), collapse="+") })
    rownames(measures) = gsub("TR.", "", rownames(measures))
    rownames(measures)[rownames(measures)==""] = as.character(data_file[,"ID.type"])[rownames(measures)==""]
    colnames(measures) = gsub("DV.", "", colnames(measures))

    return(measures[order(rownames(measures)),])
}
checkMIDAS <- function(data_file) {
    if (!any(grepl("^DV", colnames(data_file)))) { stop("This is not a MIDAS data, the mandatory 'DV' field is missing") }
    if (!any(grepl("^ID", colnames(data_file)))) { stop("This is not a MIDAS data or the field 'ID:type' is missing") }
    if (!any(grepl("^TR", colnames(data_file)))) { stop("This is not a MIDAS data, the mandatory 'TR' field is missing") }
}

