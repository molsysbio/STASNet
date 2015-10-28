########################### import_export.R ##########################
# Functions use to import and export models

# Custom paste functions
pastecoma <- function(...) {
    return(paste(..., sep=","))
}
pastetab <- function(...) {
    return(paste(..., sep="\t"))
}

#' Exports the model in a file 
#'
#' Export an MRAmodel object in a .mra file
#' @param model_description An MRAmodel object
#' @param file_name Name of the output file
#' @return Nothing
#' @export
#' @seealso importModel, rebuildModel
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
exportModel <- function(model_description, file_name="model") {
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
    for (i in 1:length(design$measured_nodes)) {
        writeLines(paste0("MN ", model_description$structure$names[1 + design$measured_nodes[i]], " ", model_description$data$unstim_data[1, i]), handle)
    }
    # CV matrix
    for (r in 1:nrow(model_description$cv)) {
        writeLines(paste0("CV ", paste0(model_description$cv[r,], collapse=" ")), handle)
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
    model_description = MRAmodel()
# Fields not covered :
    model_description$bestfit = NA # Used to indicate that it is an imported model without data
# --------------------

    if (!grepl(".mra", file_name)) {
        warn("This file does not have the expected .mra extension. Trying to extract a model anyway...")
    }

    file = readLines(file_name)
    lnb = 1
    if (!grepl("^[NH]", file[lnb])) {
        stop("This is not a valid mra model file.")
    }
    # Model name (cell line, network, ...)
    if (grepl("^H", file[lnb])) {
        model_description$name = gsub("^H[A-Z]?( |\t)", "", file[lnb])
        lnb = lnb + 1
    } else {
        model_description$name = ""
        model_description$infos = ""
    }
    while (grepl("^H", file[lnb])) {
        model_description$infos = c(model_description$infos, gsub("^H[A-Z]?( |\t)", "", file[lnb]))
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
    model_description$basal = basal_activity

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
#        print("Adjacency list not implemented yet")
#    } else { # List of links
#        if (grepl("list|LIST", file[lnb])) { lnb = lnb + 1; }
#        print("Link list not implemented yet")
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
    model_description$structure = getModelStructure(links_list)

    if (!grepl("^P", file[lnb])) {
        stop("This mra file is not valid, the parameters lines should start with a P")
    }
    # Collect the values of the parameters with the limit cases provided by the profile likelihood if available
    model_description$parameters = c()
    model_description$lower_values = c()
    model_description$upper_values = c()
    model_description$param_range = list()
    id = 1
    while (grepl("^P", file[lnb])) {
        line = unlist(strsplit(file[lnb], " +|\t|;")) # PV fitted_value lower_value upper_value
        model_description$parameters = c( model_description$parameters, as.numeric(line[2]) )
        if (length(line) > 2) {
            # NA will be introduced if there is no limit
            model_description$lower_values = c(model_description$lower_values, suppressWarnings(as.numeric(line[3])) )
            model_description$upper_values = c(model_description$upper_values, suppressWarnings(as.numeric(line[4])) )
        }
        lnb = lnb + 1
        # Parameters sets for the extreme values of the confidence interval for the parameter
        if (grepl("^PL|^PU", file[lnb])) {
            model_description$param_range[[id]] = list()
        }
        while (grepl("^PL|^PU", file[lnb])) {
            line = unlist(strsplit(file[lnb], " +|\t|;")) # One parameter set
            line = suppressWarnings(as.numeric( line[2:length(line)] ))
            if (grepl("^PL", file[lnb])) {
                model_description$param_range[[id]]$low_set = line
            }
            if (grepl("^PU", file[lnb])) {
                model_description$param_range[[id]]$high_set = line
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
    expDes = getExperimentalDesign(model_description$structure, stim_nodes, inhib_nodes, measured_nodes, stimuli, inhibitions, basal_activity)
    model_description$design = expDes
    model_description$model = new(fitmodel:::Model)
    model_description$model$setModel( expDes, model_description$structure )

    # Get the unstimulated data
    model_description$data = new(fitmodel:::Data)
    model_description$data$set_unstim_data( matrix(rep(unstim_data, each = nrow(stimuli)), nrow = nrow(stimuli)) )

    # Get the cv values if they are present
    cv_values = c()
    while(grepl("^CV", file[lnb])) {
        line = unlist(strsplit(file[lnb], " +|\t|;"))
        lnb = lnb + 1

        cv_values = rbind(cv_values, line[2:length(line)])
    }
    model_description$cv = cv_values

    return(model_description)
}

