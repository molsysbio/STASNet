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
    writeLines(paste0("LOG ", as.character(model_description$use_log)), handle)
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
    if (length(model_description$data$stim_data) > 1) {
        for (rr in 1:nrow(model_description$data$stim_data)) {
            writeLines(paste0("SD ", paste0(model_description$data$stim_data[rr,], collapse=",")), handle)
        }
    }
    if (length(model_description$data$error) > 1) {
        for (rr in 1:nrow(model_description$data$error)) {
            writeLines(paste0("ER ", paste0(model_description$data$error[rr,], collapse=",")), handle)
        }
    }
    if (length(model_description$data$scale) > 1) {
        for (rr in 1:nrow(model_description$data$scale)) {
            writeLines(paste0("SC ", paste0(model_description$data$scale[rr,], collapse=",")), handle)
        }
    }
    # CV matrix
    if (length(model_description$cv) > 1) {
        for (r in 1:nrow(model_description$cv)) {
            writeLines(paste0("CV ", paste0(model_description$cv[r,], collapse=" ")), handle)
        }
    }

    close(handle)
}

#' Import model from a file
#'
#' Import an MRAmodel object from a .mra file or from an .mra file that has been read in with readLines.
#' @param file_name Name of the .mra file
#' @param file R object from an .mra file that has been read in by readLines
#' @return An MRAmodel object
#' @export
#' @seealso exportModel, rebuildModel
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
importModel <- function(file_name=NULL,file=NULL) {

  if (is.null(file_name)){
    if(is.null(file)){stop("no input was given either 'file_name' or 'file' required")
    }
  }else{
    if (!grepl(".mra", file_name)) {
      warning("This file does not have the expected .mra extension. Trying to extract a model anyway...")
    }
    
    file = readLines(file_name)
  }
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
    use_log = FALSE
    if (grepl("^LOG", file[lnb])) {
        use_log = gsub("^LOG( |\t)", "", file[lnb])
        use_log = as.logical( use_log )
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
        unused_readouts = unlist(strsplit(unused_readouts, " |\t"))
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
    structure = STASNet:::getModelStructure(links_list, nodes)

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
        if (length(line)>1){
        inhib_nodes = line[2:length(line)]
        }else{
        inhib_nodes = integer(0)  
        }
        lnb = lnb + 1
    } else {
        stop("This mra file is not valid, the experimental design lines should be Inhibited Nodes, Inhibition matrix, Stimulated nodes, Stimulation matrix, Measured nodes (with unstimulated value)")
    }
    ## Inhibitions for each measurement
    inhibitions = c()
    while (grepl("^I", file[lnb])) {
      line = unlist(strsplit(file[lnb], " +|\t|;"))
      if (length(line)>1){
        line = as.numeric(line[2:length(line)])
        inhibitions = rbind(inhibitions, line)
      }else if (is.null(nrow(inhibitions))){  
        inhibitions = matrix(as.numeric(NA),ncol=0,nrow=1) 
      }else{
        inhibitions = rbind(inhibitions, matrix(as.numeric(NA),ncol=0,nrow=1))  
      }
      lnb = lnb + 1
    }
    ## Index of the inhibited nodes
    if (grepl("^SN", file[lnb])) {
        line = unlist(strsplit(file[lnb], " +|\t|;"))
        if (length(line)>1){
        stim_nodes = line[2:length(line)]
        }else{
        stim_nodes = integer(0)  
        }
        lnb = lnb + 1
    }
    ## Stimuli for each measurement
    stimuli = c()
    while (grepl("^S", file[lnb])) {
      line = unlist(strsplit(file[lnb], " +|\t|;"))
      if (length(line)>1){
        line = as.numeric(line[2:length(line)])
        stimuli = rbind(stimuli, line)
      }else if (is.null(nrow(stimuli))){  
        stimuli = matrix(as.numeric(NA),ncol=0,nrow=1) 
      }else{
        stimuli = rbind(stimuli, matrix(as.numeric(NA),ncol=0,nrow=1))  
      }
      lnb = lnb + 1
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
    design = STASNet:::getExperimentalDesign(structure, stim_nodes, inhib_nodes, measured_nodes, stimuli, inhibitions, basal_activity)
    model = new(STASNet:::Model)
    model$setModel( design, structure, use_log )

    # Get the unstimulated data
    data = new(STASNet:::Data)
    data$set_unstim_data( matrix(rep(unstim_data, each = nrow(stimuli)), nrow = nrow(stimuli)) )

    # Get the data values
    stim_data = c()
    while(grepl("^SD ", file[lnb])) {
        line = unlist(strsplit(file[lnb], " +|\t|;|,"))
        lnb = lnb + 1
        stim_data = rbind( stim_data, as.numeric(line[2:length(line)]) )
    }
    if (!is.null(stim_data)) {
        data$set_stim_data(stim_data)
    }
    error = c()
    while(grepl("^ER ", file[lnb])) {
        line = unlist(strsplit(file[lnb], " +|\t|;|,"))
        lnb = lnb + 1
        error = rbind( error, as.numeric(line[2:length(line)]) )
    }
    if (!is.null(error)) {
        data$set_error(error)
    }
    scale = c()
    while(grepl("^SC ", file[lnb])) {
        line = unlist(strsplit(file[lnb], " +|\t|;|,"))
        lnb = lnb + 1
        scale = rbind( scale, as.numeric(line[2:length(line)]) )
    }
    if (!is.null(scale)) {
        data$set_scale(scale)
    }
    # Get the cv values if they are present
    cv_values = c()
    while(grepl("^CV ", file[lnb])) {
        line = unlist(strsplit(file[lnb], " +|\t|;"))
        lnb = lnb + 1

        cv_values = rbind( cv_values, as.numeric(line[2:length(line)]) )
        colnames(cv_values) = structure$names[design$measured_nodes+1]
    }
    cv = cv_values
# TODO import the data, and calculate the base fit

    model_description = MRAmodel(model, design, structure, basal, data, cv, parameters, bestfit, name, infos, param_range, lower_values, upper_values, unused_perturbations, unused_readouts, min_cv, default_cv, use_log)
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
    data_file = extractMIDAS(fname)

    measures = data.matrix(data_file[grepl("^DV", colnames(data_file))])
    treatments = data.matrix(data_file[grepl("^TR", colnames(data_file))])
    rownames(measures) = sapply(1:nrow(measures), function(rr) { paste0(paste0("", colnames(treatments)[as.logical(treatments[rr,])]), collapse="+") })
    rownames(measures) = gsub("TR.", "", rownames(measures))
    rownames(measures)[rownames(measures)==""] = as.character(data_file[,"ID.type",drop=TRUE])[rownames(measures)==""]
    colnames(measures) = gsub("DV.", "", colnames(measures))

    return(measures[order(rownames(measures)),,drop=FALSE])
}
#' @param data_file A matrix to be checked for MIDAS compliance
#' @param handler Function to use to handle the error message. Should be 'stop', 'warning', 'print' or 'message'
checkMIDAS <- function(data_file, handler=stop) {
    if (!any(grepl("^DV", colnames(data_file)))) { handler("This is not a MIDAS data, the mandatory 'DV' field is missing") }
    if (!any(grepl("^ID", colnames(data_file)))) { handler("This is not a MIDAS data or the field 'ID:type' is missing") }
    if (!any(grepl("^TR", colnames(data_file)))) { handler("This is not a MIDAS data, the mandatory 'TR' field is missing") }
}

#' Compute the log-fold change to control
#'
#' Compute the log-fold change to control, remove the blank and control lines
#' @export
controlFC <- function(idata) {
    midas_format = FALSE
    if ( any(grepl("^ID.type$", colnames(idata))) ) {
        midas_format = TRUE
        full_datas = idata[,grepl("^DV", colnames(idata))]
        control_selection = which(grepl("c$|control", idata[,"ID.type"]))
        blank_selection = which(idata[,"ID.type"]=="blank")
    } else {
        full_datas = idata
        control_selection = grep("^c$|control", rownames(full_datas))
        blank_selection = which(rownames(full_datas)=="blank")
    }
    control=full_datas[control_selection,,drop=FALSE]
    if (!is.null(nrow(control)) && nrow(control) > 0) {
        control_line = colMeans(control)
    } else {
        control_line = rep(1, ncol(full_datas))
    }
    blank=full_datas[blank_selection,]
    if (length(blank_selection) > 0 || length(control_selection) > 0) {
        datas=full_datas[-c( blank_selection, control_selection ),]
    } else {
        datas = full_datas
    }
    control = matrix(rep(control_line, nrow(datas)), nrow=nrow(datas), byrow=T)
    fcdata = log( datas / control )
    if (midas_format) {
        idata = idata[-c(control_selection, blank_selection),]
        idata[grepl("^DV", colnames(idata))] = fcdata
        fcdata = idata
    }
    return(fcdata)
}

#' Write a MIDAS file from a human readable dataset
#'
#' Write a MIDAS file from a human readable dataset
#' @param data Matrix where the rownames are the treatment separated by + signs and the column names are the readouts
#' @export
midasFromData <- function(data, fname) {
    treatments = unique(unlist(lapply( rownames(data), function(tt) { unlist(strsplit(tt, "\\+")) } )))
    control_blank = which(grepl( "control|cntrl|^c$|blank", tolower(treatments) ))
    if (length(control_blank) > 0) { treatments = treatments[-control_blank] }
    midas_data = matrix(nrow=0, ncol=2+length(treatments)+ncol(data))

    colnames(midas_data) = c( "ID:type", paste0("TR:", treatments), "DA:ALL", paste0("DV:", colnames(data)) )
    for (rr in 1:nrow(data)) {
        if (rownames(data)[rr] == "blank") {
            midas_data = rbind( midas_data, c("blank", rep(0,length(treatments)+1), data[rr,]) )
        } else if ( grepl("^c$|control|cntrl", tolower(rownames(data)[rr])) ) {
            midas_data = rbind( midas_data, c("control", rep(0,length(treatments)+1), data[rr,]) )
        } else {
            midas_data = rbind( midas_data, c("t", rep(0,length(treatments)+1), data[rr,]) )
            combination = strsplit(rownames(data)[rr], "\\+")
            for (tr in combination) {
                midas_data[rr, paste0("TR:", tr)] = 1
            }
        }
    }
    write.csv(midas_data, file=fname, row.names=FALSE, quote=FALSE)
    invisible(midas_data)
}
