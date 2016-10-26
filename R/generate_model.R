getModelStructure <- function(links, struct_name="") {
  ModelStructure <- STASNet:::ModelStructure
  modelStructure =new(ModelStructure, as.character(links[,1]),as.character(links[,2]), struct_name)
  return(modelStructure);
}

getExperimentalDesign <- function(model.structure, stim.nodes, inhib.nodes, measured.nodes, stimuli, inhibitor, basal.activity = c()) {
    # Consistency checks
  if ( sum(is.na(match(as.character(stim.nodes),model.structure$names))) > 0 ) {
    stop(paste0("problem matching stim.nodes names, could not match ", stim.nodes[!stim.nodes%in%model.structure$names], collapse="\n") );
  }
  if ( sum(is.na(match(as.character(inhib.nodes),model.structure$names))) > 0 ) {
    stop(paste0("problem matching inhib.nodes names, could not match ", inhib.nodes[!inhib.nodes%in%model.structure$names], collapse="\n") );
  }
  if ( sum(is.na(match(as.character(measured.nodes),model.structure$names))) > 0 ) {
    stop(paste0("problem matching measured.nodes names, could not match ", measured.nodes[!measured.nodes%in%model.structure$names], collapse="\n") );
  }
  if ( sum(is.na(match(as.character(basal.activity),model.structure$names))) > 0 ) {
    message("Unmatched basal names :")
    message(paste0(basal.activity[is.na(match(as.character(basal.activity),model.structure$names))], collapse=" "))
  }
  if (length(basal.activity) != 0 && sum(basal.activity %in% model.structure$names) == 0 ) {
      warning("No basal names were matched!")
  }
  basal.activity = as.character(basal.activity[as.character(basal.activity) %in% model.structure$names])
  if ( !is.null(dim(stimuli)) && !is.null(dim(inhibitor)) &&  dim(stimuli)[1]!=dim(inhibitor)[1] ) {
    stop("number of experiments does not match in stimuli / inhibitor");
  } else if ( !is.null(dim(stimuli)) && dim(stimuli)[2]!=length(c(stim.nodes))) {
    stop("number stimuli does not match in stimuli / stim.nodes");
  } else if ( !is.null(dim(inhibitor)) && dim(inhibitor)[2]!=length(c(inhib.nodes)) ) {
    stop("number inhibitors does not match in inhib.nodes / inhibitor");
  }

  ExperimentalDesign <- STASNet:::ExperimentalDesign;
  expdes=new(ExperimentalDesign);

  expdes$stim_nodes=match(as.character(stim.nodes),model.structure$names)-1;
  expdes$inhib_nodes=match(as.character(inhib.nodes),model.structure$names)-1;
  expdes$measured_nodes=match(as.character(measured.nodes),model.structure$names)-1;
  expdes$basal_activity=model.structure$names %in% as.character(basal.activity);
  
  expdes$set_inhibitor(inhibitor);
  expdes$set_stimuli(stimuli);

  return(expdes)

}

