getModelStructure <- function(links) {
  ModelStructure <- fitmodel:::ModelStructure
  modelStructure =new(ModelStructure, as.character(links[,1]),as.character(links[,2]))
  return(modelStructure);
}

getExperimentalDesign <- function(model.structure, stim.nodes, inhib.nodes, measured.nodes, stimuli, inhibitor, basal.activity = c()) {
    # Consistency checks
  if ( sum(is.na(match(as.character(stim.nodes),model.structure$names))) > 0 ) {
    print(stim.nodes)
    print(model.structure$names)
    stop("problem matching stim.nodes names");
  }
  if ( sum(is.na(match(as.character(inhib.nodes),model.structure$names))) > 0 ) {
    print(inhib.nodes)
    print(model.structure$names)
    stop("problem matching inhib.nodes names");
  }
  if ( sum(is.na(match(as.character(measured.nodes),model.structure$names))) > 0 ) {
    stop("problem matching measured.nodes names");
  }
  if ( sum(is.na(match(as.character(basal.activity),model.structure$names))) > 0 ) {
    print("Unmatched names :")
    print(basal.activity[is.na(match(as.character(basal.activity),model.structure$names))])
    stop("problem matching basal.activity names");
  }
  if ( !is.null(dim(stimuli)) && !is.null(dim(inhibitor)) &&  dim(stimuli)[1]!=dim(inhibitor)[1] ) {
    stop("number of experiments does not match in stimuli / inhibitor");
  } else if ( !is.null(dim(stimuli)) && dim(stimuli)[2]!=length(c(stim.nodes))) {
    stop("number stimuli does not match in stimuli / stim.nodes");
  } else if ( !is.null(dim(inhibitor)) && dim(inhibitor)[2]!=length(c(inhib.nodes)) ) {
    stop("number inhibitors does not match in inhib.nodes / inhibitor");
  }

  ExperimentalDesign <- fitmodel:::ExperimentalDesign;
  expdes=new(ExperimentalDesign);

  expdes$stim_nodes=match(as.character(stim.nodes),model.structure$names)-1;
  expdes$inhib_nodes=match(as.character(inhib.nodes),model.structure$names)-1;
  expdes$measured_nodes=match(as.character(measured.nodes),model.structure$names)-1;
  expdes$basal_activity=model.structure$names %in% as.character(basal.activity);
  
  expdes$set_inhibitor(inhibitor);
  expdes$set_stimuli(stimuli);

  return(expdes)

}

