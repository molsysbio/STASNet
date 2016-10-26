#!/usr/bin/Rscript
# -*- coding:utf-8 -*-

# Hidden from the R installer but with the other scripts from the package
#### LIBRARIES ####
library("STASNet")

#### HELPER FUNCTIONS ####
# Print the time it took in a readable format
get_running_time <- function(init_time, text="") {
    run_time = proc.time()["elapsed"]-init_time
    run_hours = run_time %/% 3600;
    run_minutes = (run_time - 3600 * run_hours) %/% 60;
    run_seconds = run_time - 3600 * run_hours - 60 * run_minutes;
    print(paste(run_hours, "h", run_minutes, "min", run_seconds, "s", text))
}

#### SETUP ####
# Create a model from the data and fit a minimal model
# Takes relative paths as arguments in the order network data basal

data_files = ""
network = ""
basal_nodes = ""
var_files = ""
inits = 10000
var_samples=10
method = "geneticlhs"
relax = TRUE
perform_pl = FALSE
# Autodetection of the cores
cores = 0

#### READ IN ARGUMENTS ####
if (!exists("cargs")) {
  cargs = commandArgs(trailingOnly=T)
} else if (is.character(cargs)) {
  cargs = strsplit(cargs, " ")[[1]]
}

# Collect the filenames based on their extension
if (any(grepl(".tab$", cargs))) {
  network = paste0(getwd(), "/", cargs[grepl(".tab$", cargs)])
}
if (any(grepl(".data$|.csv$", cargs))) {
  data_files = paste0(getwd(), "/", cargs[grepl(".data$|.csv$", cargs)])
  data_files_name = basename(cargs[grepl(".data$|.csv$", cargs)])
  if (length(data_files)<2){
    stop("Insufficient data files to perform a model set fitting, at least 2 are required!!")
  }
  data_files=data_files[order(data_files)]
  data_files_name=data_files_name[order(data_files)]
}
if (any(grepl(".dat$", cargs))) {
  basal_nodes = paste0(getwd(), "/", cargs[grepl(".dat$", cargs)])
} 
if (any(grepl(".var$", cargs))) {
  var_files = paste0(getwd(), "/", cargs[grepl(".var$", cargs)])
  if (length(var_files)!=length(data_files)){
    stop("Please supply as many var files as data files")
  }
  var_files=var_files[order(var_files)]
}
if (any(grepl("^-i", cargs))) {
  inits = as.numeric(gsub("-i", "", cargs[grepl("^-i", cargs)]))
  if (is.na(inits)) { # If error
    inits = 1000
    print("Incorrect number of initialisations, using 1000 instead")
  }
}
if (any(grepl("^-c", cargs))) {
  cores = as.numeric(gsub("-c", "", cargs[grepl("^-c", cargs)]))
  if (is.na(cores)) { # If error
    cores = 1
    print("Incorrect number of cores (use 1 instead)")
  }
} 
if (any(grepl("^-s", cargs))) {
  var_samples = as.numeric(gsub("^-s", "", cargs[grepl("^-s", cargs)]))
  if (is.na(var_samples)) {
    var_samples = 10
    print("Incorrect number of samples, performing with 10")
  }
} 
if  (any(grepl("^-m", cargs)))
  method = gsub("^-m", "", cargs[grepl("^-s", cargs)])

if (any(grepl("^-nr", cargs)))
  relax = FALSE

if (any(cargs == "--nopl"))
  perform_pl = FALSE

if (cores == 0) {
  cores = detectCores() - 1;
}

# sanity checks
if (network == "") {
  stop("A network structure (adjacency list, .tab) file is required.")
}
if (data_files[1] == ""){ 
  stop("data files (.data) arguments required")
}
if (basal_nodes == ""){
  stop("A basal activity (.dat) list file is required")
}

# Extract the name and the number of initialisations
power = c("", "k", "M", "G", "T", "P", "Y");
power_init = floor(log(inits, base=1000))
conditions = paste0( paste0(gsub("(_MIDAS)?.(csv|data)", "", basename(data_files_name)),collapse="_"), "_", gsub(".tab", "", basename(network)), "_", inits%/%(1000^power_init), power[1+power_init]);
conditions = gsub(" ", "_", conditions)
folder = paste0( "run_", conditions, "_", Sys.Date(), "/" )
dir.create(folder)

#### 1 FIT MODEL SET ####
init_time = proc.time()["elapsed"];
#pdf(paste0(folder, "distribution_", conditions, ".pdf"))
modelset = createModelSet(model_links = network,
                          basal_file = basal_nodes,
                          csv_files = data_files,
                          var_files = var_files,
                          nb_cores = cores,
                          inits = inits,
                          perform_plots = perform_pl,
                          method = method)
#dev.off()
get_running_time(init_time, paste("to build the model with", inits, "initialisations."))
writeLines(modelset$infos[2])

#### 2 FIND VARIABLE PARAMETERS ####
if (relax){
modelset=addVariableParameters(modelset = modelset,
                               nb_cores = cores,
                              max_iterations=0,
                              nb_samples=var_samples,
                              accuracy=0.95)
  }

#### 3 PLOT RESULTS ####
modelgroup=extractSubmodels(modelset)

for ( ii in 1:length(modelgroup$names)){
model=modelgroup$models[[ii]]
# 3.1 plot graph
pdf(paste0(folder,"graph_",model$name,".pdf"))
plotModelGraph(model)
dev.off()

# 3.2 plot acuracy heatmap
mat=model$data$stim_data
pdf(paste0(folder, "accuracy_heatmap_", model$name, ".pdf"),onefile=T,width =5+ncol(mat)/3,height=4+nrow(mat)/6)
plotModelAccuracy(model)
plotModelScores(model, main=paste0("Global R = ", model$bestfitscore))
dev.off()

printParameters(model)

# 3.3 Perform the profile likelihood
#if (perform_pl) {
#    profiles = profileLikelihood(model, nb_steps, nb_cores=min(cores, length(model$parameters)));
#    model = addPLinfos(model, profiles);
#    get_running_time(init_time, paste("to run the program with", nb_steps, "points for the profile likelihood."));
#    niplotPL(profiles, data_name=conditions, folder=folder)
#}

# 3.4 Plot the simulation for all combinations of inhibitors 
#pdf(paste0("combos_", conditions, ".pdf"))
#plotModelSimulation(model, getCombinationMatrix(c("MEKi", "GSK3ABi", "IGF", "TGFA", "PI3Ki")))
#dev.off()
# Plot the simulated conditions
pdf(paste0(folder, "model_simulation_", model$name, ".pdf"))
plotModelSimulation( model )
dev.off()

#if (reduction) {
## Reduce the model and see what changed
#    print("Reduction of the model...")
#    reduced_model = selectMinimalModel(model)
## Profile likelihood on the reduced model
#    reduced_profiles = profileLikelihood(reduced_model, nb_steps, nb_cores=min(cores, length(reduced_model$parameters)));
#    reduced_model = addPLinfos(reduced_model, reduced_profiles)
#    exportModel(reduced_model, paste0(folder, "reduced_", conditions, ".mra"));
#    niplotPL(reduced_profiles, data_name=paste0("reduced_", data_name))
## Plot the simulated conditions
#    pdf(paste0(folder, "reduced_all_", conditions, ".pdf"))
#    plotModelSimulation( simulateModel(reduced_model) )
#    dev.off()
#
#    get_running_time(init_time, "with the model reduction");
#}
}

#### EXPORT MODEL DATA ####
#  export Model data
for (ii in 1:length(modelgroup$names))
exportModel(modelgroup$models[[ii]], paste0(folder,"../",modelgroup$models[[ii]]$name,".mra"))

print("Finished")

# Display the warnings if there are some
warnings()
