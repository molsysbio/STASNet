#!/usr/bin/Rscript
# -*- coding:utf-8 -*-

# Hidden from the R installer but with the other scripts from the package

library("methods")
library("fitmodel")
library("parallel")

# Print the time it took in a readable format
get_running_time <- function(init_time, text="") {
    run_time = proc.time()["elapsed"]-init_time
    run_hours = run_time %/% 3600;
    run_minutes = (run_time - 3600 * run_hours) %/% 60;
    run_seconds = run_time - 3600 * run_hours - 60 * run_minutes;
    print(paste(run_hours, "h", run_minutes, "min", run_seconds, "s", text))
}

#source("generate_model.R");
#source("create_model.R");

# Create a model from the data and fit a minimal model
# Takes relative paths as arguments in the order network data basal

reduction = F

data = ""
network = ""
basal_nodes = ""
variation = ""
nb_steps = 1000
inits = 1000
method = "default"
# Autodetection of the cores
cores = 0

# Collect the filenames based on their extension
for (argument in commandArgs()) {
    if (grepl(".tab$", argument)) {
        network = paste0(getwd(), "/", argument)
    } else if (grepl(".data$", argument) || grepl(".csv$", argument)) {
        data = paste0(getwd(), "/", argument)
        data_name = basename(argument)
    } else if (grepl(".dat$", argument)) {
        basal_nodes = paste0(getwd(), "/", argument)
    } else if (grepl(".var$", argument)) {
        variation = paste0(getwd(), "/", argument)
    } else if (grepl("^-i", argument)) {
        inits = as.numeric(gsub("-i", "", argument))
        if (is.na(inits)) { # If error
            inits = 1000
            print("Incorrect number of initialisation, using 1000 instead")
        }
    } else if (grepl("^-c", argument)) {
        cores = as.numeric(gsub("-c", "", argument))
        if (is.na(cores)) { # If error
            cores = 2
            stop("Incorrect number of cores (use -c#)")
        }
    } else if (grepl("^-s", argument)) {
        nb_steps = as.numeric(gsub("^-s", "", argument))
        if (is.na(nb_steps)) {
            nb_steps = 1000
            print("Incorrect number of steps, performing with 1000")
        }
    } else if (grepl("^-nr$", argument)) {
        reduction = FALSE
    } else if (grepl("^-m", argument)) {
        method = gsub("^-m", "", argument)
    }
}
if (cores == 0) {
    cores = detectCores() - 1;
}
if (network == "") {
    stop("A network structure (adjacency list, .tab) file is required.")
} else if (data == "") {
    stop("A data file (.data) is required")
} else if (basal_nodes == "") {
    stop("A basal activity (.dat) list file is required")
}

# Extract the name and the number of initialisations
power = c("", "k", "M", "G", "T", "P", "Y");
power_init = floor(log(inits, base=1000))
conditions = paste0( gsub("(_MIDAS)?.(csv|data)", "", basename(data_name)), "_", gsub(".tab", "", basename(network)), "_", inits%/%(1000^power_init), power[1+power_init]);
conditions = gsub(" ", "_", conditions)

#### Creates the model from network and basal files and fits a minimal model to the data
init_time = proc.time()["elapsed"];
pdf(paste0("distribution_", conditions, ".pdf"))
model = createModel(network, data, basal_nodes, variation, inits=inits, cores=cores, init_distribution=T, method=method);
dev.off()
get_running_time(init_time, paste("to build the model with", inits, "initialisations."))

pdf(paste0("accuracy_heatmap_", conditions, ".pdf"))
accuracyPlot(model);
dev.off()
printParameters(model)

if (method == "annealing") {
    stop("debugging annealing")
}

# Perform the profile likelihood
profiles = profileLikelihood(model, nb_steps, cores=min(cores, length(model$parameters)));
model = addPLinfos(model, profiles);
get_running_time(init_time, paste("to run the program with", nb_steps, "points for the profile likelihood."));
exportModel(model, paste0(conditions, ".mra"));
niplotPL(profiles, data_name=data_name)

# Plot the simulation for all combinations of inhibitors 
#pdf(paste0("combos_", conditions, ".pdf"))
#plotModelPrediction(model, getCombinationMatrix(c("MEKi", "GSK3ABi", "IGF", "TGFA", "PI3Ki")))
#dev.off()
# Plot the simulated conditions
pdf(paste0("all_", conditions, ".pdf"))
plotModelPrediction(model, "all")
dev.off()

if (reduction) {
# Reduce the model and see what changed
    print("Reduction of the model...")
    reduced_model = selectMinimalModel(model)
# Profile likelihood on the reduced model
    reduced_profiles = profileLikelihood(reduced_model, nb_steps, cores=min(cores, length(reduced_model$parameters)));
    reduced_model = addPLinfos(reduced_model, reduced_profiles)
    exportModel(reduced_model, paste0("reduced_", conditions, ".mra"));
    niplotPL(reduced_profiles, data_name=paste0("reduced_", data_name))
# Plot the simulated conditions
    pdf(paste0("reduced_all_", conditions, ".pdf"))
    plotModelPrediction(reduced_model, "all")
    dev.off()

    get_running_time(init_time, "with the model reduction");
}

print("Finished")



# IDEAS :
# data as last argument, possibility to give severals, in which case a comparison of the models is also done
# one argument : name of a file with network, basal_nodes, and the name of the data files

# Display the warnings if there are some
warnings()

