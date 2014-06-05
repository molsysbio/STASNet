#!/usr/bin/Rscript
# -*- coding=utf-8 -*-

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

data = ""
network = ""
basal_nodes = ""
variation = ""
nb_steps = 10000
inits = 1000
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

#### Creates the model from network and basal files and fits a minimal model to the data
init_time = proc.time()["elapsed"];
model = createModel(network, data, basal_nodes, variation, inits=inits, cores=cores);
get_running_time(init_time, paste("to build the model with", inits, "initialisations."))

power = c("", "k", "M", "G", "T", "P", "Y");
power_init = floor(log(inits, base=1000))
conditions = paste( gsub("(_MIDAS)?.(csv|data)", "", data_name), "_", gsub(".tab", "", network), "_", inits%/%(1000^power_init), power[1+power_init] , sep="" );

pdf(paste0("accuracy_heatmap_", data_name, ".pdf"))
plotModelAccuracy(model, conditions);
dev.off()
print_parameters(model)

# Perform the profile likelihood
profiles = profileLikelihood(model, nb_steps);
model = addPLinfos(model, profiles);

exportModel(model, paste0(conditions, ".mra"));

niPlotPL(profiles, data_name=data_name

get_running_time(init_time, paste("to run the program with", nb_steps, "points for the profile likelihood."));


# IDEAS :
# data as last argument, possibility to give severals, in which case a comparison of the models is also done
# one argument : name of a file with network, basal_nodes, and the name of the data files

# Display the warnings if there are some
warnings()

