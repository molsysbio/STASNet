#!/usr/bin/Rscript
# -*- coding=utf-8 -*-

library("methods")
#library("fitmodel")

source("R/generate_model.R");
source("R/create_model.R");

# Create a model from the data and fit a minimal model
# Takes relative paths as arguments in the order network data basal

data = ""
network = ""
basal_nodes = ""
variation = ""

# Collect the filenames based on their extension
for (argument in commandArgs()) {
    if (grepl(".tab$", argument)) {
        network = paste0(getwd(), "/", argument)
    } else if (grepl(".data$", argument)) {
        data = paste0(getwd(), "/", argument)
        data_name = basename(argument)
    } else if (grepl(".dat$", argument)) {
        basal_nodes = paste0(getwd(), "/", argument)
    } else if (grepl(".var$", argument)) {
        variation = paste0(getwd(), "/", argument)
    }
}

if (network == "") {
    stop("A network structure (adjacency list, .tab) file is required.")
} else if (data == "") {
    stop("A data file (.data) is required")
} else if (basal_nodes == "") {
    stop("A basal activity (.dat) list file is required")
}

#### Creates the model from network and basal files and fits a minimal model to the data
model = create_model(network, data, basal_nodes, variation);
plot_model_accuracy(model);

profiles = profile_likelihood(model, 10);

ni_pf_plot(profiles, data_name=data_name);

# IDEAS :
# data as last argument, possibility to give severals, in which case a comparison of the models is also done
# one argument : name of a file with network, basal_nodes, and the name of the data files

