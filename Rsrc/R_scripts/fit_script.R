#!/usr/bin/Rscript
# -*- coding=utf-8 -*-

source("./generate_model.R");
source("./create_model.R");
source("./randomLHS.r");

# Create a model from the data and fit a minimal model
# Takes relative paths as arguments in the order network data basal

data = "example_data/HT29_raw_data.data" 
network = "example_data/network.tab"
basal_nodes = "example_data/basal.dat"

if (length(commandArgs(TRUE)) >= 4) {
    data = paste(getwd(), commandArgs(TRUE)[3])
    network = paste(getwd(), commandArgs(TRUE)[2])
    basal_nodes = paste(getwd(), commandArgs(TRUE)[4])
}

#### Creates the model from network and basal files and fits a minimal model to the data
model = create_model(network, data, basal_nodes);

sorted_profiles = draw_profiles(model);

ni_pf_plot(sorted_profiles);

# IDEAS :
# data as last argument, possibility to give severals, in which case a comparison of the models is also done
# one argument : name of a file with network, basal_nodes, and the name of the data files

