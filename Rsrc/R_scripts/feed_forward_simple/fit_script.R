#!/usr/bin/Rscript
# -*- coding=utf-8 -*-

library("methods")

source("./generate_model.R");
source("./create_model.R");
source("./randomLHS.r");

# Create a model from the data and fit a minimal model
# Takes relative paths as arguments in the order network data basal
 
data_name = "sim_data.data"
data = data_name
network = "small_network.tab"
basal_nodes = "basal.dat"

if (length(commandArgs()) >= 8) {
    data_name = commandArgs()[7]
    data = paste(getwd(), "/", data_name, sep="")
    network = paste(getwd(), "/", commandArgs()[6], sep="")
    basal_nodes = paste(getwd(), "/", commandArgs()[8], sep="")
    print(data)
}

#### Creates the model from network and basal files and fits a minimal model to the data
model = create_model(network, data, basal_nodes);

sorted_profiles = profile_likelihood(model);

ni_pf_plot(sorted_profiles, data_name=data_name);

# IDEAS :
# data as last argument, possibility to give severals, in which case a comparison of the models is also done
# one argument : name of a file with network, basal_nodes, and the name of the data files

