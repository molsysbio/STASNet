#!/usr/bin/Rscript
#-*-coding:utf-8-*-

source("./randomLHS.r");
source("./create_model.R");
source("./generate_model.R");

# This is a local script made for debugging without altering the ~/bin code
# Works with a local create_model.R

#### Creates the model from network and basal files and fits a minimal model to the data
model = create_model("example_data/network.tab", "example_data/msb201329_ht29.tab", "example_data/basal.dat");

minimal_fit(model);

