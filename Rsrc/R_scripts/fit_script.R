#source("~/network_reverse_engineering-1.0/r_binding/fitmodel/R/generate_model.R");
source("./create_model.R");

# This is a local script made for debugging without altering the ~/bin code
# Works with a local create_model.R

#### Creates the model from network and basal files and fits a minimal model to the data

model = create_model("network.tab", "msb201329_ht29.tab", "basal.dat");

minimal_fit(model);

