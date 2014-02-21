<<<<<<< HEAD
#!/usr/bin/Rscript
#-*-coding:utf-8-*-

source("./randomLHS.r");
source("./create_model.R");
source("./generate_model.R");
=======
#source("~/network_reverse_engineering-1.0/r_binding/fitmodel/R/generate_model.R");
source("./create_model.R");
>>>>>>> 1124d250dad246e89a90e9592a757f48ab8c0d96

# This is a local script made for debugging without altering the ~/bin code
# Works with a local create_model.R

#### Creates the model from network and basal files and fits a minimal model to the data
<<<<<<< HEAD
model = create_model("example_data/network.tab", "example_data/msb201329_ht29.tab", "example_data/basal.dat");
=======

model = create_model("network.tab", "msb201329_ht29.tab", "basal.dat");
>>>>>>> 1124d250dad246e89a90e9592a757f48ab8c0d96

minimal_fit(model);

