library(fitmodel)
context("General testing of fitmodel")
DATA_FILE = "test_model_no_error_midas.csv"
VAR_FILE = ""

model = createModel("network.tab", "basal.dat", DATA_FILE, VAR_FILE, inits=1000, nb_cores=1, perform_plots=F, method="geneticlhs")
# Check that we get the fit we expect
expect_equal_to_reference(model$bestfit, "bestfit.rds")
expect_equal_to_reference(model$Rscores, "rscores.rds", tolerance=1e-5)

test_that("The model design is loaded correctly", {
    expect_equal_to_reference(model$design$inhibitor, "design_inhibitor.rds")
    expect_equal_to_reference(model$design$inhib_nodes, "design_inhib_nodes.rds")
    expect_equal_to_reference(model$design$stimuli, "design_stimuli.rds")
    expect_equal_to_reference(model$design$stim_nodes, "design_stim_nodes.rds")
    expect_equal_to_reference(model$design$measured_nodes, "design_measured_nodes.rds")
})


test_that("The model structure is loaded corretly", {
    expect_equal_to_reference(model$structure$names, "structure_names.rds")
    expect_equal_to_reference(model$structure$adjacencyMatrix, "structure_adjacencyMatrix.rds")
})

test_that("The data are loaded correctly", {
    expect_equal_to_reference(model$data$stim_data, "data_stim_data.rds")
    expect_equal_to_reference(model$data$unstim_data, "data_unstim_data.rds")
    expect_equal_to_reference(model$data$error, "data_error.rds")
})

test_that("Import-export works correctly", {
    expect_output(exportModel(model, "model.mra"), NA)
    expect_output(importModel("model.mra"), NA)
    expect_output(rebuildModel("model.mra", DATA_FILE), "Reading data from")
})
