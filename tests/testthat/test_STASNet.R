#context("General testing of STASNet")
DATA_FILE = "test_model_no_error_midas.csv"
VAR_FILE = ""

context("Model fitting accuracy")

model = createModel("network.tab", "basal.dat", DATA_FILE, VAR_FILE, inits=1000, nb_cores=2, perform_plots=F, method="geneticlhs")
test_that("The model fit is completed", {
# Check that we get the fit we expect
    expect_equal_to_reference(model$bestfit, "bestfit.rds")
    expect_equal_to_reference(model$Rscores, "rscores.rds", tolerance=1e-5)
})

test_that("The model design is loaded correctly", {
    expect_equal_to_reference(model$design$inhibitor, "design_inhibitor.rds")
    expect_equal_to_reference(model$design$inhib_nodes, "design_inhib_nodes.rds")
    expect_equal_to_reference(model$design$stimuli, "design_stimuli.rds")
    expect_equal_to_reference(model$design$stim_nodes, "design_stim_nodes.rds")
    expect_equal_to_reference(model$design$measured_nodes, "design_measured_nodes.rds")
})

test_that("The model structure is loaded correctly", {
    expect_equal_to_reference(model$structure$names, "structure_names.rds")
    expect_equal_to_reference(model$structure$adjacencyMatrix, "structure_adjacencyMatrix.rds")
})

test_that("The data are loaded correctly", {
    expect_equal_to_reference(model$data$stim_data, "data_stim_data.rds")
    expect_equal_to_reference(model$data$unstim_data, "data_unstim_data.rds")
    expect_equal_to_reference(model$data$error, "data_error.rds")
})

test_that("The computation is consistent", {
    expect_equal( sum( ((model$model$simulate(model$data, model$parameters)$prediction - model$data$stim_data) / model$data$error)^2, na.rm=T ), model$bestfit )
})

context("Model manipulations work")

test_that("Import-export works correctly", {
    expect_output(exportModel(model, "model.mra"), NA)
    expect_output(importModel("model.mra"), NA)
    expect_output(rebuildModel("model.mra", DATA_FILE), NA)
})

test_that("getCombinationMatrix works properly", {
    expect_equal(getCombinationMatrix(c("N1", "N2i", "N3i"), 1, 1), matrix(c(0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1), ncol=3, dimnames=list(NULL, c("N1", "N2i", "N3i")))) # Example case
    expect_equal(getCombinationMatrix(c("N1", "N2", "N2i", "N3i"), 1, 1), matrix(c(0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1), ncol=4, dimnames=list(NULL, c("N1", "N2", "N2i", "N3i")))) 
    expect_equal(getCombinationMatrix(c("N2i", "N3i"), 1, 1), matrix(c(1, 0, 0, 1), ncol=2, dimnames=list(NULL, c("N2i", "N3i")))) # With stimulation missing
    expect_equal(getCombinationMatrix(c("N1", "N2"), 1, 1), matrix(c(1, 0, 0, 1), ncol=2, dimnames=list(NULL, c("N1", "N2")))) # With inhibition missing
})

test_that("Simulations works properly", {
    res = createSimulation(as.matrix(read.csv("header.csv")), getCombinationMatrix(c("N1", "N2i", "N3i"), 1, 1), paste0("N", 2:4), noise=0.2, replicates=3 )
    expect_equal( res$noise_free_simulation, matrix(c(0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 10, 10, 10, 27.18281828, 27.18281828, 27.18281828, 10, 1.353352832, 10, 73.89056098, 2.82453563, 73.89056098, 10, 0.1831563888, 1.353352832, 545.9815003, 0.7978001573, 5.89499011816367), ncol=6, dimnames=list(NULL, c("TR:N1", "TR:N2i", "TR:N3i", "DV:N2", "DV:N3", "DV:N4"))), tolerance=1e-5)
    # Inhibitions only
    #res = createSimulation(as.matrix(read.csv("header.csv")), getCombinationMatrix(c("N2i", "N3i"), 1, 1), paste0("N", 2:4), noise=0.2, replicates=3 )
    #expect_output(createSimulation("header.csv", "", ""))
})
