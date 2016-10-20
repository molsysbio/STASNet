#context("General testing of STASNet")
DATA_FILE = "test_model_no_error_midas.csv"
VAR_FILE = ""

context("Model fitting accuracy")

model = suppressMessages( createModel("network.tab", "basal.dat", DATA_FILE, VAR_FILE, inits=1000, nb_cores=2, perform_plots=F, method="geneticlhs") )

test_that("All expected fields are present" ,{
    expect_equal(exists("min_cv", model), TRUE)
    expect_equal(exists("default_cv", model), TRUE)
})

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

context("Score computation functions")

test_that("computeFitScore works as expected", {
    computeFitScore(model)
    computeFitScore(model, T)
})

context("Model import-export")

test_that("Import-export works correctly", {
    expect_output(exportModel(model, "model.mra"), NA)
    expect_output(importModel("model.mra"), NA)
    expect_output(rebuildModel("model.mra", DATA_FILE), NA)
})

context("Simulation helper functions")

test_that("getCombinationMatrix works properly", {
    expect_equal(getCombinationMatrix(c("N1", "N2i", "N3i"), 1, 1), matrix(c(0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1), ncol=3, dimnames=list(NULL, c("N1", "N2i", "N3i")))) # Example case
    expect_equal(getCombinationMatrix(c("N1", "N2", "N2i", "N3i"), 1, 1), matrix(c(0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1), ncol=4, dimnames=list(NULL, c("N1", "N2", "N2i", "N3i")))) 
})
# With stimulation missing
test_that("getCombinationMatrix without stimulation", {
    expect_equal(getCombinationMatrix(c("N2i", "N3i"), 1, 1), matrix(c(0, 1, 0, 0, 0, 1), ncol=2, dimnames=list(NULL, c("N2i", "N3i"))))
})
# With inhibition missing
test_that("getCombinationMatrix without inhibition", {
    expect_equal(getCombinationMatrix(c("N1", "N2"), 1, 1), matrix(c(0, 1, 0, 0, 0, 1), ncol=2, dimnames=list(NULL, c("N1", "N2"))))
})
test_that("getCombinationMatrix raises an error if stim_combo is a string", {
    expect_error( getCombinationMatrix(c("N1", "N3i"), 1, "N2i"), "stim_combo" )
})
test_that("getCombinationMatrix raises an error if stim_combo is a string", {
    expect_error( getCombinationMatrix(c("N1", "N3i"), "N2i"), "inhib_combo" )
})

context("Toy data and reverse fitting")

test_that("Toy data can be generated" ,{
    expect_silent( createSimulation(as.matrix(read.csv("header.csv")), getCombinationMatrix(c("N1", "N2i", "N3i"), 1, 1), paste0("N", 2:4), noise=0.2, replicates=3 ) )
    .GlobalEnv$res = createSimulation(as.matrix(read.csv("header.csv")), getCombinationMatrix(c("N1", "N2i", "N3i"), 1, 1), paste0("N", 2:4), noise=0.2, replicates=3 )
    .GlobalEnv$simulated_data = res$noise_free_simulation[,4:6]
    colnames(.GlobalEnv$simulated_data) = gsub("^[A-Z]+:", "", colnames(.GlobalEnv$simulated_data))
})

test_that("createSimulation works properly", {
    expect_equal( res$noise_free_simulation, matrix(c(0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 10, 10, 10, 27.18281828, 27.18281828, 27.18281828, 10, 1.353352832, 10, 73.89056098, 2.82453563, 73.89056098, 10, 0.1831563888, 1.353352832, 545.9815003, 0.7978001573, 5.89499011816367), ncol=6, dimnames=list(NULL, c("TR:N1", "TR:N2i", "TR:N3i", "DV:N2", "DV:N3", "DV:N4"))), tolerance=1e-5)
    # Inhibitions only
    #expect_silent( createSimulation(as.matrix(read.csv("header.csv")), getCombinationMatrix(c("N2i", "N3i"), 1, 1), paste0("N", 2:4), noise=0.2, replicates=3 ) )
    #expect_output(createSimulation("header.csv", "", ""))
    expect_silent( createSimulation(as.matrix(read.csv("header.csv")), getCombinationMatrix(c("N1", "N2i", "N3i"), 1, 1), paste0("N", 2:4), noise=0, replicates=1 ) )
})

#res = createSimulation(as.matrix(read.csv("header.csv")), getCombinationMatrix(c("N1", "N2i", "N3i"), 1, 1), paste0("N", 2:4), noise=0, replicates=1 )
test_that("Noise free toy data is properly refited", {
    midas=data.frame("ID:type"=c("c", rep("t", 5)), res$noise_free_simulation, "DA:ALL"=rep(0, 6))
    .GlobalEnv$refit = suppressWarnings(createModel(res$model$structure, res$model$basal, midas, inits=100))
    expect_equal(refit$parameters, c(1.0,2.0,2.0,-1,-1))
    expect_equal(refit$bestfit, 0)
})


context("Model reduction")

test_that("Model gets reduced correctly", {
    expect_output(reduceModel(model), "")
    red_model = reduceModel(model)
})

context("Score computations functions")

# Dumb model
incomplete_model = list()
incomplete_model$data = list()
incomplete_model$data$stim_data = matrix(0, ncol=4, nrow=4)
incomplete_model$parameters = 1:6
incomplete_model$bestfit = 50
test_that("The reduced chi score is computed correctly", {
    expect_equal(computeReducedChiScore(incomplete_model)$reducedChi, 5)
})

context("Model extension")

test_that("The extension function works withtout bugs", {
    expect_output(suggestExtension(model, T, 8), "")
})


context("Model simulations")

test_that("Simulation of the data", {
    expect_silent( simulateModel(refit) ) # Easy case
    .GlobalEnv$resimulation = simulateModel(refit)
})
test_that("Simulation of a subset of the conditions", {
    expect_silent( simulateModel(refit, getCombinationMatrix(c("N1", "N2i"), 1)) )
    .GlobalEnv$simulation_sconditions = simulateModel(refit, getCombinationMatrix(c("N1", "N2i"), 1))
})
test_that("Duplicated conditions are eliminated", {
    expect_silent( simulateModel(refit, rbind(0, getCombinationMatrix(c("N1", "N2i"), 1))) )
})
test_that("The data are subseted properly", {
    expect_equal( simulation_sconditions$data, simulated_data[c(1, 4, 2, 5), ] )
})
test_that("Simulation without inhibitions and only one stimulation", {
    expect_silent( simulateModel(refit, getCombinationMatrix(c("N1"), 1)) )
})
test_that("Simulation without stimulations", {
    expect_silent( simulateModel(refit, getCombinationMatrix(c("N2i"), 1)) )
})
test_that("Simulation without perturbations", {
    expect_error( simulateModel(refit, matrix(ncol=0, nrow=0)), "No valid perturbations" )
})
test_that("Simulation without wrong perturbations", {
    expect_error( simulateModel(refit, getCombinationMatrix(c("N4", "N1i"), 1)), "No valid perturbations" )
})
test_that("Simulation without valid perturbations", {
    expect_error( simulateModel(refit, getCombinationMatrix(c("N4", "N1i"), 1)), "No valid perturbations" )
})
test_that("Simulation of a subset of the readouts", {
    expect_silent( simulateModel(refit, readouts=c("N3", "N4")) )
    .GlobalEnv$simulation_sreadouts = simulateModel(refit, readouts=c("N3", "N4"))
})
test_that("Simulation of a subset of the conditions and readouts", {
    expect_silent( simulateModel(refit, getCombinationMatrix(c("N1", "N2i"), 1), c("N3", "N4")) )
})
test_that("Message when inexistant readouts is required", {
    expect_message( simulateModel(refit, readouts=c("N3", "N5")), "is not in the network" )
})
test_that("Message when invalid readouts is required", {
    expect_message( simulateModel(refit, readouts=c("N3", "N1")), "cannot be measured" )
})
test_that("Message when inexistant stimulation is required", {
    expect_message( simulateModel(refit, getCombinationMatrix(c("N1", "N2"))), "not stimulated" )
    .GlobalEnv$nostim_sim = simulateModel(refit, getCombinationMatrix(c("N1", "N2")))
})
test_that("The final matrix is correctly reduced with false stimulation", {
    expect_equal( nrow(nostim_sim$bestfit), 2)
}) 
test_that("Message when inexistant inhibition is required", {
    expect_message( simulateModel(refit, getCombinationMatrix(c("N1", "N1i"), 1)), "not inhibited" )
    .GlobalEnv$noinhib_sim = simulateModel(refit, getCombinationMatrix(c("N1", "N1i"), 1))
})
test_that("The final matrix is correctly reduced with false inhibition", {
    expect_equal( nrow(noinhib_sim$bestfit), 2)
})
test_that("Prediction of new conditions", {
    expect_silent( simulateModel(refit, getCombinationMatrix(c("N1", "N2i", "N3i"))) )
})
