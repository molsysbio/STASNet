#context("General testing of STASNet")
#### model tests ####

context("Helper functions")

test_that("extractStructure correctly loads all nodes", {
    t_struct = extractStructure("test_structure.tab")
    expect_equal_to_reference(t_struct$names, "test_structure_names.rds")
    expect_equal_to_reference(t_struct$adjacencyMatrix, "test_structure_adj.rds")
})


context("Model")

DATA_FILE = "test_model_no_error_midas.csv"
VAR_FILE = ""

context("Model fitting accuracy")

test_that("createModel works", {
    expect_output( createModel("network.tab", "basal.dat", DATA_FILE, VAR_FILE, inits=1000, nb_cores=2, perform_plots=F, method="geneticlhs",rearrange = "bystim"), NA)
})
model = suppressMessages( createModel("network.tab", "basal.dat", DATA_FILE, VAR_FILE, inits=1000, nb_cores=2, perform_plots=F, method="geneticlhs",rearrange = "bystim") )
test_that("createModel works with unused readouts and perturbations", {
     expect_output( createModel("network.tab", "basal.dat", DATA_FILE, VAR_FILE, inits=1000, nb_cores=2, perform_plots=F, method="geneticlhs",rearrange = "bystim", unused_readout=c("node2"), unused_perturbation=c("node2i")), NA )
})
sub_model = suppressMessages( createModel("network.tab", "basal.dat", DATA_FILE, VAR_FILE, inits=1000, nb_cores=2, perform_plots=F, method="geneticlhs",rearrange = "bystim", unused_readout=c("node2"), unused_perturbation=c("node2i")) )
test_that("createModel works with log-space fitting", {
    expect_output( createModel("network.tab", "basal.dat", DATA_FILE, VAR_FILE, inits=1000, nb_cores=2, perform_plots=F, method="geneticlhs",rearrange = "bystim", data_space="log"), NA)
})
log_model = suppressMessages( createModel("network.tab", "basal.dat", DATA_FILE, VAR_FILE, inits=1000, nb_cores=2, perform_plots=F, method="geneticlhs",rearrange = "bystim", data_space="log") )
test_that("log-space models are fitted correctly (basal activity==0 and no replicates to compute the error)", {
    expect_equal(log_model$use_log, TRUE)
    expect_equal(is.na(log_model$bestfit), FALSE)
})

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


context("Score computation functions")

test_that("computeFitScore works as expected", {
    expect_silent(computeFitScore(model))
    expect_silent(computeFitScore(model, T))
})
test_that("The fit score computation for MRAmodel is consistent", {
    expect_equal( computeFitScore(model)$bestfit, model$bestfit )
})


context("Profile likelihood")

test_that("profileLikelihood works", {
     expect_message(profileLikelihood(model, 100, 0), "evaluate")
     .GlobalEnv$pl_data = suppressMessages( profileLikelihood(model, 1000, 0) )
})
test_that("Profile likelihood residuals are consistent", { # This should be consistent but it is not...
    expect_equal_to_reference(sapply(pl_data, function(xx) {xx$residuals}), "pl_data_residuals.rds", tolerance=1e-5)
})
test_that("Profile likelihood value are consistent", {
    expect_equal_to_reference(sapply(pl_data, function(xx) {xx$value}), "pl_data_value.rds", tolerance=1e-5)
})
test_that("Profile likelihood thresholds are consistent", {
    expect_equal_to_reference(sapply(pl_data, function(xx) {xx$thresholds}), "pl_data_thresholds.rds", tolerance=1e-5)
})
test_that("Profile likelihood explored are consistent", {
    expect_equal_to_reference(sapply(pl_data, function(xx) {xx$explored}), "pl_data_explored.rds", tolerance=1e-5)
})
test_that("addPLinfos works", {
    expect_silent({model=addPLinfos(model, pl_data)})
})
model = addPLinfos(model, pl_data)
test_that("The profile likelihood infos are correctly added to the model", {
    expect_equal_to_reference(model$param_range, "param_range.rds", tolerance=1e-5)
})


context("Model import-export")

test_that("Export works correctly", {
    expect_output(exportModel(model, "model.mra"), NA)
})
test_that("Import works correctly", {
    expect_output(importModel("model.mra"), NA)
})
imported_model = importModel("model.mra")
test_that("Imported model parameters is the same as the exported one", {
    expect_equal(imported_model$parameters, model$parameters)
})
test_that("Imported model cv is the same as the exported one", {
    expect_equal(imported_model$cv, model$cv)
})
test_that("Imported model param_range is the same as the exported one", {
    expect_equal(unlist(imported_model$param_range), unlist(model$param_range))
})
test_that("Imported model stim_data is the same as the exported one", {
    expect_equal(imported_model$data$stim_data, model$data$stim_data)
})
test_that("Imported model error is the same as the exported one", {
    expect_equal(imported_model$data$error, model$data$error)
})
test_that("Imported model scale is the same as the exported one", {
    expect_equal(imported_model$data$scale, model$data$scale)
})
test_that("Imported model unstim_data is the same as the exported one", {
    expect_equal(imported_model$data$unstim_data, model$data$unstim_data)
})
test_that("Import with data in the .mra file works correctly", {
    expect_output(importModel("model.mra"), NA)
})
data_import_model = importModel("model.mra")
test_that("Model imported with data from the .mra file parameters is the same as the exported one", {
    expect_equal(data_import_model$parameters, model$parameters)
})
test_that("Model imported with data from the .mra file cv is the same as the exported one", {
    expect_equal(data_import_model$cv, model$cv)
})
test_that("Model imported with data from the .mra file param_range is the same as the exported one", {
    expect_equal(unlist(data_import_model$param_range), unlist(model$param_range))
})
test_that("Model imported with data from the .mra file stim_data is the same as the exported one", {
    expect_equal(data_import_model$data$stim_data, model$data$stim_data)
})
test_that("Model imported with data from the .mra file error is the same as the exported one", {
    expect_equal(data_import_model$data$error, model$data$error)
})
test_that("Model imported with data from the .mra file scale is the same as the exported one", {
    expect_equal(data_import_model$data$scale, model$data$scale)
})
test_that("Model imported with data from the .mra file unstim_data is the same as the exported one", {
    expect_equal(data_import_model$data$unstim_data, model$data$unstim_data)
})
test_that("Export works correctly with unused readouts and perturbations", {
    expect_output(exportModel(sub_model, "sub_model.mra"), NA)
})
test_that("Import works correctly with unused readouts and perturbations", {
    expect_output(importModel("sub_model.mra"), NA)
})
test_that("Rebuild works correctly with unused readouts and perturbations", {
    expect_output(rebuildModel("sub_model.mra", DATA_FILE), NA)
})
test_that("Rebuild works correctly", {
    expect_output(rebuildModel("model.mra", DATA_FILE), NA)
    .GlobalEnv$reb_model = rebuildModel("model.mra", DATA_FILE)
})
test_that("The rebuilt model is consistent", {
    expect_equal(reb_model$bestfit, model$bestfit)
})
test_that("Rebuild works correctly with appropriate R objects as input", {
  mra_file = readLines("sub_model.mra")
  expect_output(rebuildModel(mra_file, DATA_FILE), NA)
  data_file = STASNet:::extractMIDAS(DATA_FILE)
  expect_output(rebuildModel(mra_file,data_file),NA)
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

context("Cloning model")

test_that("Model is cloned correctly", {
  expect_silent(STASNet:::cloneModel(model))
  alt_model = STASNet:::cloneModel(model)
  expect_false(capture.output(alt_model$model$.pointer) == capture.output(model$model$.pointer))
  expect_false(capture.output(alt_model$design$.pointer) == capture.output(model$design$.pointer))
  expect_false(capture.output(alt_model$structure$.pointer) == capture.output(model$structure$.pointer))
  expect_false(capture.output(alt_model$data$.pointer) == capture.output(model$data$.pointer))
  expect_equal(alt_model$model$modelRank(),model$model$modelRank())
  expect_equal(alt_model$use_log,model$use_log)
})

alt_model = STASNet:::cloneModel(model)

test_that("Cloned model is independent", {
  tmp_adj = alt_model$structure$adjacencyMatrix
  tmp_adj[4,3] = 0
  alt_model$structure$setAdjacencyMatrix(tmp_adj)
  alt_model$model$setModel(alt_model$design, alt_model$structure, FALSE)
  expect_gt(model$model$modelRank(), alt_model$model$modelRank())
})

context("Model reduction")

test_that("Model gets reduced correctly", {
    expect_message(reduceModel(model))
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

test_that("The extension function works without bugs", {
    expect_message(suggestExtension(model, T, 2))
    expect_message(suggestExtension(model, F))
})
test_that("The extension function correctly handles sample range", {
        expect_message(suggestExtension(model, T, 2, c(1,0,-1)))
        expect_error(suggestExtension(model,T,2,c()))
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
test_that("Requested readouts are duplicated", {
    expect_silent( simulateModel(refit, getCombinationMatrix(c("N1", "N2i"), 1), c("N2", "N3", "N2")) )
})
test_that("Prediction of new conditions", {
    expect_silent( simulateModel(refit, getCombinationMatrix(c("N1", "N2i", "N3i"))) )
})

#### modelset tests ####

context("ModelSet")

DATA_FILES = c("test_model_no_error_midas.csv","test_model_no_error_midas_2.csv")#,"test_model_no_error_midas_3.csv")
VAR_FILES = c()

modelset = suppressMessages(createModelSet("network.tab", "basal.dat", DATA_FILES, VAR_FILES,1,100,F, DEFAULT_CV=0.2))


context("Cloning modelset")

test_that("Modelset is cloned correctly", {
  expect_silent(STASNet:::cloneModel(modelset))
  alt_modelset = STASNet:::cloneModel(modelset)
  expect_false(capture.output(alt_modelset$model$.pointer) == capture.output(modelset$model$.pointer))
  expect_false(capture.output(alt_modelset$design$.pointer) == capture.output(modelset$design$.pointer))
  expect_false(capture.output(alt_modelset$structure$.pointer) == capture.output(modelset$structure$.pointer))
  expect_false(capture.output(alt_modelset$data$.pointer) == capture.output(modelset$data$.pointer))
  expect_equal(alt_modelset$model$modelRank(),modelset$model$modelRank())
})

alt_modelset = cloneModel(modelset)

test_that("Cloned modelset is independent", {
  tmp_adj = alt_modelset$structure$adjacencyMatrix
  tmp_adj[4,3] = 0
  alt_modelset$structure$setAdjacencyMatrix(tmp_adj)
  alt_modelset$model$setModel(alt_modelset$design, alt_modelset$structure, FALSE)
  alt_modelset$model$setNbModels(alt_modelset$nb_models)
  expect_gt(modelset$model$modelRank(), alt_modelset$model$modelRank())
})

context("ModelSet fitting accuracy")

test_that("Additional modelSet default fields are present" ,{
  expect_equal(exists("names", modelset), TRUE)
  expect_equal(exists("nb_models", modelset), TRUE)
})

test_that("The modelSet fit is reasonable", {
  expect_equal_to_reference(modelset$bestfit, "ms_bestfit.rds", tolerance=1e-5)
})

test_that("The modelSet information is loaded correctly", {
  expect_equal(modelset$nb_models, 2)
})

test_that("The modelSet structure is loaded correctly", {
  expect_equal_to_reference(modelset$structure$names, "ms_structure_names.rds")
  expect_equal_to_reference(modelset$structure$adjacencyMatrix, "ms_structure_adjacencyMatrix.rds")
})

test_that("The data are loaded correctly", {
  expect_equal_to_reference(modelset$data$stim_data, "ms_data_stim_data.rds")
  expect_equal_to_reference(modelset$data$unstim_data, "ms_data_unstim_data.rds")
  expect_equal_to_reference(modelset$data$error, "ms_data_error.rds")
})

test_that("The fit score computation for ModelSet is consistent", {
  expect_equal( computeFitScore(modelset)$bestfit, modelset$bestfit )
})

test_that("ModelSet breakup works", {
  expect_silent(extractSubmodels(modelset))
})

context("ModelSet relaxation")
test_that("parameters can be relaxed",{
        expect_message(addVariableParameters(modelset, 1, 0, 10))
        relax_modelset = suppressMessages(addVariableParameters(modelset, 1, 0, 10))
        expect_equal(relax_modelset$variable_parameters, 5) # The sqrt(2) factor is enough to make this variable link not detectable (p=0.72)
        expect_equal(computeFitScore(relax_modelset)$bestfit, relax_modelset$bestfit)
})

relax_modelset = suppressMessages(addVariableParameters(modelset, 1, 0, 10))

test_that("variable parameters are kept when cloned", {
  expect_silent(cloneModel(relax_modelset))
  tmp_modelset = cloneModel(relax_modelset)
  expect_equal(tmp_modelset$variable_parameters,relax_modelset$variable_parameters)
  expect_equal(tmp_modelset$parameters,relax_modelset$parameters)
})

test_that("exclusion of parameters from relaxation works",{
  expect_error(addVariableParameters(modelset,1,0,10,0.95,"geneticlhs",notVariable = c("notThere")))
  expect_error(addVariableParameters(modelset,1,0,10,0.95,"geneticlhs",notVariable = c(1,17)))
  no5set = addVariableParameters(modelset,1,0,10,0.95,"geneticlhs",notVariable = 5)
  expect_false(no5set$variable_parameters==5)
  no3_5set = addVariableParameters(modelset,1,0,10,0.95,"geneticlhs",notVariable = c(3,5))
  expect_equal(length(no3_5set$variable_parameters),0)
  expect_warning(addVariableParameters(modelset,1,0,10,0.95,"geneticlhs",notVariable =c(1:5)))
})


context("ModelSet rebuild")
MODEL_FILES=c("var_model1.mra","var_model2.mra")

test_that("Rebuild of variable modelset works correctly", {
  expect_error(rebuildModelSet(MODEL_FILES,DATA_FILES[1]))
  expect_error(rebuildModelSet(MODEL_FILES))
  expect_error(rebuildModelSet(MODEL_FILES[1],DATA_FILES))
  expect_error(rebuildModelSet(MODEL_FILES,DATA_FILES,DATA_FILES[1]))
  expect_output(rebuildModelSet(MODEL_FILES, DATA_FILES), NA)
  reb_modelset <- rebuildModelSet(MODEL_FILES, DATA_FILES)
  expect_equal(reb_modelset$variable_parameters,relax_modelset$variable_parameters)
})

test_that("Rebuild of modelset works with appropriate R objects as inputs", {
  celline = c("cell1","cell2")
  data_list = sapply(celline, function(x) STASNet::extractMIDAS(DATA_FILES[celline == x]),simplify = F)
  mra_list = sapply(celline, function(x) readLines(MODEL_FILES[celline == x]), simplify =F)
  expect_output(rebuildModelSet(MODEL_FILES, data_list),NA)
  expect_output(rebuildModelSet(mra_list, DATA_FILES),NA)
  expect_output(rebuildModelSet(mra_list, data_list),NA)
})

context("ModelSet extension")

test_that("modelSet with fixed parameters is extended correctly", {
    expect_message(suggestExtension(modelset,T))
    exprmat = suppressMessages(suggestExtension(modelset,T))
    expect_equal(all(as.numeric(as.character(exprmat$Res_delta))>=10^-5),T)
})

test_that("modelSet with variable parameters is extended correctly", {
  expect_message(suggestExtension(relax_modelset,T))
  exprmat = suppressMessages(suggestExtension(relax_modelset,T))
  expect_equal(all(as.numeric(as.character(exprmat$Res_delta))>=10^-5),T)
  })

# TODO test for nonidentifiable link that erroniously produces worse fit; add after fix!!!

context("ModelSet reduction")

test_that("modelSet with fixed parameters is reduced correctly", {
  expect_message(selectMinimalModel(modelset))
  red_modelset = selectMinimalModel(modelset)
})

test_that("modelSet with variable parameters is reduced correctly", {
  expect_message(selectMinimalModel(relax_modelset))
  red_modelset = selectMinimalModel(relax_modelset)
})

context("Direct paths extraction")

p1 = c("r_B_E", "r_C_B")
p2 = c("r_B_E^(-1)", "r_B_A")
p3 = c("r_C_B")

test_that("Product of direct paths is correct", {
    expect_equal(mul_path(p1, p3), c("r_B_E","r_C_B","r_C_B"))
})
test_that("Product of direct and inverted paths is correct", {
    expect_equal(mul_path(p1, p2), c("r_C_B","r_B_A"))
})
# Note: This model does not have profile likelihood information
test_that("getDirectPaths works", {
    expect_equal_to_reference(getDirectPaths(model), "model_direct_path.rds", tolerance=1e-4)
    .GlobalEnv$direct_path = getDirectPaths(model)
})
test_that("getDirectPaths works with node removed", {
    expect_equal_to_reference(getDirectPaths(model, c("node1")), "model_merged_direct_path.rds", tolerance=1e-5)
})

test_that("aggregateDirectPaths works", {
    expect_silent( aggregateDirectPaths(list(a=direct_path, b=direct_path)) )
})
test_that("aggregateDirectPaths raises an error if 'direct_paths' is not a list", {
    expect_error( aggregateDirectPaths(c("A->B")), "must be a list" )
})
test_that("aggregateDirectPaths raises an error if 'names(direct_paths)' is NULL", {
    expect_error( aggregateDirectPaths(list(direct_path, direct_path)), "must be non NULL" )
})
test_that("plotModelParameters works", {
    expect_silent( plotModelParameters(model) )
})

