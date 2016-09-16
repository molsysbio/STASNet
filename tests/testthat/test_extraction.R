context("Data extractions functions")

expect_structure <- function(structure, id) {
    expect_equal_to_reference(structure$names, paste0("structure_names_", id, ".rds"))
    expect_equal_to_reference(structure$adjacencyMatrix, paste0("structure_adjacencyMatrix_", id, ".rds"))
}
dumb_structure = rbind(c("N1", "N2"), c("N2", "N3"), c("N2", "N4"), c("N3", "N4"))
dumb_structure_bis = matrix(0, ncol=4, nrow=4, dimnames=list(paste0("N", 1:4), paste0("N", 1:4)))
dumb_structure_bis[2, 1] = 1
dumb_structure_bis[3, 2] = 1
dumb_structure_bis[4, 2] = 1
dumb_structure_bis[4, 3] = 1
test_that("extractStructure behaves as expected", {
    expect_structure(extractStructure(dumb_structure), 0) # two columns matrix
    expect_structure(extractStructure(dumb_structure_bis), 0) # n*n adjacency matrix
    colnames(dumb_structure_bis) = rownames(dumb_structure_bis) = NULL
    expect_warning(extractStructure(dumb_structure_bis)) # n*n adjacency matrix without node names
    expect_structure(suppressWarnings(extractStructure(dumb_structure_bis)), 1) # n*n adjacency matrix without node names
})

dumb_midas = matrix(0, nrow=5, ncol=6, dimnames=list(NULL, c("ID:type", "TR:N1", "TR:N2i", "DA:ALL", "DV:N3", "DV:N4")))
dumb_midas = as.data.frame(dumb_midas)
test_that("extractMIDAS behaves as expected", {
    expect_equal(extractMIDAS(dumb_midas), dumb_midas)
    expect_error(extractMIDAS(dumb_midas[,-1]), "'ID:type' is missing")
    expect_error(extractMIDAS(dumb_midas[,c(-2, -3)]), "'TR' field is missing")
    expect_error(extractMIDAS(dumb_midas[,c(-5, -6)]), "'DV' field is missing")
})

dumb_activity = c("N1", "N2")
test_that("extractBasalActivity behaves as expected", {
    expect_equal(extractBasalActivity(dumb_activity), dumb_activity)
    expect_equal(extractBasalActivity("dumb_basal.dat"), dumb_activity)
    suppressWarnings(expect_error(extractBasalActivity("does_not_exist"), "cannot open the connection"))
})

context("Model creation with R objects")
dumb_midas[,1] = c("c", "t", "t", "t", "t")
dumb_midas[,c(2,3)] = cbind(c(0,1,1,0,0), c(0,0,1,0,1)) # Perturbations
dumb_midas[,c(5,6)] = cbind(c(1, 2, 1.5, 1, 0.5), c(1, 4, 2, 1, 0.5))
no_control_midas = dumb_midas[-1,]
no_perturbations_midas = dumb_midas[,c(-2, -3)]
dumb_variation = dumb_midas
dumb_variation[,c(5,6)] = 0.1
test_that("createModel works with R objects", {
    expect_silent( suppressMessages(createModel(dumb_structure, dumb_activity, dumb_midas, dumb_variation, inits=10)) ) # With error model
    expect_silent( suppressMessages(createModel(dumb_structure, dumb_midas, dumb_variation, inits=10)) ) # Without error model
    expect_error(createModel(dumb_structure, dumb_activity, no_control_midas, inits=10)) # Missing control data
    expect_error(createModel(dumb_structure, dumb_activity, no_perturbations_midas, inits=10)) # Missing perturabtion
})

test_that("createModelSet works with R objects", {
    expect_silent( suppressMessages(createModelSet(dumb_structure, dumb_activity, list(m1=dumb_midas, m2=dumb_midas), list(m1=dumb_variation, m2=dumb_variation), model_name=c("m1", "m2"), inits=10)) )
})

context("Limit cases for createModel")

test_that("Limit cases number of initialisations behave as expected", {
    expect_error(suppressMessages(createModel(dumb_structure, dumb_activity, dumb_midas, inits=0)), "Number of initialisations") # No initialisation
    #expect_silent(suppressMessages(createModel(dumb_structure, dumb_activity, dumb_midas, inits=1))) # Only one initialisation
})
test_that("Misused arguments raise error", {
    expect_error(suppressMessages(createModel(dumb_activity, dumb_midas, inits=10))) # Argument forgotten
    expect_error(suppressMessages(createModel(dumb_activity, dumb_structure, dumb_midas, inits=10))) # Arguments inverted
})

context("No inhibition or simulations")

only_stim = dumb_midas[, -3]
only_inhib = dumb_midas[, -2]
test_that("Only inhibitions works", {
    expect_silent(suppressMessages(createModel(dumb_structure, dumb_activity, only_inhib, inits=1)))
})
test_that("Only inhibitions works in extractModelCore", {
    expect_silent(extractModelCore(dumb_structure, dumb_activity, only_inhib))
})
test_that("Only stimulations works", {
    expect_silent(suppressMessages(createModel(dumb_structure, dumb_activity, only_stim, inits=1)))
})
test_that("Deleting all inhibitions work", {
    expect_silent(suppressMessages( createModel(dumb_structure, dumb_activity, dumb_midas, inits=1, unused_perturbations=c("N2i")) ))
})
test_that("Deleting all stimulations work", {
    expect_silent(suppressMessages( createModel(dumb_structure, dumb_activity, dumb_midas, inits=1, unused_perturbations=c("N1")) ))
})

context("Model core data extraction")

test_that("extractModelCore works as expected", {
    expect_silent(suppressMessages(extractModelCore(dumb_structure, dumb_activity, dumb_midas)))
})
test_that("extractModelCore error when perturbations are missing", {
    expect_error(suppressMessages(extractModelCore(dumb_structure, dumb_activity, no_perturbations_midas)))
})

context("Helper functions")

alph=c("a", "b", "c", "d")
num = 1:4
combine = paste0(alph, num)
test_df = data.frame(list(a=num, b=alph, c=alph), row.names=combine)

test_that("Select one row using sub_data_frame", {
    expect_equal( sub_data_frame(test_df, 1), data.frame(a=test_df[1,1], b=test_df[1,2], c=test_df[1,3], row.names=combine[1]) )
#a=test_df$a[1], b=test_df$a[1], c=test_df$a[1] -> gives numbers
#a=test_d[1,1], b=test_df[1,2], c=test_df[1,3] -> gives factors
})
test_that("Select one row by name with sub_data_frame", {
    expect_equal( sub_data_frame(test_df, "a1"), data.frame(a=test_df[1,1], b=test_df[1,2], c=test_df[1,3], row.names=combine[1]) )
})
test_that("Select several rows with sub_data_frame", {
    expect_equal( sub_data_frame(test_df, c("a1", "c3")), test_df[c("a1", "c3"),] )
})
test_that("Select one colum using sub_data_frame", {
    expect_equal( sub_data_frame(test_df, cols=1), data.frame(a=num, row.names=combine) )
})
test_that("Select one colum by name with sub_data_frame", {
    expect_equal( sub_data_frame(test_df, cols="a"), data.frame(a=num, row.names=combine) )
})
test_that("Select several colum using sub_data_frame", {
    expect_equal( sub_data_frame(test_df, cols=c("a", "b")), test_df[, c("a", "b")] )
})
