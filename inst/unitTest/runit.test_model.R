

require(svUnit)
require(STASNet)


test.model.regression <- function() {
  
  model=createModel("inst/unitTest/adjacency.tab",
    "inst/unitTest/MIDAS.csv",
    "inst/unitTest/basal_activity.dat",inits=1000,method="correlation")

  checkEquals(length(model$parameters), 7)
}




test.model.random <- function() {
  
  model=createModel("inst/unitTest/adjacency.tab",
    "inst/unitTest/MIDAS.csv","inst/unitTest/basal_activity.dat",
    inits=1000,method="random")

  checkEquals(length(model$parameters), 7)
}
