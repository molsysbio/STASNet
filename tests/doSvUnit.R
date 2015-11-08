#!/usr/bin/Rscript
require(svUnit)  # Needed if run from R CMD BATCH

pkg <- "fitmodel"
require(fitmodel)

unlink("report.xml")  # Make sure we generate a new report
mypkgSuite <- svSuiteList(pkg, dirs="../inst/unitTest")  # List all our test suites
runTest(mypkgSuite, name = pkg)  # Run them...

protocol(Log(), type = "junit", file = "report.xml")  # ... and write report
