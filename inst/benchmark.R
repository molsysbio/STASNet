#!/usr/bin/Rscript
#-*- coding: utf8 -*-
suppressPackageStartupMessages(library(fitmodel))
# Benchmarking of the fitmodel package on a network structure

if (!exists("cargs")) {
    cargs = commandArgs()
} else if (is.character(cargs)) {
    cargs = strsplit(cargs, " ")[[1]]
}

for (arg in cargs) {
    if (grepl("^--net", arg) || grepl(".tab$", arg)) {
        network = gsub("^--net", "", arg)
        if (is.character(network) && length(network) == 1) {
            dname = gsub(".tab$", "", network)
        }
    } else if (grepl("^-x", arg)) {
        repetitions = as.numeric(gsub("^-x", "", arg))
    }
}
if (!exists("network")) { network = stop("A network file must be provided") }
if (!exists("repetitions") || repetitions < 1) { repetitions = 10 }

sim_title = paste0("benchmark_", repetitions, "_", dname)

true_adm = fitmodel:::readNetworkAdj(network)
true_inhibitions = -1
pdf(paste0(sim_title, ".pdf"))
for (ii in c(0.01, 0.05, 0.1, 0.2, 0.3)) {
    print(paste0("Benchmarking for sd = ", ii))
    noise_lvl = ii
    values = c()
    for (jj in 1:repetitions) {
        source("~/bin/toy_data_generation")
        mra = createModel(network, basal_file, simulation_file, nb_cores=0)
        local_response = mra$model$getLocalResponseFromParameter(mra$parameters)
        adm = local_response$local_response
        values = rbind( values, c((adm[mra$structure$adjacencyMatrix != 0]-true_adm[mra$structure$adjacencyMatrix != 0])/abs(true_adm[mra$structure$adjacencyMatrix != 0]), local_response$inhibitors - true_inhibitions) )
    }
    plot(apply( cbind(1:ncol(values)), 1, rep, nrow(values) ), values, xlab="Parameters", ylab="Relative difference to truth", xaxt="n", main=paste0("CV = ", ii, ", ", nrep, " replicates"), ylim=c(-1, 1), pch=16)
    axis(1, at=1:ncol(values), label=getParametersNames(mra), las=3)
    lines(0:(ncol(values)+1), rep(0, ncol(values)+2), col="red")
}
plotNetworkGraph(mra$structure, mra$design)
dev.off()


