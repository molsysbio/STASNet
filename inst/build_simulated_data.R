#!/usr/bin/Rscript
#-*- coding: utf8 -*-
suppressPackageStartupMessages(library(fitmodel))

if (!exists("cargs")) {
    cargs = commandArgs()
} else if (is.character(cargs)) {
    cargs = strsplit(cargs, " ")[[1]]
}

dname = "test"

for (arg in cargs) {
    if (grepl("^-n", arg) || grepl(".tab$", arg)) {
        network = gsub("^-n", "", arg)
        if (is.character(network) && length(network) == 1) {
            dname = gsub(".tab$", "", network)
        }
    } else if (grepl("^-r", arg)) {
        nrep = as.numeric(gsub("^-r", "", arg))
    } else if (grepl("^--noise", arg)) {
        noise_lvl = as.numeric(gsub("^--noise", "", arg))
    } else if (grepl("^-p", arg) || grepl("design|perturbation", arg)) {
        perturbations = gsub("^-p", "", arg)
    } else if (grepl("^-m", arg) || grepl("measured", arg)) {
        measured = gsub("^-m", "", arg)
    } else {
        print(paste0("Unknown argument: '", arg, "'"))
    }
}
if (!exists("network")) { network = stop("A network file must be provided") }
if (!exists("perturbations")) { perturbations = "" }
if (!exists("measured")) { measured = "" }
if (!exists("noise_lvl")) { noise_lvl = 0 }
if (!exists("nrep")) { nrep = 1 }

res = createSimulation( network, perturbations, measured, noise=noise_lvl, replicates=nrep )
# Add the columns necessary for MIDAS format validity to the simulation before writing it in a file
sim = res$simulation
cid = grep("^DV", colnames(sim))[1]
types = rep("t", nrow(sim))
types[apply(sim[,1:(cid-1)], 1, function(X){all(X==0)} )] = "c"
sim = cbind( types, sim[,1:(cid-1)], rep(0, nrow(sim)), sim[,cid:ncol(sim)] )
sim = rbind( c("blank", rep(0, ncol(sim)-1)), sim )
colnames(sim)[c(1, cid+1)] = c("ID:type", "DA:ALL")

sim_names = paste0(dname, "_n", noise_lvl, "_r", nrep)
write.csv(sim, paste0("simulation_", sim_names, ".csv"), row.names=F, quote=F)

write(res$model$structure$names, paste0("basal_", sim_names, ".dat"))

