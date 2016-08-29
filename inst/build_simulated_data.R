#!/usr/bin/Rscript
#-*- coding: utf8 -*-
suppressPackageStartupMessages(library(STASNet))

if (!exists("cargs")) {
    cargs = commandArgs()
} else if (is.character(cargs)) {
    cargs = unlist(strsplit(cargs, " "))
}

dname = "test"

for (arg in cargs) {
    if (grepl("^--net", arg) || grepl(".tab$", arg)) {
        network = gsub("^--net", "", arg)
        if (is.character(network) && length(network) == 1) {
            dname = gsub(".tab$", "", network)
        }
    } else if (grepl("^-r", arg)) {
        nrep = as.numeric(gsub("^-r", "", arg))
        if (is.na(nrep)) { nrep=3 }
    } else if (grepl("^-n", arg)) {
        noise_lvl = as.numeric(gsub("^-n", "", arg))
        if (is.na(noise_lvl)) { noise_lvl=0 }
    } else if (grepl("^-p", arg) || grepl("design|perturbation", arg)) {
        perturbations = gsub("^-p", "", arg)
    } else if (grepl("^-m", arg) || grepl("measured", arg)) {
        measured = gsub("^-m", "", arg)
    } else if (grepl("^-t", arg)) {
        sim_title = gsub("^-t", "", arg)
    } else if (grepl("^-", arg)) {
        print(paste0("Unknown argument: '", arg, "'"))
    }
}
if (!exists("sim_title")) { sim_title = "" }
if (!exists("network")) { network = stop("A network file must be provided") }
if (!exists("perturbations")) { perturbations = "" }
if (!exists("measured")) { measured = "" }
if (!exists("noise_lvl")) { noise_lvl = 0 }
if (!exists("nrep")) { nrep = 3 }

res = createSimulation( network, perturbations, measured, noise=noise_lvl, replicates=nrep )
# Add the columns necessary for MIDAS format validity to the simulation before writing it in a file
sim = res$simulation
cid = grep("^DV", colnames(sim))[1]
types = rep("t", nrow(sim))
types[apply(sim[,1:(cid-1)], 1, function(X){all(X==0)} )] = "c"
sim = cbind( types, sim[,1:(cid-1)], rep(0, nrow(sim)), sim[,cid:ncol(sim)] )
sim = rbind( c("blank", rep(0, ncol(sim)-1)), sim )
colnames(sim)[c(1, cid+1)] = c("ID:type", "DA:ALL")

if (sim_title != "") {
    sim_names = sim_title
} else if (noise_lvl > 0) {
    sim_names = paste0(dname, "_n", noise_lvl, "_r", nrep, "_", Sys.Date())
} else {
    sim_names = dname
}
simulation_file = paste0("simulation_", sim_names, ".csv")
write.csv(sim, simulation_file, row.names=F, quote=F)

basal_file = paste0("basal_", sim_names, ".dat")
write(res$model$structure$names, basal_file)

