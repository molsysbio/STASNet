#!/usr/bin/Rscript
# -*- coding:utf-8 -*-

# Hidden from the R installer but with the other scripts from the package

#out = capture.output(try({
# Print the time it took in a readable format
get_running_time <- function(init_time, text="") {
    run_time = proc.time()["elapsed"]-init_time
    run_hours = run_time %/% 3600;
    run_minutes = (run_time - 3600 * run_hours) %/% 60;
    run_seconds = run_time - 3600 * run_hours - 60 * run_minutes;
    message(paste(run_hours, "h", run_minutes, "min", run_seconds, "s", text))
}

# Create a model from the data and fit a minimal model
# Takes relative paths as arguments in the order network data basal

data = ""
network = ""
basal_nodes = ""
variation = ""
nb_steps = 1000
inits = 1000
method = "geneticlhs"
# Autodetection of the cores
cores = 0
limit = Inf
plot_accuracy = FALSE
show_values = FALSE
precorrelate = TRUE
reduction = FALSE
perform_pl = FALSE
mr_pl = FALSE
perf_plots = TRUE
extension = FALSE

if (!exists("cargs")) {
    cargs = commandArgs(trailingOnly=T)
} else if (is.character(cargs)) {
    cargs = strsplit(cargs, " ")[[1]]
}
unused_perturbations = c()
unused_readouts = c()
recomputing = FALSE

# Collect the filenames based on their extension
for (argument in cargs) {
    if (grepl("--help|-h", argument)) {
        message("Help for STASNet fitting script:")
        message("  Initial fitting:")
        message("    fit_script.R data.csv [data.var] network.tab basal.dat [options]")
        message("  Reloading a fitted network:")
        message("    fit_script.R fit.mra [data.csv] [data.var] [options]")
        message("The script expect a .csv file with data in MIDAS format, a .tab file with the network structure and a .dat file with the nodes with basal activity.")
        message("A .var file with the error in MIDAS format can also be provided.")
        message("    --help | -h                  Displays help")
        message("    -i<int>                      Number of initialisations")
        message("    -c<int>                      Maximum of cores to use (0 for auto-detection)")
        message("    --mr | --reduce              Apply model reduction")
        message("    --ext | --extension          Compute possible extensions to the network")
        message("    -m<string>                   Method to apply for the initialisation")
        message("    --nopl                       Disable profile likelihood")
        message("    --pl                         Enable profile likelihood")
        message("    --mrpl                       Perform model reduction and profile likelihood of the reduced model")
        message("    -s<int>                      Number of steps for the profile likelihood")
        message("    --noplots                    Cancel plot generation")
        message("    -l                           Force the limit of the gradients for the heatmaps")
        message("    --values                     Print values in the -l accuracy heatmaps, only used if -l is specified")
        message("    -v                           Activate debug")
        message("    -D<float>                    Default coefficient of variation")
        message("    -M<float>                    Minimum coefficient of variation")
        message("    -u'<string1> <string2> ...'  List of perturbations to ignore")
        message("    -d'<string1> <string2> ...'  List of readouts to ignore")
        quit()
    } else if (grepl(".tab$", argument)) {
        network = paste0(getwd(), "/", argument)
    } else if (grepl(".data$", argument) || grepl(".csv$", argument)) {
        data = paste0(getwd(), "/", argument)
        data_name = basename(argument)
    } else if (grepl(".mra$", argument)) {
        mra_file = paste0(getwd(), "/", argument)
        mra_name = gsub(".mra$", "", basename(mra_file))
        recomputing = TRUE
    } else if (grepl(".dat$", argument)) {
        basal_nodes = paste0(getwd(), "/", argument)
    } else if (grepl(".var$", argument)) {
        variation = paste0(getwd(), "/", argument)
    } else if (grepl("^-i", argument)) {
        power = c("k", "M", "G", "T", "P", "Y");
        inits = gsub("-i", "", argument)
        last_char = substr(inits, nchar(inits), nchar(inits))
        if (last_char %in% power) {
            inits = as.numeric(substr(inits, 1, nchar(inits)-1)) * 10^(3*which(power==last_char))
        } else {
            inits = as.numeric(gsub("-i", "", argument))
        }
        if (is.na(inits)) { # If error
            inits = 1000
            print("Incorrect number of initialisation, using 1000 instead")
        }
    } else if (grepl("^-c", argument)) {
        cores = as.numeric(gsub("-c", "", argument))
        if (is.na(cores)) { # If error
            cores = 2
            stop("Incorrect number of cores (use -c#)")
        }
    } else if (grepl("^-s", argument)) {
        nb_steps = as.numeric(gsub("^-s", "", argument))
        if (is.na(nb_steps)) {
            nb_steps = 1000
            print("Incorrect number of steps, performing with 1000")
        }
    } else if (grepl("^--mr$", argument) || grepl("^--reduce$", argument)) {
        reduction = TRUE
    } else if (grepl("^--ext$", argument) || grepl("^--extension$", argument)) {
        extension = TRUE
    } else if (grepl("^-m", argument)) {
        method = gsub("^-m", "", argument)
    } else if (argument == "--pl") {
        perform_pl = TRUE
    } else if (argument == "--mrpl") {
        mr_pl = TRUE
        reduction = TRUE
    } else if (argument == "--nopl") {
        perform_pl = FALSE
    } else if (argument == "--noplots" || argument == "--noplot") {
        perf_plots = FALSE
    } else if (grepl("^-l", argument)) {
        limit = as.numeric(gsub("^-l", "", argument))
        plot_accuracy = TRUE
    } else if (grepl("^--values", argument)) {
        show_values = TRUE
    } else if (grepl("^-v", argument)) {
        STASNet:::setDebug(T)
    } else if (grepl("^--npc", argument)) {
        precorrelate = FALSE
    } else if (grepl("^-D", argument)) {
        default_cv = as.numeric(gsub("^-D", "", argument))
    } else if (grepl("^-M", argument)) {
        min_cv = as.numeric(gsub("^-M", "", argument))
    } else if (grepl("^-u", argument)) {
        argument = gsub("^-u", "", argument)
        argument = gsub("\"", "", argument)
        unused_perturbations = c( unused_perturbations, unlist(strsplit(argument, " |\t")) )
    } else if (grepl("^-d", argument)) {
        argument = gsub("^-d", "", argument)
        argument = gsub("\"", "", argument)
        unused_readouts = c( unused_readouts, unlist(strsplit(argument, " |\t")) )
    } else if (grepl("^-", argument)) {
        print(paste0("Unknown argument: '", argument, "'"))
    }
}
library("STASNet")

if (cores == 0) {
    cores = detectCores() - 1;
}

if (recomputing) {
    if (data != "") {
        model = rebuildModel(mra_file, data, variation)
    } else {
        stop("Can't rebuild a model without data")
    }
    folder = "./"
    conditions = paste0(mra_name)
    init_time = proc.time()["elapsed"];
} else {
    if (network == "") {
        stop("A network structure (adjacency list, .tab) file is required.")
    } else if (data == "") {
        stop("A data file (.data) is required")
    } else if (basal_nodes == "") {
        stop("A basal activity (.dat) list file is required")
    }

    # Extract the name and the number of initialisations
    power = c("", "k", "M", "G", "T", "P", "Y");
    power_init = floor(log(inits, base=1000))
    conditions = paste0( gsub("(_MIDAS)?.(csv|data)", "", basename(data_name)), "_", gsub(".tab", "", basename(network)), "_", inits%/%(1000^power_init), power[1+power_init]);
    conditions = gsub(" ", "_", conditions)
    if (length(unused_perturbations) > 0) {
        conditions = paste0(conditions, "_no", paste0(unused_perturbations, collapse="-") )
    }
    if (length(unused_readouts) > 0) {
        conditions = paste0(conditions, "_ur", paste0(unused_readouts, collapse="-") )
    }
    if (exists("min_cv")) {
        conditions = paste0(conditions, "_mincv", min_cv)
    } else { min_cv = 0.1 }
    if (exists("default_cv")) {
        conditions = paste0(conditions, "_defaultcv", default_cv)
    } else { default_cv = 0.3 }
    folder = paste0( "run_", conditions, "_", Sys.Date(), "/" )
    dir.create(folder)

    #### Creates the model from network and basal files and fits a minimal model to the data
    init_time = proc.time()["elapsed"];
    pdf(file.path(folder, paste0("distribution_", conditions, ".pdf")))
    model = createModel(network, basal_nodes, data, variation, inits=inits, nb_cores=cores, perform_plots=perf_plots, method=method, precorrelate=precorrelate, unused_perturbations=unused_perturbations, unused_readouts=unused_readouts, MIN_CV=min_cv, DEFAULT_CV=default_cv, rearrange="no");
    dev.off()
    get_running_time(init_time, paste("to build the model with", inits, "initialisations."))

    # Plot the graph of the network in a pdf
    pdf(file.path( folder, paste0("graph_", gsub(" ", "_", gsub(".tab$", ".pdf", basename(network)) ) )))
    plotModelGraph(model)
    dev.off()

    mat=model$data$stim_data
    pdf(file.path(folder, paste0("accuracy_heatmap_", conditions, ".pdf")),onefile=T,width =5+ncol(mat)/3,height=4+nrow(mat)/6)
    plotModelAccuracy(model)
    plotModelScores(model, main=paste0("Global R = ", model$bestfitscore))
    dev.off()
    printParameters(model)
}
if (plot_accuracy) {
    mat=model$data$stim_data
    pdf(file.path(folder, paste0("l", limit, "_accuracy_heatmap_", conditions, ".pdf")),onefile=T,width =5+ncol(mat)/3,height=4+nrow(mat)/6)
    plotModelAccuracy(model, limit, show_values)
    dev.off()
}

if (method == "annealing") {
    stop("debugging annealing")
}

# Perform the profile likelihood
if (perform_pl) {
    conditions = paste0(conditions, "+pl")
    profiles = profileLikelihood(model, nb_steps, nb_cores=min(cores, length(model$parameters)));
    model = addPLinfos(model, profiles);
    get_running_time(init_time, paste("to run the program with", nb_steps, "points for the profile likelihood."));
    niplotPL(profiles, data_name=conditions, folder=folder)
}
exportModel(model, paste0(folder, conditions, ".mra"))

# Plot the simulation for all combinations of inhibitors 
#pdf(paste0("combos_", conditions, ".pdf"))
#plotModelSimulation(model, getCombinationMatrix(c("MEKi", "GSK3ABi", "IGF", "TGFA", "PI3Ki")))
#dev.off()
# Plot the simulated conditions
pdf(file.path(folder, paste0("model_prediction_", conditions, ".pdf")))
plotModelSimulation( model_description = model, with_data = TRUE, log_axis = TRUE)
dev.off()

if (reduction) {
    conditions = paste0(conditions, "+red")
    # Reduce the model and see what changed
    print("Reduction of the model...")
    reduced_model = selectMinimalModel(model)
    # Profile likelihood on the reduced model
    if (mr_pl) {
        reduced_profiles = profileLikelihood(reduced_model, nb_steps, nb_cores=min(cores, length(reduced_model$parameters)));
        reduced_model = addPLinfos(reduced_model, reduced_profiles)
        if (perf_plots) {
            niplotPL(reduced_profiles, data_name=paste0("reduced_", data_name))
        }
    }
    exportModel(reduced_model, paste0(folder, "reduced_", conditions, ".mra"));
    # Plot the simulated conditions
    pdf(file.path(folder, paste0("reduced_all_", conditions, ".pdf")))
    dev.off()

    pdf(file.path(folder, paste0("reduced_accuracy_", conditions, ".pdf")))
    plotModelAccuracy(reduced_model)
    plotModelScores(reduced_model, main=paste0("Global R = ", reduced_model$bestfitscore))
    dev.off()
    pdf(file.path(folder, paste0("reduced_graph_", conditions, ".pdf")))
    plotModelGraph(reduced_model)
    dev.off()

    get_running_time(init_time, "with the model reduction");
}

if (extension) {
    sug_ext = suggestExtension(model, T, cores)
    write.table(sug_ext, file.path(folder, paste0("extension_", conditions, ".csv")), row.names=FALSE, quote=FALSE, sep="\t")
}

print("Finished")

# IDEAS :
# data as last argument, possibility to give severals, in which case a comparison of the models is also done
# one argument : name of a file with network, basal_nodes, and the name of the data files

# Display the warnings if there are some
warnings()

#}), type="message")
#print(out)
