############################ simulation_model.R ################################"
# Functions for the simulation of the models

# Get the target simulations in a MIDAS design-like or list format
# Returns a matrix in a MIDAS measure-like format
# TODO change the control to put the not enough nodes error after the usability control
simulateModel <- function(model_description, targets, readouts = "all") {
    design = model_description$design
    nodes = model_description$structure$names

    # Get the experimental design constraints on the prediction capacity (index + 1 for the C++)
    stim_nodes = design$stim_nodes
    inhib_nodes = design$inhib_nodes
    inhibables = nodes[ 1 + unique(c(inhib_nodes)) ]
    stimulables = nodes[ 1 + unique(c(stim_nodes)) ]
    perID = 1 + unique(c(design$measured_nodes, stim_nodes, inhib_nodes))
    measurables = nodes[ 1 + unique(c(design$measured_nodes)) ]
    measID = 1 + unique(c(design$measured_nodes))

    # Get the names of the nodes in the network to stimulate and inhibit, and the matrix of the perturbation
    if (is.matrix(targets)) {
        # Already in matrix form
        inhibitors = gsub("i$", "", colnames(targets)[grepl("i$", colnames(targets))] )
        stimulators = colnames(targets)[!grepl("i$", colnames(targets))]
        target_matrix = targets
        colnames(target_matrix)[which( grepl("i$", colnames(targets)) )] = inhibitors
        colnames(target_matrix)[which( !grepl("i$", colnames(targets)) )] = stimulators
    } else if (targets[1] == "all") {
        inhibitors = model_description$structure$names[design$inhib_nodes + 1]
        stimulators = model_description$structure$names[design$stim_nodes + 1]
        target_matrix = cbind(design$inhibitor, design$stimuli)
        colnames(target_matrix) = c(inhibitors, stimulators)
    } else if (is.list(targets)) { # TODO distinguish between numeric and character
        # List of perturbation giving nodes names in vectors, TODO
        target_names = unlist(targets)
        inhibitors = unique(c( inhibitors, gsub("i$", "", target_names[grepl("i$", target_names)]) ))
        stimulators = unique(c( stimulators, target_names[!grepl("i$", target_names)] ))
        target_matrix = rep(0, length(c(stimulators, inhibitors)))
        colnames(target_matrix) = c(stimulators, inhibitors)

        for (combination in targets) {
            line = 1; # TODO
        }
    }

    # Set the new experimental design that will be used for the simulation
    # Remove the perturbations that cannot be used
    ## Set the inhibition matrices
    inhib_nodes = c()
    inhibitions = c()
    for (node in inhibitors) {
        if (!(node %in% inhibables)) {
            print(paste0("Node ", node, " is not inhibited in the network and won't be used"))
        } else {
            inhib_nodes = cbind(inhib_nodes, nodes[which(nodes == node)])
            inhibitions = cbind(inhibitions, target_matrix[, which(colnames(target_matrix) == node)])
        }
    }
    ## Set the stimuli matrices
    stim_nodes = c()
    stimulations = c()
    for (node in stimulators) {
        if (!(node %in% stimulables)) {
            print(paste0("Node ", node, " is not stimulated in the network and won't be used"))
        } else {
            stim_nodes = cbind(stim_nodes, nodes[which(nodes == node)])
            stimulations = cbind(stimulations, target_matrix[, which(colnames(target_matrix) == node)])
        }
    }
    ## Set the nodes to be measured
    measured_nodes = c()
    if (readouts == "all") {
        measured_nodes = measurables
    } else {
        for (node in readouts) {
            if (is.character(readouts)) {
                if (!(node %in% nodes)) {
                    print(paste0("The node ", node, " is not in the network."))
                } else if (!(node %in% measurables)) {
                    print(paste0("The node ", node, " cannot be measured with this model."))
                } else {
                    measured_nodes = c(measured_nodes, which(nodes == node)-1)
                }
            } else if (is.numeric(readouts)) { # Consider the R style numeration
                if (!(node %in% 1:length(nodes))) {
                    print(paste0("There are only ", length(nodes), " node in the network."))
                } else if (!(node %in% measID)) {
                    print(paste0("The node ", node, " (", nodes[node], ") cannot be measured with this model."))
                } else {
                    measured_nodes = c(measured_nodes, node)
                }
            }
        }
    }
    new_design = getExperimentalDesign(model_description$structure, stim_nodes, inhib_nodes, measured_nodes, stimulations, inhibitions, model_description$basal)

    # Set up the model and the data for the simulation
    model = new(fitmodel::Model)
    model$setModel( new_design, model_description$structure )
    new_data = new(fitmodel::Data)
    new_data$set_unstim_data(matrix( rep(model_description$data$unstim_data[1,], nrow(target_matrix)), byrow=T, nrow=nrow(target_matrix) ))

    # Compute the predictions
    prediction = list()
    colnames(target_matrix)[which( grepl("i$", colnames(targets)) )] = paste0(inhibitors, "i")
    prediction$conditions = target_matrix
    ## Use the optimal fit
    old_inhib_nodes = model_description$structure$names[1+design$inhib_nodes]
    prediction$bestfit = model$simulate(new_data, getParametersForNewDesign(model, model_description$model, model_description$parameters, old_inhib_nodes, inhib_nodes))$prediction
    colnames(prediction$bestfit) = measured_nodes
    ## Parameters sets provided by the profile likelihood
    params_sets = list()
    if (length(model_description$param_range) == length(model_description$parameters)) {
        idx = 1
        for (i in 1:length(model_description$param_range)) {
            if (!is.na(model_description$param_range[[i]]$low_set[1])) {
                params_sets[[idx]] = getParametersForNewDesign(model, model_description$model, model_description$param_range[[i]]$low_set, old_inhib_nodes, inhib_nodes)
                idx = idx+1
            }
            if (!is.na(model_description$param_range[[i]]$high_set[1])) {
                params_sets[[idx]] = getParametersForNewDesign(model, model_description$model, model_description$param_range[[i]]$high_set, old_inhib_nodes, inhib_nodes)
                idx = idx+1
            }
        }
    }
    ### Predictions for the extra parameter sets
    prediction$variants = list()
    i=1
    for (params in params_sets) {
        prediction$variants = c(prediction$variants, list(model$simulate(new_data, params)$prediction))
        colnames(prediction$variants[[i]]) = measured_nodes
        i=i+1
    }

    rm(model) # Free the memory
    return(prediction)
}

# Give a set of parameter usable by the new model from the parameters fitted in the old model
getParametersForNewDesign <- function(new_model, old_model, old_parameters, old_inhib, inhib_nodes, inhibition=-1, use_fitted_inhib=T) {
    # Get the adjacency matrix and the inhibitions values
    response = old_model$getLocalResponseFromParameter(old_parameters)
    inhib_values = c()
    for (inhibitor in inhib_nodes) {
        if (use_fitted_inhib && inhibitor %in% old_inhib) {
            inhib_values = c(inhib_values, response$inhibitors[which(old_inhib == inhibitor)])
        } else {
            # Possibility to define one inhibition for all or personnalised inhibitions
            if (length(inhibition) == length(inhib_nodes)) {
                cinh = inhibition[which(inhib_nodes == inhibitor)]
            } else {
                cinh = inhibition[1]
            }
            if (cinh > 0) {
                cinh = log2(cinh)
            }
            inhib_values = c(inhib_values, cinh)
        }
    }
    return(new_model$getParameterFromLocalResponse(response$local_response, inhib_values))
}

# Create the perturbation matrix for a set of perturbations, building all n-combinations of stimulators with all m-combinations of inhibitors. Add the cases with only stimulations and only inhibitions
getCombinationMatrix <- function (perturbations, inhib_combo = 2, stim_combo = 1, byStim=T) {
    stimulators = perturbations[!grepl("i$", perturbations)]
    if (stim_combo > length(stimulators) ) {
        stop ("Not enough stimulations to build the combinations")
    }
    inhibitors = perturbations[grepl("i$", perturbations)]
    if (inhib_combo > length(inhibitors) ) {
        stop ("Not enough inhibitions to build the combinations")
    }

    # Create the inhibition matrix
    inhib_combos = build_combo(seq(length(inhibitors)), inhib_combo, c())
    tmp = rep(0, length(inhibitors))
    inhib_matrix = c()
    for (i in 1:nrow(inhib_combos)) {
        inhib_matrix = rbind(inhib_matrix, tmp)
        inhib_matrix[i, inhib_combos[i,]] = 1
    }
    colnames(inhib_matrix) = inhibitors

    # Create the stimulation matrix
    stim_combos = build_combo(seq(length(stimulators)), stim_combo, c())
    tmp = rep(0, length(stimulators))
    stim_matrix = c()
    for (i in 1:nrow(stim_combos)) {
        stim_matrix = rbind(stim_matrix, tmp)
        stim_matrix[i, stim_combos[i,]] = 1
    }
    colnames(stim_matrix) = stimulators

    # Merge the two matrices to get all the combinations
    ## Put a line of 0 to get the inhibition or stimulation alone
    stim_matrix = rbind(rep(0, ncol(stim_matrix)), stim_matrix)
    inhib_matrix = rbind(rep(0, ncol(inhib_matrix)), inhib_matrix)
    ## Combination of perturbations and inhibitions, classified by inhibitions or by stimulations
    perturbation_matrix = c()
    if (byStim) {
        for (i in 1:nrow(stim_matrix)) {
            perturbation_matrix = rbind(perturbation_matrix, cbind( matrix(rep(stim_matrix[i,], nrow(inhib_matrix)), nrow(inhib_matrix), ncol(stim_matrix) , byrow=T), inhib_matrix ))
        }
        colnames(perturbation_matrix) = c(stimulators, inhibitors)
    } else {
        for (i in 1:nrow(inhib_matrix)) {
            perturbation_matrix = rbind(perturbation_matrix, cbind( matrix(rep(inhib_matrix[i,], nrow(stim_matrix)), nrow(stim_matrix), ncol(inhib_matrix) , byrow=T), stim_matrix ))
        }
        colnames(perturbation_matrix) = c(inhibitors, stimulators)
    }
    rownames(perturbation_matrix) = NULL

    return(perturbation_matrix)
}

# Recursively build the n choose k combinations for a set
build_combo <- function (symbols, remaining_steps, to_extend) {
    if (remaining_steps <= 0) {
        return(to_extend)
    } else if (length(symbols) < remaining_steps) {
        return(c())
    }
    final = c()
    # Add each symbol and deepen the recursion with remaining symbols beyond the selected one
    for ( index in seq(length(symbols)) ) {
        extension = cbind(to_extend, symbols[index])
        extension = build_combo(symbols[ -(1:index) ], remaining_steps-1, extension)
        final = rbind(final, extension)
    }
    return(final)
}

# Plot the predictions by the model
# One plot per measured node
# Can plot with error bars if available, and give absolute value or log-fold change
# TODO add the possibility to give log-fold change
plotModelPrediction <- function(model, targets, readouts="all", plotsPerFrame = 4, log_axis=F) {
    if (log_axis) {
        ylog = "y"
    } else {
        ylog = ""
    }
    if (targets[1] == "all") {
        # TODO return(plot_data_simulation)
    }

    prediction = simulateModel(model, targets, readouts)

    ratio = 2/3 # Display height ratio between the plot and its annotation
    layout(matrix(1:2, nrow=2, byrow=T), heights=c(ratio, 1-ratio))
    old_mar = par()$mar
    for (node in 1:ncol(prediction$bestfit)) {
        # Collects the positions of the bars
        par(mar = c(1, 6, 4, 4))
        bars = barplot(prediction$bestfit[,node], plot=F)
        limits = c(0, max(prediction$bestfit[,node]))
        if (length(prediction$variants) > 0) {
            low_var = c()
            high_var = c()
            # Collect the extreme values for each condition, and the global extremes to be sure everything gets included in the plot
            for (perturbation in 1:nrow(prediction$bestfit)) {
                variants = c()
                for (set in 1:length(prediction$variants)) {
                    variants = c(variants, prediction$variants[[set]][perturbation, node])
                }
                low_var = c(low_var, sort(variants)[1])
                # If the incaccuracy yields negative activity, we correct if log scale is used
                if (low_var <= 0 && log_axis) { low_var = 1 }
                high_var = c(high_var, sort(variants, decreasing=T)[1])
                limits = c( min(limits[1], low_var), max(limits[2], high_var) )
            }
            # Plot the bars with the errors
            barplot(prediction$bestfit[,node], ylim=limits, ylab="Activity", log=ylog, col="#008000", main=colnames(prediction$bestfit)[node])
            segments(bars, low_var, bars, high_var)
            space = abs(bars[2] - bars[1])/3
            segments(bars - space, low_var, bars + space, low_var)
            segments(bars - space, high_var, bars + space, high_var)
        } else {
            barplot(prediction$bestfit[,node], ylab="Activity", log=ylog)
            low_var=0;
        }

        # Write the conditions used
        par(mar = c(0, 6, 0, 4), xpd=NA)
        #eplot( xlim=c(0, max(bars)), ylim=c(0, ncol(prediction$conditions)) )
        #barplot(prediction$bestfit[,node], plot=F)
        for (pert in 1:ncol(prediction$conditions)) {
            line = rep("-", nrow(prediction$conditions))
            line[prediction$conditions[, pert] == 1] = "+"
            #line = c(colnames(prediction$conditions)[pert], line)
            text(bars, -pert * limits[2] * 0.9 * ratio / ncol(prediction$conditions), line)
            text(-1 + 3/nrow(prediction$conditions), min(0, low_var)-pert * limits[2] * 0.9 * ratio / ncol(prediction$conditions), colnames(prediction$conditions)[pert], pos=2)
        }
        eplot( xlim=c(0, 1), ylim=c(0, 1) )
    }
    par(mar=old_mar, xpd=T)


    # Invisibly returns the prediction
    return(invisible(prediction))

    # TODO add error multiplier
}

# TODO Function to plot the simulation with the experimental data
# can it be integrated in the other function ?
plot_data_simulation <- function(model, data, plotsPerFrame=4, log_axis=F) {
}

# Plots an empty zone, usefull to write only text
eplot <- function(xlim, ylim, ...) {
    plot(1, xlim=xlim, ylim=ylim, type="n", axes=F, xlab=NA, ylab=NA, ...)
}

