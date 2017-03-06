########################### parameters_comparison.R ########################
# Functions to analyse the parameters of the model

#' Multiply two paths (list of links with exponent)
mul_path <- function(p1, p2) {
    mpath = c(p1, p2)
    pos_links = mpath[!grepl("(-1)", mpath)]
    neg_links = gsub("\\^\\(-1\\)", "", mpath[grep("(-1)", mpath)] )
    remaining = c()
    for (link in neg_links) {
        if (link %in% pos_links) {
            pos_links = pos_links[-grep(link, pos_links)[1]] # Remove one occurence of the link
        } else {
            remaining = c(remaining, paste0(link, "^(-1)"))
        }
    }
    return(c(pos_links, remaining))
}

# Compute the confidence interval without any dependency assumption
product_ci <- function(ci1, ci2) {
    aa = ci1$value
    bb = ci2$value
    ha = ci1$hv
    if (is.null(ha)) { ha = Inf
    } else if (is.na(ha) || is.nan(ha)) { ha = Inf }
    hb = ci2$hv
    if (is.null(hb)) { hb = Inf
    } else if (is.na(hb) || is.nan(hb)) { hb = Inf }
    la = ci1$lv
    if (is.null(la)) { la = -Inf
    } else if (is.na(la) || is.nan(la)) { la = -Inf }
    lb = ci2$lv
    if (is.null(lb)) { lb = -Inf
    } else if (is.na(lb) || is.nan(lb)) { lb = -Inf }
    return(c( aa*bb, min(ha*lb, hb*la, la*lb, ha*hb), max(ha*hb, hb*la, ha*lb, la*lb) ))
}
# Compute the confidence interval of a product of parameters using the information provided by profile likelihood
pl_ci_product <- function(mra_model, p1id, p2id) {
    if (length(mra_model$upper_values)==0 || length(mra_model$lower_values)==0) {
        stop("Profile likelihood results are necessary to compute a confidence interval")
    }
    get_range_product <- function(values1, range1) {
        if (is.na(values1) || is.na(range1) || is.nan(values1) || is.nan(range1)) {
            return(NA)
        } else {
            return(values1 * range1)
        }
    }
    
    values = c( get_range_product(mra_model$upper_values[p1id], mra_model$param_range[[p1id]]$high_set[p2id]),
                get_range_product(mra_model$upper_values[p2id], mra_model$param_range[[p2id]]$high_set[p1id]),
                get_range_product(mra_model$lower_values[p1id], mra_model$param_range[[p1id]]$low_set[p2id]),
                get_range_product(mra_model$lower_values[p2id], mra_model$param_range[[p2id]]$low_set[p1id])
              )
    param_value = mra_model$parameters[p1id]*mra_model$parameters[p2id]
    if (length(which(!is.na(values))) > 0) {
        values = c(values[!is.na(values)], param_value)
    } else {
        values = NA
    }

    return(c( param_value, min(values), max(values) ))
}

#' Compute direct paths for a model
#' 
#' Compute the direct paths for a model and their confidence intervals by combining the paths containing links with negative exponents with paths containing the same link with a positive exponent.
#' @param mra_model An MRAmodel object
#' @param non_stop_nodes A list of nodes that should not begin or end a path, usefull to compare models with different basal activities. If the node is only present as a first or last node of a path, the corresponding paths will be deleted.
#' @return A list of lists indexed by the path name. Each sublist has fields path, value, hv and lv which represent respectively the path equation, the value of the path, the upper bound of the path and the lower bound of the path.
#' @export
getDirectPaths <- function(mra_model, non_stop_nodes=c()) {
    pnames = names(getParametersNames(mra_model))
    rev_links = list()
    paths = list()
    paths[[length(pnames)]] = c()
    final_paths = list()
    # Find the links that have a negative exponent and the path they belong to
    for (pid in 1:length(pnames)) {
        path = pnames[pid]
        links = unlist(strsplit(path, "\\*"))
        paths[[pid]] = list(links=links, value=mra_model$parameters[pid], hv=mra_model$upper_values[pid], lv=mra_model$lower_values[pid])
        denom_pos = grep("(-1)", links)
        if (any(denom_pos)) {
            rev_links[[length(rev_links)+1]] = list( pid=pid, path=gsub("\\^\\(-1\\)", "", links[denom_pos]) )
        } else {
            final_paths[[length(final_paths)+1]] = list(path=path, value=mra_model$parameters[pid], hv=mra_model$upper_values[pid], lv=mra_model$lower_values[pid] )
        }
    }

    for (rev in rev_links) {
        main_path = paths[[rev$pid]]
        all_positived = FALSE
        for (pid in 1:length(paths)) {
            path = paths[[pid]]
            if (all(rev$path %in% path$links)) {
                cip = pl_ci_product(mra_model, rev$pid, pid)
                pmul = mul_path(main_path$links, path$links)
                final_paths[[length(final_paths)+1]] = list( path=paste0(pmul, collapse="*"), value=cip[1], lv=cip[2], hv=cip[3] )
            }
        }
        if (!all_positived) { # Not tested !!
            if (length(rev$path) > 1) {
                for (pp in rev$path) {
                    rev_links[[length(rev_links)+1]] = list(pid=rev$pid, path=pp)
                }
            }
        }
    }
    names(final_paths) = sapply(final_paths, function(X){ simplify_path_name(paste0(X$path, collapse="*")) })

    # Merge links where we want the same nodes
    non_stop_nodes = unique(non_stop_nodes)
    if (length(non_stop_nodes) > 0) {
        merged_paths = list()
        to_merge = grepl( paste0("->", non_stop_nodes, "$|^", non_stop_nodes, "->", collapse="|"), sapply( final_paths, function(X){STASNet:::simplify_path_name(X$path)} ) )
        for (pid in which(!to_merge)) {
            merged_paths[[length(merged_paths)+1]] = final_paths[[pid]]
        }
        for (node in non_stop_nodes) {
            start_paths = c() # Paths that stop with the node
            stop_paths = c() # Paths that start with the node
            for (pid in which(to_merge)) {
                path = STASNet:::simplify_path_name( paste0( final_paths[[pid]]$path, collapse="*" ))
                if (grepl(paste0("->", node, "$"), path) ) {
                    stop_paths = c(stop_paths, pid)
                }
                if (grepl(paste0("^", node, "->"), path) ) {
                    start_paths = c(start_paths, pid)
                }
            }
            for (p1 in stop_paths) {
                for (p2 in start_paths) {
                    composition = c(final_paths[[p1]]$composition, final_paths[[p2]]$composition)
                    # We do pl based confidence interval if possible
                    if (length(composition) == 2) {
                        cip = pl_ci_product(mra_model, composition[1], composition[2])
                    } else {
                        cip = product_ci(final_paths[[p1]], final_paths[[p2]])
                    }
                    mul = mul_path(final_paths[[p1]]$path, final_paths[[p2]]$path)
                    merged_paths[[length(merged_paths)+1]] = list( path=mul, composition=composition, value=cip[1], lv=cip[2], hv=cip[3] )
                }
            }
        }
        final_paths = merged_paths
        names(final_paths) = sapply(final_paths, function(X){ simplify_path_name(paste0(X$path, collapse="*")) })
    }

    return(final_paths)
}

#' Print the direct paths with confidence interval
#' @rdname getDirectPaths
#' @param digits Number of digits to display for the confidence interval values
#' @return Hidden. The direct paths
#' @export
printDirectPaths <- function(mra_model, digits=2) {
    final_paths = getDirectPaths(mra_model)
    sapply(final_paths, function(X){ print(paste0( simplify_path_name(paste0(X$path, collapse="*")), " = ", signif(X$value, digits), " (", signif(X$lv, digits), " - ", signif(X$hv, digits), ")" )) })
    invisible(final_paths)
}

#' Aggregate the paths values for plotting
#'
#' @param direct_paths A list of direct paths (as outputed by getDirectPaths). The names of the elements are considered to be the names of the models.
#' @return A list with two entries. 'paths' is a matrix with the path and the models names as rownames, and columns 'value', 'hv' and 'lv' representing the value with confidence interval of the path. 'model_names' is the name of the models.
#' @export
aggregateDirectPaths <- function(direct_paths){
    if (!is.list(direct_paths)) { stop("'direct_paths' must be a list") }
    if (is.null(names(direct_paths))) { stop("'names(direct_paths)' must be non NULL for the aggregation") }
    aggregated_paths = c()
    for (pname in names(direct_paths)) {
        path = direct_paths[[pname]]
        agg_path = t(sapply(path, function(X){ c(X$value, ifelse(is.null(X$hv), NA, X$hv), ifelse(is.null(X$lv), NA, X$lv)) }))
        rownames(agg_path)=paste0(rownames(agg_path), " ", pname)
        colnames(agg_path)=c("value", "hv", "lv")

        aggregated_paths = rbind(aggregated_paths, agg_path)
    }
    aggregated_paths = aggregated_paths[order(rownames(aggregated_paths)),]

    return( list(paths=aggregated_paths, model_names=names(direct_paths)) )
}

#' Plot parameters from aggregated direct paths
#'
#' Plot parameters from aggregated direct paths. The margins are increased to c(10, 4, 4, 2) but can be increased more by specifying larger margins in 'par' before calling the function.
#' @param aggregated_paths Aggregated direct paths as returned by aggregateDirectPaths
#' @param lim The absolute limit of the plot y-axis, the parameters with value falling outside the range [-lim, lim] will have their value displayed at the edge of the plotting area.
#' @param repar Whether par() should be called, use FALSE for the function to work with layouts.
#' @param resetpar Whether par() should be called after the plotting to restore the previous plot parameters
#' @export
#' @seealso \code{\link{aggregateDirectPaths}}
plotParameters <- function(aggregated_paths, lim=2, repar=TRUE, resetpar=TRUE) {
    if (lim < 0) {
        lim = abs(lim)
        warning("Negative limit provided, using absolute value")
    }
    model_names = aggregated_paths$model_names
    aggregated_paths = aggregated_paths$paths
    aggregated_paths = t( apply(aggregated_paths, 1, function(rr){ rr[is.na(rr)] = rr["value"]; rr }) )

    # Get enough margin in the bottom for the long path names
    if (repar) {
        opar = par()
        par(mar=c(max(10, opar$mar[1]), max(4, opar$mar[2]), max(opar$mar[3], 4), max(opar$mar[4], 2)) + 0.1)
    }
    ymin = max(-lim, min(aggregated_paths, na.rm=TRUE))
    ymax = min(lim, max(aggregated_paths, na.rm=TRUE))
    aggregated_paths[is.na(aggregated_paths[,"lv"])|is.infinite(aggregated_paths[,"lv"]), "lv"] = 1.1 * ymin
    aggregated_paths[is.na(aggregated_paths[,"hv"])|is.infinite(aggregated_paths[,"hv"]), "hv"] = 1.1 * ymax
    # Plot the parameters with error bars
    plot(aggregated_paths[,"value"], xaxt="n", xlab="", pch=20, ylim=c(ymin, ymax), ylab="Path value")
    lines(1:nrow(aggregated_paths), rep(0, nrow(aggregated_paths)), col="grey", lty=2)
    segments(1:nrow(aggregated_paths), aggregated_paths[,"lv"], 1:nrow(aggregated_paths), aggregated_paths[,"hv"], xlab="")
    text(1:nrow(aggregated_paths), par("usr")[3] - (ymax-ymin)/18, srt = 45, adj = 1, labels = rownames(aggregated_paths), xpd = TRUE, cex=0.7)
    axis(1, 1:nrow(aggregated_paths), label=F)
    # Add text for the parameters whose value is outside the limits
    out_up = which(apply(aggregated_paths, 1, function(X) {X["value"] > lim}))
    if (length(out_up) > 0) {
        text(out_up, lim, signif(aggregated_paths[out_up, "value"], 3) )
    }
    out_down = which(apply(aggregated_paths, 1, function(X) {X["value"] < -lim}))
    if (length(out_down) > 0) {
        text(out_down, -lim, signif(aggregated_paths[out_down, "value"], 3) )
    }

    # Check sames paths for different models
    testline = cbind(rownames(aggregated_paths), c(rownames(aggregated_paths)[-1], NA))
    models_test = paste0(model_names, collapse="|")
    line_pos = which(apply(testline, 1, function(rr) { gsub(models_test, "", rr[1])!=gsub(models_test, "", rr[2]) })) + 0.5
    segments(line_pos, ymin, line_pos, ymax, col="gray")
    if (repar && resetpar) { suppressWarnings(par(opar)) }
}

#' Plot the parameters for a model or a list of models
#' @param aggregated_paths Aggregated direct paths as returned by aggregateDirectPaths
#' @param lim The absolute limit of the plot y-axis, the parameters with value falling outside the range [-lim, lim] will have their value displayed at the edge of the plotting area.
#' @param digits Number of digits to display for the confidence interval values
#' @return Invisibly the aggregated direct paths
#' @export
#' @seealso \code{\link{plotParameters}}, \code{\link{getDirectPaths}}, \code{\link{aggregateDirectPaths}}
plotModelParameters <- function(models, non_stop_nodes=c(), lim=2) {
    if ("MRAmodel" %in% class(models)) {
        model_list = list()
        model_list[["_"]] = models
        models = model_list
    }
    if (!is.list(models)) { stop("'models' must be a list") }
    if (is.null(names(models))) { stop("'names(models)' must be non NULL for the aggregation") }

    dpaths = aggregateDirectPaths(lapply(models, getDirectPaths, non_stop_nodes))
    plotParameters(dpaths, lim)
    invisible(dpaths)
}


