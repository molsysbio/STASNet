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

# Compute the confidence interval of a product of parameters using the information provided by profile likelihood
pl_ci_product <- function(mra_model, p1id, p2id) {
    if (length(mra_model$upper_values)==0 || length(mra_model$lower_values)==0) {
        stop("Profile likelihood results are necessary to compute a confidence interval")
    }
    get_range_product <- function(values1, range1) {
        print(values1)
        print(range1)
        if (is.na(values1) || is.na(range1) || is.nan(values1) || is.nan(range1)) {
            return(NA)
        } else {
            return(values1 * range1)
        }
    }
    
    print(p1id)
    print(p2id)
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
#' @value A list
#' @export
getDirectPaths <- function(mra_model) {
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
            final_paths[[length(final_paths)+1]] = list(path=path, composition=pid, value=mra_model$parameters[pid], hv=mra_model$upper_values[pid], lv=mra_model$lower_values[pid] )
        }
    }

    for (rev in rev_links) {
        main_path = paths[[rev$pid]]
        all_positived = FALSE
        for (pid in 1:length(paths)) {
            path = paths[[pid]]
            if (all(rev$path %in% path$links)) {
                cip = pl_ci_product(mra_model, rev$pid, pid)
                print(cip)
                pmul = mul_path(main_path$links, path$links)
                final_paths[[length(final_paths)+1]] = list( path=pmul, composition=c(rev$pid, pid),
                                                            value=cip[1], lv=cip[2], hv=cip[3] )
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
    return(final_paths)
}

#' Print the direct paths with confidence interval
#' @rdname getDirectPaths
#' @param digits Number of digits to display for the confidence interval values
printDirectPaths <- function(mra_model, digits=2) {
    final_paths = getDirectPaths(mra_model)
    sapply(final_paths, function(X){ print(paste0( simplify_path_name(paste0(X$path, collapse="*")), " = ", signif(X$value, digits), " (", signif(X$lv, digits), " - ", signif(X$hv, digits), ")" )) })
    invisible(final_paths)
}
