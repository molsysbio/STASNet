
###################### profile_likelihood.R ###############################
# Functions associated with the profile likelihood

#' Computes the profile likelihood and the parameter relationships of each parameter in the model
#' @param model_description An MRAmodel object
#' @param nb_points Number of points to plot the profile
#' @param nb_cores Maximum number of cores used for the calculation
#' @return Returns a list of the profile for each parameters
#' @export
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
profileLikelihood <- function(model_description, nb_points=10000, nb_cores=1) {
    # Get the information from the model description
    model = model_description$model
    model_structure = model_description$structure
    data = model_description$data

    init_params = model_description$parameters
    init_residual = model_description$bestfit
    message(paste(length(init_params), " paths to evaluate"))

    profiles = parallelPL(model, data, init_params, nb_points, nb_cores)

    # Print results
    message(paste("Residual =", init_residual))
    message(paste("Residual score =", model_description$bestfitscore))
    print_error_intervals(profiles)

    return(profiles)
}

# Parallelise the profile likelihood
parallelPL <- function(model, data, init_params, nb_points, NB_CORES=1) {
    if (NB_CORES == 0) {
        NB_CORES = detectCores()-1
    }

    profiles = mclapply(seq(length(init_params)), function(path, model, data, init_params) { profile = model$profileLikelihood(data, init_params, path, nb_points) ; message(paste0("Parameter ", path, "/", length(init_params), " decided")) ; return(profile) }, model, data, init_params, mc.cores=NB_CORES )
    for (path in 1:length(init_params)) {
        # Residuals bigger than the simultaneous threshold are useless and would extend y axis, hidding the information
        # TODO : Integrate in the plot, not here
        profiles[[path]]$residuals[path, profiles[[path]]$residuals[path,] > 1.1 * profiles[[path]]$thresholds[2]] = 1.1 * profiles[[path]]$thresholds[2]
        # Simplify the name of the paths, keep the old one
        profiles[[path]]$old_path = profiles[[path]]$path
        profiles[[path]]$path = simplify_path_name(profiles[[path]]$path)
    }
    names(profiles) = sapply(profiles, function(pp) { pp$path })
    return(profiles)
}

# Simplify the path name from the r_t_f* form to a f->t-> form
simplify_path_name <- function (path_name) {
  simplify_sub_path <- function (sub_path) {
    
    # Prepare a matrix 2*nb_nodes indexed by node names
    entries=unique(unlist(strsplit(sub_path, "\\*|^r_|_|\\^|\\(-1\\)")))
    elements=matrix("",nrow=sum(entries!=""),ncol=2)
    rownames(elements)<-entries[entries!=""]
    
    # Indicate for each node its previous and following node in the path
    for (link in unlist(strsplit(sub_path, "\\*"))){
      nodes = unlist(strsplit(link, "_|\\^"))[2:3]
      elements[nodes[1], 1] = nodes[2]
      elements[nodes[2], 2] = paste0("->", nodes[1])
    }

    # Look for the node(s) without predecessor (i.e the first node of the path(s))
    selected=which(elements[,1]=="")
    if (is.na(selected[1])){
      message(paste0(c("The following path appears to be circular: ", path_name, ". A random node is chosen as the path simplification starting point")))
      selected = 1
      elements[selected,1]=""
    }

    # Build the most simple sub path(s)
    final_paths=matrix("",ncol=length(selected),nrow=1)
    for (i in 1:length(selected)){	
      back_at_start = T
      node = gsub("->", "", elements[selected[i], 2])
      simple_sub_path = paste0(elements[node, 1], elements[selected[i],2])
      while (elements[node, 2] != "" & back_at_start) {
        node = elements[node, 2]
        simple_sub_path = paste0(simple_sub_path, node)
        node = gsub("->", "", node)
        back_at_start = unlist(strsplit(simple_sub_path,"->"))[1] != node
      }
      final_paths[i]=simple_sub_path
    }
    return(paste(final_paths[order(final_paths)],collapse="*"))
  }

  # Inhibitors do not need modification
  if (grepl("i$", path_name) || grepl("^i", path_name)) {
     return(path_name)
  }

  # First derive the sub path of the enumerator 
  links=unlist(strsplit(path_name, "\\*"))
  denom_pos=grepl("(-1)",links)
  simple_path=simplify_sub_path(links[!denom_pos])

  # If present add the sub path of the denominator 
  if (any(denom_pos)){
    simple_path=paste0(simple_path,"/",simplify_sub_path(links[denom_pos]))
  }
    return(simple_path)
}

# Print each path value with its error from the profile likelihood
print_error_intervals <- function(profiles) {
    for (i in 1:length(profiles)) {
        lidx = profiles[[i]]$lower_error_index[1]
        hidx = profiles[[i]]$upper_error_index[1]

        # Print differently if there is non identifiability
        if (profiles[[i]]$lower_pointwise && profiles[[i]]$upper_pointwise) {
            message(paste( profiles[[i]]$path, "=", profiles[[i]]$value, "(", profiles[[i]]$explored[lidx], "-", profiles[[i]]$explored[hidx], ")"))
        } else if (profiles[[i]]$lower_pointwise) {
            message(paste( profiles[[i]]$path, "=", profiles[[i]]$value, "(", profiles[[i]]$explored[lidx], "- ni )"))
        } else if (profiles[[i]]$upper_pointwise) {
            message(paste( profiles[[i]]$path, "=", profiles[[i]]$value, "( ni -", profiles[[i]]$explored[hidx], ")"))
        } else {
            message(paste( profiles[[i]]$path, "=", profiles[[i]]$value, "(non identifiable)"))
        }
    }
}

#' Add the information provided by the profile likelihood in the model_description object
#' (i.e. the limits of the confidence interval and the corresponding sets of coefficients)
#' @param model_description An MRAmodel object
#' @param profiles A list of the likelihood profiles of the coefficients of the model
#' @return An MRAmodel object with the profile likelihood information integrated to it
#' @export
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
#' @seealso \code{\link{profileLikelihood}}
addPLinfos <- function(model_description, profiles) {
    model_description$lower_values = c()
    model_description$upper_values = c()
    model_description$param_range=list()
    for (i in 1:length(profiles)) {
        pr = profiles[[i]]
        model_description$param_range[[i]] = list()
        if (pr$lower_pointwise) {
            lidx = pr$lower_error_index[1]
            model_description$lower_values = c(model_description$lower_values, pr$explored[lidx])
            model_description$param_range[[i]]$low_set = pr$residuals[,lidx]
            model_description$param_range[[i]]$low_set[i] = pr$explored[lidx]
        } else {
            model_description$lower_values = c(model_description$lower_values, NA)
            model_description$param_range[[i]]$low_set = NA
        }
        if (pr$upper_pointwise) {
            hidx = pr$upper_error_index[1]
            model_description$upper_values = c(model_description$upper_values, pr$explored[hidx])
            model_description$param_range[[i]]$high_set = pr$residuals[,hidx]
            model_description$param_range[[i]]$high_set[i] = pr$explored[hidx]
        } else {
            model_description$upper_values = c(model_description$upper_values, NA)
            model_description$param_range[[i]]$high_set = NA
        }
    }
    return(model_description)
}

#' Plots the functional relations between each non-identifiable parameter and the profile likelihood of all parameters
#' and redirect it in a pdf file.
#' @param profiles A list of the likelihood profiles of the parameters of the model
#' @param data_name Name for the output pdf file
#' @param folder Path to the folder for the output pdf file (must end with /)
#' @return Nothing
#' @export
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
#' @seealso \code{\link{profileLikelihood}}
niplotPL <- function(profiles, data_name="default", folder="./", file_plots=TRUE) {
    # Remove residuals bigger than the simultaneous threshold for the plot to prevent an extension of the y axis
    for (pmain in 1:length(profiles)) {
        # Scale the other parameters profiles
        for (pid in (1:length(profiles))[-pmain]) {
            profiles[[pmain]]$residuals[pid,] = profiles[[pmain]]$residuals[pid,] / max(abs(profiles[[pmain]]$residuals[pid,]))
        }
    }
    # Sort the profiles to output differently whether they are identifiable or not
    sorted_profiles = classify_profiles(profiles)
    i_profiles = sorted_profiles[[1]]
    ni_profiles = sorted_profiles[[2]]

# Non identifiables
    nbni = length(sorted_profiles$niid)
    message(paste(nbni, "non identifiable paths"))
    # Compute the dimension, minimal size if there are not enough non identifiables
    dimension = 7
    # Attribute color and style to each parameter
    colors = rep(cbbPalette, length.out=length(profiles))
    styles = rep(c(sapply(1:6, rep, length(cbbPalette))), length.out=length(profiles))

    if (file_plots) {
        pdf(paste0(folder, "NIplot_", data_name, ".pdf"), height=dimension, width=dimension)
    }
    eplot(c(0, 1), c(0, 1))
    legend( 0, 1, sapply(profiles, function(X){X$path}), col=colors[1:length(profiles)], lty=styles[1:length(profiles)], ncol=1 , bty = "n")
    for (plid in 1:length(profiles)) {
        profile = profiles[[plid]]
        th_diff = profile$thresholds[2]-profile$thresholds[1]
        identifiable = plid %in% sorted_profiles$iid
        if (!identifiable) {
            layout(matrix(c(1, 1, 2, 3), ncol=2))
        } else {
            layout(matrix(c(1, 1, 2, 2), ncol=2))
        }
        par( mar=c(5, 2, 1, 1)+0.1, oma=c(0, 2, 2, 0) )
        bfit = min(profile$residuals[plid,])
        xlabel = paste(profile$path, "value")
        plot( profile$explored, profile$residuals[plid,], type="l", main="Likelihood profile", ylim=c(bfit-th_diff/20, profile$thresholds[2] + th_diff/10), lwd=2, xlab=xlabel )
        title(ylab="Residual", xpd=NA)
        lines( profile$explored, rep(profile$thresholds[1], length(profile$explored)), lty=2, col="grey" )
        lines( profile$explored, rep(profile$thresholds[2], length(profile$explored)), lty=2, col="grey" )
        lines( rep(profile$value, 2), c(bfit-th_diff/30, bfit+th_diff/30), col="red" )

        # Functionnal relation profiles
        if (identifiable) {
            plot(0, type="n", xlim=range(profile$explored), ylim=c(-1.1, 1.1), yaxt="n", xlab=xlabel)
            for (pid in (1:length(profiles))[-plid]) {
                lines( profile$explored, profile$residuals[pid,], col=colors[pid], lty=styles[pid], lwd=2)
            }
            title(main="Other paths", xpd=NA)
        } else { # Separate identifiable and non identifiable profiles
            plot(0, type="n", xlim=range(profile$explored), ylim=c(-1.1, 1.1), yaxt="n", xlab=xlabel)
            for (pid in sorted_profiles$niid) {
                if (pid != profile$pathid) {
                    lines( profile$explored, profile$residuals[pid,], col=colors[pid], lty=styles[pid], lwd=2 )
                }
            }
            title(main="Other non identifiable paths", xpd=NA)
            plot(0, type="n", xlim=range(profile$explored), ylim=c(-1.1, 1.1), yaxt="n", xlab=xlabel)
            for (pid in sorted_profiles$iid) {
                lines( profile$explored, profile$residuals[pid,], col=colors[pid], lty=styles[pid], lwd=2 )
            }
            title(main="Identifiable paths", xpd=NA)
        }
        title(main=profile$path, outer=T)
    }
    if (file_plots) {
        dev.off()
    }
}

# Separates the profiles whether they are identifiables or not
classify_profiles <- function (profiles) {
    ni_profiles = list()
    i_profiles = list()
    i_index = c()
    ni_index = c()

    for (plid in 1:length(profiles)) {
        lprofile = profiles[[plid]]
        # Parameters are non identifiable if their profile likelihood does not reach the low threshold on both sides of the minimum 
        if (lprofile$lower_pointwise && lprofile$upper_pointwise) {
            i_profiles[[length(i_profiles)+1]] <- lprofile
            i_index = c(i_index, plid)
        } else {
            ni_profiles[[length(ni_profiles)+1]] <- lprofile
            ni_index = c(ni_index, plid)
        }
    }
    sorted_profiles = list(identifiables=i_profiles, non_identifiables=ni_profiles, iid=i_index, niid=ni_index)
    message("Profiles sorted")
    return(sorted_profiles)
}

#' Save profile likelihood results in a .rds file
#'
#' Save profile likelihood results in a .rds file name "prefix.rds"
#' @param profiles A profile likelihood object as produced by profileLikelihood
#' @param prefix The prefix for the .rds file.
#' @export
#' @rdname profiles_sharing
exportProfiles <- function(profiles, prefix="likelihood_profiles") {
    if (!grepl(".rds$", prefix)) {
        prefix = paste0(prefix, ".rds")
    }
    saveRDS(profiles, prefix)
}

#' Import profiles
#'
#' Import a .rds file and check that it corresponds to the output of a profile likelihood procedure
#' @param fname name of the rds file containing
#' @export
#' @rdname profiles_sharing
importProfiles <- function(fname) {
    if (!grepl(".rds$", fname)) {
        fname = paste0(fname, ".rds")
    }
    profiles = readRDS(fname)
    checkProfiles(profiles)
    return(profiles)
}

#' Check that an object is a valid profile likelihood
#'
#' Check that an object is a valid profile likelihood
#' @export
#' @rdname profiles_sharing
checkProfiles <- function(profiles) {
    if (!is.list(profiles)) {
        stop("Not a profile likelihood object (not a list)")
    }
    lapply(seq_along(profiles), function(pid) {
        pp = profiles[[pid]]
        for (field in c("residuals", "explored", "lower_pointwise", "upper_pointwise", "lower_simultaneous", "upper_simultaneous", "thresholds", "path", "pathid", "value", "lower_error_index", "upper_error_index", "old_path")) {
            if (!field %in% names(pp)) {
                stop(paste0("Not a valid profile likelihood object (field '", field, "' missing in item ", pid, ")"))
            }
        }
        if (!is.matrix(pp$residuals)) {
            stop(paste0("Invalid profile likelihood object (residuals in item ", pid, " is not a matrix)"))
        } else if (!is.character(pp$path)) {
            stop(paste0("Invalid profile likelihood object (path in item ", pid, " is not a character)"))
        } else if (!is.numeric(pp$explored)) {
            stop(paste0("Invalid profile likelihood object (explored in item ", pid, " is not a numeric)"))
        } else if (!length(pp$explored) == ncol(pp$residuals)) {
            stop(paste0("Invalid profile likelihood object (Incoherent 'explored' and 'residuals' dimensions in item ", pid, ")"))
        }
    })
}
