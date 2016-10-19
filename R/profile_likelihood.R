###################### profile_likelihood.R ###############################
# Functions associated with the profile likelihood

#' Computes the profile likelihood and the parameters relationships of each parameters in the model
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
    return(profiles)
}

# Simplify the path name from the r_t_f* form to a f->t-> form
simplify_path_name <- function (path_name) {
  simplify_sub_path <- function (sub_path) {
    
    # Prepare a matrix 2*nb_nodes indexed by node names
    entries=unique(unlist(strsplit(sub_path, "\\*|_|^r|\\^|\\(-1\\)")))
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
      stop(writeLines(c("Error:  the following path appears to be circular: ", path_name, "Consider a different network structure!!!")))
    }

    # Build the most simple sub path(s)
    final_paths=matrix("",ncol=length(selected),nrow=1)
    for (i in 1:length(selected)){	
      node = gsub("->", "", elements[selected[i], 2])
      simple_sub_path = paste0(elements[node, 1], elements[selected[i],2])
      while (elements[node, 2] != "") {
        node = elements[node, 2]
        simple_sub_path = paste0(simple_sub_path, node)
        node = gsub("->", "", node)
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
#' ie the limits of the confidence interval and the corresponding sets of parameters
#' @param model_description An MRAmodel object
#' @param profiles A list of the likelihood profiles of the parameters of the model
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

#' Plots the functionnal relations between each non identifiable parameter and the profile likelihood of all parameters
#' and redirect it in a pdf file.
#' @param profiles A list of the likelihood profiles of the parameters of the model
#' @param data_name Name for the output pdf file
#' @param folder Path to the folder for the output pdf file (must end with /)
#' @return Nothing
#' @export
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
#' @seealso \code{\link{profileLikelihood}}
niplotPL <- function(profiles, data_name="default", folder="./") {
    # Remove residuals bigger than the simultaneous threshold for the plot to prevent an extension of the y axis
    for (profile in profiles) {
        residual_limit = 1.1 * (profile$thresholds[2] - profile$thresholds[1]) + profile$thresholds[2]
        profile$residuals[ profile$pathid, profile$residuals[profile$pathid,] > residual_limit ] = residual_limit
    }
    # Sort the profiles to output differently whether they are identifiable or not
    sorted_profiles = classify_profiles(profiles)
    i_profiles = sorted_profiles[[1]]
    ni_profiles = sorted_profiles[[2]]

# Non identifiables
    nbni = length(ni_profiles)
    message(paste(nbni, "non identifiable paths"))
    # Compute the dimension, minimal size if there are not enough non identifiables
    if (nbni > 3) {
        dimension = 2 * nbni + 1
    }
    else {
        dimension = 7
    }

    pdf(paste0(folder, "NIplot_", data_name, ".pdf"), height=dimension, width=dimension)
    #limx = i_profiles[[1]]$
    if (nbni > 0) {
        margin = c(2, 2, 0.5, 0.5)
        par(mfcol=c(nbni, nbni), mar=margin, lab=c(3, 3, 4))
        for (ni in 1:nbni) {
            # Set the y margins
            if (ni == 1) {margin[2]=4;}
            else {margin[2]=2;}

            for (j in 1:nbni) {
                # Set the x margins
                if (j == nbni) {margin[1]=4;}
                else {margin[1]=2;}

                par(mar=margin)
                if (j == ni) {
                    limy = c( range(ni_profiles[[ni]]$residuals[ni_profiles[[j]]$pathid,], na.rm=T)[1], ni_profiles[[ni]]$thresholds[2] * 1.1)
                }
                else {
                    limy = range(ni_profiles[[ni]]$residuals[ni_profiles[[j]]$pathid,], na.rm=T)
                }
                #message (limy)
                ### Modification of limy if it reaches Inf, should not have to be done
                if (!is.finite(limy[1]) ) {limy[1] = 0; message(paste("Error in low lim value :", ni_profiles[[ni]]$thresholds[2]))}
                if (!is.finite(limy[2]) ) {limy[2] = 10000; message(paste("Error in high lim value :", ni_profiles[[ni]]$thresholds[2]))}
                
                plot(ni_profiles[[ni]]$explored, ni_profiles[[ni]]$residuals[ni_profiles[[j]]$pathid,], xlab="", ylab="", type="l", col=abs(ni-j)+1, ylim = limy)
                # Print labels on the right and bottom
                if (ni == 1) { title(ylab=ni_profiles[[j]]$path, las=2); }
                if (j == nbni) { title(xlab=ni_profiles[[ni]]$path); }
                if (j == ni) {
                    # Could be accelerated with two points instead of hundreds
                    lines( ni_profiles[[ni]]$explored, rep(ni_profiles[[ni]]$thresholds[1], length(ni_profiles[[ni]]$explored)), lty=2, col="grey" )
                    lines( ni_profiles[[ni]]$explored, rep(ni_profiles[[ni]]$thresholds[2], length(ni_profiles[[ni]]$explored)), lty=2, col="grey" )
                    pl_range = range(ni_profiles[[ni]]$resniuals[ni_profiles[[ni]]$pathni,])
                    lines( rep(ni_profiles[[ni]]$value, 2), c( pl_range[1]-0.1*(pl_range[1]-ni_profiles[[ni]]$threshold[1]), pl_range[1]+0.1*(pl_range[1]-ni_profiles[[ni]]$threshold[1]) ), col="red")
                }
            }
            message(paste("Non identifiable path", ni_profiles[[ni]]$path, "plotted"))
        }
    }

# Identifiables
    nbid = length(i_profiles)
    message(paste(nbid, "identifiable paths"))
    par( mfcol=c(1, 2), mar=c(3, 2, 0, 1), oma=c(0, 0, 2, 0) )
    if (nbid > 0) {
        for (id in 1:nbid) {
            # Plot the profile likelihood with the thresholds and a tick for the optimum
            plot(i_profiles[[id]]$explored, i_profiles[[id]]$residuals[i_profiles[[id]]$pathid, ], type="l", sub=paste(i_profiles[[id]]$path, "profile"))
            lines( i_profiles[[id]]$explored, rep(i_profiles[[id]]$thresholds[1], length(i_profiles[[id]]$explored)), lty=2, col="grey" )
            lines( i_profiles[[id]]$explored, rep(i_profiles[[id]]$thresholds[2], length(i_profiles[[id]]$explored)), lty=2, col="grey" )
            pl_range = range(i_profiles[[id]]$residuals[i_profiles[[id]]$pathid,])
            lines( rep(i_profiles[[id]]$value, 2), c( pl_range[1]-0.1*(pl_range[1]-i_profiles[[id]]$threshold[1]), pl_range[1]+0.1*(pl_range[1]-i_profiles[[id]]$threshold[1]) ), col="red")

            plot(1, type="n", xlim=range(i_profiles[[id]]$explored), ylim=range( i_profiles[[id]]$residuals[-i_profiles[[id]]$pathid,], na.rm=T) )
            # Plot the functionnal relation
            for (i in 1:dim(i_profiles[[id]]$residuals)[1]) {
                if (i != i_profiles[[id]]$pathid) {
                    lines(i_profiles[[id]]$explored, i_profiles[[id]]$residuals[i, ], sub="Functionnal relation", col=i)
                }
            }
            title (main=i_profiles[[id]]$path, outer=T)
            message(paste("Identifiable path", i_profiles[[id]]$path, "plotted") )
        }
    }

    dev.off()

}

# Separates the profiles whether they are identifiables or not
classify_profiles <- function (profiles) {
    ni_profiles = list()
    i_profiles = list()

    for (lprofile in profiles) {
        # Parameters are non identifiable if their profile likelihood does not reach the low threshold on both sides of the minimum 
        if (lprofile$lower_pointwise && lprofile$upper_pointwise) {
            i_profiles[[length(i_profiles)+1]] <- lprofile
        }
        else {
            ni_profiles[[length(ni_profiles)+1]] <- lprofile
        }
    }
    sorted_profiles = list(i_profiles, ni_profiles)
    message("Profiles sorted")
    return(sorted_profiles)
}

