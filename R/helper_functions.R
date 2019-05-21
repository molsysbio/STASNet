# Helper functions for the STASNet package, mostly patch stupid R behaviour

#' Return a data.frame even if the output is one-dimensional
sub_data_frame <- function(input, rows=NULL, cols=NULL) {
  output = input
  if (!is.null(rows)) {
    ocnames = colnames(output)
    if (is.character(rows)) {
      rows = rownames(output) %in% rows
    }
    rnames = rownames(output)[rows]
    output = output[rows,]
    if (is.null(nrow(output))) {
      output = as.data.frame(output, row.names=rnames)
      rownames(output) = rnames
      colnames(output) = ocnames
    }
  }

  if (!is.null(cols)) {
    ornames = rownames(output)
    if (is.character(cols)) {
      cols = colnames(output) %in% cols
    }
    cnames = colnames(output)[cols]
    output = output[,cols]
    if (is.null(ncol(output))) {
      output = as.data.frame(output)
      colnames(output) = cnames
      rownames(output) = ornames
    }
  }

  return(output)
}

#' Geometric mean
#'
#' Compute geometric mean with na.rm == TRUE
geom_mean <- function(xx, na.rm=TRUE) {
    xx=xx[xx>0]
    if (all(is.na(xx))) {
        return(NaN) # Return NaN for consistency with mean
    } else if (length(xx)>0) {
        xx = xx[!is.na(xx)]
        return( exp(sum(log(xx), na.rm=na.rm)/length(xx)) )
    } else {
        warning("No or only negative values for geom_mean !")
        return(NA)
    }
}

#' Log-normal standard deviation
#'
#' Compute the standard deviation of a set following a log-normal distribution
log_norm_sd <- function(xx, na.rm=TRUE) {
    xx = xx[xx>0]
    if (length(xx) > 0) {
        return( exp(sd(log(xx))^2-1) * exp(2*mean(log(xx))+sd(log(xx))^2) )
    } else {
        warn("Negative values for log_norm_sd !")
        return(NA)
    }
}

#' Standard deviation of the log
#'
#' Standard deviation of the log in linear space
linear_sd_log <- function(xx, na.rm=TRUE) {
    xx = xx[xx>0]
    if (all(is.na(xx))) {
        return(NaN)
    } else if (length(xx) > 0) {
        return( exp(sd(log(xx), na.rm=na.rm)) )
    } else {
        warn("Only negative values for linear_sd_log !")
        return(NA)
    }
}

#' Plot residuals vs rank
#'
#' Plot all residuals vs rank and best residuals vs rank to assess the quality of the optimisation
#' @param A list of residuals
residuals_plot <- function(residuals, model_name="default") {
    order_resid = order(residuals,na.last = T)
    order_id = order_resid[1:min(20, length(residuals))]
    if ( all(is.na(residuals)) ) {
        warning("All residuals are NAs! (no residual plot possible)")
    } else {
        plot(1:length(order_resid), residuals[order_resid], main=paste0("Best residuals ", model_name), ylab="Likelihood", xlab="rank", log="y",type="l",lwd=2)
        lines(1:length(order_id), residuals[order_id], col="red")
        if (length(order_resid) >= 100) {
            hundred_best = residuals[order_resid[1:100]]
            plot(1:100, hundred_best, main=paste0("Best 100 residuals ", model_name), ylab="Likelihood", xlab="rank", log="y",type="l",lwd=2, ylim=c(hundred_best[1], hundred_best[100]+1))
            lines(1:length(order_id), residuals[order_id], col="red")
        }
    }
}
