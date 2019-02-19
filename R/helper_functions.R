# Helper functions for the STASNet package, mostly patch stupid R behaviour

#' Return a data.frame even if the output as on line or one column
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
        return(NA)
    } else if (length(xx) > 0) {
        return( exp(sd(log(xx), na.rm=na.rm)) )
    } else {
        warn("Only negative values for linear_sd_log !")
        return(NA)
    }
}
