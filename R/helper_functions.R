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
