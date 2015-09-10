# Create and manage multiple dataset meant to be applied on the same network

#' Import multiple files 
#' @param data_list Name of the file containing the list of the files for the data without the extension
createDataSet <- function(model_links, data_list, basal_file, cores=1, inits=1000, init_distribution=F, method="default") {
    files = unique(gsub("\\..*$", "", readLines(data_list)))
    folder_files = dir()
    model_set = lsit()
    model_set[[length(files)]] = ""
    for (file in files) {
        if (paste0(file, ".var") %in% folder_files) {
            model_set[[i]] = create_model.R(model_links, paste0(file, ".csv"), basal_file, paste0(file, ".var"), cores, inits, init_distribution, method)
        } else {
            model_set[[i]] = create_model.R(model_links, paste0(file, ".csv"), basal_file, cores, inits, init_distribution, method)
        }
    }
}

#' Compare networks
#'
#' Compare the parameters of one network for different conditions (cell line, perturbations, ...)
#' The links must be the same for all the networks
#' @param files A list of .mra files
#' @return None
# TODO improve the function so that different networks with common links can be compared
compareModels <- function(files) {
  models = list()
  for (i in 1:length(files)) { models[[i]] = importModel(files[i])

  links = c()
  for (model in models) { links = cbind(links, model$parameters) }
  rownames(links) = models[[1]]$model$getParametersLinks()
  colnames(links) = files  
  med = median(abs(links))
  m = max(abs(links))
  breaks = unique(c(seq(-m, -2*med, length.out = 10), seq(-2*med, 2*med, length.out=50), seq(2*med, m, length.out=10)))
  pheatmap(links, breaks = breaks, color=colorRampPalette("deepskyblue", "black", "red")(length(breaks)-1))
  }
}

