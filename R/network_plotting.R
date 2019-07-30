################################ network_plotting.R #########################################

#' Plot a graph from an adjacency list
#'
#' @param structure A 2-columns matrix or a ModelStructure object. The network as an adjacency list, the first column is the upstream nodes, the second column the downstream nodes. Or a ModelStructure object as returned by getModelStructure.
#' @param expdes An ExperimentalDesign object. The measured, stimulated and inhibited nodes are highlighted if present. Signalling strengths are indicated in leftshifted edges and inhibitor strengths are denoted in red below the inhibited node.
#' @param local_values A list with entries 'local_response' (A weighted adjacency matrix representing the values of the links) and 'inhibitors' (A list of inhibition values) both compatible with the 'structure' input
#' @param print_values If local_values has values, whether those values should be printed on top of the arrows in the graph.
#' @param scaling Maximum parameter value
#' @param max_width Maximum line width to use.
#' @param min_width Minimum line width to use.
#' @details The width of the lines is set in ['min_width','max_width'] according to the 'scaling' parameter, with 0 corresponding to 'min_width' and >='scaling' to 'max_width'.
#' @export
#' @family Network graph
#' @author Bertram Klinger \email{bertram.klinger@@charite.de}
plotNetworkGraph <- function(structure, expdes="", local_values="", print_values=TRUE, scaling=5, max_width=3, min_width=1) {
    if (class(structure) == "matrix") {
      ss = extractStructure(structure)
      adm = ss$adjacencyMatrix
      colnames(adm) = rownames(adm) = ss$names
    } else if (class(structure) == "Rcpp_ModelStructure") {
        adm = structure$adjacencyMatrix
        colnames(adm) = rownames(adm) = structure$names
    } else {
        stop("Invalid 'structure' in plotNetworkGraph, must be an edge list or a ModelStructure")
    }
  
  len=length(rownames(adm))
  g1 <- graph::graphAM(adjMat=t(adm),edgemode="directed")
  
  # add inhibitors as pseudo nodes downstream of inhibited nodes in order to depict their strength  
  if (class(expdes) == "Rcpp_ExperimentalDesign" && any(local_values != "")){
    if (length(expdes$inhib_nodes)>0){
      for (nn in rownames(adm)[1+expdes$inhib_nodes]){
        g1 <- graph::addNode(paste0(nn,"i"),g1)
        g1 <- graph::addEdge(nn,paste0(nn,"i"),g1)
      }
    }
  }
  
  # setting of general and creation of changed properties
  graph::nodeRenderInfo(g1) <- list(shape="ellipse")
  graph::nodeRenderInfo(g1) <- list(textCol="black")
  graph::nodeRenderInfo(g1) <- list(lwd=1)
  graph::edgeRenderInfo(g1) <- list(fontsize=10)
  graph::edgeRenderInfo(g1) <- list(textCol="black")
  graph::edgeRenderInfo(g1) <- list(col="black")
  
  g1 <- Rgraphviz::layoutGraph(g1)

  # Add the experimental setup if provided
  if (class(expdes) == "Rcpp_ExperimentalDesign") {
    graph::nodeRenderInfo(g1)$fill[1+expdes$measured_nodes] = "#ffff66"
    
    if (length(expdes$inhib_nodes)>0){
      graph::nodeRenderInfo(g1)$lwd[1+expdes$inhib_nodes]=4 # Populate for perturbations
      graph::nodeRenderInfo(g1)$col[1+expdes$inhib_nodes] = "red"
      graph::nodeRenderInfo(g1)$col[(len+1):(len+length(expdes$inhib_nodes))]="white" # mask inhibitor pseudo nodes
      graph::nodeRenderInfo(g1)$textCol[(len+1):(len+length(expdes$inhib_nodes))]="white" 
    }
    
    if (length(expdes$stim_nodes)>0){
      graph::nodeRenderInfo(g1)$lwd[1+expdes$stim_nodes] = 4
      graph::nodeRenderInfo(g1)$col[1+expdes$stim_nodes] = "blue"
    }
  }
  if (local_values[1] != "") {
    # Add Edge Weights left justified
    efrom = graph::edgeRenderInfo(g1)$enamesFrom
    eto = graph::edgeRenderInfo(g1)$enamesTo
    edge_spline = graph::edgeRenderInfo(g1)$splines
    
    for (idx in which(adm!=0)) {
      vv = local_values$local_response[idx]
      afrom = colnames(adm)[ceiling(idx/len)]
      ato = rownames(adm)[ifelse(idx %% len==0,len,idx %% len)]
      cc = which(afrom==efrom & ato==eto)
      cc = ifelse(length(cc)!=0, cc, which(afrom==eto & ato==efrom)) # Link in both directions
#      graph::edgeRenderInfo(g1)$lwd[cc] = ifelse(abs(vv)<=1,1,ifelse(abs(vv)<=scaling,2,max_width))
      graph::edgeRenderInfo(g1)$lwd[cc] = ifelse( abs(vv)>scaling, max_width, min_width + (max_width-min_width) * (abs(vv)/scaling) )
      ifelse(abs(vv)<=1,1,ifelse(abs(vv)<=scaling,2,max_width))
      if (!is.na(vv)){
          if (vv < 0) { graph::edgeRenderInfo(g1)$col[cc] = "orange" }
      }
      if (print_values) {
          graph::edgeRenderInfo(g1)$label[cc] = trim_num(vv)
          
          coordMat=Rgraphviz::bezierPoints(edge_spline[[cc]][[1]]) # 11 x 2 matrix with x and y coordinates
          graph::edgeRenderInfo(g1)$labelX[cc] = coordMat[5,"x"]-ceiling(nchar(graph::edgeRenderInfo(g1)$label[cc])*10/2)
          graph::edgeRenderInfo(g1)$labelY[cc] = coordMat[5,"y"]
      }
    }
    
    # Add Inhibitor estimates
    if (length(expdes$inhib_nodes)>0){
      for (idx in 1:length(expdes$inhib_nodes)) {
        vv = local_values$inhibitors[idx]
        iname = paste0(colnames(adm)[expdes$inhib_nodes[idx]+1], "i")
        nname = colnames(adm)[expdes$inhib_nodes[idx]+1]
        cc = which(nname==efrom & iname==eto)
        graph::edgeRenderInfo(g1)$col[cc]="white" # mask inhibitor pseudo edges
        if (print_values) {
            graph::edgeRenderInfo(g1)$label[cc] = trim_num(vv)
            graph::edgeRenderInfo(g1)$textCol[cc]="red"

            coordMat = Rgraphviz::bezierPoints(edge_spline[[cc]][[1]]) # 11 x 2 matrix with x and y coordinates 
            graph::edgeRenderInfo(g1)$labelX[cc] = coordMat[2,"x"]
            graph::edgeRenderInfo(g1)$labelY[cc] = coordMat[2,"y"]
        }
        
      }
    }
  }
  
  Rgraphviz::renderGraph(g1)
# other options to use the width and height howver the coordinates are not thes same in the generated graph!!  
#  if (length(expdes$inhib_nodes)>0){
#      nodes = names(graph::nodeRenderInfo(g1)$nodeX)
#      inhib = colnames(adm)[expdes$inhib_nodes+1]      
#      cc = which(nodes %in% inhib)
#      ix = graph::nodeRenderInfo(g1)$labelX[cc] + graph::nodeRenderInfo(g1)$lWidth[cc]
#      iy = graph::nodeRenderInfo(g1)$labelY[cc] - 0.5*graph::nodeRenderInfo(g1)$height[cc]
#      iv = local_values$inhibitors  
#      if (local_values[1] != "" & length(expdes$inhib_nodes)>0){
#       text(x = ix, y = iy, labels = trim_num(iv), col="red",cex=0.6, pos=4, offset=0.5)
#      }
#   }
  invisible(g1)
  # (1) TODO MARK REMOVED LINKS, (2) ALLOW TO GIVE CLUSTERS THAT SHOULD BE KEPT IN CLOSE VICINITY 
}
