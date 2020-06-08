################################ network_plotting.R #########################################

#' Plot a graph from an adjacency list
#'
#' @param structure A 2-columns matrix or a ModelStructure object. The network as an adjacency list, the first column is the upstream nodes, the second column the downstream nodes. Or a ModelStructure object as returned by extractStructure.
#' @param expdes An ExperimentalDesign object. The measured, stimulated and inhibited nodes are highlighted if present. Signalling strengths are indicated in leftshifted edges and inhibitor strengths are denoted in red below the inhibited node.
#' @param local_values A list with entries 'local_response' (A weighted adjacency matrix representing the values of the links) and 'inhibitors' (A list of inhibition values) both compatible with the 'structure' input
#' @param print_values If local_values has values, whether those values should be printed on top of the arrows in the graph.
#' @param scaling Maximum parameter value to take into account to scale the edge width.
#' @param max_width Maximum line width to use.
#' @param min_width Minimum line width to use.
#' @param main Title to be written above the graph
#' @param sub Title to be written below the graph
#spline: man uebergibt ihm eine Egde wo der Knoten den man verschieben will in unserem Beispiel (MEK~ERK) er muss an der Stelle nach dem ~ stehen da die Funktion so aufgebaut ist das er den hinteren Wert als Knoten nimmt den er neu erstellen soll
#' @details The width of the lines is set in ['min_width','max_width'] according to the 'scaling' parameter, with 0 corresponding to 'min_width' and >='scaling' to 'max_width'.
#' @export
#' @family Network graph
#' @author Bertram Klinger \email{bertram.klinger@@charite.de}
plotNetworkGraph <- function(structure, expdes="", local_values="", print_values=FALSE, scaling=5, max_width=3, min_width=1, main="", sub="", spline="") {
  if (class(structure) == "matrix") {
    nstructure = extractStructure(structure)
    adm = nstructure$adjacencyMatrix
    colnames(adm) = rownames(adm) = nstructure$names
    names=nstructure$names
  } else if (class(structure) == "Rcpp_ModelStructure") {
    adm = structure$adjacencyMatrix
    colnames(adm) = rownames(adm) = structure$names
    names=structure$names
  } else {
    stop("Invalid 'structure' in plotNetworkGraph, must be an edge list or a ModelStructure")
  }

  # Generate the graph
  g2=graph::graphAM(adjMat = t(adm), edgemode = "directed")

  # The layout depends on the order of the edges
  # For those where we want to force which node is downstream, the upstream links have to be added to the graph before the downstream links
  # We duplicate the node and hope that the second one will end up at a better position
  if (spline[1]!="") {
    g1=layoutGraph(g2)
    for (ii in seq_along(spline)) {
      eTo=edgeRenderInfo(g1)$enamesTo[spline[ii]] # Target node
      eBis = paste0(" ", eTo, " ") # New node
      eFrom=which(eTo==edgeRenderInfo(g1)$enamesTo) # Edge to target
      nTo=edges(g2)[eTo][[1]] # Downstream of target node (GSK3, p90RSK, RAF)
      # Create a secondary downstream node and link it first to upstream, then to downstream nodes
      g2=addNode(eBis, g2)
      for (uid in seq_along(eFrom)) {
        g2=addEdge(edgeRenderInfo(g1)$enamesFrom[names(eFrom[uid])], eBis, g2) # Add all edges upstream->duplicated
      }
      for (did in seq_along(nTo)) {
        g2=addEdge(eBis, nTo[did], g2) # Add all edges duplicated->downstream
      }
      # Copy values in adjacency matrix
      tmp=which(eTo==colnames(adm))
      colnames(adm)[tmp]=rownames(adm)[tmp]=eBis
    }
  }
  #bidir_edges is a vector, in it are all edges their are in both directions
  bidir_edges=c()
  # for the dummy-nodes in the graph (for the inhibitors)
  dummy=c()
  # First LayoutGraph needed to know which edges are in both direction
  g1=layoutGraph(g2)
  for (ii in seq_along(edgeRenderInfo(g1)$direction)) {
    if (edgeRenderInfo(g1)$direction[[ii]]=="both") {
      bidir_edge_name=names(edgeRenderInfo(g1)$direction[ii])
      #bidir_edges is a vector, in it are all edges their are in both directions
      bidir_edges[length(bidir_edges)+1]=bidir_edge_name
      bidir_points=which(bidir_edge_name==names(edgeRenderInfo(g1)$enamesTo))
      #c is a variable for the from varaible for the new edge
      bidir_to=edgeRenderInfo(g1)$enamesTo[[bidir_points]]
      dummy[length(dummy)+1]=as.character(ii)
      g2=addNode(as.character(ii), g2)
      g2=addEdge(as.character(ii), bidir_to, g2)
      
      #updated the adm
      if (local_values[1]!="") {
        
        tmp=which(bidir_edge_name==names(edgeRenderInfo(g1)$enamesFrom))
        bidir_from=edgeRenderInfo(g1)$enamesFrom[[tmp]]
        adm=cbind(adm, c(0))
        colnames(adm)=c(rownames(adm), dummy[length(dummy)])
        adm=rbind(adm, c(0))
        rownames(adm)=colnames(adm)
        adm[bidir_to,dummy[length(dummy)]]=adm[bidir_to,bidir_from]
        
        for (aa in seq_along(local_values)) {
          local_values[[aa]]$local_response=cbind(local_values[[aa]]$local_response, c(0))
          local_values[[aa]]$local_response=rbind(local_values[[aa]]$local_response, c(0))
          oldto=which(as.character(bidir_to)==colnames(adm))
          oldfrom=which(as.character(bidir_from)==colnames(adm))
          newto=which(as.character(dummy[length(dummy)])==colnames(adm))
          local_values[[aa]]$local_response[oldto, newto]=local_values[[aa]]$local_response[oldto, oldfrom]
        }
        
      }
    }    
  }
  #inhib_nodes is the vector with all inhib edges
  inhib_nodes=c()
  #the edge for the inhib number of inhib nodes
  if (class(expdes) == "Rcpp_ExperimentalDesign"& local_values[1]!="") {
    for (ii in seq_along(expdes[["inhib_nodes"]])) {
      inhib_nodes[length(inhib_nodes)+1]=names[expdes[["inhib_nodes"]][ii]+1]
      g2=graph::addEdge(inhib_nodes[ii], inhib_nodes[ii], g2)
    }
  }
  #the second layoutGraph is with the new edges and the new nodes
  g1=layoutGraph(g2)
  # renderGraph(g1)
  nodeInfo=nodeRenderInfo(g1)$nodeX
  
  graph::nodeRenderInfo(g1) <- list(shape="ellipse")
  graph::nodeRenderInfo(g1) <- list(textCol="black")
  graph::nodeRenderInfo(g1) <- list(lwd=1)
  graph::edgeRenderInfo(g1) <- list(fontsize=10)
  graph::edgeRenderInfo(g1) <- list(textCol="black")
  graph::edgeRenderInfo(g1) <- list(col="black")
  
  
  bidir_id=0
  #the all edges on the od position of g1
  for (ss in seq_along(edgeRenderInfo(g1)$splines)) {
    # For all edges with both direction, make a new edge for the reverse direction
    if (edgeRenderInfo(g1)$direction[ss]=="both") {
      bidir_id=bidir_id+1
      edgeRenderInfo(g1)$splines[[bidir_edges[bidir_id]]]=edgeRenderInfo(g1)$splines[[bidir_edges[bidir_id]]]
      
      #there are new nodes and new edges then, we move the new egde 
      #gg is a variable so we know the new node in the graph, it is for later to make them transparent
      gg=which(dummy[bidir_id]==edgeRenderInfo(g1)$enamesFrom)
      gg=names(edgeRenderInfo(g1)$enamesFrom[gg])
      edge_spline = graph::edgeRenderInfo(g1)$splines
      
      bidir_edge_name=bidir_edges[bidir_id]
      for (jj in length(edgeRenderInfo(g1)$splines[[bidir_edge_name]]):1) {
        bidir_spline = edgeRenderInfo(g1)$splines[[bidir_edge_name]][[jj]]
        bidir_points=cPoints(bidir_spline)
        if (as.integer(bezierPoints(bidir_spline)[1,"x"])==as.integer(bezierPoints(bidir_spline)[2,"x"])) {
          #newcp is a List for the new cPoints for the move the old edge(in both direction) in a right position
          newcP=list()
          #create the new cPoints old edge
          for (ii in seq_along(bidir_points)) {
            newcP[[ii]]=new("xyPoint", x=as.integer(getX(bidir_points[[ii]])+10), y=as.integer(getY(bidir_points[[ii]])))
          }
          newBC=new("BezierCurve", cPoints=newcP)
          edge_spline[[bidir_edge_name]][[jj]]=newBC
          edgeRenderInfo(g1)$splines=edge_spline
          #newcp1 is a List for the new cPoints for the move the new edge in a right position
          newcP1=list()
          #new cPoints for the new egde
          for (ii in seq_along(bidir_points)) {
            newcP1[[ii]]=new("xyPoint", x=as.integer(getX(bidir_points[[ii]])-10), y=as.integer(getY(bidir_points[[ii]])))
          }
          newBC1=new("BezierCurve", cPoints=newcP1)
          edge_spline[[gg]][[jj]]=newBC1 
          edgeRenderInfo(g1)$splines=edge_spline
        }
        else if (as.integer(bezierPoints(bidir_spline)[1,"y"])==as.integer(bezierPoints(bidir_spline)[2,"y"])) {
          #newcp is a List for the new cPoints for the move the old edge(in both direction) in a right position
          newcP=list()
          #create the new cPoints old edge
          for (ii in seq_along(bidir_points)) {
            newcP[[ii]]=new("xyPoint", x=as.integer(getY(bidir_points[[ii]])+10), y=as.integer(getX(bidir_points[[ii]])))
          }
          newBC=new("BezierCurve", cPoints=newcP)
          edge_spline[[bidir_edge_name]][[jj]]=newBC
          edgeRenderInfo(g1)$splines=edge_spline
          #newcp1 is a List for the new cPoints for the move the new edge in a right position
          newcP1=list()
          #new cPoints for the new egde
          for (ii in seq_along(bidir_points)) {
            newcP1[[ii]]=new("xyPoint", x=as.integer(getY(bidir_points[[ii]])-10), y=as.integer(getX(bidir_points[[ii]])))
          }
          newBC1=new("BezierCurve", cPoints=newcP1)
          edge_spline[[gg]][[jj]]=newBC1
          edgeRenderInfo(g1)$splines=edge_spline
        }
        else {
          if ((bezierPoints(bidir_spline)[1,"y"]-bezierPoints(bidir_spline)[2,"y"]>0&&bezierPoints(bidir_spline)[1,"x"]-bezierPoints(bidir_spline)[2,"x"]>0)||(bezierPoints(bidir_spline)[1,"y"]-bezierPoints(bidir_spline)[2,"y"]<0&&bezierPoints(bidir_spline)[1,"x"]-bezierPoints(bidir_spline)[2,"x"]>0)) {
            #newcp is a List for the new cPoints for the move the old edge(in both direction) in a right position
            newcP=list()
            #create the new cPoints old edge
            for (ii in seq_along(bidir_points)) {
              newcP[[ii]]=new("xyPoint", x=as.integer(getX(bidir_points[[ii]])+7 ), y=as.integer(getY(bidir_points[[ii]])-5))
            }
            newBC=new("BezierCurve", cPoints=newcP)
            edge_spline[[bidir_edge_name]][[jj]]=newBC
            edgeRenderInfo(g1)$splines=edge_spline
            #newcp1 is a List for the new cPoints for the move the new edge in a right position
            newcP1=list()
            #new cPoints for the new egde
            for (ii in seq_along(bidir_points)) {
              newcP1[[ii]]=new("xyPoint", x=as.integer(getX(bidir_points[[ii]])-7), y=as.integer(getY(bidir_points[[ii]])+5))
            }
            newBC1=new("BezierCurve", cPoints=newcP1)
            edge_spline[[gg]][[jj]]=newBC1
            edgeRenderInfo(g1)$splines=edge_spline
          }
          else if ((bezierPoints(bidir_spline)[2,"y"]-bezierPoints(bidir_spline)[2,"y"]<0&&bezierPoints(bidir_spline)[2,"x"]-bezierPoints(bidir_spline)[2,"x"]<0 )||(bezierPoints(bidir_spline)[1,"y"]-bezierPoints(bidir_spline)[2,"y"]>0&&bezierPoints(bidir_spline)[1,"x"]-bezierPoints(bidir_spline)[2,"x"]<0)) {
            #newcp is a List for the new cPoints for the move the old edge(in both direction) in a right position
            newcP=list()
            #create the new cPoints old edge
            for (ii in seq_along(bidir_points)) {
              newcP[[ii]]=new("xyPoint", x=as.integer(getX(bidir_points[[ii]])+5), y=as.integer(getY(bidir_points[[ii]])+5))
            }
            newBC=new("BezierCurve", cPoints=newcP)
            edge_spline[[bidir_edge_name]][[jj]]=newBC
            edgeRenderInfo(g1)$splines=edge_spline
            #newcp1 is a List for the new cPoints for the move the new edge in a right position
            newcP1=list()
            #new cPoints for the new egde
            for (ii in seq_along(bidir_points)) {
              newcP1[[ii]]=new("xyPoint", x=as.integer(getX(bidir_points[[ii]])-5), y=as.integer(getY(bidir_points[[ii]])-5))
            }
            newBC1=new("BezierCurve", cPoints=newcP1)
            edge_spline[[gg]][[jj]]=newBC1
            edgeRenderInfo(g1)$splines=edge_spline
          }
        }
        #one arrow of the old edge disappear
        edgeRenderInfo(g1)$arrowhead[bidir_edges]="none"
        # Make the dummy-nodes invisible
        nodeRenderInfo(g1)$col[dummy]="transparent"
        nodeRenderInfo(g1)$textCol[dummy]="transparent"
      }
    } else {
      name=names(edgeRenderInfo(g1)$splines[ss])
      edgeRenderInfo(g1)$splines[[name]]=edgeRenderInfo(g1)$splines[[name]]
    }
  }
  # Colours and line width for the inhib, measured and stim nodes
  if (class(expdes) == "Rcpp_ExperimentalDesign") {
    graph::nodeRenderInfo(g1)$col[expdes[["inhib_nodes"]]+1]="red"
    graph::nodeRenderInfo(g1)$lwd[expdes[["inhib_nodes"]]+1]=4
    graph::nodeRenderInfo(g1)$col[expdes[["stim_nodes"]]+1]="blue"
    graph::nodeRenderInfo(g1)$lwd[expdes[["stim_nodes"]]+1]=4
    graph::nodeRenderInfo(g1)$fill[expdes[["measured_nodes"]]+1]="#ffff66"
  }
  # Write the weight of the edge
  if (local_values[1]!="") {
    efrom = graph::edgeRenderInfo(g1)$enamesFrom
    eto = graph::edgeRenderInfo(g1)$enamesTo
    edge_spline = graph::edgeRenderInfo(g1)$splines
    graph::edgeRenderInfo(g1)$label=NA
    for (ii in which(adm!=0)) {
      edge_values=c()
      rname=rownames(adm)[ifelse(ii%%nrow(adm)==0 , nrow(adm), ii%%nrow(adm))]
      cname=colnames(adm)[ceiling(ii/ncol(adm))]
      cc = which(cname==efrom & rname==eto)
      cc = ifelse(length(cc)!=0, cc, which(cname==eto & rname==efrom))
      neg_values=0 # counts the number of models where the link is negative, the color is orange if more than half are
      for (mid in seq_along(local_values)) {
        evalue=as.numeric( STASNet:::trim_num(local_values[[mid]]$local_response[ii]) )
        if (mid%%4 == 0) {
          edge_values=paste(edge_values, '
                  ')
        }
        if (mid==1 || mid%%4 == 0) {
          edge_values=paste0(edge_values, evalue)
        } else{
          edge_values=paste0(edge_values, "|", evalue)
        }
        if (print_values) {
          graph::edgeRenderInfo(g1)$label[cc]=edge_values
        } else {
          graph::edgeRenderInfo(g1)$label[cc]=""
        }
        
        if (evalue<0) {
          neg_values=neg_values+1
        }
      }
      coordMat=Rgraphviz::bezierPoints(edge_spline[[cc]][[1]])# 11 x 2 matrix with x and y coordinates
      if (length(local_values)==1) {
        graph::edgeRenderInfo(g1)$lwd[cc] = ifelse( abs(evalue)>scaling, max_width, min_width + (max_width-min_width) * (abs(evalue)/scaling) )
        graph::edgeRenderInfo(g1)$labelX[cc] = coordMat[5,"x"]-ceiling(nchar(graph::edgeRenderInfo(g1)$label[cc])*10/2)
        if (edgeRenderInfo(g1)$direction[[cc]]=="both") {
          graph::edgeRenderInfo(g1)$labelX[cc] = coordMat[5,"x"]+ceiling(nchar(graph::edgeRenderInfo(g1)$label[cc])*10/2)
        }
      } else{
        graph::edgeRenderInfo(g1)$labelX[cc] = coordMat[5,"x"]-ceiling((nchar(graph::edgeRenderInfo(g1)$label[cc])/2)*10/2)-3
        if (edgeRenderInfo(g1)$direction[[cc]]=="both") {
          graph::edgeRenderInfo(g1)$labelX[cc] = coordMat[5,"x"]+ceiling((nchar(graph::edgeRenderInfo(g1)$label[cc])/2)*10/2)+3
        }
      }
      graph::edgeRenderInfo(g1)$labelY[cc]= coordMat[5,"y"]
      if (neg_values > length(local_values)/2) {
        edgeRenderInfo(g1)$col[cc]="orange"
      }
    }
  }
    
  if (spline[1]!="") {
    # Delete the old node downstream of the recalculated spline
    for (gs in seq_along(spline)) {
      eTo=edgeRenderInfo(g1)$enamesTo[spline[gs]]
      eBis = paste0(" ", eTo, " ")
      nTo=edges(g1)[eTo]
      eFrom=which(eTo==edgeRenderInfo(g1)$enamesTo)
      nodeRenderInfo(g1)$col[eBis]=nodeRenderInfo(g1)$col[eTo]
      nodeRenderInfo(g1)$fill[eBis]=nodeRenderInfo(g1)$fill[eTo]
      for (did in seq_along(nTo)) {
        g1=removeEdge(eTo, nTo[[1]][did], g1)
      }
      for (uid in seq_along(eFrom)) {
        g1=removeEdge(edgeRenderInfo(g1)$enamesFrom[names(eFrom[uid])], eTo, g1)
      }
      g1=removeNode(eTo, g1)
      #nodeRenderInfo(g1)$label[eBis]=eTo
    }
  }
  
  #the number of inhib nodes and for the rigth position of the edges for this numbers
  #
  dummy=which(is.na(edgeRenderInfo(g1)$label))
  inh_id=0
  for (dummy_names in names(dummy)) {
    to_rm_edges=c()
    if (spline[1]!="") {
      for (g in seq_along(spline)) {
        eTo = edgeRenderInfo(g1)$enamesTo[spline[g]]
        to_rm_edges = which(eTo==edgeRenderInfo(g1)$enamesFrom[dummy_names] | eTo==edgeRenderInfo(g1)$enamesTo[dummy_names])
        to_rm_edges[g]=ifelse(length(to_rm_edges)==0, 0, 1)
      }
      to_rm_edges = which(1==to_rm_edges)
    }
    #if for the edges for example MEK~MEK
    if (length(to_rm_edges)==0) {
      inh_id = inh_id + 1
      #number is a variable for the difference of the position from old to new
      number=which(inhib_nodes[inh_id]==names(nodeRenderInfo(g1)$nodeX))
      number=nodeRenderInfo(g1)$nodeX[number]-nodeInfo[number]
      #is for the move of the edge in the rigth position
      for (jj in length(edgeRenderInfo(g1)$splines[[dummy_names]]):1) {
        bidir_points = cPoints( edgeRenderInfo(g1)$splines[[dummy_names]][[jj]] )
        #newcp is a List for the new cPoints for the move the old edge(in both direction) in a right position
        newcP=list()
        #create the new cPoints old edge
        for (kk in seq_along(bidir_points)) {
          newcP[[kk]]=new("xyPoint", x=as.integer(getX(bidir_points[[kk]])+number), y=as.integer(getY(bidir_points[[kk]])))
          newBC=new("BezierCurve", cPoints=newcP)
        }
        newBC=new("BezierCurve", cPoints=newcP)
        edge_spline[[dummy_names]][[jj]]=newBC
        edgeRenderInfo(g1)$splines=edge_spline
      }
      graph::edgeRenderInfo(g1)$textCol[dummy_names]="red"
      graph::edgeRenderInfo(g1)$lwd[dummy_names]=0
      if (print_values) {
        # Label edge with the parameters for all models
        edgeRenderInfo(g1)$label[dummy_names] = paste0(sapply(seq_along(local_values), function(mid) { STASNet:::trim_num(local_values[[mid]]$inhibitors[inh_id]) }), collapse="|")
      }
      coordMat=Rgraphviz::bezierPoints(edge_spline[[dummy_names]][[1]]) # 11 x 2 matrix with x and y coordinates
      graph::edgeRenderInfo(g1)$labelX[dummy_names] = coordMat[5,"x"]
      graph::edgeRenderInfo(g1)$labelY[dummy_names]= coordMat[5,"y"]
      edgeRenderInfo(g1)$col[dummy_names]="transparent"
    }
  }
  graphRenderInfo(g1) <- list(main=main, sub=sub)
  renderGraph(g1)
  invisible(g1)
}
    
