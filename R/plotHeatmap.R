########################### plotHeatmap.R ###########################
# function to generate heatmap with cutoff and colormap scale function

# A color blind palette
#                Green      Orange    Sky Blue    Blue      Vermillion  Purple    Yellow      Black
cbbPalette <- c("#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#000000")

#' Plot a customized heatmap with a symmetrical scale that is not stretched by extreme values
#' @param mat a numeric matrix that should be plotted
#' @param lim a single number indicating the maximal range [-lim lim] for the color map
#' @param fixedRange logical indicating whether the range of lim should at all means be kept (ensures comparability with other heatmaps)
#' @param stripOut numeric between [0 1] to determine the range of colors (excluding outliers from distorting the colormap) 
#' @param main a character string representing the title to be passed
#' @param col colorRampPalette object indicating the colormap
#' @param textCol string denoting the color of inset text
#' @param sig_numbers number of significant parameters to display in the heatmap
#' @return Nothing
#' @author Bertram Klinger \email{bertram.klinger@@charite.de}
plotHeatmap <- function(mat,main = "",lim = Inf,fixedRange = F,stripOut=0.05,col = colorRampPalette(c("deepskyblue","white","red1")),textCol = "gray10", sig_numbers=2){
  # helper functions to generate the breaks. When data contain only one sign: 0...+-limit, otherwise -limit...+limit 
  
  # cutoff and transformation of colour
  m = mat
  if (stripOut >= 0.5) { stop("Cannot strip more than 50% of the data to generate the color scale") } 
  if (!is.numeric(lim)) {stop("lim should be numeric")}
  if (!fixedRange) {
      lowLim = max(-lim,quantile(mat,stripOut, na.rm=T))
      upLim = min(lim,quantile(mat,1-stripOut, na.rm=T))
  } else {
      lowLim = -lim
      upLim = lim
  }
  m[m < lowLim] = lowLim
  m[m > upLim] = upLim
  breaks=define_breaks(m,lim,fixedRange)
  
  # linearized matrix order (1,1) (2,1) (3,1) (4,1)... for text inset
  ref=c(t(as.matrix(mat[nrow(mat):1,]))) 

  # Generate heatmap with textual inset
  p<- levelplot(x = t(as.matrix(m[nrow(m):1,])),
                col.regions = col,
                at = breaks,
                colorkey = list(space = "left"),
                aspect = "fill",
                xlab = "",
                ylab = "",
                main = main,
                scales = list(alternating = 2,tck = c(0,1),x = list(rot = 90)),
                panel=function(...){
                      arg<-list(...)
                      panel.levelplot(...)
                      panel.text(arg$x,
                                 arg$y,
                                 signif(ref,sig_numbers),
                                 col = textCol,
                                 cex = 0.8)})
  print(p) # Plot the heatmap
}

define_breaks <- function(m,lim = Inf,fixedRange = F) {
  if (!fixedRange) {
      limit = min(lim,(max(abs(m),na.rm = T)))
      return(seq(-1.1*(limit)*ifelse(min(m,na.rm=T)<0,1,0),1.1*limit*ifelse(max(m,na.rm=T)>0,1,0),length.out=22))
  } else {
      if (is.infinite(lim) || is.nan(lim) || is.na(lim)) {
          stop("'lim' is invalid, cannot generate breaks within a fixed range")
      }
      return(seq(-1.1*lim,1.1*lim,length.out=22))  
  }
}
