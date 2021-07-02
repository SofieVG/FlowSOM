#' PlotFlowSOM 
#' 
#' Base layer to plot a FlowSOM result
#' 
#' Base layer of the FlowSOM plot, where you can choose layout (MST, grid or 
#' coordinates of your own choosing), background colors and node size. Can
#' then be extended by e.g. \code{\link{AddStars}}, \code{\link{AddLabels}},
#' \code{\link{AddPies}}, ...
#' 
#' @param fsom              FlowSOM object, as created by \code{\link{FlowSOM}}
#' @param view              Preferred view, options: "MST", "grid" or "matrix" 
#'                          with a matrix/dataframe consisting of coordinates 
#'                          given in coords. Default = "MST"  
#' @param nodeSizes         A vector containing nodesizes. These will 
#'                          automatically be scaled between 0 and maxNodeSize 
#'                          and transformed with a sqrt. Default = fsom$MST$sizes
#' @param maxNodeSize       Determines the maximum nodesize. Default is 1.
#' @param refNodeSize       Reference for nodesize against which the nodeSizes 
#'                          will be scaled. Default = max(nodeSizes)
#' @param equalNodeSize     If \code{TRUE}, the nodes will be equal to 
#'                          maxNodeSize. If \code{FALSE} (default), the nodes  
#'                          will be scaled to the number of cells in each 
#'                          cluster
#' @param backgroundValues  Values to be used for background coloring, either
#'                          numerical values or something that can be made into
#'                          a factor (e.g. a clustering)
#' @param backgroundColors  Colorpalette to be used for the background coloring.
#'                          Can be either a function or an array specifying
#'                          colors.
#' @param backgroundLim     Only used when backgroundValues are numerical. 
#'                          Defaults to min and max of the backgroundValues.
#' @param title             Title of the plot
#'                          
#' @return A ggplot object with the base layer of a FlowSOM plot
#' 
#' @seealso \code{\link{PlotStars}}, \code{\link{PlotVariable}},
#' \code{\link{PlotMarker}}, \code{\link{PlotLabels}}, 
#' \code{\link{PlotNumbers}}, \code{\link{PlotPies}}, 
#' \code{\link{QueryStarPlot}}, \code{\link{PlotSD}}
#' 
#' @examples 
#' # Locate file on file system
#' fcs_file <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#'
#' # Build FlowSOM model
#' flowSOM.res <- FlowSOM(fcs_file, 
#'                        scale = TRUE,
#'                        compensate = TRUE, 
#'                        transform = TRUE,
#'                        toTransform = 8:18, 
#'                        colsToUse = c(9, 12, 14:18),
#'                        nClus = 10,
#'                        seed = 1)
#'                        
#' # Plot with background coloring
#' PlotFlowSOM(flowSOM.res,
#'             backgroundValues = flowSOM.res$metaclustering)
#'
#' # Own layout
#' mfis <- GetClusterMFIs(flowSOM.res)[,GetChannels(flowSOM.res, c("CD3", "CD4"))]
#' PlotFlowSOM(flowSOM.res,
#'             view = mfis,
#'             maxNodeSize = 0.1,
#'             backgroundValues = flowSOM.res$metaclustering)
#'             
#' # Adapted node sizes
#' PlotFlowSOM(flowSOM.res,
#'               nodeSizes = 1:100,
#'               view = "grid")
#' 
#' @importFrom ggplot2 ggplot coord_fixed theme_void ggtitle
#' 
#' @export
PlotFlowSOM <- function(fsom,
                        view = "MST",
                        nodeSizes = fsom$map$pctgs,
                        maxNodeSize = 1,
                        refNodeSize = max(nodeSizes),
                        equalNodeSize = FALSE,
                        backgroundValues = NULL,
                        backgroundColors = NULL,
                        backgroundLim = NULL,
                        title = NULL)
{
  #----Initialization----
  requireNamespace("ggplot2")
  fsom <- UpdateFlowSOM(fsom)
  
  nNodes <- NClusters(fsom)
  isEmpty <- fsom$map$pctgs == 0
  
  #----Warnings----
  if (length(nodeSizes) != nNodes){
    stop(paste0("Length of 'nodeSizes' should be equal to number of clusters ", 
                "in FlowSOM object (", nNodes, " clusters and ",
                length(nodeSizes), " nodesizes)."))
  }
  
  if (length(backgroundValues) != nNodes && !is.null(backgroundValues)){
    stop(paste0("Length of 'backgroundValues' should be equal to number of ",
                "clusters in FlowSOM object (", nNodes, " clusters and ",
                length(backgroundValues), " backgroundValues)."))
  }
  
  #---- Layout----
  layout <- ParseLayout(fsom, view)
  if(is.matrix(view) || is.data.frame(view)) view <- "matrix"
  
  #---- Nodesize----
  autoNodeSize <- AutoMaxNodeSize(layout = layout, 
                                 overlap = ifelse(view %in% c("grid"),
                                                  -0.3, 1)) 
  maxNodeSize <- autoNodeSize * maxNodeSize
  
  if (equalNodeSize){
    scaledNodeSize <- rep(maxNodeSize, nNodes)
  } else {
    scaledNodeSize <- ParseNodeSize(nodeSizes, maxNodeSize, refNodeSize)
  }
  
  if (any(isEmpty)) {scaledNodeSize[isEmpty] <- min(maxNodeSize, 0.05)}
  
  #----GGplot----
  plot_df <- data.frame(x = layout$x,
                        y = layout$y,
                        size = scaledNodeSize,
                        bg_size = scaledNodeSize * 1.5)
  
  p <- ggplot2::ggplot(plot_df) 
  
  # Add background
  if (!is.null(backgroundValues)){
    p <- AddBackground(p,
                       backgroundValues = backgroundValues, 
                       backgroundColors = backgroundColors,
                       backgroundLim = backgroundLim)
  }
  
  # Add MST
  if (view == "MST"){
    p <- AddMST(p, fsom)
  }
  
  # Add nodes
  p <- AddNodes(p = p,
                values = as.character(isEmpty),
                colorPalette = c("TRUE" = "grey", "FALSE" = "white"),
                showLegend = FALSE)
  
  # Fix plotlayout
  p <- p + ggplot2::coord_fixed() + ggplot2::theme_void() 
  
  if (!is.null(title)){
    p <- p + ggplot2::ggtitle(title)
  }
  return(p)
}

#' FlowSOM default colors
#' 
#' @param n Number of colors to generate
#' 
#' @return array of n colors
#' 
#' @importFrom grDevices colorRampPalette
#' @export
FlowSOM_colors <- function(n){
  grDevices::colorRampPalette(c("#00007F", 
                                "blue", 
                                "#007FFF", 
                                "cyan", 
                                "#7FFF7F", 
                                "yellow", 
                                "#FF7F00", 
                                "red", 
                                "#7F0000"))(n)
}

#' PlotStars 
#' 
#' Plot star charts
#' 
#' Plot FlowSOM grid or tree, where each node is represented by 
#' a star chart indicating median marker values
#' 
#' @param fsom FlowSOM object, as generated by \code{\link{BuildMST}}
#' @param markers Markers to plot (will be parsed by GetChannels)
#' @param colorPalette ColorPalette to use
#' @param list_insteadof_ggarrange If FALSE (default), the plot and the legend
#'                                 are combined by ggarrange. If TRUE, the
#'                                 seperate elements are returned in a list,
#'                                 to allow further customization.
#' @param ...  Additional arguments to pass to \code{\link{PlotFlowSOM}}
#' 
#' @return Nothing is returned. A plot is drawn in which each node is 
#' represented by a star chart indicating the median fluorescence intensities.
#' Resets the layout back to 1 plot at the end.
#
#' @seealso \code{\link{PlotMarker}}, \code{\link{PlotVariable}},
#' \code{\link{PlotFlowSOM}}, \code{\link{PlotLabels}}, 
#' \code{\link{PlotNumbers}}, \code{\link{PlotPies}}, 
#' \code{\link{QueryStarPlot}}, \code{\link{PlotSD}}
#' 
#' @examples
#' # Read from file, build self-organizing map and minimal spanning tree
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' flowSOM.res <- FlowSOM(fileName, compensate = TRUE, transform = TRUE,
#'                        scale = TRUE, colsToUse = c(9, 12, 14:18))
#' 
#' # Plot stars indicating the MFI of the cells present in the nodes
#' PlotStars(flowSOM.res)
#' 
#' @importFrom ggpubr get_legend ggarrange
#' @importFrom ggforce geom_circle
#' 
#' @export
PlotStars <- function(fsom,
                      markers = fsom$map$colsUsed,
                      colorPalette = FlowSOM_colors,
                      list_insteadof_ggarrange = FALSE,
                      ...){
  fsom <- UpdateFlowSOM(fsom)
  
  channels <- GetChannels(fsom, markers)
  
  p <- PlotFlowSOM(fsom = fsom, 
                   ...)
  
  if(!is.null(names(colorPalette))) {
    names(colorPalette) <- GetChannels(fsom, names(colorPalette))
  }
  p <- AddStars(p = p, 
                fsom = fsom, 
                markers = channels, 
                colorPalette = colorPalette)
  
  if(!is.null(names(colorPalette))) {
    names(colorPalette) <- fsom$prettyColnames[GetChannels(fsom, 
                                                           names(colorPalette))]
  }
  l1 <- PlotStarLegend(fsom$prettyColnames[channels], 
                       colorPalette)
  
  l2 <- ggpubr::get_legend(p)
  
  if (list_insteadof_ggarrange){
    p <- p + ggplot2::theme(legend.position = "none")
    l2 <- ggpubr::as_ggplot(l2)
    return(list(tree = p, 
                starLegend = l1, 
                backgroundLegend = l2))
  } else {
    p <- ggpubr::ggarrange(p,
                           ggpubr::ggarrange(l1, l2,
                                             ncol = 1),
                           NULL,
                           ncol = 3, widths = c(3, 1, 0.3),
                           legend = "none")
    return(p)
  }
}


#' PlotPies 
#' 
#' Plot comparison with other clustering
#' 
#' Plot FlowSOM grid or tree, with pies indicating another clustering
#'    or manual gating result
#'
#' @param fsom         FlowSOM object, as generated by \code{\link{FlowSOM}}
#' @param cellTypes    Array of factors indicating the celltypes
#' @param colorPalette Color palette to use.
#' @param ...          Additional arguments to pass to \code{\link{PlotFlowSOM}}
#' 
#' @return Ggplot plot
#' 
#' @seealso \code{\link{PlotStars}}, \code{\link{PlotVariable}},
#' \code{\link{PlotFlowSOM}}, \code{\link{PlotLabels}}, 
#' \code{\link{PlotNumbers}}, \code{\link{PlotMarker}},
#'  \code{\link{QueryStarPlot}}, \code{\link{PlotSD}}
#' 
#' @examples 
#' # Identify the files
#' fcs_file <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' wsp_file <- system.file("extdata", "gating.wsp", package = "FlowSOM")
#' 
#' # Specify the cell types of interest for assigning one label per cell
#' cellTypes <- c("B cells",
#'                 "gd T cells", "CD4 T cells", "CD8 T cells",
#'                 "NK cells", "NK T cells")
#'
#' # Parse the FlowJo workspace     
#' gatingResult <- GetFlowJoLabels(fcs_file, wsp_file,
#'                                 cellTypes = cellTypes)
#'
#' # Check the number of cells assigned to each gate
#' colSums(gatingResult$matrix)
#' 
#' # Build a FlowSOM tree
#' flowSOM.res <- FlowSOM(fcs_file,
#'                        scale = TRUE, 
#'                        compensate = TRUE, 
#'                        transform = TRUE,
#'                        toTransform = 8:18, 
#'                        colsToUse = c(9, 12, 14:18),
#'                        nClus = 10,
#'                        seed = 1)
#'    
#'  # Plot pies indicating the percentage of cell types present in the nodes
#'  PlotPies(flowSOM.res,
#'           gatingResult$manual,
#'           backgroundValues = flowSOM.res$metaclustering)
#'           
#' @importFrom grDevices colorRampPalette
#' 
#' @export
PlotPies <- function(fsom,
                     cellTypes,
                     colorPalette = grDevices::colorRampPalette(c("white",
                                                                  "#00007F", 
                                                                  "blue", 
                                                                  "#007FFF", 
                                                                  "cyan", 
                                                                  "#7FFF7F", 
                                                                  "yellow", 
                                                                  "#FF7F00", 
                                                                  "red", 
                                                                  "#7F0000")),
                     ...){
  if(length(cellTypes) != nrow(fsom$data)){
    stop("There are ", nrow(fsom$data), " cells, while you provided ",
         length(cellTypes), " labels. These numbers should match.")
  }
  p <- PlotFlowSOM(fsom = fsom, 
                   ...)
  
  p <- AddPies(p, 
               fsom = fsom,
               colorPalette = colorPalette,
               cellLabels = cellTypes)
  return(p)
}

#' PlotLabels 
#' 
#' Plot labels for each cluster
#' 
#' Plot FlowSOM grid or tree, with in each node a label. 
#' Especially useful to show metacluster numbers
#'
#' @param fsom        FlowSOM object, as generated by \code{\link{FlowSOM}}
#' @param labels      A vector of labels for every node.
#' @param maxNodeSize Determines the maximum nodesize. Default is 0.
#' @param textSize    Size for geom_text. Default (=3.88) is from geom_text.
#' @param textColor   Color for geom_text. Default = black.
#' @param ...         Additional arguments to pass to \code{\link{PlotFlowSOM}}
#' 
#' @return Nothing is returned. A plot is drawn in which each node is 
#' represented by a label.
#' 
#' 
#' @seealso \code{\link{PlotStars}}, \code{\link{PlotVariable}},
#' \code{\link{PlotFlowSOM}}, \code{\link{PlotMarker}}, 
#' \code{\link{PlotNumbers}}, \code{\link{PlotPies}}, 
#' \code{\link{QueryStarPlot}}, \code{\link{PlotSD}}
#' 
#' @examples 
#' # Read from file, build self-organizing map and minimal spanning tree
#' fileName <- system.file("extdata", "68983.fcs", package="FlowSOM")
#' ff <- flowCore::read.FCS(fileName)
#' ff <- flowCore::compensate(ff, flowCore::keyword(ff)[["SPILL"]])
#' ff <- flowCore::transform(ff,
#'          flowCore::transformList(colnames(flowCore::keyword(ff)[["SPILL"]]),
#'                                 flowCore::logicleTransform()))
#' flowSOM.res <- FlowSOM(ff,
#'                        scale = TRUE,
#'                        colsToUse = c(9, 12, 14:18),
#'                        nClus = 10,
#'                        seed = 1)
#' 
#' # Plot the node IDs
#' PlotLabels( flowSOM.res, 
#'             flowSOM.res$metaclustering)
#' 
#' @importFrom grDevices colorRampPalette
#' 
#' @export
PlotLabels <- function(fsom,
                       labels, 
                       maxNodeSize = 0,
                       textSize = 3.88,
                       textColor = "black",
                       ...){
  fsom <- UpdateFlowSOM(fsom)
  
  if(length(labels) != NClusters(fsom)) {
    stop(paste0("Length of 'labels' should be equal to number of clusters in ", 
                "FlowSOM object (", NClusters(fsom), " clusters, ", 
                length(labels), " labels)."))
  }
  
  p <- PlotFlowSOM(fsom = fsom,
                   maxNodeSize = maxNodeSize, 
                   ...)
  p <- AddLabels(p,
                 labels = labels,
                 textSize = textSize,
                 textColor = textColor)
  return(p)
}

#' PlotNumbers 
#' 
#' Plot cluster ids for each cluster
#' 
#' Plot FlowSOM grid or tree, with in each node the cluster id.
#'
#' @param fsom        FlowSOM object
#' @param level       Character string, should be either "clusters" or 
#'                    "metaclusters"
#' @param maxNodeSize Determines the maximum nodesize. Default is 0. 
#'                    See \code{\link{PlotFlowSOM}} for more options.
#' @param ...         Additional arguments to pass to \code{\link{PlotLabels}} 
#'                    and to \code{\link{PlotFlowSOM}} 
#' 
#' @return Nothing is returned. A plot is drawn in which each node is 
#' labeled by its cluster id.
#' 
#' @seealso \code{\link{PlotStars}}, \code{\link{PlotVariable}},
#' \code{\link{PlotFlowSOM}}, \code{\link{PlotLabels}}, 
#' \code{\link{PlotMarker}}, \code{\link{PlotPies}},
#'  \code{\link{QueryStarPlot}}, \code{\link{PlotSD}}
#' 
#' @examples 
#' # Read from file, build self-organizing map and minimal spanning tree
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' ff <- flowCore::read.FCS(fileName)
#' ff <- flowCore::compensate(ff, flowCore::keyword(ff)[["SPILL"]])
#' ff <- flowCore::transform(ff, flowCore::estimateLogicle(ff,
#'                                                flowCore::colnames(ff)[8:18]))
#' flowSOM.res <- FlowSOM(ff,
#'                        scale = TRUE,
#'                        colsToUse = c(9, 12, 14:18),
#'                        nClus = 10,
#'                        seed = 1)
#' 
#' # Plot the node IDs
#' PlotNumbers(flowSOM.res)
#' PlotNumbers(flowSOM.res, "metaclusters")
#' 
#' PlotNumbers(flowSOM.res,
#'             view = "grid")
#' 
#' PlotNumbers(flowSOM.res,
#'             maxNodeSize = 1,
#'             equalNodeSize = TRUE)
#' 
#' @export
PlotNumbers <- function(fsom,
                        level = "clusters",
                        maxNodeSize = 0,
                        ...){
  if (level == "clusters"){
    numbers <- seq_len(NClusters(fsom))
  } else if (level == "metaclusters") {
    numbers <- fsom$metaclustering
  } else {
    stop("level should be \"clusters\" or \"metaclusters\"")
  }
  p <- PlotLabels(fsom =  fsom,
                  labels = numbers, 
                  maxNodeSize = maxNodeSize,
                  ...)
  return(p)
}

#' PlotMarker 
#' 
#' Plot comparison with other clustering
#' 
#' Plot FlowSOM grid or tree, coloured by node values for a specific marker
#'
#' @param fsom         FlowSOM object
#' @param marker       A vector of markers/channels to plot.
#' @param refMarkers  Is used to determine relative scale of the marker 
#'                     that will be plotted. Default are all markers used in the
#'                     clustering.
#' @param title        A vector with custom titles for the plot. Default is 
#'                     the marker name.
#' @param colorPalette Colorpalette to use. Can be a function or a vector.
#' @param lim          Limits for the scale
#' @param ...     Additional arguments to pass to \code{\link{PlotFlowSOM}},
#'                e.g. view, backgroundValues, equalNodeSize ...
#' 
#' @return A ggplot figure is returned in which every cluster is colored
#'         according to the MFI value for the specified marker
#' 
#' @seealso \code{\link{PlotStars}}, \code{\link{PlotVariable}}
#' 
#' @examples 
#' # Build FlowSOM model
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' flowSOM.res <- FlowSOM(fileName, 
#'                        compensate = TRUE, transform = TRUE, scale = FALSE,
#'                        colsToUse = c(9, 12, 14:18),
#'                        nClus = 10,
#'                        seed = 1)
#' # Plot one marker
#' PlotMarker(flowSOM.res, 
#'            "CD19")
#'            
#' PlotMarker(flowSOM.res, 
#'            "CD19",
#'            colorPalette = c("grey", "red"))
#'            
#' # Plot all markers
#' PlotMarker(flowSOM.res,
#'            c(9, 12, 14:18))
#' 
#' # Use specific limits if the ones from the columns used for clustering
#' # are not relevant for  your marker of choice
#' PlotMarker(flowSOM.res, 
#'            "FSC-A",
#'             lim = c(55000, 130000))
#' 
#' # Example with additional FlowSOM plotting options
#' PlotMarker(flowSOM.res, 
#'            "CD19",
#'            view = "grid",
#'            equalNodeSize = TRUE,
#'            backgroundValues = 1:100 == 27,
#'            backgroundColors = c("white", "red"))
#'            
#'            
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 ggtitle
#' 
#' @export
PlotMarker <- function(fsom,
                       marker,
                       refMarkers = fsom$map$colsUsed,
                       title = GetMarkers(fsom, marker),
                       colorPalette = FlowSOM_colors,
                       lim = NULL,
                       ...){
  
  fsom <- UpdateFlowSOM(fsom)
  
  # Get median values to visualise
  mfis <- GetClusterMFIs(fsom)
  
  # Get column names
  channels <- GetChannels(fsom, marker)

  # Compute limits based on all reference markers
  ref_channels <- GetChannels(fsom, refMarkers)
  if (is.null(lim)) lim <- c(min(mfis[, ref_channels]), 
                             max(mfis[, ref_channels]))
  
  plotList <- lapply(seq_along(channels), function(channelI){
    
    # Use MFI values as variable to plot
    p <- PlotVariable(fsom, 
                      variable = mfis[, channels[channelI]],
                      variableName = "MFI",
                      colorPalette = colorPalette,
                      lim = lim,
                      ...)
    
    # Add title
    if (is.na(title[channelI])) title[channelI] <- ""
    p <- p + ggplot2::ggtitle(title[channelI])
  })

  p <- ggpubr::ggarrange(plotlist = plotList, common.legend = TRUE, 
                         legend = "right")
  
  return(p)
}

#' PlotVariable 
#' 
#' Plot a variable for all nodes
#' 
#' Plot FlowSOM grid or tree, coloured by node values given in variable
#'
#' @param fsom         FlowSOM object
#' @param variable     A vector containing a value for every cluster
#' @param variableName Label to show on the legend
#' @param colorPalette Colorpalette to use. Can be a function or a vector.
#' @param lim          Limits for the scale
#' @param ...     Additional arguments to pass to \code{\link{PlotFlowSOM}},
#'                e.g. view, backgroundValues, equalNodeSize ...
#'         
#' @seealso \code{\link{PlotStars}}, \code{\link{QueryStarPlot}},
#' \code{\link{PlotFlowSOM}}, \code{\link{PlotLabels}}, 
#' \code{\link{PlotNumbers}}, \code{\link{PlotMarker}}, 
#' \code{\link{PlotPies}}, \code{\link{PlotSD}}
#'         
#' @examples 
#' # Build FlowSOM model
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' flowSOM.res <- FlowSOM(fileName, 
#'                        compensate = TRUE, transform = TRUE, scale = FALSE,
#'                        colsToUse = c(9, 12, 14:18),
#'                        nClus = 10,
#'                        seed = 1)
#'                        
#' # Plot some random values
#' rand <- runif(flowSOM.res$map$nNodes)
#' PlotVariable(flowSOM.res, 
#'              variable = rand,
#'              variableName = "Random")
#'    
#' @export
PlotVariable <- function(fsom,
                         variable,
                         variableName = "",
                         colorPalette = FlowSOM_colors,
                         lim = NULL,
                         ...){
  fsom <- UpdateFlowSOM(fsom)
  if (length(variable) != fsom$map$nNodes){
    stop(paste0("Length of 'variable' should be equal to number of clusters in", 
                " FlowSOM object (", fsom$map$nNodes, " clusters and ", 
                length(variable), " variables)."))
  }
  
  # Base plot
  p <- PlotFlowSOM(fsom = fsom,
                   ...)
  p <- AddNodes(p = p, 
                values = variable,
                colorPalette = colorPalette,
                lim = lim,
                label = variableName)
  
  return(p)
}

#' PlotDimRed 
#' 
#' Plot a dimensionality reduction
#' 
#' Plot a dimensionality reduction of fsom$data
#'
#' @param fsom            FlowSOM object, as generated by \code{\link{BuildMST}}
#' @param colsToUse       The columns used for the dimensionality reduction. 
#'                        Default = fsom$map$colsUsed.
#' @param colorBy         Defines how the dimensionality reduction will be 
#'                        colored. Can be "metaclusters" (default), "clusters"
#'                        or a marker/channel/index.
#' @param cTotal          The total amount of cells to be used in the 
#'                        dimensionality reduction. Default is all the cells.
#' @param dimred          A dimensionality reduction function. 
#'                        Default = Rtsne::Rtsne. Alternatively, a data.frame or
#'                        matrix with either equal number of rows to the
#'                        fsom or an OriginalID column. Recommended to put
#'                        cTotal to NULL when providing a matrix (or ensuring
#'                        that the dimred corresponds to subsampling the
#'                        flowSOM data for cTotal cells with the same seed).
#' @param extractLayout   A function to extract the coordinates from the results
#'                        of the dimred default = function(dimred){dimred$Y}.
#' @param label           If label = TRUE (default), labels are added to plot.
#' @param returnLayout    If TRUE, this function returns a dataframe with 
#'                        the layout of dimred and the original IDs and the 
#'                        plot. Default = FALSE.                     
#' @param seed            A seed for reproducibility.
#' @param title           A title for the plot. 
#' @param ...             Additional arguments to pass to dimred.                  
#' 
#' @return A dimensionality reduction plot made in ggplot2
#'         
#' @examples 
#'    file <- system.file("extdata", "68983.fcs", package="FlowSOM")
#'    flowSOM.res <- FlowSOM(file, compensate = TRUE, transform = TRUE, 
#'                   scale = TRUE,
#'                   colsToUse = c(9, 12, 14:18), nClus = 10, silent = FALSE,
#'                   xdim = 7, ydim = 7)
#'    PlotDimRed(flowSOM.res, cTotal = 5000, seed = 1, title = "t-SNE")
#'    PlotDimRed(flowSOM.res, cTotal = 5000, colorBy = "CD3", seed = 1, 
#'               title = "t-SNE")
#'    
#' @import     ggplot2
#' @importFrom Rtsne Rtsne
#' @importFrom scattermore geom_scattermore
#' @importFrom tidyr pivot_longer
#'        
#' @export
PlotDimRed <- function(fsom,
                       colsToUse = fsom$map$colsUsed,
                       colorBy = "metaclusters",
                       cTotal = NULL,
                       dimred = Rtsne::Rtsne,
                       extractLayout  = function(dimred){dimred$Y},
                       label = TRUE,
                       returnLayout = FALSE,
                       seed = NULL,
                       title = NULL,
                       ...){
  dimred_data <- fsom$data
  if (length(colorBy) == 1 && colorBy == "metaclusters") {
    if (is.null(fsom$metaclustering)) stop("No metaclustering present")
    dimred_col <- as.data.frame(GetMetaclusters(fsom))
  } else if (length(colorBy) == 1 && colorBy == "clusters") {
    dimred_col <- as.data.frame(as.factor(GetClusters(fsom)))
  } else if (all(colorBy %in% colnames(dimred_data)) ||
             all(colorBy %in% GetMarkers(fsom, colnames(dimred_data))) ||
             all(colorBy %in% seq_len(ncol(dimred_data)))){
    dimred_col <- fsom$data[, GetChannels(fsom, colorBy), drop = FALSE]
    colnames(dimred_col) <- GetMarkers(fsom, colnames(dimred_col))
    colorBy <- "marker"
  } else stop(paste0("colorBy should be \"metaclusters\", \"clusters\" or a ",
                     "vector of channels, markers or indices"))
  if (!is.null(colsToUse)) dimred_data <- dimred_data[, GetChannels(fsom, 
                                                                    colsToUse)]
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(cTotal) && cTotal < nrow(dimred_data)) {
    downsample <- sample(seq_len(nrow(dimred_data)), cTotal)
    dimred_data <- dimred_data[downsample, , drop = FALSE]
    dimred_col <- dimred_col[downsample, , drop = FALSE]
  } else {
    downsample <- seq_len(nrow(dimred_data))
  }
  if (is.function(dimred)){
    dimred_res <- dimred(dimred_data, ...)
    dimred_layout <- as.data.frame(extractLayout(dimred_res))
    if (nrow(dimred_layout) == 0 && ncol(dimred_layout) != 2) {
      stop("Please use the right extraction function in extractLayout")
    }
  } else if((is.matrix(dimred) | is.data.frame(dimred)) & 
            (nrow(dimred) == nrow(dimred_data) | 
             any(colnames(dimred) == "Original_ID"))){
    id_col <- which(colnames(dimred) == "Original_ID")
    dimred_layout <- as.data.frame(dimred[,-id_col])
    if("Original_ID" %in% colnames(dimred)){
      dimred_data <- dimred_data[dimred[,"Original_ID"], , drop = FALSE]
      dimred_col <- dimred_col[dimred[,"Original_ID"], , drop = FALSE]
    }
  } else stop("dimred should be a dimensionality reduction method or matrix")
  
  colnames(dimred_layout) <- c("dimred_1", "dimred_2")
  dimred_plot <- cbind(dimred_layout, dimred_col)
  
  if (colorBy == "marker"){
    dimred_plot <- dimred_plot %>% tidyr::pivot_longer(3:ncol(dimred_plot), 
                                                       names_to = "markers")
    p <- ggplot2::ggplot(dimred_plot) + 
      scattermore::geom_scattermore(ggplot2::aes(x = .data$dimred_1, 
                                                 y = .data$dimred_2, 
                                                 col = .data$value), 
                                    pointsize = 1) +
      ggplot2::facet_wrap(~markers) +
      ggplot2::theme_minimal() +
      ggplot2::coord_fixed() +
      ggplot2::scale_color_gradientn(colours = FlowSOM_colors(9))
  } else {
    colnames(dimred_plot) <- c("dimred_1", "dimred_2", "colors")
    
    median_x <- tapply(dimred_plot[,"dimred_1"], dimred_plot[,"colors"], median)
    median_y <- tapply(dimred_plot[,"dimred_2"], dimred_plot[,"colors"], median)
    
    p <- ggplot2::ggplot(dimred_plot) +
      scattermore::geom_scattermore(ggplot2::aes(x = .data$dimred_1, 
                                                 y = .data$dimred_2, 
                                                 col = .data$colors), 
                                    pointsize = 1) +
      ggplot2::theme_minimal() +
      ggplot2::coord_fixed() 
    if (label){
      p <- p + ggrepel::geom_label_repel(aes(x = .data$x, 
                                             y = .data$y, 
                                             label = .data$label, 
                                             color = .data$label),
                                         data = data.frame(x = median_x,
                                                           y = median_y,
                                                           label = names(median_x)),
                                         segment.color = "grey", force = 20, 
                                         segment.size = 0.2, point.padding = 0.5)+
        labs(col = colorBy)
    }
  }
  if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
  if (returnLayout) {
    if (!is.null(cTotal)){
      dimred_layout <- data.frame(dimred_layout, "Original_ID" = downsample)
    }
    return(list("plot" = p, "layout" = dimred_layout))
  } else return(p)
}

#' QueryStarPlot 
#' 
#' Query a certain cell type
#' 
#' Identify nodes in the tree which resemble a certain profile of "high"
#' or "low" marker expressions.
#'
#' @param fsom            FlowSOM object, as generated by \code{\link{BuildMST}}
#' @param query           Array containing "high" or "low" for the specified 
#'                        column names of the FlowSOM data.
#' @param plot            If true, a plot with a gradient of scores for the 
#'                        nodes is shown.
#' @param colorPalette    ColorPalette to be used for colors for "stars", 
#'                        "pies" or "marker". Can be either a function or an 
#'                        array specifying colors.
#' @param backgroundColors Color to use for nodes with a high score in the plot.
#'                         Default is red.
#' @param ...         Additional arguments to pass to \code{\link{PlotFlowSOM}}
#' 
#' @return A list, containing the ids of the selected nodes, the individual 
#'         scores for all nodes and the scores for each marker for each node
#'         
#' @seealso \code{\link{PlotStars}}, \code{\link{PlotVariable}},
#' \code{\link{PlotFlowSOM}}, \code{\link{PlotLabels}}, 
#' \code{\link{PlotNumbers}}, \code{\link{PlotMarker}}, 
#' \code{\link{PlotPies}}, \code{\link{PlotSD}}
#'         
#' @examples 
#'    file <- system.file("extdata", "68983.fcs", package="FlowSOM")
#'    flowSOM.res <- FlowSOM(file, compensate = TRUE, transform = TRUE, 
#'                   scale = TRUE, colsToUse = c(9, 12, 14:18), nClus = 10, 
#'                   silent = FALSE, xdim = 7, ydim = 7)
#'    query <- c("CD3" = "high", #CD3
#'               "CD4" = "low", #TCRb
#'               "CD8" = "high") #CD8
#'    query_res <- QueryStarPlot(flowSOM.res, query, equalNodeSize = TRUE)
#'    
#'    cellTypes <- factor(rep("Unlabeled", 49), 
#'                        levels = c("Unlabeled", "CD8 T cells"))
#'    cellTypes[query_res$selected] <- "CD8 T cells"
#'    PlotStars(flowSOM.res,
#'              backgroundValues = cellTypes,
#'              backgroundColors = c("#FFFFFF00", "#ca0020aa"))
#'    
#' @importFrom grDevices colorRampPalette
#' @importFrom ggpubr ggarrange
#'        
#' @export
QueryStarPlot <- function(fsom,
                          query,
                          plot = TRUE,
                          colorPalette = FlowSOM_colors,
                          backgroundColors = "#CA0020",
                          ...){
  fsom <- UpdateFlowSOM(fsom)
  names(query) <- GetChannels(fsom, names(query))
  markers <- names(query)
  lowMarkers <- which(query == "low")
  query <- ParseQuery(fsom, query)
  if (plot){
    fP <- PlotStars(fsom = fsom,
                    markers = markers,
                    backgroundColors = 
                      grDevices::colorRampPalette(c(rep("#FFFFFF00", 3), 
                                                    backgroundColors)),
                    backgroundValues = query$nodeScores,
                    list_insteadof_ggarrange = TRUE,
                    ...)
    nMarkers <- length(markers)
    starHeight <- rep(1, nMarkers)
    starHeight[lowMarkers] <- 0
    l <- PlotStarLegend(fsom$prettyColnames[markers], 
                        colors = colorPalette(nMarkers), 
                        starHeight = starHeight)
    p <- ggpubr::ggarrange(fP$tree, 
                           ggpubr::ggarrange(l, fP$backgroundLegend, ncol = 1), 
                           NULL,
                           ncol = 3, 
                           widths = c(3, 1, 0.3), 
                           legend = "none")
    query$plot <- p
  }
  return(query)
}

#' QueryMultiple
#'
#' Function which takes a named list of multiple cell types, where every item is
#' a named vector with values "high"/"low" and the names correspond to the
#' markers or channels (e.g. as generated by parse_markertable).
#'
#' @param fsom       FlowSOM object
#' @param cellTypes Description of the cell types. Named list, with one named
#'                   vector per cell type containing "high"/"low" values
#' @param plotFile   Path to a pdf file to save the plots (default is 
#'                   queryMultiple.pdf). If \code{NULL}, no plots will be 
#'                   generated
#' @param ...        Additional arguments to pass to \code{\link{QueryStarPlot}}
#'
#' @return A label for every FlowSOM cluster (Unknown or one of the celltype
#'         names of the list, if selected by QueryStarPlot)
#' 
#' @examples
#'    file <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#'    ff <- flowCore::read.FCS(file)
#'    # Use the wrapper function to build a flowSOM object (saved in flowSOM.res)
#'    # and a metaclustering (saved in flowSOM.res[["metaclustering"]])
#'    flowSOM.res <- FlowSOM(ff, compensate = TRUE, transform = TRUE, scale = TRUE,
#'                   colsToUse = c(9, 12, 14:18), nClus = 10, silent = FALSE,
#'                   xdim = 7, ydim = 7)
#'    cellTypes <- list("CD8 T cells" = c("PE-Cy7-A" = "high",
#'                                         "APC-Cy7-A" = "high",
#'                                         "Pacific Blue-A" = "high"),
#'                        "B cells" = c("PE-Cy5-A" = "high"),
#'                        "NK cells" = c("PE-A" = "high",
#'                                       "PE-Cy7-A" = "low",
#'                                       "APC-Cy7-A" = "low"))
#'    query_res <- QueryMultiple(flowSOM.res, cellTypes, "query_multiple.pdf")
#'    
#' @export
QueryMultiple <- function(fsom,
                          cellTypes,
                          plotFile = "queryMultiple.pdf",
                          ...){
  fsom <- UpdateFlowSOM(fsom)
  labels <- rep("Unlabeled", fsom$map$nNodes)
  plotList <- list()
  plot <- ifelse(is.null(plotFile), FALSE, TRUE)
  for(cell_type in names(cellTypes)){
    query <- cellTypes[[cell_type]]
    query_res <- QueryStarPlot(fsom, 
                               equalNodeSize = TRUE,
                               query = query,
                               title = cell_type,
                               plot = plot,
                               ...)
    if (plot) plotList[[cell_type]] <- query_res$plot
    labels[query_res$selected] <- cell_type
  }
  if (plot){
    pdf(plotFile, useDingbats = FALSE)
    for (p in plotList){
      print(p)
    }
    dev.off()
  }
  return(labels)
}

#' PlotSD 
#' 
#' Plot FlowSOM grid or tree, coloured by standard deviaton.
#'
#' @param fsom      FlowSOM object, as generated by \code{\link{FlowSOM}}
#' @param marker    If a marker/channel is given, the sd for this marker is 
#'                  shown. Otherwise, the maximum ratio is used.
#' @param ...       Additional arguments to pass to \code{\link{PlotFlowSOM}}
#' 
#' @return Nothing is returned. A plot is drawn in which each node 
#'         is coloured depending on its standard deviation
#'         
#' @seealso \code{\link{PlotStars}}, \code{\link{PlotVariable}},
#' \code{\link{PlotFlowSOM}}, \code{\link{PlotLabels}}, 
#' \code{\link{PlotNumbers}}, \code{\link{PlotMarker}}, 
#' \code{\link{PlotPies}}, \code{\link{QueryStarPlot}}
#'         
#' @examples 
#' # Read from file, build self-organizing map and minimal spanning tree
#' fileName <- system.file("extdata", "68983.fcs", package  = "FlowSOM")
#' flowSOM.res <- ReadInput(fileName, compensate  = TRUE, transform  = TRUE,
#'                          scale  = TRUE)
#' flowSOM.res <- BuildSOM(flowSOM.res, colsToUse  = c(9, 12, 14:18))
#' flowSOM.res <- BuildMST(flowSOM.res)
#' 
#' PlotSD(flowSOM.res)
#' 
#' @importFrom grDevices colorRampPalette
#' 
#' @export
PlotSD <- function(fsom,
                   marker = NULL,
                   ...){
  fsom <- UpdateFlowSOM(fsom)
  if(!is.null(marker)) marker <- GetChannels(fsom, marker)
  SDs <- ParseSD(fsom, marker)
  PlotVariable(fsom = fsom, 
               variable = SDs,
               variableName = "SD", 
               ...)
}

#' PlotStarLegend 
#' 
#' Plots star legend
#' 
#' Function makes the legend of the FlowSOM star plot
#' 
#' @param markers           Vector of markers used in legend
#' @param colors            ColorPalette for the legend. Can be a vector or a 
#'                          function.
#' @param starHeight        Star height. Default = 1.
#'                          
#' @return Returns nothing, but plots a legend for FlowSOM star plot
#' 
#' @seealso \code{\link{PlotFlowSOM}}
#' 
#' @examples 
#' PlotStarLegend(c("CD3", "CD4", "CD8"),
#'                FlowSOM_colors(3))
#' 
#' @import ggplot2
#' @importFrom dplyr count
#' @importFrom grDevices colorRampPalette
#' 
#' @export
PlotStarLegend <- function(markers, 
                           colors, 
                           starHeight = 1){
  
  requireNamespace("ggplot2")
  markers <- factor(markers, levels = markers)
  nMarkers <- length(markers)
  circularCoords <- seq(from = 2 * pi / (nMarkers * 2), 
                        by = 2 * pi / nMarkers, 
                        length.out = nMarkers)
  dfSegments <- data.frame(Markers = markers,
                           x = sin(circularCoords),
                           y = cos(circularCoords),
                           xend = 1.1 *  ifelse(sin(circularCoords) >= 0, 
                                                1, -1),
                           yend = NA)
  nLeftRight <- dplyr::count(dfSegments, .data$x >= 0)
  if (nrow(nLeftRight) != 1){
    by <- ifelse(nMarkers <= 8, 1, 0.65)
    left <- seq(from = 0, by = by, length.out = nLeftRight[1, 2])
    right <- -seq(from = 0, by = by, length.out = nLeftRight[2, 2])
    dfSegments[which(dfSegments$x < 0), ]$yend <- left - mean(left)
    dfSegments[which(dfSegments$x >= 0), ]$yend <- right - mean(right)
  } else {
    dfSegments$yend <- -2
  }
  horizontalLines <- data.frame(Markers = dfSegments$Markers,
                                x = dfSegments$xend, y = dfSegments$yend,
                                xend = dfSegments$xend + 
                                  ifelse(dfSegments$xend >= 0, 0.5, -0.5),
                                yend = dfSegments$yend, 
                                stringsAsFactors = FALSE)
  dfSegments <- rbind(dfSegments, horizontalLines)
  dfLabels <- data.frame(x = horizontalLines$xend + 
                           ifelse(horizontalLines$xend >= 0, 0.3, -0.3),
                         y = horizontalLines$yend)
  markers_tmp <- rep(1, nMarkers)
  names(markers_tmp) <- markers
  dfStar <- ParseArcs(0, 0, markers_tmp, starHeight)
  dfStar$Markers <- factor(dfStar$Marker, levels = markers)
  l <- ggplot2::ggplot() + 
    ggplot2::xlim(c(-6, 6)) + 
    ggplot2::coord_fixed(clip = "off") + 
    ggplot2::theme_void() + 
    ggplot2::theme(legend.position = "none")
  
  l <- AddStarsPies(l, 
                    dfStar, 
                    colors)
  
  l <- AddScale(p = l,
                values = markers,
                colors = colors,
                showLegend = FALSE,
                type = "color")
  l <- l + 
    ggplot2::geom_segment(data = dfSegments, 
                          ggplot2::aes(x = .data$x, y = .data$y, 
                                       xend = .data$xend, 
                                       yend = .data$yend, 
                                       color = .data$Markers), 
                          size = 0.6)
  
  l <- AddLabels(l, 
                 labels = markers, 
                 layout = dfLabels, 
                 hjust = ifelse(horizontalLines$xend >= 0, 0, 1))
  return(l)
}

#' ParseEdges 
#' 
#' Parses edges
#' 
#' Function that parses the graph edges of the FlowSOM object into a dataframe
#' 
#' @param fsom  FlowSOM object, as generated by \code{\link{FlowSOM}}
#'                          
#' @return A dataframe consisting of start and end coordinates of edges
#' 
#' @seealso \code{\link{PlotFlowSOM}}, \code{\link{ParseNodeSize}},
#' \code{\link{ParseArcs}}, \code{\link{ParseQuery}},
#' \code{\link{ParseSD}}, \code{\link{AddMST}}
#' 
#' @importFrom igraph as_edgelist
#' 
ParseEdges <- function(fsom){
  edgeList <- as.data.frame(igraph::as_edgelist(fsom$MST$g),
                            stringsAsFactors = FALSE)
  coords <- fsom$MST$l
  segmentPlot <- lapply(seq_len(nrow(edgeList)), function(row_id){
    node_ids <- as.numeric(edgeList[row_id, ])
    row <- c(coords[node_ids[1], 1],
             coords[node_ids[1], 2],
             coords[node_ids[2], 1],
             coords[node_ids[2], 2])
    return(row)
  })
  segmentPlot <- do.call(rbind, segmentPlot)
  colnames(segmentPlot) <- c("x", "y", "xend", "yend")
  return(as.data.frame(segmentPlot))
}

#' ParseLayout
#' 
#' @param fsom FlowSOM object
#' @param layout "MST", "grid" or a matrix/dataframe with 2 columns and 1 row
#'               per cluster
#'               
#' @return dataframe with 2 columns and 1 row per cluster
#' 
ParseLayout <- function(fsom, layout){
  if (is.matrix(layout) || is.data.frame(layout)){
    if (nrow(layout) == NClusters(fsom) && 
        ncol(layout) == 2){
      layout <- as.data.frame(layout)
    } else {
      stop(paste0("Number of rows should be equal to number of clusters (", 
                  NClusters(fsom), " clusters) and number of columns should be equal ",
                  "to two."))
    }
  } else if (layout == "grid"){
    layout <- as.data.frame(fsom$map$grid)
  } else if (layout == "MST"){
    layout <- as.data.frame(fsom$MST$l)
  } else {
    stop(paste0("The layout should be 'MST', 'grid' or a matrix/dataframe ",
                "with 2 columns and 1 row per cluster."))
  }
  colnames(layout) <- c("x", "y")
  return(layout)
}

#' ParseNodeSize 
#' 
#' Parses node size
#' 
#' Function that parses the mapping of the FlowSOM object into node sizes 
#' relative to the abundances of cells per cluster
#' 
#' Scales node size relative to the abundances of cells per cluster
#' 
#' @param nodeSizes         A vector with nodesizes
#' @param maxNodeSize       Determines the maximum nodesize. 
#' @param refNodeSize       Reference for nodesize against which the nodeSizes 
#'                          will be scaled. Default = max(nodeSizes)
#'                          
#' @return A vector is returned consisting of node sizes
#' 
#' @seealso \code{\link{PlotFlowSOM}}, \code{\link{ParseEdges}}, 
#' \code{\link{AutoMaxNodeSize}}, \code{\link{ParseArcs}},
#' \code{\link{ParseQuery}}, \code{\link{ParseSD}}
#' 
#' @importFrom dplyr count mutate pull
#' 
ParseNodeSize <- function(nodeSizes, maxNodeSize, refNodeSize){
  if(any(is.na(nodeSizes))) nodeSizes[is.na(nodeSizes)] <- 0
  nNodes <- length(nodeSizes)
  if (length(unique(nodeSizes)) == 1){
    return(rep(maxNodeSize, nNodes))
  }
  scaledNodeSize <- data.frame(size = nodeSizes) %>% 
    dplyr::mutate(scaled = (.data$size / refNodeSize) * maxNodeSize ^ 2) %>%
    dplyr::mutate(sqrt_scaled = sqrt(.data$scaled)) %>% 
    dplyr::pull(.data$sqrt_scaled)
  return(scaledNodeSize)
}

#' ParseArcs 
#' 
#' Parses stars
#' 
#' Function that parses the FlowSOM object into a dataframe for the star values 
#' for ggplot
#' 
#' @param x              x coordinate of node
#' @param y              y coordinate of node
#' @param arcValues      A named vector with the frequency of how the node 
#'                       should be divided
#' @param arcHeights     The heights of the arcs
#'                          
#' @return A dataframe ready to use with ggplot, consisting of the coordinates 
#'         of centers, the radius and angles of the star values
#' 
#' @seealso \code{\link{PlotFlowSOM}}, \code{\link{ParseEdges}}, 
#' \code{\link{ParseNodeSize}},
#' \code{\link{ParseQuery}}, \code{\link{ParseSD}}
#' 
ParseArcs <- function(x, y, arcValues, arcHeights){
  markers <- names(arcValues)
  arcValues <- c(0, (arcValues / sum(arcValues)) * (2 * pi))
  arcValues <- cumsum(arcValues)
  resDf <- data.frame(Markers = markers, x0 = x, y0 = y, 
                      start = arcValues[-length(arcValues)], 
                      end = arcValues[-1], value = arcHeights)
  return(resDf)
}

#' ParseQuery 
#' 
#' Parses query
#' 
#' Identify nodes in the tree which resemble a certain profile of "high" 
#' or "low" marker expressions.
#' 
#' @param fsom   FlowSOM object, as generated by \code{\link{FlowSOM}}
#' @param query  Array containing "high" or "low" for the specified 
#'               column names of the FlowSOM data
#'                          
#' @return A list, containing the ids of the selected nodes, the individual  
#'         scores for all nodes and the scores for each marker for each node
#' 
#' @seealso \code{\link{PlotFlowSOM}}, \code{\link{ParseEdges}}, 
#' \code{\link{ParseNodeSize}}, \code{\link{ParseArcs}},
#' \code{\link{QueryStarPlot}}, \code{\link{ParseSD}}
#' 
ParseQuery <- function(fsom, query){
  scores <- matrix(NA, ncol = length(query), nrow = nrow(fsom$map$medianValues),
                   dimnames = list(NULL, names(query)))
  for (marker in names(query)) {
    data <- fsom$map$medianValues[, marker]
    if (query[marker] == "high") {
      scores[, marker] = (max(data, na.rm = TRUE) - data) ^ 2
    } else if (query[marker] == "low") {
      scores[, marker] = (data - min(data, na.rm = TRUE)) ^ 2
    } else {
      stop("Only high or low marker expressions are supported at the moment")
    }
  }
  scores <- apply(scores, 2, function(x) {
    1 - ((x - min(x, na.rm = TRUE))/ (max(x, na.rm = TRUE) - 
                                        min(x, na.rm = TRUE)))
  })
  nodeScores <- apply(scores, 1, mean)
  nodeScores[is.na(nodeScores)] <- 0
  cutoff <- max(nodeScores) * 0.95
  scoresCutoff <- nodeScores > cutoff
  return(list(selected = which(scoresCutoff), nodeScores = nodeScores, 
              fullScores = scores))
}

#' ParseSD 
#'  
#' Parses SD in FlowSOM object
#' 
#' Calculates the  standard deviaton of a FlowSOM object
#'
#' @param fsom    FlowSOM object, as generated by \code{\link{FlowSOM}}
#' @param marker  If a marker is given, the sd for this marker is shown.
#'                Otherwise, the maximum ratio is used.
#'                          
#' @return A vector containing the SDs
#' 
#' @seealso \code{\link{PlotFlowSOM}}, \code{\link{ParseEdges}}, 
#' \code{\link{ParseNodeSize}}, \code{\link{ParseArcs}},
#' \code{\link{ParseQuery}}, \code{\link{PlotSD}}
#' 
#' @importFrom stats median
#' 
ParseSD <-function (fsom, marker = NULL){
  stdevs <- fsom$map$sdValues[, fsom$map$colsUsed]
  stdev_medians <- apply(stdevs, 2, stats::median)
  if (is.null(marker)) {
    variable <- apply(stdevs, 1, function(x) {
      max(x/stdev_medians)
    })
  } else {
    variable <- stdevs[, marker]/stdev_medians[marker]
  }
  return(variable)
}

#' ScaleStarHeights 
#'
#' Scales starheights
#' 
#' Function that scales the star values between 0 and the node size
#' 
#' @param data              MedianValues of relevant markers extraced from 
#'                          FlowSOM object 
#' @param nodeSizes    A vector that is returned from \code{ParseNodeSize}
#'                          
#' @return A dataframe consisting of the scaled values of the stars. 
#'         The stars are scaled between 0 and the maximum of all stars
#' 
#' @seealso \code{\link{PlotFlowSOM}}, \code{\link{ParseNodeSize}}, 
#' \code{\link{AutoMaxNodeSize}}
#' 
ScaleStarHeights <- function(data, nodeSizes){
  nNodes <- length(nodeSizes)
  maxAllNodes <- max(data)
  minAllNodes <- min(data)
  for (i in 1:(nNodes)){
    maxRowNodeSize <- nodeSizes[i]
    scaledRow <- ((data[i, ] - minAllNodes) * maxRowNodeSize) / 
      (maxAllNodes - minAllNodes)
    data[i, ] <- scaledRow
  }
  return(data)
}

#' AutoMaxNodeSize 
#' 
#' Calculate node size
#' 
#' Function that calculates the minimum distance between the nodes
#' to use this to adapt the maxNodeSize for better plotting
#' 
#' @param layout            Coordinates of nodes
#' @param overlap           Parameter that determines how much overlap there 
#'                          will be. If negative the nodes will be smaller
#'                          
#' @return Returns the maxNodeSize with some overlap
#' 
#' @seealso \code{\link{PlotFlowSOM}}, \code{\link{ScaleStarHeights}}, 
#' \code{\link{ParseNodeSize}}
#' 
#' @importFrom stats dist
#' 
AutoMaxNodeSize <- function(layout, overlap){
  overlap <- 1 + overlap
  minDistance <- min(stats::dist(layout[, c(1, 2)]))
  return(minDistance / 2 * overlap)
}

#' UpdateFlowSOM 
#' 
#' Update old FlowSOM object to a new one and checks if it is a flowSOM object
#' 
#' Determines whether or not the fsom input is of class "FlowSOM" and returns  
#' the FlowSOM object and metaclustering object inside fsom
#' 
#' @param fsom  FlowSOM object, as generated by \code{\link{BuildMST}} or
#'              \code{\link{FlowSOM}}
#'                          
#' @return A FlowSOM object
#' 
#' @seealso \code{\link{PlotFlowSOM}}
#' 
#' @importFrom dplyr group_by summarise_all select
#' @importFrom stats  median
#' 
#' @export
UpdateFlowSOM <- function(fsom){
  if (is(fsom,"list") && !is.null(fsom$FlowSOM)) {
    fsom$FlowSOM$metaclustering <- fsom$metaclustering
    fsom <- fsom$FlowSOM
  }
  if (!is(fsom,"FlowSOM")) {
    stop("fsom should be a FlowSOM object.")
  }
  fsom$prettyColnames <- gsub("\\((.*)\\)", "<\\1>", fsom$prettyColnames)
  if (is.null(fsom$map$pctgs)){
    pctgs <- rep(0, fsom$map$nNodes)
    names(pctgs) <- as.character(seq_len(fsom$map$nNodes))
    pctgs_tmp <- table(fsom$map$mapping[, 1]) / nrow(fsom$map$mapping)
    pctgs[names(pctgs_tmp)] <- pctgs_tmp
    fsom$map$pctgs <- pctgs
  } 
  if (is.null(fsom$map$nMetaclusters)){
    fsom$map$nMetaclusters <- length(levels(fsom$metaclustering))
  }
  if (is.null(fsom$map$metaclusterMFIs) && !is.null(fsom$metaclustering)){
    fsom$map$metaclusterMFIs <- 
      data.frame(fsom$data, 
                 mcl = fsom$metaclustering[fsom$map$mapping[, 1]],
                 check.names = FALSE) %>% 
      dplyr::group_by(.data$mcl, .drop = FALSE) %>% 
      dplyr::summarise_all(stats::median) %>% 
      dplyr::select(-.data$mcl) %>% 
      data.frame(row.names = levels(fsom$metaclustering),
                 check.names = FALSE)
  }
  return(fsom)
}

#' AddMST 
#' 
#' Function plots the MST
#' 
#' @param p     GGplot object
#' @param fsom  FlowSOM object, as generated by \code{\link{FlowSOM}}
#'                          
#' @return Returns nothing, but plots the MST for FlowSOM MST view
#' 
#' @seealso \code{\link{PlotFlowSOM}}, \code{\link{ParseEdges}}, 
#' \code{\link{AddStarsPies}}, \code{\link{AddLabels}}, \code{\link{AddNodes}},
#' \code{\link{AddBackground}}, \code{\link{AddPies}}, \code{\link{AddStars}}
#' 
#' @import ggplot2
#' 
#' @export
AddMST <- function(p, 
                   fsom){
  requireNamespace("ggplot2")
  fsom <- UpdateFlowSOM(fsom)
  edges <- ParseEdges(fsom)
  p <- p + ggplot2::geom_segment(data = edges, 
                                 ggplot2::aes(x = .data$x, 
                                              y = .data$y, 
                                              xend = .data$xend, 
                                              yend = .data$yend),
                                 size = 0.2)
  return(p)
}

#' AddStarsPies 
#' 
#' Function plots stars or pies 
#' 
#' @param p            ggplot object
#' @param arcs         Dataframe that contains all the data for the plotting 
#'                     the pies or stars
#' @param colorPalette A vector of colors or a color function
#' @param showLegend   Boolean on whether to show the legend
#'                          
#' @return Returns nothing, but plots the stars or pies
#' 
#' @seealso \code{\link{PlotFlowSOM}}, \code{\link{AddLabels}},
#' \code{\link{AddNodes}}, \code{\link{AddBackground}}, \code{\link{AddPies}},
#' \code{\link{AddStars}}, \code{\link{ParseArcs}}, \code{\link{PlotStars}}
#' \code{\link{PlotPies}}
#' 
#' @import ggplot2
#' @importFrom ggforce geom_arc_bar
#' 
#' @export
AddStarsPies <- function(p, arcs, colorPalette, showLegend = TRUE){
  requireNamespace("ggplot2")
  
  p <- AddScale(p = p, 
                values = arcs$Markers, 
                colors = colorPalette, 
                showLegend = showLegend)
  
  p <- p + ggforce::geom_arc_bar(data = arcs, 
                                 ggplot2::aes(x0 = .data$x0, 
                                              y0 = .data$y0, 
                                              r0 = 0, 
                                              r = .data$value, 
                                              start = .data$start, 
                                              end = .data$end, 
                                              fill = .data$Markers),
                                 size = 0.2)
  
  
  return(p)
}

#' AddLabels 
#' 
#' @param p      ggplot object
#' @param labels Labels to be added to each node
#' @param hjust  Horizontal adjust for labels. Default is centered.
#' @param layout Dataframe with x and y columns. If null, the dataframe from
#'               the ggplot object will be reused.
#' @param textSize  Size for geom_text. Default (=3.88) is from geom_text.  
#' @param textColor Color for geom_text. Default = black.   
#' @param ...    Additional parameters to pass to geom_text
#'                          
#' @return Returns the ggplot object with labels added
#' 
#' @seealso \code{\link{PlotLabels}}, \code{\link{PlotNumbers}}
#' 
#' @import ggplot2
#' 
#' @export
AddLabels <- function(p, labels, hjust = 0.5, layout = NULL, 
                      textSize = 3.88, textColor = "black", ...){
  requireNamespace("ggplot2")
  
  if(is.null(layout)) layout <- ggplot2::ggplot_build(p)$plot$data
  
  p <- p + ggplot2::geom_text(data = layout,
                              ggplot2::aes(x = .data$x, 
                                           y = .data$y, 
                                           label = labels),
                              fontface = "bold", 
                              hjust = hjust,
                              size = textSize,
                              col = textColor,
                              ...)
  return(p)
}

#' AddNodes 
#' 
#' Function plots the nodes
#' 
#' @param p            GGplot object
#' @param nodeInfo    Dataframe with for every node an x, y and size value,
#'                     if null the dataframe from the ggplot object will be 
#'                     reused.
#' @param values       Values used for coloring the nodes. Default = NULL,
#'                     in which case all nodes are filled in fillColor.
#' @param lim          The limits of the color scale, not used if values = NULL.
#' @param colorPalette Colorpalette for color in nodes, not used if values = 
#'                     NULL. A vector of colors or a color function.
#' @param fillColor    Fixed fill for node colors, default = white.
#' @param showLegend   Boolean, default = TRUE.
#' @param label        Title for the legend.
#' @param ...          Additional arguments to pass to geom_circle
#'                          
#' @return Returns nothing, but plots the nodes
#' 
#' @seealso \code{\link{PlotFlowSOM}}, \code{\link{PlotMarker}},
#' \code{\link{PlotVariable}}, \code{\link{AddLabels}}, 
#' \code{\link{AddBackground}}, \code{\link{AddPies}}, 
#' \code{\link{AddStars}}, \code{\link{AddStarsPies}}
#' 
#' @import ggplot2
#' @importFrom ggplot2 ggplot_build
#' @importFrom ggforce geom_circle
#' @importFrom grDevices colorRampPalette
#' 
#' @export
AddNodes <- function(p,
                     nodeInfo = NULL,
                     values = NULL, 
                     lim = NULL,
                     colorPalette = NULL, 
                     fillColor = "white",
                     showLegend = TRUE,
                     label = "",
                     ...){
  requireNamespace("ggplot2")
  
  if(is.null(nodeInfo)) nodeInfo <- ggplot2::ggplot_build(p)$plot$data
  
  if (is.null(values)){
    p <- p + ggforce::geom_circle(data = nodeInfo,
                                  ggplot2::aes(x0 = .data$x,
                                               y0 = .data$y,
                                               r =  .data$size),
                                  fill = fillColor,
                                  size = 0.2,
                                  ...) 
  } else {
    
    if(is.character(values) | is.logical(values)) values <- factor(values)
    
    p <- AddScale(p, 
                  values = values,
                  colors = colorPalette,
                  limits = lim,
                  labelLegend = label,
                  showLegend = showLegend)
    
    p <- p + ggforce::geom_circle(ggplot2::aes(x0 = .data$x, 
                                               y0 = .data$y, 
                                               r =  .data$size,
                                               fill = values),
                                  size = 0.2,
                                  ...)
  }
  return(p)
}

#' AddBackground 
#' 
#' Function plots the background
#' 
#' @param p                 GGplot object
#' @param backgroundValues  Vector of values to be plotted as background for 
#'                          the nodes
#' @param backgroundColors  Colorpalette for backgroundcolors 
#' @param backgroundLim     Background limits (can be used to ensure consistent
#'                          colorpalette between plots). If NULL (default), will
#'                          be automatically adapted to the data.
#'                          
#'                          
#' @return Returns nothing, but plots the background
#' 
#' @seealso \code{\link{PlotFlowSOM}}, \code{\link{AddLabels}},
#' \code{\link{AddNodes}}, \code{\link{AddPies}}, \code{\link{AddStars}}
#' 
#' @import ggplot2
#' @importFrom ggforce geom_circle
#' @importFrom ggnewscale new_scale
#' @importFrom grDevices colorRampPalette
#' 
#' @export
AddBackground <- function(p,
                          backgroundValues, 
                          backgroundColors = NULL, 
                          backgroundLim = NULL){
  requireNamespace("ggplot2")
  
  if(is.character(backgroundValues)) {
    backgroundValues <- factor(backgroundValues)
  }
  
  p <- AddScale(p,
                backgroundValues,
                backgroundColors,
                backgroundLim,
                labelLegend = "background")
  
  p <- p + ggforce::geom_circle(ggplot2::aes(x0 = .data$x, 
                                             y0 = .data$y, 
                                             r =  .data$bg_size,
                                             fill = backgroundValues),
                                col = NA, 
                                alpha = 0.4)
  
  return(p)
}

#' AddScale
#' 
#' @param p           ggplot object
#' @param values      Values used for the fill
#' @param colors      Colors to use (can be a vector or a function)
#' @param limits      Limits to use in the scale
#' @param showLegend  Boolean on whether to show the legend
#' @param labelLegend Label to show as title of the legend
#' @param type        fill (default) or color
#' 
#' @return ggplot object with scale added
AddScale <- function(p, 
                     values = NULL,
                     colors = NULL,
                     limits = NULL,
                     showLegend = TRUE,
                     labelLegend = "",
                     type = "fill"){
  
  p <- p + ggnewscale::new_scale(type)
  
  if(is.character(values) | is.logical(values)) values <- factor(values)
  
  if (!is.null(limits) &&  is.null(colors)){
    args <- list()
    args[[type]] <- limits
    p <- p + do.call(ggplot2::lims, args)
  }
  
  if(!is.null(colors)){
    # Discrete values
    if (is.factor(values)){ 
      
      # If a function: make into the right amount of values
      if(is.function(colors)){
        colors <- colors(length(levels(values)))
      }
      
      # Check color names vs backgroundValues
      colorNames <- names(colors)
      levels <- levels(values)
      
      # If no names available, use levels
      if(is.null(colorNames)){
        colorNames <- levels[seq_along(colors)]
      }
      
      # If any names missing, show warning
      notAvailable <- ! levels %in% colorNames
      if(any(notAvailable)){
        warning("You did not provide a color for ",
                paste(levels[notAvailable], collapse =", "))
        colors <- c(colors,
                    rep("#FFFFFF", sum(notAvailable)))
        names(colors) <- c(colorNames,
                           levels[notAvailable])
      }
      scale_function <- eval(parse(text = paste0("ggplot2::scale_",
                                                 type,
                                                 "_manual")))
      if(showLegend) {
        p <- p + scale_function(values = colors)
      } else {
        p <- p + scale_function(values = colors,
                                guide = "none")
      }
      
    } else if (is.numeric(values)){ # Continuous values
      
      if (is.function(colors)){
        colors <- colors(100)
      } else {
        colors <- grDevices::colorRampPalette(colors)(100)
      }
      
      scale_function <- eval(parse(text = paste0("ggplot2::scale_",
                                                 type,
                                                 "_gradientn")))
      if(showLegend) {
        p <- p + scale_function(colors = colors, 
                                limits = limits)
      } else {
        p <- p + scale_function(colors = colors, 
                                limits = limits,
                                guide = "none")
      }
    }
  }
  
  p <- p + ggplot2::labs(fill = labelLegend)
  return(p)
}

#' AddPies 
#' 
#' Function plots the pies
#' 
#' @param p             GGplot object
#' @param fsom          FlowSOM object, as generated by \code{\link{BuildMST}}
#' @param cellLabels   Array of factors indicating the cell labels
#' @param layout        Coordinates of nodes. Uses dataframe of the ggplot 
#'                      object if NULL.
#' @param colorPalette  ColorPalette to be used for colors. Can be either a 
#'                      function or an array specifying colors.
#'                          
#' @return ggplot object with the pies added
#'
#' @seealso \code{\link{PlotFlowSOM}}, \code{\link{AddLabels}},
#' \code{\link{AddNodes}}, \code{\link{AddBackground}}, \code{\link{PlotPies}},
#' \code{\link{AddStars}}, \code{\link{ParseArcs}}
#' 
#' @importFrom grDevices colorRampPalette
#' @importFrom dplyr filter
#' 
#' @export
AddPies <- function(p, 
                    fsom, 
                    cellLabels, 
                    layout = NULL, 
                    colorPalette = NULL){
  
  fsom <- UpdateFlowSOM(fsom)
  
  if (is.null(layout)) layout <- ggplot2::ggplot_build(p)$plot$data
  
  if (!is.factor(cellLabels)) cellLabels <- as.factor(cellLabels)
  
  nNodes <- NClusters(fsom)
  count <- table(factor(GetClusters(fsom), levels = seq_len(nNodes)), 
                 factor(cellLabels, levels = c(levels(cellLabels), "Empty")))
  count[rowSums(count) == 0, "Empty"] <- 1
  
  pieValues <- lapply(seq_len(nNodes), function(cl){
    nodeData <- ParseArcs(layout[cl, "x"], 
                          layout[cl, "y"], 
                          count[cl, ], 
                          layout[cl, "size"])
    return(nodeData)
  })
  pieValues <- do.call(rbind, pieValues)
  
  pieValues$Markers <- factor(pieValues$Markers, 
                              levels = c(levels(cellLabels), "Empty"))
  p <- AddStarsPies(p, pieValues, colorPalette)
  return(p)
}

#' AddStars 
#' 
#' Function plots the stars
#' 
#' @param p             GGplot object
#' @param fsom          FlowSOM object, as generated by \code{\link{BuildMST}}
#' @param markers       Determines which markers to plot. 
#'                      Default = "fsom$map$colsUsed"
#' @param colorPalette  ColorPalette to be used for colors. Can be either a 
#'                      function or an array specifying colors.
#'                          
#' @return ggplot object with the stars added
#' 
#' @seealso \code{\link{PlotFlowSOM}}, \code{\link{AddLabels}},
#' \code{\link{AddNodes}}, \code{\link{AddBackground}}, \code{\link{PlotStars}},
#' \code{\link{AddPies}}, \code{\link{ParseArcs}}
#' 
#' @importFrom grDevices colorRampPalette
#' 
#' @export
AddStars <- function(p, 
                     fsom, 
                     markers = fsom$map$colsUsed,
                     colorPalette = NULL){
  requireNamespace("ggplot2")
  
  markers <- GetChannels(fsom, markers)
  
  fsom <- UpdateFlowSOM(fsom)
  nNodes <- NClusters(fsom)
  nMarkers <- length(markers)
  isEmpty <- which(fsom$map$pctgs == 0)
  nodeInfo <- ggplot2::ggplot_build(p)$plot$data
  nodeInfo$size[isEmpty] <- 0
  l1 <- NULL
  data <- fsom$map$medianValues[, markers, drop = FALSE]
  data[is.na(data)] <- 0
  data <- ScaleStarHeights(data, 
                           nodeInfo$size)
  markers_tmp <- rep(1, ncol(data))
  names(markers_tmp) <- colnames(data)
  starValues <- lapply(seq_len(nNodes), function(cl){
    nodeData <- ParseArcs(x = nodeInfo$x[cl], 
                          y = nodeInfo$y[cl], 
                          arcValues = markers_tmp, 
                          arcHeights = data[cl, ])
    return(nodeData)
  })
  starValues <- do.call(rbind, starValues)
  starValues$Markers <- factor(starValues$Markers, 
                               levels = colnames(data))
  p <- AddStarsPies(p, starValues, colorPalette, showLegend = FALSE)
  
  return(p)
}

#' PlotManualBars 
#' 
#' Function to plot the manual labels per FlowSOM (meta)cluster in a barplot
#' 
#' @param fsom          FlowSOM object, as generated by \code{\link{FlowSOM}} or 
#'                      by \code{\link{NewData}}. The clusters and metaclsuters 
#'                      will be plotted in the order of the factor levels.
#' @param fcs           Fcs file that should be mapped on the FlowSOM object.
#'                      Default is NULL. 
#' @param manualVector Vector with cell labels, e.g. obtained by manual gating
#' @param manualOrder  Optional vector with unique cell labels to fix in which
#'                      order the cell labels should be shown
#' @param colors        Optional color vector, should have the same length as
#'                      the number of unque cell labels 
#' @param list_insteadof_plots If \code{FALSE} (default), it returns multiple 
#'                             plots. If \code{TRUE}, it returns a list of 
#'                             ggplot objects
#'                          
#' @return Either a plot or a ggplot objects list is returned.
#' 
#' @examples 
#' # Identify the files
#' fcs_file <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' wsp_file <- system.file("extdata", "gating.wsp", package = "FlowSOM")
#' 
#' # Specify the cell types of interest for assigning one label per cell
#' cellTypes <- c("B cells",
#'                 "gd T cells", "CD4 T cells", "CD8 T cells",
#'                 "NK cells", "NK T cells")
#'                 
#' # Parse the FlowJo workspace   
#' library(flowWorkspace)             
#' gatingResult <- GetFlowJoLabels(fcs_file, wsp_file,
#'                                 cellTypes = cellTypes)               
#' 
#' # Build a FlowSOM object
#' flowSOM.res <- FlowSOM(fcs_file,
#'                        scale = TRUE,
#'                        compensate = TRUE, 
#'                        transform = TRUE,
#'                        toTransform = 8:18, 
#'                        colsToUse = c(9, 12, 14:18),
#'                        nClus = 10,
#'                        seed = 1)
#' 
#' # Make the barplot of the manual labels
#' pdf("PlotManualBars.pdf")
#' PlotManualBars(fsom = flowSOM.res,
#'                fcs = fcs_file,
#'                manualVector = gatingResult$manual,
#'                manualOrder = c(cellTypes, "Unlabeled"),
#'                colors = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", 
#'                           "#619CFF", "#F564E3", "#D3D3D3"))
#' dev.off()
#' 
#' @import     ggplot2
#' 
#' @export
#' 
PlotManualBars <- function(fsom, fcs = NULL, 
                           manualVector, manualOrder = NULL, 
                           colors = NULL, list_insteadof_plots = FALSE){
  #----Warnings----
  if (!is.null(manualOrder) & 
      !length(manualOrder) == length(unique(manualVector))){
    stop(paste0("Length of the manualOrder vector (", length(manualOrder), 
                ") should have the same length as the number of unique labels",
                "in manualVector (", length(unique(manualVector)), ")."))
  } 
  if (!is.null(colors) & !length(colors) >= length(unique(manualVector))){
    stop(paste0("Length of the colors vector (", length(colors),") should have",
                "at least the same length as the number of unique labels in ",
                "manualVector (", length(unique(manualVector)), ")."))
  } 
  
  #----Map fcs on fsom----
  C_counts <- matrix(0,
                     nrow = 1,
                     ncol = fsom$map$nNodes,
                     dimnames = 
                       list("fcs",
                            paste0("C", seq_len(fsom$map$nNodes))))
  if (!is.null(fcs)){
    fsom_tmp <- NewData(fsom, fcs)
  } else {
    fsom_tmp <- fsom
  }
  
  cluster_table <- table(GetClusters(fsom_tmp))
  
  C_counts["fcs", paste0("C", names(cluster_table))] <- cluster_table
  C_perc <- 100 * prop.table(C_counts, margin = 1)
  MC_counts <- t(apply(C_counts,
                       1,
                       function(x){
                         tapply(x, fsom$metaclustering, sum)
                       }))
  colnames(MC_counts) <- paste0("MC", colnames(MC_counts))
  MC_perc <- 100 * prop.table(MC_counts, margin  = 1)
  
  df <- data.frame("Manual" = manualVector, 
                   "MC" = GetMetaclusters(fsom_tmp),
                   "C" = GetClusters(fsom_tmp), 
                   stringsAsFactors = FALSE)
  if(is.null(manualOrder)) {
    manualOrder <- unique(manualVector)
  }
  
  #----Relative barplots MC----
  df_s <- data.frame(table(df[, 1:2]))
  p1 <- ggplot2::ggplot(data = 
                          transform(df_s,
                                    MC = factor(df_s$MC,
                                                levels = levels(fsom$metaclustering)),
                                    Manual = factor(df_s$Manual,
                                                    levels = manualOrder)),
                        ggplot2::aes(fill = .data$Manual, 
                                     y = .data$Freq, 
                                     x = .data$MC)) +
    ggplot2::geom_bar(position = "stack", stat = "identity") +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("Manual labels per MC") +
    ggplot2::ylab("size") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank()) 
  
  if(!is.null(colors)) {
    p1 <- p1 + ggplot2::scale_fill_manual(values = colors)
  }
  
  #----MC composition----
  df_s <- data.frame(table(df[, 1:2]))
  for (mc in as.character(seq_len(NMetaclusters(fsom)))){
    df_s$Freq[df_s$MC == mc] <- 
      df_s$Freq[df_s$MC == mc] / sum(df_s$Freq[df_s$MC == mc])
  }
  p2 <- ggplot2::ggplot(data = 
                          transform(df_s,
                                    MC = factor(df_s$MC,
                                                levels = levels(fsom$metaclustering)),
                                    Manual = factor(df_s$Manual,
                                                    levels = manualOrder)), 
                        ggplot2::aes(fill = .data$Manual, 
                                     y = .data$Freq, 
                                     x = .data$MC)) +
    ggplot2::geom_bar(position = "stack", stat = "identity") +
    ggplot2::theme_minimal() +
    ggplot2::ylab("") +
    ggplot2::ggtitle("Manual labels per MC") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank())
  if(!is.null(colors)) {
    p2 <- p2 + ggplot2::scale_fill_manual(values = colors)
  }
  
  #----Relative barplots C----
  df_s <- data.frame(table(df[, c(1, 3)]))
  df_s$MC <- fsom$metaclustering[as.numeric(levels(df_s$C))[df_s$C]]
  p3 <- ggplot2::ggplot(data =
                          transform(df_s,
                                    C = factor(df_s$C,
                                               levels = seq_len(NClusters(fsom))),
                                    MC = factor(df_s$MC,
                                                levels = levels(fsom$metaclustering)),
                                    Manual = factor(df_s$Manual,
                                                    levels = manualOrder)),
                        ggplot2::aes(fill = .data$Manual, 
                                     y = .data$Freq, 
                                     x = .data$C)) +
    ggplot2::geom_bar(position = "stack", stat = "identity") +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("Manual labels per C") +
    ggplot2::facet_wrap(~.data$MC, scales = "free_x") +
    ggplot2::ylab("size") +
    ggplot2::theme(axis.text.y = ggplot2:: element_blank(), 
                   strip.text = ggplot2::element_text(colour = "grey"))
  if(!is.null(colors)) {
    p3 <- p3 + ggplot2::scale_fill_manual(values = colors)
  }
  
  #----C composition----
  df_s <- data.frame(table(df[, c(1, 3)]))
  for (c in seq_len(NClusters(fsom))){
    df_s$Freq[df_s$C == c] <- df_s$Freq[df_s$C == c]/sum(df_s$Freq[df_s$C == c])
  }
  df_s$MC <- fsom$metaclustering[as.numeric(levels(df_s$C))[df_s$C]]
  p4 <- ggplot2::ggplot(data = 
                          transform(df_s,
                                    C = factor(df_s$C,
                                               levels = seq_len(
                                                 NClusters(fsom))),
                                    MC = factor(df_s$MC,
                                                levels = levels(fsom$metaclustering)),
                                    Manual = factor(df_s$Manual,
                                                    levels = manualOrder)), 
                        ggplot2::aes(fill = .data$Manual, 
                                     y = .data$Freq, 
                                     x = .data$C)) +
    ggplot2::geom_bar(position = "stack", stat = "identity") +
    ggplot2::theme_minimal() +
    ggplot2::ylab("") +
    ggplot2::ggtitle("Main labels per C") +
    ggplot2::facet_wrap(~.data$MC, scales = "free_x") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank())
  if(!is.null(colors)) {
    p4 <- p4 + ggplot2::scale_fill_manual(values = colors)
  }
  
  #----Print plot or return list----
  if (!list_insteadof_plots) {
    print(p1)
    print(p2)
    print(p3)
    print(p4)
  } else {
    return(list(p1, p2, p3, p4))
  }
}

#' gg_color_hue 
#' 
#' Helper function to get the ggplot colors
#' @param n Number of colors
#' 
#' @return array with hexadecimal color values
#' 
#' @importFrom grDevices hcl
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Plot2DScatters 
#' 
#' Function to draw 2D scatter plots of FlowSOM (meta)clusters
#' 
#' Plot multiple 2D scatter plots in a png file. A subset of fsom$data is 
#' plotted in grey, and those of the selected clusters and metaclusters are 
#' plotted in color.
#' 
#' @param fsom              FlowSOM object, as created by \code{\link{FlowSOM}}
#' @param channelpairs      List in which each element is a pair of channel
#'                          or marker names
#' @param clusters          Vector or list (to combine multiple clusters in one
#'                          plot) with clusters of interest
#' @param metaclusters      Vector or list (to combine multiple metaclusters in 
#'                          one plot) with metaclusters of interest
#' @param maxBgPoints     Maximum number of background cells to plot
#' @param sizeBgPoints    Size of the background cells
#' @param maxPoints        Maximum number of (meta)cluster cells to plot
#' @param sizePoints       Size of the (meta)cluster cells
#' @param xLim             A vector of a lower and upper limit of the x-axis
#' @param yLim             A vector of a lower and upper limit of the y-axis
#' @param density           Default is \code{TRUE} to color the (meta)cluster 
#'                          points according to density. Set to \code{FALSE} to 
#'                          use a plain color
#' @param centers           Default is \code{TRUE} to show the cluster centers
#' @param color             Colors for all the cells in the selected nodes 
#'                          (ordered list). If  \code{NULL} the default ggplot 
#'                          colors, indexed by metacluster number are used.
#' @param plotFile          If a filepath for a png is given (default = 
#'                          2DScatterPlots.png), the plots will be plotted in 
#'                          the corresponding png file. If \code{NULL}, a list 
#'                          of ggplot objects will be returned
#'                          
#'                          
#' @return If \code{plot} is \code{TRUE}, nothing is returned and a plot is 
#'         drawn in which background cells are plotted in grey and the cells of
#'         the selected nodes in color. If \code{plot} is \code{FALSE}, a ggplot
#'          objects list is returned.
#' 
#' 
#' @examples 
#' # Identify the files
#' fcs <- flowCore::read.FCS(system.file("extdata", "68983.fcs", 
#'                                       package = "FlowSOM"))
#' 
#' # Build a FlowSOM object
#' flowSOM.res <- FlowSOM(fcs, 
#'                        scale = TRUE,
#'                        compensate = TRUE, 
#'                        transform = TRUE,
#'                        toTransform = 8:18, 
#'                        colsToUse = c(9, 12, 14:18),
#'                        nClus = 10,
#'                        seed = 1)
#' 
#' # Make the 2D scatter plots of the clusters and metaclusters of interest
#' Plot2DScatters(fsom = flowSOM.res,
#'                channelpairs = list(c("PE-Cy7-A", "PE-Cy5-A"),
#'                                    c("PE-Texas Red-A", "Pacific Blue-A")),
#'                clusters = c(1, 48, 49, 82, 95),
#'                metaclusters = list(c(1, 4), 9),
#'                density = FALSE)
#'                
#' Plot2DScatters(fsom = flowSOM.res,
#'                channelpairs = list(c("PE-Texas Red-A", "Pacific Blue-A")),
#'                metaclusters = list(c(1, 4)),
#'                density = FALSE,
#'                color = list(c("red", "green")))
#' 
#' @import     ggplot2
#' @importFrom colorRamps matlab.like2
#' @importFrom ggpointdensity geom_pointdensity
#' @importFrom grDevices png dev.off
#' @importFrom ggpubr ggarrange
#' 
#' @export
#' 
Plot2DScatters <- function(fsom, 
                           channelpairs, 
                           clusters = NULL, 
                           metaclusters = NULL, 
                           maxBgPoints = 3000, 
                           sizeBgPoints = 0.5,
                           maxPoints = 1000, 
                           sizePoints = 0.5,
                           xLim = NULL,
                           yLim = NULL,
                           density = TRUE, 
                           centers = TRUE, 
                           color = NULL,
                           plotFile = "2DScatterPlots.png"){
  if(!is.null(fsom$metaclustering)){
    metacluster <- as.numeric(fsom$metaclustering)
    nMetaclusters <- NMetaclusters(fsom)
  } else {
    metacluster <- rep(1, NClusters(fsom))
    nMetaclusters <- 1
  }
  
  medianValues <- fsom$map$medianValues
  cellCluster <- GetClusters(fsom)
  
  #----Warnings----
  if (density == FALSE & 
      !is.null(color) & 
      !((length(clusters) + length(metaclusters)) == length(color))){
    stop("Length of color list should be equal to the joined length ",
         "of the clusters and metaclusters")
  } 
  
  #----Join the clusters and metaclusters of interest----
  i <- sample(nrow(fsom$data), 
              min(nrow(fsom$data), maxBgPoints))
  # loop over all subsets at once
  subsets <- list("Cluster" = as.list(clusters),
                  "Metacluster" = as.list(metaclusters)) 
  plots_list <- list()
  color_n <- 0
  for (group in names(subsets)){
    for (subset in subsets[[group]]){
      color_n <- color_n + 1
      n <- as.numeric(unlist(subset))
      for (channelpair in channelpairs){
        
        channelpair <- GetChannels(fsom, channelpair)
        #----background dataframe----
        df_bg <- data.frame(fsom$data[i, c(channelpair[1], channelpair[2])]) 
        colnames(df_bg) <- c("m1", "m2")
        
        #----(meta)cluster dataframe----
        if(group == "Cluster") {
          # dataframe with cluster's points
          df_ss <- data.frame(fsom$data[cellCluster %in% n,
                                        c(channelpair[1], channelpair[2]),
                                        drop = FALSE],
                              "Population" = cellCluster[cellCluster %in% n])
          # dataframe with centers
          df_c <- data.frame(medianValues[n, 
                                          c(channelpair[1], channelpair[2]),
                                          drop = FALSE])
          col <- gg_color_hue(nMetaclusters)[metacluster[n]]
        } else {
          df_ss <- data.frame(fsom$data[which(metacluster[cellCluster] %in% n), 
                                        c(channelpair[1], channelpair[2]),
                                        drop = FALSE],
                              "Population" =
                                metacluster[cellCluster[ 
                                  which(metacluster[cellCluster] %in% n)]])
          df_c <- data.frame(medianValues[which(metacluster
                                                [1:nrow(medianValues)] %in% n), 
                                          c(channelpair[1], channelpair[2]),
                                          drop = FALSE])
          col <- gg_color_hue(nMetaclusters)[n]
        }
        df_ss <- data.frame(df_ss[sample(nrow(df_ss), 
                                         min(nrow(df_ss), maxPoints)), ])
        df_ss$Population <- factor(df_ss$Population, levels = subset)
        colnames(df_ss)[1:2] <- c("m1", "m2")
        colnames(df_c) <- c("m1", "m2")
        
        #----ggplots----
        s <- ""
        if (length(subset) > 1){
          s <- "s"
        }
        
        p <- ggplot2::ggplot(data = df_ss, 
                             ggplot2::aes(x = .data$m1, 
                                          y = .data$m2)) +
          ggplot2::geom_point(data = df_bg, 
                              colour = "grey", 
                              size = sizeBgPoints) + # background dot plot
          ggplot2::theme_classic() +
          ggplot2::ggtitle(paste0(group, s, " ", paste(ifelse(group == "Cluster", 
                                                              subset,
                                                              levels(fsom$metaclustering)[unlist(subset)]), 
                                                       collapse = ", "))) +
          ggplot2::xlab(GetMarkers(fsom, channelpair[1])) +
          ggplot2::ylab(GetMarkers(fsom, channelpair[2])) + 
          ggplot2::theme(legend.position = "none")
        if (!is.null(xLim)) p <- p + ggplot2::xlim(xLim)
        if (!is.null(yLim)) p <- p + ggplot2::ylim(yLim)
        
        if (density) {
          # subset density plot
          p <- p  + ggpointdensity::geom_pointdensity(size = sizePoints) + 
            ggplot2::scale_color_gradientn(colors = 
                                             colorRamps::matlab.like2(10))
        } else {
          # if no colors are given, the default colors of ggplot are used
          if (is.null(color)){ 
            p <- p + ggplot2::geom_point(ggplot2::aes(color = .data$Population),
                                         size = sizePoints) +
              ggplot2::scale_color_manual(values = col) 
          } else {
            #subset plot, metacluster colors
            p <- p + ggplot2::geom_point(ggplot2::aes(color = .data$Population),
                                         size = sizePoints) +
              ggplot2::scale_color_manual(values=color[[color_n]])
            # color_n <- color_n + 1
          }
          
        }
        
        #----add cluster centers----
        if (centers) { 
          p <- p + ggplot2::geom_point(data = df_c, 
                                       shape = 21, 
                                       fill = "white", 
                                       color = "black", 
                                       size = 3)
        }
        plots_list[[length(plots_list) + 1]] <- p
      }
    }
  }
  if (!is.null(plotFile)) {
    grDevices::png(plotFile, 
                   width = 400 * length(channelpairs), 
                   height = 400 * (length(clusters) + length(metaclusters)))
    print(ggpubr::ggarrange(plotlist = plots_list, 
                            ncol = length(channelpairs), 
                            nrow = (length(clusters) + length(metaclusters)),
                            common.legend = FALSE))
    grDevices::dev.off()
  } else {
    return(plots_list)
  }
} 

#' FlowSOMmary 
#' 
#' This functions plots a summary of a flowSOM object. It includes a table of 
#' (meta)cluster data, the flowSOM trees and grid view, the (meta)cluster 
#' labels, the markers expression, the file distribution if present,
#' the cluster per metacluster percentage, a t-SNE plot,
#' and the MFI per metacluster.
#' 
#' @param fsom          FlowSOM object, as generated by \code{\link{FlowSOM}}
#' @param plotFile      Name of the pdf file that will be generated (default is 
#'                      FlowSOMmary.pdf). If \code{NULL}, a list of ggplots will 
#'                      be returned.
#'                          
#' @return Returns a summary of the FlowSOM object
#' 
#' 
#' @examples 
#' # Identify the files
#' fcs <- flowCore::read.FCS(system.file("extdata", "68983.fcs", 
#'                                       package = "FlowSOM"))
#' 
#' # Build a FlowSOM object
#' flowSOM.res <- FlowSOM(fcs, 
#'                        scale = TRUE,
#'                        compensate = TRUE, 
#'                        transform = TRUE,
#'                        toTransform = 8:18, 
#'                        colsToUse = c(9, 12, 14:18),
#'                        nClus = 10,
#'                        seed = 1)
#'                        
#' FlowSOMmary(flowSOM.res)
#' 
#' @import ggplot2
#' @importFrom ggpubr ggarrange ttheme ggtexttable
#' @importFrom dplyr count mutate group_by filter pull arrange
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices pdf
#' @importFrom scattermore geom_scattermore
#' 
#' @export
FlowSOMmary <- function(fsom, plotFile = "FlowSOMmary.pdf"){
  
  #----Initializing ----
  fsom <- UpdateFlowSOM(fsom)
  metaclustersPresent <- !is.null(fsom$metaclustering)
  if (metaclustersPresent){
    mfis <- GetMetaclusterMFIs(fsom, colsUsed = TRUE)
    metaclusters <- GetMetaclusters(fsom)
    nMetaclusters <- NMetaclusters(fsom)
  }
  clusters <- GetClusters(fsom)
  nNodes <- seq_len(NClusters(fsom))
  filePresent <- "File" %in% colnames(fsom$data)
  plotList <- list()
  
  #----Plot fsom trees and grids----
  message("Plot FlowSOM trees")
  
  for(view in c("MST", "grid")){
    
    plotList[[paste0("stars_", view)]] <- 
      PlotStars(fsom, view = view, backgroundValues = fsom$metaclustering, 
                title = paste0("FlowSOM ", view))
    
    if (metaclustersPresent){
      p2.1 <- PlotFlowSOM(fsom, view = view, title = "FlowSOM Clusters",
                          equalNodeSize = TRUE) %>% 
        AddNodes(values = fsom$metaclustering, 
                 showLegend = TRUE,
                 label = "Metaclusters") %>% 
        AddLabels(labels = seq_len(NClusters(fsom)))
      
      p2.2 <- PlotFlowSOM(fsom, view = view, equalNodeSize = TRUE,
                          title = "FlowSOM Metaclusters") %>% 
        AddNodes(values = fsom$metaclustering, showLegend = TRUE,
                 label = "Metaclusters") %>% 
        AddLabels(labels = as.numeric(fsom$metaclustering))
      
    } else {
      p2.1 <- PlotNumbers(fsom, view = view, title = "FlowSOM Clusters", 
                          maxNodeSize = "auto", equalNodeSize = TRUE)
      p2.2 <- NULL
    }
    plotList[[paste0("labels_",view)]] <- 
      ggpubr::ggarrange(p2.1, p2.2, 
                        common.legend = TRUE, legend = "right")
  }
  
  #----Plot Markers----
  plotList[["p5"]] <- PlotMarker(fsom, marker = fsom$map$colsUsed, 
                                 refMarkers = fsom$map$colsUsed, 
                                 equalNodeSize = TRUE)
  
  #----File distribution----
  if (filePresent){
    message("Plot file distribution")
    p6 <- PlotPies(fsom, cellTypes = factor(fsom$data[, "File"]), 
                   equalNodeSize = TRUE, view = "grid", 
                   title = "File distribution per cluster", 
                   colorPalette = FlowSOM_colors)
    
    filesI <- as.character(unique(fsom$data[, "File"]))
    expectedDistr <- rep(1, length(filesI))
    names(expectedDistr) <- filesI
    arcsDf <- ParseArcs(0, 0, expectedDistr, 0.7)
    arcsDf$Markers <- factor(arcsDf$Markers, levels = filesI)
    plotList[["p6"]] <- AddStarsPies(p6, arcsDf, colorPalette = FlowSOM_colors(
      length(filesI) + 1), showLegend = FALSE) +
      ggplot2::geom_text(ggplot2::aes(x = 0, y = -1, 
                                      label = "Expected distribution"))
  }
  
  #----t-SNE----
  message("Calculate t-SNE")
  dimred_res <- PlotDimRed(fsom, cTotal = 5000, colorBy = fsom$map$colsUsed,
                           seed = 1, returnLayout = TRUE,
                           title = paste0("t-SNE with markers used in FlowSOM ",
                                          "call (perplexity = 30, cells = 5000)"))
  if (metaclustersPresent){
    plotList[["p7"]] <- PlotDimRed(fsom, dimred = dimred_res$layout, seed = 1,
                                   title = paste0("t-SNE with markers used in FlowSOM ",
                                                  "call (perplexity = 30, cells = 5000)"))
  }
  plotList[["p8"]] <- dimred_res$plot
  
  #----Cluster Per Metacluster----
  if (metaclustersPresent){
    
    message("Plot cluster per metacluster distibution")
    
    clusterPerMetacluster <- data.frame(metaclusters, clusters) %>%
      dplyr::group_by(metaclusters) %>%
      dplyr::count(clusters) %>%
      dplyr::mutate(Percentage = .data$n / sum(.data$n) * 100) %>%
      dplyr::mutate(LabelPos = cumsum(.data$Percentage) - 1) %>%
      as.data.frame()
    
    clusterPerMetacluster$clusters <- factor(clusterPerMetacluster$cluster)
    plotList[["p9"]] <- ggplot2::ggplot(clusterPerMetacluster,
                                        ggplot2::aes(x = metaclusters)) +
      ggplot2::geom_bar(ggplot2::aes(y = .data$Percentage, fill = metaclusters),
                        stat = "identity", col = "black") +
      ggplot2::geom_text(ggplot2::aes(y = .data$LabelPos, label = clusters)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.minor.x = ggplot2::element_blank(),
                     legend.position = "none") +
      ggplot2::ggtitle("Percentages clusters per metacluster")
    
    #----Median expression per metacluster----
    
    message("Plot heatmap")
    colnames(mfis) <- fsom$prettyColnames[colnames(mfis)]
    rownames(mfis) <- levels(metaclusters)
    mfis_scaled <- scale(mfis)
    mfis_scaled[is.na(mfis_scaled)] <- 0
    plotList[["empty"]] <- ggplot2::ggplot() + ggplot2::theme_minimal()
    plotList[["p10"]] <-
      pheatmap::pheatmap(mfis_scaled, scale = "none",
                         display_numbers = round(mfis, 2),
                         main = "Median expression per metacluster",
                         silent = TRUE)
  }
  
  #----Table----
  message("Make tables")
  
  datatable1 <- data.frame("Total number of cells" = nrow(fsom$data), 
                           check.names = FALSE)
  rownames(datatable1) <- "FlowSOMmary"
  if (metaclustersPresent) {
    datatable1[, "Total number of metaclusters"] <- nMetaclusters
  }

  markersInFlowSOM <- split(fsom$prettyColnames[fsom$map$colsUsed], 
                            rep(seq_len(ceiling(length(fsom$prettyColnames[
                              fsom$map$colsUsed])/5)), 
                                each = 5)[seq_len(length(fsom$prettyColnames[
                                  fsom$map$colsUsed]))]) %>% 
    sapply(paste, collapse =", ") %>% 
    paste(collapse = ",\n")
  
  datatable1 <- cbind(datatable1, 
                      "Total number of clusters" = fsom$map$nNodes,
                      "Markers used for FlowSOM" = markersInFlowSOM)
  
  datatable1 <- format(datatable1, digits = 2)
  t1 <- ggpubr::ggtexttable(t(datatable1), theme = ggpubr::ttheme("minimal"))
  
  if (metaclustersPresent){
    freqMetaclusters <- data.frame(metaclusters = 
                                     as.character(metaclusters)) %>%
      dplyr::count(.data$metaclusters) %>%
      dplyr::mutate(percentage = .data$n / sum(.data$n) * 100) %>%
      as.data.frame()
    datatable2 <- data.frame("Metacluster" = levels(metaclusters),
                             "Number of cells" = 0,
                             "Percentage of cells" = 0,
                             "Number of clusters" = 0,
                             "Clusters" = "",
                             check.names = FALSE)
    rownames(datatable2) <- datatable2$Metacluster
    datatable2[freqMetaclusters[,"metaclusters"], "Number of cells"] <- 
      freqMetaclusters[,"n"]
    datatable2[freqMetaclusters[,"metaclusters"], "Percentage of cells"] <- 
      freqMetaclusters[,"percentage"]
    datatable2[, "Number of clusters"] <- 
      sapply(levels(metaclusters), function(x){
        which(fsom$metaclustering == x) %>% length()})
    datatable2[, "Clusters"] <- 
      sapply(levels(metaclusters), function(x){
        cl_selected <- which(fsom$metaclustering == x)
        clusters_string <- split(cl_selected, 
                                 rep(seq_len(ceiling(length(cl_selected)/25)), 
                                     each = 25)[seq_len(length(cl_selected))]) %>% 
          sapply(paste, collapse =", ") %>% 
          paste(collapse = ",\n")
        return(clusters_string)
      })
    
    datatable2 <- format(datatable2, digits = 2)
    split_datatable2 <- split(datatable2, rep(seq_len(
      ceiling(nrow(datatable2) / 30)), each = 30, 
      length.out=nrow(datatable2)))
  }
  
  freqClusters <- as.data.frame(clusters) %>%
    dplyr::count(.data$clusters) %>%
    dplyr::mutate(freq = .data$n / sum(.data$n) * 100)
  
  datatable3 <- data.frame("Cluster" = factor(nNodes),
                           "Number of cells" = 0,
                           "Percentage of cells" = 0,
                           check.names = FALSE)
  datatable3[freqClusters[,"clusters"], "Number of cells"] <- 
    freqClusters[,"n"]
  datatable3[freqClusters[,"clusters"], "Percentage of cells"] <- 
    freqClusters[,"freq"]
  
  if (metaclustersPresent){
    datatable3 <-  cbind(datatable3, 
                         "Belongs to metacluster" = fsom$metaclustering,
                         "Percentage in metacluster" = 0)
    datatable3[as.numeric(as.character(clusterPerMetacluster[,"clusters"])), 
               "Percentage in metacluster"] <-
      clusterPerMetacluster[, "Percentage"]
  }
  
  datatable3 <- format(datatable3, digits = 2)
  split_datatable3 <- split(datatable3, rep(seq_len(
    ceiling(nrow(datatable3) / 30)), each = 30, 
    length.out = nrow(datatable3)))
  #----Printing----
  message("Printing")
  if (!is.null(plotFile)){
    grDevices::pdf(plotFile, width = 17, height = 10)
    print(t1)
    if (metaclustersPresent){
      for (table2 in split_datatable2) {
        print(ggpubr::ggtexttable(table2, theme = ggpubr::ttheme("minimal"), 
                                  rows = NULL))
      }
    }
    for (table3 in split_datatable3){
      print(ggpubr::ggtexttable(table3, theme = ggpubr::ttheme("minimal"), 
                                rows = NULL))
    }
    for (plot in plotList){
      print(plot)
    }
    dev.off()
  } else {
    return(plotList)
  }
}

#' AddAnnotation
#' 
#' Add annotation to a FlowSOM plot
#' 
#' @param p       Plot to add annotation to. When using type = "stars", please 
#'                use plot = FALSE in \code{\link{PlotFlowSOM}}.
#' @param fsom    FlowSOM object that goes with the plot
#' @param layout  Layout to use. Default fsom$MST$l
#' @param cl      Clusters to annotate. If cl = "all", all clusters will be 
#'                annotated
#' @param mcl     Metaclusters to annotate. If mcl = "all", all metaclusters 
#'                will be annotated
#' @param clCustomLabels Cluster labels to use for annotation
#' @param mclCustomLabels Metacluster labels to use for annotation
#' @param hjust   Label adjustment. Default = 0.5
#' @param ...     Arguments passed to geom_text_repel
#' 
#' @return The updated plot
#' 
#' @examples
#' 
#' #' # Identify the files
#' fcs <- flowCore::read.FCS(system.file("extdata", "68983.fcs", 
#'                                       package = "FlowSOM"))
#' # Build a FlowSOM object
#' flowSOM.res <- FlowSOM(fcs, 
#'                        scale = TRUE,
#'                        compensate = TRUE, 
#'                        transform = TRUE,
#'                        toTransform = 8:18, 
#'                        colsToUse = c(9, 12, 14:18),
#'                        nClus = 10,
#'                        seed = 1)
#'                        
#' p <- PlotStars(flowSOM.res, backgroundValues = flowSOM.res$metaclustering,
#'                list_insteadof_ggarrange = TRUE)
#' AddAnnotation(p, flowSOM.res, cl = c(1, 2), mcl = c(3, 4))
#' 
#' @importFrom ggrepel geom_text_repel
#' 
#' @export
AddAnnotation <- function(p, fsom, layout = fsom$MST$l, cl = NULL, mcl = NULL, 
                          clCustomLabels = NULL, mclCustomLabels = NULL, 
                          hjust = 0.5, ...){
  listOrGGplot <- "tree" %in% names(p)
  if (listOrGGplot) {
    l1 <- p[["starLegend"]]
    l2 <- p [["backgroundLegend"]]
    p <- p[["tree"]]
  }
  fsom <- UpdateFlowSOM(fsom)
  if (!is.data.frame(layout)) layout <- as.data.frame(layout)
  colnames(layout) <- c("V1", "V2")
  clusterLabels <- data.frame(layout, 
                              cl = as.character(seq_len(fsom$map$nNodes)), 
                              metacl = as.character(fsom$metaclustering))
  labels <- list()
  if (!is.null(mcl)){
    metaclusterLabels <- clusterLabels %>% 
      dplyr::group_by(.data$metacl) %>% 
      dplyr::filter(dplyr::row_number() == 1) %>% 
      as.data.frame()
    if (!"all" %in% mcl){
      metaclusterLabels <- metaclusterLabels[
        which(metaclusterLabels$metacl %in% mcl), ]
    }
    if (!is.null(mclCustomLabels)){
      metaclusterLabels$label <- mclCustomLabels
    } else {
      metaclusterLabels$label <- paste0("MCL", metaclusterLabels$metacl)
    }
    labels[["metaclusters"]] <- metaclusterLabels
  }
  if (!is.null(cl)){
    if (!"all" %in% cl){
      clusterLabels <- clusterLabels[which(clusterLabels$cl %in% cl), ]
    }
    if (!is.null(clCustomLabels)){
      clusterLabels$label <- clCustomLabels
    } else {
      clusterLabels$label <- paste0("CL", clusterLabels$cl)
    }
    labels[["clusters"]] <- clusterLabels
  }
  labels <- do.call(rbind, labels)
  p <- p + ggrepel::geom_text_repel(data = labels, 
                                    ggplot2::aes(x = .data$V1, y = .data$V2, 
                                                 label = .data$label), 
                                    segment.color = "grey", force = 20, 
                                    segment.size = 0.2, point.padding = 0.5,
                                    ...)
  if (listOrGGplot){
    p <- ggpubr::ggarrange(p, ggpubr::ggarrange(l1, l2, ncol = 1), NULL,
                           ncol = 3, widths = c(3, 1, 0.3), legend = "none")
  }
  return(p)
}
