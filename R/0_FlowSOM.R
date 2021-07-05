#' Run the FlowSOM algorithm
#'
#' Method to run general FlowSOM workflow. 
#' Will scale the data and uses consensus meta-clustering by default.
#'
#' @param input         a flowFrame, a flowSet or an array of paths to files or 
#'                      directories
#' @param pattern       if input is an array of file- or directorynames, select 
#'                      only files containing pattern
#' @param compensate    logical, does the data need to be compensated
#' @param spillover     spillover matrix to compensate with
#'                      If NULL and compensate = TRUE, we will look for $SPILL 
#'                      description in fcs file.
#' @param transform     logical, does the data need to be transformed with the
#'                      transformation given in \code{transformFunction}.
#' @param toTransform   column names or indices that need to be transformed.
#'                      Will be ignored if \code{transformList} is given.
#'                      If \code{NULL} and transform = \code{TRUE}, column names
#'                      of \code{$SPILL} description in fcs file will be used.
#' @param transformFunction Defaults to logicleTransform()
#' @param transformList transformList to apply on the samples.
#' @param scale         logical, does the data needs to be rescaled. 
#'                      Default = FALSE
#' @param scaled.center see \code{\link{scale}}
#' @param scaled.scale  see \code{\link{scale}}
#' @param silent        if \code{TRUE}, no progress updates will be printed
#' @param colsToUse     Markers, channels or indices to use for building the SOM. 
#'                      Default (NULL) is all the columns used to build the 
#'                      FlowSOM object.
#' @param importance    array with numeric values. Parameters will be scaled 
#'                      according to importance
#' @param nClus         Exact number of clusters for meta-clustering. 
#'                      Ignored if maxMeta is specified.
#'                      Default = 10.
#' @param maxMeta       Maximum number of clusters to try out for 
#'                      meta-clustering. If \code{NULL} (default), only one 
#'                      option will be computed (\code{nClus}).
#' @param seed          Set a seed for reproducible results
#' @param ...           options to pass on to the SOM function 
#'                      (xdim, ydim, rlen, mst, alpha, radius, init, distf)
#'
#' @return A \code{list} with two items: the first is the flowSOM object 
#'         containing all information (see the vignette for more detailed 
#'         information about this object), the second is the metaclustering of 
#'         the nodes of the grid. This is a wrapper function for 
#'         \code{\link{ReadInput}}, \code{\link{BuildSOM}}, 
#'         \code{\link{BuildMST}} and \code{\link{MetaClustering}}. 
#'         Executing them separately may provide more options.
#'
#' @seealso \code{\link{scale}}, 
#'          \code{\link{ReadInput}}, 
#'          \code{\link{BuildSOM}},
#'          \code{\link{BuildMST}}, 
#'          \code{\link{MetaClustering}}
#' @examples
#' # Read from file
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' flowSOM.res <- FlowSOM(fileName, compensate = TRUE, transform = TRUE,
#'                       scale = TRUE, colsToUse = c(9, 12, 14:18), nClus = 10)
#' # Or read from flowFrame object
#' ff <- flowCore::read.FCS(fileName)
#' ff <- flowCore::compensate(ff, flowCore::keyword(ff)[["SPILL"]])
#' ff <- flowCore::transform(ff,
#'          flowCore::transformList(colnames(flowCore::keyword(ff)[["SPILL"]]),
#'                                 flowCore::logicleTransform()))
#' flowSOM.res <- FlowSOM(ff, 
#'                        scale = TRUE, 
#'                        colsToUse = c(9, 12, 14:18), 
#'                        nClus = 10)
#' 
#' # Plot results
#' PlotStars(flowSOM.res,
#'           backgroundValues = flowSOM.res$metaclustering)
#' 
#' # Get metaclustering per cell
#' flowSOM.clustering <- GetMetaclusters(flowSOM.res)
#' 
#' 
#' 
#' @importFrom BiocGenerics colnames
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus
#' @importFrom flowCore read.FCS compensate transform logicleTransform exprs 
#'             transformList write.FCS 'exprs<-' keyword fr_append_cols
#' @importFrom flowWorkspace gs_get_pop_paths gh_pop_get_indices gh_pop_get_data
#'             gs_get_leaf_nodes
#' @importFrom CytoML open_flowjo_xml flowjo_to_gatingset
#' @importFrom igraph graph.adjacency minimum.spanning.tree layout.kamada.kawai
#'             plot.igraph add.vertex.shape get.edges shortest.paths E V 'V<-'
#'             igraph.shape.noclip
#' @importFrom RColorBrewer brewer.pal
#' @importFrom rlang .data
#' @importFrom stats prcomp
#' @importFrom utils capture.output
#' @importFrom XML xmlToList xmlParse
#' @importFrom dplyr group_by summarise_all select
#' @importFrom stats median
#' 
#' @export
FlowSOM <- function(input, pattern = ".fcs", 
                    compensate = FALSE, spillover = NULL, 
                    transform = FALSE, toTransform = NULL, 
                    transformFunction = flowCore::logicleTransform(), 
                    transformList = NULL, scale = FALSE, 
                    scaled.center = TRUE, scaled.scale = TRUE, silent = TRUE, 
                    colsToUse = NULL, nClus = 10, maxMeta = NULL, importance = NULL, 
                    seed = NULL, ...){
  # Method to run general FlowSOM workflow. 
  # Will scale the data and uses consensus meta-clustering by default.
  #
  # Args:
  #    input: dirName, fileName, array of fileNames, flowFrame or 
  #           array of flowFrames
  #    colsToUse: column names or indices to use for building the SOM
  #    maxMeta: maximum number of clusters for meta-clustering
  #
  # Returns:
  #    list with the FlowSOM object and an array with final clusterlabels
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  t <- system.time(fsom <- ReadInput(input, pattern = pattern, 
                                     compensate = compensate, 
                                     spillover = spillover, 
                                     transform = transform, 
                                     toTransform = toTransform, 
                                     transformFunction = transformFunction, 
                                     transformList = transformList,
                                     scale = scale,
                                     scaled.center = scaled.center, 
                                     scaled.scale = scaled.scale, 
                                     silent = silent))
  if(!silent) message(t[3], "\n")
  t <- system.time(fsom <- BuildSOM(fsom, colsToUse, silent = silent, 
                                    importance = importance, ...))
  if(!silent) message(t[3], "\n")
  t <- system.time(fsom <- BuildMST(fsom, silent = silent))
  if(!silent) message(t[3], "\n")
  if(is.null(maxMeta)){
    t <- system.time(cl <- as.factor(
      metaClustering_consensus(fsom$map$codes, nClus, seed = seed)))
  } else {
    t <- system.time(cl <- as.factor(MetaClustering(fsom$map$codes,
                                                    "metaClustering_consensus", 
                                                    maxMeta,
                                                    seed = seed)))
  }
  fsom$map$nMetaclusters <- length(levels(cl))
  fsom$metaclustering <- cl
  fsom <- UpdateDerivedValues(fsom)
  fsom$info$parameters <- match.call()
  fsom$info$date <- as.character(Sys.time())
  fsom$info$version <- as.character(utils::packageVersion("FlowSOM"))
  if(!silent) message(t[3], "\n")
  return(fsom)
}

#' Print FlowSOM object
#' 
#' @param x FlowSOM object to print information about
#' @param ...  Further arguments, not used
#' @examples 
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' flowSOM.res <- FlowSOM(fileName, compensate = TRUE, transform = TRUE,
#'                       scale = TRUE, colsToUse = c(9, 12, 14:18), nClus = 10)
#' print(flowSOM.res)
#' 
#' @export
print.FlowSOM <- function(x, ...){
  if(!is.null(x$map)){
    cat("FlowSOM model trained on", nrow(x$data), "cells and", 
        length(x$map$colsUsed), "markers, \n using a",
       paste0(x$map$xdim,"x",x$map$ydim), paste0("grid (",NClusters(x)), "clusters) and",
        NMetaclusters(x), "metaclusters.")
  
    cat("\n\nMarkers used: ", paste(x$prettyColnames[x$map$colsUsed], collapse =", "))
  } else {
    cat("FlowSOM model to train on", nrow(x$data), "cells.")
  }
  
  if(!is.null(x$metaclustering)){
    cat("\n\nMetacluster cell count:\n")
    counts <- GetCounts(x)
    print(counts)
  }
  
  if(!is.null(x$outliers)){
    n_outliers <- sum(x$outliers$Number_of_outliers)
    n_mad <- round((x$outliers[1,"Threshold"] - 
                      x$outliers[1,"Median_distance"]) / 
                     x$outliers[1,"Median_absolute_deviation"])
    cat("\n", n_outliers, paste0("cells (",
                                   round(100*n_outliers/nrow(x$data),2),
                                   "%)"),
        "are further than", n_mad,"MAD from their cluster center.")
  }
}


#' Aggregate multiple fcs files together
#' 
#' Aggregate multiple fcs files to analyze them simultaneously. 
#' A new fcs file is written, which contains about \code{cTotal} cells,
#' with \code{ceiling(cTotal/nFiles)} cells from each file. Two new columns
#' are added: a column indicating the original file by index, and a noisy 
#' version of this for better plotting opportunities (index plus or minus a 
#' value between 0 and 0.1).
#' 
#' @param fileNames   Character vector containing full paths to the fcs files
#'                    to aggregate
#' @param cTotal      Total number of cells to write to the output file
#' @param channels    Channels/markers to keep in the aggregate. Default NULL 
#'                    takes all channels of the first file.
#' @param writeOutput Whether to write the resulting flowframe to a file. 
#'                    Default FALSE
#' @param outputFile  Full path to output file. Default "aggregate.fcs"
#' @param keepOrder If TRUE, the random subsample will be ordered in the same
#'                  way as they were originally ordered in the file. Default =
#'                  FALSE.
#' @param silent If FALSE, prints an update every time it starts processing a
#'               new file. Default = FALSE. 
#' @param ...     Additional arguments to pass to read.FCS
#'                  
#' @return This function does not return anything, but will write a file with
#'         about \code{cTotal} cells to \code{outputFile}
#'
#' @seealso \code{\link{ceiling}}
#'
#' @examples
#' # Define filename
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' # This example will sample 2 times 500 cells.
#' ff_new <- AggregateFlowFrames(c(fileName, fileName), 1000)
#' 
#' @importFrom flowCore read.FCS fr_append_cols keyword colnames markernames 
#'             exprs write.FCS
#' @importFrom stats rnorm
#' 
#' @export
AggregateFlowFrames <- function(fileNames, 
                                cTotal,
                                channels = NULL,
                                writeOutput = FALSE, 
                                outputFile  = "aggregate.fcs",
                                keepOrder = FALSE, 
                                silent = FALSE,
                                ...){
  
  # Compute number of cells per file
  nFiles <- length(fileNames)
  cFile <- ceiling(cTotal/nFiles)
  
  flowFrame <- NULL
  diffNumberChannels <- FALSE
  diffMarkers <- FALSE
  
  for(i in seq_len(nFiles)){
    if(!silent) {message("Reading ", fileNames[i])}
    f <- flowCore::read.FCS(fileNames[i], ...)
    # Random sampling
    ids <- sample(seq_len(nrow(f)), min(nrow(f), cFile))
    
    if(keepOrder) ids <- sort(ids)
    
    colnames <- c("File", "File_scattered", "Original_ID")
    prev_agg <- length(grep("File[0-9]*$", colnames(f)))
    if(prev_agg > 0){
      colnames[c(1, 2)] <- paste0(colnames[c(1, 2)], prev_agg + 1)
    }
    prev_ids <- length(grep("Original_ID[0-9]*$", colnames(f)))
    if(prev_ids > 0){
      colnames[3] <- paste0(colnames[3], prev_ids + 1)
    }
    
    file_ids <- rep(i, min(nrow(f), cFile))
    
    m <- cbind(file_ids,
               file_ids + stats::rnorm(length(file_ids), 0, 0.1),
               ids)
    colnames(m) <- colnames
    
    f <- flowCore::fr_append_cols(f[ids, ], m)
    
    if(is.null(flowFrame)){
      if(is.null(channels)){
        channels <- colnames(f)
        flowFrame <- f
      } else {
        channels <- GetChannels(f, channels)
        flowFrame <- f[, c(channels, colnames(m)), drop = FALSE]
      }
      flowCore::keyword(flowFrame)[["$FIL"]] <- basename(outputFile)
      flowCore::keyword(flowFrame)[["FILENAME"]] <- basename(outputFile)
    } else {
      cols_f <- flowCore::colnames(f)
      cols_flowFrame <- flowCore::colnames(flowFrame)
      commonCols <- intersect(cols_f, cols_flowFrame)
      
      if (length(commonCols) == 0) stop("No common channels between files")
      if (!diffNumberChannels && 
          length(cols_flowFrame) != length(commonCols)){
        diffNumberChannels <- TRUE
      }
      
      if (!diffMarkers && 
          any(!flowCore::markernames(f)[commonCols] %in% 
              flowCore::markernames(flowFrame)[commonCols])){
        diffMarkers <- TRUE
      }
  
      flowCore::exprs(flowFrame) <- 
        rbind(flowCore::exprs(flowFrame)[, commonCols, drop = FALSE], 
              flowCore::exprs(f)[, commonCols, drop = FALSE])
      
    }
  }
  
  if (diffNumberChannels){
    warning("Files do not contain the same number of channels/markers")
  }
  
  if (diffMarkers){ 
    warning("Files do not contain the same markers")
  }
  
  if(writeOutput){
    flowCore::write.FCS(flowFrame, filename = outputFile)
  }
  
  return(flowFrame)
}

#' PlotFileScatters
#' 
#' Make a scatter plot per channel for all provided files
#'
#' @param  input      Either a flowSet, a flowFrame (output from the 
#'                    \code{\link{AggregateFlowFrames}} function) or a 
#'                    vector of paths pointing to fcs files
#' @param  channels   Vector of channels or markers that need to be plotted, 
#'                    if NULL (default), all channels from the input will be 
#'                    plotted
#' @param  yMargin    Optional parameter to specify the margins of the 
#'                    y-axis
#' @param  yLabel     Determines the label of the y-axis. Can be "marker" and\\or
#'                    "channel". Default = "marker".
#' @param  quantiles  If provided (default NULL), a numeric vector with values
#'                    between 0 and 1. These quantiles are indicated on the plot
#' @param  names      Optional parameter to provide filenames. If \code{NULL} 
#'                    (default), the filenames will be numbers. Duplicated 
#'                    filenames will be made unique.
#' @param  groups     Optional parameter to specify groups of files, should have
#'                    the same length as the \code{input}. Id \code{NULL} 
#'                    (default), all files will be plotted in the same color
#' @param  color      Optional parameter to provide colors. Should have the same
#'                    lengths as the number of groups (or 1 if \code{groups} is 
#'                    \code{NULL})
#' @param  legend     Logical parameter to specify whether the group levels 
#'                    should be displayed. Default is \code{FALSE}
#' @param  maxPoints Total number of data points that will be plotted per 
#'                    channel, default is 50000
#' @param  ncol       Number of columns in the final plot, optional
#' @param  nrow       Number of rows in the final plot, optional
#' @param silent      If FALSE, prints an update every time it starts processing 
#'                    a new file. Default = FALSE. 
#' @param  plotFile   Path to png file, default is "FileScatters.png". If 
#'                    \code{NULL}, the output will be a list of ggplots 
#' 
#' @return List of ggplot objects if \code{plot} is \code{FALSE}, 
#'         otherwise \code{filePlot} with plot is created.
#'         
#' @examples 
#' # Preprocessing
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' ff <- flowCore::read.FCS(fileName)
#' ff <- flowCore::compensate(ff, flowCore::keyword(ff)[["SPILL"]])
#' ff <- flowCore::transform(ff,
#'          flowCore::transformList(colnames(flowCore::keyword(ff)[["SPILL"]]),
#'                                 flowCore::logicleTransform()))
#' 
#' flowCore::write.FCS(ff[1:1000, ], file = "ff_tmp1.fcs")
#' flowCore::write.FCS(ff[1001:2000, ], file = "ff_tmp2.fcs")
#' flowCore::write.FCS(ff[2001:3000, ], file = "ff_tmp3.fcs")
#'  
#' # Make plot
#' PlotFileScatters(input = c("ff_tmp1.fcs", "ff_tmp2.fcs", "ff_tmp3.fcs"),
#'                  channels = c("Pacific Blue-A", 
#'                               "Alexa Fluor 700-A", 
#'                               "PE-Cy7-A"), 
#'                  maxPoints = 1000)
#' 
#' @import ggplot2
#' @importFrom methods is
#' @importFrom flowCore fsApply exprs
#' @importFrom dplyr tibble group_by summarise
#' @importFrom ggpubr annotate_figure ggarrange text_grob
#' @importFrom stats quantile
#'  
#' @export 
PlotFileScatters <- function(input, 
                             channels = NULL, 
                             yMargin = NULL, 
                             yLabel = c("marker"),
                             quantiles = NULL,
                             names = NULL,
                             groups = NULL, 
                             color = NULL, 
                             legend = FALSE,
                             maxPoints = 50000, 
                             ncol = NULL, 
                             nrow = NULL,
                             silent = FALSE,
                             plotFile = "FileScatters.png"){
  
  #----Warnings----
  if (!is.null(color) & !is.null(groups) & 
      length(unique(groups)) != length(color)){
    stop("Color vector length should be equal to the number of groups.")
  }
  
  if (!is.null(color) & is.null(groups) & length(color) != 1){
    stop("Color vector is too long for only 1 group.")
  }
  
  #---Read in data---
  if (is(input, "flowSet")) {
    data <- flowCore::fsApply(input, function(ff) {
      flowCore::exprs(ff)
    })
    cell_counts <- flowCore::fsApply(input, function(ff) {
      nrow(ff)
    })
    file_values <- unlist(sapply(seq_len(length(cell_counts)), 
                                 function(i) {
                                   rep(i, cell_counts[i])
                                 }))
    ff <- input[[1]]
  } else if (is(input, "flowFrame")) {
    ff <- input
    data <- flowCore::exprs(ff)
    file_values <- data[, "File"]
    input <- unique(file_values)
  } else {
    channels <- GetChannels(read.FCS(input[1]), channels)
    ff <- AggregateFlowFrames(input,
                              cTotal = maxPoints, 
                              channels = channels,
                              silent = silent)
    data <- ff@exprs
    file_values <- data[, "File"]
  }
  
  subset <- sample(seq_len(nrow(data)), min(maxPoints, nrow(data)))
  if (is.null(channels)) {
    data <- data[subset, , drop = FALSE] 
  } else {
    data <- data[subset, channels, drop = FALSE]
  }
  file_values <- file_values[subset]
  channels <- colnames(data)
  
  #----Additional warnings---
  if (!is.null(names) & length(unique(file_values)) != length(names)){
    stop("Names vector should have same length as number of files.")
  }
  if (!is.null(groups) & length(unique(file_values)) != length(groups)){
    stop("Groups vector should have same length as number of files.")
  }
  
  #----Organize file names and groups----
  if (is.null(names)) { # if no names are provided, the files will be numbered
    names <- as.character(seq_len(length(input)))
  }
  if (any(duplicated(names))){
    names <- make.unique(names)
  }
  if (is.null(groups)) { # if there are no groups, all files will be labeled "1"
    groups <- rep("1", length(unique(file_values)))
  }
  
  #----Generate plots----
  plots_list <- list()
  for (channel in channels) {
    if ("marker" %in% yLabel && length(yLabel) == 1) {
      yLabs <- GetMarkers(ff, channel)
    } else if ("channel" %in% yLabel && length(yLabel) == 1){
      yLabs <- channel
    } else if (all(c("channel", "marker") %in% yLabel && length(yLabel) == 2)){
      yLabs <- paste0(GetMarkers(ff, channel), " (", channel, ")")
    } else stop("yLabel should be \"marker\" and\\or \"channel\"")
    df <- data.frame("intensity" = data[, channel],
                     "names" = factor(names[file_values], 
                                      levels = unique(names)),
                     "group" = factor(groups[file_values], 
                                      levels = unique(groups)))
    p <- ggplot2::ggplot(df, ggplot2::aes(.data$names, .data$intensity)) +
      ggplot2::geom_jitter(position = position_jitter(width = 0.1), alpha = 0.5, 
                  ggplot2::aes(colour = .data$group), shape = ".") +
      ggplot2::ylab(yLabs) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                         vjust = 0.5)) +
      ggplot2::guides(colour = ggplot2::guide_legend(
        override.aes = list(size = 5, shape = 15, alpha = 1)))
    if (!is.null(color)) { # if manual colors are provided
      p <- p + ggplot2::scale_color_manual(values = color)
    }
    
    if (!is.null(yMargin)) { # if y margins are provided
      p <- p + ggplot2::ylim(yMargin)
    }
    
    if (!legend) { # if you don't want a legend on the plot
      p <- p + ggplot2::theme(legend.position = "none")
    }
    
    if(!is.null(quantiles)){
      my_quantile <- function(x, quantiles) {
        dplyr::tibble(intensity = stats::quantile(x, quantiles), 
                      quantile = quantiles)
      }
      
      quantile_intensities <- df %>%
        dplyr::group_by(names) %>% 
        dplyr::summarise(my_quantile(.data$intensity, quantiles))
      p <- p + ggplot2::geom_point(ggplot2::aes(x = .data$names, 
                                                y = .data$intensity), 
                          col = "black", 
                          shape = 3, #95,
                          size = 3,
                          data = quantile_intensities)
    }
    
    plots_list[[length(plots_list) + 1]] <- p
  }
  
  #----Return plots----
  if (!is.null(plotFile)) {
    if (is.null(nrow) & is.null(ncol)) {
      nrow <- floor(sqrt(length(channels)))
      ncol <- ceiling(length(channels) / nrow)
    } else if (!is.null(nrow) & !is.null(ncol)) {
      if(nrow * ncol < length(channels)) (stop("Too few rows/cols to make plot"))
    } else if (is.null(nrow)) {
      nrow <- ceiling(length(channels) / ncol)
    } else {
      ncol <- ceiling(length(channels) / nrow)
    }
    png(plotFile, 
        width = ncol * (60 + 15 * length(unique(file_values))), 
        height = 250 * nrow)
    p <- ggpubr::annotate_figure(ggarrange(plotlist = plots_list,
                                           common.legend = legend, 
                                           ncol = ncol, nrow = nrow),
                                 bottom = ggpubr::text_grob("Files"))
    print(p)  
    dev.off()
  } else {
    return(plots_list)
  }
}

#' Process a flowjo workspace file
#'
#' Reads a flowjo workspace file using the \code{\link{flowWorkspace}} library 
#' and returns a list with a matrix containing gating results and a vector with 
#' a label for each cell from a set of specified gates
#'
#' @param files       The fcs files of interest
#' @param wspFile    The FlowJo wsp file to read
#' @param group       The FlowJo group to parse. Default "All Samples".
#' @param cellTypes  Cell types to use for final labeling the cells. Should
#'                    correspond with a subset of the gate names in FlowJo.
#' @param getData    If true, flowframes are returned as well.
#' @param ...         Extra arguments to pass to CytoML::flowjo_to_gatingset
#'
#' @return This function returns a list, which for every file contains a list
#' in which the first element ("matrix") is a matrix containing filtering 
#' results for each specified gate and the second element ("manual") is a vector
#' which assigns one label to each cell. If only one file is given, only one
#' list is returned instead of a list of lists.
#'
#' @seealso \code{\link{PlotPies}}
#'
#' @examples
#'
#' # Identify the files
#' fcs_file <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' wspFile <- system.file("extdata", "gating.wsp", package = "FlowSOM")
#' 
#' # Specify the cell types of interest for assigning one label per cell
#' cellTypes <- c("B cells",
#'                 "gd T cells", "CD4 T cells", "CD8 T cells",
#'                 "NK cells", "NK T cells")
#'
#' # Parse the FlowJo workspace   
#' gatingResult <- GetFlowJoLabels(fcs_file, wspFile,
#'                                 cellTypes = cellTypes,
#'                                 getData = TRUE)
#'
#' # Check the number of cells assigned to each gate
#' colSums(gatingResult$matrix)
#' 
#' # Build a FlowSOM tree
#' flowSOM.res <- FlowSOM(gatingResult$flowFrame,
#'                        colsToUse = c(9, 12, 14:18),
#'                        nClus = 10,
#'                        seed = 1)
#'    
#'  # Plot pies indicating the percentage of cell types present in the nodes
#'  PlotPies(flowSOM.res,
#'           gatingResult$manual,
#'           backgroundValues = flowSOM.res$metaclustering)
#'
#' @export
GetFlowJoLabels <- function(files,
                            wspFile,
                            group = "All Samples",
                            cellTypes = NULL,
                            getData = FALSE,
                            ...) {
  
  ws <- CytoML::open_flowjo_xml(wspFile)
  gates <- CytoML::flowjo_to_gatingset(ws, 
                                       name = group,
                                       ...)
  
  
  files_in_wsp <- flowWorkspace::sampleNames(gates)
  counts <- as.numeric(gsub(".*_([0-9]*)$", "\\1", files_in_wsp)) 
  files_in_wsp <- gsub("_[0-9]*$", "", files_in_wsp)
  result <- list()
  for(file in files){
    print(paste0("Processing ", file))
    file_id <- grep(paste0("^\\Q", basename(file), "\\E$"), 
                    files_in_wsp)
    if(length(file_id) == 0) {stop("File ", basename(file), 
                                   " not found. Files available: \n",
                                   paste0(files_in_wsp, "\n"))}
    gate_names <- flowWorkspace::gs_get_pop_paths(gates, path = "auto")
    
    gatingMatrix <- matrix(NA,
                           nrow = counts[file_id],
                           ncol = length(gate_names),
                           dimnames = list(NULL,
                                           gate_names))
    for(gate in gate_names){
      gatingMatrix[, gate] <- 
        flowWorkspace::gh_pop_get_indices(gates[[file_id]], gate)
    }
    
    if(is.null(cellTypes)){
      cellTypes <- flowWorkspace::gs_get_leaf_nodes(gates,
                                                    path = "auto")
    }
    
    manual <- ManualVector(gatingMatrix, cellTypes)
    
    result[[file]] <- list("matrix" = gatingMatrix,
                           "manual" = manual)
    
    if (getData) {
      result[[file]]$flowFrame <- 
        flowWorkspace::gh_pop_get_data(gates[[file_id]])
    }
  }
  
  if (length(files) == 1){
    result <- result[[1]]
  }
  
  #flowWorkspace::gs_cleanup_temp(gates) 
  # Commenting this out might give issues with files in the temp folder not 
  # being cleaned up, but otherwise the getData=TRUE won't work anymore
  
  return(result)
}

#' Summarise the gating matrix into one vector, only including the cell types of
#' interest
#'
#' Extract the compensated and transformed data and all gate labels.
#'
#' @param manualMatrix Matrix containing boolean values, indicating for every
#'                      gate (column) whether the cell (row) is part of it or not.
#' @param cellTypes Cell types to use in the summary vector. All others will be
#'                   ignored and cells which do not fall in one of these gates
#'                   will get the label "Unknown". Order is important!
#'
#' @return A factor with one label for every cell
#'
#' @export
ManualVector <- function(manualMatrix, cellTypes){
  
  if(is.list(manualMatrix)){ 
    manualMatrix <- do.call(rbind, manualMatrix) 
  }
  
  manual <- rep("Unlabeled", nrow(manualMatrix))
  for(cellType in cellTypes){
    manual[manualMatrix[, cellType]] <- cellType
  }
  manual <- factor(manual, levels=c("Unlabeled", cellTypes))
  return(manual)
}

#' GetChannels
#' 
#' Get channel names for an array of markers, given a flowframe or a FlowSOM
#' object. As available in "name". \code{\link{grep}} is used to look for the 
#' markers. Other regex can be added.
#' 
#' @param object  The flowFrame or the FlowSOM object of interest 
#' @param markers Vector with markers or channels of interest. Also accepts the
#'                index of the marker found in the object.
#' @param exact   If TRUE (default), the grep pattern will be extended to
#'                start with ^\\\\Q and end with \\\\E$, so only exact matches 
#'                are possible.
#'                  
#' @return Corresponding channel names
#'
#' @seealso \code{\link{GetMarkers}}
#'
#' @examples
#' 
#'    # Read the flowFrame
#'    fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#'    ff <- flowCore::read.FCS(fileName)
#'    GetChannels(ff, c("FSC-A", "CD3", "FITC-A"))
#'    GetMarkers(ff, c("FSC-A", "CD3", "FITC-A"))
#'
#' @export
GetChannels <- function(object, markers, exact = TRUE) { 
  if (is(object, "flowFrame")) {
    object_channels <- unname(flowCore::parameters(object)@data[["name"]])
    object_markers <- unname(flowCore::parameters(object)@data[["desc"]])
  } else if (is(object, "FlowSOM")) {
    object_channels <- names(object$prettyColnames)
    object_markers <- unname(gsub(" <.*", "", object$prettyColnames))
  } else {
    stop("Object should be of class flowFrame or FlowSOM")
  }
  
  if (is.logical(markers)) markers <- which(markers)
  
  channelnames <- c()
  for (marker in markers){
    if(is.numeric(marker)) {
      iChannel <- marker
    } else {
      if(exact) marker <- paste0("^\\Q", marker, "\\E$")
      iChannel <- grep(marker, object_markers)
    }
    if (length(iChannel) != 0){
      for (i in iChannel){
        channel <- object_channels[iChannel]
        names(channel) <- object_markers[iChannel]
        channelnames <- c(channelnames, channel)
      }
    } else {
      iChannel <- grep(marker, object_channels)
      if (length(iChannel) != 0){
        channel <- object_channels[iChannel]
        names(channel) <- channel
        channelnames <- c(channelnames, channel)
      } else {
        stop(paste("Marker", marker, "could not be found"))
      }
    }
  }
  return(channelnames)
}

#' GetMarkers
#' 
#' Get marker names for an array of channels, given a flowframe or a FlowSOM 
#' object. As available in "desc". If this is NA, defaults to channel name. 
#' \code{\link{grep}} is used to look for the markers. Other regex can be added.
#' 
#' @param object   The flowFrame or the FlowSOM object of interest 
#' @param channels Vector with markers or channels of interest. Also accepts the
#'                 index of the channel in the object.
#' @param exact   If TRUE (default), the grep pattern will be extended to
#'                start with ^\\\\Q and end with \\\\E$, so only exact matches 
#'                are possible.
#'                                  
#' @return Corresponding marker names
#'
#' @seealso \code{\link{GetChannels}}
#'
#' @examples
#' 
#'    # Read the flowFrame
#'    fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#'    ff <- flowCore::read.FCS(fileName)
#'    GetChannels(ff, c("FSC-A", "CD3", "FITC-A"))
#'    GetMarkers(ff, c("FSC-A", "CD3", "FITC-A"))
#'
#' @export
GetMarkers <- function(object, channels, exact = TRUE) { 
  if (is(object, "flowFrame")) {
    object_channels <- unname(flowCore::parameters(object)@data[["name"]])
    object_markers <- unname(flowCore::parameters(object)@data[["desc"]])
  } else if (is(object, "FlowSOM")) {
    object_channels <- names(object$prettyColnames)
    object_markers <- unname(gsub(" <.*", "", object$prettyColnames))
  } else {
    stop("Object should be of class flowFrame or FlowSOM")
  }
  
  if (is.logical(channels)) channels <- which(channels)
  
  markernames <- c()
  for (channel in channels){
    if (is.numeric(channel)) {
      iMarker <- channel
    } else {
      if (exact) channel <- paste0("^\\Q", channel, "\\E$")
      iMarker <- grep(channel, object_channels)
    }
    if (length(iMarker) != 0){
      for (i in iMarker){
        marker <- object_markers[i]
        if (is.na(marker)) marker <- object_channels[i]
        names(marker) <- object_channels[i]
        markernames <- c(markernames, marker)
      }
    } else {
      iMarker <- grep(channel, object_markers)
      if (length(iMarker) != 0){
        marker <- object_markers[iMarker]
        names(marker) <- marker
        markernames <- c(markernames, marker)
      } else {
        stop(paste("Channel", channel, "could not be found"))
      }
    }
  }
  return(markernames)
}
