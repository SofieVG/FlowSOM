#' Build a self-organizing map
#' 
#' Build a SOM based on the data contained in the FlowSOM object
#'
#' @param fsom      FlowSOM object containing the data, as constructed by 
#'                  the \code{\link{ReadInput}} function
#' @param colsToUse Markers, channels or indices to use for building the SOM
#' @param silent    if \code{TRUE}, no progress updates will be printed
#' @param outlierMAD Number of MAD when a cell is considered an outlier.
#'                   See also \code{\link{TestOutliers}}
#' @param ...       options to pass on to the SOM function (xdim, ydim, rlen, 
#'                  mst, alpha, radius, init, distf, importance)
#'
#' @return FlowSOM object containing the SOM result, which can be used as input
#'         for the \code{\link{BuildMST}} function
#'         
#' @seealso \code{\link{ReadInput}}, \code{\link{BuildMST}}
#' 
#' @references This code is strongly based on the \code{kohonen} package.
#'             R. Wehrens and L.M.C. Buydens, Self- and Super-organising Maps 
#'             in R: the kohonen package J. Stat. Softw., 21(5), 2007
#' 
#' @examples
#' 
#' # Read from file
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' flowSOM.res <- ReadInput(fileName, compensate = TRUE, transform = TRUE,
#'                          scale = TRUE)
#' 
#' # Build the Self-Organizing Map
#' # E.g. with gridsize 5x5, presenting the dataset 20 times, 
#' # no use of MST in neighbourhood calculations in between
#' flowSOM.res <- BuildSOM(flowSOM.res, colsToUse = c(9, 12, 14:18),
#'                         xdim = 5, ydim = 5, rlen = 20)
#' 
#' # Build the minimal spanning tree and apply metaclustering
#' flowSOM.res <- BuildMST(flowSOM.res)
#' metacl <- MetaClustering(flowSOM.res$map$codes,
#'                          "metaClustering_consensus", max = 10)
#' 
#' @export
BuildSOM <- function(fsom, 
                     colsToUse = NULL, 
                     silent = FALSE, 
                     outlierMAD = 4,
                     ...){
  if(!"data" %in% names(fsom)){
    stop("Please run the ReadInput function first!")
  }
  
  if(!silent) message("Building SOM\n")
  
  if(is.null(colsToUse)){
    colsToUse <- seq_len(ncol(fsom$data))
  }
  
  colsToUse <- GetChannels(fsom, colsToUse)
  
  fsom$map <- SOM(fsom$data[, colsToUse], silent = silent, ...)
  fsom$map$colsUsed <- colsToUse
  fsom <- UpdateDerivedValues(fsom)
  
  fsom$outliers <- TestOutliers(fsom,
                                madAllowed = outlierMAD)
  return(fsom)
}

UpdateDerivedValues <- function(fsom){
  fsom$map$medianValues <-
    t(sapply(seq_len(fsom$map$nNodes), function(i) {
      apply(subset(fsom$data, fsom$map$mapping[, 1] == i), 
            2, 
            stats::median)
    }))
  fsom$map$medianValues[is.nan(fsom$map$medianValues)] <- NA 
  colnames(fsom$map$medianValues) <- colnames(fsom$data)
  
  fsom$map$cvValues <-
    t(sapply(seq_len(fsom$map$nNodes), function(i) {
      apply(subset(fsom$data, fsom$map$mapping[, 1] == i),
            2,
            function(y){
              if(length(y) > 0 && mean(y) != 0){
                stats::sd(y)/mean(y)
              } else {
                NA
              }})
    }))
  fsom$map$cvValues[is.nan(fsom$map$cvValues)] <- NA 
  colnames(fsom$map$medianValues) <- colnames(fsom$data)
  
  fsom$map$sdValues <-
    t(sapply(seq_len(fsom$map$nNodes), function(i) {
      apply(subset(fsom$data, fsom$map$mapping[, 1] == i), 2, stats::sd)
    }))
  fsom$map$sdValues[is.nan(fsom$map$sdValues)] <- 0 
  colnames(fsom$map$sdValues) <- colnames(fsom$data)
  
  pctgs <- rep(0, fsom$map$nNodes)
  names(pctgs) <- as.character(seq_len(fsom$map$nNodes))
  pctgs_tmp <- table(fsom$map$mapping[, 1]) / nrow(fsom$map$mapping)
  pctgs[names(pctgs_tmp)] <- pctgs_tmp
  fsom$map$pctgs <- pctgs
  
  
  if(!is.null(fsom$metaclustering)){
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

#' Build a self-organizing map
#'
#' @param data  Matrix containing the training data
#' @param xdim  Width of the grid
#' @param ydim  Hight of the grid
#' @param rlen  Number of times to loop over the training data for each MST
#' @param mst   Number of times to build an MST
#' @param alpha Start and end learning rate
#' @param radius Start and end radius
#' @param init  Initialize cluster centers in a non-random way
#' @param initf Use the given initialization function if init == T
#'              (default: Initialize_KWSP)
#' @param distf Distance function (1 = manhattan, 2 = euclidean, 3 = chebyshev, 
#'              4 = cosine)
#' @param silent If FALSE, print status updates
#' @param codes Cluster centers to start with
#' @param importance array with numeric values. Parameters will be scaled 
#'                   according to importance
#'
#' @return A list containing all parameter settings and results
#'         
#' @seealso \code{\link{BuildSOM}}
#' 
#' @references This code is strongly based on the \code{kohonen} package.
#'             R. Wehrens and L.M.C. Buydens, Self- and Super-organising Maps 
#'             in R: the kohonen package J. Stat. Softw., 21(5), 2007
#' @useDynLib FlowSOM, .registration = TRUE
#' @export

SOM <- function (data, xdim = 10, ydim = 10, rlen = 10, mst = 1, 
                 alpha = c(0.05, 0.01),
                 radius = stats::quantile(nhbrdist, 0.67) * c(1, 0), 
                 init = FALSE, initf = Initialize_KWSP, distf = 2, 
                 silent = FALSE,
                 codes = NULL, importance = NULL){
  if (!is.null(codes)){
    if((ncol(codes) != ncol(data)) | (nrow(codes) != xdim * ydim)){
      stop("If codes is not NULL, it should have the same number of 
             columns as the data and the number of rows should correspond with 
             xdim*ydim")
    }
  }
  
  if(!is.null(importance)){
    data <- data * rep(importance, each = nrow(data))
  }
  
  if (is.null(colnames(data))) {
    colnames(data) <- as.character(seq_len(ncol(data)))
  }
  # Initialize the grid
  grid <- expand.grid(seq_len(xdim), seq_len(ydim))
  nCodes <- nrow(grid)
  if(is.null(codes)){
    if(init){
      codes <- initf(data, xdim, ydim)
      message("Initialization ready\n")
    } else {
      codes <- data[sample(1:nrow(data), nCodes, replace = FALSE), , 
                    drop = FALSE]
    }
  }
  
  # Initialize the neighbourhood
  nhbrdist <- as.matrix(stats::dist(grid, method = "maximum"))
  
  # Initialize the radius
  if(mst == 1){
    radius <- list(radius)
    alpha <- list(alpha)
  } else {
    radius <- seq(radius[1], radius[2], length.out = mst+1)
    radius <- lapply(1:mst, function(i){c(radius[i], radius[i+1])})
    alpha <- seq(alpha[1], alpha[2], length.out = mst+1)
    alpha <- lapply(1:mst, function(i){c(alpha[i], alpha[i+1])})
  }
  
  # Compute the SOM
  for(i in seq_len(mst)){
    res <- .C("C_SOM", data = as.double(data), 
              codes = as.double(codes), 
              nhbrdist = as.double(nhbrdist), 
              alpha = as.double(alpha[[i]]), 
              radius = as.double(radius[[i]]), 
              xdists = double(nCodes), 
              n = as.integer(nrow(data)), 
              px = as.integer(ncol(data)), 
              ncodes = as.integer(nCodes), 
              rlen = as.integer(rlen), 
              distf = as.integer(distf))
    
    codes <- matrix(res$codes, nrow(codes), ncol(codes))
    colnames(codes) <- colnames(data)
    nhbrdist <- Dist.MST(codes)
  }
  
  if(!silent) message("Mapping data to SOM\n")
  mapping <- MapDataToCodes(codes, data)
  
  return(list(xdim = xdim, ydim = ydim, rlen = rlen, mst = mst, alpha = alpha,
              radius = radius, init = init, distf = distf,
              grid = grid, codes = codes, mapping = mapping, nNodes = nCodes))
}



#' Assign nearest node to each datapoint
#
#' @param codes matrix with nodes of the SOM
#' @param newdata datapoints to assign
#' @param distf Distance function (1 = manhattan, 2 = euclidean, 3 = chebyshev, 
#'              4 = cosine)
#' 
#' @return Array with nearest node id for each datapoint
#' 
MapDataToCodes <- function (codes, newdata, distf = 2) {
  
  nnCodes <- .C("C_mapDataToCodes", 
                as.double(newdata[, colnames(codes)]), 
                as.double(codes),
                as.integer(nrow(codes)),
                as.integer(nrow(newdata)),
                as.integer(ncol(codes)),
                nnCodes = integer(nrow(newdata)),
                nnDists = double(nrow(newdata)), 
                distf = as.integer(distf))
  return(cbind(nnCodes$nnCodes, nnCodes$nnDists))
}

#' Select k well spread points from X
#' @param   X matrix in which each row represents a point
#' @param   xdim x dimension of the grid
#' @param   ydim y dimension of the grid
#'
#' @return  array containing the selected selected rows
#' 
#' @examples 
#' 
#' points <- matrix(1:1000, ncol = 10)
#' selection <- Initialize_KWSP(points, 3, 3)
#' 
#' @export
Initialize_KWSP <- function(X, xdim, ydim){
  k <- xdim * ydim
  
  # Start with a random point
  selected <- numeric(k)
  selected[1] <- sample(1:nrow(X), 1)
  dists <- apply(X, 1, function(x)sum((x - X[selected[1], ]) ^ 2))
  
  for(i in seq_len(k - 1) + 1){ #2:k
    # Add point wich is furthest away from all previously selected points
    selected[i] <- which.max(dists)
    # Update distances
    dists <- pmin(dists, 
                  apply(X, 1, function(x)sum((x - X[selected[i], ]) ^ 2)))
  }
  
  return(X[selected, ])
}

#' Create a grid from first 2 PCA components
#' @param   data matrix in which each row represents a point
#' @param   xdim x dimension of the grid
#' @param   ydim y dimension of the grid
#'
#' @return  array containing the selected selected rows
#' 
#' @examples 
#' 
#' points <- matrix(1:1000, ncol = 10)
#' selection <- Initialize_PCA(points, 3, 3)
#' 
#' @export
Initialize_PCA <- function(data, xdim, ydim){
  pca <- stats::prcomp(data, rank. = 2, retx = FALSE)
  # scale out to 5-times standard deviation, 
  # which should cover the data nicely
  sdev_scale <- 5 
  ax1 <- t(matrix(pca$rotation[, 1] * sdev_scale * pca$sdev,
                  nrow = ncol(data),
                  ncol = xdim * ydim)) *
    (2 * rep(c(1:xdim) - 1, times = ydim) / (xdim - 1) - 1)
  ax2 <- t(matrix(pca$rotation[, 2] * sdev_scale * pca$sdev,
                  nrow = ncol(data),
                  ncol = xdim * ydim)) *
    (2 * rep(c(1:ydim) - 1, each = xdim) / (ydim - 1) - 1)
  
  return(t(matrix(pca$center, 
                  nrow = ncol(data), 
                  ncol = xdim * ydim)) + 
           ax1 + ax2)
}

#' Calculate mean weighted cluster purity
#'
#' @param realClusters      array with real cluster values
#' @param predictedClusters array with predicted cluster values
#' @param weighted          logical. Should the mean be weighted
#'                          depending on the number of poins in the predicted 
#'                          clusters
#'                           
#' @return Mean purity score, worst score, number of clusters with score < 0.75
#' @examples
#' # Generate some random data as an example
#' realClusters <- sample(1:5, 100, replace = TRUE)
#' predictedClusters <- sample(1:6, 100, replace = TRUE)
#' 
#' # Calculate the FMeasure
#' Purity(realClusters, predictedClusters)
#' @export
Purity <- function(realClusters, predictedClusters, weighted = TRUE){
  
  t <- table(predictedClusters, realClusters)
  maxPercentages <- apply(t/rowSums(t), 1, max)
  
  if(weighted)
    weightedPercentages <- maxPercentages * rowSums(t)/sum(t)
  else 
    weightedPercentages <- maxPercentages/nrow(t)
  
  return(c(sum(weightedPercentages), min(maxPercentages), 
           sum(maxPercentages<0.75)))
}

#' Calculate distance matrix using a minimal spanning tree neighbourhood
#'
#' @param  X matrix in which each row represents a point
#'
#' @return Distance matrix
Dist.MST <- function(X){
  adjacency <- stats::dist(X, method = "euclidean")
  fullGraph <- igraph::graph.adjacency(as.matrix(adjacency), 
                                       mode = "undirected", 
                                       weighted = TRUE)
  mst <- igraph::minimum.spanning.tree(fullGraph)
  return(igraph::shortest.paths(mst, v = igraph::V(mst), to = igraph::V(mst), 
                                weights = NA))
}

#' Write FlowSOM clustering results to the original FCS files
#'
#' @param  fsom              FlowSOM object as generated by BuildSOM
#' @param  originalFiles     FCS files that should be extended
#' @param  preprocessedFiles FCS files that correspond to the input of FlowSOM,
#'                           If NULL (default), the originalFiles are used.
#' @param  selectionColumn   Column of the FCS file indicating the original cell
#'                           ids. If NULL (default), no selection is made.
#' @param  silent            If FALSE (default), print some extra output
#' @param  outputDir         Directory to save the fcs files. Default to the
#'                           current working directory (".")
#' @param  suffix            Suffix added to the filename. Default _FlowSOM.fcs
#'
#' @return Saves the extended fcs file as [originalName]_FlowSOM.fcs
#'
#' @examples
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' flowSOM.res <- FlowSOM(fileName, compensate = TRUE, transform = TRUE,
#'                       scale = TRUE, colsToUse = c(9, 12, 14:18), nClus = 10)
#' SaveClustersToFCS(flowSOM.res, fileName)
#' 
#' @importFrom stats runif
#' @export 
SaveClustersToFCS <- function(fsom, 
                              originalFiles, 
                              preprocessedFiles = NULL, 
                              selectionColumn = NULL, 
                              silent = FALSE,
                              outputDir = ".",
                              suffix = "_FlowSOM.fcs"){
  if (length(preprocessedFiles) != 0 & length(preprocessedFiles)!= length(originalFiles)){ 
    stop("The vector of preprocessedFiles should be the same length as ",
         "the originalFiles when provided.")
  }
  
  for(i in seq_along(originalFiles)){
    if(!silent){message("Mapping ", originalFiles[i])}
    
    ff_o <- flowCore::read.FCS(originalFiles[i])
    if(!is.null(preprocessedFiles)){
      ff <- flowCore::read.FCS(preprocessedFiles[i])
    } else {
      ff <- ff_o
    }
    
    if(!is.null(selectionColumn)){
      s <- flowCore::exprs(ff)[, selectionColumn]
    } else {
      s <- seq_len(nrow(ff_o))    
    }
    
    # Map the data
    fsom_f <- NewData(fsom, ff)
    
    # Put on corresponding indices
    newCols <- ifelse(is.null(fsom$metaclustering), 3, 4)
    m <- matrix(0, nrow = nrow(ff_o), ncol = newCols)
    
    m[s, 1] <- GetClusters(fsom_f)
    
    spread <- min(dist(fsom_f$MST$l))/3
    m[s, 2] <- fsom_f$MST$l[,1][GetClusters(fsom_f)] + 
      stats::runif(length(s), min = -spread, max = spread)
    m[s, 3] <- fsom_f$MST$l[,2][GetClusters(fsom_f)] + 
      stats::runif(length(s), min = -spread, max = spread)
    if(!is.null(fsom_f$metaclustering)){
      m[s, 4] <- GetMetaclusters(fsom_f)
    }
    colnames(m) <- c("FlowSOM_cluster", "FlowSOM_x", "FlowSOM_y", "FlowSOM_meta")[seq_len(newCols)]
    
    # Save as fcs file
    ff_o <- flowCore::fr_append_cols(ff_o, m)
    outputFile <- file.path(outputDir,
                            gsub("\\.fcs$",
                                 suffix,
                                 basename(originalFiles[i])))
    flowCore::write.FCS(ff_o, 
                        filename = outputFile)  
    
    if(!silent){message("Result written to ", outputFile)}
  }
}

#' Get cluster label for all individual cells
#'
#' @param  fsom             FlowSOM object as generated by the FlowSOM function
#'                          or the BuildSOM function
#'                          
#' @return vector label for every cell
#' @examples 
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' flowSOM.res <- FlowSOM(fileName, compensate = TRUE, transform = TRUE,
#'                       scale = TRUE, colsToUse = c(9, 12, 14:18), nClus = 10)
#' cluster_labels <- GetClusters(flowSOM.res)
#'
#' @export 
GetClusters <- function(fsom) {
  fsom <- UpdateFlowSOM(fsom)
  return(fsom$map$mapping[, 1])
}

#' Get metacluster label for all individual cells
#'
#' @param  fsom             FlowSOM object as generated by the FlowSOM function
#'                          or the BuildSOM function
#' @param  meta             Metacluster label for each FlowSOM cluster. If this
#'                          is NULL, the fsom argument should be as generated by
#'                          the FlowSOM function, and fsom$metaclustering will
#'                          be used.                          
#' @return vector label for every cell
#' @examples 
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' flowSOM.res <- FlowSOM(fileName, compensate = TRUE, transform = TRUE,
#'                       scale = TRUE, colsToUse = c(9, 12, 14:18), nClus = 10)
#' metacluster_labels <- GetMetaclusters(flowSOM.res)
#' metacluster_labels <- GetMetaclusters(flowSOM.res,
#'                                       meta = flowSOM.res$metaclustering)
#'
#' @export 
GetMetaclusters <- function(fsom, meta = NULL){
  fsom <- UpdateFlowSOM(fsom)
  if(is.null(meta)) meta <- fsom$metaclustering
  return(meta[GetClusters(fsom)])
} 

#' Get MFI values for all clusters
#'
#' @param  fsom             FlowSOM object as generated by the FlowSOM function
#'                          or the BuildSOM function
#' @param  colsUsed         logical. Should report only the columns used to 
#'                          build the SOM. Default = FALSE.
#' @param  prettyColnames   logical. Should report pretty column names instead
#'                          of standard column names. Default = FALSE.
#'                          
#' @return Matrix with median values for each marker
#'
#' @examples 
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' flowSOM.res <- FlowSOM(fileName, compensate = TRUE, transform = TRUE,
#'                       scale = TRUE, colsToUse = c(9, 12, 14:18), nClus = 10)
#' mfis <- GetClusterMFIs(flowSOM.res)
#' @export 
GetClusterMFIs <- function(fsom, colsUsed = FALSE, prettyColnames = FALSE){
  fsom <- UpdateFlowSOM(fsom)
  MFIs <- fsom$map$medianValues
  rownames(MFIs) <- seq_len(nrow(MFIs))
  if(is.null(fsom$map$colsUsed)) colsUsed <- FALSE
  if(is.null(fsom$prettyColnames)) prettyColnames <- FALSE
  if(colsUsed && !prettyColnames){
    MFIs <- MFIs[, fsom$map$colsUsed]
  } else if(!colsUsed && prettyColnames) {
    colnames(MFIs) <- fsom$prettyColnames
  } else if(colsUsed && prettyColnames) {
    MFIs <- MFIs[, fsom$map$colsUsed]
    colnames(MFIs) <- fsom$prettyColnames[fsom$map$colsUsed]
  }
  return(MFIs)
}

#' Get CV values for all clusters
#'
#' @param  fsom             FlowSOM object as generated by the FlowSOM function
#'                          or the BuildSOM function
#'                          
#' @return Matrix with coefficient of variation values for each marker
#' 
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' flowSOM.res <- FlowSOM(fileName, compensate = TRUE, transform = TRUE,
#'                       scale = TRUE, colsToUse = c(9, 12, 14:18), nClus = 10)
#' cvs <- GetClusterCVs(flowSOM.res)
#'
#' @export
GetClusterCVs <- function(fsom){
  fsom <- UpdateFlowSOM(fsom)
  return(fsom$map$cvValues)
}

#' GetFeatures
#' 
#' Map fcs files on an existing FlowSOM object
#'
#' @param  fsom             FlowSOM object as generated by the FlowSOM function
#'                          or the BuildSOM function
#' @param  files            Either a vector of fcs files or paths to fcs files  
#' @param  level            Level(s) of interest. Default is c("clusters",
#'                          "metaclusters"), but can also be only one of them
#' @param  type             Type of features to extract. Default is "counts", 
#'                          can be a vector of "counts", "percentages" and/or 
#'                          "MFIs"         
#' @param  MFI              Vector with channels / markers for which the MFI 
#'                          values must be returned when "MFIs" is in \code{type}
#' @param  filenames        An optional vector with filenames that will be used
#'                          as rownames in the count matrices. If NULL (default)
#'                          either the paths will be used or a numerical vector.
#' @param  silent           Logical. If \code{TRUE}, print progress messages.
#'                          Default = \code{FALSE}.
#' 
#' @return matrix with features per population - type combination
#'         
#' @examples 
#'  # Build FlowSom result
#'  fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#'  ff <- flowCore::read.FCS(fileName)
#'  ff <- flowCore::compensate(ff, flowCore::keyword(ff)[["SPILL"]])
#'  ff <- flowCore::transform(ff,
#'          flowCore::transformList(colnames(flowCore::keyword(ff)[["SPILL"]]),
#'                                 flowCore::logicleTransform()))
#'  flowSOM.res <- FlowSOM(ff[1:1000, ], 
#'                         scale = TRUE, 
#'                         colsToUse = c(9, 12, 14:18),
#'                         nClus = 10)
#'    
#'  # Map new data
#'  counts <- GetFeatures(fsom = flowSOM.res, 
#'                        level = "clusters",
#'                        files = c(ff[1001:2000, ], ff[2001:3000, ]))
#'  features <- GetFeatures(fsom = flowSOM.res, 
#'                          files = c(ff[1001:2000, ], ff[2001:3000, ]),
#'                          type = c("counts", "percentages", "MFIs"), 
#'                          MFI = "APC-A", 
#'                          filenames = c("ff_1001-2000", "ff_2001-3000"))
#'
#' @export 
GetFeatures <- function(fsom, 
                        files, 
                        level = c("clusters", "metaclusters"),
                        type = "counts", 
                        MFI = NULL, 
                        filenames = NULL, 
                        silent = FALSE) {
  nclus <- NClusters(fsom)
  nfiles <- length(files)
  i <- 0
  
  #----Warnings----
  if (!is.null(filenames) & length(filenames) != nfiles){
    stop("Filenames vector should have same length as files vector.")
  }
  
  if (sum(level %in% c("clusters", "metaclusters")) != length(level)){
    stop("level should be \"clusters\" and/or \"metaclusters\".")
  }
  
  if (sum(type %in% c("counts", "percentages", "MFIs")) != length(type)){
    stop("level should be \"counts\", \"percentages\" and/or \"MFIs\".")
  }
  
  if ("MFIs" %in% type & is.null(MFI)){
    stop("Please provide channel names for MFI calculation")
  }
  
  matrices <- list()
  
  #----Prepare matrices----
  if (is.null(filenames)) {
    if (is.character(files)) {
      filenames <- files
    } else {
      filenames <- as.character(seq_len(length(files)))
    }
  }
  C_counts <- matrix(data = 0,
                     nrow = nfiles,
                     ncol = nclus,
                     dimnames = list(filenames, 
                                     paste0("C", seq_len(nclus))))
  C_outliers <- matrix(data = 0,
                       nrow = nfiles,
                       ncol = nclus,
                       dimnames = list(filenames, 
                                       paste0("C", seq_len(nclus))))
  
  
  if ("MFIs" %in% type) {
    nmetaclus <- NMetaclusters(fsom)
    MFI <- GetChannels(fsom, MFI)
    nmarker <- length(MFI)
    C_MFIs <- matrix(NA,
                     nrow = nfiles,
                     ncol = nmarker * nclus,
                     dimnames = list(filenames,
                                     paste0(rep(paste0("C", seq_len(nclus)), 
                                                each = nmarker), 
                                            " ", fsom$prettyColnames[MFI])))
    MC_MFIs <- matrix(NA,
                      nrow = nfiles,
                      ncol = nmarker * nmetaclus,
                      dimnames = list(filenames,
                                      paste0(rep(paste0("MC", 
                                                        seq_len(nmetaclus)), 
                                                 each = nmarker), 
                                             " ", fsom$prettyColnames[MFI])))
  }
  
  #----Loop over files----
  for (file in files){
    i <- i + 1
    if (isFALSE(silent)){
      message(paste0("Mapping file ", i, " of ", nfiles, "."))
    }
    fsom_tmp <- suppressWarnings(NewData(fsom = fsom,
                                         input = file,
                                         silent = silent))
    
    counts_t <- table(GetClusters(fsom_tmp))
    C_counts[i, paste0("C", names(counts_t))] <- counts_t
    outliers_t <- fsom_tmp$outliers[, "Number_of_outliers", drop = FALSE]
    if (nrow(outliers_t) != 0) {
      C_outliers[i, paste0("C", rownames(outliers_t))] <- 
        outliers_t$Number_of_outliers
    }
    
    if ("MFIs" %in% type){
      if ("clusters" %in% level){
        C_MFIs[i, ] <- as.vector(t(GetClusterMFIs(fsom_tmp)[, MFI]))
      }
      if ("metaclusters" %in% level){
        MC_MFIs[i, ] <- as.vector(t(GetMetaclusterMFIs(fsom_tmp)[, MFI]))
      }
    }
  }
  
  #----Add matrices to list----
  if ("clusters" %in% level){
    if ("counts" %in% type){
      C_counts_tmp <- C_counts
      attr(C_counts_tmp, "outliers") <- C_outliers
      matrices[["cluster_counts"]] <- C_counts_tmp
    }
    if ("percentages" %in% type){
      C_pctgs <- prop.table(C_counts, margin = 1)
      colnames(C_pctgs) <- paste0("%", colnames(C_pctgs))
      matrices[["cluster_percentages"]] <- C_pctgs
    }
    if ("MFIs" %in% type){
      matrices[["cluster_MFIs"]] <- C_MFIs
    }
  }
  
  if ("metaclusters" %in% level){
    MC_counts <- t(apply(C_counts,
                         1,
                         function(x){
                           tapply(x, fsom$metaclustering, sum)
                         }))
    MC_counts[is.na(MC_counts)] <- 0
    colnames(MC_counts) <- paste0("MC", colnames(MC_counts))
    
    if ("counts" %in% type){
      matrices[["metacluster_counts"]] <- MC_counts
    }
    if ("percentages" %in% type){
      MC_pctgs <- prop.table(MC_counts, margin = 1)
      colnames(MC_pctgs) <- paste0("%", colnames(MC_pctgs))
      matrices[["metacluster_percentages"]] <- MC_pctgs
    }
    if ("MFIs" %in% type){
      matrices[["metacluster_MFIs"]] <- MC_MFIs
    }
  }
  
  return(matrices)
} 

#' GroupStats 
#' 
#' Calculate statistics between 2 groups based on the \code{\link{GetFeatures}}
#' output
#' 
#' @param features Feature matrix as generated by \code{\link{GetFeatures}}, 
#'                 e.g. a percentages matrix
#' @param groups   Named list with file or patient IDs per group (should match
#'                 with the rownames of the \code{matrix}).
#'                          
#' @return Matrix with the medians per group, the p-values (the raw, Benjamini 
#' Hochberg corrected one and the -log10) that resulted from a Wilcox test and 
#' the fold and log10 fold changes between the medians of the 2 groups
#' 
#' @examples 
#' # Build FlowSom result
#'  fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#'  ff <- flowCore::read.FCS(fileName)
#'  ff <- flowCore::compensate(ff, flowCore::keyword(ff)[["SPILL"]])
#'  ff <- flowCore::transform(ff,
#'          flowCore::transformList(colnames(flowCore::keyword(ff)[["SPILL"]]),
#'                                 flowCore::logicleTransform()))
#'  flowSOM.res <- FlowSOM(ff, scale = TRUE, colsToUse = c(9, 12, 14:18),
#'                           nClus = 10)
#'    
#' # Create new data
#' # To illustrate the output, we here generate new fcs files (with more 
#' # cells in metaclusters 1 and 9).
#' # In practice you would not generate any new file but use your different
#' # files from your different groups
#'  flowCore::write.FCS(ff[sample(1:nrow(ff), 1000), ], file = "ff_tmp1.fcs")
#'  flowCore::write.FCS(ff[sample(1:nrow(ff), 1000), ], file = "ff_tmp2.fcs")
#'  flowCore::write.FCS(ff[sample(1:nrow(ff), 1000), ], file = "ff_tmp3.fcs")
#'  ff_tmp <- ff[c(1:1000,
#'                 which(flowSOM.res$map$mapping[, 1] %in% 
#'                      which(flowSOM.res$metaclustering == 9)),
#'                 which(flowSOM.res$map$mapping[, 1] %in% 
#'                      which(flowSOM.res$metaclustering == 1))), ]
#'  flowCore::write.FCS(ff_tmp[sample(1:nrow(ff_tmp), 1000), ],
#'                      file = "ff_tmp4.fcs")
#'  flowCore::write.FCS(ff_tmp[sample(1:nrow(ff_tmp), 1000), ], 
#'                      file = "ff_tmp5.fcs")
#'  
#' # Get the count matrix
#'  percentages <- GetFeatures(fsom = flowSOM.res, 
#'                             files = c("ff_tmp1.fcs", 
#'                                       "ff_tmp2.fcs", 
#'                                       "ff_tmp3.fcs",
#'                                       "ff_tmp4.fcs", 
#'                                       "ff_tmp5.fcs"), 
#'                             type = "percentages")
#'                        
#'   
#' # Perform the statistics
#' groups <- list("Group 1" = c("ff_tmp1.fcs", "ff_tmp2.fcs", "ff_tmp3.fcs"), 
#'                "Group 2" = c("ff_tmp4.fcs", "ff_tmp5.fcs"))
#' MC_stats <- GroupStats(percentages[["metacluster_percentages"]], groups)
#' C_stats <- GroupStats(percentages[["cluster_percentages"]], groups)
#' 
#' # Process the fold changes vector
#' fold_changes <- C_stats["fold changes", ]
#' fold_changes <- factor(ifelse(fold_changes < -3, 
#'                               "Underrepresented compared to Group 1",
#'                               ifelse(fold_changes > 3, 
#'                                      "Overrepresented compared to Group 1",
#'                                       "--")), 
#'                         levels = c("--", 
#'                                    "Underrepresented compared to Group 1",
#'                                    "Overrepresented compared to Group 1"))
#' fold_changes[is.na(fold_changes)] <- "--"
#' 
#' # Show in figure
#' ## Fold change
#' gr_1 <- PlotStars(flowSOM.res, 
#'                   title = "Group 1", 
#'                   nodeSizes = C_stats["medians Group 1", ], 
#'                   list_insteadof_ggarrange = TRUE)
#' gr_2 <- PlotStars(flowSOM.res, title = "Group 2", 
#'             nodeSizes = C_stats["medians Group 2", ], 
#'             backgroundValues = fold_changes,
#'             backgroundColors = c("white", "red", "blue"), 
#'             list_insteadof_ggarrange = TRUE)
#' p <- ggpubr::ggarrange(plotlist = c(list(gr_1$tree), gr_2),
#'                        heights = c(3, 1))
#' ggplot2::ggsave("Groups_foldchanges.pdf", p, width = 10)
#' 
#' ## p values
#' p <- PlotVariable(flowSOM.res, title = "Wilcox test group 1 vs. group 2",
#' variable = C_stats["p values", ])
#' ggplot2::ggsave("Groups_pvalues.pdf", p)
#' 
#' ## volcano plot
#' p <- ggplot2::ggplot(data.frame("-log10 p values" = c(C_stats[4, ], 
#'                                                       MC_stats[4, ]), 
#'                                 "log10 fold changes" = c(C_stats[7, ],
#'                                                          MC_stats[7, ]), 
#' check.names = FALSE), ggplot2::aes(x = `log10 fold changes`, 
#'                                    y = `-log10 p values`)) +
#' ggplot2::xlim(-3, 3) +
#' ggplot2::ylim(0, 3) +
#' ggplot2::geom_point() 
#' 
#' @importFrom stats p.adjust wilcox.test
#' 
#' @export
GroupStats <- function(features, groups){
  nGroups <- lapply(groups, 
                    length)
  
  #----Warnings----
  if (length(groups) != 2 | !is.list(groups)){
    stop("Groups should be a named list with 2 groups")
  }
  
  if (!all(c(groups[[1]], groups[[2]]) %in% rownames(features))){
    stop("File or patient IDs from groups should correspond to ",
         "the features' rownames")
  }
  
  
  #----Calculate medians per group----
  i <- rep(NA, nrow(features))
  for (group in names(groups)){
    i[rownames(features) %in% groups[[group]]] <- group
  }
  
  medians <- apply(features, 2, function(x) {
    tapply(x, 
           INDEX = factor(i, 
                          levels = names(groups)), 
           stats::median, 
           na.rm = TRUE)
  })
  
  #----Calculate fold changes between groups----
  fold_change <- c()
  for (col in seq_len(ncol(medians))){
    m <- medians[, col]
    if (any(m <= 0) | any(is.na(m))){
      fold_change <- c(fold_change, NA)
    } else{
      fold_change <- c(fold_change, 
                       max(m, na.rm = TRUE) / 
                         min(m, na.rm = TRUE) * 
                         (-1) ^ which.max(m))
    }
  }
  
  logfold <- sign(fold_change)*log10(abs(fold_change))
  
  #----Perform Wilcoxon test between groups----
  pValues <- c()
  for (col in seq_len(ncol(features))){
    if(! (all(is.na(features[c(groups[[1]], groups[[2]]), col])) | 
          all(features[c(groups[[1]], groups[[2]]), col] == 0))){
      test <- stats::wilcox.test(features[groups[[1]], col], 
                                 features[groups[[2]], col],
                                 exact = FALSE)
      pValues <- c(pValues, test$p.value)
    } else {
      pValues <- c(pValues, 1)
    }
    
  }
  adjustedP <- stats::p.adjust(pValues, "BH")
  logP <- -log10(pValues)
  
  stats <- rbind(medians, pValues, logP, adjustedP, fold_change, logfold)
  rownames(stats) <- c(paste0("medians ", names(groups)), 
                       "p values",
                       "-log10 p values", 
                       "adjusted p values",
                       "fold changes", 
                       "log10 fold changes")
  return(stats)
}

#' GetCounts 
#' 
#' Get counts of number of cells in clusters or metaclusters
#' 
#' @param fsom        FlowSOM object
#' @param level       Character string, should be either "clusters" or 
#'                    "metaclusters" (default)
#' 
#' @return A named vector with the counts
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
#' GetCounts(flowSOM.res)                      
#' GetCounts(flowSOM.res, level = "clusters")
#' @export
GetCounts <- function(fsom, level = "metaclusters"){
  fsom <- UpdateFlowSOM(fsom)
  if (!is.null(fsom$metaclustering) && level == "metaclusters"){
    counts <- rep(NA, NMetaclusters(fsom))
    names(counts) <- paste("MC", levels(fsom$metaclustering))
    tmp <- table(GetMetaclusters(fsom))
    counts[paste("MC", names(tmp))] <- tmp
  } else if (level == "clusters"){
    counts <- rep(NA, NClusters(fsom))
    names(counts) <- paste("C", seq_len(NClusters(fsom)))
    tmp <- table(GetClusters(fsom))
    counts[paste("C", names(tmp))] <- tmp
  } else stop("level should be \"clusters\" or \"metaclusters\"")
  return(counts)
}

#' GetPercentages
#' 
#' Get percentages of number of cells in clusters or metaclusters
#' 
#' @param fsom        FlowSOM object
#' @param level       Character string, should be either "clusters" or 
#'                    "metaclusters" (default)
#' 
#' @return A named vector with the percentages
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
#' GetPercentages(flowSOM.res)                      
#' GetPercentages(flowSOM.res, level = "clusters")
#' @export
GetPercentages <- function(fsom, level = "metaclusters"){
  counts <- GetCounts(fsom, level = level)
  return(counts / sum(counts, na.rm = TRUE))
}
