#' MetaClustering
#'
#' Cluster data with automatic number of cluster determination for 
#' several algorithms
#'
#' @param data   Matrix containing the data to cluster
#' @param method Clustering method to use
#' @param max    Maximum number of clusters to try out
#' @param seed   Seed to pass on to given clustering method
#' @param ...    Extra parameters to pass along
#' 
#' @return Numeric array indicating cluster for each datapoint
#' @seealso   \code{\link{metaClustering_consensus}}
#'
#' @examples
#'    # Read from file, build self-organizing map and minimal spanning tree
#'    fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#'    flowSOM.res <- ReadInput(fileName, compensate = TRUE,transform = TRUE,
#'                             scale = TRUE)
#'    flowSOM.res <- BuildSOM(flowSOM.res,colsToUse = c(9, 12, 14:18))
#'    flowSOM.res <- BuildMST(flowSOM.res)
#'    
#'    # Apply metaclustering
#'    metacl <- MetaClustering(flowSOM.res$map$codes,
#'                             "metaClustering_consensus",
#'                             max = 10)
#'    
#'    # Get metaclustering per cell
#'    flowSOM.clustering <- metacl[flowSOM.res$map$mapping[, 1]]    
#'
#' @export
MetaClustering <- function(data, method, max = 20, seed = NULL, ...){
  res <- DetermineNumberOfClusters(data, max, method, seed = seed, ...)
  method <- get(method)
  method(data, k = res, seed = seed)
}

DetermineNumberOfClusters <- function(data, max, method, plot = FALSE,
                                      smooth = 0.2, seed = NULL, ...){
  # Try out a clustering algorithm for several numbers of clusters and 
  # select optimal
  #
  # Args:
  #     data:     Matrix containing the data to cluster
  #     max:        Maximum number of clusters to try
  #     method: Clustering method to use
  #     plot:     Whether to plot the results for different k
  #     smooth: Smoothing option to find elbow: 
  #             0: no smoothing, 1: maximal smoothing
  #     seed:   Seed to pass on to given method
  #
  # Returns:
  #     Optimal number of clusters
  if(method ==    "metaClustering_consensus"){
    results <- consensus(data,max, seed, ...)
    res <- rep(0,max)
    res[1] <- SSE(data,rep(1,nrow(data)))
    for(i in 2:max){
      c <- results[[i]]$consensusClass
      res[i] <- SSE(data, c)
    }
  } else {
    method <- get(method)
    res <- rep(0,max)
    for(i in 1:max){
      c <- method(data, k = i,...)
      res[i] <- SSE(data, c)
    }
  }
  
  for(i in 2:(max - 1)){
    res[i] <- (1 - smooth) * res[i] + 
      (smooth / 2) * res[i - 1] + 
      (smooth / 2) * res[i + 1]
  }
  
  if(plot) plot(1:max, res, type = "b", xlab = "Number of Clusters", 
                ylab = "Within groups sum of squares")
  findElbow(res)
}

findElbow <- function(data){
  n <- length(data)    
  data <- as.data.frame(cbind(1:n,data))
  colnames(data) <- c("X","Y")
  
  min_r <- Inf
  optimal <- 1
  for(i in 2:(n-1)){
    f1 <- stats::lm(Y~X,data[1:(i-1),])
    f2 <- stats::lm(Y~X,data[i:n,])
    r <- sum(abs(c(f1$residuals,f2$residuals)))
    if(r < min_r){
      min_r <- r
      optimal <-i
    }
  }
  optimal
}

#' MetaClustering
#' 
#' Cluster data using hierarchical consensus clustering with k clusters
#'
#' @param data Matrix containing the data to cluster
#' @param k    Number of clusters
#' @param seed Seed to pass to consensusClusterPlus
#' 
#' @return  Numeric array indicating cluster for each datapoint
#' @seealso \code{\link{MetaClustering}}
#' @examples
#'    # Read from file, build self-organizing map and minimal spanning tree
#'    fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#'    flowSOM.res <- ReadInput(fileName, compensate = TRUE,transform = TRUE,
#'                             scale = TRUE)
#'    flowSOM.res <- BuildSOM(flowSOM.res,colsToUse = c(9, 12, 14:18))
#'    flowSOM.res <- BuildMST(flowSOM.res)
#'    
#'    # Apply consensus metaclustering
#'    metacl <- metaClustering_consensus(flowSOM.res$map$codes, k = 10)    
#'
#' @export
metaClustering_consensus <- function(data, k = 7, seed = NULL){
  results <- suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(
    t(data),
    maxK = k, reps = 100, pItem = 0.9, pFeature = 1, 
    title = tempdir(), plot = "pdf", verbose = FALSE,
    clusterAlg = "hc", # "hc","km","kmdist","pam"
    distance = "euclidean" ,
    #"euclidean","pearson","spearman","binary","maximum","canberra","minkowski"
    seed = seed
  ))
  
  results[[k]]$consensusClass
}

consensus <- function(data, max, seed = NULL, ...){
  results <- suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(
    t(data),
    maxK = max, reps = 100, pItem = 0.9, pFeature = 1,
    title = tempdir(), plot = "pdf", verbose = FALSE,
    clusterAlg = "hc", # "hc","km","kmdist","pam"
    distance = "euclidean",
    #"euclidean","pearson","spearman","binary","maximum","canberra","minkowski"
    seed = seed
  ))
}

metaClustering_hclust <- function(data, k = 7){
  d <- stats::dist(data, method = "minkowski")
  fit <- stats::hclust(d, method = "ward.D2")
  stats::cutree(fit, k = k)
}

metaClustering_kmeans <- function(data, k = 7){
  stats::kmeans(data, centers = k)$cluster
}

metaClustering_som <- function(data, k = 7){
  s <- SOM(data,xdim = k,ydim = 1,rlen = 100)
  s$unit.classif
}

SSE <- function(data,clustering){
  if(!is(clustering, "numeric"))
    clustering <- as.numeric(as.factor(clustering))
  c_wss <- 0
  for(j in seq_along(clustering)){
    if(sum(clustering == j) > 1){
      c_wss <- c_wss + (nrow(data[clustering == j, , drop = FALSE]) - 1) *
        sum(apply(data[clustering == j, , drop = FALSE], 2, stats::var))
    }
  }
  c_wss
}

#' F measure
#' 
#' Compute the F measure between two clustering results
#'
#' @param realClusters Array containing real cluster labels for each sample
#' @param predictedClusters Array containing predicted cluster labels for each
#'                          sample
#' @param silent    Logical, if FALSE (default), print some information about 
#'                  precision and recall
#' 
#' @return  F measure score
#' @examples
#' # Generate some random data as an example
#' realClusters <- sample(1:5,100,replace = TRUE)
#' predictedClusters <- sample(1:6, 100, replace = TRUE)
#' 
#' # Calculate the FMeasure
#' FMeasure(realClusters,predictedClusters)
#' @export
FMeasure <- function(realClusters, predictedClusters,silent = FALSE){
  if (sum(predictedClusters) == 0)
    return(0);
  a <- table(realClusters, predictedClusters);
  p <- t(apply(a, 1, function(x) x / colSums(a)))
  r <- apply(a, 2, function(x) x / rowSums(a))
  if(!silent) message("Precision: ",
                      sum(apply(p, 1, max) * (rowSums(a) / sum(a))),
                      "\nRecall: ", 
                      sum(apply(r, 1, max) * (rowSums(a) / sum(a))), 
                      "\n")
  f <- 2 * r * p / (r + p)
  f[is.na(f)] <- 0
  sum(apply(f, 1, max) * (rowSums(a) / sum(a)))
}

#' GetMetaclusterMFIs
#' 
#' Compute the median fluorescence intensities for the metaclusters
#'
#' @param fsom     Result of calling the FlowSOM function
#' @param colsUsed Logical. Should report only the columns used to 
#'                 build the SOM. Default = FALSE.
#' @param prettyColnames    Logical. Should report pretty column names instead
#'                          of standard column names. Default = FALSE.
#'                          
#' @return  Metacluster MFIs
#' @examples
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' ff <- flowCore::read.FCS(fileName)
#' ff <- flowCore::compensate(ff, flowCore::keyword(ff)[["SPILL"]])
#' ff <- flowCore::transform(ff,
#'          flowCore::transformList(colnames(flowCore::keyword(ff)[["SPILL"]]),
#'                                 flowCore::logicleTransform()))
#' flowSOM.res <- FlowSOM(ff,
#'                        scale = TRUE,
#'                        colsToUse = c(9, 12, 14:18),
#'                        nClus = 10)
#' mfis <- GetMetaclusterMFIs(flowSOM.res)
#' @export
GetMetaclusterMFIs <- function(fsom, colsUsed = FALSE, prettyColnames = FALSE){
  fsom <- UpdateFlowSOM(fsom)
  MFIs <- fsom$map$metaclusterMFIs
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

#' GetMetaclusterCVs
#' 
#' Compute the coefficient of variation for the metaclusters
#'
#' @param fsom Result of calling the FlowSOM function
#' 
#' @return  Metacluster CVs
#' @examples
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' ff <- flowCore::read.FCS(fileName)
#' ff <- flowCore::compensate(ff, flowCore::keyword(ff)[["SPILL"]])
#' ff <- flowCore::transform(ff,
#'          flowCore::transformList(colnames(flowCore::keyword(ff)[["SPILL"]]),
#'                                 flowCore::logicleTransform()))
#' flowSOM.res <- FlowSOM(ff,
#'                        scale = TRUE,
#'                        colsToUse = c(9, 12, 14:18), 
#'                        nClus = 10)
#' cvs <- GetMetaclusterCVs(flowSOM.res)
#' @export
GetMetaclusterCVs <- function(fsom){
  CVs <- t(sapply(seq_along(levels(fsom$metaclustering)), 
                  function(i) {
                    apply(subset(fsom$data, 
                                 fsom$metaclustering[
                                   GetClusters(fsom)] == i),
                          2,
                          function(y){
                            if(length(y) > 0 && mean(y) != 0){
                              stats::sd(y)/mean(y)
                            } else {
                              NA
                            }})
                  }))
  
  return(CVs)
}

#' UpdateMetaclusters
#' 
#' Adapt the metacluster levels. Can be used to rename the metaclusters, split 
#' or merge existing metaclusters, add a metaclustering and/or reorder the levels 
#' of the metaclustering. 
#'
#' @param fsom Result of calling the FlowSOM function.
#' @param newLabels Optional. Named vector, with the names the original 
#'                  metacluster names and the values the replacement. Can be 
#'                  used to rename or merge metaclusters.
#' @param clusterAssignment Optional. Either a named vector, with the names 
#'                          the cluster numbers (characters) or a vector of 
#'                          length NClusters(fsom). Can be used to assign 
#'                          clusters to existing or new metaclusters.
#' @param levelOrder Optional. Vector showing the preferred order of the fsom 
#'                   metacluster levels.
#'               
#' @return  Updated FlowSOM object
#' @examples
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
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
#' PlotStars(flowSOM.res, backgroundValues = flowSOM.res$metaclustering)
#' GetCounts(flowSOM.res)
#'
#' # Merge MC8 and MC9
#' flowSOM.res <- UpdateMetaclusters(flowSOM.res, newLabels = c("8" = "8+9",
#'                                                              "9" = "8+9")) 
#' PlotStars(flowSOM.res, backgroundValues = flowSOM.res$metaclustering)
#' GetCounts(flowSOM.res)
#' 
#' # Split cluster 24 from metacluster 2 and order the metacluster levels
#' flowSOM.res <- UpdateMetaclusters(flowSOM.res, 
#'                                   clusterAssignment = c("24" = "debris?"),
#'                                   levelOrder = c("debris?", as.character(c(1:7)),
#'                                                   "8+9", "10"))
#' PlotStars(flowSOM.res, backgroundValues = flowSOM.res$metaclustering)
#' PlotNumbers(flowSOM.res, level = "metaclusters")
#' 
#' GetCounts(flowSOM.res)
#' 
#' @export
UpdateMetaclusters <- function(fsom, newLabels = NULL, clusterAssignment = NULL, 
                               levelOrder = NULL){
  
  if(!is.null(fsom$metaclustering)){
    # newLabels to relabel or merge some MCs
    if(!is.null(newLabels)){
      currentLevels <- levels(fsom$metaclustering)
      newLevels <- currentLevels
      names(newLevels) <- currentLevels
      
      for(original in names(newLabels)){
        newLevels[currentLevels == original] <- newLabels[original]
      }
      
      if(any(duplicated(newLevels))){
        fsom$metaclustering <- newLevels[as.character(fsom$metaclustering)]
        fsom$metaclustering <- factor(fsom$metaclustering, 
                                      levels = unique(newLevels))
        fsom$map$nMetaclusters <- length(unique(newLevels))
      } else {
        levels(fsom$metaclustering) <- newLevels
      }
      
    }
    if(!is.null(clusterAssignment)){
      if (!is.null(names(clusterAssignment))){
        # named clusterAssignment to reassign some Cs
        currentLevels <- fsom$metaclustering
        newLevels <- as.character(currentLevels)
        names(newLevels) <- as.character(1:length(newLevels))
        
        for(original in names(clusterAssignment)){
          newLevels[original] <- clusterAssignment[original]
        }
        
        fsom$metaclustering <- factor(x = newLevels)
        fsom$map$nMetaclusters <- length(unique(newLevels))
        
      } else if (length(clusterAssignment) == NClusters(fsom)){
        # clusterAssignment of length nClusters to reassign Cs
        new <- factor(clusterAssignment)
        fsom$metaclustering <- new
        fsom$map$nMetaclusters <- length(levels(new))
      } else {
        stop("The clusterAssignment vector should be named, or the the length 
        should be equal to the number of fsom clusters.")
      }
    }
  } else if (length(clusterAssignment) == NClusters(fsom)){
    new <- factor(clusterAssignment)
    fsom$metaclustering <- new
    fsom$map$nMetaclusters <- length(levels(new))
  } else {
    stop("This FlowSOM object does not include a metaclustering and the length 
         of the clusterAssignment is not equal to the number of fsom clusters.")
  }
  if (!is.null(levelOrder)){
    fsom$metaclustering <- factor(fsom$metaclustering,
                                  levels = levelOrder)
  }
  fsom <- UpdateDerivedValues(fsom)
  return(fsom)
}
