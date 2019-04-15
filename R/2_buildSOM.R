#' Build a self-organizing map
#' 
#' Build a SOM based on the data contained in the FlowSOM object
#'
#' @param fsom      FlowSOM object containing the data, as constructed by 
#'                  the \code{\link{ReadInput}} function
#' @param colsToUse column names or indices to use for building the SOM
#' @param silent    if \code{TRUE}, no progress updates will be printed
#' @param ...       options to pass on to the SOM function (xdim, ydim, rlen, 
#'                  mst, alpha, radius, init, distf, importance)
#'
#' @return FlowSOM object containing the SOM result, which can be used as input
#'         for the \code{\link{BuildMST}} function
#'         
#' @seealso \code{\link{ReadInput}},\code{\link{BuildMST}}
#' 
#' @references This code is strongly based on the \code{kohonen} package.
#'             R. Wehrens and L.M.C. Buydens, Self- and Super-organising Maps 
#'             in R: the kohonen package J. Stat. Softw., 21(5), 2007
#' 
#' @examples
#' 
#' # Read from file
#' fileName <- system.file("extdata","lymphocytes.fcs",package="FlowSOM")
#' flowSOM.res <- ReadInput(fileName, compensate=TRUE,transform=TRUE,
#'                          scale=TRUE)
#' 
#' # Build the Self-Organizing Map
#' # E.g. with gridsize 5x5, presenting the dataset 20 times, 
#' # no use of MST in neighbourhood calculations in between
#' flowSOM.res <- BuildSOM(flowSOM.res,colsToUse=c(9,12,14:18),
#'                         xdim=5,ydim=5,rlen=20)
#' 
#' # Build the minimal spanning tree and apply metaclustering
#' flowSOM.res <- BuildMST(flowSOM.res)
#' metacl <- MetaClustering(flowSOM.res$map$codes,
#'                          "metaClustering_consensus",max=10)
#' 
#' @export
BuildSOM <- function(fsom, colsToUse=NULL, silent=FALSE, ...){
    if(!"data" %in% names(fsom)){
        stop("Please run the ReadInput function first!")
    }
    
    if(!silent) message("Building SOM\n")
    
    if(is.null(colsToUse)){
        colsToUse <- seq_len(ncol(fsom$data))
    }
    
    fsom$map <- SOM(fsom$data[, colsToUse],silent=silent, ...)
    fsom$map$colsUsed <- colsToUse
    fsom <- UpdateDerivedValues(fsom)
    fsom
}

UpdateDerivedValues <- function(fsom){
  fsom$map$medianValues <-
    t(sapply(seq_len(fsom$map$nNodes), function(i) {
      apply(subset(fsom$data, fsom$map$mapping[,1] == i),2,stats::median)
    }))
  fsom$map$medianValues[is.nan(fsom$map$medianValues)] <- NA 
  colnames(fsom$map$medianValues) <- colnames(fsom$data)
  
  fsom$map$cvValues <-
    t(sapply(seq_len(fsom$map$nNodes), function(i) {
      apply(subset(fsom$data, fsom$map$mapping[,1] == i),
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
      apply(subset(fsom$data, fsom$map$mapping[,1] == i),2,stats::sd)
    }))
  fsom$map$sdValues[is.nan(fsom$map$sdValues)] <- 0 
  colnames(fsom$map$sdValues) <- colnames(fsom$data)
  
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
#' @param initf Use the given initialization function if init==T
#'              (default: Initialize_KWSP)
#' @param distf Distance function (1=manhattan, 2=euclidean, 3=chebyshev, 
#'              4=cosine)
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

SOM <- function (data, xdim=10, ydim=10, rlen=10, mst=1, alpha=c(0.05, 0.01),
                    radius = stats::quantile(nhbrdist, 0.67) * c(1, 0), 
                    init=FALSE, initf=Initialize_KWSP, distf=2, silent=FALSE,
                    codes=NULL, importance = NULL){
    if (!is.null(codes)){
      if((ncol(codes) != ncol(data)) | (nrow(codes) != xdim * ydim)){
        stop("If codes is not NULL, it should have the same number of columns
             as the data and the number of rows should correspond with 
             xdim*ydim")
      }
    }
    
    if(!is.null(importance)){
        data <- data * rep(importance,each=nrow(data))
    }
    
    # Initialize the grid
    grid <- expand.grid(seq_len(xdim),seq_len(ydim))
    nCodes <- nrow(grid)
    if(is.null(codes)){
        if(init){
            codes <- initf(data, xdim, ydim)
            message("Initialization ready\n")
        } else {
            codes <- data[sample(1:nrow(data), nCodes, replace = FALSE), , drop = FALSE]
        }
    }
    
    # Initialize the neighbourhood
    nhbrdist <- as.matrix(stats::dist(grid, method = "maximum"))
    
    # Initialize the radius
    if(mst==1){
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
    mapping <- MapDataToCodes(codes,data)
    
    list(xdim=xdim, ydim=ydim, rlen=rlen, mst=mst, alpha=alpha,
        radius=radius, init=init, distf=distf,
        grid=grid, codes=codes, mapping=mapping, nNodes=nCodes)
}



#' Assign nearest node to each datapoint
#
#' @param codes matrix with nodes of the SOM
#' @param newdata datapoints to assign
#' @param distf Distance function (1=manhattan, 2=euclidean, 3=chebyshev, 
#'              4=cosine)
#' 
#' @return Array with nearest node id for each datapoint
#' 
MapDataToCodes <- function (codes, newdata, distf=2) {
    
    nnCodes <- .C("C_mapDataToCodes", 
                    as.double(newdata[,colnames(codes)]), 
                    as.double(codes),
                    as.integer(nrow(codes)),
                    as.integer(nrow(newdata)),
                    as.integer(ncol(codes)),
                    nnCodes = integer(nrow(newdata)),
                    nnDists = double(nrow(newdata)), 
                    distf = as.integer(distf))
    cbind(nnCodes$nnCodes, nnCodes$nnDists)
}

#' Select k well spread points from X
#' @param   X matrix in which each row represents a point
#' @param   xdim x dimension of the grid
#' @param   ydim y dimension of the grid
#'
#' @return  array containing the selected selected rows
#' @export
Initialize_KWSP <- function(X, xdim, ydim){
    k <- xdim * ydim

    # Start with a random point
    selected <- numeric(k)
    selected[1] <- sample(1:nrow(X), 1)
    dists <- apply(X, 1, function(x)sum((x-X[selected[1], ])^2))
    
    for(i in seq_len(k-1)+1){ #2:k
        # Add point wich is furthest away from all previously selected points
        selected[i] <- which.max(dists)
        # Update distances
        dists <- pmin(dists, 
                    apply(X, 1, function(x)sum((x-X[selected[i], ])^2)))
    }
    
    X[selected,]
}

#' Create a grid from first 2 PCA components
#' @param   data matrix in which each row represents a point
#' @param   xdim x dimension of the grid
#' @param   ydim y dimension of the grid
#'
#' @return  array containing the selected selected rows
#' @export
Initialize_PCA <- function(data, xdim, ydim){
    pca <- stats::prcomp(data, rank.=2, retx=F)
    # scale out to 5-times standard deviation, 
    # which should cover the data nicely
    sdev_scale <- 5 
    ax1 <- t(matrix(pca$rotation[,1] * sdev_scale * pca$sdev,
             nrow=ncol(data),
             ncol=xdim * ydim)) *
             (2 * rep(c(1:xdim) - 1, times=ydim) / (xdim - 1) - 1)
    ax2 <- t(matrix(pca$rotation[,2] * sdev_scale * pca$sdev,
             nrow=ncol(data),
             ncol=xdim * ydim)) *
             (2 * rep(c(1:ydim) - 1, each=xdim) / (ydim - 1) - 1)

    t(matrix(pca$center, nrow=ncol(data), ncol=xdim * ydim)) + ax1 + ax2
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
#' realClusters <- sample(1:5,100,replace = TRUE)
#' predictedClusters <- sample(1:6, 100, replace = TRUE)
#' 
#' # Calculate the FMeasure
#' Purity(realClusters,predictedClusters)
#' @export
Purity <- function(realClusters, predictedClusters, weighted=TRUE){
    
    t <- table(predictedClusters, realClusters)
    maxPercentages <- apply(t/rowSums(t), 1, max)
    
    if(weighted)
        weightedPercentages <- maxPercentages * rowSums(t)/sum(t)
    else 
        weightedPercentages <- maxPercentages/nrow(t)
    
    c(sum(weightedPercentages), min(maxPercentages), 
            sum(maxPercentages<0.75))
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
    igraph::shortest.paths(mst, v=igraph::V(mst), to=igraph::V(mst), 
                            weights=NA)
}


#' Calculate differences in cell counts between groups
#'
#' @param  fsom     FlowSOM object as generated by BuildSOM
#' @param  groups   List containing an array with file names for each group
#' @param  plot     Logical. If TRUE, make a starplot of each individual file
#' @param  silent   Logical. If TRUE, print progress messages
#'
#' @return Distance matrix
#'
#' @examples
#'    library(FlowSOM)
#'    set.seed(1)
#'    
#'    # Build the FlowSOM tree on the example file
#'    fileName <- system.file("extdata","lymphocytes.fcs",package="FlowSOM")
#'    flowSOM.res <- FlowSOM(fileName, compensate=TRUE,transform=TRUE,
#'                     scale=TRUE,colsToUse=c(9,12,14:18),nClus = 10)
#'    
#'    # Have a look at the resulting tree
#'    PlotStars(flowSOM.res[[1]],backgroundValues = as.factor(flowSOM.res[[2]]))
#'    
#'    # Select all cells except the branch that corresponds with automated 
#'    # cluster 7 (CD3+ TCRyd +) and write te another file for the example
#'    # In practice you would not generate any new file but use your different
#'    # files from your different groups
#'    ff <- flowCore::read.FCS(fileName)
#'    ff_tmp <- ff[flowSOM.res[[1]]$map$mapping[,1] %in% 
#'                      which(flowSOM.res[[2]] != 7),]
#'    flowCore::write.FCS(ff_tmp,file="ff_tmp.fcs")
#'    # Make an extra file without cluster 7 and double amount of cluster 10
#'    ff_tmp <- ff[c(which(flowSOM.res[[1]]$map$mapping[,1] %in% 
#'                                  which(flowSOM.res[[2]] != 7)),
#'                   which(flowSOM.res[[1]]$map$mapping[,1] %in% 
#'                                  which(flowSOM.res[[2]] == 5))),]
#'    flowCore::write.FCS(ff_tmp,file="ff_tmp2.fcs")
#'    
#'    # Compare the original file with the two new files we made
#'    groupRes <- CountGroups(flowSOM.res[[1]], 
#'                  groups=list("AllCells"=c(fileName),
#'                            "Without_ydTcells"=c("ff_tmp.fcs","ff_tmp2.fcs")))
#'    PlotGroups(flowSOM.res[[1]], groupRes)
#'    
#'    # Compare only the file with the double amount of cluster 10
#'    groupRes <- CountGroups(flowSOM.res[[1]], 
#'                  groups=list("AllCells"=c(fileName),
#'                  "Without_ydTcells"=c("ff_tmp2.fcs")))
#'    PlotGroups(flowSOM.res[[1]], groupRes)
#'
#' @export 
CountGroups <- function (fsom, groups, plot = TRUE, silent = FALSE) 
{
  if (class(groups[[1]]) == "character") {
    files <- unlist(groups)
    counts <- matrix(0, nrow = length(files), ncol = fsom$map$nNodes, 
                     dimnames = list(files, as.character(1:fsom$map$nNodes)))
    for (file in files) {
      if (!silent) {
        print(file)
      }
      ff <- flowCore::read.FCS(file)
      fsom_f <- NewData(fsom, ff)
      if (plot) {
        PlotStars(fsom_f, main = file)
      }
      tmp <- table(fsom_f$map$mapping[, 1])
      counts[file, names(tmp)] <- tmp
    }
    nGroups <- lapply(groups, 
                      length)
  }
  else {
    counts <- do.call(rbind, groups)
    nGroups <- lapply(groups, 
                      nrow)
  }
  pctgs <- t(sapply(seq_len(nrow(counts)), function(i) {
    counts[i, ]/rowSums(counts)[i]
  }))
  means <- apply(pctgs, 2, function(x) {
    tapply(x, INDEX = factor(rep(names(groups), nGroups), 
                             levels = names(groups)), mean)
  })
  means <- means + 0.00000000000000000001
  medians <- apply(pctgs, 2, function(x) {
    tapply(x, INDEX = factor(rep(names(groups), nGroups), 
                             levels = names(groups)), stats::median)
  })
  medians <- medians + 0.00000000000000000001
  means_norm <- list()
  for (group in names(groups)) {
    means_norm[[group]] <- (means[group, ] - min(means))/(max(means) - 
                                                            min(means))
  }
  list(groups = rep(names(groups), unlist(nGroups)), 
       counts = counts, pctgs = pctgs, means = means, medians = medians, 
       means_norm = means_norm)
}


#' Write FlowSOM clustering results to the original FCS files
#'
#' @param  fsom             FlowSOM object as generated by BuildSOM
#' @param  original_files   FCS files that should be extended
#' @param  pp_files         FCS files that correspond to the input of FlowSOM
#' @param  selection_files  Files indicating which cells of the original files
#'                          correspond to the input files
#' @param  silent           If FALSE (default), print some extra output
#'
#' @return Saves the extended fcs file as [originalName]_FlowSOM.fcs
#'
#' @export 
SaveClustersToFCS <- function(fsom, original_files, 
                                pp_files = original_files, 
                                selection_files=NULL, 
                                silent=FALSE){
    for(i in seq_along(original_files)){
        
        if(!silent){message("Extending ",original_files[i],
                            " using the FlowSOM",
                            "mapping of ",pp_files[i], "indexed by ",
                            selection_files[i])}
        
        ff_o <- flowCore::read.FCS(original_files[i])
        ff <- flowCore::read.FCS(pp_files[i])
        if(!is.null(selection_files)){
            s <- unlist(utils::read.table(selection_files[i]))
        } else {
            s <- seq_len(nrow(ff_o))    
        }
        
        # Map the data
        fsom_f <- NewData(fsom,ff)
        
        # Put on corresponding indices
        m <- matrix(0,nrow=nrow(ff_o),ncol=1)
        m[s,] <- fsom_f$map$mapping[,1]
        colnames(m) <- "FlowSOM"
        
        # Save as fcs file
        ff_o <- flowCore::cbind2(ff_o,m)
        flowCore::write.FCS(ff_o,filename=gsub("\\.fcs",
                                                "_FlowSOM.fcs",
                                                original_files[i]))        
    }
}

#' Get cluster label for all individual cells
#'
#' @param  fsom             FlowSOM object as generated by the FlowSOM function
#'                          or the BuildSOM function
#'                          
#' @return vector label for every cell
#'
#' @export 
GetClusters <- function(fsom) {
  if (class(fsom) == "list" & !is.null(fsom$FlowSOM)) {
    fsom <- fsom$FlowSOM
  }
  if (class(fsom) != "FlowSOM") {
    stop("fsom should be a FlowSOM object.")
  }
  return(fsom$map$mapping[,1])
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
#'
#' @export 
GetMetaclusters <- function(fsom, meta = NULL){
  
  if (class(fsom) == "list" & !is.null(fsom$FlowSOM)) {
    if (is.null(meta) & !is.null(fsom$metaclustering)) {
        meta <- fsom$metaclustering
    }
    fsom <- fsom$FlowSOM 
  }
  if (class(fsom) != "FlowSOM"){
    stop("fsom should be a FlowSOM object.")
  } 
  if(is.null(meta)){
    stop("No metaclustering found.")
  }
  
  return(meta[fsom$map$mapping[,1]])
} 

#' Get MFI values for all clusters
#'
#' @param  fsom             FlowSOM object as generated by the FlowSOM function
#'                          or the BuildSOM function
#' @param  colsUsed         logical. Should report only the columns used to 
#'                          build the SOM
#' @param  prettyColnames   logical. Should report pretty column names instead
#'                          of standard column names 
#'                          
#' @return Matrix with median values for each marker
#'
#' @export 
GetMFIs <- function(fsom, colsUsed = FALSE, prettyColnames = FALSE){
  if (class(fsom) == "list" & !is.null(fsom$FlowSOM)) {
    fsom <- fsom$FlowSOM 
  }
  if (class(fsom) != "FlowSOM") {
    stop("fsom should be a FlowSOM object.")
  }
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
#' @export
GetCVs <- function(fsom){
  if (class(fsom) == "list" & !is.null(fsom$FlowSOM)) {
    fsom <- fsom$FlowSOM 
  }
  if (class(fsom) != "FlowSOM") {
    stop("fsom should be a FlowSOM object.")
  }
  return(fsom$map$cvValues)
}

# #' Find peaks and valleys in one-dimensional data
# #' 
# #' @param data                 array containing the data points
# #' @param minDensityThreshold  Only counts peaks which density > threshold
# #' @param ...                  Other parameters to be passed on to density
# #' 
# #' @return A list containing the density result, peaks and valleys
# PeaksAndValleys <- function(data, minDensityThreshold = 0.05,...){
#     dens <- stats::density(data,...)
#     secondDerivative <- diff(sign(diff(dens$y)))
#     peaks <- which(secondDerivative==-2)
#     peaks <- peaks[dens$y[peaks] > minDensityThreshold]
#     
#     # Split the y-values for each peak. The values are assigned to the
#     # next peak, the values higher than the largest peak are assigned 
#     # to group 0
#     tmp <- split(dens$y,
#                 rep(c(peaks,0),
#                     c(peaks[1],diff(peaks),length(dens$y)-max(peaks))))
#     # Find the index of the minimum in each group + the previous peak
#     # = the total index of the valleys
#     valleys <- unlist(lapply(tmp[-c(1,2)],
#                             function(l){which.min(l)[1]}))+
#         peaks[-length(peaks)]
#     list(dens=dens, peaks=peaks, valleys=valleys)
# }


##' Score valleys
##' 
##' Score the valleys of a density function on their width and the height 
##' differences with the surrounding peaks
##' 
##' @param dens    density result, as returned by PeaksAndValleys
##' @param peaks   peaks, as returend by PeaksAndValleys
##' @param valleys valleys, as returned by PeaksAndValleys
##' @param plot    if TRUE, plot the density function, indicating the best
##'                valley in red, the other valleys in transparant red and the
##'                peaks in transparant black
##'  
##' @return scores for all valleys
##' 
##' @export
# ValleyScores <- function(dens, peaks, valleys, plot=FALSE){
#     peakx <- dens$x[peaks]
#     peaky <- dens$y[peaks]
#     valleyx <- dens$x[valleys]
#     valleyy <- dens$y[valleys]
#     
#     valleyscores <- rep(NA,length(valleys))
#     if(length(valleys) > 0){
#         for(i in seq_along(valleys)){
#             height1 <- peaky[i] - valleyy[i]
#             height2 <- peaky[i+1] - valleyy[i]
#             width <- peakx[i+1] - peakx[i]
#             valleyscores[i] <- width * (height1 * height2)
#         }
#     }
#     if(plot){
#         plot(dens,main=max(valleyscores))
#         abline(v=peakx,col="#00000044")
#         abline(v=valleyx,col="#FF000044")
#         abline(v=valleyx[which.max(valleyscores)],col="#FF0000")
#     }
#     return(valleyscores)
# }

# Work in progress
# AssessQuality <- function(fsom, tresh=0.01){
#     scores <- rep(NA, nrow(fsom$map$codes))
#     scores_var <- rep(NA, nrow(fsom$map$codes))
#     for(nodeid in  seq_len(nrow(fsom$map$codes))){
#         node <- fsom$data[fsom$map$mapping[,1] == nodeid, fsom$map$colsUsed]
#         node_scores <- apply(node,2,function(x){
#             pv <- PeaksAndValleys(x)
#             if(length(pv$valleys) > 0){
#                 valleyscores <- ValleyScores(pv$dens, pv$peaks, pv$valleys)
#                 1/(1+sum(valleyscores > tresh))
#             } else {
#                 1
#             }
#         })
#         plot(density(node[,which.min(node_scores)]), main=min(node_scores))
#         scores[nodeid] <- min(node_scores)
#         scores_var[nodeid] <- names(node_scores)[which.min(node_scores)]
#     }
#     print(cbind(scores,scores_var))
#     PlotVariable(UpdateNodeSize(fsom,reset=T),scores,MST=1)
# }
# AssessQuality <- function(fsom, tresh=0.01){
#     scores <- rep(NA, nrow(fsom$map$codes))
#     scores_var <- rep(NA, nrow(fsom$map$codes))
#     for(nodeid in  seq_len(nrow(fsom$map$codes))){
#         node <- fsom$data[fsom$map$mapping[,1] == nodeid, fsom$map$colsUsed]
#         node_scores <- apply(node,2,function(x){
#             pv <- PeaksAndValleys(x,adjust=2)
#             peakValues <- sort(pv$dens$y[pv$peaks],decreasing = T)
#             if(length(pv$valleys) > 0){
#                 peakValues[1] / peakValues[2] 
#             } else {
#                 +Inf#peakValues[1]
#             }
#         })
#         plot(density(node[,which.min(node_scores)]), 
#              main=format(min(node_scores),digits=3),xlim=c(-3,5))
#         scores[nodeid] <- min(node_scores)
#         scores_var[nodeid] <- names(node_scores)[which.min(node_scores)]
#     }
#     print(cbind(scores,scores_var))
#     PlotVariable(UpdateNodeSize(fsom,reset=T),scores,view="MST")
# }
# 
# ValleyScores <- function(dens, peaks, valleys, plot=FALSE){
#     peakx <- dens$x[peaks]
#     peaky <- dens$y[peaks]
#     valleyx <- dens$x[valleys]
#     valleyy <- dens$y[valleys]
#     
#     valleyscores <- rep(NA,length(valleys))
#     if(length(valleys) > 0){
#         for(i in seq_along(valleys)){
#             height1 <- peaky[i] - valleyy[i]
#             height2 <- peaky[i+1] - valleyy[i]
#             width <- peakx[i+1] - peakx[i]
#             valleyscores[i] <- height1 * height2 #+ width
#         }
#     }
#     if(plot){
#         plot(dens,main=max(valleyscores))
#         abline(v=peakx,col="#00000044")
#         abline(v=valleyx,col="#FF000044")
#         abline(v=valleyx[which.max(valleyscores)],col="#FF0000")
#     }
#     return(valleyscores)
# }
