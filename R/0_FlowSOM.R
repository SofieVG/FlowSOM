# FlowSOM object
# List containing the following
# 
# after ReadInput:
#     data: matrix containing all the concatenated data files
#     metaData: a list, containing start and end indices for each file
#     compensate: logical, is the data compensated
#     spillover: spillover matrix the data is compensated with
#     transform: logical, is the data transformed with a logicle transform
#     toTransform: column names or indices are transformed
#     scale: logical, is the data rescaled
#     scaled.center: parameter used to rescale
#     scaled.scale: parameter used to rescale


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
#'                      If NULL and compensate=TRUE, we will look for $SPILL 
#'                      description in fcs file.
#' @param transform     logical, does the data need to be transformed with a
#'                      logicle transform
#' @param toTransform   column names or indices that need to be transformed.
#'                      If \code{NULL} and transform = \code{TRUE}, column
#'                      names of \code{$SPILL} description in fcs file will
#'                      be used.
#' @param transformFunction Defaults to logicleTransform()
#' @param scale         logical, does the data needs to be rescaled
#' @param scaled.center see \code{\link{scale}}
#' @param scaled.scale  see \code{\link{scale}}
#' @param silent        if \code{TRUE}, no progress updates will be printed
#' @param colsToUse     column names or indices to use for building the SOM
#' @param importance    array with numeric values. Parameters will be scaled 
#'                      according to importance
#' @param nClus         Exact number of clusters for meta-clustering. 
#'                      If \code{NULL}, several options will be tried 
#'                      (\code{1:maxMeta})
#' @param maxMeta       Maximum number of clusters to try out for 
#'                      meta-clustering. Ignored if nClus is specified
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
#' @seealso \code{\link{scale}},\code{\link{ReadInput}},\code{\link{BuildSOM}},
#'          \code{\link{BuildMST}},\code{\link{MetaClustering}}
#' @examples
#' # Read from file
#' fileName <- system.file("extdata","lymphocytes.fcs",package="FlowSOM")
#' flowSOM.res <- FlowSOM(fileName, compensate=TRUE,transform=TRUE,
#'                       scale=TRUE,colsToUse=c(9,12,14:18),maxMeta=10)
#' # Or read from flowFrame object
#' ff <- read.FCS(fileName)
#' ff <- compensate(ff,ff@@description$SPILL)
#' ff <- transform(ff,transformList(colnames(ff@@description$SPILL),
#'                                 logicleTransform()))
#' flowSOM.res <- FlowSOM(ff,scale=TRUE,colsToUse=c(9,12,14:18),maxMeta=10)
#' 
#' # Plot results
#' PlotStars(flowSOM.res[[1]])
#' 
#' # Get metaclustering per cell
#' flowSOM.clustering <- flowSOM.res[[2]][flowSOM.res[[1]]$map$mapping[,1]]
#' 
#' 
#' 
#' @importFrom flowCore read.FCS compensate transform logicleTransform exprs 
#'             transformList write.FCS 'exprs<-'
#' @importFrom igraph graph.adjacency minimum.spanning.tree layout.kamada.kawai
#'             plot.igraph add.vertex.shape get.edges shortest.paths E V 'V<-'
#'             igraph.shape.noclip
#' @importFrom tsne tsne
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus
#' @importFrom BiocGenerics colnames
#' @export
FlowSOM <- function(input, pattern=".fcs", compensate=FALSE, spillover=NULL, 
                    transform=FALSE, toTransform=NULL, 
                    transformFunction=flowCore::logicleTransform(), scale=TRUE, 
                    scaled.center=TRUE, scaled.scale=TRUE, silent=TRUE, 
                    colsToUse, nClus=NULL, maxMeta, importance=NULL, ...){
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
    
    t <- system.time(fsom <- ReadInput(input, pattern=pattern, 
                                        compensate=compensate, 
                                        spillover=spillover, 
                                        transform=transform, 
                                        toTransform=toTransform, 
                                        transformFunction = transformFunction, 
                                        scale=scale,
                                        scaled.center=scaled.center, 
                                        scaled.scale=scaled.scale, 
                                        silent=silent))
    if(!silent) message(t[3],"\n")
    t <- system.time(fsom <- BuildSOM(fsom, colsToUse, silent=silent, 
                                        importance=importance, ...))
    if(!silent) message(t[3],"\n")
    t <- system.time(fsom <- BuildMST(fsom, silent=silent))
    if(!silent) message(t[3],"\n")
    if(is.null(nClus)){
        t <- system.time(cl <- MetaClustering(fsom$map$codes,
                                        "metaClustering_consensus", maxMeta))
    } else {
        t <- system.time(cl <- metaClustering_consensus(fsom$map$codes, nClus))
    }
    if(!silent) message(t[3],"\n")
    list(fsom, cl)
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
#' @param writeOutput Whether to write the resulting flowframe to a file
#' @param outputFile  Full path to output file
#' @param writeMeta   If TRUE, files with the indices of the selected cells are
#'                    generated
#'                  
#' @return This function does not return anything, but will write a file with
#'         about \code{cTotal} cells to \code{outputFile}
#'
#' @seealso \code{\link{ceiling}}
#'
#' @examples
#' # Define filename
#' fileName <- system.file("extdata","lymphocytes.fcs",package="FlowSOM")
#' # This example will sample 2 times 500 cells.
#' ff_new <- AggregateFlowFrames(c(fileName,fileName),1000)
#' 
#' @export
AggregateFlowFrames <- function(fileNames, cTotal,
                            writeOutput = FALSE, outputFile="aggregate.fcs", 
                            writeMeta=FALSE){
    
    nFiles <- length(fileNames)
    cFile <- ceiling(cTotal/nFiles)
    
    flowFrame <- NULL
    
    for(i in seq_len(nFiles)){
        f <- flowCore::read.FCS(fileNames[i])
        c <- sample(seq_len(nrow(f)),min(nrow(f),cFile))
        if(writeMeta){
            #<path_to_outputfile>/<filename>_selected_<outputfile>.txt
            write.table(c,paste(gsub("[^/]*$","",outputFile),
                            gsub("\\.[^.]*$","",gsub(".*/","",fileNames[i])),
                            "_selected_",
                            gsub("\\.[^.]*$","",gsub(".*/","",outputFile)),
                            ".txt",sep=""))
        }
        m <- matrix(rep(i,min(nrow(f),cFile)))
        m2 <- m + rnorm(length(m),0,0.1)
        m <- cbind(m,m2)
        colnames(m) <- c("File","File_scattered")
        f <- cbind2(f[c,],m)
        if(is.null(flowFrame)){
            flowFrame <- f
            flowFrame@description$`$FIL` <- gsub(".*/","",outputFile)
            flowFrame@description$`FILENAME` <- gsub(".*/","",outputFile)
        }
        else {
            flowCore::exprs(flowFrame) <- rbind(flowCore::exprs(flowFrame), 
                                        flowCore::exprs(f))
        }
    }
    
    flowFrame@description[[
        paste("flowCore_$P",ncol(flowFrame)-1,"Rmin",sep="")]] <- 0
    flowFrame@description[[
        paste("flowCore_$P",ncol(flowFrame)-1,"Rmax",sep="")]] <- nFiles+1
    flowFrame@description[[
        paste("flowCore_$P",ncol(flowFrame),"Rmin",sep="")]] <- 0
    flowFrame@description[[
        paste("flowCore_$P",ncol(flowFrame),"Rmax",sep="")]] <- nFiles+1  
    
    flowFrame@description$FIL <- gsub(".*/","",outputFile)
    if(writeOutput){
        flowCore::write.FCS(flowFrame,filename=outputFile)
    }
    
    flowFrame
}