FlowSOM <- function(input, pattern=".fcs", compensate=FALSE, spillover=NULL, 
                    transform=FALSE, toTransform=NULL, scale=TRUE, 
                    scaled.center=TRUE, scaled.scale=TRUE, silent=TRUE, 
                    colsToUse, nClus=NULL, maxMeta,...){
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
                    compensate=compensate, spillover=spillover, 
                    transform=transform, toTransform=toTransform, scale=scale,
                    scaled.center=scaled.center, scaled.scale=scaled.scale, 
                    silent=silent))
    if(!silent) message(t[3],"\n")
    t <- system.time(fsom <- BuildSOM(fsom, colsToUse, silent=silent, ...))
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

AggregateFlowFrames <- function(fileNames, outputFile, cTotal,writeMeta=F){
    
    nFiles <- length(fileNames)
    cFile <- ceiling(cTotal/nFiles)
    
    flowFrame <- NULL
    
    for(i in seq_len(nFiles)){
        f <- read.FCS(fileNames[i])
        c <- sample(1:nrow(f),min(nrow(f),cFile))
        if(writeMeta){
            write.table(c,paste(gsub("[^/]*$","",outputFile),
                                gsub("\\.[^.]*$","",gsub(".*/","",fileNames[i])),
                                "_selected.txt"))
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
            exprs(flowFrame) <- rbind(exprs(flowFrame), exprs(f))
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
    
    write.FCS(flowFrame,filename=outputFile)
}