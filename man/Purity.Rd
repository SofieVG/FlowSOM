\name{Purity}
\alias{Purity}
\title{Purity of a clustering result}
\description{Compute mean weighted purity of a clustering result in comparison
with true cluster labels}
\usage{Purity(realClusters, predictedClusters, weighted=TRUE)}
\arguments{
    \item{realClusters}{array with real cluster values}
    \item{predictedClusters}{array with predicted cluster values}
    \item{weighted}{logical. Should the mean be weighted depending on the
    number of points in the predicted clusters}
}
\value{Weighted mean purity value}
\seealso{\code{\link{FMeasure}}}
\examples{
    # Read from file, build self-organizing map and minimal spanning tree
    fileName <- system.file("extdata","lymphocytes.fcs",package="FlowSOM")
    flowSOM.res <- ReadInput(fileName, compensate=TRUE,transform=TRUE,
                            scale=TRUE)
    flowSOM.res <- BuildSOM(flowSOM.res,colsToUse=c(9,12,14:18))
    flowSOM.res <- BuildMST(flowSOM.res)
    
    # Apply metaclustering
    metacl <- MetaClustering(flowSOM.res$map$codes,"metaClustering_consensus",
                            max=10)
    
    # Get the manually gated labels using a gatingML file 
    # and the flowUtils library
    ff <- read.FCS(fileName)
    ff_c <- compensate(ff,ff@description$SPILL)
    colnames(ff_c)[8:18] <- paste("Comp-",colnames(ff_c)[8:18],sep="")
    library(flowUtils)
    flowEnv <- new.env()
    gatingFile <- system.file("extdata","manualGating.xml", 
                            package="FlowSOM")
    read.gatingML(gatingFile, flowEnv) 
    filterList <- list( "B cells" = flowEnv$ID52300206,
                    "ab T cells" = flowEnv$ID785879196,
                    "yd T cells" = flowEnv$ID188379411,
                    "NK cells" = flowEnv$ID1229333490,
                    "NKT cells" = flowEnv$ID275096433
                    )
    results <- list()
    for(cellType in names(filterList)){
        results[[cellType]] <- filter(ff_c,filterList[[cellType]])@subSet
    }
    manual <- rep("Unknown",nrow(ff))
    for(celltype in names(results)){
        manual[results[[celltype]]] <- celltype
    }
    
    # Test the purity of the nodes
    Purity(manual,metacl[flowSOM.res$map$mapping[,1]])
}