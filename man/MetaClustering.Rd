\name{MetaClustering}
\alias{MetaClustering}
\title{MetaClustering}
\description{Cluster data with automatic number of cluster determination for 
several algorithms}
\usage{MetaClustering(data,method,max=20,...)}
\arguments{
    \item{data}{Matrix containing the data to cluster}
    \item{method}{Clustering method to use}
    \item{max}{Maximum number of clusters to try out}
    \item{...}{Extra parameters to pass along}
}
\value{Numeric array indicating cluster for each datapoint}
\seealso{
    \code{\link{metaClustering_consensus}}
}
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
    
    # Get metaclustering per cell
    flowSOM.clustering <- metacl[flowSOM.res$map$mapping[,1]]    
}