\name{metaClustering_consensus}
\alias{metaClustering_consensus}
\title{MetaClustering}
\description{Cluster data using hierarchical consensus clustering 
with k clusters}
\usage{metaClustering_consensus(data, k=7)}
\arguments{
    \item{data}{Matrix containing the data to cluster}
    \item{k}{Number of clusters}
}
\value{Numeric array indicating cluster for each datapoint}
\seealso{
    \code{\link{MetaClustering}}
}
\examples{
    # Read from file, build self-organizing map and minimal spanning tree
    fileName <- system.file("extdata","lymphocytes.fcs",package="FlowSOM")
    flowSOM.res <- ReadInput(fileName, compensate=TRUE,transform=TRUE,
                            scale=TRUE)
    flowSOM.res <- BuildSOM(flowSOM.res,colsToUse=c(9,12,14:18))
    flowSOM.res <- BuildMST(flowSOM.res)
    
    # Apply consensus metaclustering
    metacl <- metaClustering_consensus(flowSOM.res$map$codes,k=10)    
}