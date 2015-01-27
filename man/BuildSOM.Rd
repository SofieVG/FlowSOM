\name{BuildSOM}
\alias{BuildSOM}
\title{Build a self-organizing map}
\description{Build a SOM based on the data contained in the FlowSOM object}
\usage{BuildSOM(fsom, colsToUse=NULL, silent=FALSE, ...)}
\arguments{
    \item{fsom}{FlowSOM object containing the data, as constructed by 
    the \code{\link{ReadInput}} function}
    \item{colsToUse}{column names or indices to use for building the SOM}
    \item{silent}{if \code{TRUE}, no progress updates will be printed}
    \item{...}{options to pass on to the SOM function (xdim, ydim, rlen, mst, 
    alpha, radius, init, distf)}
}
\value{FlowSOM object containing the SOM result, which can be used as input
for the \code{\link{BuildMST}} function}
\seealso{
    \code{\link{ReadInput}},\code{\link{BuildMST}}
}
\references{This code is strongly based on the \code{kohonen} package.
R. Wehrens and L.M.C. Buydens, Self- and Super-organising Maps in R: 
the kohonen package J. Stat. Softw., 21(5), 2007}
\examples{
    # Read from file
    fileName <- system.file("extdata","lymphocytes.fcs",package="FlowSOM")
    flowSOM.res <- ReadInput(fileName, compensate=TRUE,transform=TRUE,
                            scale=TRUE)
    
    # Build the Self-Organizing Map
    # E.g. with gridsize 5x5, presenting the dataset 20 times, 
    # no use of MST in neighbourhood calculations in between
    flowSOM.res <- BuildSOM(flowSOM.res,colsToUse=c(9,12,14:18),
                            xdim=5,ydim=5,rlen=20)
    
    # Build the minimal spanning tree and apply metaclustering
    flowSOM.res <- BuildMST(flowSOM.res)
    metacl <- MetaClustering(flowSOM.res$map$codes,
                            "metaClustering_consensus",max=10)
}