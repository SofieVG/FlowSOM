\name{ReadInput}
\alias{ReadInput}
\title{Read fcs-files or flowframes}
\description{Take some input and return a matrix with preprocessed data 
(compensated, transformed, scaled)}
\usage{ReadInput(input, pattern=".fcs", compensate=FALSE, spillover=NULL, 
                transform=FALSE, toTransform=NULL, scale=FALSE, 
                scaled.center=TRUE, scaled.scale=TRUE, silent=FALSE)}
\arguments{
    \item{input}{a flowFrame, a flowSet or an array of paths to files 
    or directories}
    \item{pattern}{if input is an array of file- or directorynames, 
    select only files containing pattern}
    \item{compensate}{logical, does the data need to be compensated}
    \item{spillover}{spillover matrix to compensate with
                        If \code{NULL} and compensate=\code{TRUE}, we will
                        look for \code{$SPILL} description in fcs file.}
    \item{transform}{logical, does the data need to be transformed with 
    a logicle transform}
    \item{toTransform}{column names or indices that need to be transformed.
                        If \code{NULL} and transform=\code{TRUE}, column names
                        of \code{$SPILL} description in fcs file will be used.}
    \item{scale}{logical, does the data needs to be rescaled}
    \item{scaled.center}{see \code{\link{scale}}}
    \item{scaled.scale}{see \code{\link{scale}}}
    \item{silent}{if \code{TRUE}, no progress updates will be printed}
}
\value{FlowSOM object containing the data, which can be used as input for the
BuildSOM function}

\seealso{
    \code{\link{scale}},\code{\link{BuildSOM}}
}

\examples{
    # Read from file
    fileName <- system.file("extdata","lymphocytes.fcs",package="FlowSOM")
    flowSOM.res <- ReadInput(fileName, compensate=TRUE,transform=TRUE,
                            scale=TRUE)
    
    # Or read from flowFrame object
    ff <- read.FCS(fileName)
    ff <- compensate(ff,ff@description$SPILL)
    ff <- transform(ff,transformList(colnames(ff@description$SPILL),
                    logicleTransform()))
    flowSOM.res <- ReadInput(ff,scale=TRUE)
    
    # Build the self-organizing map and the minimal spanning tree
    flowSOM.res <- BuildSOM(flowSOM.res,colsToUse=c(9,12,14:18))
    flowSOM.res <- BuildMST(flowSOM.res)
    
    # Apply metaclustering
    metacl <- MetaClustering(flowSOM.res$map$codes,
                            "metaClustering_consensus",max=10)
    
    # Get metaclustering per cell
    flowSOM.clustering <- metacl[flowSOM.res$map$mapping[,1]]    
}