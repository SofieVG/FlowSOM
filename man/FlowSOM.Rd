\name{FlowSOM}
\alias{FlowSOM}
\title{Run the FlowSOM algorithm}
\description{Method to run general FlowSOM workflow. Will scale the data and 
uses consensus meta-clustering by default.}
\usage{FlowSOM(input, pattern=".fcs", compensate=FALSE, spillover=NULL,
                transform=FALSE, toTransform=NULL, scale=TRUE, 
                scaled.center=TRUE, scaled.scale=TRUE, silent=TRUE, colsToUse,
                nClus=NULL, maxMeta,...)}
\arguments{
    \item{input}{a flowFrame, a flowSet or an array of paths to files or 
    directories}
    \item{pattern}{if input is an array of file- or directorynames, select 
    only files containing pattern}
    \item{compensate}{logical, does the data need to be compensated}
    \item{spillover}{spillover matrix to compensate with
                    If NULL and compensate=TRUE, we will look for $SPILL 
                    description in fcs file.}
    \item{transform}{logical, does the data need to be transformed with a
    logicle transform}
    \item{toTransform}{column names or indices that need to be transformed.
                        If \code{NULL} and transform = \code{TRUE}, column
                        names of \code{$SPILL} description in fcs file will
                        be used.}
    \item{scale}{logical, does the data needs to be rescaled}
    \item{scaled.center}{see \code{\link{scale}}}
    \item{scaled.scale}{see \code{\link{scale}}}
    \item{silent}{if \code{TRUE}, no progress updates will be printed}
    \item{colsToUse}{column names or indices to use for building the SOM}
    \item{nClus}{Exact number of clusters for meta-clustering. If \code{NULL},
    several options will be tried (\code{1:maxMeta})}
    \item{maxMeta}{Maximum number of clusters to try out for meta-clustering. 
    Ignored if nClus is specified}
    \item{...}{options to pass on to the SOM function (xdim, ydim, rlen, mst, 
    alpha, radius, init, distf)}
}
\value{
    A \code{list} with two items: the first is the flowSOM object containing
    all information (see the vignette for more detailed information about this 
    object), the second is the metaclustering of the nodes of the grid. This
    is a wrapper function for \code{\link{ReadInput}}, \code{\link{BuildSOM}},
    \code{\link{BuildMST}} and \code{\link{MetaClustering}}. Executing them
    separately may provide more options.
}

\seealso{
    \code{\link{scale}},\code{\link{ReadInput}},\code{\link{BuildSOM}},
    \code{\link{BuildMST}},\code{\link{MetaClustering}}
}

\examples{
    # Read from file
    fileName <- system.file("extdata","lymphocytes.fcs",package="FlowSOM")
    flowSOM.res <- FlowSOM(fileName, compensate=TRUE,transform=TRUE,
                    scale=TRUE,colsToUse=c(9,12,14:18),maxMeta=10)
    
    # Or read from flowFrame object
    ff <- read.FCS(fileName)
    ff <- compensate(ff,ff@description$SPILL)
    ff <- transform(ff,transformList(colnames(ff@description$SPILL),
                    logicleTransform()))
    flowSOM.res <- FlowSOM(ff,scale=TRUE,colsToUse=c(9,12,14:18),maxMeta=10)
    
    # Get metaclustering per cell
    flowSOM.clustering <- flowSOM.res[[2]][flowSOM.res[[1]]$map$mapping[,1]]
}