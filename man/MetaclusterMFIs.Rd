% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/6_deprecated.R
\name{MetaclusterMFIs}
\alias{MetaclusterMFIs}
\title{MetaclusterMFIs}
\usage{
MetaclusterMFIs(fsom)
}
\arguments{
\item{fsom}{Result of calling the FlowSOM function}
}
\value{
Metacluster MFIs
}
\description{
Compute the median fluorescence intensities for the metaclusters
}
\examples{
fileName <- system.file("extdata", "68983.fcs", package="FlowSOM")
ff <- flowCore::read.FCS(fileName)
ff <- flowCore::compensate(ff,ff@description$SPILL)
ff <- flowCore::transform(ff,
         flowCore::transformList(colnames(ff@description$SPILL),
                                flowCore::logicleTransform()))
flowSOM.res <- FlowSOM(ff,scale=TRUE,colsToUse=c(9,12,14:18),maxMeta=10)
mfis <- GetMetaclusterMFIs(flowSOM.res)
}
