% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_buildSOM.R
\name{GetFeatures}
\alias{GetFeatures}
\title{GetFeatures}
\usage{
GetFeatures(
  fsom,
  files,
  level = c("clusters", "metaclusters"),
  type = "counts",
  MFI = NULL,
  filenames = NULL,
  silent = FALSE
)
}
\arguments{
\item{fsom}{FlowSOM object as generated by the FlowSOM function
or the BuildSOM function}

\item{files}{Either a vector of fcs files or paths to fcs files}

\item{level}{Level(s) of interest. Default is c("clusters",
"metaclusters"), but can also be only one of them}

\item{type}{Type of features to extract. Default is "counts", 
can be a vector of "counts", "percentages" and/or 
"MFIs"}

\item{MFI}{Vector with channels / markers for which the MFI 
values must be returned when "MFIs" is in \code{type}}

\item{filenames}{An optional vector with filenames that will be used
as rownames in the count matrices. If NULL (default)
either the paths will be used or a numerical vector.}

\item{silent}{Logical. If \code{TRUE}, print progress messages.
Default = \code{FALSE}.}
}
\value{
matrix with features per population - type combination
}
\description{
Map fcs files on an existing FlowSOM object
}
\examples{
 # Build FlowSom result
 fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
 ff <- flowCore::read.FCS(fileName)
 ff <- flowCore::compensate(ff, flowCore::keyword(ff)[["SPILL"]])
 ff <- flowCore::transform(ff,
         flowCore::transformList(colnames(flowCore::keyword(ff)[["SPILL"]]),
                                flowCore::logicleTransform()))
 flowSOM.res <- FlowSOM(ff[1:1000, ], 
                        scale = TRUE, 
                        colsToUse = c(9, 12, 14:18),
                        nClus = 10)
   
 # Map new data
 counts <- GetFeatures(fsom = flowSOM.res, 
                       level = "clusters",
                       files = c(ff[1001:2000, ], ff[2001:3000, ]))
 features <- GetFeatures(fsom = flowSOM.res, 
                         files = c(ff[1001:2000, ], ff[2001:3000, ]),
                         type = c("counts", "percentages", "MFIs"), 
                         MFI = "APC-A", 
                         filenames = c("ff_1001-2000", "ff_2001-3000"))

}
