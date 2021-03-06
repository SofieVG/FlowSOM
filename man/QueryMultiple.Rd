% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5_plotFunctions.R
\name{QueryMultiple}
\alias{QueryMultiple}
\title{QueryMultiple}
\usage{
QueryMultiple(fsom, cellTypes, plotFile = "queryMultiple.pdf", ...)
}
\arguments{
\item{fsom}{FlowSOM object}

\item{cellTypes}{Description of the cell types. Named list, with one named
vector per cell type containing "high"/"low" values}

\item{plotFile}{Path to a pdf file to save the plots (default is 
queryMultiple.pdf). If \code{NULL}, no plots will be 
generated}

\item{...}{Additional arguments to pass to \code{\link{QueryStarPlot}}}
}
\value{
A label for every FlowSOM cluster (Unknown or one of the celltype
        names of the list, if selected by QueryStarPlot)
}
\description{
Function which takes a named list of multiple cell types, where every item is
a named vector with values "high"/"low" and the names correspond to the
markers or channels (e.g. as generated by parse_markertable).
}
\examples{
   file <- system.file("extdata", "68983.fcs", package = "FlowSOM")
   ff <- flowCore::read.FCS(file)
   # Use the wrapper function to build a flowSOM object (saved in flowSOM.res)
   # and a metaclustering (saved in flowSOM.res[["metaclustering"]])
   flowSOM.res <- FlowSOM(ff, compensate = TRUE, transform = TRUE, scale = TRUE,
                  colsToUse = c(9, 12, 14:18), nClus = 10, silent = FALSE,
                  xdim = 7, ydim = 7)
   cellTypes <- list("CD8 T cells" = c("PE-Cy7-A" = "high",
                                        "APC-Cy7-A" = "high",
                                        "Pacific Blue-A" = "high"),
                       "B cells" = c("PE-Cy5-A" = "high"),
                       "NK cells" = c("PE-A" = "high",
                                      "PE-Cy7-A" = "low",
                                      "APC-Cy7-A" = "low"))
   query_res <- QueryMultiple(flowSOM.res, cellTypes, "query_multiple.pdf")
   
}
