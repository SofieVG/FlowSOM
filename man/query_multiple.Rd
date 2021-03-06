% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/6_deprecated.R
\name{query_multiple}
\alias{query_multiple}
\title{query_multiple}
\usage{
query_multiple(fsom, cell_types, pdf_name = "query_multiple.pdf", ...)
}
\arguments{
\item{fsom}{FlowSOM object}

\item{cell_types}{Description of the cell types. Named list, with one named
vector per cell type containing "high"/"low" values}

\item{pdf_name}{Path to a pdf file to save figures}

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
   file <- system.file("extdata", "68983.fcs", package="FlowSOM")
   ff <- flowCore::read.FCS(file)
   # Use the wrapper function to build a flowSOM object (saved in flowSOM.res)
   # and a metaclustering (saved in flowSOM.res[["metaclustering"]])
   flowSOM.res <- FlowSOM(ff,compensate = TRUE, transform = TRUE,scale = TRUE,
                  colsToUse = c(9,12,14:18), nClus = 10, silent = FALSE,
                  xdim=7, ydim=7)
   cell_types <- list("CD8 T cells" = c("PE-Cy7-A" = "high",
                                        "APC-Cy7-A" = "high",
                                        "Pacific Blue-A" = "high"),
                       "B cells" = c("PE-Cy5-A" = "high"),
                       "NK cells" = c("PE-A" = "high",
                                      "PE-Cy7-A" = "low",
                                      "APC-Cy7-A" = "low"))
   query_res <- QueryMultiple(flowSOM.res, cell_types, "query_multiple.pdf")
   
}
\seealso{
\code{\link{QueryStarPlot}}
}
