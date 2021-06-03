#' Read fcs-files or flowframes
#' 
#' Take some input and return FlowSOM object containing a matrix with 
#' the preprocessed data (compensated, transformed, scaled)
#'
#' @param input         a flowFrame, a flowSet or an array of paths to files 
#'                      or directories
#' @param pattern       if input is an array of file- or directorynames, 
#'                      select only files containing pattern
#' @param compensate    logical, does the data need to be compensated
#' @param spillover     spillover matrix to compensate with
#'                      If \code{NULL} and compensate = \code{TRUE}, we will
#'                      look for \code{$SPILL} description in fcs file.
#' @param transform     logical, does the data need to be transformed
#' @param toTransform   column names or indices that need to be transformed.
#'                      Will be ignored if \code{transformList} is given.
#'                      If \code{NULL} and transform = \code{TRUE}, column names
#'                      of \code{$SPILL} description in fcs file will be used.
#' @param transformFunction Defaults to logicleTransform()
#' @param transformList transformList to apply on the samples.
#' @param scale         logical, does the data needs to be rescaled
#' @param scaled.center see \code{\link{scale}}
#' @param scaled.scale  see \code{\link{scale}}
#' @param silent        if \code{TRUE}, no progress updates will be printed. 
#'                      Default = \code{FALSE}
#'
#' @return FlowSOM object containing the data, which can be used as input
#' for the BuildSOM function
#'
#' @seealso \code{\link{scale}}, \code{\link{BuildSOM}}
#' 
#' @examples
#' # Read from file
#' fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
#' flowSOM.res <- ReadInput(fileName, compensate = TRUE, transform = TRUE,
#'                          scale = TRUE)
#' 
#' # Or read from flowFrame object
#' ff <- flowCore::read.FCS(fileName)
#' ff <- flowCore::compensate(ff, flowCore::keyword(ff)[["SPILL"]])
#' ff <- flowCore::transform(ff,
#'          flowCore::transformList(colnames(flowCore::keyword(ff)[["SPILL"]]),
#'                                 flowCore::logicleTransform()))
#' flowSOM.res <- ReadInput(ff, scale = TRUE)
#' 
#' # Build the self-organizing map and the minimal spanning tree
#' flowSOM.res <- BuildSOM(flowSOM.res, colsToUse = c(9, 12, 14:18))
#' flowSOM.res <- BuildMST(flowSOM.res)
#' 
#' # Apply metaclustering
#' metacl <- MetaClustering(flowSOM.res$map$codes,
#'                          "metaClustering_consensus", max = 10)
#' 
#' # Get metaclustering per cell
#' flowSOM.clustering <- metacl[flowSOM.res$map$mapping[, 1]]    
#'  
#' @export
ReadInput <- function(input, pattern = ".fcs", 
                      compensate = FALSE, 
                      spillover = NULL,
                      transform = FALSE, 
                      toTransform = NULL, 
                      transformFunction = flowCore::logicleTransform(),
                      transformList = NULL,
                      scale = FALSE, 
                      scaled.center = TRUE, 
                      scaled.scale = TRUE, 
                      silent = FALSE){
  
  fsom <- list(pattern = pattern, 
               compensate = compensate, 
               spillover = spillover,
               transform = transform, 
               toTransform = toTransform,
               transformFunction = transformFunction, 
               transformList = transformList,
               scale = scale)
  class(fsom) <- "FlowSOM"
  
  if(is(input, "flowFrame")){
    fsom <- AddFlowFrame(fsom, input)        
  } else if(is(input, "flowSet")){
    for(i in seq_along(input)){
      fsom <- AddFlowFrame(fsom, input[[i]])    
    }
  } else if(is(input, "matrix")){
    flowFrame <- flowCore::flowFrame(input)
    fsom <- AddFlowFrame(fsom, flowFrame)
  } else if(is(input, "character")){    
    # Replace all directories by the files they contain
    toAdd <- NULL
    toRemove <- NULL
    for(i in seq_along(input)){
      if(file.info(input[i])$isdir){
        toAdd <- c(toAdd, list.files(input[i], pattern = pattern,
                                     full.names = TRUE))
        toRemove <- c(toRemove, i)
      }
    }
    if(!is.null(toRemove)){
      input <- c(input[-toRemove], toAdd)
    }
    
    # Select all files corresponding to the pattern
    input <- grep(pattern, input, value = TRUE)
    
    # Read all files
    if(length(input) > 0){
      for(i in seq_along(input)){
        if(file.exists(input[i])){
          if(!silent) message("Reading file ", input[i], "\n")
          if (tools::file_ext(input[i]) == "csv") {
            flowFrame <- flowFrame(utils::read.table(input[i]))
          } else { #else if(tools::file_ext(input[i]) == "fcs"){
            flowFrame <- suppressWarnings(
              flowCore::read.FCS(input[i]))
          }
          fsom <- AddFlowFrame(fsom, flowFrame)
        }
      }
    } else {
      stop("No files containing the pattern are found.")
    }
  } else {
    stop(paste("Inputs of type", class(input), "are not supported. 
                    Please supply either a FlowFrame, a FlowSet or an array
                    of valid paths to files or directories."))
  }
  
  if(scale){
    if(!silent) message("Scaling the data\n")
    if(!all(is.null(names(scaled.center)))){
      cols_to_scale <- intersect(colnames(fsom$data), names(scaled.center))
      sc <- scale(x = fsom$data[, cols_to_scale], 
                  center = scaled.center[cols_to_scale], 
                  scale = scaled.scale[cols_to_scale])
      fsom$data[, cols_to_scale] <- sc
      fsom$scaled.center <- attr(sc, "scaled:center")
      fsom$scaled.scale <- attr(sc, "scaled:scale")
    } else {
      fsom$data <- scale(x = fsom$data, center = TRUE, scale = TRUE)
      
      fsom$scaled.center <- attr(fsom$data, "scaled:center")
      attr(fsom$data, "scaled:center") <- NULL
      fsom$scaled.scale <- attr(fsom$data, "scaled:scale") 
      attr(fsom$data, "scaled:scale") <- NULL
    }
    
  }
  
  return(fsom)
}

#' Add a flowFrame to the data variable of the FlowSOM object
#'
#' @param fsom      FlowSOM object, as constructed by the ReadInput function
#' @param flowFrame flowFrame to add to the FlowSOM object
#'
#' @return FlowSOM object with data added
#'
#' @seealso \code{\link{ReadInput}}
AddFlowFrame <- function(fsom, flowFrame){
  # Compensation
  if(fsom$compensate){
    if(is.null(fsom$spillover)){
      if(!is.null(flowFrame@description$SPILL)){
        spillover <- flowFrame@description$SPILL    
      } else if (!is.null(flowFrame@description$`$SPILLOVER`)){
        if(is(flowFrame@description$`$SPILLOVER`, "matrix")){
          spillover <- flowFrame@description$`$SPILLOVER`
          flowFrame@description$SPILL <- spillover
        } else {
          spilloverStr <- strsplit(
            flowFrame@description$`$SPILLOVER`,
            ", ")[[1]]
          n <- as.numeric(spilloverStr[1])
          spillover <- t(
            matrix(as.numeric(spilloverStr[(n+2):length(spilloverStr)]), 
                   ncol = n))
          colnames(spillover) <- spilloverStr[2:(n+1)]
          flowFrame@description$SPILL <- spillover
        }
      } else {
        stop("No compensation matrix found")
      }
    } else {
      spillover <- fsom$spillover
    }
    flowFrame <- flowCore::compensate(flowFrame, spillover)
  }
  
  # Transform
  if(fsom$transform){
    if(!is.null(fsom$transformList)){
      
      flowFrame <- flowCore::transform(flowFrame, 
                                       fsom$transformList)
    } else {
      if(is.null(fsom$toTransform)){ 
        fsom$toTransform <- colnames(flowFrame@description$SPILL)
      } else{ 
        fsom$toTransform <- 
          colnames(flowCore::exprs(flowFrame)[, fsom$toTransform])
      }
      flowFrame <- 
        flowCore::transform(flowFrame, 
                            flowCore::transformList(fsom$toTransform,
                                                    fsom$transformFunction))
    }
  }
  
  # Save pretty names for nicer visualisation later on
  if(is.null(fsom$prettyColnames)){
    n <- flowCore::parameters(flowFrame)@data[, "name"]
    d <- flowCore::parameters(flowFrame)@data[, "desc"]
    d[is.na(d)] <- n[is.na(d)]
    if(any(grepl("#", d))){
      # Support for hashtag notation: 
      # antibody#fluorochrome -> antibody (fluorochrome)
      fsom$prettyColnames <- gsub("#(.*)$", " <\\1>", d)
    } else {
      fsom$prettyColnames <- paste(d, " <", n, ">", sep = "")
    }
    names(fsom$prettyColnames) <- colnames(flowCore::exprs(flowFrame))
  }
  
  # Add the data to the matrix
  f <- flowCore::exprs(flowFrame)
  attr(f, "ranges") <- NULL
  name <- flowCore::keyword(flowFrame,"FIL")[[1]]
  if(is.null(name)) name <- flowFrame@description$`$FIL`
  if(is.null(name)) name <- length(fsom$metaData)+1
  if(is.null(fsom$data)){ 
    fsom$data <- f 
    fsom$metaData <- list()
    fsom$metaData[[name]] <- c(1, nrow(fsom$data))
  } else {
    fsom$data <- rbind(fsom$data, f)
    fsom$metaData[[name]] <- c(nrow(fsom$data) - nrow(f) + 1, 
                               nrow(fsom$data))
  }
  
  return(fsom)
}