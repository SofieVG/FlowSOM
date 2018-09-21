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
#'                      If \code{NULL} and compensate=\code{TRUE}, we will
#'                      look for \code{$SPILL} description in fcs file.
#' @param transform     logical, does the data need to be transformed
#' @param toTransform   column names or indices that need to be transformed.
#'                      If \code{NULL} and transform=\code{TRUE}, column names
#'                      of \code{$SPILL} description in fcs file will be used.
#' @param transformFunction Defaults to logicleTransform()
#' @param scale         logical, does the data needs to be rescaled
#' @param scaled.center see \code{\link{scale}}
#' @param scaled.scale  see \code{\link{scale}}
#' @param silent        if \code{TRUE}, no progress updates will be printed
#'
#' @return FlowSOM object containing the data, which can be used as input
#' for the BuildSOM function
#'
#' @seealso \code{\link{scale}},\code{\link{BuildSOM}}
#' 
#' @examples
#' # Read from file
#' fileName <- system.file("extdata","lymphocytes.fcs",package="FlowSOM")
#' flowSOM.res <- ReadInput(fileName, compensate=TRUE,transform=TRUE,
#'                          scale=TRUE)
#' 
#' # Or read from flowFrame object
#' ff <- flowCore::read.FCS(fileName)
#' ff <- flowCore::compensate(ff,ff@@description$SPILL)
#' ff <- flowCore::transform(ff,
#'                  flowCore::transformList(colnames(ff@@description$SPILL),
#'                                  flowCore::logicleTransform()))
#' flowSOM.res <- ReadInput(ff,scale=TRUE)
#' 
#' # Build the self-organizing map and the minimal spanning tree
#' flowSOM.res <- BuildSOM(flowSOM.res,colsToUse=c(9,12,14:18))
#' flowSOM.res <- BuildMST(flowSOM.res)
#' 
#' # Apply metaclustering
#' metacl <- MetaClustering(flowSOM.res$map$codes,
#'                          "metaClustering_consensus",max=10)
#' 
#' # Get metaclustering per cell
#' flowSOM.clustering <- metacl[flowSOM.res$map$mapping[,1]]    
#'  
#' @export
ReadInput <- function(input, pattern=".fcs", 
                        compensate=FALSE, 
                        spillover=NULL,
                        transform=FALSE, 
                        toTransform=NULL, 
                        transformFunction = flowCore::logicleTransform(),
                        scale=FALSE, 
                        scaled.center=TRUE, 
                        scaled.scale=TRUE, 
                        silent=FALSE){
    
    fsom <- list(pattern=pattern, compensate=compensate, spillover=spillover,
                transform=transform, toTransform=toTransform,
                transformFunction = transformFunction, scale=scale)
    class(fsom) <- "FlowSOM"
    
    if(class(input) == "flowFrame"){
        fsom <- AddFlowFrame(fsom, input)        
    } else if(class(input) == "flowSet"){
        for(i in seq_along(input)){
            fsom <- AddFlowFrame(fsom, input[[i]])    
        }
    } else if(class(input) == "character"){    
        # Replace all directories by the files they contain
        toAdd <- NULL
        toRemove <- NULL
        for(i in seq_along(input)){
            if(file.info(input[i])$isdir){
                toAdd <- c(toAdd, list.files(input[i], pattern=pattern,
                                            full.names=TRUE))
                toRemove <- c(toRemove, i)
            }
        }
        if(!is.null(toRemove)){
            input <- c(input[-toRemove], toAdd)
        }
        
        # Select all files corresponding to the pattern
        input <- grep(pattern, input, value=TRUE)
        
        # Read all files
        if(length(input) > 0){
            for(i in seq_along(input)){
                if(file.exists(input[i])){
                    if(!silent) message("Reading file ", input[i],"\n")
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
        fsom$data <- scale(fsom$data, scaled.center, scaled.scale)
        fsom$scaled.center <- attr(fsom$data, "scaled:center")
        attr(fsom$data, "scaled:center") <- NULL
        fsom$scaled.scale <- attr(fsom$data, "scaled:scale") 
        attr(fsom$data, "scaled:scale") <- NULL
    }
    
    fsom
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
                fsom$spillover <- flowFrame@description$SPILL    
            } else if (!is.null(flowFrame@description$`$SPILLOVER`)){
                if(class(flowFrame@description$`$SPILLOVER`)=="matrix"){
                    fsom$spillover = flowFrame@description$`$SPILLOVER`
                    flowFrame@description$SPILL = fsom$spillover
                } else {
                    spilloverStr <- strsplit(
                        flowFrame@description$`$SPILLOVER`,
                        ",")[[1]]
                    n <- as.numeric(spilloverStr[1])
                    fsom$spillover <- t(matrix(as.numeric(spilloverStr[(n+2):
                                                length(spilloverStr)]),ncol=n))
                    colnames(fsom$spillover) <- spilloverStr[2:(n+1)]
                    flowFrame@description$SPILL <- fsom$spillover
                }
            } else {
                stop("No compensation matrix found")
            }
        }
        flowFrame <- flowCore::compensate(flowFrame, fsom$spillover)
    }
    
    # Transform
    if(fsom$transform){
        if(is.null(fsom$toTransform)){ 
            fsom$toTransform <- colnames(flowFrame@description$SPILL)
        } else{ 
            fsom$toTransform <- colnames(flowCore::exprs(flowFrame)[,
                                                            fsom$toTransform])
        }
        flowFrame <- flowCore::transform(flowFrame, 
            flowCore::transformList(fsom$toTransform,
                                    fsom$transformFunction))
    }
    
    # Save pretty names for nicer visualisation later on
    if(is.null(fsom$prettyColnames)){
        n <- flowFrame@parameters@data[, "name"]
        d <- flowFrame@parameters@data[, "desc"]
        d[is.na(d)] <- n[is.na(d)]
        if(any(grepl("#",d))){
            # Support for hashtag notation: 
            # antibody#fluorochrome -> antibody (fluorochrome)
            fsom$prettyColnames <- gsub("#(.*)$"," (\\1)",d)
        } else {
            fsom$prettyColnames <- paste(d, " <", n, ">", sep="")
        }
        names(fsom$prettyColnames) <- colnames(flowCore::exprs(flowFrame))
    }
  
    # Add the data to the matrix
    f <- flowCore::exprs(flowFrame)
    attr(f, "ranges") <- NULL
    name <- flowFrame@description$FIL
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
    
    fsom
}