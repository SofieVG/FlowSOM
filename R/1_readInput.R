ReadInput <- function(input, pattern=".fcs", compensate=FALSE, spillover=NULL,
                        transform=FALSE, toTransform=NULL, scale=FALSE, 
                        scaled.center=TRUE, scaled.scale=TRUE, silent=FALSE){
    # Take all possible inputs and return a matrix with preprocessed data
    # (compensated, transformed, scaled)
    #
    # Args:
    #     input:      a flowFrame, a flowSet or an array of paths to files 
    #                 or directories
    #     pattern:    if input is an array of file- or directorynames, select
    #                 only files containing pattern
    #     compensate: logical, does the data need to be compensated
    #     spillover:  spillover matrix to compensate with. 
    #                 If NULL and compensate=TRUE, we will look for $SPILL 
    #                 description in fcs file.
    #     transform:  logical, does the data need to be transformed with 
    #                 a logicle transform
    #     toTransform:column names or indices that need to be transformed.
    #                 If NULL and transform=TRUE, column names of $SPILL 
    #                 description in fcs file will be used.
    #     scale:      logical, does the data needs to be rescaled
    #     scaled.center: see scale
    #     scaled.scale:  see scale
    #     silent:     if T, no progress updates will be printed
    #
    # Returns:
    #     FlowSOM object containing the preprocessed data
    
    fsom <- list(pattern=pattern, compensate=compensate, spillover=spillover,
            transform=transform, toTransform=toTransform, scale=scale)
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
                        flowFrame <- flowFrame(read.table(input[i]))
                    } else { #else if(tools::file_ext(input[i]) == "fcs"){
                        flowFrame <- suppressWarnings(read.FCS(input[i]))
                    }
                    fsom <- AddFlowFrame(fsom, flowFrame)
                }
            }
        } else {
                stop("No files containing the pattern are found.")
        }
    } else {
        stop(paste("Inputs of type", class(input), "are not supported. 
        Please supply either a FlowFrame, a FlowSet or an array of valid 
        paths to files or directories."))
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

AddFlowFrame <- function(fsom, flowFrame){
    # Add a flowFrame to the data variable of the FlowSOM object
    #
    # Args:
    #     fsom:      FlowSOM object, as constructed by the ReadInput function
    #     flowFrame: flowFrame to add to the FlowSOM object
    #
    # Returns:
    #     FlowSOM object with data added
    
    # Compensation
    if(fsom$compensate){
        if(is.null(fsom$spillover)){
            if(!is.null(flowFrame@description$SPILL)){
                fsom$spillover <- flowFrame@description$SPILL    
            } else if (!is.null(flowFrame@description$`$SPILLOVER`)){
                spilloverStr <- strsplit(flowFrame@description$`$SPILLOVER`,
                                        ",")[[1]]
                n <- as.numeric(spilloverStr[1])
                fsom$spillover <- t(matrix(as.numeric(spilloverStr[(n+2):
                                            length(spilloverStr)]),ncol=n))
                colnames(fsom$spillover) <- spilloverStr[2:(n+1)]
                flowFrame@description$SPILL <- fsom$spillover
            } else {
                stop("No compensation matrix found")
            }
        }
        flowFrame <- compensate(flowFrame, fsom$spillover)
    }
    
    # Logicle transform
    if(fsom$transform){
        if(is.null(fsom$toTransform)){ 
            fsom$toTransform <- colnames(flowFrame@description$SPILL)
        } else{ 
            fsom$toTransform <- colnames(exprs(flowFrame)[,fsom$toTransform])
        }
        flowFrame <- transform(flowFrame, transformList(fsom$toTransform,
                                logicleTransform()))
    }
    
    # Save pretty names for nicer visualisation later on
    n <- flowFrame@parameters@data[, "name"]
    d <- flowFrame@parameters@data[, "desc"]
    d[is.na(d)] <- n[is.na(d)]
    fsom$prettyColnames <- paste(d, " <", n, ">", sep="")
    names(fsom$prettyColnames) <- colnames(exprs(flowFrame))
    
    # Add the data to the matrix
    f <- exprs(flowFrame)
    attr(f, "ranges") <- NULL
    if(is.null(fsom$data)){ 
        fsom$data <- f 
        fsom$metaData <- list()
        fsom$metaData[[flowFrame@description$FIL]] <- c(1, nrow(fsom$data))
    } else {
        fsom$data <- rbind(fsom$data, f)
        fsom$metaData[[flowFrame@description$FIL]] <- c(nrow(fsom$data) - 
                                                nrow(f) + 1, nrow(fsom$data))
    }
    fsom
}