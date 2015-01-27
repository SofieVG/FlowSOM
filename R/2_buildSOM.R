BuildSOM <- function(fsom, colsToUse=NULL, silent=FALSE, ...){
    # Build a SOM based on the data contained in the FlowSOM object
    #
    # Args:
    #     fsom:      FlowSOM object containing the data, as constructed
    #                by the ReadInput function
    #     colsToUse: column names or indices to use for building the SOM
    #     silent:    if TRUE, no progress updates will be printed
    #     ...:       options to pass on to the SOM function
    #
    # Returns:
    #     FlowSOM object containing the SOM
    
    if(!silent) message("Building SOM\n")
    
    if(is.null(colsToUse)){
        colsToUse <- 1:ncol(fsom$data)
    }

    fsom$map <- SOM(fsom$data[, colsToUse],silent=silent, ...)
    fsom$map$colsUsed <- colsToUse
    fsom$map$meanValues <-
        t(sapply(seq_len(fsom$map$xdim * fsom$map$ydim), function(i) {
            colMeans(subset(fsom$data, fsom$map$mapping[,1] == i))
        }))
    fsom$map$meanValues[is.nan(fsom$map$meanValues)] <- 0 
    fsom
}


SOM <- function (data, xdim=10, ydim=10, rlen=10, mst=1, alpha=c(0.05, 0.01),
                    radius = quantile(nhbrdist, 0.67) * c(1, 0), 
                    init=FALSE, distf=2, silent=FALSE, codes=NULL) 
{
    # Train a Self Organizing Map, code strongly based on the R kohonen package
    # However, next to rectangular grids, MST neighbourhoods can be computed in
    # between
    #
    # Args:
    #     data:   matrix containing the datapoints to train the SOM on
    #     xdim:   width of the grid
    #     ydim:   height of the grid
    #     rlen:   number of times the dataset is presented during training
    #     mst:    number of times a minimal spanning tree is build
    #             In total, the training set will be presented rlen*mst times
    #     alpha:  values to determine the start and end learning rate
    #     radius: values for the start and end neighbourhood radius
    #     init:   logical. Use smarter (but slower) initialization
    #     distf:  distance function to be used. 1=manhattan, 2=euclidean, 
    #             3=chebyshev
    #     silent: if T, no progress updates will be printed
    #
    # Returns:
    #     A list containing all parameter settings and results
    
    # Initialize the grid
    grid <- expand.grid(1:xdim,1:ydim)
    nCodes <- nrow(grid)
    if(is.null(codes)){
        if(init){
            starters <- Initialize(data, nCodes)
            message("Initialization ready\n")
        } else {
            starters <- sample(1:nrow(data), nCodes, replace = FALSE)        
        }
        codes <- data[starters, , drop = FALSE]
    }
    
    # Initialize the neighbourhood
    nhbrdist <- as.matrix(stats::dist(grid, method = "maximum"))
    
    # Initialize the radius
    if(mst==1){
        radius <- list(radius)
        alpha <- list(alpha)
    } else {
        radius <- seq(radius[1], radius[2], length.out = mst+1)
        radius <- lapply(1:mst, function(i){c(radius[i], radius[i+1])})
        alpha <- seq(alpha[1], alpha[2], length.out = mst+1)
        alpha <- lapply(1:mst, function(i){c(alpha[i], alpha[i+1])})
    }

    # Compute the SOM
    for(i in seq_len(mst)){
        res <- .C("C_SOM", data = as.double(data), 
                        codes = as.double(codes), 
                        nhbrdist = as.double(nhbrdist), 
                        alpha = as.double(alpha[[i]]), 
                        radius = as.double(radius[[i]]), 
                        xdists = double(nCodes), 
                        n = as.integer(nrow(data)), 
                        px = as.integer(ncol(data)), 
                        ncodes = as.integer(nCodes), 
                        rlen = as.integer(rlen), 
                        distf = as.integer(distf))

        codes <- matrix(res$codes, nrow(codes), ncol(codes))
        colnames(codes) <- colnames(data)
        nhbrdist <- Dist.MST(codes)
    }
    
    if(!silent) message("Mapping data to SOM\n")
    mapping <- MapDataToCodes(codes,data)
    
    list(xdim=xdim, ydim=ydim, rlen=rlen, mst=mst, alpha=alpha,
        radius=radius, init=init, distf=distf,
        grid=grid, codes=codes, mapping=mapping)
}


MapDataToCodes <- function (codes, newdata) 
{
    # Assign nearest node to each datapoint
    #
    # Args:
    #     codes: matrix with nodes of the SOM
    #     newdata: datapoints to assign
    #
    # Returns:
    #     Array with nearest node id for each datapoint
    
    nnCodes <- .C("C_mapDataToCodes", 
                    as.double(newdata[,colnames(codes)]), 
                    as.double(codes),
                    as.integer(nrow(codes)),
                    as.integer(nrow(newdata)),
                    as.integer(ncol(codes)),
                    nnCodes = integer(nrow(newdata)),
                    nnDists = double(nrow(newdata)))
    cbind(nnCodes$nnCodes, nnCodes$nnDists)
}

Initialize <- function(X, k){
    # Select k well spread points from X
    #
    # Args:
    #     X: matrix in which each row represents a point
    #     k: number of points to choose
    #
    # Returns:
    #     array containing indices of selected rows
    
    
    # Start with a random point
    selected <- numeric(k)
    selected[1] <- sample(1:nrow(X), 1)
    dists <- apply(X, 1, function(x)sum((x-X[selected[1], ])^2))
    
    for(i in seq_len(k-1)+1){ #2:k
        # Add point wich is furthest away from all previously selected points
        selected[i] <- which.max(dists)
        # Update distances
        dists <- pmin(dists, 
                    apply(X, 1, function(x)sum((x-X[selected[i], ])^2)))
    }
    
    selected
}

Purity <- function(realClusters, predictedClusters, weighted=TRUE){
    # Calculate mean weighted cluster purity
    #
    # realClusters: array with real cluster values
    # predictedClusters: array with predicted cluster values
    # weighted:     logical. Should the mean be weighted
    #                depending on the number of poins in the predicted clusters
    
    t <- table(predictedClusters, realClusters)
    maxPercentages <- apply(t/rowSums(t), 1, max)
    
    if(weighted)    weightedPercentages <- maxPercentages * rowSums(t)/sum(t)
    else weightedPercentages <- maxPercentages/nrow(t)
    
    c(sum(weightedPercentages), min(maxPercentages), 
            sum(maxPercentages<0.75))
}

Dist.MST <- function(X){
    # Calculate distance matrix using a minimal spanning tree neighbourhood
    #
    # Args:
    #     X: matrix in which each row represents a point
    #
    # Returns:
    #     Distance matrix
    adjacency <- dist(X, method = "euclidean")
    fullGraph <- graph.adjacency(as.matrix(adjacency), mode = "undirected", 
    weighted = TRUE)
    mst <- minimum.spanning.tree(fullGraph)
    shortest.paths(mst, v=V(mst), to=V(mst), weights=NA)
}