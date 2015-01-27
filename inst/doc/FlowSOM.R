### R code from vignette source 'FlowSOM.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: FlowSOM.Rnw:80-91
###################################################
set.seed(42)
library(flowCore)
library(FlowSOM)

fileName <- system.file("extdata","lymphocytes.fcs",
                        package="FlowSOM")
fSOM <- ReadInput(fileName,compensate = TRUE,transform = TRUE, 
                    toTransform=c(8:18),scale = TRUE)

ff <- read.FCS(fileName)
fSOM <- ReadInput(ff,compensate = TRUE,transform = TRUE, scale = TRUE)


###################################################
### code chunk number 3: FlowSOM.Rnw:101-102
###################################################
str(fSOM)


###################################################
### code chunk number 4: FlowSOM.Rnw:122-124
###################################################
fSOM <- BuildSOM(fSOM,colsToUse = c(9,12,14:18))
str(fSOM$map)


###################################################
### code chunk number 5: FlowSOM.Rnw:136-138
###################################################
fSOM <- BuildMST(fSOM)
str(fSOM$MST)


###################################################
### code chunk number 6: FlowSOM.Rnw:144-145
###################################################
PlotStars(fSOM)


###################################################
### code chunk number 7: FlowSOM.Rnw:150-153
###################################################
fSOM <- UpdateNodeSize(fSOM, reset=TRUE)
PlotStars(fSOM,MST=FALSE)
fSOM <- UpdateNodeSize(fSOM)


###################################################
### code chunk number 8: FlowSOM.Rnw:157-186
###################################################
library(flowUtils)
flowEnv <- new.env()
ff_c <- compensate(ff,ff@description$SPILL)
colnames(ff_c)[8:18] <- paste("Comp-",colnames(ff_c)[8:18],sep="")
gatingFile <- system.file("extdata","manualGating.xml", 
                        package="FlowSOM")
read.gatingML(gatingFile, flowEnv) 
filterList <- list( "B cells" = flowEnv$ID52300206,
                    "ab T cells" = flowEnv$ID785879196,
                    "yd T cells" = flowEnv$ID188379411,
                    "NK cells" = flowEnv$ID1229333490,
                    "NKT cells" = flowEnv$ID275096433
                )

results <- list()
for(cellType in names(filterList)){
    results[[cellType]] <- filter(ff_c,filterList[[cellType]])@subSet
}

manual <- rep("Unknown",nrow(ff))
for(celltype in names(results)){
    manual[results[[celltype]]] <- celltype
}
# Use a factor to define order of the cell types
manual <- factor(manual,levels = c("Unknown","B cells",
                                    "ab T cells","yd T cells", 
                                    "NK cells","NKT cells"))

PlotPies(fSOM,cellTypes=manual)


###################################################
### code chunk number 9: FlowSOM.Rnw:200-202
###################################################
metaClustering <- metaClustering_consensus(fSOM$map$codes,k=7)
PlotPies(fSOM,cellTypes=manual,clusters = metaClustering)


###################################################
### code chunk number 10: FlowSOM.Rnw:206-207
###################################################
metaClustering_perCell <- metaClustering[fSOM$map$mapping[,1]]


