devtools::document()
devtools::check()
devtools::build()


devtools::install_github("saeyslab/FlowSOM")

library(FlowSOM)
fileName <- system.file("extdata","lymphocytes.fcs",package="FlowSOM")
flowSOM.res <- FlowSOM(fileName, compensate=TRUE,transform=TRUE,
                       scale=TRUE,colsToUse=c(9,12,14:18),maxMeta=10)