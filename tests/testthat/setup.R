require(testthat)
require(testthat)
require(stats)

require(DGEobj)
require(DGEtools)
require(limma)
require(edgeR)

DGEobj1 <- readRDS(system.file("exampleObj.RDS", package = "DGEobj", mustWork = TRUE))
