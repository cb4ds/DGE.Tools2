require(testthat)
require(stats)

require(DGEobj)
require(DGEtools)
require(limma)
require(edgeR)

t_obj1 <- readRDS(system.file("exampleObj.RDS", package = "DGEobj", mustWork = TRUE))
