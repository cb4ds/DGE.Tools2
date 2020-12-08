context("DGEtools - tests for ggplotMDS.R functions")


test_that("ggplotMDS.R: ggplotMDS()", {
    skip_if(is.null(t_obj1$DGEList))
    skip_if(is.null(t_obj1$design$Barcode))

    mds_plot <- ggplotMDS(DGEdata = t_obj1,
                          colorBy = t_obj1$design$Barcode)
    expect_length(mds_plot, 2)
    expect_named(mds_plot, c("plot", "mdsobj"))
    expect_type(mds_plot, "list")
    expect_s3_class(mds_plot$plot, c("gg", "ggplot"))

    # testing parameter size and shape both (also colorName, shapeName, sizeName)
    mds_plot <- ggplotMDS(DGEdata   = t_obj1,
                          colorBy   = t_obj1$design$Barcode,
                          sizeBy    = t_obj1$design$Dose..mpk.,
                          shapeBy   = t_obj1$design$Species,
                          colorName = "color title",
                          shapeName = "shape title",
                          sizeName  = "size title")
    expect_s3_class(mds_plot$plot, c("gg", "ggplot"))
    # testing parameter size.
    mds_plot <- ggplotMDS(DGEdata = t_obj1,
                          colorBy = t_obj1$design$Barcode,
                          sizeBy  = t_obj1$design$Dose..mpk.)
    expect_s3_class(mds_plot$plot, c("gg", "ggplot"))
    # testing parameter shape.
    mds_plot <- ggplotMDS(DGEdata = t_obj1,
                          colorBy = t_obj1$design$Barcode,
                          shapeBy = t_obj1$design$Species)
    expect_s3_class(mds_plot$plot, c("gg", "ggplot"))

    # testing for discrete color values
    mds_plot <- ggplotMDS(DGEdata = t_obj1,
                          colorBy = t_obj1$design$Species,
                          colors  = c("green"))
    expect_s3_class(mds_plot$plot, c("gg", "ggplot"))

    # testing intercept lines
    mds_plot <- ggplotMDS(DGEdata = t_obj1,
                          colorBy = t_obj1$design$Species,
                          hlineIntercept = 0.25,
                          vlineIntercept = 0.25)
    expect_s3_class(mds_plot$plot, c("gg", "ggplot"))

    # testing assert statements
    expect_error(ggplotMDS(DGEdata = 1:10),
                 regexp = "DGEdata must be of class 'DGEList', 'DGEobj', or 'matrix'.")
    expect_error(ggplotMDS(DGEdata = t_obj1),
                 regexp = "colorBy must be specified and should be the length of the number of columns in DGEdata.")
    expect_error(ggplotMDS(DGEdata = t_obj1,
                           colorBy = t_obj1$design$Barcode,
                           shapeBy = 1:10),
                 regexp = "shapeBy should be the length of the number of columns in DGEdata.")
    expect_error(ggplotMDS(DGEdata = t_obj1,
                           colorBy = t_obj1$design$Barcode,
                           sizeBy  = 1:10),
                 regexp = "sizeBy should be the length of the number of columns in DGEdata.")
})

test_that("ggplotMDS.R: MDS_var_explained()", {
    skip_if(is.null(t_obj1$DGEList))
    skip_if(is.null(t_obj1$design$Barcode))

    mds_plot <- ggplotMDS(DGEdata = t_obj1,
                          colorBy = t_obj1$design$Barcode)
    var_result <- MDS_var_explained(mds_plot$mdsobj)
    expect_length(var_result, 3)
    expect_named(var_result, c("varexp", "cumvar", "var_explained"))
    expect_type(var_result, "list")
    expect_s3_class(var_result$varexp, c("gg", "ggplot"))
    expect_s3_class(var_result$cumvar, c("gg", "ggplot"))

    # testing mds matrix
    var_result <- MDS_var_explained(t_obj1$counts)
    expect_length(var_result, 3)
    expect_named(var_result, c("varexp", "cumvar", "var_explained"))
    expect_type(var_result, "list")
    expect_s3_class(var_result$varexp, c("gg", "ggplot"))
    expect_s3_class(var_result$cumvar, c("gg", "ggplot"))

    # testing assert statements
    expect_error(MDS_var_explained(), regexp = "mds is required and must be specified.")
})
