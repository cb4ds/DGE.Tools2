context("DGEtools - tests for checkSex.R functions")


test_that("checkSex.R: checkSex()", {
    skip_if_not("Chromosome" %in% names(DGEobj::getItem(t_obj1, "geneData_orig")))

    # testing for species rat
    sex_determination_plot <- checkSex(dgeObj  = t_obj1,
                                       species = "RAT")
    expect_s3_class(sex_determination_plot, c("gg", "ggplot"))

    # testing for species otherthen rat (mouse and human)
    sex_determination_plot <- checkSex(dgeObj  = t_obj1,
                                       species = "human")
    expect_s3_class(sex_determination_plot, c("gg", "ggplot"))

    # testing for filtered count matrix
    sex_determination_plot <- checkSex(dgeObj  = t_obj1,
                                       species = "human",
                                       orig    = FALSE)
    expect_s3_class(sex_determination_plot, c("gg", "ggplot"))

    # testing for labelCol and sexCol
    ##  with labelCol
    sex_determination_plot <- checkSex(dgeObj     = t_obj1,
                                       species    = "rat",
                                       labelCol   = "SampleName",
                                       showLabels = TRUE)
    expect_s3_class(sex_determination_plot, c("gg", "ggplot"))

    ##  with sexCol
    t_obj1$design_orig$Sex <- c("Male", "Female")
    sex_determination_plot <- checkSex(dgeObj     = t_obj1,
                                       species    = "rat",
                                       sexCol     = "Sex",
                                       showLabels = TRUE)
    expect_s3_class(sex_determination_plot, c("gg", "ggplot"))

    ##  with both sexCol and labelCol
    sex_determination_plot <- checkSex(dgeObj     = t_obj1,
                                       species    = "rat",
                                       sexCol     = "Sex",
                                       labelCol   = "SampleName",
                                       showLabels = TRUE)
    expect_s3_class(sex_determination_plot, c("gg", "ggplot"))
})
