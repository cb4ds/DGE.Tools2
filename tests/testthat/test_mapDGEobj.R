context("DGEtools tests for mapDGEobj.R functions")


test_that('mapDGEobj.R: mapDGEobj()', {
    map_DGEobj <- mapDGEobj(t_obj1)
    expect_s3_class(map_DGEobj, "igraph")

    expect_error(mapDGEobj(obj),
                 regexp = "object 'obj' not found")
    expect_error(mapDGEobj("XYZ"),
                 regexp = "dgeObj must be of class 'DGEobj'.")
})
