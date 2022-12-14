brain <- SeuratData::LoadData("stxBrain", type = "posterior1")

test_that("Throws an error when there is no pearson.pvalue metadata column", {
  expect_error(visualizeCells(brain))
  expect_equal(geterrmessage(), "visiumData does not include pearson.pvalue metadata column. Ensure visiumData has been run through the findCorrelatedCells function before using the visualizeCells function")
})
