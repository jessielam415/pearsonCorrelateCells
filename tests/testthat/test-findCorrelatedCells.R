brain <- SeuratData::LoadData("stxBrain", type = "posterior1")

test_that("Wrong number of comparator columns throws error", {
  testDf <- data.frame(firstColumn  = c("Malat1", "Plp1", "Ptgds"),
                    secondColumn = c(100, 200, 300),
                    thirdColumn = c("Malat11Info", "Plp1Info",
                                    "Plp1Info")
  )
  expect_error(findCorrelatedCells(comparator = testDf, visiumData = brain,
                                          assay="Spatial"))
  expect_equal(geterrmessage(), "Comparator must contain 2 columns")
})

test_that("Wrong comparator column type raises throws error", {
  testDf <- data.frame (firstColumn  = c(1, 2, 3),
                        secondColumn = c(100, 200, 300))
  expect_error(findCorrelatedCells(comparator = testDf, visiumData = brain,
                                   assay="Spatial"))
  expect_equal(geterrmessage(),
               "First column of comparator must contain gene names of type character")
  testDf1 <- data.frame (firstColumn  = c("Malat1", "Plp1", "Ptgds"),
                        secondColumn = c("100", "200", "300"))
  expect_error(findCorrelatedCells(comparator = testDf1, visiumData = brain,
                                   assay="Spatial"))
  expect_equal(geterrmessage(),
               "Second column of comparator must contain gene expression values of type numeric")
})

test_that("Less than 2 genes in common between selected comparator genes and visiumData throws error", {
  testDf <- data.frame (firstColumn  = c("Dbi", "Ptgds", "Gene not in visiumData"),
                        secondColumn = c(100, 200, 300))
  expect_error(findCorrelatedCells(comparator = testDf, visiumData = brain,
                                   assay="Spatial"))
  expect_equal(geterrmessage(),
               "The 3 selected genes in the comparator and do not have enough genes in common with visiumData for pearson correlations. Need at least 3 common genes")
})
