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

test_that("Having topGenes paramater be larger than total number of comparator rows throws warning", {
  testDf <- data.frame(firstColumn  = c("Malat1", "Plp1", "Ptgds"),
                         secondColumn = c(100, 200, 300))
  func <- expect_warning(findCorrelatedCells(comparator = testDf, visiumData = brain,
                                     assay="Spatial"))
  expect_equal(func$message, "The number of rows in comparator is less than inputted topGenes. Using 3 genes for correlation instead")
})
