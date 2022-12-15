#' Computes the Pearson correlation between the gene expression of a
#' comparator dataframe and each cell spot in 10X Visium spatial data.
#'
#' A function that computes the pearson correlation of the gene expression of
#' the genes in a comparator dataframe and each spot in 10X
#' Visium spatial dataset. The alternative hypothesis is that there is a
#' positive correlation between the spot and comparator dataframe.
#'
#' @param comparator A dataframe of gene expression data with two columns. First
#' column contains gene names and second column contains corresponding gene
#' expression values. The dataframe can have no or any column names.
#' @param visiumData A 10X Visium spatial dataset loaded as a Seurat object
#' @param assay The assay in the Seurat object to use for correlation with
#' comparator, uses the Spatial assay by default
#' @param slot The specific assay slot to use, uses data slot (i.e. normalized
#' data matrix) by default.
#' @param topGenes The number of genes in the comparator to use for pearson
#' correlation with spatial data, by default uses all genes in comparator. By
#' entering an argument to topGenes, the function will select the topGenes
#' number of rows with highest gene expression values.
#'
#' @return The inputted visiumData Seurat object with additional metadata fields
#' "pearson.pvalue" and "pearson.estimate" for each Visium spot.
#' "pearson.pvalue" stores the correlation coefficient
#' "pearson.pvalue" stores the p-value of the correlation
#'
#' @examples
#' \dontrun{
#' # Using the 10X Genomics Visium Mouse Brain Dataset from SeuratData
#' SeuratData::InstallData("stxBrain")
#' brain <- SeuratData::LoadData("stxBrain", type = "posterior1")
#' # Using topOligoGenes dataset in package
#' brainWithPearson <- findCorrelatedCells(comparator = topOligoGenes, visiumData = brain)
#' # View metadata
#' brainWithPearson[[]]
#' }
#'
#' @export
#' @import Seurat tibble dplyr magrittr
#' @importFrom stats cor
#'
#' @references
#' Add references here
findCorrelatedCells <- function(comparator, visiumData, assay="Spatial",
                           slot="data", topGenes=nrow(comparator)) {
  # Check the dimension of comparator
  if (dim(comparator)[2] != 2) {
    stop("Comparator must contain 2 columns")
  }

  # Check data type of comparator columns
  if (typeof(comparator[[1]]) != "character") {
    stop("First column of comparator must contain gene names of type character")
  }
  if (!is.numeric(comparator[[2]])) {
    stop("Second column of comparator must contain gene expression values of type numeric")
  }

  # Process comparator data
  comparatorCopy <- comparator
  colnames(comparatorCopy) <- c("geneName", "comparatorExprValue")
  topComparatorGenes <- comparatorCopy %>%
    dplyr::arrange(desc(comparatorExprValue)) %>%
    dplyr::slice(1:topGenes)

  # Check if there is enough common genes in comparator and visiumData to carry
  # out pearson correlations
  if (length(intersect(rownames(x = visiumData), topComparatorGenes$geneName))
      < 3) {
    stop(sprintf("The %s selected genes in the comparator and do not have enough genes in common with visiumData for pearson correlations. Need at least 3 common genes",
                 topGenes))
  }

  # Process visium assay slot
  visiumAssaySlot <- Seurat::GetAssayData(object=visiumData, slot=slot, assay=assay)
  visiumAssaySlot <- as.data.frame(visiumAssaySlot)
  visiumAssaySlot <- tibble::rownames_to_column(visiumAssaySlot, "geneName")

  # Get correlations between comparator and each Visium spot
  pValues <- numeric(ncol(x = visiumData))
  estimate <- numeric(ncol(x = visiumData))
  for (i in 2:ncol(visiumAssaySlot)) {
    spotGeneExpr <- data.frame(geneName = visiumAssaySlot$geneName,
                               spotExprValue = visiumAssaySlot[[i]])
    comparatorAndSpotDF <- merge(topComparatorGenes, spotGeneExpr,
                                 by = "geneName")
    corrTestResult <- cor.test(comparatorAndSpotDF$comparatorExprValue,
                               comparatorAndSpotDF$spotExprValue,
                               alternative="greater")
    pValues[i-1] <- corrTestResult$p.value
    estimate[i-1] <- corrTestResult$estimate
  }

  # Add correlation information to Seurat metadata
  names(pValues) <- colnames(x = visiumData)
  names(estimate) <- colnames(x = visiumData)
  visiumData <- Seurat::AddMetaData(object = visiumData, metadata=pValues,
                                    col.name = 'pearson.pvalue')
  visiumData <- Seurat::AddMetaData(object = visiumData, metadata=estimate,
                                    col.name = 'pearson.estimate')
  return(visiumData)
}

#' Visualize the significance computed Pearson correlations in each visium cell
#' spot
#'
#' A function that visualizes the pearson correlation significance of each
#' Visium spot, overlayed on the Visium image
#'
#' @param visiumData A 10X Visium spatial dataset loaded as a Seurat object that
#' has been processed using the findCorrelatedCells function
#'
#' @return A plot showing the pearson correlation significance of each
#' Visium spot overlayed on the Visium image
#'
#' @examples
#' \dontrun{
#' # Using the 10X Genomics Visium Mouse Brain Dataset from SeuratData
#' SeuratData::InstallData("stxBrain")
#' brain <- SeuratData::LoadData("stxBrain", type = "posterior1")
#' # Using topOligoGenes dataset in package
#' brainWithPearson <- findCorrelatedCells(comparator = topOligoGenes,
#' visiumData = brain)
#' visualizeCells(brainWithPearson)
#' }
#'
#' @export
#' @import Seurat
#' @importFrom stats cor
#'
#' @references
#' Add references here
visualizeCells <- function(visiumData) {
  # Check that visiumData has metadata column pearson.pvalue
  `%!in%` <- Negate(`%in%`)
  if ("pearson.pvalue" %!in% colnames(visiumData[[]])) {
    stop("visiumData does not include pearson.pvalue metadata column. Ensure visiumData has been run through the findCorrelatedCells function before using the visualizeCells function")
  }
  visiumDataCopy <- visiumData
  significant <- (visiumDataCopy[[]]$pearson.pvalue <= 0.05)
  visiumDataCopy <- Seurat::AddMetaData(object = visiumDataCopy, metadata=significant, col.name = 'Significant')
  Seurat::SpatialPlot(object = visiumDataCopy, group.by = "Significant",  alpha = c(0.8, 1))
}
