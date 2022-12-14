#' Computes the Pearson correlation between the gene expression of a
#' comparator dataframe and each spot in 10X Visium spatial data. Visualize
#' the pearson correlation significance of each Visium spot.
#'
#'
#' A function that computes the pearson correlation of the gene expression of
#' the top expressed genes in a comparator dataframe and each spot in 10X
#' Visium spatial dataset. The alternative hypothesis is that there is a
#' positive correlation between the spot and comparator dataframe.
#'
#' @param comparator A dataframe of gene expression data, with one column
#' containing expression values. The rownames of the dataframe has to be set to
#' the gene names corresponding to the expression values
#' @param visiumData A 10X Visium spatial dataset loaded as a Seurat object
#' @param assay The assay in the Seurat object to use for correlation with
#' comparator, uses the Spatial assay by default
#' @param slot The specific assay data to use, uses normalized data matrix
#' (i.e. data slot) by default.
#' @param topGenes The number of top genes in the comparator to use for pearson
#' correlation with spatial data, by default 25. If the comparator has less than
#' 25 genes, by default all genes in comparator are used for correlation with
#' spatial data.
#'
#' @return The inputted spatialDataset Seurat object modified metadata fields
#' "pearson.pvalue" and "pearson.estimate". "pearson.pvalue" stores the coefficient
#' of the pearson correlation between the comparator and the corresponding
#' Visium spot. "pearson.pvalue" that corresponds to the p-value of the pearson
#' correlation between the gene expression of the comparator and each Visium
#' spot
#'
#'
#' @export
#' @import Seurat tibble dplyr magrittr
#' @importFrom stats cor
#'
#' @references
#' Add references here
findCorrelatedCells <- function(comparator, visiumData, assay="Spatial",
                           slot="data", topGenes=25) {
  # Check the dimension of comparator
  if (dim(comparator)[2] != 2) {
    stop("Comparator must contain 2 columns")
  }

  # Check data type of comparator columns
  if (typeof(comparator[1]) != "character") {
    stop("First column of comparator must contain gene names of type character")
  }
  if (typeof(comparator[2]) != "numeric") {
    stop("Second column of comparator must contain gene names of type numeric")
  }

  # Check if topGenes is less than or equal to total comparator rows
  # If not, throw a warning and adjust topGenes to the total number of
  # comparator rows
  if (nrow(comparator) < topGenes) {
    warning(sprintf("The number of rows in comparator is less than inputted
            topGenes. Using %s genes for correlation instead", nrow(comparator)))
    topGenes <- nrow(comparator)
  }

  # Process comparator data
  comparatorCopy <- comparator
  colnames(comparatorCopy) <- c("geneName", "comparatorExprValue")
  topComparatorGenes <- comparatorCopy %>%
    dplyr::arrange(desc(comparatorExprValue)) %>%
    dplyr::slice(1:topGenes)

  # Process visium assay slot
  visiumAssaySlot <- Seurat::GetAssayData(object=visiumData, slot=slot, assay=assay)
  visiumAssaySlot <- as.data.frame(visiumAssaySlot)
  visiumAssaySlot <- tibble::rownames_to_column(visiumAssaySlot, "geneName")

  # Get correlations
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

#' A function that visualizes the pearson correlation significance of each
#' Visium spot overlayed on the Visium image
#'
#' @param visiumData A 10X Visium spatial dataset loaded as a Seurat object that
#' has been processed using the findCorrelatedCells function
#'
#' @return A plot showing the pearson correlation significance of each
#' Visium spot overlayed on the Visium image
#'
#' examples will be added
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
  if ("pearson.pvalue" %!in% visiumData[[]]) {
    stop("visiumData does not include pearson.pvalue metadata column. Ensure
         visiumData has been run through the findCorrelatedCells function before
         using the visualizeCells function")
  }
  visiumDataCopy <- visiumData
  significant <- (visiumDataCopy[[]]$pearson.pvalue <= 0.05)
  visiumDataCopy <- Seurat::AddMetaData(object = visiumDataCopy, metadata=significant, col.name = 'Significant')
  SpatialPlot(object = visiumDataCopy, group.by = "Significant",  alpha = c(0.8, 1))
}
