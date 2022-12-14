#' Top 45 Oligodendrocyte markers obtained from Allen Brain Atlas
#'
#' @description
#' A dataset generated from extracting oligodendrocyte gene expression from the
#' reference scRNA-seq dataset of ~14,000 adult mouse cortical cell taxonomy
#' from the Allen Institute. The code used to generate this dataset is available
#' in topOligoGenes.R under the data-raw directory
#'
#' @format ## `topOligoGenes`
#' A data frame with 45 rows and 2 columns:
#' \describe{
#'   \item{gene}{Gene name}
#'   \item{exp}{Gene expression value}
#' }
#' @source The reference scRNA-seq dataset is from
#' <https://www.nature.com/articles/nn.4216>. The data is processed as an rds
#' object and can be downloaded from
#' <https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds>.
"topOligoGenes"
