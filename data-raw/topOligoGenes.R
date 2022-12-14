library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tibble)
library(magrittr)

# The allen_cortex.rds is not included in this package as it is over
# 1 GB in size. It can be downloaded from
# <https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds>.
# To run this code, download the file and place it in the same directory as this
# script
allen_reference <- readRDS("allen_cortex.rds")
allen_reference <- SCTransform(allen_reference, ncells = 3000,
                               verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
oligo <- subset(x = allen_reference, subset = subclass == "Oligo")
averageExpressionOligo <- AverageExpression(oligo, group.by = "subclass")
avgExpressionOligoSCT <- as.data.frame(averageExpressionOligo[["SCT"]])
avgExpressionOligoSCT <- tibble::rownames_to_column(avgExpressionOligoSCT,
                                                    "gene")
colnames(avgExpressionOligoSCT) <- c("gene", "exp")
topOligoGenes <- avgExpressionOligoSCT %>%
  dplyr::arrange(desc(exp)) %>%
  dplyr::slice(1:45)

# usethis::use_data(topOligoGenes, overwrite = TRUE)
