#' Gene expression for 260 genes in 887 breast samples
#'
#' A data matrix containing a subset of the TCGA breast cancer dataset, with the gene level expression estimates 
#' in log2 transcripts per million for 887 breast samples.
#'
#' @format A data matrix with 260 rows and 887 columns.
#' @source \url{https://cancergenome.nih.gov/}
"breastTCGA_Mat"

#' Normal or Tumor status of breast samples
#'
#' A factor indicating whether 887 breast samples in breastTCGA_Mat are tumor or matched normal.
#'
#' @format A Factor of length 887 of levels NORMAL and TUMOR.
#' @source \url{https://cancergenome.nih.gov/}
"breastTCGA_Group"

#' ER positive or negative status of breast tumor samples
#'
#' A factor indicating whether 887 breast samples in breastTCGA_Mat are ER positive or ER negative.
#' The matched normals have empty values.
#'
#' @format A Factor of length 887 of levels Negative and Positive (with 111 missing values for the normals).
#' @source \url{https://cancergenome.nih.gov/}
"breastTCGA_ER"

#' Cancer Hallmark gene sets from the MSigDB collection
#'
#' A subset of the cancer hallmarks functional gene sets from the MSigDB collection.
#'
#' @format A list of length 10, with the hallmark gene set name, each a character vector of gene symbols.
#' @source \url{https://http://software.broadinstitute.org/gsea/msigdb/}
"msigdb_Hallmarks"
