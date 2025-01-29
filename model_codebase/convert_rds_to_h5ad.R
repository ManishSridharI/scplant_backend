#!/usr/bin/Rscript --vanilla
rm(list = ls())

if (!dir.exists("packrat")){
  # Initializes the current working directory as a Packrat project. 
  packrat::init(
    infer.dependencies = FALSE,
    enter = TRUE,
    restart = FALSE
  )

  # Install packages under the project.
  install.packages("BiocManager", repos = "https://cloud.r-project.org/");
  install.packages("devtools", repos = "https://cloud.r-project.org/");
  install.packages("Matrix", repos = "https://cloud.r-project.org/");
  install.packages('dplyr', repos = "https://cloud.r-project.org/");
  install.packages("tidyr", repos = "https://cloud.r-project.org/");
  install.packages("tibble", repos = "https://cloud.r-project.org/");
  install.packages("stringr", repos = "https://cloud.r-project.org/");
  install.packages("ggplot2", repos = "https://cloud.r-project.org/");
  install.packages("argparse", repos = "https://cloud.r-project.org/");
  install.packages("Seurat", repos = "https://cloud.r-project.org/");
  install.packages("reticulate", repos = "https://cloud.r-project.org/");
  install.packages("anndata", repos = "https://cloud.r-project.org/");
}


rm(list = ls())

set.seed(1)

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(argparse)

library(Matrix)
library(Seurat)
library(reticulate)
library(anndata)

sc <- reticulate::import("scanpy")
scipy <- reticulate::import("scipy")


##################################################
# Constants/Variables
##################################################


##################################################
# Argparse
##################################################

parser <- argparse::ArgumentParser()

parser$add_argument("-i", "--input", type = "character", help = "Input rds file")
parser$add_argument("-o", "--output", type = "character", help = "Output h5ad file")

args <- parser$parse_args()


input_file_path <- tryCatch({
  file.path(args$input)
}, error = function(e) {
  return(NULL)
})

output_file_path <- tryCatch({
  file.path(args$output)
}, error = function(e) {
  return(NULL)
})


if (is.null(input_file_path)) {
  print("Input file path is invalid")
  quit()
}

if (is.null(output_file_path)) {
  print("Output file path is invalid")
  quit()
}


##################################################
# Output folder
##################################################

output_path <- dirname(output_file_path)

if (!dir.exists(output_path)) {
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  if (!dir.exists(output_path)) {
    quit(status = 1)
  }
}


##################################################
# Input data
##################################################

dat <- readRDS(
  file = input_file_path
)


##################################################
# Process data
##################################################

# Join layers in RNA assay
dat[["RNA"]] <- tryCatch({
  JoinLayers(dat[["RNA"]])
}, error = function(e) {
  return(dat[["RNA"]])
})

matr <- GetAssayData(object = dat,  layer = "counts")

matr <- Matrix(matr, sparse = TRUE)

Matrix::writeMM(matr, file = file.path(output_path, "counts.mtx"))

x <- r_to_py(scipy$io$mmread(file.path(output_path, "counts.mtx")))

x <- x$transpose()$tocsr()

adata_seurat <- sc$AnnData(
    X = x,
    obs = dat[[]],
    var = GetAssay(dat)[[]]
)

adata_seurat$obsm[["X_umap"]] <- Embeddings(dat, "umap")
adata_seurat$obsm[["X_tsne"]] <- Embeddings(dat, "tsne")
adata_seurat$obsm[["X_pca"]] <- Embeddings(dat, "pca")


##################################################
# Write data
##################################################

adata_seurat$write_h5ad(output_file_path)

