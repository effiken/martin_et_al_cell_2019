#  Copied from scrattch package by Lucas Graybuck
#
#  https://github.com/AllenInstitute/scrattch.io/blob/master/R/write_10x.R
#
#' Write a dgCMatrix to an h5 file similar to cellRanger format
#'
#' @param mat a dgCMatrix to write.
#' @param cols_are Whether columns are "gene_names" or "sample_names". If "gene_names", mat will be transposed to match 10X conventions.
#' @param h5_target The target .h5 file for output.
#' @param ref_name Reference name for storing the data
#' @param gene_ids If available, ENSEMBL IDs for each gene
#'
#' @return an .h5 file with these mappings from dgCMatrix -> .h5:
#' \itemize{
#'   \item colnames(mat) -> /ref_name/barcodes (after transposition if cols_are == "gene_names")
#'   \item rownames(mat) -> /ref_name/gene_names (after transposition if cols_are == "gene_names")
#'   \item mat@x -> /ref_name/data
#'   \item mat@i -> /ref_name/indices
#'   \item mat@p -> /ref_name/indptr
#'   \item gene_ids -> /ref_name/gene
#' }
#'
write_dgCMatrix_h5 <- function(mat,
                               cols_are = "gene_names",
                               h5_target,
                               ref_name = "mm10-1.2.0_premrna",
                               gene_ids = NULL) {
  
  #library(Matrix)
  
  if(grepl("gene",cols_are)) {
    mat <- Matrix::t(mat)
  }
  
  # Create target file
  rhdf5::h5createFile(h5_target)
  # Create data group
  rhdf5::h5createGroup(h5_target,
                       ref_name)
  
  # Store sample ids (barcodes) and gene names
  rhdf5::h5write(colnames(mat),
                 h5_target,
                 paste0("/",ref_name,"/barcodes"))
  rhdf5::h5write(rownames(mat),
                 h5_target,
                 paste0("/",ref_name,"/gene_names"))
  
  if(is.null(gene_ids)) {
    gene_ids <- rownames(mat)
  }
  
  rhdf5::h5write(gene_ids,
                 h5_target,
                 paste0("/",ref_name,"/gene"))
  
  # Store dimensions as shape
  rhdf5::h5write(dim(mat),
                 h5_target,
                 paste0("/",ref_name,"/shape"))
  
  # Store values from mat@x as data
  rhdf5::h5createDataset(h5_target,
                         paste0("/",ref_name,"/data"),
                         dims = length(mat@x),
                         storage.mode = "integer",
                         chunk = 1000,
                         level = 4)
  rhdf5::h5write(mat@x,
                 h5_target,
                 paste0("/",ref_name,"/data"))
  
  # Store row indices from mat@i as indices
  rhdf5::h5createDataset(h5_target,
                         paste0("/",ref_name,"/indices"),
                         dims = length(mat@i),
                         storage.mode = "integer",
                         chunk = 1000,
                         level = 4)
  rhdf5::h5write(mat@i,
                 h5_target,
                 paste0("/",ref_name,"/indices"))
  
  # Store column pointers from mat@p as indptr
  rhdf5::h5write(mat@p,
                 h5_target,
                 paste0("/",ref_name,"/indptr"))
  
}








#' Read a whole sparse matrix directly from a 10X-style .h5 file
#'
#' @param h5 .h5 file to read.
#' @param target The sparse matrix to read within the .h5 file.
#'
read_10x_dgCMatrix <- function(h5,
                               target) {
  #library(Matrix)
  
  root <- rhdf5::H5Fopen(h5)
  
  i_path <- paste0(target,"/indices")
  p_path <- paste0(target,"/indptr")
  x_path <- paste0(target,"/data")
  dims_path <- paste0(target,"/shape")
  
  print("Reading indices")
  i <- read_tome_vector(root, i_path)
  print("Reading pointers")
  p <- read_tome_vector(root, p_path)
  print("Reading values")
  x <- read_tome_vector(root, x_path)
  print("Reading dimensions")
  dims <- read_tome_vector(root, dims_path)
  
  H5Fclose(root)
  
  print("Assembling dgCMatrix")
  Matrix::sparseMatrix(i = i,
                       p = p,
                       x = x,
                       index1 = FALSE,
                       dims = dims)
  
}


