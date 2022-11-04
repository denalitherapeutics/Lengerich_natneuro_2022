#' Read in 10x output files from an AWS S3 bucket and create a Seurat Object
#'
#' @param bucket Scalar character, the AWS S3 bucket containing the Cell Ranger
#' outputs (matrix.mtx.gz, features.tsv.gz, and barcodes.tsv.gz)
#' @param region Sclalar character, the region for the S3 bucket (usually
#' "us-east-1" for raw data files)
#' @param project Scalar character or `NULL`.
#' @param min.cells Scalar count, the minimum number of cells a feature must be
#' detected in to be retained.
#' @param min.features Scalar count, the minimum number of features detected in
#' a cell to be retained.
#' @import aws.s3
#' @import Seurat
#' @import dplyr
#' @import Matrix
#' @importFrom checkmate assert_count, assert_string
#' @return Returns a Seurat object
#' @export
s3.read10X <- function(bucket, region, project = NULL, min.cells = 0,
                       min.features = 0){
  checkmate::assert_string(bucket)
  checkmate::assert_string(region)
  checkmate::assert_string(project, null.ok = TRUE)
  checkmate::assert_count(min.cells)
  checkmate::assert_count(min.features)
  matrix <- aws.s3::s3read_using(FUN = readMM,
                                 bucket = bucket,
                                 object = "matrix.mtx.gz",
                                 opts = list(region = region))


  feature.names <- aws.s3::s3read_using(FUN = read.delim,
                                        bucket = bucket,
                                        object = "features.tsv.gz",
                                        header = FALSE,
                                        stringsAsFactors = FALSE,
                                        opts = list(region = region))
  barcode.names <- aws.s3::s3read_using(FUN = read.delim,
                                        bucket = bucket,
                                        object = "barcodes.tsv.gz",
                                        header = FALSE,
                                        stringsAsFactors = FALSE,
                                        opts = list(region = region))
  colnames(matrix) <- barcode.names$V1
  rownames(matrix) <- feature.names$V2

  if (is.null(project)) {
    object <- Seurat::CreateSeuratObject(counts = matrix,
                                         min.cells = min.cells,
                                         min.features = min.features)
  } else {
    object <- Seurat::CreateSeuratObject(counts = matrix,
                                         project = project,
                                         min.cells = min.cells,
                                         min.features = min.features)
  }
}
