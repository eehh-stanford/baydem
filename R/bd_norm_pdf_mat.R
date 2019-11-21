#' @title Normalize the rows of an input density matrix
#'
#' @description The input density matrix fMat has dimensions S x G,
#' where S is the number of samples and G is the number of grid points
#' at which fMat was evaluated.
#' Normalize each row to integrate to 1 using the input grid spacing dy.
#'
#' @param fMat The input density matrix with dimensions S x G
#' @param dy The grid spacing
#'
#' @return The density matrix with normalized rows
#'
#' @export
bd_norm_pdf_mat <- function(fMat, dy) {
  S <- nrow(fMat)
  for (s in 1:S) {
    fMat[s, ] <- fMat[s, ] / sum(fMat[s, ]) / dy
  }
  return(fMat)
}
