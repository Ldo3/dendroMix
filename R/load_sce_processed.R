#' Load the processed SCE dataset used in dendroMix examples
#'
#' This function downloads and loads a large processed single-cell
#' dataset used in the vignette and examples. The file is hosted
#' externally so the dendroMix package remains lightweight.
#'
#' The dataset contains a PCA-reduced matrix (41,159 × 10) derived
#' from a larger single-cell experiment object.
#'
#' @return A matrix.
#' @examples
#' \dontrun{
#'   sce <- load_sce_processed()
#'   dim(sce)
#' }
#' @export
load_sce_processed <- function() {
     
     url <- "https://github.com/Ldo3/dendroMix/releases/download/v0.1.0/sce_processed.qs"
     
     # Temporary storage
     dest <- tempfile(fileext = ".qs")
     
     # Download the file
     utils::download.file(url, dest, mode = "wb")
     
     message("Downloaded processed SCE dataset to: ", dest)
     
     # Load the data
     qs::qread(dest)
}
