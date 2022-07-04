#' Retrieve and format reads
#'
#' @description
#' `loadReads()` applies some formatting on BAM files before further dependency
#' analysis - it offers to remove:
#'
#' * reads longer than a provided length
#' * reads not mapping any maturation sites
#' * reads mapping outside of a given range
#' * reads mapping multiple locations
#'
#' @param path `String` containing the path to the BAM files.
#' @param max_length `Integer` indicating the maximal authorized length for a
#' read. Reads longer than this value will be removed. Default to `Inf`.
#'
#' @return A (list of) `GRanges` object with all the reads.
#'
#' @examples
#'
#' loadReads(1, 1)
#' loadReads(10, 1)
#'
loadReads <- function(path, max_length) {
  x + y
}
