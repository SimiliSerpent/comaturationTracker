#' Retrieve and format reads
#'
#' @description
#' `loadReads()` applies some formatting on BAM files before dependency can be
#' further analyzed. It takes as input the path to a folder with all the reads
#' of one condition inside, and consider each BAM file to be one replicate. It
#' offers to remove:
#'
#' * reads not mapping any maturation sites
#' * reads mapping multiple locations
#' * reads longer than a provided length
#' * reads mapping outside of a given range
#'
#' @param path_to_BAM `String` containing the path to the BAM files folder. Put
#' all the replicates in the same folder (along with the corresponding .bai
#' files).
#' @param path_to_maturation_sites `String` containing the path to the folder
#' with all maturation sites files, *e.g.* in gff file format.
#'
#' @param max_length `Integer` indicating the maximal authorized length for a
#' read. If specified, reads longer than this value will be removed.
#' @param genome_range Optionnal `Numeric` vector of length 2. If specified,
#' reads mapping outside this region will be removed.
#' @param remove_duplicates `Logical` indicating whether or not reads mapping
#' multiple locations will be removed. Default to `FALSE`.
#' @param verbose `Integer` indicating the level of wordiness:
#'   * 0 = muted execution
#'   * 1 = some information about the execution is displayed in the console
#'
#' @return A list of `GRanges` objects with all the reads for one condition
#' with one given replicate in one object.
#'
#' @examples
#'
#' WT_reads <- loadReads(path_to_BAM, path_to_maturation_sites, max_length=2^13,
#'                       genome_range=c(1, 200000), remove_duplicates=TRUE,
#'                       verbose=1)
#'
#' @export
loadReads <- function(path_to_BAM,
                      path_to_maturation_sites,
                      max_length,
                      genome_range,
                      remove_duplicates=FALSE,
                      verbose=0) {

  # Retrieve file names
  if (verbose > 0) {writeLines("Looking for bam files...")}
  bam_files <- list.files(path_to_BAM,
                          pattern=".bam$")
  if (verbose > 0) {writeLines(
    paste0(length(bam_files), " replicates found.")
  )}

  # Retrieve all reads
  scan_param <- Rsamtools::ScanBamParam(
    what=c("seq", "qname", "cigar", "pos", "strand")
  )
  if (verbose > 0) {writeLines(
    "done.\nBuilding Granges object for replicate from file:"
  )}
  reads <- lapply(1:length(bam_files), FUN = function(replicate) {
    if (verbose > 0) {writeLines(
      paste0(replicate, "/", length(bam_files), "...")
    )}
    methods::as(GenomicAlignments::readGAlignments(
      file=paste0(path_to_BAM, bam_files[[replicate]]),
      param=scan_param
    ),
    "GRanges")
  })
  if (verbose > 0) {writeLines("done.")}

  # Remove reads out of the expected range
  if (!missing(genome_range)) {
    reads <- lapply(1:length(bam_files), FUN = function(replicate) {
      # Check the `seqnames` attribute
      seqname <- unique(as.vector(GenomeInfoDb::seqnames(reads[[replicate]])))
      if (length(seqname) > 1) {
        stop(paste0("Provided reads have non unique sequence names for file: ",
                    path_to_BAM,
                    bam_files[[replicate]]))
      }
      # Remove the out of bounds reads
      IRanges::subsetByOverlaps(
        reads[[replicate]],
        GenomicRanges::GRanges(seqnames=seqname,
                               ranges=IRanges::IRanges(
                                 genome_range[[1]]:genome_range[[2]])
                               )
      )
    })
  }



  return()
}
















