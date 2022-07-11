#' Retrieve and format reads
#'
#' @description
#' `loadReads()` applies some formatting on BAM files before dependency can be
#' further analyzed. It takes as input the path to a folder with all the reads
#' of one condition inside, and consider each BAM file to be one replicate. It
#' offers to remove:
#'
#' * reads not mapping any maturation sites
#' * reads not mapping a given chromosome
#' * reads longer than a provided length
#' * reads mapping outside of a given range
#' * reads mapping multiple locations
#'
#' @param path_to_BAM `String` containing the path to the BAM files folder. Put
#' all the replicates in the same folder (along with the corresponding .bai
#' files).
#' @param path_to_maturation_sites Optionnal `String` containing the path to
#' the folder with all maturation sites files, *e.g.* in gff file format.
#' @param seqname `String` containing the name of the chromosome of interest.
#' @param max_length Optionnal `Integer` indicating the maximal authorized
#' length for a read. If specified, strictly longer reads will be removed.
#' @param genome_range Optionnal `Numeric` vector of length 2. If specified,
#' reads mapping outside this region will be removed.
#' @param remove_duplicates `Logical` indicating whether or not reads mapping
#' multiple locations will be removed. Default to `TRUE`.
#' @param verbose `Integer` indicating the level of wordiness:
#'   * 0 = muted execution
#'   * 1 = some information about the execution is displayed in the console
#'
#' @return A list of `GRanges` objects with all the reads for one condition
#' with one given replicate in one element.
#'
#' @examples
#'
#' path_to_BAM              <- system.file("extdata",
#'                                         "bam_files",
#'                                         package = "comaturationTracker")
#'
#' path_to_maturation_sites <- system.file("extdata",
#'                                         "sites_files",
#'                                         package = "comaturationTracker")
#'
#' WT_reads                 <- loadReads(path_to_BAM,
#'                                       path_to_maturation_sites,
#'                                       seqname="Pt",
#'                                       max_length=2^13,
#'                                       genome_range=c(1, 200000),
#'                                       remove_duplicates=TRUE,
#'                                       verbose=1)
#'
#' @export
loadReads <- function(path_to_BAM,
                      path_to_maturation_sites=NULL,
                      seqname=NULL,
                      max_length=NULL,
                      genome_range=NULL,
                      remove_duplicates=TRUE,
                      verbose=0) {

  # Check the `seqname` attribute
  if (is.null(seqname)) {
    stop("Please provide seqname.")
  }

  # Retrieve file names
  if (verbose > 0) {writeLines("Looking for bam files...")}
  bam_files <- list.files(path_to_BAM,
                          pattern=".bam$")
  if (verbose > 0) {writeLines(
    paste0(length(bam_files), " replicates found.")
  )}

  # Set up ranges and chromosome criteria
  if (is.null(genome_range)) {
    genome_range <- c(1, 536870912)
  }
  whole_range <- IRanges::IRangesList(IRanges::IRanges(genome_range[[1]],
                                                       genome_range[[2]]))
  names(whole_range) <- seqname
  scan_param <- Rsamtools::ScanBamParam(
    which=whole_range,
    what=c("seq", "qname", "cigar", "pos", "strand")
  )

  # Set up maturation sites
  if (!is.null(path_to_maturation_sites)) {
    mat_files <- list.files(path_to_maturation_sites)
    if (length(mat_files) > 0) {
      mat_sites <- c(
        rtracklayer::import(paste0(path_to_maturation_sites,
                                   "/",
                                   mat_files[[1]]))
      )
      if (length(mat_files) > 1) {
        for (file_ind in 2:length(mat_files)) {
          mat_sites <- c(mat_sites,
                         rtracklayer::import(paste0(path_to_maturation_sites,
                                                    "/",
                                                    mat_files[[file_ind]])))
        }
      }
    }
    GenomeInfoDb::seqlevels(mat_sites) <- seqname
  }

  if (verbose > 0) {writeLines(
    "Building Granges object for replicate:"
  )}
  all_reads <- lapply(1:length(bam_files), FUN = function(replicate) {
    if (verbose > 0) {writeLines(
      paste0(replicate, "/", length(bam_files), "...")
    )}

    # Retrieve reads
    reads <- methods::as(GenomicAlignments::readGAlignments(
      file=paste0(path_to_BAM, "/", bam_files[[replicate]]),
      param=scan_param
    ), "GRanges")

    # Remove reads not overlapping the maturation events sites
    if (!is.null(path_to_maturation_sites) && length(mat_files) > 0) {
      reads <- IRanges::subsetByOverlaps(reads, mat_sites)
    }

    # Remove reads exceeding a given range
    if (!is.null(max_length)) {
      reads <- reads[reads@ranges@width <= max_length]
    }

    # Remove reads mapping multiple locations
    if (remove_duplicates) {
      unique_qnames <- names(which(table(reads$qname) == 1))
      reads <- reads[reads$qname %in% unique_qnames]
    }

  })
  if (verbose > 0) {writeLines("done.")}

  return(all_reads)
}
















