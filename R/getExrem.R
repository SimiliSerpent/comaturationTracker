#' Store strand, extremities and length of a set of reads in a data frame
#'
#' @description
#' `getExtrem()` takes as input a set of reads (as output by [loadReads()]) and
#' returns a list of data frames giving informations about the reads:
#' * strand
#' * start and end positions
#' * length
#'
#' @inheritParams getStates
#'
#' @return A `List` of data frames with reads informations.
#'
#' @examples
#'
#' path_to_all_reads <- system.file("extdata",
#'                                  "loadReads_output_ex.rds",
#'                                  package = "comaturationTracker")
#'
#' all_reads         <- readRDS(path_to_all_reads)
#'
#' all_extremities   <- getExtrem(all_reads)
#'
#' @export
getExtrem <- function(all_reads) {

  all_strand <- lapply(all_reads,
                       FUN=function(reads) {
                         as.character(GenomicRanges::strand(reads))
                       })
  all_start  <- lapply(all_reads,
                       FUN=function(reads) {
                         as.integer(GenomicRanges::start(reads))
                       })
  all_end    <- lapply(all_reads,
                       FUN=function(reads) {
                         as.integer(GenomicRanges::end(reads))
                       })
  all_extremities <- lapply(1:length(all_reads),
                            FUN=function(rep_ind) {
                              data.frame(strand = all_strand[[rep_ind]],
                                         start = all_start[[rep_ind]],
                                         end = all_end[[rep_ind]],
                                         length = (all_end[[rep_ind]] -
                                                     all_start[[rep_ind]]))
                            })

  return(all_extremities)
}
