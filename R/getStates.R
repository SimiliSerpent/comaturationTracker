#' Retrieve the maturation states of a set of reads for a set of maturations
#'
#' @description
#' `getStates()` takes as input a set of reads (as output by [loadReads()]) and
#' a set of maturation events (e.g. in `gff` format) and returns a data frame
#' giving the maturation state of each read at all the maturation events sites.
#'
#' @param all_reads `List` containing `GRanges` objects (one replicate in one
#' object) - ideally output by [loadReads()].
#' @param path_to_maturation_sites `String` containing the path to the folder
#' with all maturation sites files, *e.g.* in gff file format. The type must
#' be `editing` or `intron`.
#'
#' @return A `List` of data frames with booleans indicating whether or not
#' reads are matured for each site (one data frame per replicate).
#'
#' @examples
#'
#' path_to_all_reads        <- system.file("extdata",
#'                                         "loadReads_output_ex.rds",
#'                                         package = "comaturationTracker")
#'
#' all_reads                <- readRDS(path_to_all_reads)
#'
#' path_to_maturation_sites <- system.file("extdata",
#'                                         "sites_files",
#'                                         package = "comaturationTracker")
#'
#' maturation_states        <- getStates(all_reads,
#'                                      path_to_maturation_sites)
#'
#' @export
getStates <- function(all_reads,
                     path_to_maturation_sites) {

  # Run through each replicate...
  df_with_all_events <- lapply(all_reads, function(reads) {

    # Store reads attributes in a data frame for speed
    df_4_speed <- data.frame(
      sequences = as.character(reads$seq),
      cigar     = reads$cigar,
      start     = GenomicRanges::start(reads),
      strand    = as.character(GenomicRanges::strand(reads))
    )

    # Set up maturation sites - along with a data frame for speed
    mat_sites <- getMatSites(path_to_maturation_sites)
    mat_df    <- data.frame(start = GenomicRanges::start(mat_sites),
                            end   = GenomicRanges::end(mat_sites),
                            id    = paste0(mat_sites$type, "_",
                                           mat_sites$file_id, "_",
                                           mat_sites$site_id),
                            type  = mat_sites$type)

    # Get the sites mapped by each read
    overlaps <- GenomicAlignments::findOverlaps(reads, mat_sites)

    # Run through each replicate...
    do.call(rbind, lapply(1:length(reads), function(read_ind) {

      # Retrieve the read infos for fast accession
      df_4_speed_single_read   <- df_4_speed[read_ind,]
      overlapping_sites_ind    <- overlaps[overlaps@from==read_ind,]@to
      basePosition_single_read <- getBasePos(df_4_speed_single_read$cigar,
                                             df_4_speed_single_read$start)

      # Get the state of maturation at each site mapped by the read
      read_events_state <- lapply(overlapping_sites_ind, function(site_ind) {
        if (mat_df$type[[site_ind]] == "editing") {
          state <- info.edited(df_4_speed_single_read$sequences,
                               mat_df$start[[site_ind]],
                               basePosition_single_read)
        } else if (mat_df$type[[site_ind]] == "intron") {
          state <- info.spliced(mat_df[site_ind,],
                                basePosition_single_read)
        }
        state <- convert.info(state,
                              as.character(df_4_speed_single_read$strand))
        state
      })
      names(read_events_state)     <- mat_df$id[overlapping_sites_ind]
      read_all_events_state        <- character(dim(mat_df)[1])
      names(read_all_events_state) <- mat_df$id
      read_all_events_state[match(names(read_events_state),
                                  mat_df$id)] <- unlist(read_events_state)
      read_all_events_state
    }))
  })

  return(df_with_all_events)
}

# @brief This functions builds a list with type and positions of all genome-
# aligned base of a read
# function stolen from
# https://forgemia.inra.fr/guillem.rigaill/nanopore_chloro/-/blob/master/...
#   ...first_draft/Code_extract_utils.R
#
# @param cigar: character - cigar of the read
# @param start: integer - start position of the alignment on the genome
#
# @return all_query_positions: list of positions for each aligned base of the
# read
getBasePos <- function(cigar, start) {
  occurences      <- as.integer(stringr::str_extract_all(cigar, "\\d+")[[1]])
  align_events    <- stringr::str_extract_all(cigar, "[A-Z]+")[[1]]

  notH         <- which(align_events != "H") ## not present in the query
  occurences   <- occurences[notH]
  align_events <- align_events[notH]

  element.length          <- list(M=1, I=0, S=0, D=1, N=1)
  all_cigar_events        <- rep(align_events, occurences)
  all_cigar_events_length <- unlist(element.length[all_cigar_events])
  all_cigar_positions     <- cumsum(all_cigar_events_length) + start - 1
  all_cigar_positions[all_cigar_events == "S" | all_cigar_events == "I"] <- -1
  all_query_positions     <- all_cigar_positions[all_cigar_events != "D" &
                                                 all_cigar_events != "N"]
  return(all_query_positions)
}

# @brief This functions retrieves the value of the base of the read at the
# position of the edition site
# function stolen from
# https://forgemia.inra.fr/guillem.rigaill/nanopore_chloro/-/blob/master/...
#   ...first_draft/Code_extract_utils.R
#
# @param sequence: character - sequence of the read
# @param pos_site: integer - start position of the site on the genome
# @param basePosition: list - output of the getBasePos function for
# the same read
#
# @return data.frame object: data frame with value of the base at the site
# position
info.edited <- function(sequence, pos_site, basePosition){

  posToRead <- which(basePosition == pos_site)
  if(length(posToRead) == 0){
    return(data.frame(type="edition", read="not-read"))
  } else {
    if(length(posToRead) == 1){
      return(data.frame(type="edition",
                        read=substr(sequence, posToRead, posToRead)))
    } else {
      return(data.frame(type="edition", read="er::mapped-several-times"))
    }
  }
}

# @brief This functions retrieves information on the state of the read at the
# position of the intron
# function stolen from
# https://forgemia.inra.fr/guillem.rigaill/nanopore_chloro/-/blob/master/...
#   ...first_draft/Code_extract_utils.R
#
# @param target_splicing_site: DataFrame - row of mat_df describing the
# intron
# @param basePosition: integer - start position of the site on the genome
#
# @return data.frame object: data frame with information on the read at the
# intron location
info.spliced <- function(target_splicing_site, basePosition){
  ## excluding -1
  bP <- basePosition[basePosition > -1]
  ## finding the genomic min position in basePosition (-1 excluded)
  pos_min <- min(bP)
  ## finding the genomic max position in basePosition (-1 excluded)
  pos_max <- max(bP)

  ## check whether the code should depend on the strand
  is.after.min <- (pos_min <= target_splicing_site$start)
  is.before.max <- (pos_max >= target_splicing_site$end)

  prc.intron.read <- length(which(bP >= target_splicing_site$start
                                  & bP <= target_splicing_site$end))/
    (target_splicing_site$end-target_splicing_site$start)
  return(data.frame(type="intron",
                    prc.base=prc.intron.read,
                    is.after.min=is.after.min,
                    is.before.max=is.before.max))
}

# @brief This functions retrieves the maturation state of a read for a given
# maturation event
# function stolen from
# https://forgemia.inra.fr/guillem.rigaill/nanopore_chloro/-/blob/master/...
#   ...first_draft/Code_extract_utils.R
#
# @param x: DataFrame - output of info.edited or info.spliced
# @param strand: character - strand of the read
# @param thrs.intron: numeric - maximal fraction of the intron found in the
# read to consider the transcript has been spliced (default: 0.1)
#
# @return y: state of the maturation
convert.info <- function(x, strand, thrs.intron=0.1){

  if(x$type == "edition"){
    if(strand=="+"){
      if(x$read == "A") y <- "Err"
      if(x$read == "G") y <- "Err"
      if(x$read == "C") y <- "False"
      if(x$read == "T") y <- "True"
    }
    if(strand=="-"){
      if(x$read == "A") y <- "True"
      if(x$read == "G") y <- "False"
      if(x$read == "C") y <- "Err"
      if(x$read == "T") y <- "Err"
    }
    if(x$read == "not-read") y <- "Err"
    if(x$read == "er::mapped-several-times") y <- "Err"
  }
  ## INTRON
  if(x$type == "intron"){
    if( !(x$is.after.min & x$is.before.max) ) { # CHECK THAT SOME BASE ARE
                                                # BEFORE AND AFTER THE INTRONS
      y <- "Err"
    } else {
      if(x$prc.base <= thrs.intron ){  ## BELOW THE THRESHOLD => SPLICED
        y <- "True"
      } else {   ## ABOVE THE THRESHOLD => NOT-SPLICED
        y <- "False"
      }
    }
  }
  return(y)
}














