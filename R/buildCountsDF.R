#' Compute the DESeq2 input counts table
#'
#' @description
#' `buildCountsDF()` takes as input a data frame with informations about reads
#' maturation state at each site (as output by [getStates()]) and returns a data
#' frame ready to be fed to DESeq2 to obtain comaturations.
#'
#' One row in the output data frame corresponds to a pair of maturation events.
#'
#' One column in the output data frame corresponds to a combination of the
#' maturation state at each site and a replicate.
#'
#' @inheritParams loadReads
#' @param df_with_all_events `List` containing data frames (one replicate in one
#' data frame) with the maturation state of each read at each maturation site
#' - ideally output by [getStates()].
#'
#' @return A `Data.frame` with the number of reads in each replicate having a
#' given maturation state for all pair of sites.
#'
#' @examples
#'
#' path_to_all_df       <- system.file("extdata",
#'                                     "getStates_output_ex.rds",
#'                                     package = "comaturationTracker")
#'
#' df_with_all_events   <- readRDS(path_to_all_df)
#'
#' buildCountsDF_output <- buildCountsDF(df_with_all_events)
#' counts               <- buildCountsDF_output[[1]]
#' conditions           <- buildCountsDF_output[[2]]
#'
#' @export
buildCountsDF <- function(df_with_all_events,
                          verbose=0) {

  nb_rep <- length(df_with_all_events)

  # Building names of columns for counts data frame
  state_names <- expand.grid(list(c("AB", "Ab", "aB", "ab"),
                                  c(1:nb_rep)))
  state_names <- paste0(state_names$Var1,
                        state_names$Var2)

  # Building names of rows for counts data frame
  nb_evt      <- ncol(df_with_all_events[[1]])
  all_combi   <- t(utils::combn(nb_evt, 2))
  pairs_names <- paste0(all_combi[,1], " vs ", all_combi[,2])

  # Build counts data frame
  pairs_of_sites <- utils::combn(nb_evt, 2, simplify=F)
  allowed_levels <- c("True", "False", "Err")

  if (verbose > 0) {
    rep_ind <- 0
    writeLines("Building counts data frame for replicate")
  }

  counts <- do.call(
    cbind,
    lapply(
      df_with_all_events,
      FUN = function(replicate_states) {
        if (verbose > 0) {
          rep_ind <<- rep_ind + 1
          writeLines(paste0(rep_ind, "/", nb_rep, "..."))
        }

        do.call(
          rbind,
          lapply(
            pairs_of_sites,
            FUN = function(pair_of_sites) {
              ctg_table <- table(
                factor(replicate_states[,pair_of_sites[1]],
                       levels=allowed_levels),
                factor(replicate_states[,pair_of_sites[2]],
                       levels=allowed_levels)
              )
              data.frame(ctg_table[1,1],
                         ctg_table[1,2],
                         ctg_table[2,1],
                         ctg_table[2,2])
            }
          )
        )
      }
    )
  )

  colnames(counts) <- state_names
  row.names(counts) <- pairs_names
  if (verbose > 0) {writeLines("done.")}

  # Build conditions data frame
  if (verbose > 0) {writeLines("Building conditions data frame...")}
  conditions <- data.frame(A=factor(c("True","True","False","False")),
                           B=factor(c("True","False","True","False")))
  conditions <- conditions[rep(1:4, nb_rep),]
  row.names(conditions) <- state_names
  conditions <- cbind(
    conditions,
    do.call(
      cbind,
      lapply(1:nb_rep, FUN = function(rep_ind) {
        states         <- rep(factor("False"), 4*nb_rep)
        levels(states) <- c("False", "True")
        states[((4*(rep_ind-1))+c(1:4))] <- factor("True")
        states           <- data.frame(states)
        colnames(states) <- paste0("R", rep_ind)
        states
      })
    )
  )
  if (verbose > 0) {writeLines("done.")}

  return(list(counts, conditions))
}
