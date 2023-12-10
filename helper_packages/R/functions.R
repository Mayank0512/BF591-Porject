library(tidyverse)
#' Extract timepoint from sample name
#'
#' @param x A character vector of sample names
#' @return A character vector of timepoints
#' @export
timepoint_from_sample <- function(x) {
  return(substring(x, 2, regexpr("_", x) - 1))
}

#' Extract replicate from sample name
#'
#' @param x A character vector of sample names
#' @return A character vector of replicates
#' @export
sample_replicate <- function(x) {
  return(substring(x, regexpr("_", x) + 1))
}

#' Generate meta-information from sample names
#'
#' @param sample_names A character vector of sample names
#' @return A tibble containing sample_name, timepoint, and replicate columns
#' @export
meta_info_from_labels <- function(sample_names) {
  replicate <- sapply(sample_names, FUN = sample_replicate)
  timepoint <- sapply(sample_names, FUN = timepoint_from_sample)
  final <- tibble(
    sample_name = sample_names,
    timepoint = timepoint,
    replicate = replicate
  )
  return(final)
}

