################################################################################
#' @title ramdata: Simulated Choice Data
#'
#' @description The file contains a standard choice data of 9,000 observations.
#'   There are five alternatives in the grand set.
#'
#' See \code{\link{revealPref}} for revealed preference analysis, and \code{\link{revealAtte}}
#'   for revealed attention. \code{\link{sumData}} is a low-level function that computes summary
#'   statistics, and \code{\link{genMat}} generates constraint matrices subject to given preferences.
#'
#' @format
#' \describe{
#'   \item{menu}{Numeric matrix of 0s and 1s, choice problems (1 indicates an alternative in the choice problem and 0 otherwise).}
#'   \item{choice}{Numeric matrix of 0s and 1s, choices (1 indicates an alternative being chosen).}
#' }
#'
#' @docType data
#' @name ramdata
#' @aliases menu choice
NULL
