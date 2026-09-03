################################################################################
#' @title ramchoice: Revealed Preference and Attention Analysis in Random Limited Attention Models
#'
#' @description Preferences and attention are important for understanding
#'   decision making, conducting welfare analysis, and providing robust policy
#'   recommendations. Decision makers may not pay full attention to all
#'   available alternatives, however, which can invalidate standard revealed
#'   preference analysis.
#'
#'   This package implements identification, estimation, inference, and
#'   specification procedures for the Random Attention Model of Cattaneo, Ma,
#'   Masatlioglu, and Suleymanov (2020; \doi{10.1086/706861}) and the
#'   Attention Overload Model of
#'   \href{https://arxiv.org/abs/2110.10650}{Cattaneo, Cheung, Ma, and
#'   Masatlioglu (2026)}.
#'
#'   The principal RAM and homogeneous-AOM interfaces are
#'   \code{\link{revealPref}}, \code{\link{ramTest}},
#'   \code{\link{revealAtte}}, \code{\link{revealPrefModel}},
#'   \code{\link{aomModel}}, \code{\link{aomTest}}, and
#'   \code{\link{aomIdentify}}. The heterogeneous list-based AOM interfaces are
#'   \code{\link{hlaoModel}}, \code{\link{hlaoTest}},
#'   \code{\link{hlaoNoPITest}}, \code{\link{hlaoEvent}}, and
#'   \code{\link{hlaoRankings}}. Data preparation and simulation utilities
#'   include \code{\link{sumData}}, \code{\link{genMat}},
#'   \code{\link{logitAtte}}, and \code{\link{logitSimu}}. The legacy
#'   \code{\link{rAtte}} interface and simulated \code{\link{ramdata}} dataset
#'   are retained for compatibility and illustration.
#'
#' @references
#' M. D. Cattaneo, X. Ma, Y. Masatlioglu, and E. Suleymanov (2020).
#' A Random Attention Model. \emph{Journal of Political Economy} 128(7):
#' 2796--2836. \doi{10.1086/706861}
#'
#' M. D. Cattaneo, P. H. Y. Cheung, X. Ma, and Y. Masatlioglu (2026).
#' \href{https://arxiv.org/abs/2110.10650}{Attention Overload}. Working paper.
#'
#' @author
#' Matias D. Cattaneo (maintainer), Princeton University.
#' \email{matias.d.cattaneo@gmail.com}.
#'
#' Paul Cheung, University of Maryland. \email{hycheung@umd.edu}
#'
#' Xinwei Ma, University of California San Diego. \email{x1ma@ucsd.edu}
#'
#' Yusufcan Masatlioglu, University of Maryland. \email{yusufcan@umd.edu}
#'
#' Elchin Suleymanov, Purdue University. \email{esuleyma@purdue.edu}
#'
#' @importFrom stats quantile
#' @importFrom MASS mvrnorm
#' @importFrom stats rmultinom
#' @importFrom utils combn
#'
#' @aliases ramchoice-package
"_PACKAGE"
