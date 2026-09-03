################################################################################
#' @title Generate Summary Statistics
#'
#' @description \code{sumData} generates summary statistics. Given a collection of
#'   choice problems and corresponding choices, \code{sumData} calculates the
#'   number of occurrences of each choice problem, as well as the empirical choice
#'   probabilities.
#'
#' This function is embedded in \code{\link{revealPref}}.
#'
#' @param menu Numeric matrix of 0s and 1s, the collection of choice problems.
#' @param choice Numeric matrix of 0s and 1s, the collection of choices.
#'
#' @return
#' \item{sumMenu}{Summary of choice problems, with repetitions removed.}
#' \item{sumProb}{Estimated choice probabilities as sample averages for different choice problems.}
#' \item{sumN}{Effective sample size for each choice problem.}
#' \item{sumMsize}{Size of each choice problem.}
#' \item{sumProbVec}{Estimated choice probabilities as sample averages, collapsed into a column vector.}
#' \item{Sigma}{Estimated variance-covariance matrix for the choice rule, scaled by relative sample sizes.}
#'
#' @references
#' M. D. Cattaneo, X. Ma, Y. Masatlioglu, and E. Suleymanov (2020). A Random Attention Model. \emph{Journal of Political Economy} 128(7): 2796--2836. \doi{10.1086/706861}
#'
#' M. D. Cattaneo, P. H. Y. Cheung, X. Ma, and Y. Masatlioglu (2026). \href{https://arxiv.org/abs/2110.10650}{Attention Overload}. Working paper.
#'
#' @author
#' Matias D. Cattaneo (maintainer), Princeton University. \email{matias.d.cattaneo@gmail.com}.
#'
#' Paul Cheung, University of Maryland. \email{hycheung@umd.edu}
#'
#' Xinwei Ma, University of California San Diego. \email{x1ma@ucsd.edu}
#'
#' Yusufcan Masatlioglu, University of Maryland. \email{yusufcan@umd.edu}
#'
#' Elchin Suleymanov, Purdue University. \email{esuleyma@purdue.edu}
#'
#' @examples
#' # Load data
#' data(ramdata)
#'
#' # Generate summary statistics
#' summaryStats <- sumData(ramdata$menu, ramdata$choice)
#' nrow(summaryStats$sumMenu)
#' min(summaryStats$sumN)
#'
#' summaryStats$sumMenu[1, ]
#' summaryStats$sumProb[1, ]
#' summaryStats$sumN[1]
#'
#' @export
sumData <- function(menu, choice) {

  ################################################################################
  # Generate Summary Statistics
  ################################################################################

  sumMenu  <- matrix(NA, ncol=ncol(menu), nrow=nrow(menu)) # collection of distinct menus
  sumProb  <- matrix(NA, ncol=ncol(menu), nrow=nrow(menu)) # collection of estimated choice probabilities
  sumN     <- matrix(NA, ncol=1, nrow=nrow(menu)) # collection of effective sample sizes
  sumMsize <- matrix(NA, ncol=1, nrow=nrow(menu)) # collection of menu sizes

  j <- 1
  while (nrow(menu) > 0) {
    # add new line to sumMenu, for a new menu
    sumMenu[j, ] <- menu[1, ]

    # add new line to sumMsize for the size of the new menu
    sumMsize[j, ] <- sum(menu[1, ])

    # enumerate and detect the same menus
    if (nrow(menu) >= 2) {
      i_menu <- apply(menu, MARGIN=1, FUN=function(x) all(x == menu[1, ]))
    } else {
      i_menu <- TRUE
    }

    # determine the effective sample size
    sumN[j, ] <- sum(i_menu)

    # calculate choice probabilities
    sumProb[j, ] <- apply(choice[i_menu, , drop=FALSE], MARGIN=2, FUN=mean)

    # delete rows used
    menu <- menu[!i_menu, , drop=FALSE]
    choice <- choice[!i_menu, , drop=FALSE]

    j <- j + 1
  }

  j <- j - 1
  sumMenu  <- sumMenu[1:j, ]
  sumProb  <- sumProb[1:j, ]
  sumN     <- sumN[1:j, ]
  sumMsize <- sumMsize[1:j, ]


  sumProbVec  <- matrix(0, nrow=sum(sumMsize), ncol=1)
  Sigma       <- matrix(0, nrow=sum(sumMsize), ncol=sum(sumMsize))

  j = 0
  for (i in 1:nrow(sumProb)) {
    temp <- sumProb[i, sumMenu[i, ] == 1]
    sumProbVec[(j+1):(j+sumMsize[i])] <- temp
    Sigma[(j+1):(j+sumMsize[i]), (j+1):(j+sumMsize[i])] <-
      (diag(temp) - tcrossprod(temp)) / sumN[i] * sum(sumN)

    j <- j + sumMsize[i]
  }

  return (list(sumMenu=sumMenu, sumProb=sumProb, sumN=sumN, sumMsize=sumMsize,
               sumProbVec=sumProbVec, Sigma=Sigma))

}

# Restriction matrices are often reused across Monte Carlo replications. Keep
# only the most recent construction so repeated calls do not grow package state.
.genMatCache <- new.env(parent = emptyenv())

################################################################################
#' @title Generate Constraint Matrices
#'
#' @description \code{genMat} generates constraint matrices for a range of preference orderings according to
#'   (i) the monotonic attention assumption proposed by Cattaneo, Ma, Masatlioglu, and Suleymanov (2020),
#'   (ii) the attention overload assumption proposed by Cattaneo, Cheung, Ma, and Masatlioglu (2021),
#'   and (iii) the attentive-at-binaries restriction.
#'
#' This function is embedded in \code{\link{revealPref}}.
#'
#' @param sumMenu Numeric matrix, summary of choice problems, returned by \code{\link{sumData}}.
#' @param sumMsize Numeric matrix, summary of choice problem sizes, returned by \code{\link{sumData}}.
#' @param pref_list Numeric matrix, each row corresponds to one preference. For example, \code{c(2, 3, 1)} means
#'   2 is preferred to 3 and to 1. When set to \code{NULL}, the default, \code{c(1, 2, 3, ...)},
#'   will be used.
#' @param RAM Boolean, whether the restrictions implied by the random attention model of
#'   Cattaneo, Ma, Masatlioglu, and Suleymanov (2020) should be incorporated, that is, their monotonic attention assumption (default is \code{TRUE}).
#' @param AOM Boolean, whether the restrictions implied by the attention overload model of
#'   Cattaneo, Cheung, Ma, and Masatlioglu (2021) should be incorporated, that is, their attention overload assumption (default is \code{TRUE}).
#' @param limDataCorr Boolean, whether assuming limited data (default is \code{TRUE}). When set to
#'   \code{FALSE}, will assume all choice problems are observed. This option only applies when \code{RAM} is set to \code{TRUE}.
#' @param attBinary Numeric, between 1/2 and 1 (default is \code{1}), whether additional restrictions (on the attention rule)
#'   should be imposed for binary choice problems (i.e., attentive at binaries).
#'
#' @return
#' \item{R}{Matrices of constraints, stacked vertically.}
#' \item{ConstN}{The number of constraints for each preference, used to extract from \code{R}
#'   individual matrices of constraints.}
#'
#' @references
#' M. D. Cattaneo, X. Ma, Y. Masatlioglu, and E. Suleymanov (2020). A Random Attention Model. \emph{Journal of Political Economy} 128(7): 2796--2836. \doi{10.1086/706861}
#'
#' M. D. Cattaneo, P. H. Y. Cheung, X. Ma, and Y. Masatlioglu (2026). \href{https://arxiv.org/abs/2110.10650}{Attention Overload}. Working paper.
#'
#' @author
#' Matias D. Cattaneo (maintainer), Princeton University. \email{matias.d.cattaneo@gmail.com}.
#'
#' Paul Cheung, University of Maryland. \email{hycheung@umd.edu}
#'
#' Xinwei Ma, University of California San Diego. \email{x1ma@ucsd.edu}
#'
#' Yusufcan Masatlioglu, University of Maryland. \email{yusufcan@umd.edu}
#'
#' Elchin Suleymanov, Purdue University. \email{esuleyma@purdue.edu}
#'
#' @examples
#' # Load data
#' data(ramdata)
#'
#' # Generate summary statistics
#' summaryStats <- sumData(ramdata$menu, ramdata$choice)
#'
#' # Generate constraint matrices
#' constraints <- genMat(summaryStats$sumMenu, summaryStats$sumMsize)
#' constraints$ConstN
#' constraints$R[1:10, 1:10]
#'
#' @export
genMat <- function(sumMenu, sumMsize, pref_list = NULL, RAM = TRUE, AOM = TRUE, limDataCorr = TRUE, attBinary = 1) {

  # initializing preference, for the default
  if (length(as.vector(pref_list)) == 0) {
    pref_list <- matrix(1:ncol(sumMenu), nrow=1)
  }

  ################################################################################
  # Default
  ################################################################################

  if (nrow(sumMenu) <= 1) {
    return(list(R=matrix(NA, nrow=0, ncol=sum(sumMsize)), constN=0))
  }

  ################################################################################
  # Method specification
  ################################################################################

  # option RAM
  if (length(RAM) == 0) {
    RAM <- TRUE
  } else if (length(RAM) > 1 | !RAM[1]%in%c(TRUE, FALSE)) {
    stop("Option RAM incorrectly specified.\n")
  }

  # option AOM
  if (length(AOM) == 0) {
    AOM <- TRUE
  } else if (length(AOM) > 1 | !AOM[1]%in%c(TRUE, FALSE)) {
    stop("Option AOM incorrectly specified.\n")
  }

  # check if both RAM and AOM are FALSE
  if (!RAM & !AOM) {
    stop("At least one option, RAM or AOM, has to be used.\n")
  }

  cache_key <- list(
    sumMenu = sumMenu,
    sumMsize = sumMsize,
    pref_list = pref_list,
    RAM = RAM,
    AOM = AOM,
    limDataCorr = limDataCorr,
    attBinary = attBinary
  )
  if (exists("key", envir = .genMatCache, inherits = FALSE) &&
      identical(cache_key, .genMatCache$key)) {
    return(.genMatCache$value)
  }

  ################################################################################
  # Construction
  ################################################################################

  K <- ncol(sumMenu) # dimension of the problem

  R <- matrix(0, nrow=0, ncol=sum(sumMsize))

  ConstN <- rep(NA, nrow(pref_list))

  temp <- rep(0, sum(sumMsize))

  for (pref_index in 1:nrow(pref_list)) {
    R_temp <- matrix(0, nrow=0, ncol=sum(sumMsize))
    pref <- pref_list[pref_index, ] # current preference

    ### goal: construct constraints of the form
    # mu(a|S) - mu(a|T) <= 0
    ###

    # index of menus to be enumerated in the first layer, S
    index_menu_1 <- (1:length(sumMsize))[sumMsize > min(sumMsize)]
    for (i_menu_1 in index_menu_1) {
      menu_1 <- sumMenu[i_menu_1, ] # current menu, S

      # index of menus to be enumerated in the second layer / subsets of current menu, T
      # Case: without limited data correction, i.e. assume full observability; no AOM
      if (RAM & (!AOM) & (!limDataCorr)) {
        index_menu_2 <- (1:length(sumMsize))[apply(sumMenu, MARGIN=1, FUN=function(x) {all((menu_1 - x)>=0) & (sum(menu_1) - sum(x) == 1)})]
      } else {
      # Case: with limited data correction, i.e. assume full observability
        index_menu_2 <- (1:length(sumMsize))[apply(sumMenu, MARGIN=1, FUN=function(x) {all((menu_1 - x)>=0) & (sum(menu_1) > sum(x))})]
      }

      for (i_menu_2 in index_menu_2) {
        menu_2 <- sumMenu[i_menu_2, ] # a subset, T

        if (RAM) {
          if (limDataCorr | (sum(menu_1) - sum(menu_2) == 1)) {
            # decide elements a < S-T in S and T
            i_element <- K # index for the lower contour element a
            while (menu_1[pref[i_element]] == menu_2[pref[i_element]]) { # if the index never enters S - T, so that remains in lower contour set
              if (menu_1[pref[i_element]] == 1) { # if in S (so also in T)
                pos_1 <- sum(sumMsize[1:i_menu_1]) - sumMsize[i_menu_1] + sum(menu_1 & ((1:K) <= pref[i_element]))
                pos_2 <- sum(sumMsize[1:i_menu_2]) - sumMsize[i_menu_2] + sum(menu_2 & ((1:K) <= pref[i_element]))
                temp <- temp * 0
                temp[pos_1] <- 1; temp[pos_2] <- -1
                R_temp <- rbind(R_temp, temp)
              }
              i_element <- i_element - 1
            }
          }
        }

        if (AOM) {
          pos_2 <- c()
          for (i_element in 1:K) {
            if (menu_2[pref[i_element]] == 1) { # so that this element is in T
              pos_1 <- sum(sumMsize[1:i_menu_1]) - sumMsize[i_menu_1] + sum(menu_1 & ((1:K) <= pref[i_element]))
              pos_2 <- c(pos_2, sum(sumMsize[1:i_menu_2]) - sumMsize[i_menu_2] + sum(menu_2 & ((1:K) <= pref[i_element])))
              if (length(pos_2) < sum(menu_2)) { # to rule out the degenerate case where the entire T is an upper countour set
                temp <- temp * 0
                temp[pos_1] <- 1; temp[pos_2] <- -1
                R_temp <- rbind(R_temp, temp)
              }
            }
          }
        }
      }
    }

    # now add constraints for attentive at binaries
    index_menu_1 <- (1:length(sumMsize))[sumMsize == 2]
    if (length(index_menu_1) > 0 & attBinary < 1) {
      for (i_menu_1 in index_menu_1) { # enumerate all menus of size 2
        menu_1 <- sumMenu[i_menu_1, ] # current menu
        pos_1 <- sum(sumMsize[1:i_menu_1]) - 1
        pos_2 <- sum(sumMsize[1:i_menu_1])
        temp <- temp * 0
        if (which(pref == which(menu_1 == 1)[1]) < which(pref == which(menu_1 == 1)[2])) {
          temp[pos_1] <- -1; temp[pos_2] <- (1-attBinary)/attBinary
        } else {
          temp[pos_1] <- (1-attBinary)/attBinary; temp[pos_2] <- -1
        }
        R_temp <- rbind(R_temp, temp)
      }
    }

    # combine them
    R <- rbind(R, R_temp)
    ConstN[pref_index] <- nrow(R_temp)
  }

  rownames(R) <- NULL
  result <- list(R=R, ConstN=ConstN)
  .genMatCache$key <- cache_key
  .genMatCache$value <- result
  return(result)
}

################################################################################
#' @title Compute Choice Probabilities and Attention Frequencies for the Logit Attention Rule
#'
#' @description \code{logitAtte} computes choice probabilities and attention frequencies for the logit attention rule
#'   considered by Brady and Rehbeck (2016). To be specific, for a choice problem \code{S} and its subset \code{T}, the attention that \code{T}
#'   attracts is assumed to be proportional to its size: \code{|T|^a}, where \code{a} is a parameter that one can specify. It will be assumed that
#'   the first alternative is the most preferred, and that the last alternative is the least preferred.
#'
#'   This function is useful for replicating the simulation results in Cattaneo, Ma, Masatlioglu, and Suleymanov (2020; \doi{10.1086/706861})
#'   and \href{https://arxiv.org/abs/2110.10650}{Cattaneo, Cheung, Ma, and Masatlioglu (2026)}.
#'
#' @param mSize Positive integer, size of the choice problem.
#' @param a Numeric, the parameter of the logit attention rule.
#'
#' @return
#' \item{choiceProb}{The vector of choice probabilities.}
#' \item{atteFreq}{The attention frequency.}
#'
#' @references
#' R. L. Brady and J. Rehbeck (2016). Menu-Dependent Stochastic Feasibility. \emph{Econometrica} 84(3): 1203-1223. \doi{10.3982/ECTA12694}
#'
#' M. D. Cattaneo, X. Ma, Y. Masatlioglu, and E. Suleymanov (2020). A Random Attention Model. \emph{Journal of Political Economy} 128(7): 2796--2836. \doi{10.1086/706861}
#'
#' M. D. Cattaneo, P. H. Y. Cheung, X. Ma, and Y. Masatlioglu (2026). \href{https://arxiv.org/abs/2110.10650}{Attention Overload}. Working paper.
#'
#' @author
#' Matias D. Cattaneo (maintainer), Princeton University. \email{matias.d.cattaneo@gmail.com}.
#'
#' Paul Cheung, University of Maryland. \email{hycheung@umd.edu}
#'
#' Xinwei Ma, University of California San Diego. \email{x1ma@ucsd.edu}
#'
#' Yusufcan Masatlioglu, University of Maryland. \email{yusufcan@umd.edu}
#'
#' Elchin Suleymanov, Purdue University. \email{esuleyma@purdue.edu}
#'
#' @examples
#' logitAtte(mSize = 5, a = 2)
#'
#' @export
logitAtte <- function(mSize = NULL, a = NULL) {

  ################################################################################
  # Some error handling
  ################################################################################

  if (length(mSize) == 0) {
    mSize <- 5
  } else {
    mSize <- as.integer(mSize)
    if (length(mSize) > 1 | mSize[1] <= 0) {
      stop("Option mSize incorrectly specified.\n")
    }
  }

  if (length(a) == 0) {
    a <- 2
  } else {
    if (length(a) > 1 | !is.numeric(a[1])) {
      stop("Option a incorrectly specified.\n")
    }
  }

  # find the (unscaled) choice probabilities
  choiceProb <- rep(0, mSize)
  for (i in 1:mSize) { # position of the alternative
    for (j in 1:(mSize-i+1)) { # size of the consideration set
      choiceProb[i] <- choiceProb[i] + choose(mSize-i, j-1) * (j^a)
    }
  }

  # find the (unscaled) attention frequencies
  atteFreq <- 0
  total <- 0
  for (j in 1:mSize) { # size of the consideration set
    atteFreq <- atteFreq + choose(mSize-1, j-1) * (j^a)
    total <- total + choose(mSize, j) * (j^a)
  }

  choiceProb <- choiceProb / total
  atteFreq <- atteFreq / total

  return(list(choiceProb=choiceProb, atteFreq=atteFreq))
}

################################################################################
#' @title Choice Data Simulation Following the Logit Attention Rule
#'
#' @description \code{logitSimu} simulates choice data according to the logit attention rule
#'   considered by Brady and Rehbeck (2016). To be specific, for a choice problem \code{S} and its subset \code{T}, the attention that \code{T}
#'   attracts is assumed to be proportional to its size: \code{|T|^a}, where \code{a} is a parameter that one can specify. It will be assumed that
#'   the first alternative is the most preferred, and that the last alternative is the least preferred.
#'
#'   This function is useful for replicating the simulation results in Cattaneo, Ma, Masatlioglu, and Suleymanov (2020; \doi{10.1086/706861})
#'   and \href{https://arxiv.org/abs/2110.10650}{Cattaneo, Cheung, Ma, and Masatlioglu (2026)}.
#'
#' @param n Positive integer, the effective sample size for each choice problem.
#' @param uSize Positive integer, total number of alternatives.
#' @param mSize Positive integer, size of the choice problem.
#' @param a Numeric, the parameter of the logit attention rule.
#'
#' @return
#' \item{menu}{The choice problems.}
#' \item{choice}{The simulated choices.}
#'
#' @references
#' R. L. Brady and J. Rehbeck (2016). Menu-Dependent Stochastic Feasibility. \emph{Econometrica} 84(3): 1203-1223. \doi{10.3982/ECTA12694}
#'
#' M. D. Cattaneo, X. Ma, Y. Masatlioglu, and E. Suleymanov (2020). A Random Attention Model. \emph{Journal of Political Economy} 128(7): 2796--2836. \doi{10.1086/706861}
#'
#' M. D. Cattaneo, P. H. Y. Cheung, X. Ma, and Y. Masatlioglu (2026). \href{https://arxiv.org/abs/2110.10650}{Attention Overload}. Working paper.
#'
#' @author
#' Matias D. Cattaneo (maintainer), Princeton University. \email{matias.d.cattaneo@gmail.com}.
#'
#' Paul Cheung, University of Maryland. \email{hycheung@umd.edu}
#'
#' Xinwei Ma, University of California San Diego. \email{x1ma@ucsd.edu}
#'
#' Yusufcan Masatlioglu, University of Maryland. \email{yusufcan@umd.edu}
#'
#' Elchin Suleymanov, Purdue University. \email{esuleyma@purdue.edu}
#'
#' @examples
#' set.seed(42)
#' logitSimu(n = 5, uSize = 6, mSize = 5, a = 2)
#'
#' @export
logitSimu <- function(n, uSize, mSize, a) {
  ################################################################################
  # Some error handling
  ################################################################################
  if (length(n) == 0) {
    n <- 20
  } else {
    n <- as.integer(n)
    if (length(n) > 1 | n[1] <= 0) {
      stop("Option n incorrectly specified.\n")
    }
  }

  if (length(uSize) == 0) {
    uSize <- 6
  } else {
    uSize <- as.integer(uSize)
    if (length(uSize) > 1 | uSize[1] <= 0) {
      stop("Option uSize incorrectly specified.\n")
    }
  }

  if (length(mSize) == 0) {
    mSize <- 5
  } else {
    mSize <- as.integer(mSize)
    if (length(mSize) > 1 | mSize[1] <= 0) {
      stop("Option mSize incorrectly specified.\n")
    }
  }

  if (uSize < mSize) {
    stop("Option uSize incorrectly specified: it cannot be smaller than mSize.\n")
  }

  if (length(a) == 0) {
    a <- 2
  } else {
    if (length(a) > 1 | !is.numeric(a[1])) {
      stop("Option a incorrectly specified.\n")
    }
  }

  # determine population choice rule
  choiceProb <- logitAtte(mSize = mSize, a = a)$choiceProb

  # initialize
  allMenus <- t(combn(uSize, mSize))  # all possible subsets
  menu <- choice <- matrix(0, nrow=n*nrow(allMenus), ncol=uSize)

  for (i in 1:nrow(allMenus)) {
    for (j in 1:n) {
      menu[j+(i-1)*n, allMenus[i, ]] <- 1
      choice[j+(i-1)*n, sort(allMenus[i, ])[rmultinom(1, 1, choiceProb) == 1]] <- 1
    }
  }

  return(list(menu=menu, choice=choice))
}
