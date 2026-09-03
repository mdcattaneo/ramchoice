.ram_validate_cluster <- function(cluster, n_observations) {
  if (is.null(cluster)) {
    return(NULL)
  }
  if (length(cluster) != n_observations || anyNA(cluster)) {
    stop(
      "'cluster' must contain one nonmissing identifier per observation.",
      call. = FALSE
    )
  }
  labels <- unique(as.character(cluster))
  if (length(labels) < 2L) {
    stop("Cluster-robust inference requires at least two clusters.", call. = FALSE)
  }
  list(
    id = match(as.character(cluster), labels),
    labels = labels,
    n = length(labels)
  )
}

.ram_cluster_covariance <- function(scores) {
  n_cluster <- nrow(scores)
  crossprod(scores) * n_cluster / (n_cluster - 1L)
}

.ram_cluster_multiplier_draws <- function(scores, n_draws) {
  n_cluster <- nrow(scores)
  weights <- matrix(
    stats::rnorm(n_draws * n_cluster),
    nrow = n_draws,
    ncol = n_cluster
  )
  weights %*% scores * sqrt(n_cluster / (n_cluster - 1L))
}

.ram_row_maximum <- function(value) {
  if (!ncol(value)) {
    return(numeric(nrow(value)))
  }
  apply(value, 1L, max)
}
