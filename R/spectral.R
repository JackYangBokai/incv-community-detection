#' Spectral clustering for a Stochastic Block Model
#'
#' Performs spectral clustering on an adjacency matrix by computing the top
#' \code{k} singular vectors and applying k-means++ via \code{ClusterR::KMeans_rcpp}.
#'
#' @param A Symmetric adjacency matrix (n x n, binary or weighted).
#' @param k Number of clusters (default 2).
#' @return A list with component \code{cluster}, an integer vector of length n.
#' @export
#' @examples
#' net <- community.sim(k = 3, n = 90, n1 = 30, p = 0.5, q = 0.1)
#' cl <- SBM.spectral.clustering(net$adjacency, k = 3)
#' table(cl$cluster, net$membership)
SBM.spectral.clustering <- function(A, k = 2) {
  svd.k <- RSpectra::svds(A, k = k)
  U_k <- svd.k$u

  kpp_result <- ClusterR::KMeans_rcpp(U_k, clusters = k,
    num_init = 25, max_iters = 100,
    initializer = "kmeans++", fuzzy = FALSE)
  list(cluster = kpp_result$clusters)
}

#' Estimate SBM connection probabilities and negative log-likelihood
#'
#' Given a clustering and adjacency matrix, estimates the block probability
#' matrix and computes the negative log-likelihood.
#'
#' @param cluster Integer vector of cluster labels.
#' @param k Number of clusters.
#' @param A Adjacency matrix corresponding to the nodes in \code{cluster}.
#' @param restricted Logical. If \code{TRUE}, uses a restricted SBM with a
#'   single within-community probability \code{p} and a single
#'   between-community probability \code{q}. If \code{FALSE}, estimates a
#'   separate probability for each pair of communities.
#' @return A list with:
#'   \item{p.matrix}{k x k estimated block probability matrix.}
#'   \item{negloglike}{Negative log-likelihood of the model.}
#' @export
SBM.prob <- function(cluster, k, A, restricted = TRUE) {
  p.matrix <- matrix(0, nrow = k, ncol = k)
  negloglike <- 0

  edge.vector <- c(A)[c(upper.tri(A))]
  one <- edge.index.map(which(edge.vector == 1))
  zero <- edge.index.map(which(edge.vector == 0))

  if (restricted) {
    within.connect <- sum(cluster[one$x] == cluster[one$y])
    within.disconnect <- sum(cluster[zero$x] == cluster[zero$y])
    between.connect <- length(one$x) - within.connect
    between.disconnect <- length(zero$x) - within.disconnect

    within.total <- within.connect + within.disconnect
    between.total <- between.connect + between.disconnect
    if (within.total == 0) p <- 0
    else p <- within.connect / within.total
    if (between.total == 0) q <- 0
    else q <- between.connect / between.total
    diag(p.matrix) <- p
    p.matrix[(lower.tri(p.matrix)) | (upper.tri(p.matrix))] <- q

    negloglike <- neglog(within.connect, p) + neglog(within.disconnect, 1 - p) +
      neglog(between.connect, q) + neglog(between.disconnect, 1 - q)
  } else {
    for (i in 1:k) {
      for (j in i:k) {
        if (i == j) {
          connect <- sum((cluster[one$x] == i) & (cluster[one$y] == i))
          disconnect <- sum((cluster[zero$x] == i) & (cluster[zero$y] == i))
          total <- connect + disconnect
          if (total == 0) p <- 0
          else p <- connect / total
          p.matrix[i, i] <- p
          negloglike <- negloglike + neglog(connect, p) + neglog(disconnect, 1 - p)
        } else {
          connect <- sum((cluster[one$x] == i) & (cluster[one$y] == j)) +
            sum((cluster[one$x] == j) & (cluster[one$y] == i))
          disconnect <- sum((cluster[zero$x] == i) & (cluster[zero$y] == j)) +
            sum((cluster[zero$x] == j) & (cluster[zero$y] == i))
          total <- connect + disconnect
          if (total == 0) q <- 0
          else q <- connect / total
          p.matrix[i, j] <- q
          p.matrix[j, i] <- q
          negloglike <- negloglike + neglog(connect, q) + neglog(disconnect, 1 - q)
        }
      }
    }
  }
  list(p.matrix = p.matrix, negloglike = negloglike)
}

#' Spectral clustering for a Degree-corrected Stochastic Block Model
#'
#' Performs spectral clustering on an adjacency matrix by computing the top
#' \code{k} singular vectors and applying kGmedian via \code{Gmedian::kGmedian}.
#'
#' @param A Symmetric adjacency matrix (n x n, binary or weighted).
#' @param k Number of clusters (default 2).
#' @return A list with:
#'   \item{label_hat}{An integer vector of length n.}
#'   \item{theta_hat}{A vector with length n, estimator of normalized degree parameter.}
#' @export
DCSBM.spectral.clustering <- function(A, k=2){
  svd.k <- svds(A, k = k)
  U_k <- svd.k$u  
  D_k <- diag(svd.k$d)      
  V_k <- svd.k$v      
  
  theta_hat <- sqrt(rowSums(U_k^2))
  normalized_U <- U_k / theta_hat  
  
  normalized_U <- as.matrix(normalized_U)
  
  if (k == 1) {
    label_hat <- rep(1, nrow(normalized_U))
  } else {
    normalized_U <- as.matrix(normalized_U)
    
    kmedian_result <- kGmedian(normalized_U, ncenters = k, nstart = 25, iter.max = 100)
    label_hat <- kmedian_result$cluster
  }
  
  return(list(label_hat = label_hat, theta_hat = theta_hat))
}

#' Estimate DESBM connection probabilities
#'
#' Given a clustering and adjacency matrix, estimates the block probability
#' matrix of the DCSBM model.
#'
#' @param cluster Integer vector of cluster labels.
#' @param k Number of clusters.
#' @param A Adjacency matrix corresponding to the nodes in \code{cluster}.
#' @param eps Gap that ensures the probability is less than 1.
#' @return A list with \code{B.matrix}, a k x k estimated block probability matrix.
#' @export
DCSBM.prob <- function(cluster, k, A, theta_hat, eps = 1e-8) {
  n <- nrow(A)
  
  if (ncol(A) != n) {
    stop("A must be a square adjacency matrix.")
  }
  if (length(cluster) != n) {
    stop("length(cluster) must equal nrow(A).")
  }
  if (length(theta_hat) != n) {
    stop("length(theta_hat) must equal nrow(A).")
  }
  
  B.matrix <- matrix(0, nrow = k, ncol = k)
  
  idx <- which(upper.tri(A), arr.ind = TRUE)
  ii <- idx[, 1]
  jj <- idx[, 2]
  
  A.vec <- A[cbind(ii, jj)]
  theta.prod <- theta_hat[ii] * theta_hat[jj]
  
  for (a in 1:k) {
    for (b in a:k) {
      
      if (a == b) {
        sel <- (cluster[ii] == a) & (cluster[jj] == a)
      } else {
        sel <- ((cluster[ii] == a) & (cluster[jj] == b)) |
          ((cluster[ii] == b) & (cluster[jj] == a))
      }
      
      numerator <- sum(A.vec[sel])
      denominator <- sum(theta.prod[sel])
      
      if (denominator <= 0 || sum(sel) == 0) {
        B.hat <- 0
      } else {
        B.hat <- numerator / denominator
      }
      
      if (sum(sel) > 0) {
        max.theta.prod <- max(theta.prod[sel])
        if (max.theta.prod > 0) {
          B.upper <- (1 - eps) / max.theta.prod
          B.hat <- min(B.hat, B.upper)
        }
      }
      
      B.matrix[a, b] <- B.hat
      B.matrix[b, a] <- B.hat
    }
  }
  
  return(list(B.matrix = B.matrix))
}
