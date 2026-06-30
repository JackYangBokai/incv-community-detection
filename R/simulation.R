#' Simulate a Stochastic Block Model network
#'
#' Generates an adjacency matrix from a planted-partition SBM with \code{k}
#' communities. The first community has size \code{n1}; the remaining
#' communities share the leftover nodes roughly equally.
#'
#' @param k Number of communities (default 2).
#' @param n Total number of nodes (default 1000).
#' @param n1 Size of the smallest community (default 100).
#' @param p Within-community connection probability (default 0.3).
#' @param q Between-community connection probability (default 0.1).
#' @return A list with:
#'   \item{membership}{Integer vector of community labels (length \code{n}).}
#'   \item{adjacency}{n x n binary symmetric adjacency matrix.}
#' @export
#' @examples
#' net <- community.sim(k = 3, n = 120, n1 = 30, p = 0.5, q = 0.1)
#' table(net$membership)
community.sim <- function(k = 2, n = 1000, n1 = 100, p = 0.3, q = 0.1) {
  size <- rep(0, k)
  size[1] <- n1
  i <- 2
  while (i < k) {
    size[i] <- floor((n - n1) / (k - 1))
    i <- i + 1
  }
  size[k] <- n - sum(size)

  mem <- sample(rep(1:k, size))

  conn <- diag(rep(0, n))
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (mem[i] == mem[j]) conn[i, j] <- rbinom(1, 1, p)
      else conn[i, j] <- rbinom(1, 1, q)
      conn[j, i] <- conn[i, j]
    }
  }
  list(membership = mem, adjacency = conn)
}

#' Simulate an SBM with distance-decaying block probabilities
#'
#' Generates an SBM where block probabilities decay exponentially with
#' the distance between community indices: \code{B[k1,k2] = rho * eta^min(|k1-k2|,3)}.
#'
#' @param n Total number of nodes.
#' @param n1 Size of the first (smallest) community.
#' @param eta Decay parameter for inter-community probabilities (default 0.3).
#' @param rho Scaling factor for the block probability matrix (default 0.1).
#' @param K Number of communities (default 3).
#' @return A list with:
#'   \item{adjacency}{n x n binary symmetric adjacency matrix.}
#'   \item{membership}{Integer vector of community labels.}
#'   \item{conn}{K x K block probability matrix.}
#' @export
community.sim.sbm <- function(n, n1, eta = 0.3, rho = 0.1, K = 3) {
  B0 <- matrix(0, nrow = K, ncol = K)
  for (k1 in 1:K) {
    for (k2 in 1:K) {
      B0[k1, k2] <- eta^(min(abs(k1 - k2), 3))
    }
  }
  B0 <- B0 * rho
  Member_mat <- matrix(0, n, K)

  size <- rep(0, K)
  size[1] <- n1
  i <- 2
  while (i < K) {
    size[i] <- floor((n - n1) / (K - 1))
    i <- i + 1
  }
  size[K] <- n - sum(size)

  c0 <- sample(rep(1:K, size))
  for (i in 1:n) {
    Member_mat[i, c0[i]] <- 1
  }
  P0 <- Member_mat %*% B0 %*% t(Member_mat)

  A <- matrix(0, nrow = n, ncol = n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      A[i, j] <- rbinom(n = 1, size = 1, prob = P0[i, j])
      A[j, i] <- A[i, j]
    }
  }
  list(adjacency = A, membership = c0, conn = B0)
}

#' Simulate an DCSBM with distance-decaying block probabilities
#'
#' Generates an DCSBM where block probabilities decay exponentially with
#' the distance between community indices: \code{B[k1,k2] = rho * eta^min(|k1-k2|,3)}.
#' The degree heterogeneity parameter theta is generated from Uniform(1,10) independently,
#' and normalized for each cluster.
#'
#' @param n Total number of nodes.
#' @param n1 Size of the first (smallest) community.
#' @param eta Decay parameter for inter-community probabilities (default 0.3).
#' @param rho Scaling factor for the block probability matrix (default 0.1).
#' @param K Number of communities (default 3).
#' @return A list with:
#'   \item{adjacency}{n x n binary symmetric adjacency matrix.}
#'   \item{membership}{Integer vector of community labels.}
#'   \item{conn}{K x K block probability matrix.}
#'   \item{theta}{Vector of degree-heterogeneity parameter.}
#' @export
community.sim.dcsbm <- function(n, n1, eta=0.3, rho = 0.1, K = 3) {
  
  B0 <- matrix(0, nrow = K, ncol = K)
  for (k1 in 1:K) {
    for (k2 in 1:K) {
      B0[k1, k2] <- eta^(min(abs(k1-k2),3))
    }
  }
  B0 <- B0 * rho
  Member_mat <- matrix(0,n,K)
  
  size <- rep(0, K)
  size[1] <- n1 # Set the first community be the one with the smallest size n1
  i <- 2
  while (i<K) {
    size[i] <- floor((n-n1)/(K-1)) # The rest community has size (n-n1)/(k-1)
    i <- i + 1
  }
  size[K] <- n - sum(size) # The last community size is determine by the substraction from n
  
  c0 <- sample(rep(1:K,size))
  for (i in 1:n) {
    Member_mat[i, c0[i]] <- 1
  }
  P0 <- Member_mat %*% B0 %*% t(Member_mat)
  
  theta <- runif(n,1,10)
  for (k in 1:K) {
    theta[c0 == k] <- theta[c0 == k]/sqrt(sum(theta[c0 == k]^2)) * sqrt(sum(c0 == k))
  }
  
  P0 <- outer(theta, theta) * P0
  P0[P0 > 1] <- 1
  
  A <- matrix(0, nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      A[i, j] <- rbinom(n = 1, size = 1, prob = P0[i, j])
      A[j, i] <- A[i, j]
    }
  }
  return(list(adjacency = A, membership = c0, conn = B0, theta = theta))
}

#' Generate a fast SBM or DCBM network (sparse)
#'
#' Uses parallel edge sampling for speed; returns a sparse adjacency matrix.
#'
#' @param n Number of nodes.
#' @param K Number of communities.
#' @param B K x K block probability matrix.
#' @param psi Degree heterogeneity vector (length \code{n}); set to all 1s for SBM.
#' @param PI Community prior probabilities (length \code{K}).
#' @param ncore Number of cores for parallel computation.
#' @return A list with:
#'   \item{A}{Sparse symmetric adjacency matrix.}
#'   \item{member}{Community membership vector.}
#'   \item{psi}{Scaled degree heterogeneity vector.}
#' @export
blockmodel.gen.fast <- function(n, K, B, psi = rep(1, n),
                                PI = rep(1 / K, K), ncore = 1) {
  on.exit(gc())

  g <- sample.int(K, n, TRUE, PI)

  psi.scale <- psi
  for (kk in 1:K)
    psi.scale[g == kk] <- psi[g == kk] / max(psi[g == kk])

  stor <- do.call("rbind",
    parallel::mclapply(1:(n - 1), function(i) {
      tmp <- which(rbinom(n - i, 1,
        B[g[i], g[(i + 1):n]] * psi.scale[i] * psi.scale[(i + 1):n]) == 1)
      if (length(tmp) == 0) return(NULL)
      else return(cbind(rep(i, length(tmp)), i + tmp))
    }, mc.cores = ncore))

  A <- Matrix::sparseMatrix(stor[, 1], stor[, 2],
                            dims = c(n, n), symmetric = TRUE)
  list(A = A, member = g, psi = psi.scale)
}

#' Generate a sparse RDPG network
#'
#' @param n Number of nodes.
#' @param d Latent dimension.
#' @param sparsity.multiplier Multiplier to control network density.
#' @param ncore Number of cores for parallel computation.
#' @return A list with:
#'   \item{A}{Sparse symmetric adjacency matrix.}
#'   \item{P}{n x n probability matrix.}
#' @export
sparse.RDPG.gen <- function(n, d, sparsity.multiplier = 1, ncore = 1) {
  X <- matrix(runif(n * d), nrow = n, ncol = d)
  P <- tcrossprod(X)
  P <- P * sparsity.multiplier / max(P)
  diag(P) <- 0

  stor <- do.call("rbind",
    parallel::mclapply(1:(n - 1), function(i) {
      tmp <- which(rbinom(n - i, 1, P[i, (i + 1):n]) == 1)
      if (length(tmp) == 0) return(NULL)
      else return(cbind(rep(i, length(tmp)), i + tmp))
    }, mc.cores = ncore))

  A <- Matrix::sparseMatrix(i = stor[, 1], j = stor[, 2],
                            dims = c(n, n), symmetric = TRUE)
  list(A = A, P = P)
}

#' Generate a latent space network
#'
#' @param n Number of nodes.
#' @param d Latent dimension.
#' @param alpha Intercept parameter controlling overall density.
#' @param sparsity Sparsity multiplier.
#' @param ncore Number of cores for parallel computation.
#' @return A list with:
#'   \item{A}{Sparse symmetric adjacency matrix.}
#'   \item{Z}{n x d latent position matrix.}
#' @export
latent.gen <- function(n, d, alpha = 1, sparsity = 1, ncore = 1) {
  Z <- matrix(runif(n * d), nrow = n, ncol = d)

  stor <- do.call("rbind",
    parallel::mclapply(1:(n - 1), function(i) {
      logodds <- alpha - sapply((i + 1):n,
        function(jj) sqrt(sum((Z[i, ] - Z[jj, ])^2)))
      pp <- sparsity * exp(logodds) / (1 + exp(logodds))
      tmp <- which(rbinom(n - i, 1, pp) == 1)
      if (length(tmp) == 0) return(NULL)
      else return(cbind(rep(i, length(tmp)), i + tmp))
    }, mc.cores = ncore))

  A <- as(Matrix::sparseMatrix(i = stor[, 1], j = stor[, 2],
                               dims = c(n, n), symmetric = TRUE),
          "dMatrix")
  list(A = A, Z = Z)
}
