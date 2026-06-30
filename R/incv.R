#' Inductive Node-Splitting Cross-Validation with f-fold splitting
#'
#' Performs f-fold node-split cross-validation to select the number of
#' communities \code{k} in a Stochastic Block Model. Nodes are randomly
#' partitioned into \code{f} folds; for each fold the held-out nodes are
#' assigned to communities via the training-set spectral clustering, and the
#' held-out negative log-likelihood and MSE are computed.
#'
#' @param A Symmetric adjacency matrix (n x n, binary).
#' @param k.vec Integer vector of candidate community numbers (default 2:6).
#' @param restricted Logical. If \code{TRUE}, use a restricted SBM (one within-
#'   and one between-community probability). Default \code{TRUE}.
#' @param f Number of folds (default 10).
#' @param method Inference method for assigning held-out nodes:
#'   \code{"affinity"} (default) assigns by maximum average connection, or
#'   \code{"loss"} assigns by minimum negative log-likelihood.
#' @param p.est.type How to re-estimate the probability matrix for evaluation:
#'   1 = training only, 2 = training + testing, 3 = testing only (default 3).
#' @return A list with:
#'   \item{k.loss}{Selected \code{k} that minimises CV negative log-likelihood.}
#'   \item{k.mse}{Selected \code{k} that minimises CV MSE.}
#'   \item{cv.loss}{Average CV negative log-likelihood for each \code{k}.}
#'   \item{cv.mse}{Average CV MSE for each \code{k}.}
#' @export
#' @examples
#' set.seed(42)
#' net <- community.sim(k = 3, n = 150, n1 = 50, p = 0.5, q = 0.1)
#' result <- nscv.f.fold(net$adjacency, k.vec = 2:5, f = 5)
#' result$k.loss
Incv.f.fold <- function(A, k.vec = 2:6, restricted = TRUE, f = 10,
                        method = "affinity", p.est.type = 3) {
  n <- ncol(A)
  cv.loss <- rep(0, length(k.vec))
  loss.f <- matrix(0, nrow = f, ncol = length(k.vec))
  mse.f  <- matrix(0, nrow = f, ncol = length(k.vec))
  cv.loss <- rep(0, length(k.vec))
  node <- 1:n

  f.group <- rep(1:f, floor(n / f))
  if ((n %% f) != 0) f.group <- c(f.group, 1:(n %% f))
  f.group <- sample(f.group)

  for (i in 1:f) {
    training <- node[f.group != i]
    testing <- node[f.group == i]

    A11 <- A[training, training]
    A12 <- A[training, testing]
    A22 <- A[testing, testing]

    for (j in seq_along(k.vec)) {
      k <- k.vec[j]

      tr.cluster <- SBM.spectral.clustering(A11, k)$cluster
      p.matrix <- SBM.prob(tr.cluster, k, A11, restricted)$p.matrix

      te.cluster <- rep(0, length(testing))
      affinity <- rep(0, k)
      ni <- rep(0, k)
      mi <- rep(0, k)
      for (t in 1:k) {
        ni[t] <- length(tr.cluster[tr.cluster == t])
      }

      for (s in 1:length(testing)) {
        loss <- rep(0, k)
        for (t in 1:k) {
          mi[t] <- sum(A12[, s][tr.cluster == t])
          affinity[t] <- mi[t] / ni[t]
        }
        for (t in 1:k) {
          loss_u <- 0
          for (u in 1:k) {
            loss_u <- loss_u +
              neglog(mi[u], p.matrix[t, u]) +
              neglog(ni[u] - mi[u], 1 - p.matrix[t, u])
          }
          loss[t] <- loss_u
        }

        if (method == "affinity") {
          group <- which.max(affinity)
        } else if (method == "loss") {
          group <- which.min(loss)
        }
        te.cluster[s] <- group
      }

      if (p.est.type == 1) {
        cluster <- tr.cluster
        AA <- A11
      } else if (p.est.type == 2) {
        cluster <- c(tr.cluster, te.cluster)
        AA <- cbind(A11, A12)
      } else if (p.est.type == 3) {
        cluster <- te.cluster
        AA <- A22
      }

      p.matrix <- SBM.prob(cluster, k, AA, restricted)$p.matrix

      te.negloglike <- 0
      te.mse <- 0

      te.edge.vector <- c(A22)[c(upper.tri(A22))]
      te.one <- edge.index.map(which(te.edge.vector == 1))
      te.zero <- edge.index.map(which(te.edge.vector == 0))

      for (s in 1:k) {
        for (t in s:k) {
          if (s == t) {
            connect <- sum((te.cluster[te.one$x] == s) &
                             (te.cluster[te.one$y] == s))
            disconnect <- sum((te.cluster[te.zero$x] == s) &
                                (te.cluster[te.zero$y] == s))
            p <- p.matrix[s, s]
            te.negloglike <- te.negloglike +
              neglog(connect, p) + neglog(disconnect, 1 - p)
            te.mse <- te.mse + connect * (1 - p)^2 + disconnect * p^2
          } else {
            connect <- sum((te.cluster[te.one$x] == s) &
                             (te.cluster[te.one$y] == t)) +
              sum((te.cluster[te.one$x] == t) &
                    (te.cluster[te.one$y] == s))
            disconnect <- sum((te.cluster[te.zero$x] == s) &
                                (te.cluster[te.zero$y] == t)) +
              sum((te.cluster[te.zero$x] == t) &
                    (te.cluster[te.zero$y] == s))
            q <- p.matrix[s, t]
            te.negloglike <- te.negloglike +
              neglog(connect, q) + neglog(disconnect, 1 - q)
            te.mse <- te.mse + connect * (1 - q)^2 + disconnect * q^2
          }
        }
      }
      loss.f[i, j] <- te.negloglike
      mse.f[i, j] <- te.mse
    }
    cv.loss <- colMeans(loss.f)
    cv.mse <- colMeans(mse.f)
    k.loss <- k.vec[which.min(cv.loss)]
    k.mse  <- k.vec[which.min(cv.mse)]
  }
  list(k.loss = k.loss, k.mse = k.mse,
       cv.loss = cv.loss, cv.mse = cv.mse)
}

#' Inductive Node-Splitting Cross-Validation with random node splits
#'
#' Performs repeated random-split node-level cross-validation. At each
#' iteration a random fraction \code{split} of nodes is used for training.
#'
#' @inheritParams Incv.f.fold
#' @param split Fraction of nodes used for training (default 0.66).
#' @param ite Number of random split iterations (default 100).
#' @return A list with:
#'   \item{k.loss}{Selected \code{k} that minimises CV negative log-likelihood.}
#'   \item{k.mse}{Selected \code{k} that minimises CV MSE.}
#'   \item{cv.loss}{Average CV negative log-likelihood for each \code{k}.}
#'   \item{cv.mse}{Average CV MSE for each \code{k}.}
#' @export
Incv.random.split <- function(A, k.vec = 2:6, restricted = TRUE,
                              split = 0.66, ite = 100,
                              method = "affinity", p.est.type = 3) {
  n <- ncol(A)
  cv.loss <- rep(0, length(k.vec))
  loss.f <- matrix(0, nrow = ite, ncol = length(k.vec))
  mse.f  <- matrix(0, nrow = ite, ncol = length(k.vec))
  cv.loss <- rep(0, length(k.vec))
  node <- 1:n

  for (i in 1:ite) {
    training <- sort(sample(1:n, ceiling(n * split)))
    testing <- sort(setdiff(1:n, training))

    A11 <- A[training, training]
    A12 <- A[training, testing]
    A22 <- A[testing, testing]

    for (j in seq_along(k.vec)) {
      k <- k.vec[j]

      tr.cluster <- SBM.spectral.clustering(A11, k)$cluster
      p.matrix <- SBM.prob(tr.cluster, k, A11, restricted)$p.matrix

      te.cluster <- rep(0, length(testing))
      affinity <- rep(0, k)
      ni <- rep(0, k)
      mi <- rep(0, k)
      for (t in 1:k) {
        ni[t] <- length(tr.cluster[tr.cluster == t])
      }

      for (s in 1:length(testing)) {
        loss <- rep(0, k)
        for (t in 1:k) {
          mi[t] <- sum(A12[, s][tr.cluster == t])
          affinity[t] <- mi[t] / ni[t]
        }
        for (t in 1:k) {
          loss_u <- 0
          for (u in 1:k) {
            loss_u <- loss_u +
              neglog(mi[u], p.matrix[t, u]) +
              neglog(ni[u] - mi[u], 1 - p.matrix[t, u])
          }
          loss[t] <- loss_u
        }

        if (method == "affinity") {
          group <- which.max(affinity)
        } else if (method == "loss") {
          group <- which.min(loss)
        }
        te.cluster[s] <- group
      }

      if (p.est.type == 1) {
        cluster <- tr.cluster
        AA <- A11
      } else if (p.est.type == 2) {
        cluster <- c(tr.cluster, te.cluster)
        AA <- cbind(A11, A12)
      } else if (p.est.type == 3) {
        cluster <- te.cluster
        AA <- A22
      }
      p.matrix <- SBM.prob(cluster, k, AA, restricted)$p.matrix

      te.negloglike <- 0
      te.mse <- 0

      te.edge.vector <- c(A22)[c(upper.tri(A22))]
      te.one <- edge.index.map(which(te.edge.vector == 1))
      te.zero <- edge.index.map(which(te.edge.vector == 0))

      for (s in 1:k) {
        for (t in s:k) {
          if (s == t) {
            connect <- sum((te.cluster[te.one$x] == s) &
                             (te.cluster[te.one$y] == s))
            disconnect <- sum((te.cluster[te.zero$x] == s) &
                                (te.cluster[te.zero$y] == s))
            p <- p.matrix[s, s]
            te.negloglike <- te.negloglike +
              neglog(connect, p) + neglog(disconnect, 1 - p)
            te.mse <- te.mse + connect * (1 - p)^2 + disconnect * p^2
          } else {
            connect <- sum((te.cluster[te.one$x] == s) &
                             (te.cluster[te.one$y] == t)) +
              sum((te.cluster[te.one$x] == t) &
                    (te.cluster[te.one$y] == s))
            disconnect <- sum((te.cluster[te.zero$x] == s) &
                                (te.cluster[te.zero$y] == t)) +
              sum((te.cluster[te.zero$x] == t) &
                    (te.cluster[te.zero$y] == s))
            q <- p.matrix[s, t]
            te.negloglike <- te.negloglike +
              neglog(connect, q) + neglog(disconnect, 1 - q)
            te.mse <- te.mse + connect * (1 - q)^2 + disconnect * q^2
          }
        }
      }
      loss.f[i, j] <- te.negloglike
      mse.f[i, j] <- te.mse
    }
    cv.loss <- colMeans(loss.f)
    cv.mse <- colMeans(mse.f)
    k.loss <- k.vec[which.min(cv.loss)]
    k.mse  <- k.vec[which.min(cv.mse)]
  }
  list(k.loss = k.loss, k.mse = k.mse,
       cv.loss = cv.loss, cv.mse = cv.mse)
}

#' Inductive Node-Splitting Cross-Validation with f-fold splitting for DCSBM
#'
#' Performs f-fold node-split cross-validation to select the number of
#' communities \code{k} in a Degree-corrected Stochastic Block Model. Nodes 
#' are randomly partitioned into \code{f} folds; for each fold the held-out 
#' nodes are assigned to communities via the training-set spectral clustering, 
#' and the held-out negative log-likelihood and MSE are computed.
#'
#' @param A Symmetric adjacency matrix (n x n, binary).
#' @param k.vec Integer vector of candidate community numbers (default 2:6).
#' @param f Number of folds (default 10).
#' @return A list with:
#'   \item{k.loss}{Selected \code{k} that minimises CV negative log-likelihood.}
#'   \item{k.mse}{Selected \code{k} that minimises CV MSE.}
#'   \item{cv.loss}{Average CV negative log-likelihood for each \code{k}.}
#'   \item{cv.mse}{Average CV MSE for each \code{k}.}
#' @export
Incv.f.fold.dcsbm <- function(A, k.vec = 2:6, f=10) {
  n <- ncol(A)
  cv.loss <- rep(0, length(k.vec))
  loss.f <- matrix(0, nrow=f, ncol=length(k.vec))
  mse.f  <- matrix(0, nrow=f, ncol=length(k.vec))
  node <- 1:n
  
  f.group <- rep(1:f, floor(n/f))
  if ((n%%f) != 0) f.group <- c(f.group, 1:(n%%f))
  f.group <- sample(f.group) # Get a random assignment of nodes to the f folds.
  
  for(i in 1:f) {
    
    # Split training and testing set
    training <- node[f.group != i]
    testing <- node[f.group == i]  
    
    # Split the adjacency matrix
    A11 <- A[training, training]
    A12 <- A[training, testing]
    A22 <- A[testing, testing]
    
    for (j in (1:length(k.vec))) {
      
      k <- k.vec[j]
      spec.res <- DCSBM.spectral.clustering(A11, k)
      tr.cluster <- spec.res$label_hat
      tr.theta <- spec.res$theta_hat
      B.matrix <- DCSBM.prob(tr.cluster, k, A11, tr.theta)$B.matrix
      
      
      te.cluster <- rep(0, length(testing))
      te.theta <- rep(0,length(testing))
      affinity <- rep(0,k)
      ni <- rep(0,k)
      mi <- rep(0,k)
      for (t in 1:k) {
        ni[t] <- sum(tr.theta[tr.cluster==t])
      }
      
      for (s in 1:length(testing)) {
        loss <- rep(0,k)
        for (t in 1:k) {
          mi[t] <-  sum(A12[,s][tr.cluster==t])
          affinity[t] <- mi[t]/ni[t]
        }
        group <- which.max(affinity)
        te.cluster[s] <- group
        
        mean.thetaP <- 0
        for (t in 1:k){
          mean.thetaP <- mean.thetaP + ni[t]*B.matrix[t,group]
        }
        if (mean.thetaP <= 1e-8) {
          te.theta[s] <- 0
        } else {
          te.theta[s] <- sum(mi) / mean.thetaP
        }
      }
      
      idx22 <- which(upper.tri(A22), arr.ind = TRUE)
      ii <- idx22[, 1]
      jj <- idx22[, 2]
      
      A22.vec <- A22[cbind(ii, jj)]
      
      prob.vec <- te.theta[ii] *
        te.theta[jj] *
        B.matrix[cbind(te.cluster[ii], te.cluster[jj])]
      
      prob.vec <- pmin(pmax(prob.vec, 1e-8), 1 - 1e-8)
      
      te.negloglike <- -sum(
        A22.vec * log(prob.vec) +
          (1 - A22.vec) * log(1 - prob.vec)
      )
      
      te.mse <- sum((A22.vec - prob.vec)^2)
      loss.f[i, j] <- te.negloglike
      mse.f[i, j] <- te.mse
    }
  }
  
  cv.loss <- colMeans(loss.f)
  cv.mse <- colMeans(mse.f)
  k.loss <- k.vec[which.min(cv.loss)]
  k.mse  <- k.vec[which.min(cv.mse)]
  
  return(list(
    k.loss = k.loss,
    k.mse = k.mse,
    cv.loss = cv.loss,
    cv.mse = cv.mse
  ))
}
