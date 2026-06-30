#' INCVCommunityDetection: Inductive Node-Splitting Cross-Validation for Community Detection
#'
#' Implements Inductive Node-Splitting Cross-Validation (INCV) for selecting the
#' number of communities in stochastic block models, along with competing
#' methods (CROISSANT, ECV, NCV).
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{Incv.f.fold}}}{F-fold node-split cross-validation (INCV).}
#'   \item{\code{\link{Incv.random.split}}}{Random-split cross-validation (INCV).}
#'   \item{\code{\link{Incv.f.fold.dcsbm}}}{F-fold node-split cross-validation (INCV) for DCSBM.}
#'   \item{\code{\link{croissant.blockmodel}}}{CROISSANT for SBM/DCBM selection.}
#'   \item{\code{\link{ECV.for.blockmodel}}}{Edge cross-validation.}
#'   \item{\code{\link{NCV.for.blockmodel}}}{Node cross-validation.}
#' }
#'
#' @docType package
#' @name INCVCommunityDetection-package
#' @keywords internal
"_PACKAGE"

#' @import Matrix
#' @importFrom RSpectra svds
#' @importFrom ClusterR KMeans_rcpp
#' @importFrom Gmedian kGmedian
#' @importFrom irlba irlba partial_eigen
#' @importFrom parallel mclapply
#' @importFrom cluster pam
#' @importFrom Rfast auc
#' @importFrom data.table data.table
#' @importFrom IMIFA Procrustes
#' @importFrom stats kmeans rbinom runif setNames
#' @importFrom methods as
NULL
