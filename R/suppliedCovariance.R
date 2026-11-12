##' Pseudo-inverse square root of a covariance, with the rank cut on the correlation scale
##'
##' A wrapper for \code{XtX_pseudoinv_sqrt(mat, mat.is.XtX = TRUE)} that
##' differs from calling it directly in two ways, both about which singular
##' values count as zero (#146).
##'
##' First, the matrix is put on the correlation scale before the cut and the
##' result scaled back, so the rank does not depend on the units of the
##' variables: a count in the hundreds beside a proportion would otherwise
##' have a variance ratio that a relative cut reads as a null direction.
##'
##' Second, \code{XtX_pseudoinv_sqrt()} compares square roots of singular
##' values against its \code{tol}, so it is given the square root of
##' \code{tol}; singular values of the correlation matrix below \code{tol}
##' times the largest then count as zero. With its default the cut sits at
##' machine precision, where the rank of a singular covariance, such as one
##' over a factor's full set of indicator columns, varies with the BLAS.
##'
##' The result \code{V} satisfies \code{V \%*\% t(V)} = a generalized
##' inverse of \code{mat}, which gives the same quadratic form as the
##' Moore-Penrose inverse for any vector in the column space of \code{mat}.
##' Variables with zero variance get zero rows.
##'
##' @param mat A symmetric positive semi-definite matrix.
##' @param tol Relative tolerance on the singular values of the correlation
##'   matrix.
##' @return As \code{XtX_pseudoinv_sqrt()}: a matrix with one row per column
##'   of \code{mat} and its rank as the \code{r} attribute.
##' @keywords internal
.pseudoinvSqrtOfCovariance <- function(mat, tol = .Machine$double.eps^0.5) {
  if (!nrow(mat)) return(XtX_pseudoinv_sqrt(mat, mat.is.XtX = TRUE))
  sds <- sqrt(pmax(diag(mat), 0))
  ok <- sds > 0
  inv_sd <- ifelse(ok, 1 / sds, 0)
  ans <- XtX_pseudoinv_sqrt(mat * tcrossprod(inv_sd), mat.is.XtX = TRUE, tol = sqrt(tol))
  r <- attr(ans, "r")
  ans <- ans * inv_sd
  attr(ans, "r") <- r
  ans
}

##' Arrange supplied covariances by stratification
##'
##' @param covariances \code{NULL}, or a named list of k by k by strata
##'   arrays (or k by k matrices), one per stratification to be tested
##'   against supplied covariances.
##' @param stratifications The names of the stratifications being tested.
##' @return A list parallel to \code{stratifications}, \code{NULL} where
##'   none was supplied.
##' @keywords internal
.matchCovariances <- function(covariances, stratifications) {
  out <- setNames(vector("list", length(stratifications)), stratifications)
  if (is.null(covariances) || (is.list(covariances) && !length(covariances))) return(out)
  if (!is.list(covariances) || is.null(names(covariances)) || any(names(covariances) == "") ||
      anyDuplicated(names(covariances))) {
    stop("covariances must be a list with distinct names, those of the stratifications it applies to (",
         paste(stratifications, collapse = ", "), ")")
  }
  unknown <- setdiff(names(covariances), stratifications)
  if (length(unknown)) {
    stop("covariances given for stratification(s) not being tested: ",
         paste(unknown, collapse = ", "), ". Those tested are: ",
         paste(stratifications, collapse = ", "))
  }
  for (s in names(covariances)) {
    cv <- covariances[[s]]
    d <- dim(cv)
    if (!is.numeric(cv) || !(length(d) %in% 2:3) || d[1] != d[2]) {
      stop("covariances[[\"", s, "\"]] must be a k by k by strata array or a k by k matrix")
    }
    dn <- dimnames(cv)
    if (is.null(dn) || is.null(dn[[1]]) || !identical(dn[[1]], dn[[2]]) || anyDuplicated(dn[[1]])) {
      stop("covariances[[\"", s, "\"]] needs matching, distinct names on its first two dimensions")
    }
    if (length(d) == 3 && (is.null(dn[[3]]) || anyDuplicated(dn[[3]]))) {
      stop("covariances[[\"", s, "\"]] needs its third dimension named for the strata")
    }
    if (any(is.infinite(cv))) stop("covariances[[\"", s, "\"]] has infinite entries")
    if (!isTRUE(all.equal(cv, aperm(cv, c(2, 1, 3)[seq_along(d)]), check.attributes = FALSE))) {
      stop("covariances[[\"", s, "\"]] must be symmetric in its first two dimensions")
    }
    out[[s]] <- cv
  }
  out
}

##' Hansen and Bowers (2008) inferentials against supplied within-stratum covariances
##'
##' As \code{\link{HB08}}, except that the covariance of the covariates
##' within each stratum is supplied instead of being calculated from the
##' sample. The null covariance of the adjusted differences is the sum over
##' strata of the stratum's supplied covariance times the squared stratum
##' weight ratio over \eqn{1/n_{tb} + 1/n_{cb}}, the stratum's treated and
##' control counts; strata with one unit, or lacking either group, contribute
##' nothing. Supplying each stratum's own sample covariance reproduces
##' \code{HB08} when no covariate value is missing and unit weights are 1.
##' Where values were imputed, a stratum's contribution to an entry is
##' scaled by the share of its units observed on both variables.
##'
##' @param alignedcovs A \code{CovsAlignedToADesign} object.
##' @param covariance A k by k by strata array, named for the covariate
##'   columns and the strata, or a k by k matrix used for every stratum.
##' @param tol Relative tolerance for the rank of the null covariance.
##' @return A list as \code{\link{HB08}} returns.
##' @keywords internal
HB08_supplied <- function(alignedcovs, covariance, tol = .Machine$double.eps^0.5) {
  zz <- as.numeric(alignedcovs@Z)
  S <- alignedcovs@StrataMatrix
  Covs <- alignedcovs@Covariates
  NMcov <- alignedcovs@NotMissing
  n_ <- nrow(Covs)
  s_ <- ncol(S)
  wtr <- alignedcovs@StrataWeightRatio
  stopifnot(length(wtr) == n_, identical(dim(NMcov), dim(Covs)))
  vars <- colnames(Covs)

  ## which columns the supplied covariances cover
  cvars <- dimnames(covariance)[[1]]
  covariate_cols <- vars[!startsWith(vars, "(")]
  lacking <- setdiff(covariate_cols, cvars)
  if (length(lacking)) {
    stop("covariate column(s) not in the supplied covariances: ",
         paste(lacking, collapse = ", "),
         ". They must be over the columns expandCovariates() produces.")
  }
  tested <- intersect(cvars, vars)

  ## the stratum each unit is in, and each stratum's layer of the array
  stratum <- sparseToVec(S %*% seq_len(s_))
  if (length(dim(covariance)) == 3) {
    levs <- alignedcovs@StrataLevels
    layer <- match(levs, dimnames(covariance)[[3]])
  } else {
    layer <- NULL
  }

  ## stratum counts and scalars
  n <- sparseToVec(t(S) %*% rep(1, n_))
  n1 <- sparseToVec(t(S) %*% zz)
  n0 <- n - n1
  n1_over_n <- sparseToVec(S %*% ifelse(n > 0, n1 / n, 0))
  wtr_b <- ifelse(n > 0, sparseToVec(t(S) %*% wtr) / n, 0)
  usable <- n > 1 & n0 > 0 & n1 > 0
  c_b <- ifelse(usable, wtr_b^2 / (1 / n0 + 1 / n1), 0)
  used <- which(c_b > 0)

  ## only the strata that contribute need a layer, and a usable one
  if (!is.null(layer)) {
    if (anyNA(layer[used])) {
      absent <- levs[used][is.na(layer[used])]
      stop("no supplied covariance for stratum level(s): ",
           paste(utils::head(absent), collapse = ", "), if (length(absent) > 6) ", ..." else "")
    }
    sub <- covariance[tested, tested, layer[used], drop = FALSE]
  } else {
    sub <- covariance[tested, tested, drop = FALSE]
  }
  if (anyNA(sub)) stop("the supplied covariances have missing entries for strata that are tested")
  negative <- if (is.null(layer)) {
    any(diag(sub) < 0)
  } else {
    any(apply(sub, 3, function(m) any(diag(matrix(m, nrow(sub))) < 0)))
  }
  if (negative) {
    stop("the supplied covariances have a negative variance on a diagonal")
  }

  ## adjusted differences, as in HB08
  x_tilde <- Covs * wtr
  ssn <- drop(crossprod(zz - n1_over_n, x_tilde))
  names(ssn) <- vars

  ## null covariance: sum over strata of c_b times the stratum's covariance
  k <- length(tested)
  block <- matrix(0, k, k, dimnames = list(tested, tested))
  sigma_of <- if (is.null(layer)) {
    function(j) sub
  } else {
    function(j) matrix(sub[, , j], k, k)
  }
  if (length(used)) {
    block[] <- if (is.null(layer)) sum(c_b) * sub else matrix(sub, k * k, length(used)) %*% c_b[used]
    ## strata where a tested value was imputed: their entries are scaled by
    ## the share of units observed on both variables
    NMt <- NMcov[, tested, drop = FALSE]
    touched <- intersect(unique(stratum[rowSums(NMt != 1) > 0]), used)
    if (length(touched)) {
      rows_of <- split(seq_len(n_), stratum)
      for (b in touched) {
        share <- crossprod(NMt[rows_of[[as.character(b)]], , drop = FALSE]) / n[b]
        block <- block + c_b[b] * (share - 1) * sigma_of(match(b, used))
      }
    }
  }
  tcov <- matrix(0, length(vars), length(vars), dimnames = list(vars, vars))
  tcov[tested, tested] <- block
  ssvar <- diag(tcov)
  zero_variance <- ssvar <= .Machine$double.eps
  zstat <- ifelse(zero_variance, NA_real_, ssn / sqrt(ssvar))
  p <- 2 * pnorm(abs(zstat), lower.tail = FALSE)

  keep <- tested[!zero_variance[tested]]
  inv <- .pseudoinvSqrtOfCovariance(block[keep, keep, drop = FALSE], tol)
  check_for_degenerate(inv, n_, s_)
  DF <- attr(inv, "r")
  csq <- if (length(keep) && DF > 0) sum(drop(crossprod(ssn[keep], inv))^2) else NA_real_

  list(z = zstat, p = p, Msq = csq, DF = DF,
       adj.diff.of.totals = ssn, tcov = tcov)
}
