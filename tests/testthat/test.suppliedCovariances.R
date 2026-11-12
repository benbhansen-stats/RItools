################################################################################
# balanceTest() generic: formula and design methods; supplied covariances
################################################################################
library("testthat")
context("balanceTest methods and supplied covariances")

make_data <- function(n = 240, seed = 20261029) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- 0.6 * x1 + rnorm(n)
  x3 <- rnorm(n)
  g <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  s <- factor(rep(1:4, length.out = n))
  z <- rbinom(n, 1, 0.5)
  data.frame(z = z, x1 = x1, x2 = x2, x3 = x3, g = g, s = s,
             stringsAsFactors = FALSE)
}

dat <- make_data()
Sigma <- cov(expandCovariates(~ x1 + x2 + g, dat))          # singular: ga + gb + gc = 1
SigmaN <- cov(expandCovariates(~ x1 + x2 + x3, dat))        # full rank

test_that("expandCovariates gives the columns balanceTest tests", {
  x <- expandCovariates(~ x1 + x2 + g, dat)
  expect_equal(colnames(x), c("x1", "x2", "ga", "gb", "gc"))
  expect_equal(unname(x[, "x1"]), dat$x1)
  expect_equal(unname(x[, "ga"]), as.numeric(dat$g == "a"))
  bt <- balanceTest(z ~ x1 + x2 + g, data = dat)
  vars <- dimnames(bt$results)[["vars"]]
  expect_equal(colnames(x), vars[!startsWith(vars, "(")])
  expect_equal(attr(x, "originals"), c(1L, 2L, 3L, 3L, 3L))

  ## response, strata() and cluster() terms are ignored
  expect_equal(expandCovariates(z ~ x1 + g + strata(s), dat),
               expandCovariates(~ x1 + g, dat))
  expect_equal(expandCovariates(~ x1 + cluster(s) + g, dat),
               expandCovariates(~ x1 + g, dat))

  ## character columns become factors; NAs are kept
  d2 <- dat
  d2$g <- as.character(d2$g)
  d2$x1[3] <- NA
  x2 <- expandCovariates(~ x1 + g, d2)
  expect_equal(colnames(x2), c("x1", "ga", "gb", "gc"))
  expect_true(is.na(x2[3, "x1"]))
})

test_that("balanceTest is a generic and the formula interface is unchanged", {
  expect_true(is.function(balanceTest.formula))
  bt <- balanceTest(z ~ x1 + x2 + g + strata(s), data = dat)
  expect_s3_class(bt, "balancetest")
  expect_equal(bt, balanceTest(z ~ x1 + x2 + g + strata(s), dat))
  expect_equal(bt, balanceTest(data = dat, z ~ x1 + x2 + g + strata(s)))
  expect_equal(bt, balanceTest(fmla = z ~ x1 + x2 + g + strata(s), data = dat))
  ## an argument neither method knows is an error, not silently dropped
  expect_error(balanceTest(z ~ x1, data = dat, p.adjust.methods = "none"), "p.adjust.methods")
  ## the design method gives the formula method's answer, less the formula
  d <- dat; d$`(weights)` <- 1
  des <- RItools:::makeDesigns(z ~ x1 + x2 + g + strata(s), d)
  bd <- balanceTest(des)
  expect_s3_class(bd, "balancetest")
  expect_equal(bd$overall, bt$overall)
  expect_equal(bd$results, bt$results, check.attributes = FALSE)
  expect_null(attr(bd, "fmla"))
  expect_error(balanceTest(des, data = dat), "carries its own data")
  expect_error(balanceTest(des, subset = dat$x1 > 0), "subset")
  expect_equal(dimnames(bt$results)[["strata"]], c("s", "--"))
  ## named arguments still reach the formula method
  d <- dat
  d$w <- 1
  bt2 <- balanceTest(z ~ x1 + x2, data = d, p.adjust.method = "none", unit.weights = w)
  expect_s3_class(bt2, "balancetest")
  expect_equal(unlist(bt2$overall["--", ]),
               unlist(balanceTest(z ~ x1 + x2, data = d, p.adjust.method = "none")$overall["--", ]))
})

test_that("the rank cut depends on neither rounding nor the variables' units (#146)", {
  ## rounding-level noise on a singular covariance does not raise its rank
  noisy <- Sigma + 1e-13 * tcrossprod(rep(1, nrow(Sigma)))
  expect_equal(attr(RItools:::.pseudoinvSqrtOfCovariance(noisy), "r"), 4)
  ## the formula method's omnibus df with a factor is the rank on both rows
  bt <- balanceTest(z ~ x1 + x2 + g + strata(s), data = dat)
  expect_equal(bt$overall[, "df"], c(4, 4))

  ## a count in the thousands beside a proportion: both directions are kept
  ## and the statistic does not depend on the units
  set.seed(146)
  d <- dat
  d$enroll <- round(1000 + 800 * d$x1 + 300 * rnorm(nrow(d)))
  d$prop <- plogis(d$x2) / 20
  bt <- balanceTest(z ~ enroll + prop, data = d)
  expect_equal(bt$overall["--", "df"], 2)
  d$enroll_k <- d$enroll / 1e6
  bt_k <- balanceTest(z ~ enroll_k + prop, data = d)
  expect_equal(bt_k$overall["--", "chisquare"], bt$overall["--", "chisquare"])
  S <- cov(expandCovariates(~ enroll + prop, d))
  expect_gt(S["enroll", "enroll"] / S["prop", "prop"], 1e8)
  bs <- balanceTest(z ~ enroll + prop, data = d, covariances = list("--" = S))
  expect_equal(bs$overall["--", "df"], 2)
  expect_equal(bs$overall["--", "chisquare"], bt$overall["--", "chisquare"])
  ## the wrapper returns a generalized inverse: V V' S V V' = V V'
  V <- RItools:::.pseudoinvSqrtOfCovariance(Sigma)
  G <- tcrossprod(V)
  expect_equal(Sigma %*% G %*% Sigma, Sigma, check.attributes = FALSE)
})

test_that("the formula method's model frame holds the formula's columns (#146)", {
  ## unit.weights and subset expressions built from the caller's variables
  ## work inside a function, where they used to be looked up in the global
  ## environment only
  f <- function(d) {
    w <- rep(2, nrow(d))
    keep <- d$x1 > -10
    list(balanceTest(z ~ x1 + x2 + strata(s), data = d, unit.weights = w),
         balanceTest(z ~ x1 + x2 + strata(s), data = d, subset = keep),
         balanceTest(z ~ x1 + x2 + strata(s), data = d, unit.weights = rep(1, nrow(d))))
  }
  res <- f(dat)
  ref <- balanceTest(z ~ x1 + x2 + strata(s), data = dat)
  for (r in res) {
    expect_equal(r$overall[, "chisquare"], ref$overall[, "chisquare"])
    expect_equal(r$results[, "adj.diff", ], ref$results[, "adj.diff", ])
  }
  ## a subset given by a local vector takes effect
  g <- function(d) { keep <- d$g != "c"; balanceTest(z ~ x1 + x2, data = d, subset = keep) }
  expect_equal(g(dat)$results[, "adj.diff", "--"],
               balanceTest(z ~ x1 + x2, data = dat[dat$g != "c", ])$results[, "adj.diff", "--"])
  ## a term that uses an object that is not a column, of another length
  brks <- c(-Inf, 0, Inf)
  bt <- balanceTest(z ~ cut(x1, brks) + x2, data = dat)
  expect_equal(sum(grepl("^cut", dimnames(bt$results)[["vars"]])), 2)
  ## `.` still means every column, and a column the formula does not name
  ## does not matter
  d2 <- dat[, c("z", "x1", "x2")]
  expect_equal(balanceTest(z ~ ., data = d2)$overall, balanceTest(z ~ x1 + x2, data = d2)$overall)
  d3 <- dat; d3$junk <- list(1)[rep(1, nrow(d3))]
  expect_equal(balanceTest(z ~ x1 + x2, data = d3)$overall, balanceTest(z ~ x1 + x2, data = dat)$overall)
  ## a column whose name is not a syntactic name
  d4 <- dat; d4$`x 3` <- d4$x3; d4$"a`b" <- d4$x2
  expect_equal(unname(balanceTest(z ~ `x 3`, data = d4)$results[1, "adj.diff", "--"]),
               unname(balanceTest(z ~ x3, data = dat)$results[1, "adj.diff", "--"]))
  expect_s3_class(balanceTest(z ~ ., data = d4[, c("z", "x1", "a`b")]), "balancetest")
  ## the comparison issue #74 once left commented out: `subset` against
  ## subsetting the data beforehand. The descriptives agree. The inferentials
  ## are not expected to, since `subset` sets weights to zero and those units
  ## still count (see #74); the model frame was not the reason it failed.
  data(nuclearplants)
  xb4 <- balanceTest(pr ~ . - pt, data = nuclearplants, subset = pt == 0)
  xb5 <- balanceTest(pr ~ . - pt, data = subset(nuclearplants, pt == 0))
  for (stat in c("Treatment", "Control", "adj.diff", "std.diff")) {
    expect_equal(xb4$results[, stat, "--"], xb5$results[, stat, "--"], info = stat)
  }
  ## and `subset` may name a column the formula does not use
  expect_equal(balanceTest(pr ~ date + cap, data = nuclearplants, subset = pt == 0)$results[, "adj.diff", "--"],
               balanceTest(pr ~ date + cap, data = subset(nuclearplants, pt == 0))$results[, "adj.diff", "--"])
  ## no `data` at all: the variables come from where the formula was written
  h <- function() {
    zz <- dat$z; a <- dat$x1; b <- dat$x2; blk <- dat$s
    brks <- c(-Inf, 0, Inf)
    list(balanceTest(zz ~ a + b + strata(blk)),
         balanceTest(zz ~ cut(a, brks) + b))
  }
  nodata <- h()
  ref <- balanceTest(z ~ x1 + x2 + strata(s), data = dat)
  expect_equal(unname(nodata[[1]]$overall[, "chisquare"]), unname(ref$overall[, "chisquare"]))
  expect_equal(unname(nodata[[1]]$results[, "adj.diff", "blk"]), unname(ref$results[, "adj.diff", "s"]))
  expect_s3_class(nodata[[2]], "balancetest")
  ## a formula with no environment still works
  fm <- z ~ x1 + x2; environment(fm) <- NULL
  expect_equal(balanceTest(fm, data = dat)$overall, balanceTest(z ~ x1 + x2, data = dat)$overall)
})

## each stratum's own sample covariance, as a k by k by strata array
strata_covs <- function(fm, data, by) {
  x <- expandCovariates(fm, data)
  lev <- levels(factor(by))
  arr <- array(0, c(ncol(x), ncol(x), length(lev)), dimnames = list(colnames(x), colnames(x), lev))
  for (l in lev) arr[, , l] <- cov(x[by == l, , drop = FALSE])
  arr
}

test_that("supplying each stratum's sample covariance reproduces HB08 exactly", {
  fm <- ~ x1 + x2 + x3
  arr <- strata_covs(fm, dat, dat$s)
  one <- cov(expandCovariates(fm, dat))
  bt <- balanceTest(z ~ x1 + x2 + x3 + strata(s), data = dat)
  bs <- balanceTest(z ~ x1 + x2 + x3 + strata(s), data = dat,
                    covariances = list(s = arr, "--" = one))
  expect_equal(attr(bs, "covariance"), "supplied")
  expect_equal(dimnames(bs$results), dimnames(bt$results))
  expect_equal(bs$overall, bt$overall)
  for (stat in c("adj.diff", "z", "p", "Treatment", "Control", "std.diff")) {
    expect_equal(bs$results[, stat, ], bt$results[, stat, ], info = stat)
  }
  ## a stratification with no entry is tested in the usual way
  bs1 <- balanceTest(z ~ x1 + x2 + x3 + strata(s), data = dat, covariances = list(s = arr))
  expect_equal(bs1$overall, bt$overall)
  ## and the covariances matter: doubling one stratification's halves its chi-square only
  bs2 <- balanceTest(z ~ x1 + x2 + x3 + strata(s), data = dat, covariances = list(s = 2 * arr))
  expect_equal(bs2$overall["s", "chisquare"], bt$overall["s", "chisquare"] / 2)
  expect_equal(bs2$overall["--", "chisquare"], bt$overall["--", "chisquare"])
  ## the design method takes the same argument
  d <- dat; d$`(weights)` <- 1
  des <- RItools:::makeDesigns(z ~ x1 + x2 + x3 + strata(s), d)
  expect_equal(balanceTest(des, covariances = list(s = arr, "--" = one))$overall, bt$overall)
})

test_that("several stratifications each take their own array; a matrix serves every stratum", {
  d <- dat
  d$s2 <- factor(rep(c("p", "q", "r"), length.out = nrow(d)))
  fm <- ~ x1 + x2 + x3
  covs <- list(s = strata_covs(fm, d, d$s), s2 = strata_covs(fm, d, d$s2))
  bt <- balanceTest(z ~ x1 + x2 + x3 + strata(s) + strata(s2) - 1, data = d)
  bs <- balanceTest(z ~ x1 + x2 + x3 + strata(s) + strata(s2) - 1, data = d, covariances = covs)
  expect_equal(bs$overall, bt$overall)
  ## layers are matched to strata by name, whatever their order
  shuffled <- covs
  shuffled$s <- covs$s[, , rev(dimnames(covs$s)[[3]])]
  expect_equal(balanceTest(z ~ x1 + x2 + x3 + strata(s) + strata(s2) - 1, data = d,
                           covariances = shuffled)$overall, bt$overall)
  ## one matrix for every stratum: the null covariance is that matrix times the design scalar
  one <- cov(expandCovariates(fm, d))
  bm <- balanceTest(z ~ x1 + x2 + x3 + strata(s) - 1, data = d, covariances = list(s = one))
  same <- array(one, c(dim(one), nlevels(d$s)), dimnames = c(dimnames(one), list(levels(d$s))))
  ba <- balanceTest(z ~ x1 + x2 + x3 + strata(s) - 1, data = d, covariances = list(s = same))
  expect_equal(bm$overall, ba$overall)
  tc <- attr(bm$overall, "tcov")[["s"]][colnames(one), colnames(one)]
  expect_equal(tc / tc[1, 1], one / one[1, 1])
})

test_that("a factor's columns and a subset of the variables", {
  fm <- ~ x1 + x2 + g
  arr <- strata_covs(fm, dat, dat$s)
  bt <- balanceTest(z ~ x1 + x2 + g + strata(s), data = dat)
  bs <- balanceTest(z ~ x1 + x2 + g + strata(s), data = dat, covariances = list(s = arr))
  expect_equal(bs$overall["s", "df"], 4)       # the indicators sum to one
  expect_equal(bs$overall["s", ], bt$overall["s", ], check.attributes = FALSE)
  ## covariances over more columns than are tested: the block is taken by name
  big <- strata_covs(~ x1 + x2 + x3 + g, dat, dat$s)
  expect_equal(balanceTest(z ~ x1 + x2 + g + strata(s), data = dat, covariances = list(s = big))$overall["s", "chisquare"],
               bt$overall["s", "chisquare"])
})

test_that("malformed covariances are refused with a reason", {
  fm <- ~ x1 + x2 + x3
  arr <- strata_covs(fm, dat, dat$s)
  f <- function(cv) balanceTest(z ~ x1 + x2 + x3 + strata(s), data = dat, covariances = cv)
  expect_error(f(arr), "list with distinct names")
  expect_error(f(list(nope = arr)), "not being tested: nope")
  expect_error(f(list(s = arr[, 1:2, ])), "k by k")
  expect_error(f(list(s = unname(arr))), "names")
  expect_error(f(list(s = arr[c("x1", "x2"), c("x1", "x2"), ])), "x3")
  expect_error(f(list(s = arr[, , 1:3])), "no supplied covariance for stratum level\\(s\\): 4")
  bad <- arr; bad[1, 1, 1] <- NA
  expect_error(f(list(s = bad)), "missing")
  asym <- arr; asym[1, 2, 1] <- asym[1, 2, 1] + 1
  expect_error(f(list(s = asym)), "symmetric")
  expect_error(f(list(s = arr, s = arr)), "distinct names")
  ## an empty list is no covariances at all
  expect_equal(f(list())$overall, balanceTest(z ~ x1 + x2 + x3 + strata(s), data = dat)$overall)
  ## combinations the supplied covariances cannot describe
  expect_error(balanceTest(z ~ x1 + x2 + x3 + strata(s), data = dat, subset = x1 > -1,
                           covariances = list(s = arr)), "Subset the data instead")
  expect_error(balanceTest(z ~ x1 + x2 + x3 + strata(s), data = dat,
                           post.alignment.transform = function(x, w) rank(x),
                           covariances = list(s = arr)), "post.alignment.transform")
})

test_that("a stratum that contributes nothing needs no usable layer", {
  d <- rbind(dat, transform(dat[1, ], s = "lone", z = 1))
  d$s <- factor(d$s)
  fm <- ~ x1 + x2 + x3
  arr <- strata_covs(fm, d, d$s)                 # the singleton's layer is all NA
  expect_true(all(is.na(arr[, , "lone"])))
  bt <- suppressWarnings(balanceTest(z ~ x1 + x2 + x3 + strata(s) - 1, data = d))
  bs <- suppressWarnings(balanceTest(z ~ x1 + x2 + x3 + strata(s) - 1, data = d, covariances = list(s = arr)))
  expect_equal(bs$overall, bt$overall)
  ## and it may be left out of the array altogether
  bs2 <- suppressWarnings(balanceTest(z ~ x1 + x2 + x3 + strata(s) - 1, data = d,
                                      covariances = list(s = arr[, , levels(dat$s)])))
  expect_equal(bs2$overall, bt$overall)
})

test_that("clusters run through the supplied-covariance calculation", {
  d <- dat
  d$s <- factor(rep(1:4, each = nrow(d) / 4))
  d$cl <- factor(rep(seq_len(nrow(d) / 4), each = 4))
  d$z <- as.numeric(as.integer(d$cl) %% 2 == 0)
  one <- cov(expandCovariates(~ x1 + x2, d))
  bs <- balanceTest(z ~ x1 + x2 + strata(s) + cluster(cl), data = d, covariances = list(s = one, "--" = one))
  expect_true(all(is.finite(bs$overall[, "chisquare"])))
  expect_equal(bs$overall[, "df"], c(2, 2))
})

test_that("missing covariate values reduce the null variance as the sample calculation does", {
  set.seed(3)
  n <- 6000
  big <- data.frame(z = rbinom(n, 1, 0.5), x1 = rnorm(n), s = factor(rep(1:6, length.out = n)))
  big$x2 <- 0.6 * big$x1 + rnorm(n)
  Strue <- matrix(c(1, 0.6, 0.6, 1.36), 2, 2, dimnames = list(c("x1", "x2"), c("x1", "x2")))
  arr <- array(Strue, c(2, 2, 6), dimnames = c(dimnames(Strue), list(levels(big$s))))
  full <- balanceTest(z ~ x1 + x2 + strata(s) - 1, data = big, covariances = list(s = arr))
  big$x1[sample(n, 0.3 * n)] <- NA
  miss <- balanceTest(z ~ x1 + x2 + strata(s) - 1, data = big, covariances = list(s = arr))
  ref <- balanceTest(z ~ x1 + x2 + strata(s) - 1, data = big)
  pick <- function(b) attr(b$overall, "tcov")[["s"]][c("x1", "x2"), c("x1", "x2")]
  expect_true(all(abs(pick(miss) / pick(ref) - 1) < 0.1))
  expect_equal(pick(miss)["x1", "x1"] / pick(full)["x1", "x1"], 0.7, tolerance = 0.03)
  expect_equal(pick(miss)["x2", "x2"], pick(full)["x2", "x2"])
  vars <- dimnames(miss$results)[["vars"]]
  expect_true(all(is.na(miss$results[startsWith(vars, "("), "z", "s"])))
  expect_equal(miss$overall["s", "df"], 2)
})
