##' @rdname balanceTest
##' @order 3
##' @method balanceTest DesignOptions
##' @export
balanceTest.DesignOptions <- function(fmla, data = NULL,
                                      stratum.weights = harmonic_times_mean_weight,
                                      include.NA.flags = TRUE,
                                      covariate.scales = setNames(numeric(0), character(0)),
                                      post.alignment.transform = NULL,
                                      inferentials.calculator = HB08,
                                      p.adjust.method = "holm",
                                      covariances = NULL,
                                      ...) {
  .stopOnUnusedArguments(match.call(expand.dots = FALSE)$...)
  if (!is.null(data)) {
    stop("a design carries its own data; `data` is for the formula method")
  }
  stopifnot(is.null(post.alignment.transform) || is.function(post.alignment.transform))
  design <- fmla
  if (is.null(p.adjust.method)) {
    p.adjust.method <- "none"
  }

  ## we used to allow select report options, not any more.
  report <- c("adj.means","adj.mean.diffs","chisquare.test", "std.diffs","z.scores","p.values")

  ## Which of the NM cols to look at for a given variable's NM info
  NMpatterns <-   c("_any Xs recorded_", colnames(design@NotMissing))[1L+design@NM.Covariates]
  NMpatterns <- paste0("(",NMpatterns,")")

  aggDesign       <- aggregateDesigns(design)
  ## (Creation of stratum weightings for use in
  ##  descriptives calculations would go here, if
  ## we wanted to allow departures from the ETT default.
  ## Something like `design@Sweights <- DesignWeights(aggDesign, <...>)`.)
  descriptives    <- designToDescriptives(design, covariate.scales)
  ## The notmissing indicators can't themselves be missing, so no role for
  ## another notmissing indicator to record where its missingness pattern.
  ## Instead we just record a blank for them:
  NMpatterns <- c(NMpatterns, rep("", dim(descriptives)[1]-length(NMpatterns)))

  # these weights govern inferential but not descriptive calculations

  aggDesign <- as(aggDesign, "StratumWeightedDesignOptions")
  aggDesign@Sweights <-
      DesignWeights(aggDesign, stratum.weights)

  strataAligned <- sapply(colnames(aggDesign@StrataFrame),
                          alignDesignsByStrata,
                          design=aggDesign,
                          post.align.transform=post.alignment.transform,
                          simplify = FALSE, USE.NAMES = TRUE)

  origvars <- strataAligned[[1]]@OriginalVariables #to include NotMissing columns

  ## One set of inferentials per stratification: the usual calculator, or the
  ## supplied-covariance one where a covariance array was given for it.
  covariances <- .matchCovariances(covariances, names(strataAligned))
  if (any(!vapply(covariances, is.null, logical(1)))) {
    if (!is.null(post.alignment.transform)) {
      stop("covariances cannot be combined with post.alignment.transform: ",
           "they describe the untransformed covariates")
    }
    if (any(design@UnitWeights == 0)) {
      stop("covariances cannot be combined with zero unit weights, which is how `subset` ",
           "is implemented: units of weight zero would still be counted. Subset the data instead.")
    }
  }
  tmp <- mapply(function(aligned, covariance) {
    if (is.null(covariance)) inferentials.calculator(aligned) else HB08_supplied(aligned, covariance)
  }, strataAligned, covariances, SIMPLIFY = FALSE)
  names(tmp) <- colnames(aggDesign@StrataFrame)

  ans <- list()

  # append the z and p to the "descriptives" array (making it somewhat misnamed)
  tmp.z <- as.data.frame(lapply(tmp, function(tt) { tt$z }))
  tmp.p <- as.data.frame(lapply(tmp, function(tt) { tt$p }))
  nstats.previous <- dim(descriptives)[2]
  descriptives <- abind(descriptives, along = 2, tmp.z, tmp.p, use.first.dimnames = TRUE)
  names(dimnames(descriptives)) <- c("vars", "stat", "strata")

  dimnames(descriptives)[[2]][nstats.previous + 1:2] <- c("z", "p")

  inferentials <- do.call(rbind, lapply(tmp, function(s) {
    data.frame(s$Msq, s$DF, pchisq(s$Msq, df = s$DF, lower.tail = FALSE))
  }))
  colnames(inferentials) <- c("chisquare", "df", "p.value")

  tcovs <- lapply(tmp, function(r) {
    tcov <- r$tcov
    dimnames(tcov) <-
      rep(list(dimnames(descriptives)[["vars"]]),2)
    tcov
  })

  # Deal with summaries of not-missing indicators that only ever take the value T
  nmvars <- identify_NM_vars(dimnames(descriptives)[["vars"]])
  # next line assumes every "stat" not in the given list is a mean
  # over a group assigned to some treatment condition.
  group_mean_labs <- setdiff(dimnames(descriptives)[["stat"]],
                             c("std.diff", "adj.diff", "pooled.sd", "z", "p"))
  if (length(nmvars) & length(group_mean_labs)) #cf #111
  {
	  groupmeans <- descriptives[nmvars, group_mean_labs,,drop=FALSE]
	  bad <- apply(abs(groupmeans - 1) < sqrt(.Machine$double.eps), 1, all)
	  toremove <- match(nmvars[bad], dimnames(descriptives)[["vars"]])
	  if(length(toremove)>0){ ## if toremove=integer(0) then it drops all vars from descriptives
		  descriptives <- descriptives[-toremove,,,drop=FALSE]
	    tcovs <- lapply(tcovs,
        function(mat){mat[-toremove, -toremove, drop=FALSE]}
        )
      origvars <- origvars[-toremove]

      ## also trim covars to NM patterns lookup table:
      NMpatterns <- NMpatterns[-toremove]
      ## If the NM pattern var that a covar got mapped to
      ## was removed, we take the covar to have no missings
      ## anywhere, and report "" as the name of its nonmissing
      ## indicator variable:
      strings_to_remove <- dimnames(descriptives)[["vars"]][toremove]
      NMpatterns[ NMpatterns%in% strings_to_remove] <- ""
	  }
  }

  for (s in 1L:dim(descriptives)[3]) {
  descriptives[, "p", s] <- 
    p.adjust(descriptives[, "p", s], method = p.adjust.method)
  }
  attr(descriptives, "NMpatterns") <- NMpatterns
  attr(descriptives, "originals") <- origvars
  attr(descriptives, "term.labels") <- design@TermLabels
# Now `attr(d, "term.labels")[attr(d, "originals")]`
# associates variables to terms of fmla
  attr(descriptives, "include.NA.flags") <- 
    include.NA.flags # hint to print and plot methods

# the meat of our xbal object
  ans$overall <- inferentials
  attr(ans$overall, "tcov") <- tcovs
  ans$results <- descriptives

  if (!is.null(covariances) && any(!vapply(covariances, is.null, logical(1)))) {
    attr(ans, "covariance") <- "supplied"
  }
  attr(ans, "report") <- report # hinting to our summary method later
  class(ans) <- c("balancetest", "xbal", "list")
  ans
}

