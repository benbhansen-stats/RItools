##' xBalance helper function
##'
##' Make engine
##' @param ss ss
##' @param zz zz
##' @param mm mm
##' @param report report
##' @param swt swt
##' @param s.p s.p
##' @param normalize.weights normalize.weights
##' @param zzname zzname
##' @param post.align.trans post.align.trans
##' @param p.adjust.method Method to adjust P.
##' @return List
balanceTestEngine <- function(ss,zz,mm,report, swt, s.p, normalize.weights, zzname, post.align.trans, p.adjust.method) {
  ##ss is strata, zz is treatment, mm is the model matrix defined by the formula and data input to xBalance, swt is stratum weights, s.p. is the pooled sd, normalize.weights is logical (for creation of stratum weights)

  cnms <-
    c(
      if ('adj.means'%in%report) c("Control", "Treatment") else character(0), ##c("Tx.eq.0","Tx.eq.1") else character(0),
      if ('adj.mean.diffs'%in%report) 'adj.diff' else character(0),
      if ('std.diffs'%in%report) 'std.diff' else character(0),
      if ('z.scores'%in%report) 'z' else character(0),
      if ('p.values'%in%report) 'p' else 'p'#character(0) turns out that it may be useful to have p-values in the object whether or not they are requested for printing
      )


  ##Set up a data frame to contain the results
  ans <-
    if (length(setdiff(report,"chisquare.test"))) {
      as.data.frame(matrix(0,dim(mm)[2], length(cnms),
			   dimnames= list(dimnames(mm)[[2]], cnms)))
    } else {
      as.data.frame(matrix(0,0,0))
    }

  if (!length(zz)) {
    return(list(dfr=as.data.frame(matrix(0,0, length(cnms),
					 dimnames=list(character(0), cnms))),
		chisq=c(pre.chisquare=0,pre.df=0,
			post.chisquare=0,post.df=0)))
  }


  ##Number of strata
  nlev <- nlevels(ss)

  ### Calculate post.difference
  ZtH <- unsplit(tapply(zz,ss,mean),ss) ##proportion treated (zz==1) in strata s.
  ssn <- drop(crossprod(zz-ZtH, mm*swt$wtratio)) ##weighted sum of mm in treated (using centered treatment)
  wtsum <- sum(unlist(tapply(zz,ss,function(x){var(x)*(length(x)-1)}))) ## h=(m_t*m_c)/m
  post.diff <- ssn/(wtsum) ##diff of means
  if ('adj.mean.diffs'%in%report) ans[['adj.diff']] <- post.diff
  if ('std.diffs' %in% report) ans[['std.diff']] <- post.diff/s.p

  ### Calculate post.Tx.eq.0, post.Tx.eq.1 --- now called "the treatment var"=0 and =1
  if ("adj.means"%in%report) 	{
    postwt0 <- unsplit(swt$sweights/tapply(zz<=0, ss, sum),
		       ss[zz<=0], drop=TRUE)
    ans[["Control"]] <- apply(mm[zz<=0,,drop=FALSE]*postwt0, 2,sum)
    ans[["Treatment"]] <- ans[["Treatment"]] + post.diff
  }


  msmn <- xBalance.make.stratum.mean.matrix(ss, mm)

  tmat <- (mm - msmn)
  ##dv is sample variance of treatment by stratum
  dv <- unsplit(tapply(zz,ss,var),ss)
  ssvar <- apply(dv*swt$wtratio^2*tmat*tmat, 2, sum) ## for 1 column in  mm, sum(tmat*tmat)/(nrow(tmat)-1)==var(mm) and sum(dv*(mm-mean(mm))^2)=ssvar or wtsum*var(mm)

  if (!is.null(post.align.trans)) {
    # Transform the columns of tmat using the function in post.align.trans
    tmat.new <- apply(tmat, 2, post.align.trans)
    # Ensure that post.align.trans wasn't something that changes the size of tmat (e.g. mean).
    # It would crash later anyway, but this is more informative
    if (is.null(dim(tmat.new)) || !all(dim(tmat) == dim(tmat.new))) {
      stop("Invalid post.alignment.transform given")
    }
    ## recenter on stratum means
    tmat <- tmat.new
    msmn <- xBalance.make.stratum.mean.matrix(ss, tmat)
    tmat <- tmat - msmn
    tmat <- tmat *swt$wtratio
    # Recalculate on transformed data the difference of treated sums and their null expectations
    # (NB: since tmat has just been recentered,
    # crossprod(zz,tmat) is the same as crossprod(zz-ZtH,tmat))
    ssn <- drop(crossprod(zz, tmat))
    ssvar <- apply(dv*tmat*tmat, 2, sum)
  } else {
      tmat <- tmat *swt$wtratio
  }


    if ('z.scores' %in% report) {
    ans['z'] <- ifelse(ssvar<=.Machine$double.eps,0, ssn/sqrt(ssvar))
  }
  if (any(c("adj.means","adj.mean.diffs","std.diffs","z.scores","p.values") %in% report)) {
    ##always produce a pvalue to use to create signif stars.
    ans['p'] <- p.adjust(
      ifelse(ssvar <= .Machine$double.eps,
             1,
             2 * pnorm(abs(ssn/sqrt(ssvar)), lower.tail=FALSE)),
      p.adjust.method)
  }

  if ("chisquare.test" %in% report && any(ssvar > .Machine$double.eps)) {
    # nu=0 stops calculation of the U matrix, which is not used.
    pst.svd <- try ( svd(tmat*sqrt(dv), nu=0) )
    if (inherits(pst.svd,'try-error')) {
      pst.svd<-propack.svd(tmat*sqrt(dv))
    }
    ##	pst.svd <- svd(tmat*sqrt(dv))
    Positive <- pst.svd$d > max(sqrt(.Machine$double.eps)*pst.svd$d[1], 0)
    Positive[is.na(Positive)]<-FALSE # JB Note: Can we imagine a situation in which we dont want to do this?
    if (all(Positive)) { ## is this faster? { ytl <- sweep(pst.svd$v,2,1/pst.svd$d,"*") }
      ytl <- pst.svd$v *
      matrix(1/pst.svd$d, nrow=dim(mm)[2],ncol=length(pst.svd$d), byrow=T)
    } else if (!any(Positive)) {
      ytl <- array(0, dim(mm)[2:1] )
    } else  {
      ytl <- pst.svd$v[, Positive, drop = FALSE] *
      matrix(1/pst.svd$d[Positive],ncol=sum(Positive),nrow=dim(mm)[2],byrow=TRUE)
    }

    mvz <- drop(crossprod(zz, tmat)%*%ytl)

    csq <- drop(crossprod(mvz))
    DF <- sum(Positive)
    tcov <- crossprod(sqrt(dv) * tmat * (1 / wtsum))

  } else {
    csq <- DF <- 0
    tcov <- matrix(0, ncol = length(ssvar), nrow = length(ssvar),
                   dimnames = list(names(ssvar), names(ssvar)))
  }

  list(dfr   = ans,
       chisq = c('chiquare' = csq, 'df' = DF),
       tcov  = tcov) # test statistic covariance matrix
}



xBalance.make.stratum.mean.matrix <- function(ss, mm) {
  stopifnot(inherits(ss, "factor")) # just in case a numeric variable is passed in.

  post.nobs <- dim(mm)[1]
  nlev <- nlevels(ss)

  # for this matrix, finding the indices of the rows is easy, as there is only one
  # item per row, and there post.nobs number of rows.
  tR <- new("matrix.csr",
            ja = as.integer(as.integer(ss)),
            ia = as.integer(1:(post.nobs+ 1)),
            ra = unsplit(1/tapply(ss,ss,length),ss),
            dimension = c(post.nobs,nlev))

  # With many items per row, we need to break the list of strata
  # down into indices of where each row starts and ends
  # e.g. Say ss = 0 0 0 0 1 1 1 1 1 1 1 2 2 2 2 the row indices would be
  # 1 5 12 16 (where 16 is the start of the non existant 4th row)

  L <- new("matrix.csr", #ifelse(oldver,"tripletMatrix","dgTMatrix"),
           ia = as.integer(1:(post.nobs + 1)),
           ja = as.integer(as.integer(ss)),
           ra = rep(1,length(ss)),
           dimension = c(post.nobs,nlev))


  msmn <- t(tR) %*% mm
  msmn <- L %*% msmn

  msmn <- as.matrix(msmn)

  return(msmn)
}
