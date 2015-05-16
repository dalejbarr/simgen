#' Fit model with error-trapping
#' 
#' Try to fit an lmer model, catching any warnings about convergence
#' 
#' @param tf.formula The model formula
#' @param tf.data The data frame
#' @param ... Arguments to be passed onto the fitting function \code{\link{lmer}}
#' @return A list with:
#' \item{value}{fitted model object from lmer}
#' \item{converged}{whether or not the model converged}
#' @seealso \code{\link{tryUpdate}}
#' @examples
#' 
#' # random intercepts only
#' ff <- tryFit(Resp~Cond + (1+Cond|SubjID)+(1+Cond|ItemID), mkDf())
#'
#' @export tryFit
tryFit <- function(tf.formula, tf.data, ...) {
    converged <- TRUE
    w.handler <- function(w) {
        converged <- FALSE
        invokeRestart("muffleWarning")
    }
    arg.list <- c(list(formula=tf.formula, data=tf.data), list(...))
    list(value=withCallingHandlers(tryCatch(
             do.call(lmer, arg.list),
             error=function(e) e),
             warning=w.handler),
         converged=converged)
}


#' Update model with error-trapping
#' 
#' Try to update an lmer model, catching any errors/warnings
#' 
#' @param tf.formula The model formula
#' @param tf.model The data frame
#' @param ... Further arguments to be passed onto the function \code{\link{update}}
#' @return A list with:
#' \item{value}{fitted model object from lmer}
#' \item{converged}{whether or not the model converged}
#' @seealso \code{\link{tryFit}}
#' @examples
#' mod1 <- tryFit(Resp~Cond+(1+Cond|SubjID)+(1+Cond|ItemID), mkDf())
#' mod2 <- tryUpdate(.~.-Cond, mod1)
#' lrTest(mod1, mod2)
#'
#' @export tryUpdate
tryUpdate <- function(tu.formula, tu.model, ...) {
    converged <- TRUE
    w.handler <- function(w) {
        converged <- FALSE
        invokeRestart("muffleWarning")
    }
    m1 <- tu.model$value
    arg.list <- c(list(object=m1, formula=tu.formula), list(...))
    list(value=withCallingHandlers(tryCatch(
             do.call(update, arg.list),
             error=function(e) e),
             warning=w.handler),
         converged=converged)
}


#' Proportion of runs achieving statistical significance
#'
#' @param x Vector or matrix where each row is a single run
#' @param alpha alpha level (default=.05)
#' @return A vector with the proportion of significant runs.
#' @seealso \code{\link{pconv}}
#' @examples
#'
#' nmc <- 10
#' pmx <- cbind(randParams(genParamRanges(), nmc, 1001), seed=mkSeeds(nmc, 1001))
#'
#' res <- mcRun("fitanova", mcr.fnArgs=list(wsbi=FALSE),
#'              mcr.datFn="mkDf", mcr.datArgs=list(wsbi=FALSE),
#'              mcr.varying=pmx)
#' 
#' psig(res[,c("p1","p2","pmf","pmax")])
#' 
#' @export psig
psig <- function(x, alpha=.05) {
    res <- NULL
    if (is.vector(x)) {
        res <- sum(x<=alpha)/sum(!is.na(x))
    } else {
        if (!is.matrix(x) && !is.data.frame(x)) {
            stop("argument 'x' must be a vector of a matrix: was ", class(x))
        } else {}
        res <- apply(x, 2, function(xx) {sum(xx<=alpha)/sum(!is.na(xx))})
    }
    return(res)
}


#' Proportion of runs achieving convergence
#'
#' @param x Vector or matrix where each row is a single run
#' @return A vector with the proportion of converging runs.
#' @seealso \code{\link{psig}}
#' @examples
#'
#' nmc <- 10
#' pmx <- cbind(randParams(genParamRanges(), nmc, 1001), seed=mkSeeds(nmc, 1001))
#'
#' res <- mcRun("fitanova", mcr.fnArgs=list(wsbi=FALSE),
#'              mcr.datFn="mkDf", mcr.datArgs=list(wsbi=FALSE),
#'              mcr.varying=pmx)
#' 
#' pconv(res[,c("p1","p2","pmf","pmax")])
#' 
#' @export pconv
pconv <- function(x) { # probability of convergence
    res <- NULL
    if (is.vector(x)) {
        res <- sum(!is.na(x))/length(x)
    } else {
        if (!is.matrix(x) && !is.data.frame(x)) {
            stop("argument 'x' must be a vector of a matrix: was ", class(x))
        } else {}
        res <- apply(x, 2, function(xx) {sum(!is.na(xx))/length(xx)})
    }
    return(res)
}


#' Make random seeds for Monte Carlo runs
#'
#' Creates a vector of unique random seeds for reproducible generation of datasets.
#'
#' @param nmc Number of Monte Carlo runs (default 1000)
#' 
#' @param firstseed Seed for generating seeds (meta-seed?)
#' 
#' @return A vector of seeds up to .Machine$integer.max
#' 
#' @examples
#'
#' nmc <- 10
#' pmx <- cbind(randParams(genParamRanges(), nmc, 1001), seed=mkSeeds(nmc, 1001))
#' @export mkSeeds
mkSeeds <- function(nmc=1000, firstseed=NULL) {
    randSeed <- function(n=1) {
                                        # seeds must all be unique!
        seeds <- c()
        nremaining <- n
        nTries <- 1
        nMaxTries <- 1000
        while (nremaining & (nTries < nMaxTries)) {
            seeds <- unique(c(seeds, round((.Machine$integer.max-1)*(runif(nremaining, min=0, max=1)),0)))
            nremaining <- n-length(seeds)
            nTries <- nTries + 1
        }
        if (nTries == nMaxTries) {
            stop("couldn't create enough unique random seeds after 1000 tries!")
        } else {}
        return(seeds)
    }
    if (!is.null(firstseed)) {
        set.seed(firstseed)
    } else {}
    randSeed(nmc)    
}


#' Deviation-Coded Contrast Matrices
#'
#' Return a matrix of deviation-coded contrasts.
#'
#' @param n a vector of levels for a factor, or the number of levels.
#' @param base an integer specifying which group is considered the
#' baseline group. Ignored if ‘contrasts’ is ‘FALSE’.
#' @param contrasts a logical indicating whether contrasts should be computed.
#'
#' @export 
contr.dev <- function(n, base = 1, contrasts = TRUE) {
    if (length(n) <= 1L) {
        if (is.numeric(n) && (length(n) == 1L) && (n > 1L))
            levels <- seq_len(n)
        else stop("not enough degrees of freedom to define contrasts")
    } else {
        levels <- n
    }
    #mx <- apply(contr.treatment(n), 2, function(x) {x-mean(x)})
    ctreat <- contr.treatment(levels, base, contrasts)
    mx <- apply(ctreat, 2, scale, scale = FALSE)
    dimnames(mx) <- dimnames(ctreat)
    mx    
}

# NOT EXPORTED!
getLmer.pValue <- function(m1,m2) {
    pval <- c(chisq=NA, df=NA, p=NA)
    if (m1$converged && m2$converged) {
        df1 <- length(fixef(m1$value))+sum(unlist(lapply(VarCorr(m1$value), function(x) {dim(x)[1]})))
        df2 <- length(fixef(m2$value))+sum(unlist(lapply(VarCorr(m2$value), function(x) {dim(x)[1]})))
        pval["df"] <- abs(df1-df2)
        pval["chisq"] <- abs((-2*logLik(m1$value))-(-2*logLik(m2$value)))
        pval["p"] <- as.numeric(pchisq(pval["chisq"], pval["df"], lower.tail=FALSE))
    } else {}
    return(pval)
}

#' Likelihood-Ratio Test for Fitted mer Objects
#'
#' @param m1 Model to be compared (list object from \code{\link{tryFit}})
#' @param m2 Model to be compared (list object from \code{\link{tryFit}})
#'
#' @return A list with:
#' \item{Number of Estimated Parameters}{A matrix;}
#' \item{Likelihood-Ratio Test}{Output from the test, or \code{NA} if either or both models did not converge.}
#' @seealso \code{\link{tryFit}}
#' @examples
#' dat <- mkDf()
#' m1 <- tryFit(Resp~Cond+(1+Cond|SubjID)+(1+Cond|ItemID), dat, REML=FALSE)
#' m2 <- tryFit(Resp~(1+Cond|SubjID)+(1+Cond|ItemID), dat, REML=FALSE)
#' lrTest(m1, m2)
#' @export lrTest
lrTest <- function(m1, m2) {
    m1v <- m1$value; m2v <- m2$value
    # calculate no. estimated parameters
    m1d <- m1v@dims
    m2d <- m2v@dims
    m1.df <- c(rand=m1d["nt"], fixed=m1d["p"], link=m1d["s"])
    m1.df <- c(m1.df, tot=sum(m1.df))
    m2.df <- c(rand=m2d["nt"], fixed=m2d["p"], link=m2d["s"])
    m2.df <- c(m2.df, tot=sum(m2.df))
    df.mx <- rbind(m1=m1.df, m2=m2.df)
    colnames(df.mx) <- c("random","fixed","link","total")
    # get estimated log likelihood (deviance)
    if (m1$converged && m2$converged) {
        dv.m1 <- c(m1=-2*logLik(m1v),
                   m2=-2*logLik(m2v))
        dv.m1 <- c(dv.m1, chisq=as.numeric(abs(dv.m1["m2"]-dv.m1["m1"])))
        diff.df <- abs(df.mx["m2","total"]-df.mx["m1","total"])
        dv.m1 <- c(dv.m1, df=diff.df)
        dv.m1 <- c(dv.m1, p=as.numeric(pchisq(abs(dv.m1["chisq"]), diff.df, lower.tail=FALSE)))
        res <- data.frame(deviance1=round(dv.m1["m1"],3), deviance2=round(dv.m1["m2"],3),
                          chisq=round(dv.m1["chisq"],3), df=as.integer(dv.m1["df"]), p=dv.m1["p"])
    } else {
        res <- data.frame(deviance1=NA, deviance2=NA, chisq=NA, df=NA, p=NA)
    }
    list(`Number of Estimated Parameters`=df.mx,
         `Likelihood-Ratio Test`=res)
}

.getLmer.chisq <- function(m1,m2) {
    pval <- NA
    if (m1$converged && m2$converged) {
        chisq.val <- (-2*logLik(m1$value))-(-2*logLik(m2$value))
        pval <- as.numeric(pchisq(abs(chisq.val), 1, lower.tail=FALSE))
    } else {}
    return(pval)
}

#' Fit mixed-models and obtain p-values from model comparison
#'
#' Compare a set of models to a single 'base' model.  The base model
#' is fit using \code{\link{tryFit}}, and all subsequent models using
#' \code{\link{tryUpdate}}.  Models are compared using the
#' likelihood-ratio test.
#'
#' This function is useful for testing multiple IVs in complex designs
#' (e.g., all main effects and interactions in a factorial design).
#'
#' @param mcr.data Data to be analyzed
#' @param lrc.mods List containing model formulae.  The first element
#' specifies the 'base' model that will be compared to all subsequent
#' models (specified, in turn, by the remaining elements of the list).
#' The base model is fit by a call to \code{tryFit}, while the
#' remaining models are fit by \code{tryUpdate}, and thus should have
#' formulae as defined in \code{\link{update.formula}}.
#' @param lrc.modCompare A vector of mode 'character'.  The first
#' element names the 'base' model listed in \code{lrc.mods}, and the
#' remaining elements name the models to be compared to 'base'.  If
#' NULL (default) then automatically compares \code{lrc.mods[1]}
#' to all subsequent models in \code{lrc.mods}.
#' @param lrc.method Currently, no effect.  This argument has been
#' reserved for future development
#' @param lrc.cluster A processing cluster (typically derived from
#' \code{\link{makeCluster}} in the \code{parallel} package.  NULL
#' means run on a single core.
#' @param ... Arguments to be passed along to \code{\link{tryFit}}
#' (and thus to \code{lmer}).
#' @return A vector with the p-values from tests specified by
#' \code{lrc.modCompare}. If a test could not be performed (e.g.,
#' because one or both of the corresponding models did not converge),
#' \code{NA} is returned.
#' @seealso \code{\link{tryFit}}, \code{\link{tryUpdate}}, \code{\link{update}}, \code{\link{update.formula}}
#'
#' @examples
#' dat <- mkDf()
#' 
#' lrCompare(dat,
#'           lrc.mods=c(max=Resp~Cond+(1+Cond|SubjID)+(1+Cond|ItemID),
#'                      max.noC=.~.-Cond,
#'           lrc.modCompare=c("max","max.noC"),
#'           REML=FALSE, na.action=na.omit) # pass to tryFit/lmer
#' 
#' lrCompare(dat,  # same effect as above, with modCompare set to NULL
#'           lrc.mods=c(max=    Resp~Cond+(1+Cond|SubjID)+(1+Cond|ItemID),
#'                      max.noC=Resp~     (1+Cond|SubjID)+(1+Cond|ItemID)),
#'           lrc.modCompare=NULL, # compare all to first
#'           REML=FALSE, na.action=na.omit) # pass to tryFit/lmer
#'
#' dat <- mkDf.facMixedAB(nitem=24)
#' 
#' lrCompare(dat,
#'           lrc.mods=list(
#'               base=Y~Ad*Bd+(1+Bd|SubjID)+(1+Ad*Bd|ItemID),
#'               A=.~.-Ad, B=.~.-Bd, AB=.~.-Ad:Bd),
#'           REML=FALSE, na.action=na.omit)
#' @export lrCompare
lrCompare <- function(mcr.data, lrc.mods, lrc.modCompare=NULL, lrc.method="update",
                      lrc.cluster=NULL, ...) {
    if (length(lrc.mods) < 2) {
        stop("'lrc.mods' must contain at least two model formulae")
    } else {}
    if (is.null(names(lrc.mods))) { # no names; create them
        if (!is.null(lrc.modCompare)) {
            stop("names provided in 'lrc.modCompare', but elements of 'lrc.mods' were not named")
        } else {}
        names(lrc.mods) <- c("base", paste("M", 1:(length(lrc.mods)-1), sep=""))
    } else {}
    if (is.null(lrc.modCompare)) {
        lrc.modCompare <- c(names(lrc.mods)[1], names(lrc.mods)[2:length(lrc.mods)])
    } else {}
    # error checking...
    modsNotFit <- setdiff(lrc.modCompare, names(lrc.mods))
    if (length(modsNotFit)!=0) {
        stop("model(s) '", paste(lrc.modCompare, collapse=", "), "' not in 'lrc.mods'")
    } else {}    
    base.mod <- do.call(tryFit, c(list(tf.formula=lrc.mods[[1]], tf.data=mcr.data), list(...)))
    if (is.null(lrc.cluster)) {
        modres <- lapply(lrc.mods[-1], tryUpdate, tu.model=base.mod)
    } else {
        initializeCluster(lrc.cluster)
        modres <- parallel::parLapplyLB(lrc.cluster, lrc.mods[-1], tryUpdate, tu.model=base.mod)
    }
    ff <- do.call("rbind", lapply(lrc.modCompare[-1], function(y) {
        lrTest(base.mod, modres[[y]])[["Likelihood-Ratio Test"]]
    }))
    rownames(ff) <- paste(lrc.modCompare[1], lrc.modCompare[-1], sep=" vs ")
    return(ff)
}

#' Generate population parameters from data-generating parameters
#'
#' Takes a list of ranges for data-generating parameters and randomly
#' instantiates parameters for a series of populations, drawing each
#' value from a uniform distribution.
#'
#' @param plist Named list of parameters (e.g., from
#' \code{\link{genParamRanges}}), with each element of the list a
#' vector that defines the parameter's range (min, max).  If the
#' vector has only a single element, that parameter is treated as
#' constant.
#' @param nmc Number of populations to create
#' @param firstseed Seed the random number generator before generating
#' the populations (unless \code{NULL}).
#' @return A matrix where each row has the parameter values for a population
#'
#' @examples
#' randParams(genParamRanges(), 2, 1001)
#' randParams(genParamRanges(), 2, 1001) # same result as above
#' randParams(genParamRanges(), 10) # different
#' 
#' @export randParams
randParams <- function(plist, nmc=1, firstseed=NULL) {
    if (!is.null(firstseed)) {
        set.seed(firstseed)
    } else {}
    if (is.vector(plist)) {
        plist <- as.list(plist)
    } else {}
    if (is.null(names(plist))) {
        names(plist) <- paste("V", 1:length(plist), sep="")
    } else {}
    ff <- lapply(names(plist), function(nx) {
        x <- plist[[nx]]
        if (is.list(x)) {
            rge <- x[[sample(1:length(x),1)]] # randomly choose ranges with eq prob
        } else {
            rge <- x
        }
        if (length(rge)==1) {
            res <- rep(rge, nmc)
        } else {
            res <- runif(nmc, rge[1], rge[2])
        }
    })
    matrix(unlist(ff), nrow=nmc, dimnames=list(NULL, names(plist)))
}


#' Find valid parameters for a variance-covariance matrix
#'
#' To generate data from a multivariate normal distribution using
#' \code{mvrnorm} in the \code{MASS} package, it is necessary to
#' specify a square variance-covariance matrix that is symmetric and
#' 'positive definite.'
#'
#' This function takes data-generating parameter ranges (variances and
#' covariances) and randomly generates values for a single population.
#' Use a brute-force method, so to avoid hanging, set \code{maxtries}
#' to a reasonable value.
#'
#' A 3x3 matrix corresponds to 3 variances on the diagonal, and 6
#' variances (forming the upper or lower triangle).  Stated generally,
#' an NxN matrix requires \code{sum(1:(N-1))} covariances.
#'
#' @param vars (named) list of ranges of variances; if not named,
#' given names V1..VN
#' @param corrs (named) list of correlations (range -1 to 1); if not
#' named, given names C1..CN.  List them in the following order: take
#' the upper triangle of the matrix, list all elements in the first
#' row (left to right), then the second row, etc.
#' @param maxtries the number of attempts before giving up
#' @return A vector of variances and correlations
#' @examples
#' pars <- findValidVCovParams(vars=list(A=c(1,3),B=c(1,3), C=c(0,2)),
#'                             corrs=list(AB=c(-.5,.5),AC=c(-.5,.5),BC=c(-.5,.5)))
#'
#' formVCovMx(pars[1:3], pars[4:6])
#' @seealso \code{\link{formVCovMx}}
#' @export findValidVCovParams
findValidVCovParams <- function(vars, corrs=as.list(rep(0,sum(1:(length(vars)-1)))),
                                maxtries=10000) {
    if (is.vector(vars)) {
        vars <- as.list(vars)
    } else {}
    if (is.null(names(vars))) {
        names(vars) <- paste("V", 1:length(vars), sep="")
    } else {}
    if (is.null(names(corrs))) {
        names(corrs) <- paste("C", 1:length(corrs), sep="")
    } else {}
    lapply(corrs, function(x) {
        lapply(x, function(y) {
            if (y<(-1) || y>1) {
                stop("correlations must be between -1 and 1")
            } else {}
        })
    })
    ntries <- 0
    done <- FALSE
    # use brute force to find a definite positive matrix
    while (!done) {
        if (!is.null(maxtries)) {
            if (ntries >= maxtries) {
                stop("'maxtries' reached without finding a positive-definite matrix")
            } else {}
        } else {}
        parms <- randParams(c(vars, corrs))
        vcov.mx <- formVCovMx(parms[1:length(vars)], parms[-(1:length(vars))])
       # test whether the matrix is positive-definite
       # if so, all eigenvalues should be positive
        evals <- eigen(vcov.mx, TRUE)$values 
        if(!any(evals < 0)) {
            done <- TRUE
        } else {}
        ntries <- ntries + 1
    }
    return(parms[1,])
}


#' Form a symmetric variance-covariance matrix
#'
#' Takes a vector of variances and correlations, calculates
#' covariances and forms a symmetric matrix (for use in generating
#' data from a multivariate normal distribution)
#'
#' @param vars vector of variances
#' @param corrs vector of correlations
#' @return A symmetric variance-covariance matrix
#' @examples
#' pars <- findValidVCovParams(vars=list(A=c(1,3),B=c(1,3), C=c(0,2)),
#'                             corrs=list(AB=c(-.5,.5),AC=c(-.5,.5),BC=c(-.5,.5)))
#'
#' formVCovMx(pars[1:3], pars[4:6])
#' @export formVCovMx
formVCovMx <- function(vars, corrs) {
    if (length(corrs)!=sum(1:(length(vars)-1))) {
        stop("need ", sum(1:(length(vars)-1)), " correlation parameters for ",
             length(corrs), " variances")
    } else {}
    if (is.list(vars)) {
        vv <- as.numeric(vars)
        if (is.null(names(vars))) {
            names(vv) <- paste("V", 1:length(vars), sep="")
        } else {
            names(vv) <- names(vars)
        }
        vars <- vv
    } else {}
    if (is.list(corrs)) {
        cc <- as.numeric(corrs)
        if (is.null(names(corrs))) {
            names(cc) <- paste("C", 1:length(corrs), sep="")
        } else {
            names(cc) <- names(corrs)
        }
        corrs <- cc
    } else {}    
    mx <- matrix(nrow=length(vars), ncol=length(vars))
    k <- 1
    for (i in 1:length(vars)) {
        for (j in i:length(vars)) {
            if (i==j) { # is a variance
                mx[i,j] <- vars[i]
            } else {
                mx[i,j] <- sqrt(vars[i])*sqrt(vars[j])*corrs[k]
                mx[j,i] <- sqrt(vars[i])*sqrt(vars[j])*corrs[k]
                k <- k+1
            }
        }
    }
    return(mx)
}

#' Load in package simgen on all available clusters
#'
#' This function is called by \code{\link{mcRun}} to load the package
#' functions across all clusters.  Typically users should not need to
#' call this function.
#'
#' @param cl cluster object (created e.g. by
#' \code{\link{makeCluster}}; see package parallel for details.
#' 
#' @export initializeCluster
initializeCluster <- function(cl) {
    invisible(parallel::clusterCall(cl, function(x) {library(simgen)}))
}
