#' Monte Carlo Simulation of Simultaneous Generalization to Subjects and Items
#' 
#' Routines for running Monte Carlo simulations testing the generalizability of
#' various ways of analyzing datasets with crossed random effects.  This
#' accompanies a paper by Barr, Levy, Scheepers, & Tily (under revision).  See
#' http://idiom.ucsd.edu/~rlevy/papers/barr-levy-scheepers-tily-8-aug-2011-submitted.pdf.
#' 
#' \tabular{ll}{ Package: \tab simgen\cr Type: \tab Package\cr Version: \tab
#' 1.5\cr Date: \tab 2012-02-20\cr License: \tab GPL (>= 2)\cr LazyLoad: \tab
#' yes\cr } To get started, see \code{\link{createParamMx}},
#' \code{\link{mkDf}}, and \code{\link{mcRun}}.
#' 
#' @name simgen-package
#' @aliases simgen-package simgen
#' @docType package
#' @author Dale J. Barr <dale.barr@@glasgow.ac.uk>
#' @keywords package



#' Fit model with error-trapping
#' 
#' Try to fit an lmer model, catching any errors/warnings
#' 
#' 
#' @param mf The model formula
#' @param xdat The data frame
#' @return %% ~Describe the value returned %% If it is a LIST, use
#' 
#' %% ...
#' @return
#' \item{value} fitted model object from lmer
#' \item{converged} whether or not the model converged
#' @seealso \code{\link{fitlmer}}, \code{\link{modSpace}}
#' @examples
#' 
#' # random intercepts only
#' ff <- tryFit(Resp~Cond + (1+Cond|SubjID)+(1+Cond|ItemID), mkDf())
#'
#' @importFrom lme4 lmer
#' @export tryFit
tryFit <- function(tf.formula, tf.data, ...) {
    converged <- TRUE
    w.handler <- function(w) {
        converged <- FALSE
        invokeRestart("muffleWarning")
    }
    arg.list <- c(list(formula=tf.formula, data=tf.data), list(...))
    list(value=withCallingHandlers(tryCatch(
             do.call(lme4:::lmer, arg.list),
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


#' Matrix for Deviation-Coded Predictors
#'
#' Return a matrix of deviation-coded predictors (analogous to
#' \code{\link{contr.sum}}).
#'
#' @param n A vector of levels for a factor, or the number of levels.
#' @export contr.deviation
contr.deviation <- function(n) {
    if (length(n)<=1L) {
        if (is.numeric(n) && length(n) == 1L && n > 1L)
            levels <- seq_len(n)
        else stop("not enough degrees of freedom to define contrasts")
    } else {}
    apply(contr.treatment(n), 2, function(x) {x-mean(x)})
}



#' Get p-value from model comparison of lmer objects
#'
#' Derive p-value by comparing two models fitted by
#' \code{\link{tryFit}}.
#'
#' @param m1 A model fitted by \code{\link{tryFit}}
#' @param m2 A model fitted by \code{\link{tryFit}}
#' @return A p-value or \code{NA} if either (or both) models did not converge
#' @examples
#' d1 <- mkDf()
#' 
#' m1 <- tryFit(Resp ~ Cond + (1 + Cond | SubjID) + (1 + Cond | ItemID), d1, REML=FALSE,
#'               na.action=na.omit)
#' m2 <- tryFit(Resp ~ (1 + Cond | SubjID) + (1 + Cond | ItemID), d1, REML=FALSE,
#'               na.action=na.omit)
#'
#' getLmer.pValue(m1, m2)
#' 
#' @importFrom lme4 lmer
#' @export getLmer.pValue
getLmer.pValue <- function(m1,m2) {
    pval <- NA
    if (m1$converged && m2$converged) {
        chisq.val <- (-2*lme4::logLik(m1$value))-(-2*lme4::logLik(m2$value))
        pval <- as.numeric(pchisq(abs(chisq.val), 1, lower.tail=FALSE))
    } else {}
    return(pval)
}

#' Fit mixed-models and obtain p-values from model comparison
#'
#' A set of model formula are fed to \code{lrCompare} as named
#' elements in a list.  These models are fit using
#' \code{\link{tryFit}}.  Models are compared as specified in the
#' argument \code{lrc.modCompare} using LLR (likelihood-ratio) tests.
#' This function is useful for performing multiple tests on relatively
#' complex designs (e.g., all main effects and interactions in a
#' factorial design).
#'
#' @param mcr.data Data to be analyzed
#' @param lrc.mods List of (named) model formulae
#' @param lrc.modCompare List of comparisons to perform (named)
#' @param ... Arguments to be passed along to \code{\link{tryFit}}
#' (and thus to \code{lmer}).
#' @return A vector with the p-values from tests specified by
#' \code{lrc.modCompare}. If a test could not be performed (e.g.,
#' because one or both of the corresponding models did not converge),
#' \code{NA} is returned.
#' 
#' @export lrCompare
lrCompare <- function(mcr.data, lrc.mods, lrc.modCompare, ...) {
    extraArgs <- c(list(tf.data=mcr.data), list(...))
    modres <- do.call(lapply, c(list(X=lrc.mods, FUN=tryFit), extraArgs))
    unlist(lapply(lrc.modCompare, function(x) {
        getLmer.pValue(modres[[x[[1]]]], modres[[x[[2]]]])
    }))
}


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
    ntries <- 0
    done <- FALSE
    # use brute force to find a definite positive matrix
    while (!done && (ntries < maxtries)) {
        parms <- randParams(c(vars, corrs))
        vcov.mx <- formVCovMx(parms[1:length(vars)], parms[-(1:length(vars))])
       # test whether the matrix is positive-definite
       # if so, all eigenvalues should be positive
        evals <- eigen(vcov.mx, TRUE)$values 
        if(sum(evals>0)==dim(vcov.mx)[1]) {
            done <- TRUE
        } else {}
        ntries <- ntries + 1
    }
    return(parms)
}

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


#' @importFrom parallel clusterCall
#' @export initializeCluster
initializeCluster <- function(cl) {
    invisible(clusterCall(cl, function(x) {library(simgen)}))
}
