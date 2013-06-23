#' Data-generating parameter ranges
#'
#' Generates parameter ranges for one-factor design.  These parameters
#' are used for generating datasets.
#'
#' @return A list containing ranges for the following parameters, default values in parentheses:
#' \item{int}{intercept (-3 to 3)}
#' \item{eff}{the effect of the IV (.8)}
#' \item{err}{error variance (1 to 3)}
#' \item{miss}{proportion missing data (0 to 5)}
#' \item{pMin}{lower bound on rate of missing data (0)}
#' \item{pMax}{upper bound on rate of missing data (.8)}
#' \item{t00}{by-subject intercept variance (0 to 3)}
#' \item{t11}{by-subject slope variance (0 to 3)}
#' \item{rsub}{by-subject random intercept/slope correlation (-.8 to .8)}
#' \item{w00}{by-item intercept variance (0 to 3)}
#' \item{w11}{by-item slope variance (0 to 3)}
#' \item{ritm}{by-item random intercept/slope correlation (-.8 to .8)}
#' @seealso \code{\link{genParamRanges.facMixedAB}}, \code{\link{randParams}}, \code{\link{mkDf}}
#' 
#' @export genParamRanges
genParamRanges <- function() {
    list(
        int=c(-3,3),       # range of intercept value, continuous simulations
        eff=.8,            # set to 0 to test Type I error rate
        err=c(0,3),         # range for error variance
        miss=c(0,.05),     # proportion of missing data
        pMin=0,                # lower bound on condition/cluster-level rate of  missing data
        pMax=0.8,              # lower bound on condition/cluster-level rate of  missing data
        t00=c(0, 3),         # tau_00 is the subject variance for the intercept
        t11=c(0, 3),         # tau_11 is the subject variance for the slope
        rsub=c(-.8, .8), # range of the by-subject intercept/slope correlation
        w00=c(0, 3),         # by-item intercept variance
        w11=c(0, 3),         # by-item slope variance
        ritm=c(-.8, .8) # by-item intercept/slope correlation
        )
}


#' Data-generating parameter ranges
#'
#' Generates parameter ranges for two-way factorial design, with A
#' between-subjects and B within-subjects, and A,B within-items.
#' 
#' @return A list containing ranges for the following parameters, default values in parentheses:
#' 
#' \item{t_00}{by-subject random intercept variance (1 to 3)}
#' \item{t_11}{by-subject random slope variance (1 to 3)}
#' \item{r_01}{by-subject random intercept/slope correlation (-.9,.9)}
#' \item{evar}{error variance (fixed at 3)}
#' \item{mu}{intercept}
#' \item{A}{effect of A}
#' \item{B}{effect of B}
#' \item{AB}{effect of AB}
#' \item{w_00}{by-item random intercept (1 to 3)}
#' \item{w_11}{by-item random slope variance for A (1 to 3)}
#' \item{w_22}{by-item random slope variance for B (1 to 3)}
#' \item{w_33}{by-item random slope variance for AB (1 to 3)}
#' \item{w_01}{by-item random intercept/A-slope variance (-.9,.9)}
#' \item{w_02}{by-item random intercept/B-slope variance (-.9,.9)}
#' \item{w_03}{by-item random intercept/AB-slope variance (-.9,.9)}
#' \item{w_12}{by-item random A/B-slope variance (-.9,.9)}
#' \item{w_13}{by-item random A/AB-slope variance (-.9,.9)}
#' \item{w_23}{by-item random B/AB-slope variance (-.9,.9)}
#' 
#' @seealso \code{\link{genParamRanges.facMixedAB}}, \code{\link{randParams}}, \code{\link{mkDf}}
#'
#' @export genParamRanges.facMixedAB
genParamRanges.facMixedAB <- function() {
    params1 <- list(t_00=c(0,3), # random intercept variance
                    t_11=c(0,3), # random slope variance
                    r_01=c(-.9,.9), # random intercept/slope correlation
                    evar=3, # error variance
                    mu=list(c(-2,-1), c(1,2)),
                    A=list(c(-2,-1), c(1,2)),
                    B=list(c(-2,-1), c(1,2)),
                    AB=list(c(-2,-1), c(1,2)))
    params2 <- list(w_00=c(0,3), # by-item int variance
                    w_11=c(0,3), # by-item slope var (A)
                    w_22=c(0,3), # by-item slope var (B)
                    w_33=c(0,3), # by-item slope var (AB)
                    w_01=c(-.9,.9), # by-item int/A corr
                    w_02=c(-.9,.9), # by-item int/B corr
                    w_03=c(-.9,.9), # by-item int/AB corr
                    w_12=c(-.9,.9), # by-item A/B corr
                    w_13=c(-.9,.9), # by-item A/AB corr
                    w_23=c(-.9,.9)) # by-item B/AB corr
    c(params1, params2)
}

#' Make matrix of data-generating parameters for mixed factorial design
#'
#' @param param.list list of DGP ranges (e.g., from \code{\link{genParamRanges.facMixedAB}})
#' @param nmc number of populations to generate
#' @param firstseed random number seed to use before generating populations
#' @return A matrix, each row of which corresponds to parameters for a single population
#' @seealso \code{\link{genParamRanges.facMixedAB}}, \code{\link{mkDf.facMixedAB}}
#' @examples
#' mkParamMx.facMixedAB(genParamRanges.facMixedAB(), 10)
#' @export mkParamMx.facMixedAB
mkParamMx.facMixedAB <- function(param.list, nmc=1000, firstseed=NULL) {
    if (!is.null(firstseed)) {
        set.seed(firstseed)
    } else {}
    item.ix <- grep("^w_", names(param.list))
    param.noitem <- param.list[-item.ix]
    param.item <- param.list[item.ix]
    mx <- randParams(param.noitem, nmc)
    item.mx <- t(replicate(nmc,
                           findValidVCovParams(param.item[c("w_00", "w_11", "w_22", "w_33")],
                                               param.item[c("w_01", "w_02", "w_03", "w_12", "w_13", "w_23")]),
                           simplify="matrix"))
    colnames(item.mx) <- names(param.item)
    seed.mx <- mkSeeds(nmc)
    cbind(mx, item.mx, seed=seed.mx)
}

#' Generate simulated data from a two-way mixed factorial design
#'
#' Generates a data set corresponding to a 2-way fixed factorial design.
#'
#' If \code{nitem} is defined, then it will generate item effects as
#' well as subject effects.  If \code{nreps} is defined, will generate
#' subject effects only.
#'
#' @param mcr.params Data-generating parameters
#' @param nsubj Number of subjects (must be a multiple of 4)
#' @param nitem Number of items (must be a multiple of 2, or NULL if \code{nreps} is defined)
#' @param nreps Number of replicated observations per cell per subject (must be NULL if \code{nitem} is defined)
#' @param showRfx Whether to show the random effects in the output
#' @return A data frame
#' @seealso \code{\link{genParamRanges.facMixedAB}}, \code{\link{randParams}}, \code{\link{mkParamMx.facMixedAB}}
#' @examples
#' mkDf.facMixedAB(nitem=24)
#' mkDf.facMixedAB(nreps=2)
#' 
#' #' @export mkDf.facMixedAB
mkDf.facMixedAB <- function(mcr.params=mkParamMx.facMixedAB(genParamRanges.facMixedAB(), 1)[1,],
                            nsubj=24, nitem=NULL, nreps=NULL, showRfx=FALSE) {
    if (is.null(nitem) && is.null(nreps)) {
        stop("need to define either 'nitem' or 'nreps' argument. see ?mkDf.facMixedAB")
    } else {}
    if (!is.null(nitem) && !is.null(nreps)) {
        stop("both 'nitem' and 'nreps' were defined; only one allowed. see ?mkDf.facMixedAB")
    } else {}
    if (!is.null(nitem)) {
        if ((nitem %% 2) != 0) {
            stop("nitem must be a multiple of 2")
        } else {}
    } else {}
    if ((nsubj %% 4) != 0) {
        stop("nsubj must be a multiple of 4")
    } else {}
    # mcr.params is a vector of parameters for data generation
    x <- mcr.params
    if (!is.list(x)) {
        x <- as.list(x)
    } else {}
    set.seed(x[["seed"]])
  
    library(MASS)
    ranef.subj <- mvrnorm(nsubj, mu=c(0,0),
                          Sigma=formVCovMx(x[c("t_00","t_11")], x["r_01"]))
    if (is.null(nitem)) {
        srs <- rep(ranef.subj[,2], each=2)*c(-.5,.5) # subject random slope
        dat <- data.frame(SubjID=factor(rep(1:nsubj,each=2*nreps)),
                          A=factor(rep(c("A1","A2"),each=2*nreps)),
                          B=factor(rep(c("B1","B2"),each=nreps)))
        dat$s.i <- rep(ranef.subj[,1], each=2*nreps)
        dat$s.s <- rep(srs, each=nreps)
        dat$err <- rnorm(nrow(dat), sd=sqrt(x[["evar"]]))
        dat$Y <- x[["mu"]] + # fixed intercept
            rep(c(-x[["A"]],x[["A"]])/2, each=2*nreps) + # fixed effect of A
                rep(c(-x[["B"]],x[["B"]])/2, each=nreps) +   # fixed effect of B
                    rep(c(x[["AB"]],-x[["AB"]],-x[["AB"]],x[["AB"]])/4, each=nreps) + # fixed interaction
                        dat$s.i + # random intercept
                            dat$s.s + # random slope
                                dat$err # error variance
        dat$At <- 1*(dat$A=="A2")
        dat$Bt <- 1*(dat$B=="B2")
        dat$Ad <- dat$At-mean(dat$At)
        dat$Bd <- dat$Bt-mean(dat$Bt)
        if (!showRfx) {
            dat <- dat[,setdiff(colnames(dat), c("s.i","s.s","err"))]
        } else {}
    } else { # data structure with item effects
        subj <- data.frame(SubjID=factor(1:nsubj), ListID=factor(1:2),
                           A=factor(rep(c("A1","A2"),each=nsubj/2)),
                           s.i=ranef.subj[,1], s.B=ranef.subj[,2])
        list.info <- data.frame(ListID=factor(rep(1:2,each=nitem)),
                                ItemID=factor(1:nitem),
                                B=factor(rep(c("B1","B2","B2","B1"),each=nitem/2)))
        ranef.item <- cbind(ItemID=factor(1:nitem),
                            as.data.frame(mvrnorm(nitem, mu=rep(0,4),
                                                  Sigma=formVCovMx(x[c("w_00","w_11","w_22","w_33")],
                                                      x[c("w_01","w_02","w_03","w_12","w_13","w_23")]))))
        colnames(ranef.item) <- c("ItemID","i.i","i.A","i.B","i.AB")
        dat <- merge(subj, merge(list.info, ranef.item, by="ItemID"))
        dat$At <- 1*(dat$A=="A2")
        dat$Bt <- 1*(dat$B=="B2")
        dat$Ad <- dat$At-mean(dat$At)
        dat$Bd <- dat$Bt-mean(dat$Bt)
        front.cols <- c("SubjID","ListID","ItemID","A","B","At","Bt","Ad","Bd")
        dat <- dat[order(dat$SubjID, dat$ItemID), c(front.cols, setdiff(colnames(dat), front.cols))]
        rownames(dat) <- NULL
        dat$rfx <- with(dat, s.i+Bd*s.B + i.i+Ad*i.A+Bd*i.B+Ad*Bd*i.AB)
        dat$ffx <- with(dat, x[["mu"]] + x[["A"]]*Ad + x[["B"]]*Bd + x[["AB"]]*Ad*Bd)
        dat$Y <- dat$ffx + dat$rfx + rnorm(nrow(dat), sd=sqrt(x[["evar"]]))
        dat <- dat[,c(front.cols, "Y")]
    }
  
    return(dat)
}
