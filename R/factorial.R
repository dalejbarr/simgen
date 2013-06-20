#' Data-generating parameter ranges
#'
#' Generates parameter ranges for one-factor design.  These parameters
#' are used for generating datasets.
#'
#' @return A list containing ranges for the following parameters, default values in parentheses:
#' @returnItem int intercept (-3 to 3)
#' @returnItem eff the effect of the IV (.8)
#' @returnItem err error variance (1 to 3)
#' @returnItem miss proportion missing data (0 to 5)
#' @returnItem pMin lower bound on rate of missing data (0)
#' @returnItem pMax upper bound on rate of missing data (.8)
#' @returnItem t00 by-subject intercept variance (0 to 3)
#' @returnItem t11 by-subject slope variance (0 to 3)
#' @returnItem rsub by-subject random intercept/slope correlation (-.8 to .8)
#' @returnItem w00 by-item intercept variance (0 to 3)
#' @returnItem w11 by-item slope variance (0 to 3)
#' @returnItem ritm by-item random intercept/slope correlation (-.8 to .8)
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

#' @export genParamRanges.facMixedAB
genParamRanges.facMixedAB <- function() {
    params1 <- list(t_00=c(1,3), # random intercept variance
                    t_11=c(1,3), # random slope variance
                    r_01=c(-.9,.9), # random intercept/slope correlation
                    evar=3, # error variance
                    mu=list(c(-2,-1), c(1,2)),
                    A=list(c(-2,-1), c(1,2)),
                    B=list(c(-2,-1), c(1,2)),
                    AB=list(c(-2,-1), c(1,2)))
    params2 <- list(w_00=c(1,3), # by-item int variance
                    w_11=c(1,3), # by-item slope var (A)
                    w_22=c(1,3), # by-item slope var (B)
                    w_33=c(1,3), # by-item slope var (AB)
                    w_01=c(-.9,.9), # by-item int/A corr
                    w_02=c(-.9,.9), # by-item int/B corr
                    w_03=c(-.9,.9), # by-item int/AB corr
                    w_12=c(-.9,.9), # by-item A/B corr
                    w_13=c(-.9,.9), # by-item A/AB corr
                    w_23=c(-.9,.9)) # by-item B/AB corr
    c(params1, params2)
}

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


#' @export fitlmer.facMixedAB
fitlmer.facMixedAB <- function(mcr.data, treat=TRUE, items=FALSE) {
    # models to fit
    if (!items) {
        models <- list(noB1=Y~A*B+(1|SubjID),
                       noB2=Y~A+B+(1|SubjID),
                       max1=Y~A*B+(1+B|SubjID),
                       max2=Y~A+B+(1+B|SubjID))
    } else {
        models <- list(noB1=Y~A*B+(1|SubjID)+(1|ItemID),
                       noB2=Y~A+B+(1|SubjID),
                       max1=Y~A*B+(1+B|SubjID),
                       max2=Y~A+B+(1+B|SubjID))
    }
    modCompare <- list(noB.t=list("noB1","noB2"), max.t=list("max1","max2"))
    # with treatment coding
    if (treat) {
        contrasts(mcr.data$A) <- contr.treatment(levels(mcr.data$A))
        contrasts(mcr.data$B) <- contr.treatment(levels(mcr.data$B))
        treat.p <- unlist(lrCompare(mcr.data, models, modCompare,
                                    REML=FALSE, na.option=na.omit))
        names(treat.p) <- names(modCompare)
    } else {}
    # with deviation coding
    contrasts(mcr.data$A) <- contr.deviation(levels(mcr.data$A))
    contrasts(mcr.data$B) <- contr.deviation(levels(mcr.data$B))
    names(modCompare) <- c("noB.d", "max.d")
    dev.p <- unlist(lrCompare(mcr.data, models, modCompare,
                              REML=FALSE, na.option=na.omit))
    names(dev.p) <- names(modCompare)
    dat.aov <- summary(aov(Y~A*B + Error(SubjID/B), mcr.data))
    if (treat) {
        res <- c(treat.p, dev.p,
                 aov=dat.aov[["Error: SubjID:B"]][[1]]$`Pr(>F)`[2])
    } else {
        res <- c(dev.p,
                 aov=dat.aov[["Error: SubjID:B"]][[1]]$`Pr(>F)`[2])
    }
    return(res)
}


#' @export mkDf.facMixedAB
mkDf.facMixedAB <- function(mcr.params=mkParamMx(defaultGenParamRanges.facMixedAB(), 1)[1,],
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
