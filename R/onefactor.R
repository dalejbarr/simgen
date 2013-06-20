.reassembleStepwise <- function(mx, crit) {
  forward <- grep("forward", colnames(mx))
  backward <- grep("backward", colnames(mx))
  tlist <- list()
  if (length(forward)>0) {
    ff.forward <- aperm(array(c(as.matrix(mx)), dim=c(nrow(mx),ncol(mx[,forward])/6, 6)), c(1,2,3))
    dimnames(ff.forward) <- list(1:nrow(mx), c(.01,.05,seq(.1,.8,.1)),
                                 c("fm","t","chi","pt","pchi","esteff"))
    tlist <- list(forward=ff.forward)
  } else {}
  if (length(backward)>0) {
    ff.backward <- aperm(array(c(as.matrix(mx)), dim=c(nrow(mx),ncol(mx[,backward])/6, 6)), c(1,2,3))
    dimnames(ff.backward) <- list(1:nrow(mx), c(.01,.05,seq(.1,.8,.1)),
                                  c("fm","t","chi","pt","pchi","esteff"))
    tlist[["backward"]] <- ff.backward
  } else {}
  if (length(tlist)==1) {
    result <- tlist[[1]]
  } else {
    result <- tlist
  }
  return(result)
}

#' @importFrom lme4 VarCorr
.modCompare <- function(mods) {
  getDf <- function(mm) {
    1+sum(c(unlist(lapply(lme4:::VarCorr(mm), function(y) {sum(1:dim(y)[1])})), length(fixef(mm))))
  }
  mod.inf <- rbind(unlist(lapply(mods, deviance)),
                   unlist(lapply(mods, getDf)))
                                        # compare all models sequentially
  ff <- lapply(2:ncol(mod.inf), function(x) {
    chi.obs=abs(mod.inf[1,x-1]-mod.inf[1,x])
    df=abs(mod.inf[2,x-1]-mod.inf[2,x])
    res=c(chi.obs=chi.obs, df=df, p=pchisq(chi.obs, df, lower.tail=FALSE))
    names(res) <- c("chi.obs", "df","p")
    return(res)
  })
  cmp.mx <- matrix(unlist(ff), nrow=3)
  rownames(cmp.mx) <- names(ff[[1]])
  colnames(cmp.mx) <- lapply(2:ncol(mod.inf), function(x) {paste(colnames(mod.inf)[(x-1):x], collapse=".")})
  return(cmp.mx)
}

.stepwiseFit <- function(mf, mcr.data, crit=c(.01,.05,seq(.1,.8,.1))) {
                                        # mf: list of model formulae in order of testing
                                        # xd: the data frame
                                        # forward: in forward selection (TRUE), testing proceeds
                                        #   while p<crit; in backwards, testing proceeds while p>crit
                                        # crit: critical value for taking a step
  xd <- mcr.data
  if (missing(mf)) {
    stop("missing argument 'mf' for stepwise fitting function (see ?fitstepwise)")
  } else {}
  
  compareModels <- function(cmp.mx, crit, forward) {
    if (forward) {
      mod.sig <- lapply(crit, function(x) {
        ff <- (1:ncol(cmp.mx))[cmp.mx["p",]>x] # stop when p > crit
        if (length(ff)==0) {
          ff1 <- ncol(cmp.mx)+1
        } else {
          ff1 <- min(ff)
        }
        return(ff1)
      })
      mod.win.ix <- cvg.ix[unlist(lapply(mod.sig, max))]
    } else { # backward
      mod.sig <- lapply(crit, function(x) {
        ff <- (1:ncol(cmp.mx))[cmp.mx["p",]<=x] # stop when p <= crit
        if (length(ff)==0) {
          ff1 <- 1
        } else {
          ff1 <- max(ff)+1
        }
        return(ff1)
      })
      mod.win.ix <- cvg.ix[unlist(lapply(mod.sig, min))]
    }
    mod.win <- lapply(mod.win.ix, function(x) {
      ff <- list(mods[[x]])
      names(ff) <- names(mfits)[cvg.ix[x]]
      return(ff)
    })
    names(mod.win) <- crit
    return(mod.win)    
  }
  tot.mods <- length(mf) # total models to test
  if (tot.mods < 2) {
    stop("need 2 or more models to be tested")
  } else {}
  mfits <- lapply(mf, tryFit.default, xdat=xd)
  cvg.ix <- (1:length(mfits))[unlist(lapply(mfits, function(x) {x$converged}))]
  mods <- lapply(mfits[cvg.ix], function(x) {x$value})
  if (length(mods) < 2) {
    if (length(mods) == 0) {
      mod.win <- lapply(crit, function(x) {
        ff <- list(NULL)
        names(ff) <- "null"
        return(ff)
      })
    } else {
      mod.win <- lapply(crit, function(x) {
        ff <- list(mods[[1]])
        names(ff) <- names(mfits)[cvg.ix]
        return(ff)
      })
    }
    names(mod.win) <- crit
    result <- list(forward=mod.win, backward=mod.win, cmp.mx=NULL)
  } else {
    result <- list(forward=NULL, backward=NULL, cmp.mx=NULL)
    cmp.mx <- .modCompare(mods)
    result$forward <- compareModels(cmp.mx, crit, forward=TRUE)
    result$backward <- compareModels(cmp.mx, crit, forward=FALSE)
    result$cmp.mx <- cmp.mx
  }
  return(result)
}

.mod2Code <- function(mod) {
  mf <- .modSpace(FALSE)
  mf2 <- c(names(c(mf[1], mf[[2]][1:2], mf[3])), "nofit")
  gg <- (1:5)[names(mod)==mf2]
  if (length(gg)==0) {
    mf2 <- c(names(c(.modSpace(FALSE), recursive=TRUE)), "nofit")
    gg <- (1:5)[names(mod)==mf2]
  } else {}
  return(gg)
}

#' @importFrom Matrix diag
#' @importFrom lme4 vcov
.modInfo <- function(mod, xd, wsbi) {
                                        # do comparison model
  mf2 <- .modSpace(wsbi)[[names(mod)]]
  if (is.null(mf2)) {
    mf2 <- c(.modSpace(wsbi), recursive=TRUE)[[names(mod)]]
    if (is.null(mf2)) {
      mf2 <- .modSpace(wsbi)[[2]][[names(mod)]]
    } else {}    
  } else {}
  mf3 <- as.formula(paste(mf2[2], "~", mf2[3], "-Cond", sep=""))
  if ( (mod2 <- tryFit.default(mf3, xd))$converged ) {
    m1 <- mod2[["value"]]; m2 <- mod[[1]]
    mcmp <- anova(m1, m2)
    chi.obs <- mcmp$`Chisq`[2]
    chi.p <- mcmp$`Pr(>Chisq)`[2]
  } else {
    chi.obs <- NA
    chi.p <- NA
  }
  t.obs <- (fixef(mod[[1]])[2])/(sqrt(Matrix:::diag(lme4:::vcov(mod[[1]])))[2])
  v1 <- c(fm=.mod2Code(mod), t.obs[1],
          chi.obs,
          2*(1-pnorm(abs(t.obs))),
          chi.p,
          eff=fixef(mod[[1]])[2])
  names(v1) <- c("fm", "t","chi","pt","pchi","esteff")
  return(v1)
}

#' @importFrom Matrix diag
#' @importFrom lme4 vcov
.modInfo1 <- function(mf2, xd, mod) {
  mf3 <- as.formula(paste(mf2[2], "~", mf2[3], "-Cond", sep=""))
  if ( (mod2 <- tryFit.default(mf3, xd))$converged ) {
    m1 <- mod2[["value"]]; m2 <- mod[[1]]
    mcmp <- anova(m1, m2)
    chi.obs <- mcmp$`Chisq`[2]
    chi.p <- mcmp$`Pr(>Chisq)`[2]
  } else {
    chi.obs <- NA
    chi.p <- NA
  }
  t.obs <- (fixef(mod[[1]])[2])/(sqrt(Matrix:::diag(lme4:::vcov(mod[[1]])))[2])
  v1 <- c(fm=mod$converged, t.obs[1],
          chi.obs,
          2*(1-pnorm(abs(t.obs))),
          chi.p,
          eff=fixef(mod[[1]])[2])
  names(v1) <- c("fm", "t","chi","pt","pchi","esteff")
  return(v1)
}

.multiModInfo <- function(mfits, xd, wsbi) {
  allMods <- unlist(lapply(mfits, function(x) {names(x)}))
  modTypes <- unique(allMods)
  names(modTypes) <- modTypes
  mtFits <- lapply(modTypes, function(x) {mfits[[min((1:length(allMods))[allMods==x])]]})
  mtInfo <- lapply(mtFits, .modInfo, xd=xd, wsbi=wsbi)
  ff <- mtInfo[allMods]
  ff.mx <- matrix(unlist(ff), ncol=length(ff[[1]]), byrow=TRUE,
                  dimnames=list(names(mfits), names(ff[[1]])))
  return(ff.mx)
}

.modSpace <- function(wsbi) {
  if (wsbi) {
    mf <- list(max=Resp ~ Cond + (1 + Cond | SubjID) + (1 | ItemID),
               min=Resp ~ Cond + (1 | SubjID) + (1 | ItemID))
  } else {
    mf <- list(max=Resp ~ Cond + (1 + Cond | SubjID) + (1 + Cond | ItemID),
               mid=list(srs=Resp ~ Cond + (1 + Cond | SubjID) + (1 | ItemID),
                 irs=Resp ~ Cond + (1 | SubjID) + (1 + Cond | ItemID)),
               min=Resp ~ Cond + (1 | SubjID) + (1 | ItemID))
  }
  return(mf)
}

.tryFit.default <- function(mf, xdat) {
    tryFit(mf, xdat, na.action=na.omit, REML=FALSE)
}

#' Create Population Parameter Matrix for One-Factor Design
#' 
#' Create a matrix of population parameters for Monte Carlo simulation
#' with one-factor design.
#' 
#' Each row of the matrix returned from this function represents the
#' population parameters for a given experiment.  As of version 1.6 of
#' simgen, This function has been superceded by
#' \code{\link{genParamRanges}} and \code{\link{randParams}}, but is
#' preserved for backwards compatibility.
#' 
#' @param nexp Number of experiments to run (default 100000).
#' @param simparam.env An \code{\link{environment}} object containing the
#' parameter ranges; see \code{\link{getParamRanges}} (default) for the
#' appropriate format.
#' @param firstseed First seed to start off generation of the matrix.  Included
#' for replicability of results.  A parameter matrix of a given size with a
#' given seed will always be identical.
#' @param h0 Status of the null hypothesis: TRUE or FALSE.
#' @param outfile Name of a file where the matrix will be stored; if NULL
#' (default), no data will be stored.
#' @return A matrix, with the following columns: %% ~Describe the value
#' returned %% If it is a LIST, use
#' 
#' %% ...
#' @returnItem int (grand mean) intercept
#' @returnItem eff treatment effect
#' @returnItem err error variance
#' @returnItem miss proportion of response values missing
#' @returnItem t00 by-subject intercept variance
#' @returnItem t11 by-subject slope variance
#' @returnItem rsub by-subject intercept/slope correlation
#' @returnItem w00 by-item intercept variance
#' @returnItem w11 by-item slope variance
#' @returnItem ritm by-item intercept/slope correlation
#' @returnItem seed random number seed for creating the dataframe (see
#' \code{\link{mkDf}})
#' @seealso \code{\link{getParamRanges}}, \code{\link{mkDf}}
#' @examples
#' 
#' # using defaults
#' createParamMx(10)
#' 
#' # now let's change one of the ranges
#' p.env <- getParamRanges()
#' get("t11.range", env=p.env)
#' assign("t11.range", c(5,10), env=p.env)
#' createParamMx(10, simparam.env=p.env)
#' 
#' @export createParamMx
createParamMx <- function(nexp=100000,    # no. of simulated experiments for which to generate population parameters
                          simparam.env=getParamRanges(),
                                          # environment containing ranges of population params (usually param.env)
                          firstseed=666,  # seed to start things off (this function generates further seeds)
                          h0=TRUE,        # null hypothesis TRUE or FALSE
                          outfile=NULL) { # name of file to write the matrix to
  #   Generates a matrix of parameters.  Each row are the population parameters
  #   used to generate data from a single hypothetical "experiment."
  warning("this function is only available for backwards compatibility with simgen_1.54; please consider using genParamRanges and randParams instead")
  if (!is.null(firstseed)) {
    set.seed(firstseed)
  } else {}
  param.mx <- matrix(nrow=nexp, ncol=13)
  colnames(param.mx) <- c("int", "eff", "err","miss","pMin","pMax","t00","t11","rsub","w00","w11","ritm","seed")
  param.mx[, "int"] <- runif(nexp, min=simparam.env$icept.range[1],
                             max=simparam.env$icept.range[2])
  param.mx[, "eff"] <- ifelse(h0, simparam.env$slope["h0"], simparam.env$slope["h1"])
  param.mx[, "err"] <- runif(nexp, min=simparam.env$evar.range[1],
                             max=simparam.env$evar.range[2])
  param.mx[, "miss"] <- runif(nexp, min=simparam.env$pmissing.range[1],
                              max=simparam.env$pmissing.range[2])
  param.mx[, "pMin"] <- simparam.env$pMin
  param.mx[, "pMax"] <- simparam.env$pMax
  param.mx[, "t00"] <- runif(nexp, min=simparam.env$t00.range[1],
                        max=simparam.env$t00.range[2])
  param.mx[, "t11"] <- runif(nexp, min=simparam.env$t11.range[1],
                        max=simparam.env$t11.range[2])
  param.mx[, "rsub"] <- runif(nexp, min=simparam.env$r01.subj.range[1],
                        max=simparam.env$r01.subj.range[2])
  param.mx[, "w00"] <- runif(nexp, min=simparam.env$w00.range[1],
                        max=simparam.env$w00.range[2])
  param.mx[, "w11"] <- runif(nexp, min=simparam.env$w11.range[1],
                        max=simparam.env$w11.range[2])
  param.mx[, "ritm"] <- runif(nexp, min=simparam.env$r01.item.range[1],
                        max=simparam.env$r01.item.range[2])
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
  param.mx[, "seed"] <- randSeed(nexp)

  if (!is.null(outfile)) {
    save(param.mx, file=outfile)
    return(invisible(param.mx))
  } else {}
  
  return(param.mx)
}


#' Data-generating parameter ranges for one-factor design
#' 
#' Get the default list of parameter ranges used in the Monte Carlo
#' simulations for one-factor design.  These parameter ranges can then
#' be modified by the user.
#' 
#' @return An environment object, containing the following variables:
#' @returnItem evar.range range of error variance (default [0,3])
#' @returnItem icept.range range of intercept (default [-3, 3])
#' @returnItem pmissing.range range of proportion missing obs (default [0,
#' .05])
#' @returnItem r01.item.range range of by-subject intercept/slope correlation
#' (default [-.8, .8])
#' @returnItem r01.subj.range range of by-item intercept/slope correlation
#' (default [-.8, .8])
#' @returnItem slope treatment effect, when h0 TRUE, 0; when h0 FALSE, .8
#' @returnItem t00.range range of by-subject random intercept variance (default
#' [0, 3]
#' @returnItem t11.range range of by-subject random slope variance (default [0,
#' 3])
#' @returnItem w00.range range of by-item random intercept variance (default
#' [0, 3])
#' @returnItem w11.range range of by-item random slope variance (default [0,
#' 3])
#' @seealso \code{\link{createParamMx}}
#' @examples
#' 
#' # display all the parameter ranges
#' print(as.list(getParamRanges()))
#' 
#' # how to alter one of the values
#' p.env <- getParamRanges()
#' p.env$pmissing.range <- c(0, .25)
#' assign("pmissing.range", c(0, .25), env=p.env) # alternate way
#' print(as.list(p.env))
#' 
#' @export getParamRanges
getParamRanges <- function() {
  param.env <- new.env()
  with(param.env, {
    pmissing.range <- c(0,.05)  # proportion of missing data
    pMin <- 0                   # lower bound on condition/cluster-level rate of  missing data
    pMax <- 0.8                   # lower bound on condition/cluster-level rate of  missing data
    icept.range <- c(-3,3)       # range of intercept value, continuous simulations
    slope <- c(h0=0, h1=.8)      # the treatment effect (H0 true; H0 false    )
    evar.range <- c(0,3)         # range for error variance
    t00.range <- c(0, 3)         # tau_00 is the subject variance for the intercept
    t11.range <- c(0, 3)         # tau_11 is the subject variance for the slope
    r01.subj.range <- c(-.8, .8) # range of the by-subject intercept/slope correlation
    w00.range <- c(0, 3)         # by-item intercept variance
    w11.range <- c(0, 3)         # by-item slope variance
    r01.item.range <- c(-.8, .8) # by-item intercept/slope correlation
  })
  return(param.env)
}

#' Run ANOVA analysis
#' 
#' Given a dataset \code{xd}, calculates F1, F2, F'min, and F1+F2, and
#' associated p-values (one-factor design only).
#' 
#' 
#' @param mcr.data A dataframe formatted as described in \code{\link{mkDf}}.
#' @param wsbi Whether the design is between-items (TRUE) or within-items
#' (FALSE).
#' @return A vector, with:
#' @returnItem F1 the by-subjects F value
#' @returnItem F2 the by-items F value
#' @returnItem minF the min F' value
#' @returnItem p1 the p-value associated with F1
#' @returnItem p2 the p-value associated with F2
#' @returnItem pmf the p-value associated with min F'
#' @returnItem pmax the p-value associated with F1+F2 ( =max(c(p1,p2)) )
#' @seealso \code{\link{mkDf}}
#' @examples
#' 
#' nmc <- 10
#' pmx <- cbind(randParams(genParamRanges(), nmc, 1001), seed=mkSeeds(nmc, 1001))
#' 
#' # between-items dataset
#' x.bi <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=TRUE)
#' 
#' # within-items dataset
#' x.wi <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=FALSE)
#' 
#' fitanova(x.bi, wsbi=TRUE)
#' fitanova(x.wi, wsbi=FALSE)
#' 
#' @export fitanova
fitanova <- function(mcr.data, wsbi) {
  xd <- mcr.data
  xd$C <- factor(xd$Cond)
  xd$S <- factor(xd$SubjID)
  xd$I <- factor(xd$ItemID)
  F1.all <- summary(aov(Resp ~ C + Error(S/C), data=xd))[[2]][[1]]
  F1 <- F1.all$`F value`[1]
  F1.df.denom <- F1.all$`Df`[2]
  p1 <- F1.all$`Pr(>F)`[1]
  if (wsbi) {
    F2.all <- summary(aov(Resp ~ C + Error(I), data=xd))[[1]][[1]]
  } else {
    F2.all <- summary(aov(Resp ~ C + Error(I/C), data=xd))[[2]][[1]]
  }
  F2 <- F2.all$`F value`[1]
  F2.df.denom <- F2.all$`Df`[2]
  p2 <- F2.all$`Pr(>F)`[1]
  minF <- (F1*F2) / (F1 + F2)
  minF.df.denom <- ( (F1+F2)^2 )  /   ( (F1^2)/F2.df.denom + (F2^2)/F1.df.denom ) 
  v1 <- c(F1=F1, F2=F2, minF=minF,
          p1=p1, p2=p2,  pmf=pf(minF, 1, minF.df.denom, lower.tail=FALSE),
          pmax=max(c(p1,p2)))
  return(v1)
}



#' Run mixed-effects model (\code{\link{lmer}}) with MCMC p-value
#' 
#' Uses \code{\link{lmer}} to fit a random-intercepts only model, and then uses
#' \code{\link{mcmcsamp}} to obtain p-values (one-factor design only).
#' 
#' Tries to fit the most complex model requested, and if that fails to
#' converge, then tries progressively simpler models, before calculating mcmc
#' p-value.  If no model converges, returns NAs.  The MCMC procedure is based
#' on Baayen's \code{pvals.fnc} in package \code{languageR}.
#' 
#' @param mcr.data A dataframe formatted as described in \code{\link{mkDf}}.
#' @param wsbi Whether the design is between-items (TRUE) or within-items
#' (FALSE); has no effect because the model is random-intercepts only, but was
#' included for consistency with \code{\link{fitlmer}} and
#' \code{\link{fitanova}}.
#' @param nmcmc Number of Markov-Chain Monte Carlo simulations (default =
#' 10000).
#' @return A single element vector with the mcmc p-value.
#' @note This function ONLY fits random intercept models.
#' @seealso \code{\link{fitlmer}}
#' @examples
#' 
#' nmc <- 10
#' pmx <- cbind(randParams(genParamRanges(), nmc, 1001), seed=mkSeeds(nmc, 1001))
#' 
#' # between-items dataset
#' x.bi <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=TRUE)
#' 
#' fitlmer.mcmc(x.bi, wsbi=FALSE, nmcmc=1000) # maximal
#' 
#' @export fitlmer.mcmc
fitlmer.mcmc <- function(mcr.data, wsbi, nmcmc=10000) {
  xd <- mcr.data
  p.mcmc <- NA
  library(lme4, quietly=TRUE)
  mf <- Resp ~ Cond + (1 | SubjID) + (1 | ItemID)
  xd.lmer <- tryCatch(lmer(mf,
                           family=gaussian, data=xd, na.action=na.omit, REML=FALSE),
                      warning = function(w) {return (NULL)},
                      error = function(e) {return (NULL)})
  mcmc <- try(lme4::mcmcsamp(xd.lmer, n = nmcmc, silent=TRUE))
  if (!is(mcmc, "try-error")) {
    mcmcfixef <- t(mcmc@fixef)
    nr <- nrow(mcmcfixef)
    prop <- colSums(mcmcfixef > 0)/nr
    ans <- 2 * pmax(0.5/nr, pmin(prop, 1 - prop))
    names(ans) <- colnames(mcmcfixef)
    p.mcmc <- ans["Cond"]
  } else {}
  names(p.mcmc) <- "pmcmc"
  return(p.mcmc)
}



#' Run mixed-effects model using (\code{\link{lmer}})
#' 
#' Runs either a random-intercepts only model or a maximal
#' random-effects model (by-subject and by-item random intercepts,
#' by-subject random slope, and by-item random slope for within-item
#' design); NB: this function is for one-factor design only.
#' 
#' \code{fitlmer} will attempt to fit the model specified by the user, and will
#' progressively simplify the model as needed to get it to converge.  If no
#' model converges, it returns NAs.  \code{fitlmer} performs a likelihood ratio
#' test for the treatment effect, as well as returns a p-value for the
#' t-statistic using an approximation from the normal distribution.
#' 
#' @param mcr.data A dataframe formatted as described in \code{\link{mkDf}}.
#' @param ri.only Whether the random effects specification is to be
#' random-intercepts only (TRUE) or maximal random-effects (FALSE).
#' @param wsbi Whether the design is between-items (TRUE) or within-items
#' (FALSE).
#' @return
#' @returnItem fm Code for the model that converged: (1) dropped one slope; (2)
#' dropped two slopes; (3) main model did not converge; (4) comparision model
#' for likelihood ratio test did not converge.
#' @returnItem t t-statistic for the treatment effect
#' @returnItem chi chi-square statistic for the likelihood ratio test (1 df)
#' @returnItem pt p-value for the t-statistic (normal distribution)
#' @returnItem pchi p-value for the chi-square statistic
#' @seealso \code{\link{mkDf}}, \code{\link{fitlmer.mcmc}}
#' @examples
#' 
#' nmc <- 10
#' pmx <- cbind(randParams(genParamRanges(), nmc, 1001), seed=mkSeeds(nmc, 1001))
#' 
#' # between-items dataset
#' x.bi <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=TRUE)
#' 
#' # within-items dataset
#' x.wi <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=FALSE)
#' 
#' # maximal model
#' fitlmer(x.bi, wsbi=TRUE)
#' 
#' fitlmer(x.wi, wsbi=FALSE, ri.only=FALSE) # maximal
#' fitlmer(x.wi, wsbi=FALSE, ri.only=TRUE) # random intercepts only
#' 
#' @export fitlmer
fitlmer <- function(mcr.data, ri.only=FALSE, wsbi=FALSE) {
  xd <- mcr.data
  v1 <- c(3,rep(NA,4))
  names(v1) <- c("fm","t","chi","pt","pchi")
  p.mcmc <- NA  
  fullModel <- 0
                                        # fullModel codes:
                                        # 0 : alles gut
                                        # 1 : dropped one slope
                                        # 2 : dropped two slopes
                                        # 3 : didn't converge
                                        # 4 : comparison model without cond didn't converge
  library(lme4, quietly=TRUE)
  if (ri.only) {
    mf <- Resp ~ Cond + (1 | SubjID) + (1 | ItemID)
    xd.lmer <- tryCatch(lmer(mf,
                             family=gaussian, data=xd, na.action=na.omit, REML=FALSE),
                        warning = function(w) {return (NULL)},
                        error = function(e) {return (NULL)})
  } else {
                                        # rirs model try full model
    if (!wsbi) {
      mf <- Resp ~ Cond + (1 + Cond | SubjID) + (1 + Cond | ItemID)
    } else {
      mf <- Resp ~ Cond + (1 + Cond | SubjID) + (1 | ItemID)
    }
    xd.lmer <- tryCatch(lmer(mf,
                             family=gaussian, data=xd, na.action=na.omit, REML=FALSE),
                        warning = function(w) {return (NULL)},
                        error = function(e) {return (NULL)})
    if (is.null(xd.lmer)) {
      if (wsbi) {
        fullModel <- 1
        mf <- Resp ~ Cond + (1 | SubjID) + (1 | ItemID)
        xd.lmer <- tryCatch(lmer(mf,
                                 family=gaussian, data=xd, na.action=na.omit, REML=FALSE),
                            warning = function(w) {return (NULL)},
                            error = function(e) {return (NULL)})        
      } else {
                                        # turn off warnings and get partially converged model
        mywarn <- getOption("warn")
        options(warn=-1)
        xd.lmer.1 <- lmer(mf, family=gaussian, data=xd, na.action=na.omit, REML=FALSE)
        options(warn=mywarn)
        
                                        # choose one of the two random slopes to drop
        if (lme4:::VarCorr(xd.lmer.1)$SubjID[2,2] > lme4:::VarCorr(xd.lmer.1)$ItemID[2,2]) {
                                        # drop Item random slope
          fullModel <- 1
          mf <- Resp ~ Cond + (1 + Cond | SubjID) + (1 | ItemID)
          xd.lmer <- tryCatch(lmer(mf,
                                   family=gaussian, data=xd, na.action=na.omit, REML=FALSE),
                              warning = function(w) {return (NULL)},
                              error = function(e) {return (NULL)})
          if (is.null(xd.lmer)) {
            fullModel <- 2
            mf <- Resp ~ Cond + (1 | SubjID) + (1 | ItemID)          
            xd.lmer <- tryCatch(lmer(mf,
                                     family=gaussian, data=xd, na.action=na.omit, REML=FALSE),
                                warning = function(w) {return (NULL)},
                                error = function(e) {return (NULL)})
          } else {}
        } else {
                                        # drop Subj random slope
          fullModel <- 1
          mf <- Resp ~ Cond + (1 | SubjID) + (1 + Cond | ItemID)
          xd.lmer <- tryCatch(lmer(mf,
                                   family=gaussian, data=xd, na.action=na.omit, REML=FALSE),
                              warning = function(w) {return (NULL)},
                              error = function(e) {return (NULL)})
          if (is.null(xd.lmer)) {
            fullModel <- 2
                                        # drop both random slopes
            mf <- Resp ~ Cond + (1 | SubjID) + (1 | ItemID)          
            xd.lmer <- tryCatch(lmer(mf,
                                     family=gaussian, data=xd, na.action=na.omit, REML=FALSE),
                                warning = function(w) {return (NULL)},
                                error = function(e) {return (NULL)})
          } else {}
        }
      } # end non-convergence proc for wswi
    } else {} # end handling of non-converged model
  } # end random slopes
  if (is.null(xd.lmer)) {
    return(v1)
  } else {}
  mf2 <- as.formula(paste(deparse(mf),"-Cond"))
  xd.lmer.2 <- tryCatch(lmer(mf2,
                             family=gaussian, data=xd, na.action=na.omit, REML=FALSE),
                        warning=function(w) {return(NULL)},
                        error=function(e) {return(NULL)})
  if (is.null(xd.lmer.2)) {
    fullModel <- 4
    ts.chi <- NA
    p.chi <- NA
  } else {
    ts.chi <- deviance(xd.lmer.2)-deviance(xd.lmer)
    p.chi <- pchisq(abs(ts.chi), 1, lower.tail=F)
  }
  
  ts.tval <- abs(lme4:::fixef(xd.lmer)[2]/sqrt(Matrix:::diag(lme4:::vcov(xd.lmer))[2]))
  p.t <- 2*(1-pnorm(ts.tval))
  v1 <- c(fm=fullModel, t=ts.tval, chi=ts.chi, pt=p.t, pchi=p.chi)
  names(v1) <- c("fm", "t","chi","pt","pchi")
  return(v1)
}



#' Run mixed-effects model using (\code{\link{lmer}})
#' 
#' Fits a LMEM with by-subject (and by-item, when appropriate) random slopes,
#' and no random intercepts (one-factor design only).
#' 
#' 
#' @param mcr.data A dataframe formatted as described in \code{\link{mkDf}}.
#' @param wsbi Whether the design is between-items (TRUE) or within-items
#' (FALSE).
#' @return
#' @returnItem fm Code for the model that converged: (1) dropped one slope; (2)
#' dropped two slopes; (3) main model did not converge; (4) comparision model
#' for likelihood ratio test did not converge.
#' @returnItem t t-statistic for the treatment effect
#' @returnItem chi chi-square statistic for the likelihood ratio test (1 df)
#' @returnItem pt p-value for the t-statistic (normal distribution)
#' @returnItem pchi p-value for the chi-square statistic
#' @seealso \code{\link{mkDf}}, \code{\link{fitlmer.mcmc}}
#' @examples
#' 
#' nmc <- 10
#' pmx <- cbind(randParams(genParamRanges(), nmc, 1001), seed=mkSeeds(nmc, 1001))
#' 
#' # between-items dataset
#' x.bi <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=TRUE)
#' 
#' # within-items dataset
#' x.wi <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=FALSE)
#' 
#' # maximal model
#' fitrsonly(x.bi, wsbi=TRUE)
#' fitrsonly(x.wi, wsbi=FALSE)
#' 
#' @export fitrsonly
fitrsonly <- function(mcr.data, wsbi) {
  xd <- mcr.data
  if (wsbi) {
    mf <- Resp ~ Cond + (0 + Cond | SubjID) + (1 | ItemID)
  } else {
    mf <- Resp ~ Cond + (0 + Cond | SubjID) + (0 + Cond | ItemID)
  }
  mod <- .tryFit.default(mf, xd)
  return(.modInfo1(mf, xd, mod))
}



#' Run mixed-effects model using (\code{\link{lmer}})
#' 
#' Fits a LMEM with random slopes and random intercepts but no random
#' correlations (one-factor design only).
#' 
#' 
#' @param mcr.data A dataframe formatted as described in \code{\link{mkDf}}.
#' @param wsbi Whether the design is between-items (TRUE) or within-items
#' (FALSE).
#' @return
#' @returnItem fm Code for the model that converged: (1) dropped one slope; (2)
#' dropped two slopes; (3) main model did not converge; (4) comparision model
#' for likelihood ratio test did not converge.
#' @returnItem t t-statistic for the treatment effect
#' @returnItem chi chi-square statistic for the likelihood ratio test (1 df)
#' @returnItem pt p-value for the t-statistic (normal distribution)
#' @returnItem pchi p-value for the chi-square statistic
#' @seealso \code{\link{mkDf}}, \code{\link{fitlmer.mcmc}}
#' @examples
#' 
#' nmc <- 10
#' pmx <- cbind(randParams(genParamRanges(), nmc, 1001), seed=mkSeeds(nmc, 1001))
#' 
#' # between-items dataset
#' x.bi <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=TRUE)
#' 
#' # within-items dataset
#' x.wi <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=FALSE)
#' 
#' # maximal model
#' fitnocorr(x.bi, wsbi=TRUE)
#' fitnocorr(x.wi, wsbi=FALSE)
#' 
#' @export fitnocorr
fitnocorr <- function(mcr.data, wsbi=FALSE) {
  xd <- mcr.data
  if (wsbi) {
    mf <- Resp ~ Cond + (1 | SubjID) + (0 + Cond | SubjID) + (1 | ItemID)
  } else {
    mf <- Resp ~ Cond + (1 | SubjID) + (0 + Cond | SubjID) + (1 | ItemID) + (0 + Cond | ItemID)
  }
  mod <- .tryFit.default(mf, xd)
  return(.modInfo1(mf, xd, mod))
}



#' Run mixed-effects model (\code{\link{lmer}}) with MCMC p-value
#' 
#' Uses \code{\link{lmer}} to fit a random-intercepts only model, and then uses
#' \code{\link{mcmcsamp}} to obtain p-values.
#' 
#' If the model does not converge, returns \code{NA}.  The MCMC procedure is
#' based on Baayen's \code{pvals.fnc} in package \code{languageR}.
#' 
#' @param mcr.data A dataframe formatted as described in \code{\link{mkDf}}.
#' @param wsbi Whether the design is between-items (TRUE) or within-items
#' (FALSE).
#' @param nmcmc Number of Markov-Chain Monte Carlo simulations (default =
#' 10000).
#' @return A single element vector with the mcmc p-value.
#' @note This function ONLY fits models with independent random intercepts and
#' slopes.
#' @seealso \code{\link{fitnocorr}}
#' @examples
#' 
#' nmc <- 10
#' pmx <- cbind(randParams(genParamRanges(), nmc, 1001), seed=mkSeeds(nmc, 1001))
#' 
#' # between-items dataset
#' x.bi <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=TRUE)
#' 
#' # NB: small number of MCMC runs so that the example runs quickly
#' # increase the number of runs for stable results
#' fitnocorr.mcmc(x.bi, wsbi=TRUE, nmcmc=1000) 
#' 
#' @export fitnocorr.mcmc
fitnocorr.mcmc <- function(mcr.data, wsbi=FALSE, nmcmc=10000) {
  xd <- mcr.data
  p.mcmc <- NA
  library(lme4, quietly=TRUE)
  if (wsbi) {
    mf <- Resp ~ Cond + (1 | SubjID) + (0 + Cond | SubjID) + (1 | ItemID)
  } else {
    mf <- Resp ~ Cond + (1 | SubjID) + (0 + Cond | SubjID) + (1 | ItemID) + (0 + Cond | ItemID)
  }
  mod <- .tryFit.default(mf, xd)
  if (mod$converged) {
    mcmc <- try(lme4::mcmcsamp(mod$value, n = nmcmc, silent=TRUE))
    mcmcfixef <- t(mcmc@fixef)
    nr <- nrow(mcmcfixef)
    prop <- colSums(mcmcfixef > 0)/nr
    ans <- 2 * pmax(0.5/nr, pmin(prop, 1 - prop))
    names(ans) <- colnames(mcmcfixef)
    p.mcmc <- ans["Cond"]    
  } else {    
  }
  names(p.mcmc) <- "pmcmc"
  return(p.mcmc)
}



#' Fit lmer model using model selection
#' 
#' Analyzes data using \code{\link{lmer}}, using model selection to test
#' significance of random slope terms in the model (likelihood ratio tests).
#' Does forward and backward selection, starting with subject slope or item
#' slope first (one-factor design only)
#' 
#' 
#' @param mcr.data A dataframe formatted as described in \code{\link{mkDf}}.
#' Named thus to interface with the function \code{\link{mcRun}}.
#' @param wsbi Whether the design is between-items (TRUE) or within-items
#' (FALSE).
#' @param mf List of the models to be tested, in decreasing order of complexity
#' (see \code{\link{modSpace}} for examples of generating such a list).
#' @param crit alpha level for each likelihood-ratio test of slope variance.
#' @return A single row of a dataframe with number of fields depending on
#' \code{crit}.  Stepwise lmer models output six values for each alpha level
#' (i.e., level of \code{crit}.  These values are:
#' 
#' The values are given twice, once from each direction (forward or backward).
#' Thus, if there are two values of crit, there will be 2 (direction) x 6
#' (value) x 2 (levels of crit) = 24 values in each row of the dataframe.  To
#' assemble a dataframe of results from a file into a three-dimensional array,
#' see \code{\link{reassembleStepwiseFile}}.
#' @returnItem fm the model that was selected.  4 means that no model
#' converged.  For forward stepping models, 3 = maximal model was selected; 2 =
#' model includes only first slope; 1 = model is random intercept only.  For
#' backward stepping models, 1 = maximal model, 2 = model minus one slope, 3 =
#' random intercept model.
#' @returnItem t t-statistic for the treatment effect
#' @returnItem chi chi-square statistic for the likelihood ratio test (1 df)
#' @returnItem pt p-value for the t-statistic (normal distribution)
#' @returnItem pchi p-value for the chi-square statistic
#' @seealso \code{\link{fitlmer}}, \code{\link{modSpace}},
#' \code{\link{fitstepwise.bestpath}}, \code{\link{reassembleStepwiseFile}}
#' @examples
#' 
#' nmc <- 10
#' pmx <- cbind(randParams(genParamRanges(), nmc, 1001), seed=mkSeeds(nmc, 1001))
#' 
#' x8 <- mkDf(nsubj=24, nitem=24, pmx[8,], wsbi=FALSE)
#' mf <- modSpace(wsbi=FALSE) # figure out set of possible models
#' mf.sfirst <- c(mf[3], mf[[2]][1], mf[1]) # test subject slope first
#' mf.ifirst <- c(mf[3], mf[[2]][2], mf[1]) # test items slope first
#' 
#' # forward, subj first
#' fitstepwise(x8, wsbi=FALSE, mf=mf.sfirst, crit=.05)
#' 
#' # forward, item first
#' fitstepwise(x8, wsbi=FALSE, mf=mf.ifirst, crit=.05)
#' 
#' @export fitstepwise
fitstepwise <- function(mcr.data, wsbi, mf, crit=c(.01,.05,seq(.1,.8,.1))) {
  xd <- mcr.data
  mfits <- .stepwiseFit(mf, xd, crit)
  gg <- lapply(mfits[c("forward","backward")], .multiModInfo, xd=xd, wsbi=wsbi)
  return(unlist(gg))
}



#' Fit lmer model using model selection
#' 
#' Analyzes data using \code{\link{lmer}}, using model selection to
#' test significance of random slope terms in the model (likelihood
#' ratio tests).  Does forward or backward selection.  Unlike
#' \code{\link{fitstepwise}}, before taking a step forward (or
#' backward) both possible slopes are tested, and the "best" (most
#' conservative) step is taken.  NB: one-factor design only.
#' 
#' It only makes sense to run this with data from a within-subject/within-items
#' design (since for within-subject/between-items data, there is only one slope
#' to be tested, and the procedure is therefore equivalent to a simple stepwise
#' algorithm.
#' 
#' Note that for purposes of efficiency, forward and backward models are not
#' run simultaneously as with \code{\link{fitstepwise}}.
#' 
#' @param mcr.data A dataframe formatted as described in \code{\link{mkDf}}.
#' Named thus to interface with the function \code{\link{mcRun}}.
#' @param forward Should a forward model (TRUE) or backward model (FALSE) be
#' run?
#' @param crit alpha level for each likelihood-ratio test of slope variance.
#' @return A single row of a dataframe with number of fields depending on
#' \code{crit}.  Stepwise lmer models output six values for each alpha level
#' (i.e., level of \code{crit}.  These values are:
#' 
#' The values are given once for the specified direction (forward or backward).
#' Thus, if there are two values of crit, there will be 6 (value) x 2 (levels
#' of crit) = 12 values in each row of the dataframe.  To assemble a dataframe
#' of results from a file into a three-dimensional array, see
#' \code{\link{reassembleStepwiseFile}}.
#' @returnItem fm the model that was selected.  4 means that no model
#' converged.  For forward stepping models, 3 = maximal model was selected; 2 =
#' model includes only first slope; 1 = model is random intercept only.  For
#' backward stepping models, 1 = maximal model, 2 = model minus one slope, 3 =
#' random intercept model.
#' @returnItem t t-statistic for the treatment effect
#' @returnItem chi chi-square statistic for the likelihood ratio test (1 df)
#' @returnItem pt p-value for the t-statistic (normal distribution)
#' @returnItem pchi p-value for the chi-square statistic
#' @seealso \code{\link{fitlmer}}, \code{\link{modSpace}},
#' \code{\link{fitstepwise}}, \code{\link{reassembleStepwiseFile}}
#' @examples
#' 
#' nmc <- 10
#' pmx <- cbind(randParams(genParamRanges(), nmc, 1001), seed=mkSeeds(nmc, 1001))
#' 
#' x8 <- mkDf(nsubj=24, nitem=24, pmx[8,], wsbi=FALSE)
#' 
#' # forward
#' fitstepwise.bestpath(mcr.data=x8, forward=TRUE, crit=.05)
#' 
#' # backward
#' fitstepwise.bestpath(mcr.data=x8, forward=FALSE, crit=.05)
#' 
#' @export fitstepwise.bestpath
fitstepwise.bestpath <- function(mcr.data, forward, crit=c(.01,.05,seq(.1,.8,.1))) {
  xd <- mcr.data
  library(lme4, quietly=TRUE)

  mf <- .modSpace(FALSE)
  mf <- c(mf$min, mf$mid$srs, mf$mid$irs, mf$max)
  names(mf) <- c("min","srs","irs","max")
  mods <- lapply(mf[c("srs","irs")], .tryFit.default, xdat=xd)
  if (forward) {
    mods <- c(list(min=.tryFit.default(mf[["min"]], xdat=xd)), mods)
    cvg.lx <- unlist(lapply(mods, function(x) {x$converged}))
    if (sum(cvg.lx) < 3) { # at least one failed to converge
      if (sum(cvg.lx)==0) { # NONE converged
        ff <- .tryFit.default(mf[["max"]], xdat=xd)
        if (ff$converged) {
          ff.inf <- .modInfo(list(max=ff$value), xd, wsbi=FALSE)
          return(c(matrix(apply(t(ff.inf), 1, rep, length(crit)), nrow=length(crit), byrow=TRUE)))
        } else {
          return(rep(NA, length(crit)*6))
        }
      } else if (sum(cvg.lx)==1) { # one converged
        mf <- mf[c(names(cvg.lx)[cvg.lx], "max")] # the converged plus the maximal
      } else { # two converged
        if ("min" %in% (names(cvg.lx)[cvg.lx])) {
          mf <- mf[c(names(cvg.lx)[cvg.lx], "max")] # the min plus the converged
        } else { # the two "mid" models converged
                                        # pick based on which one has larger SD
          srs.sd <- attr(lme4:::VarCorr(mods[["srs"]]$value)$SubjID, "stddev")[2]
          irs.sd <- attr(lme4:::VarCorr(mods[["irs"]]$value)$ItemID, "stddev")[2]
          if (srs.sd > irs.sd) {
            mf <- mf[c("srs","max")] #c(mf$mid$srs, mf$max)
          } else {
            mf <- mf[c("irs","max")] #c(mf$mid$irs, mf$max) #mf[c("irs", "max")]
          }
        }
      }
    } else { # all three converged
      minvsrs <- .modCompare(list(mods[["min"]]$value, mods[["srs"]]$value))
      minvirs <- .modCompare(list(mods[["min"]]$value, mods[["irs"]]$value))
      if (minvsrs["p",1] < minvirs["p",1]) {
        mf <- mf[c("min","srs","max")] #c(mf$min, mf$mid$srs, mf$max) #mf[c("min", "srs", "max")]
      } else {
        mf <- mf[c("min","irs","max")] #c(mf$min, mf$mid$irs, mf$max) #mf[c("min", "irs", "max")]
      }
    }    
  } else { # backward model
    mods <- c(list(max=.tryFit.default(mf[["max"]], xdat=xd)), mods)
    cvg.lx <- unlist(lapply(mods, function(x) {x$converged}))
    if (sum(cvg.lx) < 3) { # at least one failed to converge
      if (sum(cvg.lx)==0) { # NONE converged
        ff <- .tryFit.default(mf[["min"]], xdat=xd)
        if (ff$converged) {
          ff.inf <- .modInfo(list(max=ff$value), xd, wsbi=FALSE)
          return(c(matrix(apply(t(ff.inf), 1, rep, 10), nrow=10, byrow=TRUE)))
        } else {
          return(rep(NA, 60))
        }
      } else if (sum(cvg.lx)==1) { # one converged
        mf <- mf[c("min", names(cvg.lx)[cvg.lx])] # the min plus the converged
      } else { # two converged
        if ("max" %in% (names(cvg.lx)[cvg.lx])) {
          mf <- mf[c("min", names(cvg.lx)[cvg.lx])] # the min plus the converged
        } else { # the two "mid" models converged
                                        # pick based on which one has smaller SD
          srs.sd <- attr(lme4:::VarCorr(mods[["srs"]]$value)$SubjID, "stddev")[2]
          irs.sd <- attr(lme4:::VarCorr(mods[["irs"]]$value)$ItemID, "stddev")[2]
          if (srs.sd < irs.sd) {
            mf <- mf[c("min","srs")] # drop out the larger one (mid.srs has no item slope)
          } else {
            mf <- mf[c("min","irs")] # drop out the larger one (mid.irs has no subj slope)
          }
        }
      }
    } else { # all three converged
      maxvsrs <- .modCompare(list(mods[["max"]]$value, mods[["srs"]]$value))
      maxvirs <- .modCompare(list(mods[["max"]]$value, mods[["irs"]]$value))
      if (maxvsrs["p",1] < maxvirs["p",1]) {
        mf <- mf[c("min","irs","max")] #c(mf$min, mf$mid$irs, mf$max) #mf[c("min", "mid.irs", "max")]
      } else {
        mf <- mf[c("min","srs","max")] #c(mf$min, mf$mid$srs, mf$max) #mf[c("min", "mid.srs", "max")]
      }
    }    
  }
  mfits <- .stepwiseFit(mf, xd, crit)
  gg <- lapply(mfits[c("forward","backward")], .multiModInfo, xd=xd, wsbi=FALSE)
  if (forward) {
    gg <- gg["forward"]
  } else {
    gg <- gg["backward"]
  }
  return(unlist(gg))
}


#' Make a dataframe with simulated data given a set of population parameters
#' 
#' Generates data for a simulated experiment given population parameters and
#' information about the design (one-factor design).
#' 
#' Note that for between-items designs, the variable \code{w11} in the
#' parameter vector is simply ignored.  The default value for \code{missMeth}
#' of 'random' will randomly remove between 0-5\% of the observations.  The
#' option 'randomBig' will randomly remove 10-80\%; 'bysubj' will randomly
#' remove 10-80\% of each subject's data; 'bycond' will randomly remove 10-80\%
#' of each condition's data; and 'bysubjcond' will randomly remove 10=80\% of
#' each subject/condition combination's data.
#' 
#' @param nsubj number of subjects in the experiment
#' @param nitem number of items in the experiment
#' @param mcr.params population parameters (typically a row from the matrix
#' generated by \code{\link{createParamMx}}
#' @param wsbi whether the design is between-items (TRUE) or within-items
#' (FALSE)
#' @param missMeth method used for generating missing data (default is
#' 'random'; alternatives are 'none', 'randomBig', 'bysubj', 'bycond',
#' bysubjcond'; see Details)
#' @param rigen use a 'random-intercepts-only' generative model (default FALSE)
#' @return A dataframe, with fields:
#' @returnItem SubjID a factor, identifying subject number
#' @returnItem ItemID a factor, identifying item number
#' @returnItem Cond treatment condition, deviation coded (-.5, .5)
#' @returnItem Resp response variable
#' @seealso \code{\link{genParamRanges}}, \code{\link{mkDf.facMixedAB}}
#' @examples
#' 
#' nmc <- 10
#' pmx <- cbind(randParams(genParamRanges(), nmc, 1001), seed=mkSeeds(nmc, 1001))
#' 
#' x.df <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=FALSE)
#' 
#' ## set rate of missing observations manually
#' pmx[,"miss"] <- 0.4
#' x.df2 <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=FALSE)
#' 
#' ## by-condition missing observation rates
#' pmx[,"pMin"] <- 0.1
#' pmx[,"pMax"] <- 0.8
#' x.df3 <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=FALSE,missMeth="bycond")
#' 
#' ## by-subject missing observation rates
#' x.df4 <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=FALSE,missMeth="bysubj")
#' 
#' ## by-subject/condition pair missing observation rates
#' x.df5 <- mkDf(nsubj=24, nitem=24, mcr.params=pmx[1,], wsbi=FALSE,missMeth="bysubjcond")
#' 
#' 
#' 
#' @export mkDf
#' @importFrom MASS mvrnorm
mkDf <- function(nsubj=24,    # number of subjects
                  nitem=24,    # number of items
                  wsbi=FALSE,     # if design is between items (TRUE) or within (FALSE)
                  mcr.params=createParamMx(1,firstseed=sample(1:1000000,1))[1,],
                                     # vector of parameters for generation (see createParamMx)
                  missMeth="random", # method for generating missing obs
                  rigen=FALSE) { # random-intercepts-only generative model
  #   This function creates simulated data given a row of parameters and
  #   no. subjects and no. items

  mkSigma <- function(v1, v2, r) {
    # make a variance-covariance matrix
    # from variance 1 (intercept), variance 2 (slope), and correlation
    cv <- r*sqrt(v1)*sqrt(v2)
    return(matrix(c(v1,cv,cv,v2),ncol=2))
  }

  set.seed(mcr.params[["seed"]])
  
  if (((nitem %% 2) != 0) || (nitem < 4)) {
    stop("number of items must be >= 4 and a factor of 2")
  } else {}
    
  exp.list <- expand.grid(list(ItemID=factor(1:nitem), SubjID=factor(1:nsubj)))[,2:1]    

  #library(MASS, quietly=TRUE)
  # assign condition variable and create item random-effects depending on design
  if (wsbi) {
    exp.list$Cond <- rep(c(-.5,.5),times=nitem/2)
  } else {
    exp.list$Cond <- c(rep(c(-.5,.5),times=nitem/2), rep(c(.5,-.5),times=nitem/2))
  }

  # subject random effects
  subj <- cbind(SubjID=factor(1:nsubj),
                mvrnorm(nsubj, mu=c(0,0),
                        Sigma=mkSigma(mcr.params[["t00"]], mcr.params[["t11"]], mcr.params[["rsub"]])))
  colnames(subj) <- c("SubjID","ri.subj","rs.subj")

  # item random effects
  item <- cbind(ItemID=factor(1:nitem),
                mvrnorm(nitem, mu=c(0,0),
                        Sigma=mkSigma(mcr.params[["w00"]], mcr.params[["w11"]], mcr.params[["ritm"]])) )
  colnames(item) <- c("ItemID","ri.item","rs.item")

  
  x <- merge(subj, merge(exp.list, item))[,c("SubjID","ItemID","Cond",
                                             colnames(subj)[-1],colnames(item)[-1])]
                                          # make it pretty
  x <- x[order(x$SubjID, x$ItemID),]; rownames(x) <- 1:nrow(x)

  x$err <- rnorm(nrow(x), mean=0, sd=sqrt(mcr.params[["err"]]))  # unexpl variance
  
  # calculate response variable
  x$Resp <- mcr.params[["int"]] + x$ri.subj + x$ri.item + # intercept
    mcr.params[["eff"]]*x$Cond + x$err
  if (!rigen) { # add in the slope, except for RI-only version
    x$Resp <- x$Resp + x$rs.subj*x$Cond # slope
  } else {}

  # for wswi designs, add in the random item slope
  if (!wsbi) {
    if (!rigen) {
      x$Resp <- x$Resp + x$rs.item*x$Cond
    } else {}
  } else {
    x <- x[,setdiff(colnames(x), "rs.item")]
  }


  ### missing data
  
  ## random deletion of some proportion of observations -- this is NOT
  ## stochastic deletion of each observation with some proportion p
  dropRand <- function(nsubj, nitem, pobs) {
    ##cat("Dropping random...\n")
    nmiss <- round((nsubj*nitem)*pobs)
    idx <- sample(c(rep(FALSE, nsubj*nitem-nmiss), rep(TRUE, nmiss)))
    return(idx)
  }

  dropRandBig <- function(pMin, pMax) {
    drops <- rbinom(nrow(x), 1, runif(1,pMin,pMax))
    return(drops>0)
  }
  
  ## condition-specific drop rates
  dropByCond <- function(pMin,pMax) {
    condRates <- runif(2,pMin,pMax)
    drops <- rbinom(dim(x)[1], 1, condRates[(x$Cond > 0) + 1])
    tmp <- drops * 1:length(drops)
    return(tmp[tmp>0])
  }

  ## subject-specific drop rates
  dropBySubj <- function(pMin,pMax) {
    subjRates <- runif(nsubj,pMin,pMax)
    drops <- rbinom(dim(x)[1], 1, subjRates[x$SubjID])
    tmp <- drops * 1:length(drops)
    return(tmp[tmp>0])
  }

  ## subj X condition specific drop rates
  dropBySubjCond <- function(pMin,pMax) {
    rates <- matrix(runif(2*nsubj,pMin,pMax),nsubj,2)
    drops <- rbinom(dim(x)[1], 1, rates[cbind(x$SubjID,(x$Cond > 0) + 1)])
    tmp <- drops * 1:length(drops)
    return(tmp[tmp>0])
  }

  dropNone <- function(nsubj,nitem) {
    return(rep(FALSE,nsubj*nitem))
  }

  # now randomly delete observations
  missing.lx <- switch(missMeth,
                       none=dropNone(nsubj, nitem),
                       random=dropRand(nsubj, nitem, mcr.params[["miss"]]),
                       randomBig=dropRandBig(mcr.params[["pMin"]],mcr.params[["pMax"]]),
                       bycond=dropByCond(mcr.params[["pMin"]],mcr.params[["pMax"]]),
                       bysubj=dropBySubj(mcr.params[["pMin"]],mcr.params[["pMax"]]),
                       bysubjcond=dropBySubjCond(mcr.params[["pMin"]],mcr.params[["pMax"]]))
  if (is.null(missing.lx)) {
    stop("invalid argument '", missMeth, "' for 'missMeth', must be one of 'none','random','randomBig','bycond','bysubj','bysubjcond'")
  } else {}
  x[missing.lx,"Resp"] <- NA   # delete them

  x <- x[,c("SubjID","ItemID","Cond","Resp")]
  x$SubjID <- factor(x$SubjID)
  return(x)
}


#' Parse and assemble data in a stepwise file into an array.
#' 
#' Convert file containing values from stepwise lmer models to a
#' three-dimensional array in the workspace (one-factor design only).
#' 
#' 
#' @param fname Name of the file containing output from \code{\link{mcRun}}
#' calling either \code{\link{fitstepwise}} or
#' \code{\link{fitstepwise.bestpath}}.
#' @return Returns either a three dimension array (if file is the result of
#' \code{\link{fitstepwise.bestpath}}) or a list of two three dimensional
#' arrays if \code{\link{fitstepwise}} was used, one for forward and one for
#' backward selection.  Each array has the following dimensions:
#' @returnItem run the Monte Carlo run
#' @returnItem crit the alpha level for selection
#' @returnItem params the resulting values from individual calls to
#' \code{fitstepwise} or \code{fitstepwise.bestpath}
#' @seealso \code{\link{fitstepwise}}, \code{\link{fitstepwise.bestpath}}.
#' @examples
#' 
#' nmc <- 10 # number of monte carlo simulations
#' create the parameter matrix
#' pmx <- cbind(randParams(genParamRanges(), nmc, 1001), seed=mkSeeds(nmc, 1001))
#'
#' # mf contains model formulae for the various models to fit
#' mf <- list(max=Resp ~ Cond + (1 + Cond | SubjID) + (1 + Cond | ItemID),
#'             mid=srs=Resp ~ Cond + (1 + Cond | SubjID) + (1 | ItemID),
#'             min=Resp ~ Cond + (1 | SubjID) + (1 | ItemID))
#' 
#' mcRun("fitstepwise", mcr.outfile="test.txt",
#'       mcr.xdatFnc="mkDf", mcr.varying=pmx, mf=mf, wsbi=FALSE)
#' 
#' ff <- reassembleStepwiseFile("test.txt")
#' ff$forward[,2,]  # alpha-level=.05
#' 
#' @export reassembleStepwiseFile
reassembleStepwiseFile <- function(fname) {
  mx <- read.csv(fname, header=TRUE)
  return(.reassembleStepwise(mx))
}
