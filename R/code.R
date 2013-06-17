# Matrix, lme4, MASS

.multicoreExists <- function() {
  return(is.element("multicore", installed.packages()[,1]))
}

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
  mf <- modSpace(FALSE)
  mf2 <- c(names(c(mf[1], mf[[2]][1:2], mf[3])), "nofit")
  gg <- (1:5)[names(mod)==mf2]
  if (length(gg)==0) {
    mf2 <- c(names(c(modSpace(FALSE), recursive=TRUE)), "nofit")
    gg <- (1:5)[names(mod)==mf2]
  } else {}
  return(gg)
}

.modInfo <- function(mod, xd, wsbi) {
                                        # do comparison model
  mf2 <- modSpace(wsbi)[[names(mod)]]
  if (is.null(mf2)) {
    mf2 <- c(modSpace(wsbi), recursive=TRUE)[[names(mod)]]
    if (is.null(mf2)) {
      mf2 <- modSpace(wsbi)[[2]][[names(mod)]]
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

#### EXPORTED FUNCTIONS

modSpace <- function(wsbi) {
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

tryFit.default <- function(mf, xdat) {
    tryFit(mf, xdat, na.action=na.omit, REML=FALSE)
}

### FUNCTIONS ABOVE HAVEN'T BEEN UPDATED YET
### FUNCTIONS BELOW DON'T NEED UPDATING

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

createParamMx <- function(nexp=100000,    # no. of simulated experiments for which to generate population parameters
                          simparam.env=getParamRanges(),
                                          # environment containing ranges of population params (usually param.env)
                          firstseed=666,  # seed to start things off (this function generates further seeds)
                          h0=TRUE,        # null hypothesis TRUE or FALSE
                          outfile=NULL) { # name of file to write the matrix to
  #   Generates a matrix of parameters.  Each row are the population parameters
  #   used to generate data from a single hypothetical "experiment."
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

loessPred <- function(pval, paramX, paramY,
                      minx=0.01,maxx=2.99,miny=.01,maxy=2.99,
                      svec=rep(TRUE, length(paramX)),span=.9) {
                                        # calculate predictions (Type I error or Power)
                                        # get rid of NAs, and use svec to select out observations
  lvec <- (!is.na(pval)) & svec
  ff <- pval[lvec]
  p.sig <- ifelse(ff<=.05,1,0)
  x <- paramX[lvec]; y <- paramY[lvec]
  p.map <- loess(p.sig ~ x * y, span=span)
  return(p.map)
}

heatmap2 <- function(lofit, labelX, labelY, title,
                     mxseq=seq(.01,2.99,length.out=300),
                     myseq=seq(.01,2.99,length.out=300),
                     powermap=FALSE) {
  p.grid <- expand.grid(x=mxseq,y=myseq)
  lpred.mx <- predict(lofit, newdata=p.grid)
  rownames(lpred.mx) <- mxseq
  colnames(lpred.mx) <- myseq
                                        # color scheme for heatmap
  makergb <- function(x,ming,maxg) {
    rge.x <- max(x)-min(x)
    px <- (x-min(x))/rge.x
    rge.g <- abs(ming-maxg)
    pg <- ming + px*rge.g*(ifelse(ming>maxg,-1,1))
    return(unlist(lapply(pg, function(y) {rgb(y,y,y)})))
  }
  if (!powermap) {
    mycol1 <- makergb(seq(.05,.09,.01), .95, .70)
    mycol2 <- makergb(seq(.10,.20,.01), .70, .50)
#    mycol2 <- makergb(seq(.10,.20,.01), .70, .40)    
    mycol3 <- makergb(seq(.21,.70,.01), .50, .25)
#    mycol3 <- makergb(seq(.21,.70,.01), .40, .15)
    mycol <- c(mycol1,mycol2,mycol3)
    mybreaks <- c(seq(.05,.70,.01),2.0)
  } else {
    mycol1 <- makergb(seq(.05,.09,.01), .95, .90)
    mycol2 <- makergb(seq(.10,.20,.01), .90, .70)
    mycol3 <- makergb(seq(.21,.70,.01), .70, .25)
    mycol <- c(mycol1,mycol2,mycol3)
    mybreaks <- c(seq(.05,.70,.01),2.0)
  }
  lab2expr <- function(x) { # map label to expression
    ee <- x
    if (!is.expression(ee)) {
      if (x=="t00") {
        ee <- expression({tau[0][0]}^2)
      } else if (x == "t11") {
        ee <- expression({tau[1][1]}^2)
      } else if (x == "w00") {
        ee <- expression({omega[0][0]}^2)
      } else if (x == "w11") {
        ee <- expression({omega[1][1]}^2)
      } else if (x == "subj_logratio") {
        ee <- expression(paste("log (", {tau[1][1]}/{tau[0][0]}, ")", sep=""))
      } else if (x == "item_logratio") {
        ee <- expression(paste("log (", {omega[1][1]}/{omega[0][0]}, ")", sep=""))
      } else if (x == "subj_cr") {
        ee <- expression(paste("log (", {tau[1][1]}/{omega[0][0]}, ")", sep=""))
      } else if (x == "item_cr") {
        ee <- expression(paste("log (", {omega[1][1]}/{tau[0][0]}, ")", sep=""))
      }
    } else {}
    return(ee)
  } 
  e1 <- lab2expr(labelX); e2 <- lab2expr(labelY) 
  
  mxseq <- as.numeric(rownames(lpred.mx))
  myseq <- as.numeric(colnames(lpred.mx))
  minx <- floor(min(mxseq)); maxx <- ceiling(max(mxseq))
  miny <- floor(min(myseq)); maxy <- ceiling(max(myseq))
  image(z=lpred.mx,
        x=mxseq, y=myseq,
        xlab=e1,ylab=e2,
        breaks=mybreaks, col=mycol,
        main=title, cex.lab=1.4,
        xaxt='n',yaxt='n',useRaster=TRUE)
  
  contour(lpred.mx, x=mxseq,y=myseq, add=TRUE, levels=c(seq(0,.2,by=.01),seq(.25,.5,by=.05),seq(.6,1,.1)))
  axis(1,at=round(seq(minx,maxx,length.out=7),2))
  axis(2,at=round(seq(miny,maxy,length.out=7),2))  
}

### FUNCTIONS BELOW HAVE BEEN UPDATED

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

fitrsonly <- function(mcr.data, wsbi) {
  xd <- mcr.data
  if (wsbi) {
    mf <- Resp ~ Cond + (0 + Cond | SubjID) + (1 | ItemID)
  } else {
    mf <- Resp ~ Cond + (0 + Cond | SubjID) + (0 + Cond | ItemID)
  }
  mod <- tryFit.default(mf, xd)
  return(.modInfo1(mf, xd, mod))
}

fitnocorr <- function(mcr.data, wsbi=FALSE) {
  xd <- mcr.data
  if (wsbi) {
    mf <- Resp ~ Cond + (1 | SubjID) + (0 + Cond | SubjID) + (1 | ItemID)
  } else {
    mf <- Resp ~ Cond + (1 | SubjID) + (0 + Cond | SubjID) + (1 | ItemID) + (0 + Cond | ItemID)
  }
  mod <- tryFit.default(mf, xd)
  return(.modInfo1(mf, xd, mod))
}

fitnocorr.mcmc <- function(mcr.data, wsbi=FALSE, nmcmc=10000) {
  xd <- mcr.data
  p.mcmc <- NA
  library(lme4, quietly=TRUE)
  if (wsbi) {
    mf <- Resp ~ Cond + (1 | SubjID) + (0 + Cond | SubjID) + (1 | ItemID)
  } else {
    mf <- Resp ~ Cond + (1 | SubjID) + (0 + Cond | SubjID) + (1 | ItemID) + (0 + Cond | ItemID)
  }
  mod <- tryFit.default(mf, xd)
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

fitstepwise <- function(mcr.data, wsbi, mf, crit=c(.01,.05,seq(.1,.8,.1))) {
  xd <- mcr.data
  mfits <- .stepwiseFit(mf, xd, crit)
  gg <- lapply(mfits[c("forward","backward")], .multiModInfo, xd=xd, wsbi=wsbi)
  return(unlist(gg))
}

fitstepwise.bestpath <- function(mcr.data, forward, crit=c(.01,.05,seq(.1,.8,.1))) {
  xd <- mcr.data
  library(lme4, quietly=TRUE)

  mf <- modSpace(FALSE)
  mf <- c(mf$min, mf$mid$srs, mf$mid$irs, mf$max)
  names(mf) <- c("min","srs","irs","max")
  mods <- lapply(mf[c("srs","irs")], tryFit.default, xdat=xd)
  if (forward) {
    mods <- c(list(min=tryFit.default(mf[["min"]], xdat=xd)), mods)
    cvg.lx <- unlist(lapply(mods, function(x) {x$converged}))
    if (sum(cvg.lx) < 3) { # at least one failed to converge
      if (sum(cvg.lx)==0) { # NONE converged
        ff <- tryFit.default(mf[["max"]], xdat=xd)
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
    mods <- c(list(max=tryFit.default(mf[["max"]], xdat=xd)), mods)
    cvg.lx <- unlist(lapply(mods, function(x) {x$converged}))
    if (sum(cvg.lx) < 3) { # at least one failed to converge
      if (sum(cvg.lx)==0) { # NONE converged
        ff <- tryFit.default(mf[["min"]], xdat=xd)
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
                MASS:::mvrnorm(nsubj, mu=c(0,0),
                        Sigma=mkSigma(mcr.params[["t00"]], mcr.params[["t11"]], mcr.params[["rsub"]])))
  colnames(subj) <- c("SubjID","ri.subj","rs.subj")

  # item random effects
  item <- cbind(ItemID=factor(1:nitem),
                MASS:::mvrnorm(nitem, mu=c(0,0),
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

mcRun <- function(mcr.FUN,  # name of function to be applied
                  mcr.cluster=NULL, # parallel cluster; NULL for single cluster
                   mcr.outfile=tempfile(fileext=".csv"), # name of text file to write to
                   mcr.xdatFnc=NULL,
                   mcr.xdat=NULL,
                   mcr.constant=NULL,
                   mcr.varying=NULL,
                   mcr.LoadOnExit=TRUE,
                   mcr.nCores=NULL,
                   mcr.reportInt=100,
                  ...) { # user parameters
    statusUpdate <- function(i, elapsed) {
        fmtstr1 <- paste("%",nchar(as.character(nEpochs)),"d",
                         sep="")
        fmtstr1.2 <- paste("%", nchar(as.character(nrow(mcr.varying))),"d", sep="")
        fmtstr2 <- paste(fmtstr1, "/", nEpochs, " (", fmtstr1.2,
                         "/", nrow(mcr.varying),
                         ") ", ":", sep="")
        stmp1 <- sprintf(fmtstr2, i, ix.b[i], nrow(mcr.varying))
        
        e2 <- elapsed[1:i]
        efac <- mean(e2)
        cat(stmp1, sprintf("%1.3fs/sweep, ", efac))
        totsecs <- efac*(nEpochs-i)
        cat(sprintf("%02dh:", floor(totsecs / 3600)))
        if (totsecs >= 3600) {
            remn <- totsecs %% 3600
        } else {
            remn <- totsecs
        }
        cat(sprintf("%02dm:", floor(remn/60)))
        cat(sprintf("%02ds", round(remn %% 60)), "left\n")
        flush.console()
    }

                                        # basic error checking
    if (missing(mcr.FUN)) {
        stop("need to supply 'FUN' argument to mcRun")
    } else {}
    
    if (is.function(mcr.FUN)) {
        mcr.FUN <- as.character(substitute(mcr.FUN))
    } else {}

    fbDataArgument <- c("mcr.data") %in% names(formals(mcr.FUN)) # do we need to pass data to the function?
    if (fbDataArgument && is.null(mcr.xdatFnc) && is.null(mcr.xdat)) {
        stop("function '", mcr.FUN, "' needs data (has 'mcr.data' argument), but mcRun arguments 'mcr.xdat' and 'mcr.xdatFnc' were not defined!")
    } else {}

    if (!fbDataArgument && is.null(mcr.xdatFnc) && is.null(mcr.xdat)) {
        warning("data arguments supplied to mcRun ('mcr.xdat' and/or 'mcr.xdatFnc'), but target function '", mcr.FUN, "' has no mcr.data argument defined!")
    } else {}
  
    if (!is.null(mcr.xdatFnc) && !is.null(mcr.xdat)) {
        stop("can only define one of 'xdatFnc' OR 'xdat', NOT both! (see ?mcRun for details)")
    } else {}

    if (is.function(mcr.xdatFnc)) {
        mcr.xdatFnc <- as.character(substitute(mcr.xdatFnc))
    } else {}

    if (!is.null(mcr.constant)) {
        if (!is.list(mcr.constant)) {
            stop("mcr.constant must be a list")
        } else {}
    } else {}
  
    # for 'varying' arguments that are matrices, make a data frame
    if (is.matrix(mcr.varying)) {
        mcr.varying <- as.data.frame(mcr.varying)
    } else if (is.list(mcr.varying) && !is.data.frame(mcr.varying)) {
    # for 'varying' arguments that are lists,
    # compile a data frame    
        ff <- lapply(mcr.varying, function(x) {
            if (is.list(x)) {
                as.data.frame(x)
            } else {
                x
            }
        })
        mcr.varying <- do.call("rbind", ff)
    } else {}

    # prepare arguments for the data generation function
    # if there is no function, create one that just returns mcr.xdat
    if (!is.null(mcr.xdat)) {
        .mcr.mkDat <- function() {return(mcr.xdat)}
        mcr.xdatFnc <- ".mcr.mkDat"
    } else {  # data generation function
        dparams <- list(...)[intersect(names(formals(mcr.xdatFnc)), names(list(...)))]
    }

    # prepare arguments for the target function
    fbParamArgument <- c("mcr.params") %in% names(formals(mcr.FUN)) # do we need to pass data to the function?
    fparams <- list(...)[intersect(names(formals(mcr.FUN)), names(list(...)))]

    # write the function that we will call repeatedly
    .mcr.doOnce <- function(theseParams) {
        if (fbDataArgument) { # create data if necessary
            xdat <- do.call(mcr.xdatFnc, c(dparams, list(mcr.params=c(mcr.constant, theseParams))))
            fparams2 <- c(fparams, list(mcr.data=xdat))
        } else {
            xdat <- NULL
            fparams2 <- fparams
        }
        if (fbParamArgument) {
            fparams3 <- c(fparams2, list(mcr.params=c(mcr.constant, theseParams)))
        } else {
            fparams3 <- fparams2
        }
        res <- do.call(mcr.FUN, fparams3)
        return(res)
    }
  
    # calculate how many epochs ("sweeps") to perform
    fullEpochs <- floor(nrow(mcr.varying)/mcr.reportInt)
    remainderEpochs <- ((nrow(mcr.varying) %% mcr.reportInt)>0)*1
    nEpochs <- fullEpochs + remainderEpochs
    ix.a <- (0:(nEpochs-1))*mcr.reportInt+1
    ix.b <- ix.a + c(rep(mcr.reportInt, fullEpochs), rep(nrow(mcr.varying)%%mcr.reportInt, remainderEpochs))-1
  
    # set up the result matrix
    ff <- .mcr.doOnce(as.list(mcr.varying[1,]))
    mcr.colnames <- names(ff)
    mcr.nelements <- length(ff)
    mcr.orig <- ff

    # open the file
    mcr.con <- file(mcr.outfile, "w", blocking=FALSE)
    if (length(mcr.colnames)==0) {
        mcr.colnames <- paste("V", 1:mcr.nelements, sep="")
    } else {}
    cat(mcr.colnames, file=mcr.con, append=FALSE, sep=",")
    cat("\n", append=TRUE, file=mcr.con)
    elapsed <- rep(NA, nEpochs)

    # export the data function to cluster
    if (!is.null(mcr.cluster)) {
        parallel:::clusterExport(cl, c(mcr.xdatFnc, mcr.FUN))
    } else {}
    
    # ready to go, start running
    for (i in 1:nEpochs) {
        elapsed[i] <-
            system.time( {
                todo <- split(mcr.varying[ix.a[i]:ix.b[i],], ix.a[i]:ix.b[i])
                if (!is.null(mcr.cluster)) {
                    ff <- parallel:::parLapplyLB(mcr.cluster, todo, .mcr.doOnce)
                                        # ff <- mclapply(todo, .mcr.doOnce, mc.cores=mcr.nCores)
                } else {
                    ff <- lapply(todo, .mcr.doOnce)
                }
                for (ii in 1:length(ff)) {
                    cat(ff[[ii]], file=mcr.con, append=TRUE, sep=",")
                    cat("\n", file=mcr.con)
                }
            })["elapsed"]
        statusUpdate(i, elapsed)
    }

    close(mcr.con)
  
    cat("total time (seconds): ", sum(elapsed), "\n")

    if (mcr.LoadOnExit) {
        result <- read.csv(mcr.outfile, header=TRUE)
        return(invisible(result))
    } else {
        return(NULL)
    }
}


reassembleStepwiseFile <- function(fname) {
  mx <- read.csv(fname, header=TRUE)
  return(.reassembleStepwise(mx))
}

###########################################
# adding mixed factorial designs
# these were added for higher order designs

defaultGenParamRanges <- function() {
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

defaultGenParamRanges.facMixedAB <- function() {
    list(t_00=c(1,3), # random intercept variance
         t_11=c(1,3), # random slope variance
         r_10=c(-.9,.9), # random intercept/slope correlation
         evar=3, # error variance
         mu=list(c(-2,-1), c(1,2)),
         A=list(c(-2,-1), c(1,2)),
         B=list(c(-2,-1), c(1,2)),
         AB=list(c(-2,-1), c(1,2)))
}

mkParamMx <- function(param.list, nmc=1000, firstseed=NULL) {
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
    if ("seed" %in% names(param.list)) {
        stop("cannot use 'seed' as the name of a generative parameter")
    } else {}
    
    if (!is.null(firstseed)) {
        set.seed(firstseed)
    } else {}


    ff <- lapply(names(param.list), function(nx) {
        x <- param.list[[nx]]
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
    mx <- matrix(unlist(ff), nrow=nmc, dimnames=list(NULL, names(param.list)))
    cbind(mx, seed=randSeed(nmc))
}

contr.deviation <- function(n) {
    if (length(n)<=1L) {
        if (is.numeric(n) && length(n) == 1L && n > 1L)
            levels <- seq_len(n)
        else stop("not enough degrees of freedom to define contrasts")
    } else {}
    apply(contr.treatment(n), 2, function(x) {x-mean(x)})
}

mkDf.facMixedAB <- function(mcr.params=mkParamMx(defaultGenParamRanges.facMixedAB(), 1)[1,],
                            nsubj=24, nreps=2, deviationCoding=TRUE) {
                                        # mcr.params is a vector of parameters for data generation
  x <- mcr.params
  if (!is.list(x)) {
      x <- as.list(x)
  } else {}
  set.seed(x[["seed"]])  
  library(MASS)
  ranef <- mvrnorm(nsubj, mu=c(0,0),
                   Sigma=matrix(c(x[["t_00"]], x[["r_10"]]*sqrt(x[["t_00"]])*sqrt(x[["t_11"]]),
                     x[["r_10"]]*sqrt(x[["t_00"]])*sqrt(x[["t_11"]]), x[["t_11"]]), nrow=2))
  srs <- rep(ranef[,2], each=2)*c(-.5,.5)
  dat <- data.frame(SubjID=factor(rep(1:nsubj,each=2*nreps)),
                    A=factor(rep(c("A1","A2"),each=2*nreps)),
                    B=factor(rep(c("B1","B2"),each=nreps)))
  if (deviationCoding) {
      contrasts(dat$A) <- contr.deviation(levels(dat$A))
      contrasts(dat$B) <- contr.deviation(levels(dat$B))
  } else { # use whatever the default coding is (usually 'treatment')
  }
  dat$Y <- x[["mu"]] + # fixed intercept
    rep(c(-x[["A"]],x[["A"]])/2, each=2*nreps) + # fixed effect of A
    rep(c(-x[["B"]],x[["B"]])/2, each=nreps) +   # fixed effect of B
    rep(c(x[["AB"]],-x[["AB"]],-x[["AB"]],x[["AB"]])/4, each=nreps) + # fixed interaction
    rep(ranef[,1], each=2*nreps) + # random intercept
    rep(srs, each=nreps) + # random slope
    rnorm(nrow(dat), sd=sqrt(x[["evar"]])) # error variance
  return(dat)
}

# fitlmer.facMixedAB()

getLmer.pValue <- function(m1,m2) {
    pval <- NA
    if (m1$converged && m2$converged) {
        chisq.val <- (-2*lme4::logLik(m1$value))-(-2*lme4::logLik(m2$value))
        pval <- as.numeric(pchisq(abs(chisq.val), 1, lower.tail=FALSE))
    } else {}
    return(pval)
}

fitlmer.facMixedAB <- function(mcr.data) {
    # models to fit
    models <- list(noB1=Y~A*B+(1|SubjID),
                   noB2=Y~A+B+(1|SubjID),
                   max1=Y~A*B+(1+A*B|SubjID),
                   max2=Y~A+B+(1+A*B|SubjID))
    # with treatment coding
    contrasts(mcr.data$A) <- contr.treatment(levels(mcr.data$A))
    contrasts(mcr.data$B) <- contr.treatment(levels(mcr.data$B))
    treat.mods <- lapply(models, tryFit, tf.data=mcr.data)
    # with deviation coding
    contrasts(mcr.data$A) <- contr.deviation(levels(mcr.data$A))
    contrasts(mcr.data$B) <- contr.deviation(levels(mcr.data$B))
    dev.mods <- lapply(models, tryFit, tf.data=mcr.data)
    treat.noB.p <- getLmer.pValue(treat.mods$noB2, treat.mods$noB1)
    treat.max.p <- getLmer.pValue(treat.mods$max2, treat.mods$max1)
    dev.noB.p <- getLmer.pValue(dev.mods$noB2, dev.mods$noB1)
    dev.max.p <- getLmer.pValue(dev.mods$max2, dev.mods$max1)
    dat.aov <- summary(aov(Y~A*B + Error(SubjID/B), mcr.data))
    return(c(noB.t=treat.noB.p, max.t=treat.max.p,
                        noB.d=dev.noB.p, max.d=dev.max.p,
                        aov=dat.aov[["Error: SubjID:B"]][[1]]$`Pr(>F)`[2]))
}
