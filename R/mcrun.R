#' Monte Carlo Simulation
#' 
#' Run Monte Carlo simulations, taking advantage of multiple processing cores
#' or a computing cluster
#' 
#' The number of simulations that will be run is given by
#' \code{nrow(mcr.varying)}, the number of rows in the parameter matrix passed
#' as an argument to the function.  The results matrix will be stored in the
#' file whose path is given by \code{mcr.outfile} in Comma Separated Values
#' (CSV) format.
#' 
#' If \code{mcr.fn} is a user-supplied function, the corresponding function
#' must return a vector with all elements of mode 'numeric', and accept
#' arguments \code{mcr.data} and (optionally) \code{mcr.params}.  All arguments
#' in \code{mcr.fnList} will be passed to that function.
#' 
#' If \code{mcr.datFn} is a user-supplied function, the corresponding
#' function must return a data frame and accept the argument \code{mcr.params}.
#' All arguments in \code{mcr.datArgs} will be passed to that function.
#' 
#' @param mcr.fn A character string naming the function that performs
#' the analysis: for example, for one-factor designs this could be any
#' one of \code{\link{fitanova}}, \code{\link{fitlmer}},
#' \code{\link{fitlmer.mcmc}}, \code{\link{fitstepwise}},
#' \code{\link{fitstepwise.bestpath}}, \code{\link{fitrsonly}},
#' \code{\link{fitnocorr}}, or \code{\link{lrCompare}}.  May also be a
#' user-supplied function (see 'Details').
#' @param mcr.fnArgs Arguments to be passed from \code{mcRun} to
#' the function named by \code{mcr.fn}.
#' @param mcr.cluster A processing cluster object, typically created by a call
#' to \code{makeCluster} from the package \code{parallel}.  If NULL (default),
#' uses a single processing core.
#' @param mcr.outfile name of file to which to write the output ("comma
#' separated value" format).  Defaults to "out.csv".
#' @param mcr.datFn Name of function to be used to generate data set for each
#' Monte Carlo run (e.g., 'mkDf').  May be a user-supplied function that
#' returns a data frame.  If a the same dataset is to be used for all runs,
#' then set to 'NULL' and define \code{mcr.dat} instead.
#' @param mcr.datArgs Arguments to be passed from \code{mcRun} to the
#' function named by \code{mcr.datFn}.
#' @param mcr.dat A data frame to be used in all Monte Carlo runs.  If the
#' data are to be generated by a function, set to 'NULL' and use
#' \code{mcr.datFn} instead.
#' @param mcr.constant list of parameter values that are constant over all
#' runs.
#' @param mcr.varying A matrix, data frame or list of lists containing the
#' parameters to be varied over runs, with each row corresponding to a single
#' run.  This determines the number of Monte Carlo runs (see 'Details').
#' @param mcr.LoadOnExit whether the data should be loaded from the file upon
#' completion into an R object, and passed on as the return value from the
#' function.
#' @param mcr.reportInt Interval at which to give status updates on progress
#' (default = every 100 runs)
#' @return If mcr.LoadOnExit is true, loads the data from the CSV file and
#' returns it to the calling function.
#' @examples
#' 
#' nmc <- 10
#' pmx <- cbind(randParams(genParamRanges(), nmc, 1001), seed=mkSeeds(nmc, 1001))
#'
#' cl <- NULL # single-core processing
#' # alternatively:
#' # cl <- parallel:::makeCluster(number.of.clusters)
#' # see the 'parallel' package
#' 
#' mx <- mcRun("fitanova", mcr.cluster=cl,
#' mcr.fnArgs=list(wsbi=TRUE), # pass this along to fitanova
#' mcr.varying=pmx, # parameters that are varying; each row is a single run
#' mcr.datFn="mkDf",  # data-generating function to call
#' mcr.datArgs=list(nitem=12, wsbi=TRUE), # values to be passed along to mkDf
#' mcr.reportInt=5) # report progress every 5 runs
#' 
#' @export mcRun
mcRun <- function (mcr.fn,
                    mcr.fnArgs = NULL,
                    mcr.cluster = NULL,
                    mcr.outfile = tempfile(fileext = ".csv"), 
                    mcr.datFn = NULL,
                    mcr.datArgs = NULL, mcr.dat = NULL,
                    mcr.constant = NULL,
                    mcr.varying = NULL,
                    mcr.LoadOnExit = TRUE, mcr.reportInt = 100) {
    statusUpdate <- function(i, elapsed, tot_epochs) {
        fmtstr1 <- paste("%", nchar(as.character(nEpochs)), "d", 
            sep = "")
        fmtstr1.2 <- paste("%", nchar(as.character(tot_epochs)), 
            "d", sep = "")
        fmtstr2 <- paste(fmtstr1, "/", nEpochs, " (", fmtstr1.2, 
            "/", tot_epochs, ") ", ":", sep = "")
        stmp1 <- sprintf(fmtstr2, i, ix.b[i], tot_epochs)
        e2 <- elapsed[1:i]
        efac <- mean(e2)
        cat(stmp1, sprintf("%1.3fs/sweep, ", efac))
        totsecs <- efac * (nEpochs - i)
        cat(sprintf("%02dh:", floor(totsecs/3600)))
        if (totsecs >= 3600) {
            remn <- totsecs%%3600
        } else {
            remn <- totsecs
        }
        cat(sprintf("%02dm:", floor(remn/60)))
        cat(sprintf("%02ds", round(remn%%60)), "left\n")
        flush.console()
    }
    if (missing(mcr.fn)) {
        stop("need to supply 'FUN' argument to mcRun")
    } else {}
    if (is.function(mcr.fn)) {
        mcr.fn <- as.character(substitute(mcr.fn))
    } else {}
    fbDataArgument <- c("mcr.data") %in% names(formals(mcr.fn))
    if (fbDataArgument && is.null(mcr.datFn) && is.null(mcr.dat)) {
        stop("function '", mcr.fn, "' needs data (has 'mcr.data' argument), but mcRun arguments 'mcr.dat' and 'mcr.datFn' ", 
            "were not defined!")
    } else {}
    if (!fbDataArgument && is.null(mcr.datFn) && is.null(mcr.dat)) {
        warning("data arguments supplied to mcRun ('mcr.dat' and/or 'mcr.datFn'), but target function '", 
            mcr.fn, "' has no mcr.data argument defined!")
    } else {}
    if (!is.null(mcr.datFn) && !is.null(mcr.dat)) {
        stop("can only define one of 'mcr.datFn' OR 'mcr.dat', NOT both! (see ?mcRun for details)")
    } else {}
   if (!is.character(mcr.datFn)) {
        if (is.function(mcr.datFn)) {
            mcr.datFn <- as.character(substitute(mcr.datFn))
        } else {}
    } else {}
    if (!is.null(mcr.constant)) {
        if (!is.list(mcr.constant)) {
            stop("mcr.constant must be a list")
        } else {}
    } else {}
    ## covert a matrix to a data frame
    if (is.matrix(mcr.varying)) {
        mcr.varying <- as.data.frame(mcr.varying)
    } else {}
    if (is.data.frame(mcr.varying)) {
        ## it's a data frame: split it up by row
        mcr.varying <- split(mcr.varying, seq_len(nrow(mcr.varying)))
    } else {}
    tot_epochs <- if (is.data.frame(mcr.varying) || is.matrix(mcr.varying))
                      nrow(mcr.varying) else length(mcr.varying)
    if (!is.null(mcr.dat)) {
        .mcr.mkDat <- function() {
            return(mcr.dat)
        }
        mcr.datFn <- ".mcr.mkDat"
    } else {
        datArgs <- formals(mcr.datFn)
        argTest <- names(mcr.datArgs) %in% names(datArgs)
        if (sum(argTest) < length(argTest) && !("..." %in% names(formals(mcr.datFn)))) {
            badArgs <- paste(names(mcr.datArgs)[!argTest], collapse = ", ")
            stop("arg(s) '", badArgs, "' not in list of formal args for function '", 
                mcr.datFn, "'")
        } else {}
    }
    fbParamArgument <- c("mcr.params") %in% names(formals(mcr.fn))
    fparams <- mcr.fnArgs
    .mcr.doOnce <- function(theseParams, fbDataArg, fbParamArg, 
        datFn, datArgs, parConst, fpar, fn) {
        if (fbDataArg) {
            xdat <- do.call(datFn, c(datArgs, list(mcr.params = c(parConst, 
                theseParams))))
            fparams2 <- c(fpar, list(mcr.data = xdat))
        } else {
            xdat <- NULL
            fparams2 <- fpar
        }
        if (fbParamArg) {
            fparams3 <- c(fparams2, list(mcr.params = c(parConst, 
                theseParams)))
        } else {
            fparams3 <- fparams2
        }
        res <- do.call(fn, fparams3)
        return(res)
    }
    fullEpochs <- floor(tot_epochs / mcr.reportInt)
    remainderEpochs <- ((tot_epochs %% mcr.reportInt) > 
        0) * 1
    nEpochs <- fullEpochs + remainderEpochs
    ix.a <- (0:(nEpochs - 1)) * mcr.reportInt + 1
    ix.b <- ix.a + c(rep(mcr.reportInt, fullEpochs), rep(tot_epochs %% mcr.reportInt, 
        remainderEpochs)) - 1
    ff <- .mcr.doOnce(as.list(mcr.varying[[1]]), fbDataArg = fbDataArgument, 
        fbParamArg = fbParamArgument, datFn = mcr.datFn, datArgs = mcr.datArgs, 
        parConst = mcr.constant, fpar = fparams, fn = mcr.fn)
    mcr.colnames <- names(ff)
    mcr.nelements <- length(ff)
    mcr.orig <- ff
    mcr.con <- file(mcr.outfile, "w", blocking = FALSE)
    if (length(mcr.colnames) == 0) {
        mcr.colnames <- paste("V", 1:mcr.nelements, sep = "")
    }
    else {
    }
    cat(mcr.colnames, file = mcr.con, append = FALSE, sep = ",")
    cat("\n", append = TRUE, file = mcr.con)
    elapsed <- rep(NA, nEpochs)
    if (!is.null(mcr.cluster)) {
        parallel::clusterExport(mcr.cluster, c(mcr.datFn, mcr.fn))
    }
    else {
    }
    for (i in 1:nEpochs) {
        elapsed[i] <- system.time({
            todo <- mcr.varying[ix.a[i]:ix.b[i]]
            if (!is.null(mcr.cluster)) {
                ff <- parallel:::parLapplyLB(mcr.cluster, todo, 
                  .mcr.doOnce, fbDataArg = fbDataArgument, fbParamArg = fbParamArgument, 
                  datFn = mcr.datFn, datArgs = mcr.datArgs, parConst = mcr.constant, 
                  fpar = fparams, fn = mcr.fn)
            }
            else {
                ff <- lapply(todo, .mcr.doOnce, fbDataArg = fbDataArgument, 
                  fbParamArg = fbParamArgument, datFn = mcr.datFn, 
                  datArgs = mcr.datArgs, parConst = mcr.constant, 
                  fpar = fparams, fn = mcr.fn)
            }
            for (ii in 1:length(ff)) {
                cat(ff[[ii]], file = mcr.con, append = TRUE, 
                  sep = ",")
                cat("\n", file = mcr.con)
            }
        })["elapsed"]
        statusUpdate(i, elapsed, tot_epochs)
    }
    close(mcr.con)
    cat("total time (seconds): ", sum(elapsed), "\n")
    if (mcr.LoadOnExit) {
        result <- read.csv(mcr.outfile, header = TRUE)
        return(invisible(result))
    }
    else {
        return(NULL)
    }
}
