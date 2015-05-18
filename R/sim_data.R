#' Calculate marginal and cell counts for factorially-designed data
#'
#' @param iv_names Names of independent variables in data.frame given by \code{dat}.
#' @param dat A data frame
#' @param unit_names Names of the fields containing sampling units (subjects, items)
#' @return a list, the elements of which have the marginal/cell counts for each factor in the design
#' @export
fac_counts <- function(iv_names, dat, unit_names = c("subj_id", "item_id")) {    
    fac_info <- attr(terms(as.formula(paste0("~", paste(iv_names, collapse = "*")))), "factors")
    rfx <- sapply(unit_names, function(this_unit) {
        ## figure out how many things are replicated by unit, how many times
        rep_mx <- xtabs(paste0("~", paste(iv_names, collapse = "+"), "+", this_unit), dat)
        lvec <- apply(fac_info, 2, function(x) {
            ix <- seq_along(x)[as.logical(x)]
            ## create margin table
            marg_mx <- apply(rep_mx, c(ix, length(dim(rep_mx))), sum)
            mmx <- apply(marg_mx, length(dim(marg_mx)), c)
            ## as.logical(prod(apply(mmx, 2, function(xx) all(xx > 1))))
        })
        ## res <- fac_info[, lvec, drop = FALSE]
        ## keep_term <- rep(TRUE, ncol(res))
        ## try to simplify the formula
    }, simplify = FALSE)
    return(rfx)
}

#' Generate numerical deviation-coded predictors
#'
#' Add deviation-coded predictors to data frame.
#'
#' @param dat A data frame with columns containing factors to be converted.
#' @param iv_names Names of the variables to be converted.
#'
#' @return A data frame including additional deviation coded predictors.
#'
#' @examples
#' with_dev_pred(stim_lists(c(A = 3)), "A")
#' @export
with_dev_pred <- function(dat, iv_names = NULL) {
    if (is.null(iv_names)) {
        iv_names <- names(dat)
    } else {}
    mform <- as.formula(paste0("~", paste(iv_names, collapse = "+")))
    cont <- as.list(rep("contr.dev", length(iv_names)))
    names(cont) <- iv_names
    cbind(dat, model.matrix(mform, dat, contrasts.arg = cont)[, -1])
}

check_design_args <- function(design_args) {
    ## TODO check integrity of design args
    return(TRUE)
}

#' Generate trial lists from a stimulus presentation lists
#'
#' Merge stimulus presentation lists with subject data to create a
#' trial list.
#'
#' @param design_args Stimulus presentation lists (see \code{\link{stim_lists()}}).
#' @param subjects One of the following three: (1) an integer
#' specifying the desired number of subjects (must be a multiple of
#' number of stimulus lists); (2) a data frame with assignment
#' information (must include a column \code{subj_id}); or (3)
#' \code{NULL}, in which case there will be one subject per list.
#' @param seq_assign If TRUE, assignment of subjects to lists will
#' be sequential rather than random (default is FALSE)
#'
#' @return A data frame containing all trial information.
#' @export
trial_lists <- function(design_args,
                        subjects = NULL, seq_assign = FALSE) {
    sp_lists <- stim_lists(design_args)
    sp2 <- split(sp_lists, sp_lists[["list_id"]])
    if (is.null(subjects)) {
        subjects <- length(sp2)
    } else {}
    if (is.numeric(subjects)) {
        if ((subjects %% length(sp2)) != 0) {
            stop("'subjects' must be a multiple of number of lists (",
                 length(sp2), ")")
        } else {}
        list_ord = rep(seq_along(sp2), subjects / length(sp2))
        if (!seq_assign) {
            list_ord = sample(list_ord) # randomize assignment to lists
        } else {}
        subj_dat <- data.frame(subj_id = seq_len(subjects),
                               list_id = list_ord)
    } else {
        if (!is.data.frame(subjects)) {
            stop("'subjects' must be an integer, data.frame, or NULL")
        } else {}
        subj_dat <- subjects
        if (any(!("subj_id" %in% names(subj_dat)),
                !("list_id" %in% names(subj_dat)))) {
            stop("'subjects' must contain fields 'subj_id', 'list_id'")
        } else {}
    }
    res <- lapply(seq_len(nrow(subj_dat)), function(rx) {
               cbind(subj_id = subj_dat[rx, "subj_id"],
                     sp2[[subj_dat[rx, "list_id"]]])
           })
    res2 <- do.call("rbind", res)
    rownames(res2) <- NULL
    return(res2)
}

#' Get names for predictors in factorial design
#'
#' Get the names for the numerical predictors corresponding to all
#' main effects and interactions of categorical IVs in a factorial
#' design.
#'
#' @param design_args A list with experimental design information (see \code{link{stim_lists}})
#' @param design_formula A formula (default NULL, constructs from \code{design_args})
#' @param contr_type Name of formula for generating contrasts (default "contr.dev")
#' @return A character vector with names of all the terms
#' @export 
term_names <- function(design_args,
                       design_formula = NULL,
                       contr_type = "contr.dev") {
    check_design_args(design_args)
    plists <- stim_lists(design_args)
    cont <- as.list(rep("contr.dev", length(design_args[["ivs"]])))
    names(cont) <- names(design_args[["ivs"]])
    if (is.null(design_formula)) design_formula <- as.formula(paste0("~",
                                                               paste(names(design_args[["ivs"]]),
                                                                     collapse = " * ")))
    suppressWarnings(mmx <- model.matrix(design_formula, plists, contrasts.arg = cont))
    return(colnames(mmx))
}

#' Get the GLM formula for a factorially-designed experiment
#'
#' Get the formula corresponding to the general linear model for a
#' factorially designed experiment, with maximal random effects.
#' 
#' @param design_args A list with experimental design information (see \code{link{stim_lists}})
#' @param n_subj Number of subjects
#' @param dv_name Name of dependent variable; \code{NULL} (default) for one-sided formula
#' @param lmer_format Do you want the results combined as the model formula for a \code{lmer} model? (default \code{TRUE})
#' @return A formula, character string, or list, depending
#' @export
design_formula <- function(design_args,
                           n_subj = NULL,
                           dv_name = NULL,
                           lmer_format = TRUE) {
    iv_names <- names(design_args[["ivs"]])

    fixed <- paste(iv_names, collapse = " * ")

    fac_cnts <- fac_counts(iv_names, trial_lists(design_args, subjects = n_subj))
    fac_info <- attr(terms(as.formula(paste0("~", paste(iv_names, collapse = "*")))), "factors")

    rfx <- lapply(fac_cnts, function(lx) {
        lvec <- sapply(lx, function(mx) {
            as.logical(prod(apply(mx, 2, function(xx) all(xx > 1))))
        })
        res <- fac_info[, lvec, drop = FALSE]        
        keep_term <- rep(TRUE, ncol(res))
        ## try to simplify the formula
        for (cx in rev(seq_len(ncol(res))[-1])) {
            drop_term <- sapply(seq_len(cx - 1), function(ccx) {
                identical(as.logical(res[, ccx, drop = FALSE]) | as.logical(res[, cx, drop = FALSE]),
                          as.logical(res[, cx, drop = FALSE]))
            })
            keep_term[seq_len(cx - 1)] <- keep_term[seq_len(cx - 1)] & (!drop_term)
        }
        fterms <- apply(res[, keep_term, drop = FALSE], 2, function(llx) {
            paste(names(llx)[as.logical(llx)], collapse = " * ")
        })
        need_int <- any(apply(lx[[1]], 2, sum) > 1)
        fterms2 <- if (need_int) c("1", fterms) else fterms
        paste(fterms2, collapse = " + ")
    })
    
    form_list <- c(list(fixed = fixed), rfx)
    
    if (lmer_format) {
        form_str <- paste0(dv_name, " ~ ", form_list[["fixed"]], " + ",
               paste(sapply(names(form_list[-1]),
                            function(nx) paste0("(", rfx[[nx]], " | ", nx, ")")),
                     collapse = " + "))
        result <- as.formula(form_str, env = NULL)
    } else {
        result <- lapply(form_list, function(x) formula(paste0("~", x), env = NULL))
    }
    
    return(result)
}

#' Compose data from fixed and random effects, with normally distributed errors
#'
#' @param design_args List containing information about the experiment
#' design; see \code{\link{stim_lists}}
#' @param fixed vector of fixed effects
#' @param subj_rmx matrix of by-subject random effects
#' @param item_rmx matrix of by-item random effects
#' @param err_var error variance
#' @param verbose give debugging info (default = \code{FALSE})
#' @return a data frame containing simulated data
#' @export
compose_data <- function(design_args,
                         n_subj = NULL,
                         fixed = NULL,
                         subj_rmx = NULL, item_rmx = NULL,
                         err_var = 1, verbose = FALSE) {
    ## utility function for doing matrix multiplication
    multiply_mx <- function(des_mx, rfx, row_ix) {
        ## make sure all cols in rfx are represented in des_mx
        diff_cols <- setdiff(colnames(rfx), colnames(des_mx))
        if (length(diff_cols) != 0) {
            stop("column(s) '", paste(diff_cols, collapse = ", "),
                 "' not represented in terms '",
                 paste(term_names(ivs, between_subj, between_item),
                       collapse = ", "), "'")
        } else {}

        reduced_des <- des_mx[, colnames(rfx), drop = FALSE]

        t_rfx <- t(rfx)
        res_vec <- vector("numeric", length(row_ix))
        for (ix in unique(row_ix)) {
            lvec <- row_ix == ix
            res_vec[lvec] <- c(reduced_des[lvec, , drop = FALSE] %*%
                                   t_rfx[, ix, drop = FALSE])
        }
        res_vec
    }

    plists <- stim_lists(ivs, n_item = n_item, n_rep = n_rep)

    tlists <- trial_lists(plists, subjects = n_subj)

    cont <- as.list(rep("contr.dev", length(ivs)))
    names(cont) <- names(ivs)

    mmx <- model.matrix(as.formula(paste0("~", paste(names(ivs), collapse = "*"))),
                        tlists,
                        contrasts.arg = cont)

    if (is.null(fixed)) {
        fixed <- runif(ncol(mmx), -3, 3)
        names(fixed) <- colnames(mmx)
    } else {}

    ## fixed component of Y
    fix_y <- c(mmx %*% fixed) # fixed component of Y

    if (is.null(subj_rmx)) {
        stop("Autogeneration of subj_rmx not implemented yet; please define 'subj_rmx'")
    } else {}
    if (nrow(subj_rmx) != length(unique(tlists[["subj_id"]]))) {
        stop("Argument 'subj_rmx' has ", nrow(subj_rmx), " rows; needs ",
             length(unique(tlists[["subj_id"]])))
    } else {}
    sre <- multiply_mx(mmx, subj_rmx, tlists[["subj_id"]])

    if (is.null(item_rmx)) {
        stop("Autogeneration of item_rmx not implemented yet; please define 'item_rmx'")
    } else {}
    if (nrow(item_rmx) != length(unique(tlists[["item_id"]]))) {
        stop("Argument 'item_rmx' has ", nrow(item_rmx), " rows; needs ",
             length(unique(tlists[["item_id"]])))
    } else {}
    ire <- multiply_mx(mmx, item_rmx, tlists[["item_id"]])
    err <- rnorm(nrow(tlists), sd = sqrt(err_var))
    comb_mx <- matrix(nrow = nrow(tlists), ncol = 0)
    if (verbose) {
        comb_mx <- cbind(fix_y = fix_y, sre = sre, ire = ire, err = err)
    } else {}
    cbind(tlists, Y = fix_y + sre + ire + err, comb_mx)
}
