#' Derive names for a given design
#'
#' Derive names of predictor terms for a given design
#' 
#' @param ivs named list of independent variables (IVs) with each list
#' element a vector giving the levels of that IV, or a single number
#' stating the number of desired levels
#' @param between_subj character vector with names of IVs whose
#' levels are administered between subjects
#' @param between_item charater vector with names of IVs whose levels
#' are administered between items
#' @return A vector of term names
#' @seealso \code{\link{stim_lists}}
#' @export
term_names <- function(ivs, between_subj = c(), between_item = c()) {
    plists <- stim_lists(ivs, between_subj, between_item)
    cont <- as.list(rep("contr.dev", length(ivs)))
    names(cont) <- names(ivs)
    mmx <- model.matrix(as.formula(paste0("~", paste(names(ivs), collapse = "*"))),
                        plists,
                        contrasts.arg = cont)
    colnames(mmx)
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

#' Generate trial lists from a stimulus presentation lists
#'
#' Merge stimulus presentation lists with subject data to create a
#' trial list.
#'
#' @param sp_lists Stimulus presentation lists (see \code{\link{stim_lists()}}).
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
trial_lists <- function(sp_lists, subjects = NULL, seq_assign = FALSE) {
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

#' Simulate data with normally distributed errors
#'
#' @param ivs named list of independent variables (IVs) with each list
#' element a vector giving the levels of that IV, or a single number
#' stating the number of desired levels
#' @param between_subj character vector with names of IVs whose
#' levels are administered between subjects
#' @param between_item charater vector with names of IVs whose levels
#' are administered between items
#' @param n_item desired number of stimulus items; if \code{NULL},
#' will generate lists with the minimum possible number
#' @param n_rep number of repetitions of each stimulus item for each
#' participant (default 1)
#' @param fixed vector of fixed effects
#' @param subj_rmx matrix of by-subject random effects
#' @param item_rmx matrix of by-item random effects
#' @param err_var error variance
#' @param verbose give debugging info (default = \code{FALSE})
#' @return a data frame containing simulated data
#' @export
sim_norm <- function(ivs,
                     between_subj = c(),
                     between_item = c(),
                     n_subj = NULL, n_item = NULL, n_rep = 1,
                     fixed = NULL, subj_rmx = NULL, item_rmx = NULL,
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
    sre <- multiply_mx(mmx, subj_rmx, tlists[["subj_id"]])

    if (is.null(item_rmx)) {
        stop("Autogeneration of item_rmx not implemented yet; please define 'item_rmx'")
    } else {}
    ire <- multiply_mx(mmx, item_rmx, tlists[["item_id"]])
    err <- rnorm(nrow(tlists), sd = sqrt(err_var))
    comb_mx <- matrix(nrow = nrow(tlists), ncol = 0)
    if (verbose) {
        comb_mx <- cbind(fix_y = fix_y, sre = sre, ire = ire, err = err)
    } else {}
    cbind(tlists, Y = fix_y + sre + ire + err, comb_mx)
}
