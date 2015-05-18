#' Generate stimulus presentation lists
#'
#' Generates counterbalanced presentation lists for factorially
#' designed experiments involving stimulus presentation
#'
#' @param design_args A list describing the experimental design, which
#' must have an element \code{ivs}, giving a named list of independent
#' variables, with each list element a vector giving the levels of
#' that IV, or a single number stating the number of desired levels.
#' If any IVs are administered between-subjects or between-items,
#' these should be named in list elements \code{between_subj} and
#' \code{between_item} respectively.  The argument \code{design_args}
#' also can optionally include the following two elements:
#' \code{n_item}, the desired number of stimulus items, which if
#' unspecified, will result in lists with the minimum possible number
#' of items; and \code{n_rep}, the number of repetitions of each
#' stimulus item for each participant (default 1).
#' @param as_one boolean (default \code{TRUE}) specifying whether the
#' presentation lists are to be returned as a single data frame or as
#' elements in a list object
#'
#' @return a single \code{data.frame} (default) with each list
#' identified by \code{list_id} or a \code{list} of dataframes,
#' depending on the value of \code{as_one}
#'
#' @examples
#' stim_lists(list(ivs = c(A = 2, B = 2))) # 2x2 within-subjects within-item
#' 
#' stim_lists(list(c(A = 2, B = 2), n_item = 16)) # same but w/more items
#'
#' stim_lists(list(c(A = 2, B = 2), n_item = 16, n_rep = 3))
#'
#' # mixed by subject, fully within by item
#' stim_lists(list(ivs = list(group = c("adult", "child"),
#'                            task = c("easy", "hard")),
#'            between_subj = "group",
#'            n_item = 12))
#'
#' # mixed by subject, mixed by item
#' stim_lists(list(ivs = c(A = 2, B = 2),
#'            between_subj = "A",
#'            between_item = "B"))
#' @export
stim_lists <- function(design_args, 
                       as_one = TRUE) {
    fac_combine_levels <- function(vars, iv_list, dframe = TRUE) {
        row_indices <- rev(do.call("expand.grid",
                                   lapply(rev(vars),
                                          function(x) seq_along(iv_list[[x]]))))
        res <- mapply(function(x, y) x[y], iv_list[vars], row_indices)
        if (dframe) {
            as.data.frame(res, stringsAsFactors = FALSE)
        } else {
            res
        }
    }

    rotate_cells <- function(x, combine = FALSE) {
        res <- lapply(seq_len(nrow(x)),
                      function(ix) x[c(ix:nrow(x), seq_len(ix - 1)), , drop = FALSE])
        if (combine) {
            do.call("rbind", res)
        } else {
            res
        }
    }

    bs_combine <- function(dat, plists) {
        ## factorially combine across lists for between subject variables
        if (nrow(dat) == 0) {
            return(plists)
        } else {}
        res <- c(lapply(seq_len(nrow(dat)), function(rx) {
            if (length(plists) > 0) {
                lapply(plists, function(lx) {
                    cbind(dat[rep(rx, nrow(lx)), , drop = FALSE], lx)
                })
            } else {
                dat[rx, , drop = FALSE]
            }
        }))
        if (length(plists) > 0) {
            do.call("c", res)
        } else {
            res
        }
    }

    check_design_args(design_args)
    iv_names <- names(design_args[["ivs"]])
    ## check whether elements of 'ivs' are numbers and convert to char vector
    ivs2 <- lapply(iv_names, function(nx) {
        x <- design_args[["ivs"]][[nx]]
        if ((length(x) == 1) && is.numeric(x)) {
            paste0(nx, seq_len(x))
        } else {x}
    })
    names(ivs2) <- iv_names
    
    item_within <- setdiff(names(ivs2), design_args[["between_item"]])
    subj_within <- setdiff(names(ivs2), design_args[["between_subj"]])
    ## iv_levels <- sapply(ivs2, length) # IS THIS NEEDED?

    ww_fac <- intersect(item_within, subj_within)

    ww_chunks <- rotate_cells(fac_combine_levels(intersect(item_within, subj_within),
                                                 ivs2))

    wb_chunks <- fac_combine_levels(intersect(subj_within, design_args[["between_item"]]), ivs2)

    ## combine the WSWI and WSBI chunks to create the base presentation lists
    if (nrow(wb_chunks) > 0) {
        if (length(ww_chunks) > 0) {
            base_plists <- lapply(ww_chunks, function(ww) {
                cbind(wb_chunks[rep(seq_len(nrow(wb_chunks)), each = nrow(ww)), , drop = FALSE], ww)
            })
        } else {
            base_plists <- list(wb_chunks)
        }
    } else {
        base_plists <- ww_chunks
    }

    ## handle BSWI
    bswi <- fac_combine_levels(intersect(design_args[["between_subj"]],
                                         item_within), ivs2)
    bswi_lists <- bs_combine(bswi, base_plists)

    ## handle bsbi factors (if they exist)
    bsbi <- fac_combine_levels(intersect(design_args[["between_subj"]],
                                         design_args[["between_item"]]), ivs2)
    bsbi_lists <- bs_combine(bsbi, bswi_lists)
    div_fac <- if (nrow(bsbi)) nrow(bsbi) else 1
    n_item <- design_args[["n_item"]]
    if (is.null(design_args[["n_item"]])) { # dynamically choose minimum n_item
        if (length(bswi_lists) > 0) {
            n_item <- nrow(bswi_lists[[1]]) * div_fac
        } else {
            n_item <- div_fac
        }
    } else {}
    if (length(bswi_lists) > 0) {
        item_fac <- div_fac * nrow(bswi_lists[[1]])
        if ((n_item %% item_fac) != 0) {
            stop("n_item must be a factor of ", item_fac)
        } else {}
    } else {}
    
    if (length(bsbi_lists) == 0) {
        bsbi_lists <- bswi_lists
    } else {}

    if ((n_item %% div_fac) != 0) stop("n_item must be a multiple of ", div_fac)
    
    rep_times <- if (length(bswi_lists) > 0) length(bswi_lists) else 1
    it_chunks <- rep(seq_len(div_fac), each = rep_times)
    it_lists <- split(seq_len(n_item),
                      rep(seq_len(div_fac), each = n_item / div_fac))[it_chunks]

    plists <- mapply(function(x, y) {
        ix <- rep(seq_len(nrow(y)), each = length(x) / nrow(y))
        cbind(item_id = x, y[ix, , drop = FALSE])
    },
                     it_lists, bsbi_lists, SIMPLIFY = FALSE)

    n_rep <- design_args[["n_rep"]]
    if (is.null(design_args[["n_rep"]])) n_rep <- 1

    if (n_rep > 1) {
        plists <- lapply(plists, function(x) {
                             cbind(n_rep = paste0("r", rep(seq_len(n_rep), each = nrow(x))),
                                   x[rep(seq_len(nrow(x)), n_rep), , drop = FALSE])
                         })
    } else {}
    
    if (as_one) {
        res <- mapply(function(x, y) {
            cbind(list_id = x, y)
        },
                      seq_along(plists), plists, SIMPLIFY = FALSE)
        final_lists <- do.call("rbind", res)
    } else {
        final_lists <- plists
    }
    rownames(final_lists) <- NULL
    return(final_lists)
}
