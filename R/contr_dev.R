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
    ctreat <- contr.treatment(n, base, contrasts)
    mx <- apply(ctreat, 2, scale, scale = FALSE)
    dimnames(mx) <- dimnames(ctreat)
    mx
}
