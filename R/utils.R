# =============================================================================
# utils.R — Shared utility functions
# =============================================================================


#' #' Null-coalescing operator
#' #'
#' #' Returns `x` if non-NULL, otherwise `y`.
#' #'
#' #' @param x,y Values to coalesce.
#' #' @return `x` if non-NULL, else `y`.
#' #' @export
#' `%||%` <- function(x, y) {
#'  if (is.null(x)) y else x
#' }

#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


#' Safe unique with NA handling
#'
#' Returns unique values, optionally removing NAs.
#'
#' @param x Vector.
#' @param na.rm Logical. Remove NAs? Default TRUE.
#' @return Unique values.
#' @export
safe_unique <- function(x, na.rm = TRUE) {
  u <- unique(x)
  if (na.rm) u[!is.na(u)] else u
}


#' Geometric mean
#'
#' @param x Numeric vector (must be positive).
#' @param na.rm Logical. Remove NAs?
#' @return Scalar geometric mean.
#' @keywords internal
geom_mean <- function(x, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  x <- x[x > 0]
  if (length(x) == 0L) return(NA_real_)
  exp(mean(log(x)))
}


#' Model registry: parameter names for each model family
#'
#' Returns the canonical parameter names for a model.
#'
#' @param model Character. Model name.
#' @return Character vector of parameter names.
#' @export
model_params <- function(model) {
  switch(model,
    logistic4    = c("a", "b", "c", "d"),
    loglogistic4 = c("a", "b", "c", "d"),
    gompertz4    = c("a", "b", "c", "d"),
    logistic5    = c("a", "b", "c", "d", "g"),
    loglogistic5 = c("a", "b", "c", "d", "g"),
    stop("Unknown model: ", model)
  )
}


#' List all available model names
#'
#' @return Character vector of the five canonical model names.
#' @export
available_models <- function() {
  c("logistic4", "loglogistic4", "gompertz4", "logistic5", "loglogistic5")
}
