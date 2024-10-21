#' @title From the summation series of a functional basis to function value
#' @description
#' Generic function to compute function value from summation series of a functional basis.
#' @param object An object that represents a functional basis.
#' @param x point(s) to take value.
#' @return A numeric atomic vector.
#' See \code{\link{bsplineSeries2fun}} and \code{\link{FourierSeries2fun}}.
#'
#' @export
#' @author Heyang Ji
#' @details
#' When applied to \code{\link{bspline_series}} object, equivalent to \code{\link{bsplineSeries2fun}}.\cr
#' When applied to \code{\link{Fourier_series}} object, equivalent to \code{\link{FourierSeries2fun}}.
basis2fun = function(object,x) UseMethod("basis2fun")
#' @rdname basis2fun
#' @export
# basis2fun.bspline_series = function(object,x) bsplineSeries2fun(object,x)
setMethod("basis2fun",
          signature(object="bspline_series",
                    x = "numeric"),
          function(object,x) bsplineSeries2fun(object,x)
)

#' @rdname basis2fun
#' @export
# basis2fun.Fourier_series = function(object,x) FourierSeries2fun(object,x)
setMethod("basis2fun",
          signature(object="Fourier_series",
                    x = "numeric"),
          function(object,x) FourierSeries2fun(object,x)
)
