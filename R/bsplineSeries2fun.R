#' @title Compute the value of the Fourier summation series at certain points.
#'
#' @param object an object of \code{\link{bspline_series}} class.
#' @param x Value of $x$.
#'
#' @return A numeric atomic vector
#' @export
#' @author Heyang Ji
#' @importFrom splines bs
#' @examples bsb = bspline_basis(
#'             Boundary.knots = c(0,24),
#'             intercept      = TRUE,
#'             df             = NULL,
#'             degree         = 3
#' )
#' bss = bspline_series(
#'           coef = c(2,1,1.5,0.5),
#'           bspline_basis = bsb
#' )
#' bsplineSeries2fun(bss,(1:239)/10)
setGeneric("bsplineSeries2fun",
           function(object,x) standardGeneric("bsplineSeries2fun")

)
#' @rdname bsplineSeries2fun
#' @export
setMethod("bsplineSeries2fun",
          signature(object="bspline_series",
                    x = "numeric"),
          function(object,x){
            bs_basis_val = splines::bs(x,
                                       df             = object@bspline_basis@df,
                                       knots          = object@bspline_basis@knots,
                                       degree         = object@bspline_basis@degree,
                                       intercept      = object@bspline_basis@intercept,
                                       Boundary.knots = object@bspline_basis@Boundary.knots)
            as.vector(bs_basis_val%*%as.matrix(object@coef))
          })
