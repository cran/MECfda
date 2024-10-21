#' @method plot bspline_series
#' @title Plot b-splines baisi summation series.
#'
#' @param x A \code{\link{bspline_series}} object.
#' @return No return value. Generate a scatter plot.
#'
#' @author Heyang Ji
#' @export
#' @examples
#' bsb = bspline_basis(
#' Boundary.knots = c(0,24),
#' intercept      = TRUE,
#' df             = NULL,
#' degree         = 3
#' )
#' bss = bspline_series(
#' coef = c(2,1,1.5,3),
#' bspline_basis = bsb
#' )
#' plot(bss)

setMethod("plot",
          signature(x="bspline_series"),
          function(x){
            t = x@bspline_basis@Boundary.knots[1] + x@bspline_basis@Boundary.knots[2]*(0:1000/1000)
            plot(x = t, y = bsplineSeries2fun(x,t),
                 xlab = "Domain",
                 ylab = "Series value",
                 main = "B-splines Series")
          })
# plot.bspline_series = function(x){
#   t = x@bspline_basis@Boundary.knots[1] + x@bspline_basis@Boundary.knots[2]*(0:1000/1000)
#   plot(x = t, y = bsplineSeries2fun(x,t),
#        xlab = "x",
#        ylab = "Series value",
#        main = "Curve of the B-splines Series within a period")
# }
