#' @method plot numericBasis_series
#' @title Plot numeric basis function summation series.
#'
#' @param x A \code{\link{numericBasis_series}} object.
#'
#' @author Heyang Ji
#' @export
#' @examples
#' t_0 = 0
#' period = 1
#' t_points = seq(0.05,0.95,length.out = 19)
#' nb = numeric_basis(
#'   basis_function = cbind(1/2,cos(2*pi*t_points),sin(2*pi*t_points)),
#'   t_points       = t_points,
#'   t_0            = t_0,
#'   period         = period
#' )
#' ns = numericBasis_series(coef = c(0.8,1.2,1.6),numeric_basis = nb)
#' plot(ns)

setMethod("plot",
          signature(x="numericBasis_series"),
          function(x){
            t = x@numeric_basis@t_0 + x@numeric_basis@period*(0:1000/1000)
            plot(x = t, y = numericBasisSeries2fun(x,t),
                 xlab = "Domain",
                 ylab = "Series value",
                 main = "Numeric Basis Function Series")
          })
# plot.numericBasis_series = function(x){
#   t = x@numeric_basis@t_0 + x@numeric_basis@period*(0:1000/1000)
#   plot(x = t, y = numericBasisSeries2fun(x,t),
#        xlab = "Domain",
#        ylab = "Series value",
#        main = "Numeric Basis Function Series")
# }
