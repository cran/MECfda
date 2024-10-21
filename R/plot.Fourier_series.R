#' @method plot Fourier_series
#' @title Plot Fourier basis summation series.
#'
#' @param x A \code{\link{Fourier_series}} object.
#' @return No return value. Generate a scatter plot.
#'
#' @author Heyang Ji
#' @export
#' @examples
#' fsc = Fourier_series(
#' double_constant = 0.5,
#' cos = c(0,0.3),
#' sin = c(1,0.7),
#' k_cos = 1:2,
#' )
#' plot(fsc)

setMethod("plot",
          signature(x="Fourier_series"),
          function(x){
            t = x@t_0 + x@period*(0:1000/1000)
            plot(x = t, y = FourierSeries2fun(x,t),
                 xlab = "Domain",
                 ylab = "Series value",
                 main = "Fourier Series")
          })
# plot.Fourier_series = function(x){
#   t = x@t_0 + x@period*(0:1000/1000)
#   plot(x = t, y = FourierSeries2fun(x,t),
#        xlab = "x",
#        ylab = "Series value",
#        main = "Curve of the Fourier Series within a period")
# }
