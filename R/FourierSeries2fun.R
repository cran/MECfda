#' @title Compute the value of the Fourier summation series
#' @description
#' Compute the value of the Fourier summation series
#' \deqn{f(x) = \frac{a_0}{2} +
#' \sum_{k=1}^{p_a} a_k \cos{(\frac{2\pi}{T}k(x-t_0))} +
#' \sum_{k=1}^{p_b} b_k \sin{(\frac{2\pi}{T}k(x-t_0))},
#' \qquad x\in[t_0,t_0+T]}
#' at some certain point(s).
#'
#' @param object an object of \code{\link{Fourier_series}} class.
#' @param x Value of \eqn{x}.
#'
#' @return A numeric atomic vector
#' @export
#' @author Heyang Ji
#' @examples fsc = Fourier_series(
#'            double_constant = 0.5,
#'            cos = c(0,0.3),
#'            sin = c(1,0.7),
#'            k_cos = 1:2,
#'            )
#'           FourierSeries2fun(fsc,1:5)
#'
setGeneric("FourierSeries2fun",
           function(object,x) standardGeneric("FourierSeries2fun")
           # ,valueClass = "numeric"
)

#' @rdname FourierSeries2fun
#' @export
setMethod("FourierSeries2fun",
          signature(object="Fourier_series",
                    x = "numeric"),

          function(object, x){
            x1 = 2*pi*((object@k_cos) %*% t((x - object@t_0)/object@period))
            x2 = 2*pi*((object@k_sin) %*% t((x - object@t_0)/object@period))
            return(object@double_constant /2 + colSums(object@cos * cos(x1) + object@sin * sin(x2)))
          })
