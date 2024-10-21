#' @title Method of class Fourier_series to extract Fourier coefficients
#'
#' @param object an object of \code{\link{Fourier_series}} class.
#'
#' @return A list that contains the coefficients.
#' @export
#' @author Heyang Ji
#' @examples fsc = Fourier_series(
#'            double_constant = 0.5,
#'            cos = c(0,0.3),
#'            sin = c(1,0.7),
#'            k_cos = 1:2,
#'            )
#'  extractCoef(fsc)
setGeneric("extractCoef",
           function(object) standardGeneric("extractCoef"))
#' @rdname extractCoef
#' @export
setMethod("extractCoef",
          signature(object="Fourier_series"),
          function(object){
            a_0 = object@double_constant
            a_k = object@cos
            b_k = object@sin
            names(a_0) = '0'
            names(a_k) = object@k_cos
            names(b_k) = object@k_sin
            return(list(a_0 = a_0,
                        a_k = a_k,
                        b_k = b_k))
          })
