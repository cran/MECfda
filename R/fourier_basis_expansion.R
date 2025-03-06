#' @title Fourier basis expansion for functional variable data
#' @description
#' For a function \eqn{f(x), x\in\Omega}, and a basis function sequence \eqn{\{\rho_k\}_{k\in\kappa}},
#' basis expansion is to compute \eqn{\int_\Omega f(t)\rho_k(t) dt}.
#' Here we do basis expansion for all \eqn{f_i(t), t\in\Omega = [t_0,t_0+T]} in functional variable data, \eqn{i=1,\dots,n}.
#' We compute a matrix \eqn{(b_{ik})_{n\times p}}, where \eqn{b_{ik} = \int_\Omega f(t)\rho_k(t) dt}.
#' The basis used here is the Fourier basis, \deqn{\frac{1}{2},\ \cos(\frac{2\pi}{T}k[x-t_0]),\ \sin (\frac{2\pi}{T}k[x-t_0])}
#' where \eqn{x\in[t_0,t_0+T]} and \eqn{k = 1,\dots,p_f}.
#' @param object a \code{\link{functional_variable}} class object.
#' @param order_fourier_basis the order of Fourier basis, \eqn{p_f}.
#'
#' @return Returns a numeric matrix, \eqn{(b_{ik})_{n\times p}}, where \eqn{b_{ik} = \int_\Omega f(t)\rho_k(t) dt}.
#' @export
#' @author Heyang Ji
setGeneric("fourier_basis_expansion",
           function(object,order_fourier_basis) standardGeneric("fourier_basis_expansion")

)
#' @rdname fourier_basis_expansion
#' @export
setMethod("fourier_basis_expansion",
          signature(object="functional_variable",
                    order_fourier_basis = "integer"),
          function(object,order_fourier_basis){

            N = nrow(object@X)
            n_t = ncol(object@X)
            t_points = object@t_points



            {
              t_fourier = 2*pi*(((t_points - object@t_0)/object@period)%*%t(1:order_fourier_basis))
              a_k = object@period * (object@X%*%cos(t_fourier))/nrow(t_fourier)
              b_k = object@period * (object@X%*%sin(t_fourier))/nrow(t_fourier)
            }


            a_0 = object@period * rowMeans(object@X)/2
            colnames(a_k) = paste0('cos',1:ncol(a_k))
            colnames(b_k) = paste0('sin',1:ncol(b_k))
            return(cbind(a_0,a_k,b_k))
          })
