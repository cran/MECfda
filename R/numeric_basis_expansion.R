#' @title Numeric basis expansion for functional variable data
#' @description
#' For a function \eqn{f(t), t\in\Omega}, and a basis function sequence \eqn{\{\rho_k\}_{k\in\kappa}},
#' basis expansion is to compute \eqn{\int_\Omega f(t)\rho_k(t) dt}.
#' Here we do basis expansion for all \eqn{f_i(t), t\in\Omega = [t_0,t_0+T]} in functional variable data, \eqn{i=1,\dots,n}.
#' We compute a matrix \eqn{(b_{ik})_{n\times p}}, where \eqn{b_{ik} = \int_\Omega f(t)\rho_k(t) dt}.
#' The basis we use here is numerically represented by the value of basis functions
#' at some points in the domain.
#' @param object a \code{\link{functional_variable}} class object.
#' The minimum and maximum of the slot \code{t_points} should be respectively
#' equal to the slot \code{t_0} and slot \code{t_0} plus slot \code{period}.
#' @param num_basis a \code{numeric_basis} class object, representing the function basis.
#' See \code{\link{numeric_basis}}.
#' @return Returns a numeric matrix, \eqn{(b_{ik})_{n\times p}},
#' with an extra attribute \code{numeric_basis}, which is the \code{numeric_basis} object input
#' by the argument \code{num_basis}.
#' @export
#' @author Heyang Ji
setGeneric("numeric_basis_expansion",
           function(object,num_basis) standardGeneric("numeric_basis_expansion")

)
#' @rdname numeric_basis_expansion
#' @export
setMethod("numeric_basis_expansion",
          signature(object = "functional_variable",
                    num_basis = "numeric_basis"),
          function(object,num_basis){
            N = nrow(object@X)
            n_t = ncol(object@X)

            if(any(c( object@t_points!= num_basis@t_points,
                      object@t_0!= num_basis@t_0,
                      object@period!= num_basis@period))){
              stop("2 input argument don't match.")
            }

            t_points = object@t_points
            t_0 = object@t_0
            period = object@period


            ip_fpc = object@period * (object@X%*%num_basis@basis_function)/nrow(num_basis@basis_function)
            colnames(ip_fpc) = paste0('basisFunction',1:ncol(num_basis@basis_function))

            attr(ip_fpc,'numeric_basis') = num_basis
            return(ip_fpc)
          })
