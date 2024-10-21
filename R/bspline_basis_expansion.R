#' @title B-splines basis expansion for functional variable data
#' @description
#' For a function \eqn{f(t), t\in\Omega}, and a basis function sequence \eqn{\{\rho_k\}_{k\in\kappa}},
#' basis expansion is to compute \eqn{\int_\Omega f(t)\rho_k(t) dt}.
#' Here we do basis expansion for all \eqn{f_i(t), t\in\Omega = [t_0,t_0+T]} in functional variable data, \eqn{i=1,\dots,n}.
#' We compute a matrix \eqn{(b_{ik})_{n\times p}}, where \eqn{b_{ik} = \int_\Omega f(t)\rho_k(t) dt}.
#' The basis used here is the b-splines basis, \eqn{\{B_{i,p}(x)\}_{i=-p}^{k}}, \eqn{x\in[t_0,t_{k+1}]},
#' where \eqn{t_{k+1} = t_0+T} and \eqn{B_{i,p}(x)} is defined as
#' \deqn{B_{i,0}(x) = \left\{
#'   \begin{aligned}
#'   &I_{(t_i,t_{i+1}]}(x), & i = 0,1,\dots,k\\
#'   &0, &i<0\ or\ i>k
#'   \end{aligned}
#'   \right.}
#'  \deqn{B_{i,r}(x) = \frac{x - t_{i}}{t_{i+r}-t_{i}} B_{i,r-1}(x) + \frac{t_{i+r+1} - x}
#'     {t_{i+r+1} - t_{i+1}}B_{i+1,r-1}(x)}
#' @param object a \code{\link{functional_variable}} class object.
#' @param n_splines the number of splines, equal to \eqn{k+p+1}. See \code{df} in \code{\link[splines]{bs}}.
#' @param bs_degree the degree of the piecewise polynomial of the b-splines. See \code{degree} in \code{\link[splines]{bs}}.
#'
#' @return Returns a numeric matrix, \eqn{(b_{ik})_{n\times p}}, where \eqn{b_{ik} = \int_\Omega f(t)\rho_k(t) dt}
#' @export
#' @author Heyang Ji
#' @importFrom splines bs
setGeneric("bspline_basis_expansion",
           function(object,n_splines,bs_degree) standardGeneric("bspline_basis_expansion")

)
#' @rdname bspline_basis_expansion
#' @export
setMethod("bspline_basis_expansion",
          signature(object="functional_variable",
                    n_splines = "integer"),
          function(object,n_splines,bs_degree){
            N = nrow(object@X)
            n_t = ncol(object@X)
            t_points = object@t_points

            bs_basis = splines::bs(t_points,
                                   df = n_splines,
                                   degree = bs_degree,
                                   intercept = TRUE,
                                   Boundary.knots = c(object@t_0, object@t_0 + object@period))

            ip_bs = (object@X%*%bs_basis)/nrow(bs_basis)
            colnames(ip_bs) = paste0('bs',1:n_splines)
            return(ip_bs)
          })
