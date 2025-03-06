#' @title Functional principal component basis expansion for functional variable data
#' @description
#' For a function \eqn{f(t), t\in\Omega}, and a basis function sequence \eqn{\{\rho_k\}_{k\in\kappa}},
#' basis expansion is to compute \eqn{\int_\Omega f(t)\rho_k(t) dt}.
#' Here we do basis expansion for all \eqn{f_i(t), t\in\Omega = [t_0,t_0+T]} in functional variable data, \eqn{i=1,\dots,n}.
#' We compute a matrix \eqn{(b_{ik})_{n\times p}}, where \eqn{b_{ik} = \int_\Omega f(t)\rho_k(t) dt}.
#' The basis we use here is the functional principal component (FPC) basis induced by
#' the covariance function of the functional variable.
#' Suppose \eqn{K(s,t)\in L^2(\Omega\times \Omega)}, \eqn{f(t)\in L^2(\Omega)}.
#' Then \eqn{K} induces an linear operator \eqn{\mathcal{K}} by
#' \deqn{(\mathcal{K}f)(x) = \int_{\Omega} K(t,x)f(t)dt}
#'   If \eqn{\xi(\cdot)\in L^2(\Omega)} such that
#' \deqn{\mathcal{K}\xi = \lambda \xi}
#' where \eqn{\lambda\in {C}},
#' we call \eqn{\xi} an eigenfunction/eigenvector of
#' \eqn{\mathcal{K}}, and \eqn{\lambda} an eigenvalue associated with \eqn{\xi}.\cr
#' For a stochastic process \eqn{\{X(t),t\in\Omega\}}
#' we call the orthogonal basis \eqn{\{\xi_k\}_{k=1}^\infty}
#' corresponding to eigenvalues \eqn{\{\lambda_k\}_{k=1}^\infty}
#' (\eqn{\lambda_1\geq\lambda_2\geq\dots}),
#' induced by
#' \deqn{K(s,t)=\text{Cov}(X(t),X(s))}
#' a functional principal component (FPC) basis for \eqn{L^2(\Omega)}.
#' @param object a \code{\link{functional_variable}} class object.
#' The minimum and maximum of the slot \code{t_points} should be respectively
#' equal to the slot \code{t_0} and slot \code{t_0} plus slot \code{period}.
#' @param npc The number of functional principal components. See \code{npc} in \code{\link[refund]{fpca.sc}}.
#' @return Returns a numeric matrix, \eqn{(b_{ik})_{n\times p}},
#' with an extra attribute \code{numeric_basis}, which represents the FPC basis.
#' The attribute \code{numeric_basis} is a \code{numeric_basis} object. See \code{\link{numeric_basis}}.
#' The slot \code{basis_function} is also a numeric matrix, denoted as \eqn{(\zeta_{jk})_{m\times p}}
#' \deqn{b_{ik} = \int_\Omega f(t)\xi_k(t) dt}
#' \deqn{\zeta_{jk} = \xi_k(t_j)}
#' @export
#' @author Heyang Ji
#' @importFrom refund fpca.sc
#' @examples
#' n<-50; ti<-seq(0,1,length.out=101)
#' X<-t(sin(2*pi*ti)%*%t(rnorm(n,0,1)))
#' object = functional_variable(X = X, t_0 = 0, period = 1, t_points = ti)
#' a = FPC_basis_expansion(object,3L)
#' dim(a)

setGeneric("FPC_basis_expansion",
           function(object,npc) standardGeneric("FPC_basis_expansion")

)
#' @rdname FPC_basis_expansion
#' @export
setMethod("FPC_basis_expansion",
          signature(object="functional_variable",
                    npc = "integer"),
          function(object,npc){
            N = nrow(object@X)
            n_t = ncol(object@X)
            t_points = object@t_points
            t_0 = object@t_0
            period = object@period

            if(t_0!=min(t_points) | (t_0+period)!=max(t_points)){
              stop("The minimum and maximum and  of t_points should be t_0 and t_0 + period.")
            }


            fpca<- refund::fpca.sc(object@X, argvals=t_points,var=T,center=T,npc=npc)
            eigenfunctions<-fpca$efunctions
            ip_fpc = object@period * (object@X%*%eigenfunctions)/nrow(eigenfunctions)
            colnames(ip_fpc) = paste0('FPC',1:npc)
            colnames(eigenfunctions) = paste0('FPC',1:npc)

            FPC_basis = numeric_basis(basis_function = eigenfunctions,
                                      t_points = t_points)
            attr(ip_fpc,'numeric_basis') = FPC_basis
            # attributes(ip_fpc)$numeric_basis
            # attr(ip_fpc,'numeric_basis')
            return(ip_fpc)
          })

