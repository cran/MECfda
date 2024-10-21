setClassUnion('num.int.nul',c("numeric","integer","NULL"))

#' @title b-spline basis
#' @description
#' A s4 class that represents a b-spline basis \eqn{\{B_{i,p}(x)\}_{i=-p}^{k}} on the interval \eqn{[t_0,t_{k+1}]},
#' where \eqn{B_{i,p}(x)} is defined as
#' \deqn{B_{i,0}(x) = \left\{
#'   \begin{aligned}
#'   &I_{(t_i,t_{i+1}]}(x), & i = 0,1,\dots,k\\
#'   &0, &i<0\ or\ i>k
#'   \end{aligned}
#'   \right.}
#'  \deqn{B_{i,r}(x) = \frac{x - t_{i}}{t_{i+r}-t_{i}} B_{i,r-1}(x) + \frac{t_{i+r+1} - x}
#'     {t_{i+r+1} - t_{i+1}}B_{i+1,r-1}(x)}
#'       For all the discontinuity points of \eqn{B_{i,r}} (\eqn{r>0}) in the interval \eqn{(t_0,t_k)},
#'       let the value equals its limit, which means
#'  \deqn{B_{i,r}(x) = \lim_{t\to x} B_{i,r}(t)}
#'
#' @slot Boundary.knots boundary of the domain of the splines (start and end), which is \eqn{t_0} and \eqn{t_{k+1}}.
#' Default is \eqn{[0,1]}. See \code{Boundary.knots} in \code{\link[splines]{bs}}.
#' @slot knots knots of the splines, which is \eqn{(t_1,\dots,t_k)},
#' equally spaced sequence is chosen by the function automatically with equal space
#' (\eqn{t_j = t_0 + j\cdot\frac{t_{k+1}-t_0}{k+1}}) when not assigned.
#' See \code{knots} in \code{\link[splines]{bs}}.
#' @slot intercept Whether an intercept is included in the basis,
#' default value is TRUE, and must be TRUE. See \code{intercept} \code{\link[splines]{bs}}.
#' @slot df degree of freedom of the basis, which is the number of the splines, equal to \eqn{p+k+1}.
#' By default \eqn{k = 0}, and \code{df} \eqn{= p+1}. See \code{df} \code{\link[splines]{bs}}.
#' @slot degree degree of the splines, which is the degree of piecewise polynomials \eqn{p}, default value is 3.
#' See \code{degree} in \code{\link[splines]{bs}}.
#'
#' @export
#' @import methods
#' @author Heyang Ji
#' @examples bsb = bspline_basis(
#'             Boundary.knots = c(0,24),
#'             intercept      = TRUE,
#'             df             = NULL,
#'             degree         = 3
#' )
bspline_basis = setClass(
  "bspline_basis",
  slots = c(
    Boundary.knots = "numeric",
    knots          = "num.int.nul",
    # knots          = "ANY",
    intercept      = "logical",
    df             = "num.int.nul",
    # df             = "ANY",
    degree         = "integer"
  )
)

#' @title b-splines summation series.
#' @description
#' A s4 class that represents
#' the summation \eqn{\sum_{i=0}^{k}b_i B_{i,p}(x)} by a bspline_basis object
#' and coefficients \eqn{b_i} (\eqn{i = 0,\dots,k}).
#'
#' @slot coef coefficients of the b-splines, \eqn{b_i} (\eqn{i = 0,\dots,k}).
#' @slot bspline_basis a \code{\link{bspline_basis}} object,
#' represents the b-splines basis used, \eqn{\{B_{i,p}(x)\}_{i=-p}^{k}}.
#'
#' @export
#' @import methods
#' @author Heyang Ji
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
bspline_series = setClass(
  "bspline_series",
  slots = c(
    coef = "numeric",
    bspline_basis = "bspline_basis")
)

#' @title s4 class of Fourier summation series
#' @description
#' A s4 class that represents the linear combination of Fourier basis functions below:
#' \deqn{\frac{a_0}{2} +
#' \sum_{k=1}^{p_a} a_k \cos{(\frac{2\pi}{T}k(x-t_0))} +
#' \sum_{k=1}^{p_b} b_k \sin{(\frac{2\pi}{T}k(x-t_0))},
#' \qquad x\in[t_0,t_0+T]}
#'
#' @slot double_constant value of \eqn{a_0}.
#' @slot cos values of coefficients of \eqn{\cos} waves, \eqn{a_k}.
#' @slot sin values of coefficients of \eqn{\sin} waves, \eqn{b_k}.
#' @slot k_cos values of \eqn{k} corresponding to the coefficients of \eqn{\cos} waves
#' @slot k_sin values of \eqn{k} corresponding to the coefficients of \eqn{\sin} waves
#' @slot t_0 left end of the domain interval, \eqn{t_0}
#' @slot period length of the domain interval, \eqn{T}.
#' @details
#' If not assigned, \eqn{t_0 = 0}, \eqn{T = 2\pi}.
#' If not assigned, k_cos and k_sin equals 1, 2, 3, ...
#' @export
#' @import methods
#' @author Heyang Ji
#' @examples fsc = Fourier_series(
#'            double_constant = 0.5,
#'            cos = c(0,0.3),
#'            sin = c(1,0.7),
#'            k_cos = 1:2,
#'            )
Fourier_series = setClass(
  "Fourier_series",
  slots = c(
    double_constant = "numeric",
    cos             = "numeric",
    sin             = "numeric",
    k_cos           = "num.int.nul",
    k_sin           = "num.int.nul",
    # k_cos           = "ANY",
    # k_sin           = "ANY",
    t_0             = "numeric",
    period          = "numeric"
  )
)
