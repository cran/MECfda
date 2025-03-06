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
#' default value is \code{TRUE}, and must be \code{TRUE}. See \code{intercept} \code{\link[splines]{bs}}.
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

#' @title Numeric representation of a function basis
#' @description
#' A s4 class that numerically represents a basis of linear space of function. \cr
#' \eqn{\{\rho_k\}_{k=1}^\infty} denotes a basis of function linear space.
#' Some times the basis cannot be expressed analytically.
#' But we can numerically store the space by the value of
#' a finite subset of the basis functions at some certain points in the domain,
#' \eqn{\rho_k(t_j), k = 1,\dots,p, j = 1,\dots,m}.
#' The s4 class is to represent a finite sequence of functions by their values
#' at a finite sequence of points within their domain,
#' in which all the functions have the same domain and the domain is an interval.
#'
#' @slot basis_function matrix of the value of the functions,
#' \eqn{(\zeta_{jk})_{m\times p}},
#' where \eqn{\zeta_{ik} = \rho_k(t_j), j = 1,\dots,m, k = 1,\dots,p}.
#' Each row of the matrix is corresponding to a point of \eqn{t}.
#' Each column of the matrix is corresponding to a basis function.
#' @slot t_points a numeric atomic vector,
#' represents the points in the domains of the function
#' where the function values are taken.
#' The \eqn{j}th element is corresponding to \eqn{j}th row of
#' slot \code{basis_function}.
#' @slot t_0 left end of the domain interval.
#' @slot period length of the domain interval.
#' @details
#' The units of a basis of a linear space should be linearly independent.
#' But the program doesn't check the linear dependency of the basis function
#' when a \code{numeric_basis} object is initialized.
#'
#' @author Heyang Ji
#' @export
#' @import methods
#' @examples
#' t_0 = 0
#' period = 1
#' t_points = seq(0.05,0.95,length.out = 19)
#' numeric_basis(
#'   basis_function = cbind(1/2,cos(t_points),sin(t_points)),
#'   t_points       = t_points,
#'   t_0            = t_0,
#'   period         = period
#' )

numeric_basis = setClass(
  "numeric_basis",
  slots = c(
    basis_function = "matrix",
    ## matrix of the value of basis functions
    ## each column represents a basis function
    ## each row represents a value of $t$
    t_points       = "numeric",
    t_0            = "numeric",
    period         = "numeric"
  )
)


#' @title Linear combination of a sequence of basis functions represented numerically
#' @description
#' A linear combination of basis function \eqn{\{\rho_k\}_{k=1}^p},
#' \deqn{\sum_{k=1}^p c_k \rho_k(t).}
#'
#' @slot coef linear coefficient \eqn{\{c_k\}_{k=1}^p}.
#' @slot numeric_basis \eqn{\{\rho_k\}_{k=1}^p} represented by a \code{numeric_basis} object.
#' See \code{\link{numeric_basis}}.
#'
#' @author Heyang Ji
#' @export
#' @import methods
#' @examples
#' t_0 = 0
#' period = 1
#' t_points = seq(0.05,0.95,length.out = 19)
#' nb = numeric_basis(
#'   basis_function = cbind(1,cos(t_points),sin(t_points)),
#'   t_points       = t_points,
#'   t_0            = t_0,
#'   period         = period
#' )
#' ns = numericBasis_series(coef = c(0.8,1.2,1.6),numeric_basis = nb)
numericBasis_series = setClass(
  "numericBasis_series",
  slots = c(
    coef = "numeric",
    numeric_basis  = "numeric_basis"
  )
)
