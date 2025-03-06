#' @title Compute the value of the basis function summation series at certain points.
#' @description Compute the function \eqn{f(x) = \sum_{k=0}^{p}c_k \rho_{k}(x)}, \eqn{x\in\Omega}
#'
#' @param object an object of \code{\link{numericBasis_series}} class.
#' @param x Value of $x$.
#'
#' @return A numeric atomic vector
#' @export
#' @author Heyang Ji
#' @importFrom pracma interp1
#'
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
#' numericBasisSeries2fun(ns,seq(0,1,length.out = 51))
#'
setGeneric("numericBasisSeries2fun",
           function(object,x) standardGeneric("numericBasisSeries2fun")
)
#' @rdname numericBasisSeries2fun
#' @export
setMethod("numericBasisSeries2fun",
          signature(object="numericBasis_series",
                    x = "numeric"),
          function(object,x){
            num_fun = as.vector(object@numeric_basis@basis_function%*%as.matrix(object@coef))
            t_points = object@numeric_basis@t_points

            if(min(x)<object@numeric_basis@t_0 | max(x)>object@numeric_basis@t_0+object@numeric_basis@period){
              stop("x exceeds the value range.")
            }

            range0 = range(t_points)
            which0 = which(x %in% t_points)
            which1 = setdiff(which(x > range0[1] & x < range0[2]),which0)
            which2 = setdiff(which(x < range0[1]                ),which0)
            which3 = setdiff(which(                x > range0[2]),which0)

            ret = numeric(length(x))
            if(length(which0>0)){
              ret[which0] = num_fun[match(x[which0],t_points)]
            }
            if(length(which1>0)){
              ## cubic Hermite interpolation
              ret[which1] = pracma::interp1(
                x = t_points,
                y = num_fun,
                xi = x[which1],
                method = 'cubic')
            }
            ## use the two points on the edge for extrapolation
            if(length(which2>0)){
              which_lin_extra_points = order(t_points)[1:2]
              ys = num_fun[which_lin_extra_points]
              xs = t_points[which_lin_extra_points]
              ret[which2] = predict(lm(ys~xs),newdata = data.frame(xs = x[which2]))
            }
            if(length(which3>0)){
              which_lin_extra_points = order(t_points,decreasing = T)[1:2]
              ys = num_fun[which_lin_extra_points]
              xs = t_points[which_lin_extra_points]
              ret[which3] = predict(lm(ys~xs),newdata = data.frame(xs = x[which3]))
            }
            return(ret)
          })
