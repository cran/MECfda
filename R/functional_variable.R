#' @title Function-valued variable data.
#' @description
#' A s4 class that represents data of a function-valued variable.
#' The format is
#'   \eqn{f_i(t),\ t\in\Omega=[t_0,t_0 + T]}
#'   where \eqn{i} is the observation (subject) index, \eqn{t} represents the measurement (time) points.
#'
#' @slot X a matrix \eqn{(x_{ij})_{n\times m}}, where \eqn{x_{ij} = f_i(t_j)}, represents the value of \eqn{f_i(t_j)}, each row represent an observation (subject), each column is corresponding to a measurement (time) point.
#' @slot t_0  start of the domain (time period), \eqn{t_0}. Default is 0.
#' @slot period length of the domain (time period), \eqn{T}. Default is 1.
#' @slot t_points sequence of the measurement points, \eqn{(t_1,\dots,t_m)}. Default is \eqn{t_k = t_0 + \frac{(2k-1)T}{2(m+1)}}.
#'
#' @export
#' @import methods
#' @author Heyang Ji
#' @examples X = array(rnorm(12),dim = 4:3)
#' functional_variable(X=X,period = 3)
functional_variable = setClass(
  "functional_variable",
  slots = c(X = 'array',
            t_0 = 'numeric',
            period = 'numeric',
            t_points = 'numeric'
  ))
