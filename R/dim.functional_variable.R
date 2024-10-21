#' @title Extract dimensionality of functional data.
#' @description
#' Extract the dimensionality of slot X of functional_variable object.
#'
#' @param x a \code{\link{functional_variable}} object.
#'
#' @return Retruns a 2-element numeric vector.
#' @export
#' @author Heyang Ji
#' @examples fv = functional_variable(X=array(rnorm(12),dim = 4:3),period = 3)
#' dim(fv)
setMethod("dim",
          signature(x="functional_variable"),
          function(x){
            out = dim(slot(x,'X'))
            names(out) = c('subject','time_points')
            return(out)
          })
