#' @title Extract the value of coefficient parameter function
#' @description
#' Generic function to extract the value of coefficient parameter function of the covariates
#' from linear model with functional covariates at some certain points.
#' @param object An object that represents a functional covariates linear model.
#' @param FC An integer, represent the ordinal number of the functional covariate.
#' Default is 1, which is take the first functional covariate.
#' @param t_points Sequence of the measurement (time) points.
#' @param ... More arguments.
#' @return A numeric atomic vector
#' @export
#' @author Heyang Ji
fc.beta <- function(object,...) UseMethod('fc.beta')
#' @rdname fc.beta
#' @export
setMethod("fc.beta",
          signature(object="fcRegression"),
          function(object, FC = 1, t_points = NULL){
            if(is.null(t_points)){
              t_points = object$data$fc_list[[FC]]@t_points
            }else if( !(
              all(t_points>=object$data$fc_list[[FC]]@t_0) &
              all(t_points<=object$data$fc_list[[FC]]@t_0 + object$data$fc_list[[FC]]@period)
            ) ){
              stop('t_points out of boundary')
            }
            return(basis2fun(object$FC.BasisCoefficient[[FC]],t_points))
          }
)
# fc.beta.fcRegression = function(object, FC = 1, t_points = NULL){
#   if(is.null(t_points)){
#     t_points = object$data$fc_list[[FC]]@t_points
#   }else if( !(
#     all(t_points>=object$data$fc_list[[FC]]@t_0) &
#     all(t_points<=object$data$fc_list[[FC]]@t_0 + object$data$fc_list[[FC]]@period)
#   ) ){
#     stop('t_points out of boundary')
#   }
#   return(basis2fun(object$FC.BasisCoefficient[[FC]],t_points))
# }
#' @rdname fc.beta
#' @export
setMethod("fc.beta",
          signature(object="fcQR"),
          function(object, FC = 1, t_points = NULL){
            if(is.null(t_points)){
              t_points = object$data$fc_list[[FC]]@t_points
            }else if( !(
              all(t_points>=object$data$fc_list[[FC]]@t_0) &
              all(t_points<=object$data$fc_list[[FC]]@t_0 + object$data$fc_list[[FC]]@period)
            ) ){
              stop('t_points out of boundary')
            }
            return(basis2fun(object$FC.BasisCoefficient[[FC]],t_points))
          }
)
# fc.beta.fcQR = function(object, FC = 1, t_points = NULL){
#   if(is.null(t_points)){
#     t_points = object$data$fc_list[[FC]]@t_points
#   }else if( !(
#     all(t_points>=object$data$fc_list[[FC]]@t_0) &
#     all(t_points<=object$data$fc_list[[FC]]@t_0 + object$data$fc_list[[FC]]@period)
#   ) ){
#     stop('t_points out of boundary')
#   }
#   return(basis2fun(object$FC.BasisCoefficient[[FC]],t_points))
# }
#################################################################################
