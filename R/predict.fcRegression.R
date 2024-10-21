#' @method predict fcRegression
#' @title Predicted values based on fcRegression object
#' @description
#' Predicted values based on the linear model with functional covariates represented by a "fcRegression" class object.
#' @details
#' If no new data is input, will return the fitted value.
#' @param object A fcRegression class object produced by \code{\link{fcRegression}}.
#' @param newData.FC A atomic vector or a matrix or a dataframe or
#' a functional_variable class object or a list of objects above.
#' See argument FC in \code{\link{fcRegression}}.
#' @param newData.Z A dataframe or a matrix or a atomic vector.
#' See arguement Z in \code{\link{fcRegression}}.
#' @param ... Further arguments passed to or from other methods,
#' including \code{\link[stats]{predict.lm}}, \code{\link[stats]{predict.glm}},
#' \code{\link[lme4]{predict.merMod}}.
#' @return See \code{\link[stats]{predict.lm}}, \code{\link[stats]{predict.glm}},
#' \code{\link[lme4]{predict.merMod}}.
#' @export
#' @author Heyang Ji
#' @examples
#' data(MECfda.data.sim.0.0)
#' res = fcRegression(FC = MECfda.data.sim.0.0$FC, Y=MECfda.data.sim.0.0$Y, Z=MECfda.data.sim.0.0$Z,
#'                    basis.order = 5, basis.type = c('Bspline'),
#'                    formula.Z = ~ Z_1 + (1|Z_2))
#' data(MECfda.data.sim.1.0)
#' predict(object = res, newData.FC = MECfda.data.sim.1.0$FC,newData.Z = MECfda.data.sim.1.0$Z)
#' @importFrom methods hasArg


# setMethod("predict",
#           signature(object="fcRegression"),
#           function(object,newData.FC,newData.Z = NULL, ...) {
#
#           }
# )
###################################################################################
predict.fcRegression = function(object,newData.FC,newData.Z = NULL, ...) {
  if(!methods::hasArg(newData.FC)){
    if(!is.null(newData.Z)){
      stop('newData.FC must be input')
    }else{
      return(predict(object$regression_result, ...))
    }
  }

  if(is.atomic(newData.FC) | is.matrix(newData.FC) | is.data.frame(newData.FC) | any(class(newData.FC) == "functional_variable")){

    newData.fc_list = list(newData.FC)
  }else if(is.list(newData.FC)){

    newData.fc_list = newData.FC
  }else{
    stop('Incorrect format of newData.FC')
  }
  rm(newData.FC)

  if (length(newData.fc_list) != length(object$FC.BasisCoefficient)) stop("dimensionality of variables doesn't match")

  for (k in 1:length(newData.fc_list)) {
    fc = newData.fc_list[[k]]
    if (!any(class(fc) == "functional_variable")){
      if(is.atomic(fc)&(!is.matrix(fc))){
        fc = t(fc)
        fc = functional_variable(X=fc,
                                 t_0      = object$data$fc_list[[k]]@t_0,
                                 period   = object$data$fc_list[[k]]@period,
                                 t_points = object$data$fc_list[[k]]@t_points)
      }else if( is.matrix(fc) | is.data.frame(fc) ){
        fc = as.matrix(fc)
        fc = functional_variable(X=fc,
                                 t_0      = object$data$fc_list[[k]]@t_0,
                                 period   = object$data$fc_list[[k]]@period,
                                 t_points = object$data$fc_list[[k]]@t_points)
      }else{
        stop('Incorrect format of newData.FC')
      }
    }
    newData.fc_list[[k]] = fc
  }


  if(!is.null(newData.Z)){
    newData.Z = as.data.frame(newData.Z)

  }else{
    fc = newData.fc_list[[k]]
    newData.Z = as.data.frame(matrix(nrow = dim(fc)['subject'],ncol = 0))
  }





  for (fc in newData.fc_list) {
    if (dim(fc)['subject'] != nrow(newData.Z)) stop("dimensionality of variables doesn't match")
  }
  rm(fc)
  names(newData.fc_list) = names(object$data$fc_list)


  {
    basis.order = object$basis.order
    switch (object$function.basis.type,
            'Fourier' = {
              BE = NULL
              for (i in 1:length(newData.fc_list)) {
                X = newData.fc_list[[i]]
                n_k = basis.order[i]
                BE_X = fourier_basis_expansion(X,n_k)
                colnames(BE_X) = paste(names(newData.fc_list)[i],colnames(BE_X),sep = '.')
                BE = cbind(BE,BE_X)
              }
            },
            'Bspline' = {
              bs_degree = object$bs_degree
              BE = NULL
              for (i in 1:length(newData.fc_list)) {
                X = newData.fc_list[[i]]
                n_k = basis.order[i]
                BE_X = bspline_basis_expansion(X,n_k,bs_degree)
                colnames(BE_X) = paste(names(newData.fc_list)[i],colnames(BE_X),sep = '.')
                BE = cbind(BE,BE_X)
              }
            }
    )
    rm(X,i)
    data.funRegress = as.data.frame(cbind(BE,newData.Z))
    rm(BE_X)
  }
  predict(object = object$regression_result, newdata = data.funRegress, ...)
}



