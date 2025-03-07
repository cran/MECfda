#' @method predict fcQR
#' @title Predicted values based on fcQR object
#' @description
#' Predicted values based on the Quantile linear model with functional covariates represented by a "fcQR" class object.
#' @details
#' If no new data is input, will return the fitted value.
#' @param object A fcQR class object produced by \code{\link{fcQR}}.
#' @param newData.FC A atomic vector or a matrix or a dataframe or a functional_variable class object or a list of objects above.
#' See argument FC in \code{\link{fcRegression}}.
#' @param newData.Z A dataframe or a matrix or a atomic vector. See argument Z in \code{\link{fcRegression}}.
#' @param ... Further arguments passed to or from other methods \code{\link[quantreg]{predict.rq}}.
#' @return See \code{\link[quantreg]{predict.rq}}.
#' @export
#' @author Heyang Ji
#' @importFrom methods hasArg

# setMethod("predict",
#           signature(object="fcQR"),
#           function(object,newData.FC,newData.Z = NULL, ...) {
#
#           }
# )
#################################################################################
predict.fcQR = function(object,newData.FC,newData.Z = NULL, ...) {
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
            },
            'FPC' = {
              BE = NULL
              for (i in 1:length(newData.fc_list)) {
                X = newData.fc_list[[i]]
                nb = object$FC.BasisCoefficient[[i]]@numeric_basis
                BE_X = numeric_basis_expansion(X,nb)
                colnames(BE_X) = gsub('basisFunction','FPC',colnames(BE_X), fixed = TRUE)
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
