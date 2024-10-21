setMethod("initialize", signature(.Object="functional_variable"),
          function(.Object, X, t_0 = 0, period = 1, t_points = NULL){
            .Object@X = X
            .Object@t_0 = t_0
            .Object@period = period
            if (is.null(t_points)){
              n_t = ncol(X)
              .Object@t_points = .Object@t_0 + .Object@period*(1:n_t - 0.5)/n_t
            }else{
              .Object@t_points = t_points
            }
            return(.Object)
          })
