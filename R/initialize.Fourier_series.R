setMethod("initialize", signature(.Object="Fourier_series"),
          function(.Object, double_constant = 0, sin = 0, cos = 0, k_sin = NULL,k_cos = NULL,t_0 = 0, period = 2*pi){
            .Object@double_constant = double_constant
            .Object@cos = cos
            .Object@sin = sin
            if(is.null(k_cos)) k_cos = 1:length(cos)
            if(is.null(k_sin)) k_sin = 1:length(sin)
            .Object@k_cos = as.integer(k_cos)
            .Object@k_sin = as.integer(k_sin)
            .Object@t_0 = t_0
            .Object@period = period
            if (min(.Object@k_cos)<=0 | min(.Object@k_sin)<=0){
              stop("error: k_sin and k_cos must be positive integer")
            }
            if (length(.Object@cos)!=length(.Object@k_cos) |
                length(.Object@sin)!=length(.Object@k_sin)){
              stop("error: length of k_sin and k_cos must match sin and cos")
            }
            return(.Object)
          })
