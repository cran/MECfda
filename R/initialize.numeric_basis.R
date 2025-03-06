setMethod("initialize", signature(.Object="numeric_basis"),
          function(.Object, basis_function, t_points, t_0 = 0, period = 1){
            # if(nrow(basis_function)<ncol(basis_function)){
            #   stop("The matrix of basis function values should be column full rank.")
            # }
            # if(Matrix::rankMatrix(basis_function)<ncol(basis_function)){
            #   stop("The matrix of basis function values should be column full rank.")
            # }
            if(any(t_0>t_points)){
              stop("The value of t_points should be within the time domain.")
            }
            if(any(t_0+period<t_points)){
              stop("The value of t_points should be within the time domain.")
            }
            .Object@basis_function = basis_function
            .Object@t_points       = t_points
            .Object@t_0            = t_0
            .Object@period         = period
            return(.Object)
          })
