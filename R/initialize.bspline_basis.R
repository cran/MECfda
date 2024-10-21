setMethod("initialize", signature(.Object="bspline_basis"),
          function(.Object, Boundary.knots = 0:1, knots = NULL, intercept = TRUE, df = NULL, degree = 3L){
            .Object@Boundary.knots = Boundary.knots
            .Object@knots          = knots
            .Object@intercept      = intercept
            .Object@df             = df
            .Object@degree         = as.integer(degree)
            if(!is.null(.Object@df)){
              if (.Object@df - .Object@intercept < .Object@degree){
                stop("degrees of freedom (df) minus intercept (intercept) must be no less than degree of the piecewise polynomial (degree)")
              }
            }
            return(.Object)
          })
