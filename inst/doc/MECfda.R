## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(MECfda)

## ----fd-----------------------------------------------------------------------
fv = functional_variable(
  X = matrix(rnorm(10*24),10,24),
  t_0 = 0,
  period = 1,
  t_points = (0:9)/10
)
dim(fv)

## ----c1-----------------------------------------------------------------------
fsc = Fourier_series(
  double_constant = 3,
  cos = c(0,2/3),
  sin = c(1,7/5),
  k_cos = 1:2,
  k_sin = 1:2,
  t_0 = 0,
  period = 1
)
plot(fsc)
FourierSeries2fun(fsc,seq(0,1,0.05))
extractCoef(fsc)

## ----c2-----------------------------------------------------------------------
bsb = bspline_basis(
  Boundary.knots = c(0,24),
  df             = 7,
  degree         = 3
)
bss = bspline_series(
  coef = c(2,1,3/4,2/3,7/8,5/2,19/10),
  bspline_basis = bsb
)
plot(bss)
bsplineSeries2fun(bss,seq(0,24,0.5))

## ----basis2fun----------------------------------------------------------------
basis2fun(fsc,seq(0,1,0.05))
basis2fun(bss,seq(0,24,0.5))

## ----be-----------------------------------------------------------------------
data(MECfda.data.sim.0.0)
fv = MECfda.data.sim.0.0$FC[[1]]
BE.fs = fourier_basis_expansion(fv,5L)
BE.bs = bspline_basis_expansion(fv,5L,3L)

## ----shili1, eval = FALSE-----------------------------------------------------
# fcRegression(Y, FC, Z, formula.Z, family = gaussian(link = "identity"),
#              basis.type = c("Fourier", "Bspline"), basis.order = 6L,
#              bs_degree = 3)

## ----fcglmm-------------------------------------------------------------------
data(MECfda.data.sim.0.0)
res = fcRegression(FC = MECfda.data.sim.0.0$FC, 
                   Y=MECfda.data.sim.0.0$Y, 
                   Z=MECfda.data.sim.0.0$Z,
                   family = gaussian(link = "identity"),
                   basis.order = 5, basis.type = c('Bspline'),
                   formula.Z = ~ Z_1 + (1|Z_2))
t = (0:100)/100
plot(x = t, y = fc.beta(res,1,t), ylab = expression(beta[1](t)))
plot(x = t, y = fc.beta(res,2,t), ylab = expression(beta[2](t)))
data(MECfda.data.sim.1.0)
predict(object = res, newData.FC = MECfda.data.sim.1.0$FC,
        newData.Z = MECfda.data.sim.1.0$Z)

## ----shili2, eval = FALSE-----------------------------------------------------
# fcQR(Y, FC, Z, formula.Z, tau = 0.5, basis.type = c("Fourier", "Bspline"),
#      basis.order = 6L, bs_degree = 3)

## ----fcqr---------------------------------------------------------------------
data(MECfda.data.sim.0.0)
res = fcQR(FC = MECfda.data.sim.0.0$FC, 
           Y=MECfda.data.sim.0.0$Y, 
           Z=MECfda.data.sim.0.0$Z,
           tau = 0.5,
           basis.order = 5, basis.type = c('Bspline'),
           formula.Z = ~ .)
t = (0:100)/100
plot(x = t, y = fc.beta(res,1,t), ylab = expression(beta[1](t)))
plot(x = t, y = fc.beta(res,2,t), ylab = expression(beta[2](t)))
data(MECfda.data.sim.1.0)
predict(object = res, newData.FC = MECfda.data.sim.1.0$FC,
        newData.Z = MECfda.data.sim.1.0$Z)

## ----shili3, eval = FALSE-----------------------------------------------------
# ME.fcRegression_MEM(
#   data.Y,
#   data.W,
#   data.Z,
#   method = c("UP_MEM", "MP_MEM", "average"),
#   t_interval = c(0, 1),
#   t_points = NULL,
#   d = 3,
#   family.W = c("gaussian", "poisson"),
#   family.Y = "gaussian",
#   formula.Z,
#   basis.type = c("Fourier", "Bspline"),
#   basis.order = NULL,
#   bs_degree = 3,
#   smooth = FALSE,
#   silent = TRUE
# )

## ----MEM, eval = FALSE--------------------------------------------------------
# data(MECfda.data.sim.0.1)
# res = ME.fcRegression_MEM(data.Y = MECfda.data.sim.0.1$Y,
#                           data.W = MECfda.data.sim.0.1$W,
#                           data.Z = MECfda.data.sim.0.1$Z,
#                           method = 'UP_MEM',
#                           family.W = "gaussian",
#                           basis.type = 'Bspline')

## ----shili4, eval = FALSE-----------------------------------------------------
# ME.fcQR_IV.SIMEX(
#   data.Y,
#   data.W,
#   data.Z,
#   data.M,
#   tau = 0.5,
#   t_interval = c(0, 1),
#   t_points = NULL,
#   formula.Z,
#   basis.type = c("Fourier", "Bspline"),
#   basis.order = NULL,
#   bs_degree = 3
# )

## ----iv.simex, eval = FALSE---------------------------------------------------
# rm(list = ls())
# data(MECfda.data.sim.0.2)
# res = ME.fcQR_IV.SIMEX(data.Y = MECfda.data.sim.0.2$Y,
#                        data.W = MECfda.data.sim.0.2$W,
#                        data.Z = MECfda.data.sim.0.2$Z,
#                        data.M = MECfda.data.sim.0.2$M,
#                        tau = 0.5,
#                        basis.type = 'Bspline')

## ----shili5, eval = FALSE-----------------------------------------------------
# ME.fcQR_CLS(
#   data.Y,
#   data.W,
#   data.Z,
#   tau = 0.5,
#   t_interval = c(0, 1),
#   t_points = NULL,
#   grid_k,
#   grid_h,
#   degree = 45,
#   observed_X = NULL
# )

## ----cls, eval = FALSE--------------------------------------------------------
# rm(list = ls())
# data(MECfda.data.sim.0.1)
# res = ME.fcQR_CLS(data.Y = MECfda.data.sim.0.1$Y,
#                   data.W = MECfda.data.sim.0.1$W,
#                   data.Z = MECfda.data.sim.0.1$Z,
#                   tau = 0.5,
#                   grid_k = 4:7,
#                   grid_h = 1:2)

## ----shili6, eval = FALSE-----------------------------------------------------
# ME.fcLR_IV(
#   data.Y,
#   data.W,
#   data.M,
#   t_interval = c(0, 1),
#   t_points = NULL,
#   CI.bootstrap = F
# )

## ----lriv, eval = FALSE-------------------------------------------------------
# rm(list = ls())
# data(MECfda.data.sim.0.3)
# res = ME.fcLR_IV(data.Y = MECfda.data.sim.0.3$Y,
#                  data.W = MECfda.data.sim.0.3$W,
#                  data.M = MECfda.data.sim.0.3$M)

