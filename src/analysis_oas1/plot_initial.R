load("results/analysis_oas1/test_initial_nosmooth.RData")
res_nosmooth = res
load("results/analysis_oas1/test_initial_smooth.RData")
res_smooth = res

#plot the clusters
par(mfrow = c(2, 1))
plot(res_nosmooth$phi[2, ], ylim = c(0, 0.015), type = 'l', main = "unsmoothed clusters")
lines(res_nosmooth$phi[1, ], col = 2)
lines(res_nosmooth$phi[3, ], col = 4)
plot(res_smooth$phi[1, ], ylim = c(0, 0.015), type = 'l', main = "smoothed clusters")
lines(res_smooth$phi[2, ], col = 2)
lines(res_smooth$phi[3, ], col = 4)


load("data/oas1/OAS1.Robj")
comp_nosmooth = cbind(res_nosmooth$pi, Robj$g)

#compare the actual genotype with the estimated profile
est_nosmooth = res_nosmooth$pi %*% res_nosmooth$lambda
est_smooth = res_smooth$pi %*% res_smooth$lambda

#look at example profiles belonging to different genotype classes
par(mfrow = c(2, 1))
plot(est_nosmooth[48, ], type = 'l')
lines(est_nosmooth[17, ], col = 2)
lines(est_nosmooth[3, ], col = 3)
plot(est_smooth[48, ], type = 'l')
lines(est_smooth[17, ], col = 2)
lines(est_smooth[3, ], col = 3)

#smooth unsmoothed profiles post-hoc
par(mfrow = c(1, 1))
plot(ashsmooth.pois(est_nosmooth[48, ], log = FALSE), type = 'l')
lines(ashsmooth.pois(est_nosmooth[17, ], log = FALSE), col = 2)
lines(ashsmooth.pois(est_nosmooth[3, ], log = FALSE), col = 3)
