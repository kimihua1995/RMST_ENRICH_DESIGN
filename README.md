<h2>Biomerker-Guided Adaptive Enrichment Design with Threshold Detection for Clinical Trials with Time-to-Event Outcome</h2>

I) Code for Numerical Example
  1) `numerical example.R`
      - paramters set-ups
      - plot the KM curves
      - calculate true biomarker cutpoint and RMST difference
      - calculate asymptotic variance in Equation (4.4)
      - calculate critical values of q by Monte Carlo Method
      - plot the global power by required total sample size
  2) `source_numerical example.R`: source codes for the numerical example

II) Code for Smulation Study
  1) `JBS simulation.R`
     - alternative setting
     - global null setting
  2) `source_simulation.R`: source codes for the simulation study
     - `sim.data.piece()`: simulate piecewise exponential time-to-event data
     - `true.cut.piece0()`: true biomarker cutpoint from the piecewise exponential hazard model
     - `rmst.true()`: true marginal RMST within certain biomarker subgroup
     - `dat.modify()`: function for modifying censoring time
     - `my.design()`: function for simulation under alternative setting
     - `my.design.null()`: function for simulation under null setting
     - `stat1()`: RMST estimator
  3) `source_cutpoint.R`: source codes for estimating biomarker-cutpoint
     - functions for RMST regression method: `my.func_surv()`, `my.rmst2reg()`, `fit.rmst.reg()`, `find.cutpoint()`
     - functions for prediction method:
       - alternative setting: `true.cut.piece()`, `fit.piece.0()`, `fit.piece.1()`, `cut.pred1()`, `cut.pred2()`, `cut.pred3()`
       - null setting: `true.cut.piece.null()`, `fit.piece.null.0()`, `fit.piece.null.1()`, `cut.pred1.null()`, `cut.pred3.null()`

