# phenofit 0.1.3 (2018-07-25)

* fix `Init_param`   
    In previous version, parameters are only rely on "good" (w >= 0.5) values. 
    But when "good" values in a relative low ratio, wield parameters will be 
    generated.   
    In this version, when "good" values ratio is lower than 40%, all values 
    will be used to generate initial parameters.   
* Export weigths of every iteration for growing season dividing function, 
(i.e. `wHANT`, `sgfitw` and `whitsmw2`). And unified their weights updating 
strategy.
* `doubleLog.zhang` is still not as stable as others.
* `wTSM_cpp` iter parameter is ignored now.

# phenofit 0.1.5 (2018-09-19)

* shiny app `check_season` online now.
* `season` can export rough curve fitting result, even no peaks or trough found.

# phenofit 0.1.6 (2018-10-26)

* Melt parameters (i.e. `nptperyear` and `south`) into `INPUT`. `check_input`,
 `season`, `season_3y` and `curvefits` are impacted.
* Add `adj.param` parameter to `season`, which determine whether to automatically 
adjust roughn curve fitting parameters.

# phenofit 0.1.7 (2018-10-29)

* Add V-curve optimization in `season_3y` for Whittaker's parameter lambda.
* Remove `check_fit` and `upper envelope` in `wWHIT`.
* Update `v-curve`.
