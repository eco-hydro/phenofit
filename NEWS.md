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
