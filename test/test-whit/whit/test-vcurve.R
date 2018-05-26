## test V-curve
system.time(whitV(x, IsPlot = TRUE))
system.time(v_curve(INPUT$y, INPUT$w, llas = seq(1, 6, 0.1),  d = 2, show = T))
