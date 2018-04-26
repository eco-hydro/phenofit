# str <- paste("A long string of text goes here just for the purpose \n of illustrating my point Weight ")
# p <- ggplot(mtcars,aes(x=wt,y=mpg))+
#     geom_point()+
#     xlab(str, parse = T)
# p
# try(ggsave(plot=p,filename=<some file>,height=4,width=6))

par(mar = c(4, 4, 2, 0.1))
x_mean <- 1.5
x_sd <- 1.2
hist(rnorm(100, x_mean, x_sd),
     main = substitute(
         paste(X[i], " ~ N(", mu, "=", m, ", ", sigma^2, "=", s2, ")"),
         list(m = x_mean, s2 = x_sd^2)
     )
)
