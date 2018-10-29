## second solution
library(grid)
library(gridExtra)

label_value <- function (labels, multi_line = TRUE, sep = "*'~'*")
{
    out <- do.call("Map", c(list(paste, sep = sep), labels))
    out <- list(unname(unlist(out)))

    lapply(out, function(values) {
        values <- paste0("list(", values, ")")
        lapply(values, function(expr) c(parse(text = expr)))
    })
}

## when compared with other methods, iter2 was used.

itersI <- "iter2"
## Update 20180904
# 1. Roughness error,
# Rg is the Roughness of normalized Y_pred
# Rg_0 is ... of Y_pred

subfix <- "Rg_norm_by_pred"
# Rg_norm_by_obs, Rg_norm_by_pred, Rg
d <- df[meth %in% methods2, .(site, meth, type, iters, RMSE, R2 = R2, Bias, Rg = Rg_norm_by_pred)] %>%
    melt(id.vars = c("site", "meth", "type", "iters"), variable.name = "index") #, "perc"
d <- merge(st[, .(site, IGBPname)], d)

d$index %<>% factor(indice, indice_label)# <- indice# factor(indices)
# d$index %<>% mapvalues(indice, )

d2 <- d[iters == itersI] %>% dcast(., site+IGBPname+iters+index~meth, value.var = "value") %>%
    melt(c("site", "IGBPname", "iters", "index", "wWH"), variable.name = "meth")
d2[, kind:=0]

d_diff <- data.table(index = levels(d2$index) %>% factor(., .),
                     dmax = c(.05, .02, 0.02, .005))
d2 %<>% merge(d_diff)
d2[, label:= sprintf("%s~%s", meth, index)]

d2[wWH - value >  dmax, kind := 1]
d2[wWH - value < -dmax, kind := -1]
d2$kind %<>% factor(levels = c(1, 0, -1),
                    labels = c("Bigger", "Similar", "Smaller"))

colors <- scales::hue_pal()(3)
colors <- c(colors[2], "grey60", colors[1])

