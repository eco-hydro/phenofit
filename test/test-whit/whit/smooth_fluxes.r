# Apply V-curve to production data Dongdong Kong

library(readxl)
library(spam)
source('V-pack.R')

# Read the data from Excel workbook
fname = 'fluxsites_GPP.xlsx'
shts = excel_sheets(fname)
ns = length(shts)
sht = shts[2]
Data = read_excel(fname, sheet = sht)
x = Data$date
m = length(x)
# m = 400
x = x[1:m]
yraw = Data$GPP_NT[1:m]
w = rep(1, m)
y = yraw[1:m]
nay = is.na(y)
y[nay] = 0
w[nay] = 0

par(mfrow = c(2, 1), mar = c(4, 3, 2, 1), mgp = c(1.5, 0.6, 0))

# Find best lambda with V-curve
d = 2
llas = seq(3, 9, by = 0.2)
vc = v_curve(y, w = w, llas,  d = d, show = T) 

# Plot data and smooth
plot(x, yraw, type = 'l', col = 'darkgrey', xlab = '', ylab = 'GPP_NT')
lines(x, vc$z, col = 'blue', lwd = 1.5)
title(paste(sht, 'd =', d,  'log10(lambda) =', log10(vc$lambda)))

fname = paste(sht, '.pdf', sep = '')
dev.copy2pdf(file = fname)

