rm(list = ls())

folder = '/Users/rachaelshudde/Desktop/'
exact = load(paste(folder, 'exactGP100.RData', sep = ""))
nngp = load(paste(folder, 'NNGP100.Rdata', sep = ""))
wc = load(paste(folder, 'WCGP100.RData', sep = ""))

exact = exactGP100
nngp = NNGP100
wc = WCGP100

# get time means
time_exact = mean(exact[,4])
time_nngp = mean(nngp[,4])
tme_wc = mean(wc[,4])
