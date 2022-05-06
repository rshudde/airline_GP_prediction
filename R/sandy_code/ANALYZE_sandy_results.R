rm(list = ls())

# read in 100
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
time_wc = mean(wc[,4])

# read in 200

# read in 300

# read in 400

# read in 500

# make plots