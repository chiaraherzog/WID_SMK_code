source("0-source/installLibs.R")

cat('Installing libraries if not available...\n')

load("0-source/libs_CRAN.Rdata")
installLibs(libraries, type = 'CRAN')

load("0-source/libs_bioC.Rdata")
installLibs(libraries_bioC, type = 'bioc')

load("0-source/libs_gith.Rdata")
installLibs(libraries_gith, type = 'github')


cat('Done.\n')