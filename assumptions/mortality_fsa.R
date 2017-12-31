

## FTs levetidsmodel

benchmarkFT <- read.table("data/Benchmark_doedelighed2012.csv", header = T, sep = ';', dec = '.')
benchmarkR  <- read.table("data/Benchmark_levetidsforbedringer2012.csv", header = T, sep = ';', dec = '.')

# head(benchmarkFT)
# head(benchmarkR)

femaleR  <- benchmarkR$Kvinder
femaleMu <- benchmarkFT$Kvinder
maleR    <- benchmarkR$Maend
maleMu   <- benchmarkFT$Maend

# # FTs levetidsmodel
# muFT <- function(x,t,k) {
#   # x: alder
#   # t: ?r ud i fremtiden
#   # k: k?n
#   if (k == 'K') (1 - femaleR[x + 1 + t]) ^ t * femaleMu[x + 1 + t]
#   else if (k == 'M') (1 - maleR[x + 1 + t]) ^ t * maleMu[x + 1 + t]
# }


# FTs levetidsmodel
muFT <- function(x, t, k) {
  # x: alder
  # t: ?r ud i fremtiden
  # k: k?n
  
  lower <- floor(x + t)
  upper <- ceiling(x + t)
  
  if (k == 0) {
    femaleMu_temp <- ((x + t) - lower) * femaleMu[upper + 1] + (upper - (x + t)) * femaleMu[lower + 1]
    femaleR_temp <- ((x + t) - lower) * femaleR[upper + 1] + (upper - (x + t)) * femaleR[lower + 1]
    res <- (1 - femaleR_temp) ^ t * femaleMu_temp
  }
  else if (k == 1) {
    maleMu_temp <- ((x + t) - lower) * maleMu[upper + 1] + (upper - (x + t)) * maleMu[lower + 1]
    maleR_temp <- ((x + t) - lower) * maleR[upper + 1] + (upper - (x + t)) * maleR[lower + 1]
    res <- (1 - maleR_temp) ^ t * maleMu_temp
    
  }
  
  res <- ifelse(res == 0,
                (1 - maleR[x + 1 + t]) ^ t * maleMu[x + 1 + t], ## ONLY working for men at the moment!!!
                res)
  
  
  return(res)
  
}
