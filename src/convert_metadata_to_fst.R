library(fst)
library(pool)
library(RSQLite)
library(dplyr)
args = commandArgs(trailingOnly = T)
scEiaD_2020_v01 <- dbPool(drv = SQLite(), dbname = args[1], idleTimeout = 3600000)
meta_filter <- scEiaD_2020_v01 %>% tbl('metadata_filter') %>% collect
dir.create('data/', showWarnings = F)
write_fst(meta_filter, 'data/metadata_filter.fst')
