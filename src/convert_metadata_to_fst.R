library(fst)
library(pool)
library(RSQLite)
library(dplyr)
args = commandArgs(trailingOnly = T)
scEiaD_2020_v01 <- dbPool(drv = SQLite(), dbname = args[1], idleTimeout = 3600000)
meta_filter <- scEiaD_2020_v01 %>% tbl('metadata_filter') %>% collect

write_fst(meta_filter, args[2])
# write_fst(meta_filter, 'inst/app/www/metadata_filter.fst')
