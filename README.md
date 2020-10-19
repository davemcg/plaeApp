# PLatform for Analysis of scEiad

# Installation (bash)

  - You'll need [magick](https://imagemagick.org/index.php) otherwise the package install will fail
  
# Installation (in R)
  - `Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")` #otherwise some warnings will, by default, be converted to errors (https://github.com/r-lib/remotes/issues/403)
  - `remotes::install_github('davemcg/plaeApp')`
  -  Get scEiaD (run in bash/Terminal)

    wget http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/scEiaD-0-2000-counts-TabulaDroplet-batch-scVI-8-0.001-500-0.6.sqlite.gz
  -  Decompress sqlite file
  
    pigz -d -p 8 scEiaD-0-2000-counts-TabulaDroplet-batch-scVI-8-0.001-500-0.6.sqlite.gz
  - Edit absolute path to sqlite file in `inst/app/server.R`
  - Run App
