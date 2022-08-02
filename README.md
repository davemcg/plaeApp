# PLatform for Analysis of scEiad

## What is this?
The code-base for plae.nei.nih.gov

By cloning this repo and downloading (see below) the scEiaD database, you can run PLAE on your local computer.

## Warning
scEiaD is a sqlite file that is approximately **306 GB when uncompressed**. Not a typo! 

# Installation (bash)

  - You'll need [magick](https://imagemagick.org/index.php) otherwise the package install will fail
  
# Installation (in R)
  - `Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")` #otherwise some warnings will, by default, be converted to errors (https://github.com/r-lib/remotes/issues/403)
  - `remotes::install_github('davemcg/plaeApp')`
  -  Get scEiaD (run in bash/Terminal) and the metadata in [FST](https://www.fstpackage.org)

    wget http://hpc.nih.gov/~mcgaugheyd/scEiaD/2022_03_22/MOARTABLES__anthology_limmaFALSE___4000-counts-universe-study_accession-scANVIprojection-15-5-0.1-50-20__pointRelease01.sqlite.gz
    wget http://hpc.nih.gov/~mcgaugheyd/scEiaD/2022_03_22/meta_filter.fst
    
  -  Decompress sqlite file
  
    pigz -d -p 8 MOARTABLES__anthology_limmaFALSE___4000-counts-universe-study_accession-scANVIprojection-15-5-0.1-50-20__pointRelease01.sqlite.gz
  - Edit absolute path to sqlite file in `inst/app/server.R`

# Note
  - If you've cloned this repo, then you can create the `meta_filter.fst` file yourself by running `Rscript ~/path/to/plaeApp/src/convert_metadata_to_fst.R /path/to/scEiaD.sqlite /output/path/meta_filter.fst`
