# PLatform for Analysis of scEiad

# Installation (bash)

  - You'll need [magick](https://imagemagick.org/index.php) otherwise the package install will fail
  
# Installation (in R)
  - `remotes::install_github('davemcg/plaeApp')`
  -  Get scEiaD (run in bash/Terminal)

    wget http://hpc.nih.gov/~mcgaugheyd/scEiaD/scEiaD__2020_08_12__Mus_musculus_Macaca_fascicularis_Homo_sapiens-5000-counts-TabulaDroplet-batch-scVI-8-0.1-15-7.sqlite.gz
  -  Decompress sqlite file
  
    pigz -d -p 8 scEiaD__2020_08_12__Mus_musculus_Macaca_fascicularis_Homo_sapiens-5000-counts-TabulaDroplet-batch-scVI-8-0.1-15-7.sqlite.gz
  - Edit absolute path to sqlite file in `inst/app/server.R`
  - Run App
