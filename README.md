# PLatform for Analysis of scEiad

# Installation (in R)
  - `devtools::install_github('davemcg/plaeApp')`
  -  Get scEiaD (run in bash/Terminal)
  
    wget http://hpc.nih.gov/~mcgaugheyd/scEiaD/MOARTABLES__anthology_limmaFALSE___Mus_musculus_Macaca_fascicularis_Homo_sapiens-5000-counts-TabulaDroplet-batch-scVI-8-0.1-15-7.sqlite.gz
  -  Decompress sqlite file
  
    pigz -d -p 8 MOARTABLES__anthology_limmaFALSE___Mus_musculus_Macaca_fascicularis_Homo_sapiens-5000-counts-TabulaDroplet-batch-scVI-8-0.1-15-7.sqlite.gz
  - Edit absolute path to sqlite file in `inst/app/server.R`
  - Run App
