library(Cairo)
library(scattermore)
library(pals)

library(patchwork)
library(purrr)
library(pool)
library(RSQLite)
library(dplyr)
library(magick)
library(stringr)
library(ComplexHeatmap)

make_heatmap <- function(input, db, meta_filter, cat_to_color_df){
  ### this makes it a little easier to test
  # input <- list()
  # input[['heatmap_Gene']] <- c('CABP5 (ENSG00000105507)','ARR3 (ENSG00000120500)', 'AIF1L (ENSG00000126878)', 'WIF1 (ENSG00000156076)', 'RHO (ENSG00000163914)', 'ONECUT1 (ENSG00000169856)', 'AIF1 (ENSG00000204472)', 'TTR (ENSG00000118271)') #, 'POU4F2')#c('RHO', 'CRX')
  # input[['heatmap_groups']] <- c('CellType_predict')
  # input[['heatmap_organism']] <- c('Mus musculus', "Homo sapiens")
  # input[['heatmap_filter_on']] <- c('Rod','Cone','Horizontal Cell')
  # db <- scEiaD_2020_v01

  swap <- c('Cluster','CellType','CellType (Predict)')
  names(swap) <- c('cluster', 'CellType','CellType_predict')
  gene <- input$heatmap_Gene
  cat(input$gene)
  cat(input$heatmap_groups)
  cat(input$heatmap_organism)
  heatmap_data <- db %>% tbl('diff_testing') %>%
    filter(Gene %in% gene) %>%
    collect() %>%
    filter(Against == 'All') %>%
    filter(Organism  %in% input$heatmap_organism) %>%
    filter(Group == swap[input$heatmap_groups])
  cat(class(input$heatmap_filter_on))
  if (class(input$heatmap_filter_on) == 'character'){
    cat('yo')
    heatmap_data <- heatmap_data %>%
      filter(Base %in% input$heatmap_filter_on)
  }

  x <- heatmap_data %>%
    dplyr::select(Gene, log2FoldChange, Base, Organism) %>%
    mutate(Base = paste(Base, Organism, sep = '__')) %>%
    dplyr::select(-Organism) %>%
    tidyr::pivot_wider(values_from=log2FoldChange, names_from = Base)
  celltype <- colnames(x[,-1]) %>% gsub('__.*','',.)
  organism <- colnames(x[,-1]) %>% gsub('.*__','',.)
  x <- x %>% data.frame()
  row.names(x) <- x$Gene
  x <- x[,-1]
  hm_plot <- Heatmap((x),
                     name = 'log2 Fold Change',
                     column_split = celltype,
                     use_raster = TRUE,
                     column_title_rot = 90,
                     column_labels = organism)
  out <- list()
  #out$plot <- ComplexHeatmap::draw(hm_plot, padding = unit(c(0, 0.5, 1.5, 0.5), "in"), row_title = "Genes")
  out$plot <- hm_plot
  out$data <- as_tibble(x)
  out
  #ComplexHeatmap::draw(hm_plot, padding = unit(c(0, 0.5, 1.5, 0.5), "in"), row_title = "Genes")
}
