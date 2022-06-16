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
  # input[['heatmap_groups']] <- c('CellType (Predict)')
  # input[['heatmap_organism']] <- 'Mus musculus'
  # db <- scEiaD_2020_v01

  gene <- input$heatmap_Gene
  cat(input$gene)
  cat(input$heatmap_groups)
  cat(input$heatmap_organism)
  heatmap_data <- db %>% tbl('diff_testing') %>%
    filter(Gene %in% gene) %>%
    collect() %>%
    filter(Against == 'All') %>%
    filter(Organism == input$heatmap_organism) %>%
    filter(Group == input$heatmap_groups)

  x <- heatmap_data %>%
    dplyr::select(Gene, log2FoldChange, Base) %>%
    tidyr::pivot_wider(values_from=log2FoldChange, names_from = Base)
  celltype <- colnames(x[,-1])
  x <- x %>% data.frame()
  row.names(x) <- x$Gene
  x <- x[,-1]
  Heatmap(t(x), row_labels = celltype,
          name = 'log2 Fold Change',
          use_raster = TRUE,
          row_title = input$heatmap_groups,
          column_title = input$heatmap_organism)
}
