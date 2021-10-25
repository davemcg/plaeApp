make_gene_scatter_umap_plot <- function(input,
                                        db,
                                        mf,
                                        meta_filter,
                                        celltype_predict_labels,
                                        celltype_labels,
                                        tabulamuris_predict_labels,
                                        cluster_labels){
  cat(file=stderr(), paste0(Sys.time(), ' Gene Scatter Plot Call\n'))
  # input <- list()
  # input$Gene <- 'CRX (ENSG00000105392)'
  # input$pt_size <- 1
  # input$gene_scater_slider <- c(1,15)
  # db <- scEiaD_2020_v01
  gene <- input$Gene


  pt_size <- input$pt_size_gene %>% as.numeric()
  expression_range <- input$gene_scatter_slider

  p <-  db %>% tbl('counts') %>%
    filter(Gene == gene) %>%
    collect() %>%
    #mutate(counts = counts - min(counts) + 1) %>%
    left_join(., meta_filter, by = 'Barcode') %>%
    filter(!is.na(UMAP_1), !is.na(UMAP_2), !is.na(counts))
  cat(input$gene_filter_cat)
  cat(class(input$gene_filter_cat))
  if (!is_null(input$gene_filter_cat)){
    validate( need(input$gene_filter_on != '', 'Please select at least one value in "Gene Filter on" '  ))
    if (class(input$gene_filter_on) == 'character'){
      p <- p %>%
        #filter(!!as.symbol(input$gene_filter_cat) %in% input$gene_filter_on)
        filter_at(vars(all_of(input$gene_filter_cat)), all_vars(. %in% input$gene_filter_on))
    } else {
      p <- p %>%
        filter(!!as.symbol(input$gene_filter_cat) >= input$gene_filter_on[1],
               !!as.symbol(input$gene_filter_cat) <= input$gene_filter_on[2])
    }
  }


  # remove effect of super long with capping max val at 5 * SD
  max_value <- sd(p$counts) * 5
  p <- p %>% select(Barcode, UMAP_1, UMAP_2, counts) %>% filter(!is.na(counts)) %>% mutate(counts = case_when(counts > max_value ~ max_value, TRUE ~ counts)) %>%
    filter(counts >= as.numeric(expression_range[1]),
           counts < as.numeric(expression_range[2]))
  color_range <- range(p$counts)
  #color_range <- c(1,5)
  plot <- p %>% ggplot() +
    geom_scattermost(cbind(mf$UMAP_1, mf$UMAP_2), color = '#707070',
                     pointsize = 0,
                     pixels=c(1000,1000)) +
    geom_scattermost(cbind(p$UMAP_1, p$UMAP_2),
                     color = viridis::magma(100, alpha=0.2)
                     [1+99*(p$counts-color_range[1])/diff(color_range)],
                     pointsize= pt_size ,
                     pixels=c(1000,1000),
                     interpolate=FALSE) +
    geom_point(data=data.frame(x=double(0)), aes(x,x,color=x)) +
    scale_color_gradientn(  #add the manual guide for the empty aes
      limits=c(min(p$counts),max(p$counts)),
      colors=viridis::magma(100),
      name="log2(counts+1)") +
    theme_cowplot() +
    theme(axis.line = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank()) +
    annotate("text", -Inf, Inf, label = paste0(gene, ' expression'), hjust = 0, vjust = 1, size = 6)

  more <- NULL
  if ('1' %in% input$gene_label_toggle){
    more <- geom_text_repel(data = celltype_labels, bg.color = 'white',
                            aes(x = UMAP_1, y = UMAP_2, label = CellType))
  }
  if ('2' %in% input$gene_label_toggle){
    more <- geom_text_repel(data = celltype_predict_labels, bg.color = 'white',
                            aes(x = UMAP_1, y = UMAP_2, label = CellType_predict))
  }
  if ('3' %in% input$gene_label_toggle){
    more <- geom_text_repel(data = cluster_labels, bg.color = 'white',
                            aes(x = UMAP_1, y = UMAP_2, label = cluster),
                            max.iter = 20)
  }
  if ('4' %in% input$gene_label_toggle){
    more <- geom_text_repel(data = tabulamuris_predict_labels, bg.color = 'white',
                            aes(x = UMAP_1, y = UMAP_2, label = TabulaMurisCellType_predict),
                            max.iter = 20)
  }


  suppressWarnings(plot + more)
}
