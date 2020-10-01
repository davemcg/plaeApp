make_gene_scatter_umap_plot <- function(input, db, mf, meta_filter){
  cat(file=stderr(), paste0(Sys.time(), ' Gene Scatter Plot Call\n'))
  gene <- input$Gene
  tech <- input$gene_and_meta_scatter_tech
  pt_size <- input$pt_size_gene %>% as.numeric()
  expression_range <- input$gene_scatter_slider
  mf <- mf %>% filter(TechType == tech)
  p <-  db %>% tbl('cpm') %>%
    filter(Gene == gene) %>%
    collect() %>%
    mutate(cpm = cpm - min(cpm) + 1) %>%
    filter(cpm > as.numeric(expression_range[1]),
           cpm < as.numeric(expression_range[2])) %>%
    left_join(., meta_filter, by = 'Barcode') %>%
    filter(TechType == tech, !is.na(UMAP_1), !is.na(UMAP_2), !is.na(cpm))
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


  color_range <- range(p$cpm)
  plot <- p %>% ggplot() +
    geom_scattermost(cbind(mf$UMAP_1, mf$UMAP_2), color = '#D3D3D333',
                     pointsize = pt_size * 1.5,
                     pixels=c(750,750)) +
    geom_scattermost(cbind(p$UMAP_1, p$UMAP_2),
                     color = viridis::magma(100, alpha=0.3)
                     [1+99*(p$cpm-color_range[1])/diff(color_range)],
                     pointsize= pt_size,
                     pixels=c(750,750),
                     interpolate=FALSE) +
    geom_point(data=data.frame(x=double(0)), aes(x,x,color=x)) +
    scale_color_gradientn(  #add the manual guide for the empty aes
      limits=c(min(p$cpm),max(p$cpm)),
      colors=viridis::magma(100),
      name="log2(cpm+1)") +
    theme_cowplot() +
    theme(axis.line = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank()) +
    annotate("text", -Inf, Inf, label = paste0(gene, ' expression'), hjust = 0, vjust = 1, size = 6)

  suppressWarnings(plot)
}
