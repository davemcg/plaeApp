make_meta_scatter_umap_plot <- function(input, mf, meta_filter,
                                        celltype_predict_labels,
                                        celltype_labels,
                                        tabulamuris_predict_labels,
                                        cluster_labels,
                                        cat_to_color_df
){
  cat(file=stderr(), paste0(Sys.time(), ' Meta Plot Call\n'))
  # input <- list()
  # input[['meta_column']] <- 'CellType_predict'
  # input[['pt_size_meta']] <- 1
  # input[['gene_and_meta_scatter_tech']] <- 'Droplet'
  # input[['meta_column_transform']] <- 'None'
  # input[['meta_filter_cat']] <- 'CellType_predict'
  # input[['meta_filter_on']] <- 'Bipolar Cells'

  meta_column <- input$meta_column
  transform <- input$meta_column_transform

  pt_size <- input$pt_size_meta %>% as.numeric()
  filter_column <- input$meta_column
  if (transform == 'log2' && is.numeric(meta_filter[,meta_column] %>% pull(1))){
    cat('log2 time')
    meta_filter[,meta_column] <- log2(meta_filter[,meta_column] + 1)
  }
  p_data <- meta_filter %>%
    filter(!grepl('Doub|\\/Margin\\/Periocular', CellType)) %>%
    filter(!is.na(!!as.symbol(meta_column))) %>%
    sample_frac(0.3)
  # category filtering
  if (!is_null(input$meta_filter_cat)){
    validate( need(input$meta_filter_on != '', 'Please select at least one value in "Meta Filter on" '  ))
    if (class(input$meta_filter_on) == 'character'){
      p_data <- p_data %>%
        #filter(!!as.symbol(input$meta_filter_cat) %in% input$meta_filter_on)
        filter_at(vars(all_of(input$meta_filter_cat)), all_vars(. %in% input$meta_filter_on))
    } else {
      p_data <- p_data %>%
        filter(!!as.symbol(input$meta_filter_cat) >= input$meta_filter_on[1],
               !!as.symbol(input$meta_filter_cat) <= input$meta_filter_on[2])
    }
  }


  # metadata NUMERIC plot --------------
  if (is.numeric(meta_filter[,meta_column] %>% pull(1)) ){
    color_range <- range(p_data[,meta_column] %>% pull(1))
    suppressWarnings(plot <- ggplot() +
                       geom_scattermost(cbind(mf %>%
                                                filter(is.na(!!as.symbol(meta_column))) %>% pull(UMAP_1),
                                              mf %>%
                                                filter(is.na(!!as.symbol(meta_column))) %>% pull(UMAP_2)),
                                        pointsize = pt_size - 1, color = '#D3D3D333',
                                        pixels = c(1000,1000)) +
                       geom_scattermost(cbind(p_data$UMAP_1, p_data$UMAP_2),
                                        color = viridis::viridis(100, alpha=0.3)
                                        [1+99*((p_data[,meta_column] %>% pull(1))-color_range[1])/diff(color_range)],
                                        pointsize= pt_size - 1,
                                        pixels=c(1000,1000),
                                        interpolate=FALSE) +
                       geom_point(data=data.frame(x=double(0)), aes(x,x,color=x))  +
                       scale_color_gradientn(  #add the manual guide for the empty aes
                         limits=c(min(p_data[,meta_column] %>% pull(1)),
                                  max(p_data[,meta_column] %>% pull(1))),
                         colors=viridis::viridis(100),
                         name=meta_column) +
                       theme_cowplot() +
                       guides(colour = guide_legend(override.aes = list(size=7, alpha = 1))) +
                       theme(axis.line = element_blank(),
                             axis.title = element_blank(),
                             axis.ticks = element_blank(),
                             axis.text = element_blank()) +
                       annotate("text", -Inf, Inf, label = "Metadata", hjust = 0, vjust = 1, size = 6))
    # metadata CATEGORICAL plot --------------
  } else {
    cur_color_df <- cat_to_color_df %>%
      filter( meta_category %in%  meta_column,
              value %in% p_data[[meta_column]]) %>% distinct

    color_list <- cur_color_df %>% pull(color)
    # replaced join with this for speed
    k <- cur_color_df$color
    names(k) <- cur_color_df$value
    np_color <- {k[p_data[[meta_column]] ]} %>% paste0(., '33')
    names(np_color) <- NULL
    color_data <- cur_color_df  %>%
      select(value) %>%
      mutate( x=0)
    names(color_list) <- color_data$value

    suppressWarnings(plot <- ggplot() +
                       geom_scattermost(cbind(mf %>%
                                                filter(is.na(!!as.symbol(meta_column))) %>% pull(UMAP_1),
                                              mf %>%
                                                filter(is.na(!!as.symbol(meta_column))) %>% pull(UMAP_2)),
                                        pointsize = pt_size - 1, color = '#D3D3D333',
                                        pixels = c(1000,1000)) +
                       geom_scattermost(cbind(p_data$UMAP_1, p_data$UMAP_2),
                                        color = np_color ,
                                        pointsize= pt_size - 1,
                                        pixels=c(1000,1000),
                                        interpolate=FALSE) +
                       #geom_point(data=data.frame(x=double(0)), aes(x,x,color=x))  +
                       geom_point(data=color_data, aes(x,x,color=value)) +
                       scale_colour_manual(name= meta_column,
                                           values = color_list) +
                       theme_cowplot() +
                       theme(axis.line = element_blank(),
                             axis.title = element_blank(),
                             axis.ticks = element_blank(),
                             axis.text = element_blank()) +
                       guides(colour = guide_legend(override.aes = list(alpha = 1, size = 7)))
    )

  }


  more <- NULL
  if ('1' %in% input$label_toggle){
    more <- geom_text_repel(data = celltype_labels, bg.color = 'white',
                            alpha = 0.7,
                            aes(x = UMAP_1, y = UMAP_2, label = CellType))
  }
  if ('2' %in% input$label_toggle){
    more <- geom_text_repel(data = celltype_predict_labels, bg.color = 'white',
                            alpha = 0.7,
                            aes(x = UMAP_1, y = UMAP_2, label = CellType_predict))
  }
  if ('3' %in% input$label_toggle){
    more <- geom_text_repel(data = cluster_labels, bg.color = 'white',
                            alpha = 0.7,
                            aes(x = UMAP_1, y = UMAP_2, label = cluster),
                            max.iter = 20)
  }
  if ('4' %in% input$label_toggle){
    more <- geom_text_repel(data = tabulamuris_predict_labels, bg.color = 'white',
                            alpha = 0.7,
                            aes(x = UMAP_1, y = UMAP_2, label = TabulaMurisCellType_predict),
                            max.iter = 20)
  }
  col_size <- {meta_filter[[meta_column]]} %>% n_distinct
  out <- list()
  out$data <- p_data %>% select(Barcode, UMAP_1, UMAP_2, all_of(meta_column ))
  out$plot <- plot + more
  out$col_size <- col_size
  out
}
