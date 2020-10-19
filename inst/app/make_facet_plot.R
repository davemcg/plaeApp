make_facet_plot <- function(input, meta_filter){
  cat(file=stderr(), paste0(Sys.time(), ' Facet Plot Call\n'))
  # # testing
  # input <- list()
  # input$pt_size_facet = 1
  # input$facet <- 'organism'
  # input$facet_color <- 'CellType'
  facet_column <- input$facet
  color_column <- input$facet_color
  pt_size <- input$pt_size_facet %>% as.numeric()


  if (!is.null(input$facet_filter_cat)){
    validate(
      need(input$facet_filter_on != '', "Please select at least one feature to filter on")
    )}

  if (!is.null(input$facet_filter_cat)){
    gray_data <- meta_filter %>%
      filter(is.na(!!as.symbol(color_column))) %>%
      filter_at(vars(all_of(input$facet_filter_cat)), all_vars(. %in% input$facet_filter_on))
    p_data <- meta_filter %>%
      filter(!is.na(!!as.symbol(facet_column)),
             !is.na(!!as.symbol(color_column))) %>%
      filter_at(vars(all_of(input$facet_filter_cat)), all_vars(. %in% input$facet_filter_on))
  } else {
    gray_data <- meta_filter %>%
      filter(is.na(!!as.symbol(color_column)))
    p_data <- meta_filter %>%
      filter(!is.na(!!as.symbol(color_column)))
  }
  # if(input$facet_filter != '' | any(input$facet_filter != '') ){
  #   p_data <- p_data %>% filter(!!as.symbol(facet_column) %in% input$facet_filter)
  # }


  num_col = min(5, p_data %>% pull(!!as.symbol(facet_column)) %>% n_distinct )#fix the number of facet columns at 5, because it get too squished past that
  suppressWarnings(plot <- ggplot(data = p_data) +
                     geom_scattermore(data = gray_data,
                                      aes(x = UMAP_1, y = UMAP_2),
                                      color = 'gray',
                                      pointsize = pt_size,
                                      pixels = c(750,750),
                                      alpha = 0.4) +
                     geom_scattermore(aes(x = UMAP_1, y = UMAP_2,
                                          color = !!as.symbol(color_column)) ,
                                      pointsize= pt_size,
                                      pixels = c(750,750),
                                      alpha = 0.6) +
                     facet_wrap(vars(!!(as.symbol(facet_column))), ncol = num_col) +
                     scale_colour_manual(values = rep(c(pals::alphabet() %>% unname(),
                                                        pals::alphabet2() %>% unname()),
                                                      times = 20),
                                         na.value = 'gray') +
                     guides(colour = guide_legend(override.aes = list(alpha = 1, size = 7))) +
                     theme_cowplot() +
                     theme(axis.line = element_blank(),
                           axis.title = element_blank(),
                           axis.ticks = element_blank(),
                           axis.text = element_blank())
  )
  plot

}
