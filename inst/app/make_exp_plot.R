make_exp_plot <- function(input, db, meta_filter){
  cat(file=stderr(), paste0(Sys.time(), ' Exp Plot Call\n'))
  gene <- input$exp_plot_genes

  grouping_features <- input$exp_plot_groups

  if (input$exp_filter_cat != ''){
    # box_data <- db %>% tbl('grouped_stats') %>%
    #   filter(Gene %in% gene) %>%
    #   collect() %>%
    #   filter(!!as.symbol(input$exp_filter_cat) %in% input$exp_filter_on)
    #

    if (class(input$exp_filter_on) == 'character'){
      box_data <- db %>% tbl('grouped_stats') %>%
        filter(Gene %in% gene) %>%
        collect() %>%
        filter_at(vars(all_of(input$exp_filter_cat)), all_vars(. %in% input$exp_filter_on))
      #filter(!!as.symbol(input$exp_filter_cat) %in% input$exp_filter_on)
    } else {
      box_data <- db %>% tbl('grouped_stats') %>%
        filter(Gene %in% gene) %>%
        collect() %>%
        #filter(!!as.symbol(input$exp_filter_cat) %in% input$exp_filter_on) %>%
        filter(!!as.symbol(input$exp_filter_cat) >= input$exp_filter_on[1],
               !!as.symbol(input$exp_filter_cat) <= input$exp_filter_on[2])
    }

  } else {
    box_data <- db %>% tbl('grouped_stats') %>%
      filter(Gene %in% gene) %>%
      collect()
  }
  validate(
    need(input$exp_plot_groups != '', "Please select at least one grouping feature")
  )

  #cat(input)
  box_data <- box_data %>%
    #filter(!is.na(!!as.symbol(grouping_features))) %>%
    group_by_at(vars(one_of(c('Gene', input$exp_plot_facet, grouping_features)))) %>%
    summarise(cpm = sum(cpm * cell_exp_ct) / sum(cell_exp_ct),
              cell_exp_ct = sum(cell_exp_ct, na.rm = TRUE)) %>%
    full_join(., meta_filter %>%
                group_by_at(vars(one_of(grouping_features))) %>%
                summarise(Count = n())) %>%
    mutate(cell_exp_ct = ifelse(is.na(cell_exp_ct), 0, cell_exp_ct)) %>%
    mutate(`%` = round((cell_exp_ct / Count) * 100, 2),
           Expression = round(cpm * (`%` / 100), 2)) %>%
    select_at(vars(one_of(c('Gene', grouping_features, 'cell_exp_ct', 'Count', '%', 'Expression')))) %>%
    arrange(-Expression) %>%
    rename(`Cells # Detected` = cell_exp_ct,
           `Total Cells` = Count,
           `Mean CPM` = Expression,
           `% of Cells Detected` = `%`) %>%
    tidyr::drop_na()
  box_data$Group <- box_data[,c(2:(length(grouping_features)+1))] %>% tidyr::unite(x, sep = ' ') %>% pull(1)

  box_data %>%
    ggplot(aes(x=Gene, y = !!as.symbol(input$exp_plot_ylab), color = !!as.symbol(grouping_features))) +
    geom_boxplot(color = 'black', outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(aes(size = `Total Cells`), grouponX = TRUE) +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_colour_manual(values = rep(c(pals::alphabet() %>% unname()), 20)) +
    theme(legend.position="bottom") +
    facet_wrap(ncol = as.numeric(input$exp_plot_col_num), scales = 'free_x', vars(!!as.symbol(input$exp_plot_facet)))
}
