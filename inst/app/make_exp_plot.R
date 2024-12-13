make_exp_plot <- function(input, db, meta_filter){
  # TESTING VALUES
  # db <- scEiaD_2020_v01
  # input <- list()
  # input$exp_plot_facet <- 'CellType_predict'
  # input$exp_filter_cat <- 'CellType_predict'
  # input$exp_filter_on <- c('Rod' ,'Cone')
  # input$exp_plot_facet <- c('cluster','Stage')
  # input$exp_plot_genes <- c('PAX6 (ENSG00000007372)', 'CRX (ENSG00000105392)', 'NRL (ENSG00000129535)', 'RHO (ENSG00000163914)')
  # input$exp_plot_groups <- c('study_accession')
  # input$exp_plot_ylab <- 'Mean log2(Counts + 1)'
  # input$exp_plot_col_num <- 10
  # input$exp_filter_min_cell_number <- 50

  cat(file=stderr(), paste0(Sys.time(), ' Exp Plot Call\n'))
  gene <- input$exp_plot_genes

  grouping_features <- input$exp_plot_groups

  if (!is.null(input$exp_filter_cat)){
    validate(
      need(input$exp_filter_on != '', "Please select at least one feature to filter on")
    )}

  meta_filter_EXP <- meta_filter
  if (!is.null(input$exp_filter_cat)){
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
      meta_filter_EXP <- meta_filter %>%
        filter_at(vars(all_of(input$exp_filter_cat)), all_vars(. %in% input$exp_filter_on))
      #filter(!!as.symbol(input$exp_filter_cat) %in% input$exp_filter_on)
    } else {
      box_data <- db %>% tbl('grouped_stats') %>%
        filter(Gene %in% gene) %>%
        collect() %>%
        #filter(!!as.symbol(input$exp_filter_cat) %in% input$exp_filter_on) %>%
        filter(!!as.symbol(input$exp_filter_cat) >= input$exp_filter_on[1],
               !!as.symbol(input$exp_filter_cat) <= input$exp_filter_on[2])
      meta_filter_EXP <- meta_filter %>%
        filter(!!as.symbol(input$exp_filter_cat) >= input$exp_filter_on[1],
               !!as.symbol(input$exp_filter_cat) <= input$exp_filter_on[2])
    }

  } else {
    box_data <- db %>% tbl('grouped_stats') %>%
      filter(Gene %in% gene) %>%
      collect()
  }

  box_data_temp <- list()
  box_data_NA <- list()
  for (i in gene){
    #print(i)
    box_data_temp[[i]] <- box_data %>%
      filter(Gene == i) %>%
      filter(!Platform %in% c('SCRBSeq')) %>%
      mutate(Stage = factor(Stage, levels = c('Early Dev.', 'Late Dev.', 'Maturing', 'Mature'))) %>%
      group_by_at(vars(one_of(c('Gene', input$exp_plot_facet, grouping_features)))) %>%
      summarise(counts = sum(counts * cell_exp_ct) / sum(cell_exp_ct),
                cell_exp_ct = sum(cell_exp_ct, na.rm = TRUE)) %>%
      full_join(., meta_filter_EXP %>%
                  group_by_at(vars(one_of(input$exp_plot_facet, grouping_features))) %>%
                  summarise(Count = n())) %>%
      mutate(cell_exp_ct = ifelse(is.na(cell_exp_ct), 0, cell_exp_ct)) %>%
      mutate(`%` = round((cell_exp_ct / Count) * 100, 2),
             Expression = round(counts * (`%` / 100), 2)) %>%
      select_at(vars(one_of(c('Gene', grouping_features, 'cell_exp_ct', 'Count', '%', 'Expression')))) %>%
      arrange(-Expression) %>%
      dplyr::rename(`Cell # Detected` = cell_exp_ct,
             `Total Cells` = Count,
             `Mean log2(Counts + 1)` = Expression,
             `% of Cells Detected` = `%`) %>%
      mutate(`Mean log2(Counts + 1)` = case_when(is.na(`Mean log2(Counts + 1)`) ~ 0, TRUE ~ `Mean log2(Counts + 1)`)) %>%
      filter(`Total Cells` > as.integer(input$exp_filter_min_cell_number))
    # expand missing genes (nothing detected) to all genes used
    box_data_NA[[i]] <-  box_data_temp[[i]] %>% filter(is.na(Gene))
    if (nrow(box_data_NA[[i]] > 0)){
      box_data_temp[[i]] <- bind_rows(box_data_temp[[i]] %>% filter(!is.na(Gene)),
                                      box_data_NA[[i]] %>% mutate(Gene = i))
    }
  }
  box_data <- box_data_temp %>% bind_rows()
  if (nrow(box_data) == 0){
    validate(
      "Warning: no data remains after the filtering. Consider lowering the \"Minimum # Cells in Group\" threshold"
    )}

  out <- list()
  out$box_data <- box_data



  if ((input$exp_plot_flip == "Gene")){
    out$plot <- box_data %>%
      ggplot(aes(x=!!as.symbol(input$exp_plot_facet), y = !!as.symbol(input$exp_plot_ylab), color = !!as.symbol(grouping_features))) +
      geom_boxplot(color = 'black', outlier.shape = NA) +
      ggbeeswarm::geom_quasirandom(aes(size = `Total Cells`), groupOnX = TRUE) +
      cowplot::theme_cowplot() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_radius(range=c(2, 6)) +
      scale_colour_manual(values = rep(c(pals::alphabet() %>% unname()), 20)) +
      theme(legend.position="bottom") +
      facet_wrap(ncol = as.numeric(input$exp_plot_col_num), scales = 'free_y', ~Gene)
  } else {
    out$plot <- box_data %>%
      ggplot(aes(x=Gene, y = !!as.symbol(input$exp_plot_ylab), color = !!as.symbol(grouping_features))) +
      geom_boxplot(color = 'black', outlier.shape = NA) +
      ggbeeswarm::geom_quasirandom(aes(size = `Total Cells`), groupOnX = TRUE) +
      cowplot::theme_cowplot() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_radius(range=c(2, 6)) +
      scale_colour_manual(values = rep(c(pals::alphabet() %>% unname()), 20)) +
      theme(legend.position="bottom") +
      facet_wrap(ncol = as.numeric(input$exp_plot_col_num), scales = 'free_y', vars(!!as.symbol(paste(input$exp_plot_facet, collapse = ','))))
  }
  out
}
