make_temporal_plot <- function(input, db, meta_filter) {
  # gdata::keep(meta_filter, scEiaD_2020_v01, sure = T)
  # input <- list()
  # db <- scEiaD_2020_v01
  # input[['temporal_gene']] <- c('PAX6','POU4F2', 'RHO')
  # input[['temporal_group']]  <- 'CellType'
  # input[['temporal_y_val']] <- 'Mean CPM'
  cat(file=stderr(), paste0(Sys.time(), ' Temporal Plot Call\n'))
  gene <- input$temporal_gene
  grouping <- input$temporal_group
  y_val <- input$temporal_y_val
  if (grouping == 'CellType (predict)'){grouping <- 'CellType_predict'}
  if (y_val == 'Mean CPM') {y_val <- 'cpm'} else {y_val <- 'Ratio'}
  meta_data <- meta_filter %>%
    mutate(Stage = replace(Stage, Age >=17, 'Adult')) %>%
    group_by(organism, !!as.symbol(grouping), Stage) %>%
    summarise(full_count = n()) %>%
    mutate(Stage = factor(Stage, levels = c('Early', 'Late', 'Adult')))

  temporal_data <- db %>% tbl('cpm') %>%
    filter(Gene %in% gene) %>%
    collect() %>%
    mutate(cpm = cpm - min(cpm) + 1) %>%
    left_join(., meta_filter, by = 'Barcode') %>%
    mutate(Stage = replace(Stage, Age >=17, 'Adult')) %>%
    filter(!is.na(!!as.symbol(grouping)), !grepl('Doub|RPE', !!as.symbol(grouping))) %>%
    group_by(organism, !!as.symbol(grouping), Stage, Gene) %>%
    summarise(cpm = mean(cpm), count = n()) %>%
    right_join(., meta_data) %>%
    mutate(count = ifelse(is.na(count), 0, count)) %>%
    mutate(Ratio = count/full_count) %>%
    filter(!is.na(!!as.symbol(grouping)), !is.na(Gene)) %>%
    ungroup() %>%
    mutate(Stage = factor(Stage, levels = c('Early', 'Late', 'Adult')))



  suppressWarnings(temporal_human <-  temporal_data %>%
                     filter(organism == 'Homo sapiens') %>%
                     ggplot(aes(x=Stage, y = !!as.symbol(y_val), color = Gene, group=Gene)) +
                     geom_point(stat = 'identity') +
                     ggtitle('Human') +
                     geom_line() + ylab(input$temporal_y_val) +
                     cowplot::theme_cowplot() + xlab('Age (days from birth)') +
                     facet_wrap(vars(!!as.symbol(grouping))) +
                     theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
                     scale_colour_manual(values = rep(c(pals::alphabet() %>% unname()))))
  suppressWarnings(temporal_mouse <- temporal_data %>%
                     filter(organism == 'Mus musculus') %>%
                     ggplot(aes(x=Stage, y = !!as.symbol(y_val), color = Gene, group =Gene)) +
                     geom_point(stat = 'identity') +
                     ggtitle('Mouse') +
                     geom_line() +
                     ylab(input$temporal_y_val) +
                     cowplot::theme_cowplot() + xlab('Age (days from birth)') +
                     facet_wrap(vars(!!as.symbol(grouping))) +
                     theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
                     scale_colour_manual(values = rep(c(pals::alphabet() %>% unname()))))
  # draw plot
  temporal_human + temporal_mouse + plot_layout(ncol = 1)
}
