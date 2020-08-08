library(ggplot2)
library(Cairo)
library(scattermore)
library(pals)
library(ggrepel)
library(patchwork)
library(purrr)
library(pool)
library(RSQLite)
library(dplyr)
library(magick)
library(stringr)

make_dotplot <- function(input, db, meta_filter){
  ### this makes it a little easier to test
  # input <- list()
  # input[['dotplot_Gene']] <- c('RHO','WIF1', 'CABP5', 'AIF1', 'ARR3', 'ONECUT1', 'GRIK1', 'GAD1', 'POU4F2')#c('RHO', 'CRX')
  # input[['dotplot_groups']] <- 'CellType_predict'
  # input[['dotplot_filter_cat']] <- 'CellType_predict'
  # input[['dotplot_filter_on']] <- 'Macrophage'
  gene <- input$dotplot_Gene
  grouping_features <- input$dotplot_groups
  validate(
    need(input$dotplot_Gene !='', "Please select at least one gene"),
    need(input$dotplot_groups != '', "Please select at least one grouping feature")
  )

  if (input$dotplot_filter_cat != ''){
    #NOTE: dotplot_filter_cat options on app do not match database
    validate(need(input$dotplot_filter_on!= '', 'Select a value to filter on' ))
    dotplot_data <- db %>% tbl('grouped_stats') %>%
      filter(Gene %in% gene) %>%
      collect() %>%
      filter(!!as.symbol(input$dotplot_filter_cat) %in% input$dotplot_filter_on)

  } else {
    dotplot_data <- db %>% tbl('grouped_stats') %>%
      filter(Gene %in% gene) %>%
      collect()
  }

  dotplot_data <- dotplot_data %>%
    filter(Gene %in% gene) %>%
    group_by_at(vars(one_of(c('Gene', grouping_features)))) %>%
    summarise(cell_exp_ct = sum(cell_exp_ct, na.rm = TRUE),
              cpm = mean(cpm)) %>%
    collect() %>%
    tidyr::drop_na() %>%
    full_join(., meta_filter %>%
                group_by_at(vars(one_of(grouping_features))) %>%
                summarise(Count = n())) %>%
    mutate(`%` = round((cell_exp_ct / Count) * 100, 2),
           Expression = cpm * (`%` / 100)) %>%
    filter(!is.na(Count),
           `%` > 2) %>%
    select_at(vars(one_of(c('Gene', grouping_features, 'cell_exp_ct', 'Count', '%', 'Expression')))) %>%
    filter(!is.na(Gene))
  if (length(grouping_features) == 2){
    colnames(dotplot_data)[c(2,3)] <- c("Group1","Group2")
    dotplot_data <- dotplot_data %>%
      mutate(Column = paste(Group1,Group2, sep = ':'),
             Column = case_when(grepl('^\\d', Column) ~ paste0('X', Column),
                                TRUE ~ Column))
  } else {
    colnames(dotplot_data)[c(2)] <- c("Group1")
    dotplot_data <- dotplot_data %>%
      mutate(Column = Group1,
             Column = case_when(grepl('^\\d', Column) ~ paste0('X', Column),
                                TRUE ~ Column))
  }

  # cluster
  # make data square to calculate euclidean distance
  # NOTE: pivot_wider only in tidyr >=1.0.0
  mat <- dotplot_data %>%
    select(Gene,Group1, Column, Expression) %>%  # drop unused columns to faciliate widening
    tidyr::pivot_wider(names_from = Column, values_from = Expression) %>%
    mutate(id_col = paste(Gene, Group1, sep = '-' )) %>%
    select(Gene, Group1, id_col, everything()) %>%
    as.data.frame() # make df as tibbles -> matrix annoying
  row.names(mat) <- mat$id_col  # put gene in `row`
  mat <- mat %>% select(-colnames(mat)[1:3]) #drop gene column as now in rows; Has to be this way for edge case  of having only one row
  mat[is.na(mat)] <- 0

  validate(
    need(nrow(mat)!=0 | is.null(mat), 'Gene not expressed in current selection')
  )
  if(nrow(mat) ==1 | ncol(mat) ==1) {# in some cases
    h_clust <- list()
    h_clust[['labels']] <- rownames(mat)
    h_clust[['order']] <- 1:nrow(mat)
    v_clust <- list()
    v_clust[['labels']] <- colnames(mat)
    v_clust[['order']] <- 1:ncol(mat)
  }else{
    h_clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix
    v_clust <- hclust(dist(mat %>% as.matrix() %>% t())) # hclust with distance matrix
  }
  dotplot <- dotplot_data %>%
    mutate(Gene = factor(Gene, levels = {h_clust$labels[h_clust$order]} %>%
                           str_split('-') %>% sapply(function(x) x[1]) %>% unique ),
           Column = factor(Column %>% str_replace_all(':',' ') , levels = v_clust$labels[v_clust$order] %>%
                             str_replace_all(':', ' '))) %>%
    ggplot(aes(x=Column, y = Gene, size = `%`, color = Expression)) +
    geom_point() +
    cowplot::theme_cowplot() +
    scale_color_viridis_c(option = 'magma') +
    theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab('') + xlab('') +
    scale_y_discrete(position = 'right')+
    coord_flip()

  dp_legened <- get_legend(dotplot)
  dotplot <- dotplot+
    theme(axis.ticks = element_blank(), axis.text.x= element_text(angle = -45), legend.position = 'none')
  # labels
  order <- left_join(tibble::enframe(v_clust$labels[v_clust$order] , value = 'Column'),
                     dotplot_data) %>%
    mutate(Column = Column %>% str_replace_all(':', ' '))

  group1_labels <- order %>%
    ggplot() +
    geom_tile(aes(x = Column, y = 1, fill = Group1)) +
    scale_fill_manual(name = grouping_features[1],
                      values = rep(c(pals::alphabet2() %>% unname(),
                                     pals::alphabet() %>% unname()), times = 10)) +
    theme_nothing() +
    aplot::xlim2(dotplot) +
    coord_flip()
  # remove legend if same size as matrix as will be HUGE
  if ((dotplot_data$Group1 %>% unique() %>% length()) == ncol(mat)) {
    group1_legend <- NULL
  } else {
    group1_legend <- plot_grid(get_legend(group1_labels + theme(legend.position="bottom")))
  }

  if (length(grouping_features) == 2){
    group2_labels <- order %>%
      ggplot() +
      geom_tile(aes(x = Column, y = 1, fill = Group2)) +
      scale_fill_manual(name = grouping_features[2],
                        values = rep(c(pals::alphabet() %>% unname(),
                                       pals::alphabet2() %>% unname()), times = 10)) +
      theme_nothing() +
      aplot::xlim2(dotplot) +
      coord_flip()
    if ((dotplot_data$Group2 %>% unique() %>% length()) == ncol(mat)) {
      group2_legend <- NULL
    } else {
      group2_legend <- plot_grid(get_legend(group2_labels + theme(legend.position="bottom")))
    }

    # group2_labels +
    #   group1_labels +
    #   dotplot +
    #   group2_legend +
    #   group1_legend +
    #   plot_layout(ncol= 1, heights = c(0.1,0.1, 1, 0.1, 0.2))


    top <- dotplot |
      group1_labels |
      group2_labels|
      plot_spacer() |
      dp_legened |plot_spacer() |plot_layout(nrow = 1, widths = c(1,.05,.05,.05,.05,.05))
    bottom <- (group1_legend +plot_spacer())/plot_spacer()/group2_legend
    top/plot_spacer()/bottom +plot_layout(ncol = 1, heights = c(1,.05, .2))

  } else {
    # group1_labels +
    #   dotplot +
    #   group1_legend +
    #   plot_layout(nrow= 1, heights = c(0.1, 1, 0.1))
    top <- dotplot |
      group1_labels |
      plot_spacer() |
      dp_legened |plot_spacer() | plot_layout(nrow = 1, widths = c(1,.05,.05, .05,.05))
    top/plot_spacer() / group1_legend +plot_layout(ncol = 1, heights = c(1,.05, .2))
  }
}
