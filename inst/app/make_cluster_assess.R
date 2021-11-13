# cluster properties visualized with a bar plot

# meta_filter <- data.table::fread('~/data/scEiaD_v2/metadata_filter.tsv.gz') %>%
#   mutate(CellType_predict = case_when(!is.na(TabulaMurisCellType_predict) ~ 'Tabula Muris',
#                                       is.na(CellType_predict) ~ 'Unlabelled',
#                                       TRUE ~ CellType_predict))

meta_filter_TM <- fst::read_fst('~/data/scEiaD_v3/2021_11_09_meta_filter.fst') %>%
  mutate(CellType_predict = case_when(!is.na(TabulaMurisCellType_predict) ~ 'Tabula Muris',
                                      is.na(CellType_predict) ~ 'Unlabelled',
                                      TRUE ~ CellType_predict)) #%>%
  # mutate(CellType_predict = case_when(TabulaMurisCellType_predict == 'T cell' ~ 'T/NK-Cell',
  #                                     TabulaMurisCellType_predict == 'B cell' ~ 'B-Cell',
  #                                     TabulaMurisCellType_predict == 'endothelial cell' ~ 'Endothelial',
  #                                     TabulaMurisCellType_predict == 'epithelial cell' ~ 'Epithelial',
  #                                     TabulaMurisCellType_predict == 'endothelial cell' ~ 'Epithelial',
  #                                     TabulaMurisCellType_predict == 'keratinocyte' ~ 'Keratinocyte',
  #                                     TabulaMurisCellType_predict == 'blood cell' ~ 'Red Blood Cell',
  #                                     TabulaMurisCellType_predict == 'hepatocyte' ~ 'Hepatocyte',
  #                                     TabulaMurisCellType_predict == 'mesenchymal cell' ~ 'Mesenchymal',
  #                                     TabulaMurisCellType_predict == 'bladder cell' ~ 'Bladder',
  #                                     TabulaMurisCellType_predict == 'mesenchymal stem cell' ~ 'Mesenchymal (Stem)',
  #                                     TabulaMurisCellType_predict == 'bladder urothelial cell' ~ 'Bladder Urothelial',
  #                                     TabulaMurisCellType_predict == 'kidney proximal straight tubule epithelial cell' ~ 'Kidney Proximal Tubule',
  #                                     TabulaMurisCellType_predict == 'basal cell of epidermis' ~ 'Basal Cell',
  #                                     TabulaMurisCellType_predict == 'macrophage' ~ 'Macrophage',
  #                                     TabulaMurisCellType_predict == 'natural killer cell' ~ 'T/NK-Cell',
  #                                     TabulaMurisCellType_predict == 'monocyte' ~ 'Monocyte',
  #                                     TRUE ~ CellType_predict))

meta_filter_TM <- meta_filter_TM %>% mutate(Age = as.character(Age),
                                            cluster = as.character(cluster),
                                            SubCellType = tidyr::replace_na(SubCellType, 'None'),
                                      subcluster = as.character(subcluster))
map_color <- function(column, meta_filter){
  #master_colorlist <- c(pals::polychrome()[3:length(pals::polychrome())], pals::alphabet2())
  #master_colorlist <- c(pals::glasbey()[-c(3,4,8,18)],pals::alphabet2()[-c(5,7,8,9,23,24)])
  master_colorlist <- c(pals::cols25()[1:23],pals::alphabet())
  values <- meta_filter %>% pull(!!column) %>% unique %>% sort
  if(length(values) > length(master_colorlist) ){
    r= round(length(values) / length(master_colorlist)) +1
    master_colorlist <- rep(master_colorlist, r)
  }
  colors <- master_colorlist[1:length(values)]
  return(tibble(meta_category = column,value = values, color=colors))

}
categorical_columns <- c("organism","cluster","CellType","CellType_predict", "Phase")
cat_to_color_df <- lapply(categorical_columns, function(col) map_color(col, meta_filter_TM)) %>% bind_rows()



cluster_facet_maker <- function(meta_filter,
                                faceting_on = 'CellType_predict',
                                filter_on = NA,
                                perc_cutoff = 10){
  pdata <- meta_filter %>%
    ungroup() %>%
    group_by(!!as.symbol(faceting_on),  cluster) %>%
    summarise(Count = n())
  colnames(pdata)[1] <- 'value'
  pdata <- pdata %>%
    left_join(cat_to_color_df %>%
                filter(meta_category == faceting_on) %>%
                select(value, color) %>% unique(),
              by =  'value') %>%
    ungroup()
  col <- pdata$color
  names(col) <- pdata$value

  if (!is.na(filter_on)){
    pdata <- pdata %>% filter(cluster == filter_on)
  }

  ## cluster facet
  facet_plot <- pdata %>%
    group_by(cluster) %>%
    mutate(Percentage = (Count / sum(Count) * 100)) %>%
    filter(Percentage > perc_cutoff)  %>%
    mutate(cluster = factor(cluster, levels = seq(1,max(pdata$cluster %>% as.numeric()),1) %>% as.character())) %>%
    ggplot(aes(x=value, y=Count, fill = value)) +
    geom_bar(stat = 'identity') + cowplot::theme_cowplot() +
    scale_y_continuous(expand = c(0,0)) +
    coord_flip() +
    scale_fill_manual(values = col) +
    facet_wrap(~cluster, scales = 'free', ncol = 6) +
    theme(legend.position = "none") +
    xlab(faceting_on)
  facet_plot
}

cluster_facet_maker(meta_filter_TM, 'CellType_predict')

cowplot::plot_grid(
  cluster_facet_maker(meta_filter_TM, 'CellType_predict', 1),
  cluster_facet_maker(meta_filter_TM, 'Phase', 1),
  cluster_facet_maker(meta_filter_TM, 'study_accession', 1, 5),
  cluster_facet_maker(meta_filter_TM, 'organism', 1),
  align = 'hv',
  ncol = 1)


cowplot::plot_grid(
  cluster_facet_maker(meta_filter_TM, 'CellType_predict',66),
  cluster_facet_maker(meta_filter_TM, 'Phase', 66),
  cluster_facet_maker(meta_filter_TM, 'study_accession', 66, 5),
  cluster_facet_maker(meta_filter_TM, 'organism', 66),
  align = 'hv',
  ncol = 1)

