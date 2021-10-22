# cluster properties visualized with a bar plot

# meta_filter <- data.table::fread('~/data/scEiaD_v2/metadata_filter.tsv.gz') %>%
#   mutate(CellType_predict = case_when(!is.na(TabulaMurisCellType_predict) ~ 'Tabula Muris',
#                                       is.na(CellType_predict) ~ 'Unlabelled',
#                                       TRUE ~ CellType_predict))

meta_filter_TM <- data.table::fread('~/data/scEiaD_v2/metadata_filter.tsv.gz') %>%
  mutate(CellType_predict = case_when(!is.na(TabulaMurisCellType_predict) ~ 'Tabula Muris',
                                      is.na(CellType_predict) ~ 'Unlabelled',
                                      TRUE ~ CellType_predict)) %>%
  mutate(CellType_predict = case_when(TabulaMurisCellType_predict == 'T cell' ~ 'T/NK-Cell',
                                      TabulaMurisCellType_predict == 'B cell' ~ 'B-Cell',
                                      TabulaMurisCellType_predict == 'endothelial cell' ~ 'Endothelial',
                                      TabulaMurisCellType_predict == 'epithelial cell' ~ 'Epithelial',
                                      TabulaMurisCellType_predict == 'endothelial cell' ~ 'Epithelial',
                                      TabulaMurisCellType_predict == 'keratinocyte' ~ 'Keratinocyte',
                                      TabulaMurisCellType_predict == 'blood cell' ~ 'Red Blood Cell',
                                      TabulaMurisCellType_predict == 'hepatocyte' ~ 'Hepatocyte',
                                      TabulaMurisCellType_predict == 'mesenchymal cell' ~ 'Mesenchymal',
                                      TabulaMurisCellType_predict == 'bladder cell' ~ 'Bladder',
                                      TabulaMurisCellType_predict == 'mesenchymal stem cell' ~ 'Mesenchymal (Stem)',
                                      TabulaMurisCellType_predict == 'bladder urothelial cell' ~ 'Bladder Urothelial',
                                      TabulaMurisCellType_predict == 'kidney proximal straight tubule epithelial cell' ~ 'Kidney Proximal Tubule',
                                      TabulaMurisCellType_predict == 'basal cell of epidermis' ~ 'Basal Cell',
                                      TabulaMurisCellType_predict == 'macrophage' ~ 'Macrophage',
                                      TabulaMurisCellType_predict == 'natural killer cell' ~ 'T/NK-Cell',
                                      TabulaMurisCellType_predict == 'monocyte' ~ 'Monocyte',
                                      TRUE ~ CellType_predict))

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



## CellType_predict
pdata <- meta_filter_TM %>% ungroup()%>% group_by(CellType_predict,  cluster) %>%
  summarise(Count = n()) %>%
  left_join(cat_to_color_df %>%
              filter(meta_category == 'CellType_predict') %>%
              select(value, color) %>% unique(),
            by = c('CellType_predict' = 'value')) %>%
  ungroup()
col <- pdata$color
names(col) <- pdata$CellType_predict
cluster_plot <- pdata %>%
  mutate(cluster = factor(cluster, levels = seq(1,max(pdata$cluster %>% as.numeric()),1) %>% as.character())) %>%
  ggplot(aes(x=cluster, y= Count, fill = CellType_predict)) +
  geom_bar(stat = 'identity') + cowplot::theme_cowplot() +
  scale_y_continuous(expand = c(0,0)) +
  coord_flip() +
  scale_fill_manual(values = col) +
  guides(fill=guide_legend(ncol=1))

## cluster facet
facet_plot <- pdata %>%
  group_by(cluster) %>%
  mutate(Percentage = (Count / sum(Count) * 100)) %>%
  filter(Percentage > 10)  %>%
  mutate(cluster = factor(cluster, levels = seq(1,max(pdata$cluster %>% as.numeric()),1) %>% as.character())) %>%
  ggplot(aes(x=CellType_predict, y=Count, fill = CellType_predict)) +
  geom_bar(stat = 'identity') + cowplot::theme_cowplot() +
  scale_y_continuous(expand = c(0,0)) +
  coord_flip() +
  scale_fill_manual(values = col) +
  facet_wrap(~cluster, scales = 'free', ncol = 5) +
  theme(legend.position = "none")

## phase
pdata <- meta_filter_TM %>% ungroup()%>% group_by(Phase,  cluster) %>%
  summarise(Count = n()) %>%
  left_join(cat_to_color_df %>%
              filter(meta_category == 'Phase') %>%
              select(value, color) %>% unique(),
            by = c('Phase' = 'value')) %>%
  ungroup()
col <- pdata$color
names(col) <- pdata$Phase
phase_plot <- pdata %>%
  mutate(cluster = factor(cluster, levels = seq(1,max(pdata$cluster %>% as.numeric()),1) %>% as.character())) %>%
  ggplot(aes(x=cluster, y= Count, fill = Phase)) +
  geom_bar(stat = 'identity') + cowplot::theme_cowplot() +
  coord_flip() +
  scale_fill_manual(values = col)

## organism
pdata <- meta_filter_TM %>% ungroup()%>% group_by(organism,  cluster) %>%
  summarise(Count = n()) %>%
  left_join(cat_to_color_df %>%
              filter(meta_category == 'organism') %>%
              select(value, color) %>% unique(),
            by = c('organism' = 'value')) %>%
  ungroup()
col <- pdata$color
names(col) <- pdata$organism
org_plot <- pdata %>%
  mutate(cluster = factor(cluster, levels = seq(1,max(pdata$cluster %>% as.numeric()),1) %>% as.character())) %>%
  ggplot(aes(x=cluster, y= Count, fill = organism)) +
  geom_bar(stat = 'identity') + cowplot::theme_cowplot() +
  coord_flip() +
  scale_fill_manual(values = col)




#
#
# meta_filter %>% group_by(cluster, organism) %>%
#   summarise(Count = n()) %>%
#   mutate(Percent = (Count/sum(Count)) * 100) %>% mutate(value = case_when(Percent > 15.0 ~ organism, TRUE ~ '')) %>% ungroup() %>% group_by(cluster) %>% summarise(value = paste(value, collapse = ' ')) %>% mutate(value = gsub('\\s\\s+', ' ', value))
