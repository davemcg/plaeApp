# cluster properties visualized with a bar plot


## CellType_predict
pdata <- meta_filter %>% group_by(CellType_predict,  cluster) %>%
  summarise(Count = n()) %>%
  left_join(cat_to_color_df %>%
              filter(meta_category == 'CellType_predict') %>%
              select(value, color) %>% unique(),
            by = c('CellType_predict' = 'value'))
col <- pdata$color
names(col) <- pdata$CellType_predict
pdata %>%
  mutate(cluster = as.integer(cluster)) %>%
  ggplot(aes(x=cluster, y= Count, fill = CellType_predict)) +
  geom_bar(stat = 'identity') + cowplot::theme_cowplot() +
  coord_flip() +
  scale_fill_manual(values = col)


## organism

pdata <- meta_filter %>% group_by(cluster, organism) %>%
  summarise(Count = n()) %>%
  mutate(Percent = mutate(Percent = (Count/sum(Count)) * 100)) %>%
  left_join(cat_to_color_df %>%
              filter(meta_category == 'organism') %>%
              select(value, color) %>% unique(),
            by = c('organism' = 'value'))
col <- pdata$color
names(col) <- pdata$organism
pdata %>%
  mutate(cluster = as.integer(cluster)) %>%
  ggplot(aes(x=cluster, y= Count, fill = organism)) +
  geom_bar(stat = 'identity') + cowplot::theme_cowplot() +
  coord_flip() +
  scale_fill_manual(values = col)






meta_filter %>% group_by(cluster, organism) %>%
  summarise(Count = n()) %>%
  mutate(Percent = (Count/sum(Count)) * 100) %>% mutate(value = case_when(Percent > 15.0 ~ organism, TRUE ~ '')) %>% ungroup() %>% group_by(cluster) %>% summarise(value = paste(value, collapse = ' ')) %>% mutate(value = gsub('\\s\\s+', ' ', value))
