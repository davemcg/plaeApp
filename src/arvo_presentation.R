# ARVO presentation
# run "build_tables.R" to get data
t2 <- cell_info_labels %>%
  mutate(CellType = gsub('AC/HC_Precurs', 'Amacrine/Horizontal Precursors', CellType)) %>%
  #mutate(CellType = gsub("Rod Bipolar Cells", "Bipolar Cells", CellType)) %>%
  filter(!is.na(CellType),
         !is.na(study_accession),
         !CellType %in% c('Doublet', 'Doublets'),
         !grepl('RPE|Vascul', CellType)) %>%
  mutate(organism = case_when(grepl('Homo', organism) ~ 'HS',
                              grepl('Mus', organism) ~ 'MM',
                              TRUE ~ 'MF')) %>%
  group_by(CellType,organism, study_accession) %>%
  summarise(Count = n()) %>% filter(Count > 50) %>%
  summarise(Studies = length(study_accession), Count = sum(Count)) %>%
  summarise(Species = paste(organism, collapse = ', '), Studies = sum(Studies), Count = sum(Count)) %>%
  arrange(-Count)

t3 <- meta %>%
  mutate(CellType = gsub('AC/HC_Precurs', 'Amacrine/Horizontal Precursors', CellType_predict)) %>%
  #mutate(CellType = gsub("Rod Bipolar Cells", "Bipolar Cells", CellType)) %>%
  filter(!is.na(CellType),
         !is.na(study_accession)) %>%
  mutate(organism = case_when(grepl('Homo', organism) ~ 'HS',
                              grepl('Mus', organism) ~ 'MM',
                              TRUE ~ 'MF')) %>%
  group_by(CellType,organism, study_accession) %>%
  summarise(Count = n()) %>% filter(Count > 50) %>%
  summarise(Studies = length(study_accession), Count = sum(Count)) %>%
  summarise(Species = paste(organism, collapse = ', '), Studies = sum(Studies), Count = sum(Count)) %>%
  arrange(-Count)


t2 %>%
  left_join(t3, by = 'CellType') %>%
  arrange(CellType) %>% filter(!is.na(`Studies.y`)) %>%
  select(CellType, `Studies.x`,`Studies.y`) %>%
  pivot_longer(-CellType, names_to = 'Labels') %>%
  mutate(`Labels from` = case_when(Labels == 'Studies.x' ~ 'Publication', TRUE ~ 'Machine Learning')) %>%
  mutate(`Labels from` = factor(`Labels from`, levels = c('Publication','Machine Learning'))) %>%

  filter(`Labels from` == 'Publication') %>%


  ggplot(aes(x=CellType,y=value,color=`Labels from`, group = CellType)) +
  geom_point(size = 2) +
  coord_flip() +
  cowplot::theme_cowplot() +
  geom_line(color = 'Black') +
  ylab("Studies with more than 50 cells") +
  xlab("Cell Type") +
  scale_y_continuous(limits = c(0,18)) +
  scale_color_brewer(palette = 'Set1')


t2 %>%
  left_join(t3, by = 'CellType') %>%
  arrange(CellType) %>% filter(!is.na(`Studies.y`)) %>%
  select(CellType, `Studies.x`,`Studies.y`) %>%
  pivot_longer(-CellType, names_to = 'Labels') %>%
  mutate(`Labels from` = case_when(Labels == 'Studies.x' ~ 'Publication', TRUE ~ 'Machine Learning')) %>%
  mutate(`Labels from` = factor(`Labels from`, levels = c('Publication','Machine Learning'))) %>%

  ggplot(aes(x=CellType,y=value,color=`Labels from`, group = CellType)) +
  geom_point(size = 2) +
  coord_flip() +
  cowplot::theme_cowplot() +
  geom_line(color = 'Black') +
  ylab("Study count from studies with more than 50 cells") +
  xlab("Cell Type") +
  scale_y_continuous(limits = c(0,18)) +
  scale_color_brewer(palette = 'Set1')




t2 %>%
  left_join(t3, by = 'CellType') %>%
  arrange(CellType) %>% filter(!is.na(`Count.y`)) %>%
  select(CellType, `Count.x`,`Count.y`) %>%
  pivot_longer(-CellType, names_to = 'Labels') %>%
  mutate(`Labels from` = case_when(Labels == 'Count.x' ~ 'Publication', TRUE ~ 'Machine Learning')) %>%
  mutate(`Labels from` = factor(`Labels from`, levels = c('Publication','Machine Learning'))) %>%

  filter(`Labels from` == 'Publication') %>%

  ggplot(aes(x=CellType,y=value,color=`Labels from`, group = CellType)) +
  geom_point(size = 2) +
  coord_flip() +
  cowplot::theme_cowplot() +
  geom_line(color = 'Black') +
  ylab("Cell Type Count from Studies with more than 50 cells") +
  xlab("Cell Type") +
  scale_y_continuous(limits = c(0,144000), labels = scales::comma)+
  scale_color_brewer(palette = 'Set1')


t2 %>%
  left_join(t3, by = 'CellType') %>%
  arrange(CellType) %>% filter(!is.na(`Count.y`)) %>%
  select(CellType, `Count.x`,`Count.y`) %>%
  pivot_longer(-CellType, names_to = 'Labels') %>%
  mutate(`Labels from` = case_when(Labels == 'Count.x' ~ 'Publication', TRUE ~ 'Machine Learning')) %>%
  mutate(`Labels from` = factor(`Labels from`, levels = c('Publication','Machine Learning'))) %>%
  ggplot(aes(x=CellType,y=value,color=`Labels from`, group = CellType)) +
  geom_point(size = 2) +
  coord_flip() +
  cowplot::theme_cowplot() +
  geom_line(color = 'Black') +
  ylab("Cell type count from studies with more than 50 cells") +
  xlab("Cell Type") +
  scale_y_continuous(limits = c(0,144000), labels = scales::comma)+
  scale_color_brewer(palette = 'Set1')
