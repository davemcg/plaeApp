box_data <- scEiaD_2020_v01 %>% tbl('grouped_stats') %>%
  filter(organism == 'Homo sapiens') %>%
  collect()

box_data2 <- box_data %>%
  filter(!Platform %in% c('SCRBSeq')) %>%
  mutate(Stage = factor(Stage, levels = c('Early Dev.', 'Late Dev.', 'Maturing', 'Mature'))) %>%
  #filter(!is.na(!!as.symbol(grouping_features))) %>%
  group_by_at(vars(one_of(c('Gene', 'CellType_predict', 'organism', 'study_accession', 'Stage')))) %>%
  summarise(counts = sum(counts * cell_exp_ct) / sum(cell_exp_ct),
            cell_exp_ct = sum(cell_exp_ct, na.rm = TRUE)) %>%
  full_join(., meta_filter %>%
              filter(organism == 'Homo sapiens') %>%
              group_by_at(vars(one_of(c('organism', 'study_accession', 'Stage')))) %>%
              summarise(Count = n())) %>%
  mutate(cell_exp_ct = ifelse(is.na(cell_exp_ct), 0, cell_exp_ct)) %>%
  mutate(`%` = round((cell_exp_ct / Count) * 100, 2),
         Expression = round(counts * (`%` / 100), 2)) %>%
  select_at(vars(one_of(c('Gene', 'CellType_predict', 'organism', 'study_accession',
                          'Stage', 'cell_exp_ct', 'Count', '%', 'Expression'
                          )))) %>%
  arrange(-Expression) %>%
  rename(`Cell # Detected` = cell_exp_ct,
         `Total Cells` = Count,
         `Mean log2(Counts + 1)` = Expression,
         `% of Cells Detected` = `%`) %>%
  mutate(`Mean log2(Counts + 1)` = case_when(is.na(`Mean log2(Counts + 1)`) ~ 0, TRUE ~ `Mean log2(Counts + 1)`)) %>%
  filter(`Total Cells` > as.integer(100))
# expand missing genes (nothing detected) to all genes used
box_data_NA <- box_data2 %>% filter(is.na(Gene))
if (nrow(box_data_NA > 0)){
  box_data_NA_list <- list()
  for (i in gene){
    box_data_NA_list <- box_data_NA %>% mutate(Gene = i)
  }
  box_data2 <- bind_rows(box_data2 %>% filter(!is.na(Gene)),
                        box_data_NA_list %>% bind_rows())
}
scEiaD_CT_table <- box_data2 %>% filter(!is.na(Gene) )

save(scEiaD_CT_table, file = '/Volumes/ARC168/EiaD_databases/scEiaD_CT_table.Rdata')
