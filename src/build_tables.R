library(tidyverse)
library(formattable)
library(webshot)
library(htmltools)

# all cells
load('~/data/scEiaD_v3/cell_info_labelled.Rdata')
# pre mt filtering
srt <- data.table::fread('~/git/scEiaD/data/sample_run_layout_organism_tech_biosample_organ_2021_06_07.tsv')
mito <- data.table::fread('~/data/scEiaD_v3/QC.tsv.gz')
mito <- mito %>% left_join(srt %>% select(sample_accession, SX = Source) %>% unique(), by = 'sample_accession') %>% filter(SX %in% c('iPSC','Tissue'))
# load('/Volumes/data/projects/nei/mcgaughey/massive_integrated_eye_scRNA/fastMNN_umap_full.Rdata')
study_meta <- read_tsv('~/git/scEiaD/data/GEO_Study_Level_Metadata.tsv')
meta <- fst::read_fst('~/data/scEiaD_v3/meta_filter.fst')

stats <- meta %>% group_by(study_accession, Platform) %>% summarise(Counts = n())


# Colors from https://personal.sron.nl/~pault/#sec:qualitative
color_bar_factor <- formatter("span",
                              style = function(x) style(
                                display = "block",
                                color = "black",
                                border.radius = "4px",
                                background = c("#44BB99", "#EE8866", "#BBCC33", "#dc74e8")[factor(as.character(x))]))

color_bar_factor2 <- formatter("span",
                               style = function(x) style(
                                 display = "block",
                                 color = "black",
                                 border.radius = "4px",
                                 background = c("#DDDDDD", "#FFAABB")[factor(as.character(x))]))

color_bar_factor3 <- formatter("span", width = "5px",
                               style = function(x) style(
                                 display = "block",
                                 color = "black",
                                 border.radius = "4px",
                                 background = type_val[factor(as.character(x))]))

export_formattable <- function(f, file, width = "100%", height = NULL,
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,zoom = 4,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}

post <- stats %>% rename(`Post QC<br/>Count` = Counts) %>%
  select(study_accession, Platform, `Post QC<br/>Count`) %>%
  group_by(study_accession, Platform) %>%
  summarise(`Post QC<br/>Count` = sum(`Post QC<br/>Count`)) %>%
  mutate(study_accession = case_when(study_accession == 'OGVFB_Hufnagel_iPSC_RPE' ~ 'SRP329495',
                                                                                              TRUE ~ study_accession))

table01 <- meta %>%
  mutate(study_accession = case_when(study_accession == 'OGVFB_Hufnagel_iPSC_RPE' ~ 'SRP329495',
                                      TRUE ~ study_accession)) %>%
  select(-PMID, -Citation, -GSE, -Design, -Summary) %>%
  left_join(study_meta, by = 'study_accession') %>%
  #left_join(cell_info_labels %>% select(barcode = value, CellType)) %>%
  #filter(Source != 'Cell Culture', Source != 'Organoid') %>%
  filter(!grepl('Bharti',study_accession)) %>%
  mutate(PMID = as.character(PMID)) %>%
  mutate(Citation = case_when(is.na(Citation) ~ '',
                              TRUE ~ paste0(substr(Citation, 1, 30), ' ...')),
         PMID = case_when(is.na(PMID) ~ '',
                          TRUE ~ PMID),
         `SRA Accession` = study_accession,
         Labels = case_when(!is.na(CellType) ~ 'Yes',
                            TRUE ~ 'No')) %>%
  group_by(Citation, PMID, `SRA Accession`, organism, Platform) %>%
  summarise(Count = n(), Labels = unique(Labels) %>% sort() %>% tail(1)) %>%
  arrange(organism, -Count) %>%
  filter(Count > 0)

#table01$`Post QC<br/>Count`[is.na(table01$`Post QC<br/>Count`)] <- 0
formattable_01 <- table01 %>%
  mutate(PMID = case_when(!grepl('doi', PMID) ~ glue::glue('<a href = https://pubmed.ncbi.nlm.nih.gov/{PMID}>{PMID}</a>') %>% as.character(),
                          TRUE ~ '<a <href = https://doi.org/10.1101/774950>bioRxiv 774950</a>')) %>%
  formattable(., list(Count = normalize_bar("lightblue"),
                      `Post QC<br/>Count` = normalize_bar("lightblue"),
                      organism = color_bar_factor,
                      Labels = color_bar_factor2,
                      Platform = formatter("span", style = x ~ ifelse(!grepl('SMART|C1|SCRB', x),
                                                                      style(font.weight = "bold"), NA))))
formattable_01
write_file(format_table(formattable_01), path = 'inst/app/www/table_01.html')


table02 <- meta %>%
  filter(!is.na(CellType)) %>%
  mutate(CellType = gsub('AC/HC_Precurs', 'Amacrine/Horizontal Precursors', CellType)) %>%
  #mutate(CellType = gsub("Rod Bipolar Cells", "Bipolar Cells", CellType)) %>%
  filter(!is.na(CellType),
         !is.na(study_accession),
         !CellType %in% c('Doublet', 'Doublets'),
         !grepl('RPE|Vascul', CellType)) %>%
  mutate(organism = case_when(grepl('Homo', organism) ~ 'HS',
                              grepl('Mus', organism) ~ 'MM',
                              grepl('Gal', organism) ~ 'GG',
                              TRUE ~ 'MF')) %>%
  group_by(CellType,organism, study_accession) %>%
  summarise(Count = n()) %>%
  filter(Count > 50) %>%
  summarise(Studies = length(study_accession), Count = sum(Count)) %>%
  summarise(Species = paste(organism, collapse = ', '), Studies = sum(Studies), Count = sum(Count)) %>%
  arrange(-Count) %>%
  formattable(., list(Count = normalize_bar("lightblue")), width = "1")
table02
formattable_02 <- table02
write_file(format_table(formattable_02), path = 'inst/app/www/table_02.html')



table03 <- meta %>%
  mutate(CellType = gsub('AC/HC_Precurs', 'Amacrine/Horizontal Precursors', CellType_predict)) %>%
  #mutate(CellType = gsub("Rod Bipolar Cells", "Bipolar Cells", CellType)) %>%
  filter(!is.na(CellType),
         !is.na(study_accession)) %>%
  mutate(organism = case_when(grepl('Homo', organism) ~ 'HS',
                              grepl('Mus', organism) ~ 'MM',
                              grepl('Gal', organism) ~ 'GG',
                              TRUE ~ 'MF')) %>%
  group_by(CellType,organism, study_accession) %>%
  summarise(Count = n()) %>%
  filter(Count > 50) %>%
  summarise(Studies = length(study_accession), Count = sum(Count)) %>%
  summarise(Species = paste(organism, collapse = ', '), Studies = sum(Studies), Count = sum(Count)) %>%
  arrange(-Count) %>%
  formattable(., list(Count = normalize_bar("lightblue")), width = "1")
table03
formattable_03 <- table03
write_file(format_table(formattable_03), path = 'inst/app/www/table_03.html')
