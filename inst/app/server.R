# server.R
# if running locally: setwd('inst/app')
time <- Sys.time()
cat(file = stderr(), 'Server Go!\n')
#options(shiny.trace=TRUE)
options(shiny.sanitize.errors = FALSE)

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
library(shinyalert)
library(fst)

scEiaD_2020_v01 <- dbPool(drv = SQLite(), dbname ="/Volumes/McGaughey_S/data/scEiaD/MOARTABLES__anthology_limmaFALSE___5000-counts-universe-batch-scVIprojection-6-15-0.1-50-20.sqlite", idleTimeout = 3600000)
#scEiaD_2020_v01 <- dbPool(drv = SQLite(), dbname = "/data/swamyvs/plaeApp/sql_08132020.sqlite", idleTimeout = 3600000)

# # find "common" tabula muris cell type labels to move over
# meta_filter %>%
#   group_by(cluster, TabulaMurisCellType_predict) %>%
#   summarise(Count = n()) %>%
#   mutate(Percentage = (Count / sum(Count) * 100)) %>%
#   filter(Percentage > 10) %>%
#   arrange(cluster) %>%
#   data.frame() %>%
#   filter(!is.na(TabulaMurisCellType_predict))

x_dir <- 1
y_dir <- 1
# load('~/data/scEiaD_v2/n_features-5000__transform-counts__partition-universe__covariate-batch__method-scVIprojection__dims-6__epochs-15__dist-0.1__neighbors-50__knn-20__pacmap.Rdata')
meta_filter <- read_fst('~/data/scEiaD_v2//2021_10_22_meta_filter.fst') %>%
  as_tibble() %>%
  mutate(CellType_predict = case_when(!is.na(TabulaMurisCellType_predict) ~ 'Tabula Muris',
                                      is.na(CellType_predict) ~ 'Unlabelled',
                                      TRUE ~ CellType_predict)) %>%
  mutate(UMAP_a = UMAP_2 * x_dir,
         UMAP_b = UMAP_1 * y_dir) %>%
  mutate(UMAP_1 = UMAP_a, UMAP_2 = UMAP_b)

tabulamuris_predict_labels <-scEiaD_2020_v01 %>% tbl('tabulamuris_predict_labels') %>% collect %>%
  mutate(UMAP_a = UMAP_2 * x_dir,
         UMAP_b = UMAP_1 * y_dir) %>%
  mutate(UMAP_1 = UMAP_a, UMAP_2 = UMAP_b)
celltype_predict_labels <-scEiaD_2020_v01 %>% tbl('celltype_predict_labels')  %>% collect %>%
  mutate(UMAP_a = UMAP_2 * x_dir,
         UMAP_b = UMAP_1 * y_dir) %>%
  mutate(UMAP_1 = UMAP_a, UMAP_2 = UMAP_b)
celltype_labels <- scEiaD_2020_v01 %>% tbl('celltype_labels') %>% collect %>%
  mutate(UMAP_a = UMAP_2 * x_dir,
         UMAP_b = UMAP_1 * y_dir) %>%
  mutate(UMAP_1 = UMAP_a, UMAP_2 = UMAP_b)
cluster_labels <-scEiaD_2020_v01 %>% tbl('cluster_labels') %>% collect %>%
  mutate(UMAP_a = UMAP_2 * x_dir,
         UMAP_b = UMAP_1 * y_dir) %>%
  mutate(UMAP_1 = UMAP_a, UMAP_2 = UMAP_b)
mf <- meta_filter %>% sample_frac(0.2)

# generate color_mappings
categorical_columns <- c("Phase","batch","study_accession","library_layout","organism","Platform",
                         "Covariate","CellType","CellType_predict","TabulaMurisCellType","TabulaMurisCellType_predict",
                         "GSE","Summary","Design","Citation","PMID","Stage","cluster",
                         "Doublet","TechType", "SubCellType", 'subcluster', 'Age', "retina_region",
                         'Tissue','Organ', 'Source','sample_accession')
#"SubCellType" and subcluster are problems
meta_filter <- meta_filter %>% mutate(Age = as.character(Age), SubCellType = tidyr::replace_na(SubCellType, 'None'),
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

#this is to avoid a color collision within the same cluster
# sub_cluster_df <- meta_filter %>% select(cluster, subcluster) %>% distinct
# sub_cluster_map <- lapply(sub_cluster_df$cluster, function(x) sub_cluster_df %>%
#                             filter(cluster == x) %>%
#                             map_color('subcluster',.) ) %>% bind_rows

# this could work, but the Celltype  > SubCellType mapping is not 1:1
# sub_celltype_df <- meta_filter %>% select(CellType_predict, SubCellType) %>% distinct
# sub_celltype_map <- lapply(sub_celltype_df$CellType_predict, function(x) sub_celltype_df %>%
#                              filter(CellType_predict == x) %>%
#                              map_color('SubCellType',.) ) %>% bind_rows %>%
#   mutate(color = replace(color, value == 'None', '#D3D3D3'))


cat_to_color_df <- lapply(categorical_columns, function(col) map_color(col, meta_filter)) %>% bind_rows()
# %>%
#   bind_rows(sub_cluster_map)



# # attach colors to cell types
# cell_types <- meta_filter %>%
#   pull(CellType_predict) %>% unique() %>% sort()
# type_val <- setNames(pals::alphabet(n = cell_types %>% length()), cell_types)
# type_col <- scale_colour_manual(values = type_val)
# type_fill <- scale_fill_manual(values = type_val)
cat(file=stderr(), 'Data loaded in ')
cat(file=stderr(), Sys.time() - time)
cat(file=stderr(), ' seconds.\n')

shinyOptions(cache = cachem::cache_disk(file.path(dirname(tempdir()), "plae-cache")))

# site begins! ---------
shinyServer(function(input, output, session) {
  #bootstraplib::bs_themer()
  observe({
    # URL scanning to jump directly to tab/section
    query <- parseQueryString(session$clientData$url_search)
    if(!is.null(query$url)) {
      cat(query$url)
      url <- strsplit(query$url,"\"")[[1]][2]
      cat(url)
      updateTabsetPanel(session, 'nav', query$url)
    }

    # server help queries ------
    # gene plot updateSelectizeInput -------
    if (is.null(query[['Gene']])){
      updateSelectizeInput(session, 'Gene',
                           choices = scEiaD_2020_v01 %>% tbl('genes') %>% collect() %>% pull(1),
                           options = list(placeholder = 'Type to search'),
                           selected = 'CRX (ENSG00000105392)',
                           server = TRUE)
    }
    # gene plot category filtering ----
    if (is.null(query[['gene_filter_cat']])){
      updateSelectizeInput(session, 'gene_filter_cat',
                           label = 'Scatter Filter Category: ',
                           choices = meta_filter %>%
                             dplyr::select(nCount_RNA:doublet_score_scran) %>% colnames() %>% sort(),
                           selected = '',
                           server = TRUE)
    }
    observeEvent(input$gene_filter_cat, {
      if (input$gene_filter_cat == ''){
        choice = ''
      } else {
        choice = meta_filter[,input$gene_filter_cat] %>% t() %>% c() %>% unique() %>% sort()
      }
      output$gene_filter_on_dynamicUI <- renderUI({
        if (class(choice) == 'character'){
          selectizeInput('gene_filter_on', strong('Gene Select: '),
                         choices = choice, selected = NULL, multiple = TRUE)
        } else {
          shinyWidgets::setSliderColor(c("#3399ff"), c(1))
          sliderInput("gene_filter_on", label = strong("Gene Filter Range: "), min = min(choice),
                      max = max(choice), value = c(min(choice), max(choice)))
        }
      })
    })

    # meta plot updateSelectizeInput ------
    if (is.null(query[['meta_column']])){
      updateSelectizeInput(session, 'meta_column',
                           label = 'Meta Color:',
                           choices = meta_filter %>%
                             dplyr::select(nCount_RNA:doublet_score_scran) %>% colnames() %>% sort(),
                           options = list(placeholder = 'Type to search'),
                           selected = 'CellType_predict',
                           server = TRUE)
    }
    # meta plot category filtering ----
    if (is.null(query[['meta_filter_cat']])){
      updateSelectizeInput(session, 'meta_filter_cat',
                           label = 'Meta Filter Category:',
                           choices = meta_filter %>%
                             dplyr::select(nCount_RNA:doublet_score_scran) %>% colnames() %>% sort(),
                           selected = '',
                           server = TRUE)
    }

    observeEvent(input$meta_filter_cat, {
      if (input$meta_filter_cat == ''){
        choice = ''
      } else {
        choice = meta_filter[,input$meta_filter_cat] %>% t() %>% c() %>% unique() %>% sort()
      }
      output$meta_filter_on_dynamicUI <- renderUI({
        if (class(choice) == 'character'){
          selectizeInput('meta_filter_on', strong('Meta Select: '),
                         choices = choice, selected = NULL, multiple = TRUE)
        } else {
          shinyWidgets::setSliderColor(c("#3399ff"), c(1))
          sliderInput("meta_filter_on", label = strong("Meta Filter Range: "), min = min(choice),
                      max = max(choice), value = c(min(choice), max(choice)))
        }
      })

    })

    # dotplot updateSelectizeInput ----
    if (is.null(query[['dotplot_Gene']])){
      updateSelectizeInput(session, 'dotplot_Gene',
                           label = 'Genes: ',
                           choices = scEiaD_2020_v01 %>% tbl('genes') %>% collect() %>% pull(1),
                           options = list(placeholder = 'Type to search'),
                           selected = c('ARR3 (ENSG00000120500)','AIF1L (ENSG00000126878)','GAD1 (ENSG00000128683)','POU4F2 (ENSG00000151615)','WIF1 (ENSG00000156076)','RHO (ENSG00000163914)','ONECUT1 (ENSG00000169856)','GRIK1 (ENSG00000171189)','AIF1 (ENSG00000204472)'),
                           server = TRUE)
    }
    if (is.null(query[['dotplot_groups']])){
      updateSelectizeInput(session, 'dotplot_groups',
                           label = 'Group by (two max): ',
                           choices = meta_filter %>%
                             select(-Barcode) %>%
                             select_if(purrr::negate(is.numeric)) %>%
                             colnames() %>% sort(),
                           options = list(placeholder = 'Type to search',
                                          maxItems = 2),
                           selected = c('CellType_predict'),
                           server = TRUE)
    }
    if (is.null(query[['dotplot_filter_cat']])){
      updateSelectizeInput(session, 'dotplot_filter_cat',
                           label = 'Filter Category: ',
                           choices = meta_filter %>%
                             dplyr::select(nCount_RNA:doublet_score_scran) %>% colnames() %>% sort(),
                           selected = '',
                           server = TRUE)
    }
    observeEvent(input$dotplot_filter_cat, {
      if (input$dotplot_filter_cat == ''){
        choice = ''
      } else {
        choice = meta_filter[,input$dotplot_filter_cat] %>% pull(1) %>% unique() %>% sort()
      }
      updateSelectizeInput(session, 'dotplot_filter_on',
                           label = 'Select: ',
                           choices = choice,
                           server = TRUE)
    })

    # insitu updateSelectizeInput ----
    if (is.null(query[['insitu_Gene']])){
      updateSelectizeInput(session, 'insitu_Gene',
                           choices = scEiaD_2020_v01 %>% tbl('genes') %>% as_tibble() %>% pull(1),
                           options = list(placeholder = 'Type to search'),
                           selected = c('RHO (ENSG00000163914)'),
                           server = TRUE)
    }

    if (is.null(query[['insitu_filter_cat']])){
      updateSelectizeInput(session, 'insitu_filter_cat',
                           label = 'Filter category: ',
                           choices = scEiaD_2020_v01 %>% tbl('grouped_stats') %>%
                             select(-Gene, -cell_ct, -cell_exp_ct, -counts) %>% colnames() %>% sort(),
                           selected = '',
                           server = TRUE)
    }
    observeEvent(input$insitu_filter_cat, {
      if (input$insitu_filter_cat == ''){
        choice = ''
      } else {
        choice = meta_filter[,input$insitu_filter_cat] %>% pull(1) %>% unique() %>% sort()
      }
      updateSelectizeInput(session, 'insitu_filter_on',
                           label = "Select: ",
                           choices = choice,
                           server = TRUE)
    })

    #
    if (is.null(query[['grouping_features']])){
      updateSelectizeInput(session, 'grouping_features',
                           choices = scEiaD_2020_v01 %>% tbl('grouped_stats') %>%
                             select(-Gene, -cell_ct, -cell_exp_ct, -counts) %>% colnames() %>% sort(),
                           options = list(placeholder = 'Type to search'),
                           selected = c('CellType_predict', 'organism'),
                           server = TRUE)
    }

    if (is.null(query[['meta_groupings']])){
      updateSelectizeInput(session, 'meta_groupings',
                           choices = meta_filter %>%
                             select(-Barcode, -UMAP_1, -UMAP_2, -nCount_RNA, -nFeature_RNA, -percent_mt, -UMAP_a, -UMAP_b) %>%
                             colnames() %>% sort(),
                           options = list(placeholder = 'Type to search'),
                           selected = c('CellType_predict', 'organism'),
                           server = TRUE)
    }
    # exp_plot plot updateSelect -----
    if (is.null(query[['exp_plot_genes']])){
      # The java stuffs allw for comma, space, and semicolon
      # deliminted input
      # https://github.com/rstudio/shiny/issues/1663
      updateSelectizeInput(session, 'exp_plot_genes',
                           choices = scEiaD_2020_v01 %>% tbl('genes') %>% collect() %>% pull(1),
                           #options = list(placeholder = 'Type to search'),
                           selected = c('PAX6 (ENSG00000007372)','CRX (ENSG00000105392)','NRL (ENSG00000129535)','POU4F2 (ENSG00000151615)'),
                           options = list(
                             placeholder = 'Type to search',
                             splitOn = I("(function() { return /[, ;]/; })()"),
                             create = I("function(input, callback){
                                          return {
                                            value: input,
                                            text: input
                                          };
                                        }")
                           ),
                           server = TRUE)
    }
    if (is.null(query[['exp_filter_cat']])){
      updateSelectizeInput(session, 'exp_filter_cat',
                           choices = scEiaD_2020_v01 %>% tbl('grouped_stats') %>%
                             select(-Gene, -cell_ct, -cell_exp_ct, -counts) %>% colnames() %>% sort(),
                           selected = 'CellType_predict',
                           server = TRUE)
    }
    observeEvent(input$exp_filter_cat, {
      if (input$exp_filter_cat == ''){
        choice = ''
      } else {
        choice = meta_filter[,input$exp_filter_cat] %>% t() %>% c() %>% unique() %>% sort()
      }

      updateSelectizeInput(session, 'exp_filter_on',
                           choices = choice,
                           selected = c('Rod','Retinal Ganglion Cell', 'Horizontal Cell'),
                           server = TRUE)
    })
    if (is.null(query[['exp_plot_groups']])){
      updateSelectizeInput(session, 'exp_plot_groups',
                           choices = scEiaD_2020_v01 %>% tbl('grouped_stats') %>%
                             select(-Gene, -cell_ct, -cell_exp_ct, -counts) %>% colnames() %>% sort(),
                           selected = c('study_accession'),
                           server = TRUE)
    }



    # temporal plot updateSelect -----
    if (is.null(query[['temporal_gene']])){
      updateSelectizeInput(session, 'temporal_gene',
                           choices = scEiaD_2020_v01 %>% tbl('genes') %>% collect() %>% pull(1),
                           options = list(placeholder = 'Type to search'),
                           selected = c('PAX6 (ENSG00000007372)','POU4F2 (ENSG00000151615)'),
                           server = TRUE)
    }


    # facet plot updateSelect --------
    if (is.null(query[['facet']])){
      updateSelectizeInput(session, 'facet',
                           choices = meta_filter %>%
                             dplyr::select_if(is.character) %>% colnames() %>% sort(),
                           options = list(placeholder = 'Type to search'),
                           selected = 'organism',
                           server = TRUE)
    }
    #facet plot updateSelectizeInput
    if (is.null(query[['facet_filter_cat']])){
      updateSelectizeInput(session, 'facet_filter_cat',
                           choices = scEiaD_2020_v01 %>% tbl('grouped_stats') %>%
                             select(-Gene, -cell_ct, -cell_exp_ct, -counts) %>% colnames() %>% sort(),
                           selected = 'CellType_predict',
                           server = TRUE)
    }
    observeEvent(input$facet_filter_cat, {
      if (input$exp_filter_cat == ''){
        choice = ''
      } else {
        choice = meta_filter[,input$facet_filter_cat] %>% t() %>% c() %>% unique() %>% sort()
      }

      updateSelectizeInput(session, 'facet_filter_on',
                           choices = choice,
                           selected = c('Cone','Retinal Ganglion Cell', 'Horizontal Cell'),
                           server = TRUE)
    })

    if (is.null(query[['facet_color']])){
      updateSelectizeInput(session, 'facet_color',
                           choices = meta_filter %>%
                             dplyr::select_if(is.character) %>% colnames() %>% sort(),
                           options = list(placeholder = 'Type to search'),
                           selected = 'CellType_predict',
                           server = TRUE)
    }
    # diff table updateSelect ------
    if (is.null(query[['diff_gene']])){
      updateSelectizeInput(session, 'diff_gene',
                           choices = scEiaD_2020_v01 %>% tbl('wilcox_diff_AUC_genes') %>% collect() %>% pull(1),
                           options = list(placeholder = 'Type to search'),
                           selected = 'RHO (ENSG00000163914)',
                           server = TRUE)
    }
    if (is.null(query[['diff_base']])){
      group = input$search_by
      choices = scEiaD_2020_v01 %>%
        tbl('wilcox_diff_AUC_sets') %>%
        filter(Group == group) %>%
        collect() %>% filter(!grepl('Doubl', Base)) %>%
        pull(Base)
      updateSelectizeInput(session, 'diff_base',
                           choices = choices,
                           options = list(placeholder = 'Type to search'),
                           server = TRUE)
    }
    updateSelectizeInput(session, 'diff_against',
                         choices = scEiaD_2020_v01 %>%
                           tbl('wilcox_diff_AUC_sets') %>%
                           filter(Group == group) %>%
                           collect() %>% filter(!grepl('Doubl', Base)) %>%
                           pull(Base),
                         selected = '',
                         options = list(placeholder = 'Optional Filtering'),
                         server = TRUE)

    # Meta Plot modal ----------
    observeEvent(input$BUTTON_show_meta_legend, {
      # Show a modal when the button is pressed
      showModal(shinyjqui::draggableModalDialog(size = 'l', title = 'Click to Drag',
                                                plotOutput('meta_plot_legend') %>% shinycssloaders::withSpinner(type = 3, size = 0.5, color = "#3399ff", color.background = 'white'),
                                                easyClose = TRUE))
    })
    # HELP button descriptions ----------
    ## table
    dt_help_html <- "<p>The data tables in this app are reactive and searchable. The first (top right) field searches across all columns. If you need to search on individual columns, below each column is a search field that lets you search by typing. If the column is numeric, a slider appears allowing you to set a range of values to filter to. There are two ways to export information from a table:
                                                <ul>
                                                <li> Click \"copy\" and paste into Excel </li>
                                                <li> Click \"CSV\" and a comma separated file will download </li> </ul></p>
                                                     <p>The tables tend to be huge in size, so only 10 (by default) rows are shown at a time. You can click on the \"Show 10 rows\" button and select up to 100 rows to display</p>"
    observeEvent(input$data_table_help1, {
      # Show a modal when the button is pressed
      showModal(shinyjqui::draggableModalDialog(size = 'm',
                                                title = "Data Tables",
                                                HTML(dt_help_html),
                                                easyClose = TRUE))
    })
    observeEvent(input$data_table_help2, {
      # Show a modal when the button is pressed
      showModal(shinyjqui::draggableModalDialog(size = 'm',
                                                title = "Data Tables",
                                                HTML(dt_help_html),
                                                easyClose = TRUE))
    })
    observeEvent(input$data_table_help3, {
      # Show a modal when the button is pressed
      showModal(shinyjqui::draggableModalDialog(size = 'm',
                                                title = "Data Tables",
                                                HTML(dt_help_html),
                                                easyClose = TRUE))
    })
    ## umap
    observeEvent(input$umap_table_help, {
      # Show a modal when the button is pressed
      showModal(shinyjqui::draggableModalDialog(size = 'm',
                                                title = "UMAP - Tables",
                                                HTML("<p>The UMAP is a 2D Projection of a higher dimensional space which tries to bring together
                            closely related elements while maintaining the overall structure. The higher dimensional space
                            is built from the gene expression patterns of each cell. The left panel shows the expression pattern
                            of a gene in the UMAP space. You can remove cells with low or high expression of your gene of choice by using the slider in \"Filter Gene Expression\" to select a range. The right panel shows metadata associated with each cell, including
                            predicted cell type.</p>

                            Tips:
                            <ul>
                              <li>You can click and drag to set a box, then double click to zoom in!</li>
                              <li>Double click again to zoom back</li>
                            </ul>"),
                                                easyClose = TRUE))
    })
    observeEvent(input$exp_plot_help, {
      showModal(shinyjqui::draggableModalDialog(size = 's',
                                                title = "Expression Plot",
                                                HTML("<p>This highly flexible scatter plot view allows you to plot gene expression by Cell Type
                                                (published), Cell Type (our inferred cell labels), or cluster (unsupervised grouping
                                                of the single cell transcriptomes). You can control how the data is displayed by
                                                selecting how the data points are colored and you have advanced filtering
                                                functionality that lets you select what fields (e.g. Citation) to filter on.</p>

                            There are two fundamental ways to view the data:
                            <ul>
                              <li>Expression, which is log2 scaled counts</li>
                              <li>% Cells Detected, which is the proportion of cells that have any detectable transcript</li>
                            </ul>"),
                                                easyClose = TRUE))
    })
    observeEvent(input$insitu_help, {
      showModal(shinyjqui::draggableModalDialog(size = 's',
                                                title = "In Situ Projection",
                                                HTML("<p>As developmental biologists may be more experienced in viewing
                                                stained cross-sections of the retina, we have created this visualization
                                                which colors each of the major cell type (e.g. Rods, Cones) by the intensity
                                                of the expression of user selected gene. Like elsewhere in PLAE, you can
                                                filter to only show data by flexible criteria (e.g. only show gene / cell
                                                expression infrmation from mouse)</p>"),
                                                easyClose = TRUE))
    })
    observeEvent(input$facet_umap_help, {
      showModal(shinyjqui::draggableModalDialog(size = 's',
                                                title = "Facet UMAP",
                                                HTML("<p>Because some of the information you are interested in
                                                     may be overlapping, it sometimes is useful to be able to split
                                                     the UMAP visualization into separate plots by some field (e.g.
                                                     organism) of interest.</p>"),
                                                easyClose = TRUE))
    })
    observeEvent(input$dotplot_help, {
      showModal(shinyjqui::draggableModalDialog(size = 's',
                                                title = "Dotplot",
                                                HTML("<p>The dotplot visualization is highly space efficent as it displays
                                                     both amount of expression (by color intensity) and percent cells
                                                     that detect the transcript (by size of dot) with a user selected
                                                     combination of grouping features (e.g. organism and CellType).</p>"),
                                                easyClose = TRUE))
    })
    observeEvent(input$diff_testing_help, {
      showModal(shinyjqui::draggableModalDialog(size = 'l',
                                                title = "Differential Testing",
                                                HTML("<p>We have pre-computed 3 different differential expression tests.

                              <li>CellType, which are based on published cell type assignments</li>
                              <li>CellType (predict), which uses ML to project CellType labels onto (nearly) all of the cells</li>
                              <li>Cluster (droplet), which are created from leiden method on the scVI correct lower dimension space.</li>
                            </ul></p>"),
                                                easyClose = TRUE))
    })
    observeEvent(input$diff_testing_help2, {
      showModal(shinyjqui::draggableModalDialog(size = 'l',
                                                title = "PseudoBulk Design",
                                                HTML("<p>What is pseudobulk? Most scRNA-seq diff testing is done with complicated
                                                     and computationally expensive tests designed to maximize signal in relatively
                                                     homogenous data. As the scEiaD dataset includes a great many biological
                                                     replicates we have chosen to maximize the power of the replicates by summing
                                                     the counts data into groups (e.g. sum all CRX counts in labelled Rods for Clark et al.,
                                                     Lu et al., etc). This makes the data bulk (traditional) RNA-seq like in its
                                                     statistical properties, which allows us to use more established tools (like edgeR)
                                                     for the differential testing.</p>

                                                     <p>Model design for [ ] against Remaining and Pairwise [ ] against [ ]</p>
                                                     <ul>
                                                     <li>design <- model.matrix(~0+group+org_covariate) where group is CellType/Cluster and
                                                     org_covariate is Human/Mouse/Macaque</li>
                                                     </ul>
                                                     <p>Model design for Organism specific test within [ ]</p>
                                                     <ul>
                                                     <li>design <- model.matrix(~0+group) where group is Human/Mouse/Macaque and
                                                     the test is limited by contrasts to a specific CellType or Cluster</li>
                                                     </ul>
                                                     <p>The edgeR glmQLFTest is used to fit the linear model for the differential testing.</p>"),
                                                easyClose = TRUE))
    })

    ## exp plot
    observeEvent(input$exp_plot_help, {
      # Show a modal when the button is pressed
      shinyalert("Help", 'help_text', type = "info")
    })
    ## insitu
    observeEvent(input$insitu_help, {
      # Show a modal when the button is pressed
      shinyalert("Help", 'help_text', type = "info")
    })
    ## facet umap
    observeEvent(input$facet_umap_help, {
      # Show a modal when the button is pressed
      shinyalert("Help", 'help_text', type = "info")
    })
    ## temporal plot
    observeEvent(input$temporal_plot_help, {
      # Show a modal when the button is pressed
      shinyalert("Help", 'help_text', type = "info")
    })
    ## dotplot
    observeEvent(input$dotplot_help, {
      # Show a modal when the button is pressed
      shinyalert("Help", 'help_text', type = "info")
    })
    ## diff testing
    observeEvent(input$diff_testing_help, {
      # Show a modal when the button is pressed
      shinyalert("Help", 'help_text', type = "info")
    })


    # BREAK -------
    # gene scatter plot ------------

    x_range =  c(meta_filter$UMAP_1 %>% min(), meta_filter$UMAP_1 %>% max())
    y_range = c(meta_filter$UMAP_2 %>% min(), meta_filter$UMAP_2 %>% max())

    gene_scatter_ranges <- reactiveValues(x = x_range,
                                          y = y_range)
    source('make_gene_scatter_umap_plot.R')
    gene_scatter_plot <- eventReactive(input$BUTTON_draw_scatter, {
      make_gene_scatter_umap_plot(input,
                                  scEiaD_2020_v01,
                                  mf,
                                  meta_filter,
                                  celltype_predict_labels,
                                  celltype_labels,
                                  tabulamuris_predict_labels,
                                  cluster_labels)
    })

    observeEvent(input$gene_scatter_plot_dblclick, {
      brush <- input$gene_scatter_plot_brush
      if (!is.null(brush)) {
        gene_scatter_ranges$x <- c(brush$xmin, brush$xmax)
        gene_scatter_ranges$y <- c(brush$ymin, brush$ymax)

      } else {
        gene_scatter_ranges$x <- x_range
        gene_scatter_ranges$y <- y_range
      }
    })
    output$gene_scatter_plot <- renderPlot({
      gene_scatter_plot() + coord_cartesian(xlim = gene_scatter_ranges$x, ylim = gene_scatter_ranges$y)
    }) %>%
      bindCache(list(input$Gene,
                     input$pt_size_gene,
                     input$gene_scatter_slider,
                     input$gene_label_toggle,
                     input$gene_filter_cat,
                     input$gene_filter_cat_on,
                     gene_scatter_ranges$x,
                     gene_scatter_ranges$y)) %>%
      bindEvent(input$BUTTON_draw_scatter,
                input$gene_scatter_plot_dblclick)

    # gene scatter plot download ------
    output$BUTTON_download_scatter <- downloadHandler(
      filename = function() { ('plae_gene_scatter.png') },
      content = function(file) {
        ggsave(file, plot = gene_scatter_plot() + coord_cartesian(xlim = gene_scatter_ranges$x, ylim = gene_scatter_ranges$y), device = "png")
      }
    )

    # metadata plot --------------
    source('make_meta_scatter_umap_plot.R')
    meta_plot <- eventReactive(input$BUTTON_draw_meta, {
      make_meta_scatter_umap_plot(input, mf, meta_filter,
                                  celltype_predict_labels,
                                  celltype_labels,
                                  tabulamuris_predict_labels,
                                  cluster_labels,
                                  cat_to_color_df
      )
    })

    x_range =  c(meta_filter$UMAP_1 %>% min(), meta_filter$UMAP_1 %>% max())
    y_range = c(meta_filter$UMAP_2 %>% min(), meta_filter$UMAP_2 %>% max())

    meta_ranges <- reactiveValues(x = x_range,
                                  y = y_range)
    observeEvent(input$meta_plot_dblclick, {
      brush <- input$meta_plot_brush
      if (!is.null(brush)) {
        meta_ranges$x <- c(brush$xmin, brush$xmax)
        meta_ranges$y <- c(brush$ymin, brush$ymax)
      } else {
        meta_ranges$x <- x_range
        meta_ranges$y <- y_range
      }
    })

    output$meta_plot <- renderPlot({
      plot_data <- meta_plot()
      if (plot_data$col_size < 10) {
        plot_data$plot + coord_cartesian(xlim = meta_ranges$x, ylim = meta_ranges$y)
      } else {
        plot_data$plot + coord_cartesian(xlim = meta_ranges$x, ylim = meta_ranges$y) +
          theme(legend.position = 'none')
      }
    }) %>%
      bindCache(list(input$meta_column,
                     input$pt_size_meta,
                     input$label_toggle,
                     input$meta_column_transform,
                     input$meta_filter_cat,
                     input$meta_filter_on,
                     meta_ranges$x,
                     meta_ranges$y)) %>%
      bindEvent(input$BUTTON_draw_meta,
                input$meta_plot_dblclick)

    output$BUTTON_download_meta <- downloadHandler(
      filename = function() { ('plae_meta.png') },
      content = function(file) {
        if (meta_plot()$col_size < 10) {
          ggsave(file, plot = meta_plot()$plot + coord_cartesian(xlim = meta_ranges$x, ylim = meta_ranges$y), device = "png")
        } else {
          ggsave(file, plot = meta_plot()$plot + coord_cartesian(xlim = meta_ranges$x, ylim = meta_ranges$y) +
                   theme(legend.position = 'none'), device = "png")
        }
      }
    )
    output$meta_plot_legend <- renderPlot({
      plot <- meta_plot()$plot
      legend <- cowplot::get_legend(plot)
      plot_grid(legend)
    })


    # gene cluster table  --------
    gene_cluster_stats_maker <- eventReactive(input$BUTTON_make_gene_table, {
      grouping_features <- input$grouping_features
      gene <- input$Gene
      table <- scEiaD_2020_v01 %>% tbl('grouped_stats') %>%
        filter(Gene == gene) %>%
        group_by_at(vars(one_of(c('Gene', grouping_features)))) %>%
        summarise(counts = sum(counts * cell_exp_ct) / sum(cell_exp_ct),
                  cell_exp_ct = sum(cell_exp_ct, na.rm = TRUE)) %>%
        collect() %>%
        tidyr::drop_na() %>%
        full_join(., meta_filter %>%
                    group_by_at(vars(one_of(grouping_features))) %>%
                    summarise(Count = n())) %>%
        mutate(cell_exp_ct = ifelse(is.na(cell_exp_ct), 0, cell_exp_ct)) %>%
        mutate(`%` = round((cell_exp_ct / Count) * 100, 2),
               Expression = round(counts * (`%` / 100), 2)) %>%
        select_at(vars(one_of(c('Gene', grouping_features, 'cell_exp_ct', 'Count', '%', 'Expression')))) %>%
        arrange(-Expression) %>%
        rename(`Cells # Detected` = cell_exp_ct,
               `Total Cells` = Count,
               `log2(counts+1)` = Expression) %>%
        ungroup() %>%
        select(-Gene)

      table %>% DT::datatable(extensions = 'Buttons',
                              filter = list(position = 'bottom', clear = TRUE, plain = TRUE),
                              options = list(pageLength = 10, scrollX = T, searchHighlight = TRUE, dom = 'frtBip', buttons = c('pageLength','copy', 'csv')))
    })
    output$gene_cluster_stats <- DT::renderDataTable({ gene_cluster_stats_maker()})


    # metadata table -----
    metadata_stats <- eventReactive(input$BUTTON_make_meta_table, {
      grouping_features <- input$meta_groupings
      table <- meta_filter %>%
        group_by_at(vars(one_of(grouping_features))) %>%
        summarise(Count = n()) %>%
        select_at(vars(one_of(c(grouping_features, 'Count')))) %>%
        arrange(-Count) %>%
        rename(`Total Cells` = Count)

      table %>% DT::datatable(extensions = 'Buttons',
                              filter = list(position = 'bottom', clear = TRUE, plain = TRUE),
                              options = list(pageLength = 10, scrollX = T, searchHighlight = TRUE, dom = 'frtBip', buttons = c('pageLength','copy', 'csv')))
    })
    output$metadata_stats <- DT::renderDataTable({ metadata_stats()})

    # facet plot -----------
    source('make_facet_plot.R')

    facet_plot <- eventReactive(input$BUTTON_draw_filter, {make_facet_plot(input,  meta_filter)})

    output$facet_plot <- renderPlot({
      facet_plot()
    }, height = eventReactive(input$BUTTON_draw_filter, {input$facet_height %>% as.numeric()}))
    output$BUTTON_download_facet <- downloadHandler(
      filename = function() { ('plae_facet.png') },
      content = function(file) {
        ggsave(file, plot = facet_plot(),
               height = as.numeric(input$facet_height) / 50,
               width = 15,
               device = "png")
      }
    )
    ## exp_plot -----------
    source('make_exp_plot.R')
    exp_plot <- eventReactive(input$BUTTON_draw_exp_plot, {
      make_exp_plot(input, scEiaD_2020_v01, meta_filter)
    })

    output$exp_plot <- renderPlot({
      exp_plot()
    }, height = eventReactive(input$BUTTON_draw_exp_plot, {as.numeric(input$exp_plot_height)}))

    output$BUTTON_download_exp <- downloadHandler(
      filename = function() { ('plae_exp.png') },
      content = function(file) {
        ggsave(file, plot = make_exp_plot(input, scEiaD_2020_v01, meta_filter),
               height = as.numeric(input$exp_plot_height) / 50,
               width = 15,
               device = "png")
      }
    )

    # ## temporal plot -----------
    # source('make_temporal_plot.R')
    # temporal_plot <- eventReactive(input$BUTTON_draw_temporal, {
    #   make_temporal_plot(input, scEiaD_2020_v01, meta_filter)
    # })
    # output$temporal_plot <- renderPlot({
    #   temporal_plot()
    # }, height = as.numeric(input$temporal_plot_height ) )

    ## dotplot ---------
    source('make_dotplot.R')
    draw_dotplot <- eventReactive(input$BUTTON_draw_dotplot,
                                  {make_dotplot(input, scEiaD_2020_v01, meta_filter,cat_to_color_df)}
    )
    output$dotplot <- renderPlot({
      draw_dotplot()
    }, height = eventReactive(input$BUTTON_draw_dotplot, {input$dotplot_height %>% as.numeric()}))
  })
  output$BUTTON_download_dotplot <- downloadHandler(
    filename = function() { ('plae_dotplot.png') },
    content = function(file) {
      ggsave(file, plot = make_dotplot(input, scEiaD_2020_v01, meta_filter,cat_to_color_df),
             height =
               input$dotplot_height %>% as.numeric() / 50,
             width = 12,
             device = "png")
    }
  )

  # in situ ----
  # Functions used to generate in situ plots
  source('make_in_situ_plot.R')
  # Event reactives to produce image and table
  ## Reactive that generates data table
  insitu_table_maker <- eventReactive(input$BUTTON_draw_insitu, {

    full_table <- get_insitu_table(input, scEiaD_2020_v01, meta_filter)
    full_table <- full_table %>%
      rename(`Cells # Detected` = cell_exp_ct,
             `Total Cells` = Count,
             `log2(counts+1)` = Expression) %>%
      tidyr::drop_na()
    full_table %>% DT::datatable(extensions = 'Buttons',
                                 filter = list(position = 'bottom', clear = TRUE, plain = TRUE),
                                 options = list(pageLength = 10, searchHighlight = TRUE, dom = 'frtBip', buttons = c('pageLength','copy', 'csv')))
  })

  ## Reactive that generates in situ image
  make_pic <- eventReactive( input$BUTTON_draw_insitu, {
    make_insitu_plot(input, scEiaD_2020_v01, meta_filter)
  })

  # Output in situ
  ## Output table
  output$insitu_gene_stats <- DT::renderDataTable({ insitu_table_maker()})

  ## Output image
  output$insitu_img <- renderImage({
    make_pic()
  }, deleteFile = TRUE)

  ## diff table code --------
  ### diff table auc -----
  output$make_diff_table_auc <- DT::renderDataTable(server = TRUE, {
    gene <- input$diff_gene
    if (input$search_by == 'Gene'){
      # out_gene <- scEiaD_2020_v01 %>% tbl('wilcox_diff_testing') %>%
      #   filter(Gene %in% gene) #%>%
      out_auc <- scEiaD_2020_v01 %>% tbl('wilcox_diff_AUC') %>%
        filter(Gene %in% gene, Base != 'Doublet') %>%
        head(2000) %>%
        collect() %>%
        select(Base, `Tested Against`, AUC, logFC) %>%
        mutate(
          AUC = format(AUC, digits = 3) %>% as.numeric(.),
          logFC = format(logFC, digits = 3) %>% as.numeric(.))
    }
    else {
      req(input$diff_base)
      diff_base <- input$diff_base
      filter_term <- input$search_by
      if (input$diff_against == ''){
        out_auc <- scEiaD_2020_v01 %>% tbl('wilcox_diff_AUC') %>%
          filter(Base == diff_base) %>%
          head(2000) %>%
          filter(Base == diff_base,
                 Group == filter_term)
        # out_gene <- scEiaD_2020_v01 %>% tbl('wilcox_diff_testing') %>%
        #   filter(Base == diff_base) %>%
        #   filter(Group == filter_term)
      } else {
        against <- input$diff_against
        out_auc <- scEiaD_2020_v01 %>% tbl('wilcox_diff_AUC') %>%
          filter(Base == diff_base,
                 Group == filter_term,
                 `Tested Against` == against) %>%
          head(2000)
      }
      out_auc <- out_auc %>%
        select(Gene, Base, `Tested Against`, AUC, logFC) %>%
        collect() %>%
        mutate(
          AUC = format(AUC, digits = 3),
          AUC = as.numeric(AUC),
          logFC = format(logFC, digits = 3) %>% as.numeric(.))
    }
    out_auc %>%
      DT::datatable(extensions = 'Buttons',
                    caption = htmltools::tags$caption( style = 'caption-side: top; text-align: left; color:black; font-size:200% ;','Table 2: Pairwise AUC Diff Testing'),
                    filter = list(position = 'bottom', clear = TRUE, plain = TRUE),
                    options = list(pageLength = 10, scrollX = TRUE,
                                   dom = 'rtBip', buttons = c('pageLength','copy'))) %>%
      DT::formatStyle(columns = c(8), width='250px')
  })
  ### diff table gene - base level -----
  output$make_diff_table <- DT::renderDataTable(server = TRUE, {
    gene <- input$diff_gene
    if (input$search_by == 'Gene'){
      table_name = 'Table 1: Group - Gene - Base Diff Testing'
      out_gene <- scEiaD_2020_v01 %>% tbl('wilcox_diff_testing') %>%
        filter(Gene %in% gene) %>%
        select(Group, Base,  p.value, FDR,  mean_auc, mean_logFC) %>%
        collect() %>%
        mutate(Group = as.factor(Group)) %>%
        mutate(PValue = format(as.numeric(p.value), digits = 3) %>% as.numeric(),
               FDR = format(as.numeric(FDR), digits = 3) %>% as.numeric(),
               `Mean AUC` = format(mean_auc, digits = 3) %>% as.numeric(),
               `Mean logFC` = format(mean_logFC, digits = 3) %>% as.numeric()) %>%
        select(-p.value, -mean_auc, -mean_logFC)
    } else {
      table_name = 'Table 1: Group - Gene - Base Diff Testing'
      req(input$diff_base)
      diff_base <- input$diff_base
      filter_term <- input$search_by
      out_gene <- scEiaD_2020_v01 %>% tbl('wilcox_diff_testing') %>%
        filter(Base == diff_base,
               Group == filter_term) %>%
        select(Gene, Base,  p.value, FDR, mean_auc, mean_logFC) %>%
        head(2000) %>%
        collect() %>%
        mutate(PValue = format(as.numeric(p.value), digits = 3) %>% as.numeric(),
               FDR = format(as.numeric(FDR), digits = 3) %>% as.numeric(),
               `Mean AUC` = format(mean_auc, digits = 3) %>% as.numeric(),
               `Mean logFC` = format(mean_logFC, digits = 3) %>% as.numeric()) %>%
        select(-p.value, -mean_auc, -mean_logFC)
    }
    out_gene %>%
      DT::datatable(extensions = 'Buttons',
                    caption = htmltools::tags$caption( style = 'caption-side: top; text-align: left; color:black; font-size:200% ;',table_name),
                    filter = list(position = 'bottom', clear = TRUE, plain = TRUE),
                    options = list(pageLength = 10, scrollX = TRUE,
                                   dom = 'rtBip', buttons = c('pageLength','copy'))) %>%
      DT::formatStyle(columns = c(8), width='250px')
  })


  output$diff_table_download <- downloadHandler(
    filename = function() {
      paste("plae_diff_table_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      gene <- input$diff_gene
      if (input$search_by == 'Gene'){
        out <- scEiaD_2020_v01 %>% tbl('wilcox_diff_testing') %>%
          filter(Gene %in% gene) %>%
          head(2000)
      } else {
        req(input$diff_term)
        test_val <- input$diff_term
        filter_term <- input$search_by
        out <- scEiaD_2020_v01 %>% tbl('wilcox_diff_testing') %>%
          filter(Group == test_val) %>%
          head(2000) %>%
          filter(Base == filter_term)
      }
      out <- out %>%
        collect() %>%
        mutate(Group = as.factor(Group)) %>%
        mutate(FDR = format(FDR, digits = 3),
               FDR = as.numeric(FDR),
               AUC = format(AUC, digits = 3),
               AUC = as.numeric(AUC),
               PValue = format(p.value, digits = 3),
               PValue = as.numeric(PValue)) %>%
        select(-`p.value`)
      write.csv(out, file)
    }
  )

  ### haystack table -----
  output$make_haystack_table <- DT::renderDataTable(server = TRUE, {
    out <- scEiaD_2020_v01 %>%
      tbl('haystack') %>%
      filter(D_KL > 0.3, !is.na(`log.p.adj`), `log.p.adj` < 0, T.counts > 500) %>%
      left_join(scEiaD_2020_v01 %>% tbl('gene_auto_label'), by = 'Gene') %>%
      arrange(log.p.adj) %>%
      select(-log.p.vals) %>%
      collect() %>%
      mutate(D_KL = format(D_KL, digits = 3) %>%as.numeric(),
             log.p.adj = round(log.p.adj)) %>%
      rename(`CellType (Predict)` = CellType_predict)


    out %>% DT::datatable(extensions = 'Buttons',
                          filter = list(position = 'bottom', clear = TRUE, plain = TRUE),
                          options = list(pageLength = 10, scrollX = TRUE,
                                         dom = 'rtBip', buttons = c('pageLength','copy'))) %>%
      DT::formatStyle(columns = c(8), width='250px')
  })

})
