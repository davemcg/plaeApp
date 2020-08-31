# server.R
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

scEiaD_2020_v01 <- dbPool(drv = SQLite(), dbname = "~/data/massive_integrated_eye_scRNA/scEiaD__2020_08_20__Mus_musculus_Macaca_fascicularis_Homo_sapiens-5000-counts-TabulaDroplet-batch-scVI-8-0.1-15-7.sqlite", idleTimeout = 3600000)
#scEiaD_2020_v01 <- dbPool(drv = SQLite(), dbname = "/data/swamyvs/plaeApp/sql_08132020.sqlite", idleTimeout = 3600000)
meta_filter <- read_fst('www/metadata_filter.fst') %>% as_tibble()
tabulamuris_predict_labels <-scEiaD_2020_v01 %>% tbl('tabulamuris_predict_labels') %>% collect
celltype_predict_labels <-scEiaD_2020_v01 %>% tbl('celltype_predict_labels') %>% collect
celltype_labels <-scEiaD_2020_v01 %>% tbl('celltype_labels') %>% collect
cluster_labels <-scEiaD_2020_v01 %>% tbl('cluster_labels')
mf <- meta_filter %>% sample_frac(0.2)

# generate color_mappings
categorical_columns <- c("Phase","batch","study_accession","library_layout","organism","Platform",
                         "Covariate","CellType","CellType_predict","TabulaMurisCellType","TabulaMurisCellType_predict",
                         "GSE","Summary","Design","Citation","PMID","Stage","cluster",
                         "Doublet","TechType", "SubCellType", 'subcluster' )
#"SubCellType" and subcluster are problems
meta_filter <- meta_filter %>% mutate(SubCellType = tidyr::replace_na(SubCellType, 'None'))
map_color <- function(column, meta_filter){
  master_colorlist <- c(pals::alphabet(), pals::alphabet2())
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
                           selected = 'CRX',
                           server = TRUE)
    }
    # gene plot category filtering ----
    if (is.null(query[['gene_filter_cat']])){
      updateSelectizeInput(session, 'gene_filter_cat',
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
          selectizeInput('gene_filter_on', strong('Filter on: '),
                         choices = choice, selected = NULL, multiple = TRUE)
        } else {
          shinyWidgets::setSliderColor(c("#3399ff"), c(1))
          sliderInput("gene_filter_on", label = strong("Filter Range: "), min = min(choice),
                      max = max(choice), value = c(min(choice), max(choice)))
        }
      })
    })

    # meta plot updateSelectizeInput ------
    if (is.null(query[['meta_column']])){
      updateSelectizeInput(session, 'meta_column',
                           choices = meta_filter %>%
                             dplyr::select(nCount_RNA:doublet_score_scran) %>% colnames() %>% sort(),
                           options = list(placeholder = 'Type to search'),
                           selected = 'CellType_predict',
                           server = TRUE)
    }
    # meta plot category filtering ----
    if (is.null(query[['meta_filter_cat']])){
      updateSelectizeInput(session, 'meta_filter_cat',
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
          selectizeInput('meta_filter_on', strong('Filter on: '),
                         choices = choice, selected = NULL, multiple = TRUE)
        } else {
          shinyWidgets::setSliderColor(c("#3399ff"), c(1))
          sliderInput("meta_filter_on", label = strong("Filter Range: "), min = min(choice),
                      max = max(choice), value = c(min(choice), max(choice)))
        }
      })

    })

    # dotplot updateSelectizeInput ----
    if (is.null(query[['dotplot_Gene']])){
      updateSelectizeInput(session, 'dotplot_Gene',
                           choices = scEiaD_2020_v01 %>% tbl('genes') %>% collect() %>% pull(1),
                           options = list(placeholder = 'Type to search'),
                           selected = c('RHO','WIF1','CABP5', 'AIF1','AQPT4','ARR3','ONECUT1','GRIK1','GAD1','POU4F2'),
                           server = TRUE)
    }
    if (is.null(query[['dotplot_groups']])){
      updateSelectizeInput(session, 'dotplot_groups',
                           choices = meta_filter %>%
                             select(-Barcode) %>%
                             select_if(purrr::negate(is.numeric)) %>%
                             colnames() %>% sort(),
                           options = list(placeholder = 'Type to search',
                                          maxItems = 2),
                           selected = c('CellType_predict','organism'),
                           server = TRUE)
    }
    if (is.null(query[['dotplot_filter_cat']])){
      updateSelectizeInput(session, 'dotplot_filter_cat',
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
                           choices = choice,
                           server = TRUE)
    })

    # insitu updateSelectizeInput ----
    if (is.null(query[['insitu_Gene']])){
      updateSelectizeInput(session, 'insitu_Gene',
                           choices = scEiaD_2020_v01 %>% tbl('genes') %>% as_tibble() %>% pull(1),
                           options = list(placeholder = 'Type to search'),
                           selected = c('RHO'),
                           server = TRUE)
    }

    if (is.null(query[['insitu_filter_cat']])){
      updateSelectizeInput(session, 'insitu_filter_cat',
                           choices = scEiaD_2020_v01 %>% tbl('grouped_stats') %>%
                             select(-Gene, -cell_ct, -cell_exp_ct, -cpm) %>% colnames() %>% sort(),
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
                           choices = choice,
                           server = TRUE)
    })

    #
    if (is.null(query[['grouping_features']])){
      updateSelectizeInput(session, 'grouping_features',
                           choices = scEiaD_2020_v01 %>% tbl('grouped_stats') %>%
                             select(-Gene, -cell_ct, -cell_exp_ct, -cpm) %>% colnames() %>% sort(),
                           options = list(placeholder = 'Type to search'),
                           selected = c('CellType_predict'),
                           server = TRUE)
    }

    if (is.null(query[['meta_groupings']])){
      updateSelectizeInput(session, 'meta_groupings',
                           choices = meta_filter %>%
                             select(-Barcode, -UMAP_1, -UMAP_2, -nCount_RNA, -nFeature_RNA, -percent_mt) %>%
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
                           selected = c('PAX6','POU4F2','CRX','NRL'),
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
                             select(-Gene, -cell_ct, -cell_exp_ct, -cpm) %>% colnames() %>% sort(),
                           selected = 'CellType',
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
                           selected = c('Cones','Retinal Ganglion', 'Horizontal Cells'),
                           server = TRUE)
    })
    if (is.null(query[['exp_plot_groups']])){
      updateSelectizeInput(session, 'exp_plot_groups',
                           choices = scEiaD_2020_v01 %>% tbl('grouped_stats') %>%
                             select(-Gene, -cell_ct, -cell_exp_ct, -cpm) %>% colnames() %>% sort(),
                           selected = c('study_accession'),
                           server = TRUE)
    }



    # temporal plot updateSelect -----
    if (is.null(query[['temporal_gene']])){
      updateSelectizeInput(session, 'temporal_gene',
                           choices = scEiaD_2020_v01 %>% tbl('genes') %>% collect() %>% pull(1),
                           options = list(placeholder = 'Type to search'),
                           selected = c('PAX6','POU4F2'),
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
                             select(-Gene, -cell_ct, -cell_exp_ct, -cpm) %>% colnames() %>% sort(),
                           selected = 'CellType',
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
                           selected = c('Cones','Retinal Ganglion', 'Horizontal Cells'),
                           server = TRUE)
    })

    if (is.null(query[['facet_color']])){
      updateSelectizeInput(session, 'facet_color',
                           choices = meta_filter %>%
                             dplyr::select_if(is.character) %>% colnames() %>% sort(),
                           options = list(placeholder = 'Type to search'),
                           selected = 'CellType',
                           server = TRUE)
    }
    # diff table updateSelect ------
    if (is.null(query[['diff_gene']])){
      updateSelectizeInput(session, 'diff_gene',
                           choices = scEiaD_2020_v01 %>% tbl('genes') %>% collect() %>% pull(1),
                           options = list(placeholder = 'Type to search'),
                           selected = 'CRX',
                           server = TRUE)
    }
    if (is.null(query[['diff_term']])){
      term = input$search_by
      updateSelectizeInput(session, 'diff_term',
                           choices = scEiaD_2020_v01 %>%
                             tbl('PB_Test_terms') %>%
                             filter(PB_Test == term) %>%
                             collect() %>% pull(terms) %>%
                             strsplit(., '___') %>% unlist(),
                           options = list(placeholder = 'Type to search'),
                           server = TRUE)
    }
    # Meta Plot modal ----------
    observeEvent(input$BUTTON_show_meta_legend, {
      # Show a modal when the button is pressed
      showModal(shinyjqui::draggableModalDialog(size = 'l', title = 'Click to Drag',
                                                plotOutput('meta_plot_legend') %>% shinycssloaders::withSpinner(type = 3, size = 0.5, color = "#3399ff", color.background = 'white'),
                                                easyClose = TRUE))
    })
    # HELP button descriptions ----------
    ## umap
    observeEvent(input$umap_table_help, {
      # Show a modal when the button is pressed
      showModal(shinyjqui::draggableModalDialog(size = 'm',
                                                title = "UMAP - Tables",
                                                HTML("<p>The UMAP is a 2D Projection of a higher dimensional space which tries to bring together
                            closely related elements while maintaining the overall structure. The higher dimensional space
                            is built from the gene expression patterns of each cell. The left panel shows the expression pattern
                            of a gene in the UMAP space. The right panel shows metadata associated with each cell, including
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
                              <li>Expression, which is log2 scaled CPM (counts per million) </li>
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
                                                HTML("<p>We have pre-computed 12 different differential expression tests. They
                                                     can be grouped into 3 categories:
                                                     <ul>
                              <li>[ ] against Remaining, which tests [ ] against all other cells. The effect of organism is controlled
                              by giving it as a covariate in the test</li>
                              <li>Pairwise [ ] against [ ], which tests genes differentially expressed in pairwise combinations
                              (for example Rods against Cones, ignoring all other cells)</li>
                              <li>Organism specific test within [ ]. For example you can search for genes differentially expressed
                              between mouse and human WITHIN rods.</li>
                            </ul>
                            <p>[ ] is either:</p>
                            <ul>
                              <li>CellType, which are based on published cell type assignments</li>
                              <li>CellType (predict), which uses ML to project CellType labels onto (nearly) all of the cells</li>
                              <li>Cluster (droplet or well), which groups the droplet or well (e.g. 10X or SmartSeq) based cells into clusters in an
                              unsupervised manner. Well and droplet were clustered separately as the integration performance
                              was suboptimal when combining these two technologies.</li>
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
    if (input$gene_and_meta_scatter_tech == 'Droplet'){
      temp_filter <- meta_filter %>% filter(TechType == 'Droplet')
      x_range =  c(temp_filter$UMAP_1 %>% min(), temp_filter$UMAP_1 %>% max())
      y_range = c(temp_filter$UMAP_2 %>% min(), temp_filter$UMAP_2 %>% max())
    } else {
      temp_filter <- meta_filter %>% filter(TechType == 'Well')
      x_range =  c(temp_filter$UMAP_1 %>% min(), temp_filter$UMAP_1 %>% max())
      y_range = c(temp_filter$UMAP_2 %>% min(), temp_filter$UMAP_2 %>% max())
    }
    gene_scatter_ranges <- reactiveValues(x = x_range,
                                          y = y_range)
    source('make_gene_scatter_umap_plot.R')
    gene_scatter_plot <- eventReactive(input$BUTTON_draw_scatter, {
      make_gene_scatter_umap_plot(input, scEiaD_2020_v01, mf, meta_filter)
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
    })


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
    if (input$gene_and_meta_scatter_tech == 'Droplet'){
      temp_filter <- meta_filter %>% filter(TechType == 'Droplet')
      x_range =  c(temp_filter$UMAP_1 %>% min(), temp_filter$UMAP_1 %>% max())
      y_range = c(temp_filter$UMAP_2 %>% min(), temp_filter$UMAP_2 %>% max())
    } else {
      temp_filter <- meta_filter %>% filter(TechType == 'Well')
      x_range =  c(temp_filter$UMAP_1 %>% min(), temp_filter$UMAP_1 %>% max())
      y_range = c(temp_filter$UMAP_2 %>% min(), temp_filter$UMAP_2 %>% max())
    }
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
    })

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
        summarise(cpm = sum(cpm * cell_exp_ct) / sum(cell_exp_ct),
                  cell_exp_ct = sum(cell_exp_ct, na.rm = TRUE)) %>%
        collect() %>%
        tidyr::drop_na() %>%
        full_join(., meta_filter %>%
                    group_by_at(vars(one_of(grouping_features))) %>%
                    summarise(Count = n())) %>%
        mutate(cell_exp_ct = ifelse(is.na(cell_exp_ct), 0, cell_exp_ct)) %>%
        mutate(`%` = round((cell_exp_ct / Count) * 100, 2),
               Expression = round(cpm * (`%` / 100), 2)) %>%
        select_at(vars(one_of(c('Gene', grouping_features, 'cell_exp_ct', 'Count', '%', 'Expression')))) %>%
        arrange(-Expression) %>%
        rename(`Cells # Detected` = cell_exp_ct,
               `Total Cells` = Count,
               `log2(cpm+1)` = Expression)

      table %>% DT::datatable(extensions = 'Buttons',
                              filter = list(position = 'bottom', clear = FALSE),
                              options = list(pageLength = 10, dom = 'frtBip', buttons = c('pageLength','copy', 'csv')))
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
                              filter = list(position = 'bottom', clear = FALSE),
                              options = list(pageLength = 10, dom = 'frtBip', buttons = c('pageLength','copy', 'csv')))
    })
    output$metadata_stats <- DT::renderDataTable({ metadata_stats()})

    # facet plot -----------
    source('make_facet_plot.R')
    facet_plot <- eventReactive(input$BUTTON_draw_filter, {make_facet_plot(input,  meta_filter)})

    output$facet_plot <- renderPlot({
      facet_plot()
    }, height = eventReactive(input$BUTTON_draw_filter, {input$facet_height %>% as.numeric()}))

    ## exp_plot -----------
    source('make_exp_plot.R')
    exp_plot <- eventReactive(input$BUTTON_draw_exp_plot, {
      make_exp_plot(input, scEiaD_2020_v01, meta_filter)
    })

    output$exp_plot <- renderPlot({
      exp_plot()
    }, height = eventReactive(input$BUTTON_draw_exp_plot, {as.numeric(input$exp_plot_height)}))



    ## temporal plot -----------
    source('make_temporal_plot.R')
    temporal_plot <- eventReactive(input$BUTTON_draw_temporal, {
      make_temporal_plot(input, scEiaD_2020_v01, meta_filter)
    })
    output$temporal_plot <- renderPlot({
      temporal_plot()
    }, height = as.numeric(input$temporal_plot_height ) )

    ## dotplot ---------
    source('make_dotplot.R')
    draw_dotplot <- eventReactive(input$BUTTON_draw_dotplot,
                                  {make_dotplot(input, scEiaD_2020_v01, meta_filter,cat_to_color_df)}
    )
    output$dotplot <- renderPlot({
      draw_dotplot()
    }, height = eventReactive(input$BUTTON_draw_dotplot, {input$dotplot_height %>% as.numeric()}))
  })

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
             `log2(cpm+1)` = Expression) %>%
      tidyr::drop_na()
    full_table %>% DT::datatable(extensions = 'Buttons',
                                 filter = list(position = 'bottom', clear = FALSE),
                                 options = list(pageLength = 10, dom = 'frtBip', buttons = c('pageLength','copy', 'csv')))
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

  ## diff testing table ----
  # diff_table <- reactive({
  #   req(input$diff_table_select)
  #   if (input$diff_table_select == 'Cluster'){
  #     out <- 'diff_testing_A'
  #   } else if (input$diff_table_select == 'CellType') {
  #     out <- 'diff_testing_E'
  #   } else if (input$diff_table_select == 'SubCluster') {
  #     out <- 'diff_testing_G'
  #   } else {out <- 'diff_testing_C'}
  #   return(out)
  # })

  output$make_diff_table <- DT::renderDataTable(server = TRUE, {
    gene <- input$diff_gene
    #cat(diff_table)
    if (input$search_by == 'Gene'){
      out <- scEiaD_2020_v01 %>% tbl('PB_results') %>%
        filter(Gene %in% gene, FDR < 0.05, abs(logFC) > 0.5) %>%
        arrange(FDR)
    } else {
      # isolate({
      req(input$diff_term)
      test_val <- input$diff_term
      filter_term <- input$search_by
      out <- scEiaD_2020_v01 %>% tbl('PB_results') %>%
        filter(test == test_val, FDR < 0.05, abs(logFC) > 0.5) %>%
        arrange(FDR) %>%
        collect() %>%
        filter(PB_Test == filter_term)
      #})
    }
    out %>%
      # mutate(`AUC Score` = round(count / 55, 1),
      #        mean_auc = round(mean_auc, 2)) %>%
      # select(-status, -model_component, -count, -med_auc, -estimate, -test_val, -p_value) %>%
      # collect() %>%
      # filter(!grepl('Doub|RPE|Astro|Red|Vascular', term)) %>%
      collect() %>%
      select(-comparison) %>%
      mutate(PB_Test = as.factor(PB_Test)) %>%
      mutate(FDR = format(FDR, digits = 3),
             FDR = as.numeric(FDR),
             PValue = format(PValue, digits = 3),
             PValue = as.numeric(PValue)) %>%
      DT::datatable(extensions = 'Buttons',
                    filter = list(position = 'bottom', clear = FALSE),
                    options = list(pageLength = 10,
                                   dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatRound(columns = c('logFC','logCPM','F'), digits = 2) %>%
      DT::formatStyle(columns = c(8), width='250px')

  })

  # output$formattable01 <- renderFormattable({formattable_01})
  # output$formattable02 <- renderFormattable({formattable_02})
})
