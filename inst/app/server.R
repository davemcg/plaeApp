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


scEiaD_2020_v01 <- dbPool(drv = SQLite(), dbname = "~/data/massive_integrated_eye_scRNA/MOARTABLES__anthology_limmaFALSE___Mus_musculus_Macaca_fascicularis_Homo_sapiens-5000-counts-TabulaDroplet-batch-scVI-8-0.1-15-7.sqlite", idleTimeout = 3600000)

# these will be pre-processed and moved
# into the sqlite db when they are finalized
#anthology_2020_v01 <- dbPool(drv = SQLite(), dbname = "~/data/massive_integrated_eye_scRNA/MOARTABLES__anthology_limmaFALSE___Mus_musculus_Macaca_fascicularis_Homo_sapiens-2000-counts-onlyDROPLET-batch-scVI-6-0.1-500-10.sqlite", idleTimeout = 3600000)

# fancy tables
# they come from `tables.Rmd` in analysis/
# load('www/formattables.Rdata')
# filter
meta_filter <- left_join(scEiaD_2020_v01 %>% tbl('metadata_filter'),
                         scEiaD_2020_v01 %>% tbl('doublets'), by ='Barcode') %>%
  collect() %>%
  mutate(`Doublet Probability` = as.numeric(`Doublet Probability`),
         doublet_score_scran = as.numeric(doublet_score_scran)) %>%
  mutate(PMID = as.character(PMID)) %>%
  mutate(Tech = case_when(Platform %in% c('10xv2','10xv3','DropSeq') ~ 'Droplet',
                          TRUE ~ 'Well'))

# cutdown mf for plotting
mf <- meta_filter %>% sample_frac(0.2)
# get coords for cell labels
celltype_predict_labels <-
  bind_rows(meta_filter %>%
              filter(Platform %in% c('10xv2','10xv3','DropSeq')) %>%
              group_by(CellType_predict) %>%
              summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2)) %>%
              mutate(Tech = 'Droplet'),
            meta_filter %>%
              filter(!Platform %in% c('10xv2','10xv3','DropSeq')) %>%
              group_by(CellType_predict) %>%
              summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2)) %>%
              mutate(Tech = 'Well'))
celltype_labels <- meta_filter %>%
  group_by(CellType) %>%
  summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2)) %>%
  mutate(Tech = 'Droplet')
# tabulamuris_labels <- meta_filter %>%
#   group_by(TabulaMurisCellType) %>%
#   summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
tabulamuris_predict_labels <- meta_filter %>%
  group_by(TabulaMurisCellType_predict) %>%
  summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2)) %>%
  mutate(Tech = 'Droplet')
# get coords for cell labels
cluster_labels <-
  bind_rows(meta_filter %>%
              filter(Platform %in% c('10xv2','10xv3','DropSeq')) %>%
              group_by(cluster) %>% summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2)) %>%
              mutate(Tech = 'Droplet'),
            meta_filter %>%
              filter(!Platform %in% c('10xv2','10xv3','DropSeq')) %>%
              group_by(cluster) %>% summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2)) %>%
              mutate(Tech = 'Well'))
#----


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
    observeEvent(input$facet, {
      if(input$facet == ''){
        choice=''
      }else {
        choice = meta_filter[,input$facet] %>% pull(1) %>% unique() %>% sort()
      }
      updateSelectizeInput(session, 'facet_filter',
                           choices = choice,
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
    # BREAK -------
    # gene scatter plot ------------
    gene_scatter_ranges <- reactiveValues(x = c(meta_filter$UMAP_1 %>% min(), meta_filter$UMAP_1 %>% max()),
                                          y = c(meta_filter$UMAP_2 %>% min(), meta_filter$UMAP_2 %>% max()))
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
        gene_scatter_ranges$x <- c(meta_filter$UMAP_1 %>% min(), meta_filter$UMAP_1 %>% max())
        gene_scatter_ranges$y <- c(meta_filter$UMAP_2 %>% min(), meta_filter$UMAP_2 %>% max())
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
                                  cluster_labels
      )
    })

    meta_ranges <- reactiveValues(x = c(meta_filter$UMAP_1 %>% min(), meta_filter$UMAP_1 %>% max()),
                                  y = c(meta_filter$UMAP_2 %>% min(), meta_filter$UMAP_2 %>% max()))
    observeEvent(input$meta_plot_dblclick, {
      brush <- input$meta_plot_brush
      if (!is.null(brush)) {
        meta_ranges$x <- c(brush$xmin, brush$xmax)
        meta_ranges$y <- c(brush$ymin, brush$ymax)

      } else {
        meta_ranges$x <- c(meta_filter$UMAP_1 %>% min(), meta_filter$UMAP_1 %>% max())
        meta_ranges$y <- c(meta_filter$UMAP_2 %>% min(), meta_filter$UMAP_2 %>% max())
      }
    })

    output$meta_plot <- renderPlot({
      meta_plot() + coord_cartesian(xlim = meta_ranges$x, ylim = meta_ranges$y)
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

      table %>% DT::datatable(extensions = 'Buttons', rownames = F,
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

      table %>% DT::datatable(extensions = 'Buttons', rownames = F,
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
                                  {make_dotplot(input, scEiaD_2020_v01, meta_filter)}
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
    full_table %>% DT::datatable(extensions = 'Buttons', rownames = F,
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
        filter(Gene %in% gene) %>%
        arrange(FDR)
    } else {
      # isolate({
      req(input$diff_term)
      test_val <- input$diff_term
      filter_term <- input$search_by
      out <- scEiaD_2020_v01 %>% tbl('PB_results') %>%
        filter(test == test_val) %>%
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
      DT::datatable(extensions = 'Buttons', rownames = F,
                    filter = list(position = 'bottom', clear = FALSE),
                    options = list(pageLength = 10,
                                   dom = 'frtBip', buttons = c('pageLength','copy', 'csv'))) %>%
      DT::formatRound(columns = c('logFC','logCPM','F'), digits = 2) %>%
      DT::formatStyle(columns = c(8), width='250px')

  })

  # output$formattable01 <- renderFormattable({formattable_01})
  # output$formattable02 <- renderFormattable({formattable_02})
})
