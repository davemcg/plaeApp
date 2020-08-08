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

#anthology_2020_v01 <- dbPool(drv = SQLite(), dbname = "~/data/massive_integrated_eye_scRNA/MOARTABLES__anthology_limmaFALSE___Mus_musculus_Macaca_fascicularis_Homo_sapiens-2000-counts-onlyDROPLET-batch-scVI-6-0.1-500-10.sqlite", idleTimeout = 3600000)

#scEiaD_2020_v01 <- dbPool(drv = SQLite(), dbname = "~/data/massive_integrated_eye_scRNA/MOARTABLES__anthology_limmaFALSE___Mus_musculus_Macaca_fascicularis_Homo_sapiens-5000-counts-TabulaDroplet-batch-scVI-8-0.1-15-7.sqlite", idleTimeout = 3600000)
scEiaD_2020_v01 <- dbPool(drv = SQLite(), dbname = "/data/swamyvs/plaeApp/MOARTABLES__anthology_limmaFALSE___Mus_musculus_Macaca_fascicularis_Homo_sapiens-5000-counts-TabulaDroplet-batch-scVI-8-0.1-15-7.sqlite", idleTimeout = 3600000)
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
      #species <- strsplit(query$species, "\"")[[1]][2]
      updateTabsetPanel(session, 'nav', query$url)
      #updateSelectInput(session, 'species',selected = species)
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

    gene_scatter_plot <- eventReactive(input$BUTTON_draw_scatter, {
      cat(file=stderr(), paste0(Sys.time(), ' Gene Scatter Plot Call\n'))
      gene <- input$Gene
      tech <- input$gene_and_meta_scatter_tech
      pt_size <- input$pt_size_gene %>% as.numeric()
      expression_range <- input$gene_scatter_slider
      mf <- mf %>% filter(Tech == tech)
      p <-  scEiaD_2020_v01 %>% tbl('cpm') %>%
        filter(Gene == gene) %>%
        collect() %>%
        mutate(cpm = cpm - min(cpm) + 1) %>%
        filter(cpm > as.numeric(expression_range[1]),
               cpm < as.numeric(expression_range[2])) %>%
        left_join(., meta_filter, by = 'Barcode') %>%
        filter(Tech == tech, !is.na(UMAP_1), !is.na(UMAP_2), !is.na(cpm))
      cat(input$gene_filter_cat)
      cat(class(input$gene_filter_cat))
      if (!is_null(input$gene_filter_cat)){
        if (class(input$gene_filter_on) == 'character'){
          p <- p %>%
            #filter(!!as.symbol(input$gene_filter_cat) %in% input$gene_filter_on)
            filter_at(vars(all_of(input$gene_filter_cat)), all_vars(. %in% input$gene_filter_on))
        } else {
          p <- p %>%
            filter(!!as.symbol(input$gene_filter_cat) >= input$gene_filter_on[1],
                   !!as.symbol(input$gene_filter_cat) <= input$gene_filter_on[2])
        }
      }
      color_range <- range(p$cpm)
      plot <- p %>% ggplot() +
        geom_scattermost(cbind(mf$UMAP_1, mf$UMAP_2), color = '#D3D3D333',
                         pointsize = pt_size * 1.5,
                         pixels=c(750,750)) +
        geom_scattermost(cbind(p$UMAP_1, p$UMAP_2),
                         color = viridis::magma(100, alpha=0.3)
                         [1+99*(p$cpm-color_range[1])/diff(color_range)],
                         pointsize= pt_size,
                         pixels=c(750,750),
                         interpolate=FALSE) +
        geom_point(data=data.frame(x=double(0)), aes(x,x,color=x)) +
        scale_color_gradientn(  #add the manual guide for the empty aes
          limits=c(min(p$cpm),max(p$cpm)),
          colors=viridis::magma(100),
          name="log2(cpm+1)") +
        theme_cowplot() +
        theme(axis.line = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank()) +
        annotate("text", -Inf, Inf, label = paste0(gene, ' expression'), hjust = 0, vjust = 1, size = 6)

      suppressWarnings(plot)

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
    meta_plot <- eventReactive(input$BUTTON_draw_meta, {
      cat(file=stderr(), paste0(Sys.time(), ' Meta Plot Call\n'))
      meta_column <- input$meta_column
      transform <- input$meta_column_transform

      pt_size <- input$pt_size_meta %>% as.numeric()
      filter_column <- input$meta_column
      # cut down to match tech selected
      tech <- input$gene_and_meta_scatter_tech
      mf <- mf %>% filter(Tech == tech)
      meta_filter <- meta_filter %>% filter(Tech == tech)
      celltype_predict_labels <- celltype_predict_labels %>% filter(Tech == tech)
      celltype_labels <- celltype_labels %>% filter(Tech == tech)
      tabulamuris_predict_labels <- tabulamuris_predict_labels %>% filter(Tech == tech)
      cluster_labels <- cluster_labels %>% filter(Tech == tech)

      if (transform == 'log2' && is.numeric(meta_filter[,meta_column] %>% pull(1))){
        cat('log2 time')
        meta_filter[,meta_column] <- log2(meta_filter[,meta_column] + 1)
      }
      p_data <- meta_filter %>%
        filter(!grepl('Doub|\\/Margin\\/Periocular', CellType)) %>%
        filter(!is.na(!!as.symbol(meta_column)))
      # category filtering
      if (!is_null(input$meta_filter_cat)){
        if (class(input$meta_filter_on) == 'character'){
          p_data <- p_data %>%
            #filter(!!as.symbol(input$meta_filter_cat) %in% input$meta_filter_on)
            filter_at(vars(all_of(input$meta_filter_cat)), all_vars(. %in% input$meta_filter_on))
        } else {
          p_data <- p_data %>%
            filter(!!as.symbol(input$meta_filter_cat) >= input$meta_filter_on[1],
                   !!as.symbol(input$meta_filter_cat) <= input$meta_filter_on[2])
        }
      }

      # metadata NUMERIC plot --------------
      if (is.numeric(meta_filter[,meta_column] %>% pull(1)) ){
        color_range <- range(p_data[,meta_column] %>% pull(1))
        suppressWarnings(plot <- ggplot() +
                           geom_scattermost(cbind(mf %>%
                                                    filter(is.na(!!as.symbol(meta_column))) %>% pull(UMAP_1),
                                                  mf %>%
                                                    filter(is.na(!!as.symbol(meta_column))) %>% pull(UMAP_2)),
                                            pointsize = pt_size, color = '#D3D3D333',
                                            pixels = c(750,750)) +
                           geom_scattermost(cbind(p_data$UMAP_1, p_data$UMAP_2),
                                            color = viridis::viridis(100, alpha=0.3)
                                            [1+99*((p_data[,meta_column] %>% pull(1))-color_range[1])/diff(color_range)],
                                            pointsize= pt_size,
                                            pixels=c(750,750),
                                            interpolate=FALSE) +
                           geom_point(data=data.frame(x=double(0)), aes(x,x,color=x))  +
                           scale_color_gradientn(  #add the manual guide for the empty aes
                             limits=c(min(p_data[,meta_column] %>% pull(1)),
                                      max(p_data[,meta_column] %>% pull(1))),
                             colors=viridis::viridis(100),
                             name=meta_column) +
                           guides(colour = guide_legend(override.aes = list(alpha = 1))) +
                           theme_cowplot() +
                           theme(axis.line = element_blank(),
                                 axis.title = element_blank(),
                                 axis.ticks = element_blank(),
                                 axis.text = element_blank()) +
                           annotate("text", -Inf, Inf, label = "Metadata", hjust = 0, vjust = 1, size = 6))
        # metadata CATEGORICAL plot --------------
      } else {
        group <- p_data[,meta_column] %>% pull(1) %>% as.factor()
        p_color = rep(c(pals::alphabet(),pals::alphabet2()), times = 20)[group] %>% paste0(., '33') # alpha 0.33
        color_data <- group %>% levels() %>% tibble::enframe() %>% mutate(x=0)
        suppressWarnings(plot <- ggplot() +
                           geom_scattermost(cbind(mf %>%
                                                    filter(is.na(!!as.symbol(meta_column))) %>% pull(UMAP_1),
                                                  mf %>%
                                                    filter(is.na(!!as.symbol(meta_column))) %>% pull(UMAP_2)),
                                            pointsize = pt_size, color = '#D3D3D333',
                                            pixels = c(750,750)) +
                           geom_scattermost(cbind(p_data$UMAP_1, p_data$UMAP_2),
                                            color = p_color ,
                                            pointsize= pt_size,
                                            pixels=c(750,750),
                                            interpolate=FALSE) +
                           geom_point(data=color_data, aes(x,x,color=value), alpha = 0) +
                           scale_colour_manual(name= meta_column,
                                               values = rep(c(pals::alphabet() %>% unname(),
                                                              pals::alphabet2() %>% unname()),
                                                            times = 20)) +
                           guides(colour = guide_legend(override.aes = list(alpha = 1, size = 7))) +
                           theme_cowplot() +
                           theme(axis.line = element_blank(),
                                 axis.title = element_blank(),
                                 axis.ticks = element_blank(),
                                 axis.text = element_blank()) +
                           annotate("text", -Inf, Inf, label = "Metadata", hjust = 0, vjust = 1, size = 6))
      }

      more <- NULL
      if ('1' %in% input$label_toggle){
        more <- geom_text_repel(data = celltype_labels, bg.color = 'white',
                                aes(x = UMAP_1, y = UMAP_2, label = CellType))
      }
      if ('2' %in% input$label_toggle){
        more <- geom_text_repel(data = celltype_predict_labels, bg.color = 'white',
                                aes(x = UMAP_1, y = UMAP_2, label = CellType_predict))
      }
      if ('3' %in% input$label_toggle){
        more <- geom_text_repel(data = cluster_labels, bg.color = 'white',
                                aes(x = UMAP_1, y = UMAP_2, label = cluster),
                                max.iter = 20)
      }
      if ('4' %in% input$label_toggle){
        more <- geom_text_repel(data = tabulamuris_predict_labels, bg.color = 'white',
                                aes(x = UMAP_1, y = UMAP_2, label = TabulaMurisCellType_predict),
                                max.iter = 20)
      }
      if (meta_column %in% c('cluster','subcluster', 'TabulaMurisCellType', 'TabulaMurisCellType_predict')){
        suppressWarnings(plot + more + theme(legend.position = 'none'))
      } else {
        suppressWarnings(plot + more)
      }

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
    facet_plot <- eventReactive(input$BUTTON_draw_filter, {
      cat(file=stderr(), paste0(Sys.time(), ' Facet Plot Call\n'))
      facet_column <- input$facet
      color_column <- input$facet_color
      #transform <- input$facet_column_transform
      pt_size <- input$pt_size_facet %>% as.numeric()

      gray_data <- meta_filter %>%
        filter(is.na(!!as.symbol(color_column)))
      p_data <- meta_filter %>%
        filter(!is.na(!!as.symbol(facet_column)),
               !is.na(!!as.symbol(color_column)))

      suppressWarnings(plot <- ggplot(data = p_data) +
                         geom_scattermore(data = gray_data,
                                          aes(x = UMAP_1, y = UMAP_2),
                                          color = 'gray',
                                          pointsize = pt_size,
                                          pixels = c(750,750),
                                          alpha = 0.4) +
                         geom_scattermore(aes(x = UMAP_1, y = UMAP_2,
                                              color = !!as.symbol(color_column)) ,
                                          pointsize= pt_size,
                                          pixels = c(750,750),
                                          alpha = 0.6) +
                         facet_wrap(vars(!!(as.symbol(facet_column)))) +
                         scale_colour_manual(values = rep(c(pals::alphabet() %>% unname(),
                                                            pals::alphabet2() %>% unname()),
                                                          times = 20),
                                             na.value = 'gray') +
                         guides(colour = guide_legend(override.aes = list(alpha = 1, size = 7))) +
                         theme_cowplot() +
                         theme(axis.line = element_blank(),
                               axis.title = element_blank(),
                               axis.ticks = element_blank(),
                               axis.text = element_blank())
      )
      plot

    })

    output$facet_plot <- renderPlot({
      facet_plot()
    }, height = eventReactive(input$BUTTON_draw_filter, {input$facet_height %>% as.numeric()}))

    ## exp_plot -----------
    exp_plot <- eventReactive(input$BUTTON_draw_exp_plot, {
      cat(file=stderr(), paste0(Sys.time(), ' Exp Plot Call\n'))
      gene <- input$exp_plot_genes

      grouping_features <- input$exp_plot_groups

      if (input$exp_filter_cat != ''){
        # box_data <- scEiaD_2020_v01 %>% tbl('grouped_stats') %>%
        #   filter(Gene %in% gene) %>%
        #   collect() %>%
        #   filter(!!as.symbol(input$exp_filter_cat) %in% input$exp_filter_on)
        #

        if (class(input$exp_filter_on) == 'character'){
          box_data <- scEiaD_2020_v01 %>% tbl('grouped_stats') %>%
            filter(Gene %in% gene) %>%
            collect() %>%
            filter_at(vars(all_of(input$exp_filter_cat)), all_vars(. %in% input$exp_filter_on))
            #filter(!!as.symbol(input$exp_filter_cat) %in% input$exp_filter_on)
        } else {
          box_data <- scEiaD_2020_v01 %>% tbl('grouped_stats') %>%
            filter(Gene %in% gene) %>%
            collect() %>%
            #filter(!!as.symbol(input$exp_filter_cat) %in% input$exp_filter_on) %>%
            filter(!!as.symbol(input$exp_filter_cat) >= input$exp_filter_on[1],
                   !!as.symbol(input$exp_filter_cat) <= input$exp_filter_on[2])
        }

      } else {
        box_data <- scEiaD_2020_v01 %>% tbl('grouped_stats') %>%
          filter(Gene %in% gene) %>%
          collect()
      }
      validate(
        need(input$exp_plot_groups != '', "Please select at least one grouping feature")
      )

      #cat(input)
      box_data <- box_data %>%
        #filter(!is.na(!!as.symbol(grouping_features))) %>%
        group_by_at(vars(one_of(c('Gene', input$exp_plot_facet, grouping_features)))) %>%
        summarise(cpm = sum(cpm * cell_exp_ct) / sum(cell_exp_ct),
                  cell_exp_ct = sum(cell_exp_ct, na.rm = TRUE)) %>%
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
               `Mean CPM` = Expression,
               `% of Cells Detected` = `%`) %>%
        tidyr::drop_na()
      box_data$Group <- box_data[,c(2:(length(grouping_features)+1))] %>% tidyr::unite(x, sep = ' ') %>% pull(1)

      box_data %>%
        ggplot(aes(x=Gene, y = !!as.symbol(input$exp_plot_ylab), color = !!as.symbol(grouping_features))) +
        geom_boxplot(color = 'black', outlier.shape = NA) +
        ggbeeswarm::geom_quasirandom(aes(size = `Total Cells`), grouponX = TRUE) +
        cowplot::theme_cowplot() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        scale_colour_manual(values = rep(c(pals::alphabet() %>% unname()), 20)) +
        theme(legend.position="bottom") +
        facet_wrap(ncol = as.numeric(input$exp_plot_col_num), scales = 'free_x', vars(!!as.symbol(input$exp_plot_facet)))
    })

    output$exp_plot <- renderPlot({
      exp_plot()
    }, height = eventReactive(input$BUTTON_draw_exp_plot, {as.numeric(input$exp_plot_height)}))



    ## temporal plot -----------
    temporal_plot <- eventReactive(input$BUTTON_draw_temporal, {
      cat(file=stderr(), paste0(Sys.time(), ' Temporal Plot Call\n'))
      gene <- input$temporal_gene
      grouping <- input$temporal_group
      y_val <- input$temporal_y_val
      if (grouping == 'CellType (predict)'){grouping <- 'CellType_predict'}
      if (y_val == 'Mean CPM') {y_val <- 'cpm'} else {y_val <- 'Ratio'}
      meta_data <- meta_filter %>%
        group_by(organism, !!as.symbol(grouping), Age) %>%
        summarise(full_count = n())
      temporal_data <- scEiaD_2020_v01 %>% tbl('cpm') %>%
        filter(Gene %in% gene) %>%
        collect() %>%
        mutate(cpm = cpm - min(cpm) + 1) %>%
        left_join(., meta_filter, by = 'Barcode') %>%
        filter(!is.na(!!as.symbol(grouping)), !grepl('Doub|RPE', !!as.symbol(grouping))) %>%
        group_by(organism, !!as.symbol(grouping), Age, Gene) %>%
        summarise(cpm = mean(cpm), count = n()) %>%
        right_join(., meta_data) %>%
        mutate(count = ifelse(is.na(count), 0, count)) %>%
        mutate(Ratio = count/full_count) %>%
        filter(!is.na(!!as.symbol(grouping)), !is.na(Gene)) %>%
        ungroup() %>%
        mutate(Age = case_when(organism == 'Mus musculus' & Age == 1000 ~ 20,
                               organism == 'Homo sapiens' & Age == 1000 ~ 65,
                               is.na(Age) ~ 65,
                               Age == 31360 ~ 65,
                               TRUE ~ Age))


      suppressWarnings(temporal_human <-  temporal_data %>%
                         filter(organism == 'Homo sapiens') %>%
                         ggplot(aes(x=Age, y = !!as.symbol(y_val), color = Gene)) +
                         geom_point(stat = 'identity') +
                         ggtitle('Human') +
                         geom_line() + ylab(input$temporal_y_val) +
                         cowplot::theme_cowplot() + xlab('Age (days from birth)') +
                         facet_wrap(vars(!!as.symbol(grouping))) +
                         theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
                         scale_colour_manual(values = rep(c(pals::alphabet() %>% unname()))))
      suppressWarnings(temporal_mouse <- temporal_data %>%
                         filter(organism == 'Mus musculus') %>%
                         ggplot(aes(x=Age, y = !!as.symbol(y_val), color = Gene)) +
                         geom_point(stat = 'identity') +
                         ggtitle('Mouse') +
                         geom_line() + ylab(input$temporal_y_val) +
                         cowplot::theme_cowplot() + xlab('Age (days from birth)') +
                         facet_wrap(vars(!!as.symbol(grouping))) +
                         theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
                         scale_colour_manual(values = rep(c(pals::alphabet() %>% unname()))))
      # draw plot
      temporal_human + temporal_mouse + plot_layout(ncol = 1)
    })

    output$temporal_plot <- renderPlot({
      temporal_plot()
    }, height = 700)

    ## dotplot ---------
    source('make_dotplot.R')
    draw_dotplot <- eventReactive(input$BUTTON_draw_dotplot,{make_dotplot(input, scEiaD_2020_v01, meta_filter)}
                                  )
    output$dotplot <- renderPlot({
      draw_dotplot()
    }, height = eventReactive(input$BUTTON_draw_dotplot, {input$dotplot_height %>% as.numeric()}))
  })

  # in situ ----
  # Functions used to generate in situ plots

  ## Function to read in original images
  get_image <- function(file) {
    image_read(file.path(paste0('www/insitu_layers/ret_',file, ".png")))
  }

  ## Function to recolor each individual cell layer based on expression
  recolor <- function(ret_layer, color){
    if (length(color) == 0) {
      recolored_layer <- ret_layer
    } else {
      recolored_layer <- ret_layer %>% image_colorize(75,color)
    }
    return(recolored_layer)
  }

  ## Function to gather gene expression data in table
  get_insitu_table <- function() {

    ### Pull the data for the gene of interest
    gene <- input$insitu_Gene
    grouping_features <- "CellType_predict"

    if (input$insitu_filter_cat !=''){
      validate(
        need(input$insitu_filter_on != '', "Please select at least one feature to filter on")
      )}

    ### Filter expression table, if filters active
    if (input$insitu_filter_cat != ''){
      filt_cat <- input$insitu_filter_cat
      filt_on <- input$insitu_filter_on
      full_table <- scEiaD_2020_v01 %>% tbl('grouped_stats') %>%
        filter(Gene == gene) %>%
        filter(!!as.symbol(filt_cat) %in% filt_on) %>%
        group_by_at(vars(one_of(c('Gene', grouping_features)))) %>%
        summarise(cpm = sum(cpm * cell_exp_ct) / sum(cell_exp_ct),
                  cell_exp_ct = sum(cell_exp_ct, na.rm = TRUE)) %>%
        as_tibble() %>%
        tidyr::drop_na() %>%
        full_join(., meta_filter %>%
                    group_by_at(vars(one_of(grouping_features))) %>%
                    summarise(Count = n())) %>%
        mutate(cell_exp_ct = ifelse(is.na(cell_exp_ct), 0, cell_exp_ct)) %>%
        mutate(`%` = round((cell_exp_ct / Count) * 100, 2),
               Expression = round(cpm * (`%` / 100), 2)) %>%
        select_at(vars(one_of(c('Gene', grouping_features, 'cell_exp_ct', 'Count', '%', 'Expression')))) %>%
        arrange(-Expression)
    }

    ### Or make expression table without filtering if none selected
    else {
      full_table <- scEiaD_2020_v01 %>% tbl('grouped_stats') %>%
        filter(Gene == gene) %>%
        group_by_at(vars(one_of(c('Gene', grouping_features)))) %>%
        summarise(cpm = sum(cpm * cell_exp_ct) / sum(cell_exp_ct),
                  cell_exp_ct = sum(cell_exp_ct, na.rm = TRUE)) %>%
        as_tibble() %>%
        tidyr::drop_na() %>%
        full_join(., meta_filter %>%
                    group_by_at(vars(one_of(grouping_features))) %>%
                    summarise(Count = n())) %>%
        mutate(cell_exp_ct = ifelse(is.na(cell_exp_ct), 0, cell_exp_ct)) %>%
        mutate(`%` = round((cell_exp_ct / Count) * 100, 2),
               Expression = round(cpm * (`%` / 100), 2)) %>%
        select_at(vars(one_of(c('Gene', grouping_features, 'cell_exp_ct', 'Count', '%', 'Expression')))) %>%
        arrange(-Expression)
    }
  }
  # Event reactives to produce image and table
  ## Reactive that generates data table
  insitu_table_maker <- eventReactive(input$BUTTON_draw_insitu, {

    full_table <- get_insitu_table()
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

    ### Load unedited images
    amacrine <- get_image('amacrine')
    artery <- get_image('artery')
    astrocyte <- get_image('astrocyte')
    axons <- get_image('axons')
    layer_labels <- get_image('background')
    bipolar <- get_image('bipolar')
    bruch <- get_image('bruch')
    choriocap <- get_image('choriocap')
    cones <- get_image('cones')
    horizontal <- get_image('horizontal')
    cell_labels <- get_image('labels')
    melanocytes <- get_image('melanocytes')
    microglia <- get_image('microglia')
    muller <- get_image('muller')
    rgc <- get_image('rgc')
    rods <- get_image('rods')
    rpe <- get_image('rpe')
    sclera <- get_image('sclera')
    vein <- get_image('vein')

    full_table <- get_insitu_table()

    ### Create mini table with color codes
    p <- full_table %>%
      select(CellType_predict,Expression) %>%
      tidyr::drop_na() %>%
      arrange(Expression)

    ### Convert expression to color scale
    p$col <- viridis(length(p$Expression))

    ### Generate legend plot
    leg_lab <- seq(min(p$Expression),max(p$Expression),l=5)
    leg_lab[] <- lapply(leg_lab, round,2)
    legend <- image_graph(width = 300, height = 600, res = 96)
    legend_image <- as.raster(matrix(rev(p$col), ncol=1))
    plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = expression(bold("Log"[2] * "(CPM + 1)")), cex.main=1.5)
    text(x=2, y = seq(0,1,l=5), labels = leg_lab[], cex=1.5)
    rasterImage(legend_image, 0, 0, 1,1)
    dev.off()

    ### Recolor each layer based on expression
    amacrine <- recolor(amacrine, p$col[which(p$CellType_predict == "Amacrine Cells")])
    artery <- recolor(artery,  p$col[which(p$CellType_predict == "Artery")])
    astrocyte <- recolor(astrocyte,  p$col[which(p$CellType_predict == "Astrocytes")])
    axons <- recolor(axons, p$col[which(p$CellType_predict == "Axons")])
    bipolar <- recolor(bipolar, p$col[which(p$CellType_predict == "Bipolar Cells")])
    bruch <- recolor(bruch, p$col[which(p$CellType_predict == "Bruch Membrane")])
    cones <- recolor(cones, p$col[which(p$CellType_predict == "Cones")])
    choriocap <- recolor(choriocap, p$col[which(p$CellType_predict == "Choriocapillaris")])
    horizontal <- recolor(horizontal, p$col[which(p$CellType_predict == "Horizontal Cells")])
    melanocytes <- recolor(melanocytes, p$col[which(p$CellType_predict == "Melanocytes")])
    microglia <- recolor(microglia, p$col[which(p$CellType_predict == "Microglia")])
    muller <- recolor(muller, p$col[which(p$CellType_predict == "Muller Glia")])
    rpe <- recolor(rpe, p$col[which(p$CellType_predict == "RPE")])
    rgc <- recolor(rgc, p$col[which(p$CellType_predict == "Retinal Ganglion Cells")])
    rods <- recolor(rods, p$col[which(p$CellType_predict == "Rods")])
    sclera <- recolor(sclera, p$col[which(p$CellType_predict == "Sclera")])
    vein <- recolor(vein,  p$col[which(p$CellType_predict == "Vein")])

    ### Merge the recolored layers into single image
    retina_insitu <- c(layer_labels, amacrine, artery, astrocyte,choriocap, bipolar, bruch, cones, horizontal, melanocytes, microglia, muller, axons, rpe, rgc, rods, sclera, vein, cell_labels)
    ret_img <- retina_insitu %>%
      image_mosaic() %>%
      image_flatten()

    ### Append the legend to the side and write a temporary file with the complete image
    tmpfile <- image_append(c(ret_img,legend), stack=FALSE) %>%
      image_write(tempfile(fileext='png'), format = 'png')

    ### Return the location of the temporary file
    return(list(src = tmpfile,
                height = input$insitu_height,
                contentType = "image/png"))
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
