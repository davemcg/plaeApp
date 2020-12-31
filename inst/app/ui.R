print('UI Start')
print(Sys.time())

library(shiny)
library(magrittr)
library(Cairo)
library(ggplot2)
library(scattermore)
library(pals)
#library(shinythemes)
library(cowplot)
library(shinyalert)
library(magrittr)
library(DT)
# header color
## orig: #2c3e50
## new:  #6633ff
# header select color
## orig: #1a242f
## new: #400a91
# text color
## orig: #2c3e50
## new: #0f0f0f


# button and slider styling
# #3269FF

linebreaks <- function(n){HTML(strrep(br(), n))}

shinyUI(
  navbarPage('plae',id = 'nav',
             theme = 'flatly_mod.css',
             selected = 'Overview',
             navbarMenu('Viz', # UMAP ----------

                        tabPanel('UMAP - Tables',
                                 fluidPage(tags$html(lang="en"),
                                           tags$head(
                                             tags$style(HTML("
                                                    .shiny-output-error-validation {
                                                      color: #242526;
                                                      font-weight: bold;
                                                    }
                                             "))
                                           ),
                                           fluidRow(column(6, tags$h1("UMAP"))),
                                           fluidRow(
                                             # Gene Scatter  ---------------
                                             column(6,
                                                    plotOutput('gene_scatter_plot',
                                                               height = '500px',
                                                               dblclick = "gene_scatter_plot_dblclick",
                                                               brush = brushOpts(
                                                                 id = "gene_scatter_plot_brush",
                                                                 resetOnNew = TRUE)) %>%
                                                      shinycssloaders::withSpinner(type = 3, size = 0.5,
                                                                                   color = "#3269FF",
                                                                                   color.background = 'white'),
                                                    fluidRow(column(5,
                                                                    selectizeInput('Gene', strong('Scatter Gene: '),
                                                                                   choices=NULL, multiple=FALSE)),
                                                             column(5,
                                                                    selectizeInput('pt_size_gene', strong('Scatter Point Size: '),
                                                                                   choices=c(1,3,5,10),
                                                                                   selected = 1, multiple=FALSE))),
                                                    shinyWidgets::setSliderColor(c("#3269FF"), c(1)),
                                                    fluidRow(column(5,
                                                                    sliderInput("gene_scatter_slider", label = strong("Filter Gene Expression (log2(cpm + 1)): "), min = 1,
                                                                                max = 15, value = c(1, 15))
                                                    )),
                                                    fluidRow(column(5,
                                                                    selectizeInput('gene_filter_cat', label = strong('Scatter Filter Category: '),
                                                                                   choices = NULL, selected = NULL, multiple = TRUE)),
                                                             column(5,
                                                                    uiOutput('gene_filter_on_dynamicUI'))),
                                                    br(),
                                                    fluidRow(
                                                      column(5, actionButton('BUTTON_draw_scatter','
                                                                      Draw Scatter Plot', icon = icon("arrow-up"), alt =  'Button scatter up',
                                                                             style='background-color: #3269FF; color: #ffffff')),
                                                      column(5, downloadButton('BUTTON_download_scatter', 'Download Scatter Plot', alt = 'Download Scatter Plot to PNG in your browser default folder',
                                                                               style='background-color: #3269FF; color: #ffffff')))),
                                             # Meta Plot ------
                                             column(6,
                                                    plotOutput('meta_plot',
                                                               height = '500px',
                                                               dblclick = "meta_plot_dblclick",
                                                               brush = brushOpts(
                                                                 id = "meta_plot_brush",
                                                                 resetOnNew = TRUE)) %>%
                                                      shinycssloaders::withSpinner(type = 3,
                                                                                   size = 0.5,
                                                                                   color = "#3269FF",
                                                                                   color.background = 'white'),
                                                    fluidRow(column(5,
                                                                    selectizeInput('meta_column', strong('Meta Color: '),
                                                                                   choices= NULL, selected = 'CellType_predict')),
                                                             column(5,
                                                                    selectizeInput('pt_size_meta', strong('Meta Point Size: '),
                                                                                   choices=c(1,3,5),
                                                                                   selected = 1, multiple=FALSE))),
                                                    fluidRow(column(5,
                                                                    selectInput("label_toggle", label = strong("Meta Label: "),
                                                                                choices = list("None" = 0,
                                                                                               "CellType (published)" = 1,
                                                                                               "CellType (predict)" = 2,
                                                                                               "Cluster" = 3,
                                                                                               "Tabula Muris" = 4), multiple = FALSE,
                                                                                selected = 2)),
                                                             column(2,
                                                                    radioButtons('meta_column_transform',
                                                                                 label = 'Meta Numeric Transform', inline = FALSE,
                                                                                 choices = list("None" = "None", "log2" = "log2"))),
                                                             column(2, actionButton('BUTTON_show_meta_legend', 'Meta Legend', alt = 'Makes a pop up legend for the plot - will be blank is no legend for the particular plot', style='background-color: #3269FF; color: #ffffff')),

                                                    ),
                                                    fluidRow(column(5,
                                                                    selectizeInput('meta_filter_cat', strong('Meta Filter Category: '),
                                                                                   choices = NULL, selected = NULL, multiple = TRUE)),
                                                             column(5,
                                                                    uiOutput('meta_filter_on_dynamicUI'))),
                                                    fluidRow(
                                                      column(5,
                                                             actionButton('BUTTON_draw_meta',' Draw Meta Plot', icon = icon("arrow-up"), alt = 'meta plot displays at top of page',
                                                                          style='background-color: #3269FF; color: #ffffff')),
                                                      column(5, downloadButton('BUTTON_download_meta','Download Meta Plot', alt = 'Download Meta Plot to PNG in your browser default folder',
                                                                               style='background-color: #3269FF; color: #ffffff'))),
                                                    br()
                                                    # selectizeInput('meta_filter_on', strong('Filter on: '),
                                                    #                choices = NULL, selected = NULL, multiple = TRUE)))
                                             )
                                           ),
                                           fluidRow(column(6, tags$h1("Tables"))),
                                           # tags$head(
                                           #   tags$style(HTML("hr {border-top: 1px dashed #0f0f0f;}"))
                                           # ),
                                           # fluidRow(tags$hr()),
                                           fluidRow(
                                             column(6,
                                                    selectizeInput('grouping_features', strong('Gene Table Grouping(s)'),
                                                                   choices = NULL,
                                                                   multiple = TRUE),
                                                    actionButton('BUTTON_make_gene_table',' Make Gene Table', icon = icon("arrow-down"), alt = 'gene tables displays at top of page',
                                                                 style='background-color: #3269FF; color: #ffffff'),
                                                    br(), br(),
                                                    div(DT::dataTableOutput('gene_cluster_stats'), style='font-size:75%')),
                                             column(6,
                                                    selectizeInput('meta_groupings', strong('Metadata Table Groupings '),
                                                                   choices = NULL,
                                                                   multiple = TRUE),
                                                    actionButton('BUTTON_make_meta_table',' Make Meta Table', icon = icon("arrow-down"),alt = 'meta tables displays at top of page',
                                                                 style='background-color: #3269FF; color: #ffffff'),
                                                    br(), br(),
                                                    div(DT::dataTableOutput('metadata_stats'), style='font-size:75%'))
                                           ),
                                           fluidRow(column(6,
                                                           actionButton("umap_table_help", "Page Pop Up Info"),
                                                           actionButton("data_table_help1", "Data Table Pop Up Info")))),
                                 linebreaks(8),
                                 fluidRow(includeHTML("www/footer.html"))
                        ),
                        # exp_plots ------
                        tabPanel('Expression Plot',
                                 column(12,
                                        fluidRow(column(6, tags$h1('Expression Plot'))),
                                        fluidRow(
                                          column(4,
                                                 (selectizeInput('exp_plot_genes', strong('Gene(s): '),
                                                                 choices = NULL, multiple = TRUE))),
                                          column(2,
                                                 selectizeInput('exp_plot_height', strong('Plot Height: '),
                                                                choices = seq(200, 4000, by = 100),
                                                                selected = 400, multiple = FALSE)),
                                          column(3,
                                                 selectInput('exp_plot_ylab', strong('Value: '),
                                                             choices = c('Mean CPM', '% of Cells Detected')))
                                        ),
                                        fluidRow(
                                          column(3,
                                                 (selectizeInput('exp_plot_facet', strong('Facet on: '),
                                                                 choices = c('CellType','cluster','CellType_predict'), multiple = FALSE))),
                                          column(3,
                                                 (selectizeInput('exp_plot_groups', strong('Color on: '),
                                                                 choices = NULL, multiple = FALSE))),
                                          column(3,
                                                 numericInput('exp_plot_col_num', strong('Number of columns: '),
                                                              min = 1, max = 30, value = 5)),
                                        ),
                                        fluidRow(
                                          column(3, selectizeInput('exp_filter_cat', strong('Filter Category: '),
                                                                   choices = NULL, multiple = TRUE)),
                                          column(3, selectizeInput('exp_filter_on', strong('Filter On: '),
                                                                   choices = NULL, multiple = TRUE)),
                                        ),
                                        fluidRow(column(3, actionButton('BUTTON_draw_exp_plot','Draw Plot', icon = icon("arrow-down"),
                                                                        alt = 'BUTTON draw exp plot below',
                                                                        style='background-color: #3269FF; color: #ffffff')),
                                                 column(3, downloadButton('BUTTON_download_exp','
                                                                      Download Scatter Plot', alt = 'BUTTON download exp plot',
                                                                          style='background-color: #3269FF; color: #ffffff'))),
                                        br(),
                                        fluidRow(column(10, plotOutput('exp_plot') %>% shinycssloaders::withSpinner(type = 3, size = 0.5, color = "#3269FF", color.background = 'white'))),
                                        br(),
                                        actionButton("exp_plot_help", "Page Pop Up Info"),),
                                 linebreaks(160),
                                 fluidRow(includeHTML("www/footer.html"))),

                        # in situ ---------
                        tabPanel('In Situ Projection',
                                 fluidPage(
                                   column(8,
                                          fluidRow(column(6, tags$h1('In Situ Projection'))),
                                          fluidRow(
                                            column(5, selectizeInput('insitu_Gene', strong('Genes: '),
                                                                     choices=NULL, multiple=FALSE)),
                                            column(5, selectizeInput('insitu_height', strong('Plot Height: '),
                                                                     choices = seq(500, 1000, by = 100), selected = 700))
                                          ),
                                          fluidRow(
                                            column(5, selectizeInput('insitu_filter_cat', label = strong('Filter category: '),
                                                                     choices=NULL, multiple=FALSE)),
                                            column(5, selectizeInput('insitu_filter_on', label = strong('Filter on: '),
                                                                     choices=NULL, multiple=TRUE)),
                                          ),
                                          fluidRow(
                                            column(5,actionButton('BUTTON_draw_insitu','Draw In Situ Projection!', icon = icon("arrow-down"),
                                                                  alt = 'button draw insitu below',
                                                                  style='background-color: #3269FF; color: #ffffff')),
                                            column(5,radioButtons('RADIO_show_insitu_table', "Show data table?",
                                                                  choices = c("Yes"="yes", "No"="no"), selected = "yes", inline=TRUE))
                                          ),
                                          fluidRow(
                                            plotOutput('insitu_img', height = "auto")
                                          ),
                                          conditionalPanel(condition = "input.RADIO_show_insitu_table == 'yes'",
                                                           br(),
                                                           fluidRow(
                                                             div(DT::dataTableOutput('insitu_gene_stats'), style='font-size:75%'))
                                          ),
                                          actionButton("insitu_help", "Page Pop Up Info"),
                                          actionButton("data_table_help2", "Data Table Pop Up Info")
                                   )),
                                 linebreaks(120),
                                 fluidRow(includeHTML("www/footer.html"))),
                        tabPanel('Facet UMAP', # Facet UMAP ---------

                                 column(10,
                                        fluidRow(column(6, tags$h1('Facet UMAP'))),
                                        fluidRow(
                                          column(10,
                                                 fluidRow(column(5,
                                                                 selectizeInput('facet', strong('Facet On: '),
                                                                                choices=NULL, multiple=FALSE)),
                                                          column(5,
                                                                 selectizeInput('facet_color', strong('Color On: '),
                                                                                choices=NULL, multiple=FALSE)),
                                                          column(5,
                                                                 selectizeInput('facet_filter_cat', strong('Filter Category: '),
                                                                                choices = NULL, multiple=TRUE)),
                                                          column(5,
                                                                 selectizeInput('facet_filter_on', strong('Filter On: '),
                                                                                choices = NULL, multiple=TRUE)),
                                                          column(5,
                                                                 selectizeInput('pt_size_facet', strong('Point Size: '),
                                                                                choices=c(1,3,5,10),
                                                                                selected = 1, multiple=FALSE)),
                                                          column(5,
                                                                 selectizeInput('facet_height', strong('Plot Height: '),
                                                                                choices = c(100,200,300,400,600, 800),
                                                                                selected = 400, multiple = FALSE))),
                                                 fluidRow(column(5, actionButton('BUTTON_draw_filter','Draw Plot', icon = icon("arrow-down"),
                                                                                 alt = 'BUTTON draw plot below',
                                                                                 style='background-color: #3269FF; color: #ffffff')),
                                                          column(5, downloadButton('BUTTON_download_facet','Download Scatter Plot',
                                                                                   alt = 'download facet plot button',
                                                                                   style='background-color: #3269FF; color: #ffffff'))),
                                                 br(),
                                                 plotOutput('facet_plot') %>% shinycssloaders::withSpinner(type = 3, size = 0.5, color = "#3269FF", color.background = 'white'))
                                        ),
                                        br(),
                                        actionButton("facet_umap_help", "Page Pop Up Info")),
                                 linebreaks(72),
                                 fluidRow(includeHTML("www/footer.html"))),
                        # temporal plot -----
                        # tabPanel('Temporal Gene x Cell Type',
                        #          actionButton("temporal_plot_help", "Page Pop Up Info"),
                        #          column(10,
                        #                 fluidRow(
                        #                   column(10,
                        #                          fluidRow(column(3, selectizeInput('temporal_gene', strong('Gene(s): '),
                        #                                                            choices = NULL, multiple = TRUE)),
                        #                                   column(3, selectInput('temporal_group', strong('Split on: '),
                        #                                                         choices = c('CellType', 'CellType (predict)'))),
                        #                                   column(3, selectInput('temporal_y_val', strong('Value: '),
                        #                                                         choices = c('Mean CPM', 'Ratio Detected'))),
                        #                                   column(3, selectInput('temporal_plot_height',strong('Plot Height: '),
                        #                                                         selected = 1000,
                        #                                                         choices = seq(400, 2000, 200)) )))),
                        #                 fluidRow(column(5,
                        #                                 actionButton('BUTTON_draw_temporal','Draw Plot', icon = icon("arrow-down"),
                        #                                              style='background-color: #3269FF; color: #ffffff'))),
                        #                 br(), br(),
                        #                 fluidRow(column(10, plotOutput('temporal_plot')))
                        #          )),
                        tabPanel('Dotplot', # Dotplot ---------
                                 column(8,
                                        fluidRow(column(6, tags$h1('Dotplot'))),
                                        fluidRow(
                                          column(6, selectizeInput('dotplot_Gene', strong('Genes: '),
                                                                   choices=NULL, multiple=TRUE)),
                                          column(4, selectizeInput('dotplot_height', strong('Plot Height: '),
                                                                   choices = seq(400, 2000, by = 100), selected = 800))),
                                        fluidRow(
                                          column(4, selectizeInput('dotplot_groups', strong('Group by (two max): '),
                                                                   choices=NULL, multiple=TRUE)),
                                          column(4, selectizeInput('dotplot_filter_cat', strong('Filter category: '),
                                                                   choices=NULL, multiple=FALSE)),
                                          column(4, selectizeInput('dotplot_filter_on', strong('Filter on: '),
                                                                   choices=NULL, multiple=TRUE)),
                                        ),
                                        actionButton('BUTTON_draw_dotplot','Draw Dotplot!', icon = icon("arrow-down"),
                                                     alt = 'button draw dotplot below',
                                                     style='background-color: #3269FF; color: #ffffff'),
                                        downloadButton('BUTTON_download_dotplot','Download Dot Plot',
                                                       alt = 'button download dotplot',
                                                       style='background-color: #3269FF; color: #ffffff'),
                                        actionButton("dotplot_help", "Page Pop Up Info"),
                                        br(), br(),
                                        plotOutput('dotplot') %>% shinycssloaders::withSpinner(type = 3, size = 0.5, color = "#3269FF", color.background = 'white')),
                                 linebreaks(120),
                                 fluidRow(includeHTML("www/footer.html")))
             ),
             # # diff testing  tables ------------
             tabPanel('Diff Testing',
                      fluidPage(column(8,
                                       fluidRow(tags$h1('Pseudo Bulk Diff Testing')),
                                       fluidRow(
                                         selectInput('search_by', strong('Search by: '),
                                                     choices = c('Gene',
                                                                 "CellType (Predict) against Remaining",
                                                                 "CellType against Remaining",
                                                                 "Cluster against Remaining",
                                                                 "Organism against Organism within CellType",
                                                                 "Organism against Organism within CellType (Predict)",
                                                                 "Organism against Organism within Cluster",
                                                                 "Pairwise CellType (Predict) against CellType (Predict)",
                                                                 "Pairwise CellType against CellType",
                                                                 "Pairwise Cluster against Cluster"),
                                                     selected = 'Gene')
                                       )),
                                column(8,
                                       fluidRow(
                                         conditionalPanel("input.search_by == 'Gene'",
                                                          selectizeInput('diff_gene', strong('Genes: '),
                                                                         choices =  NULL,
                                                                         multiple = FALSE)),
                                         conditionalPanel("input.search_by != 'Gene'",
                                                          selectizeInput('diff_term', strong('Term: '),
                                                                         choices =  NULL,
                                                                         multiple = FALSE))
                                       )),
                                column(12,
                                       fluidRow(
                                         div(DT::DTOutput('make_diff_table'), style='font-size:75%'),
                                         downloadButton("diff_table_download","Download all results as csv"),
                                         br(), br())),
                                fluidRow(column(12,
                                                actionButton("diff_testing_help", "Page Pop Up Info"),
                                                actionButton("diff_testing_help2", "Pseudo Bulk?"),
                                                actionButton("data_table_help3", "Data Table Pop Up Info")))),
                      linebreaks(10),
                      fluidRow(includeHTML("www/footer.html"))),
             tabPanel('Data', # Data ---------
                      fluidRow(column(width = 8, offset = 1, h1("Data"))),
                      fluidRow(column(width = 8, offset = 1, "If size not given, it is less than 1 GB")),
                      fluidRow(column(width = 8, offset = 1, h2("Run plae Locally"))),
                      fluidRow(column(width = 8, offset = 1, 'If you have 200GB of free hard drive space, you can run plae on your own computer. ', tags$a(href="https://www.github.com/davemcg/plaeApp", "Installation instructions are available in our Github repository"), ' (this is the codebase for the app you are using now).')),
                      fluidRow(column(width = 8, offset = 1, h2("Seurat Objects"))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/scEiaD_droplet_seurat_v3.Rdata", "Droplet (~3.5 GB)")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/scEiaD_well_seurat_v3.Rdata", "Well")))),
                      fluidRow(column(width = 8, offset = 1, h2("AnnData Objects"))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/scEiaD_droplet_anndata.h5ad", "Droplet (~3.4 GB)")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/scEiaD_well_anndata.h5ad", "Well")))),
                      fluidRow(column(width = 8, offset = 1, h2("PseudoBulk Diff Results"))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/pseudoBulk_diff_results.tsv.gz", "All PseudoBulk Results (~4.5GB)")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/pseudoBulkTable_CellTypeagainstRemaining.tsv.gz", "CellType (Predict) against Remaining")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/pseudoBulkTable_CellType(Predict)againstRemaining.tsv.gz", "CellType against Remaining")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/pseudoBulkTable_Cluster(Droplet)againstRemaining.tsv.gz", "Cluster (Droplet) against Remaining")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/pseudoBulkTable_Cluster(Well)againstRemaining.tsv.gz", "Cluster (Well) against Remaining")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/pseudoBulkTable_OrganismagainstOrganismwithinCellType(Predict).tsv.gz", "Organism against Organism within CellType")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/pseudoBulkTable_OrganismagainstOrganismwithinCellType.tsv.gz", "Organism against Organism within CellType (Predict)")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/pseudoBulkTable_OrganismagainstOrganismwithinCluster(Droplet).tsv.gz", "Organism against Organism within Cluster (Droplet)")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/pseudoBulkTable_OrganismagainstOrganismwithinCluster(Well).tsv.gz", "Organism against Organism within Cluster (Well)")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/pseudoBulkTable_PairwiseCellTypeagainstCellType.tsv.gz", "Pairwise CellType (Predict) against CellType (Predict)")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/pseudoBulkTable_PairwiseCellType(Predict)againstCellType(Predict).tsv.gz", "Pairwise CellType against CellType")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/pseudoBulkTable_PairwiseClusteragainstCluster(Droplet).tsv.gz", "Pairwise Cluster against Cluster (Droplet) (~1.6 GB)")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/pseudoBulkTable_PairwiseClusteragainstCluster(Well).tsv.gz", "Pairwise Cluster against Cluster (Well)")))),
                      fluidRow(column(width = 8, offset = 1, h2("PseudoBulk Matrices"))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/pseudoBulk_celltype.tsv.gz", "Counts summed across CellType and batch")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/pseudoBulk_celltypePredict.tsv.gz", "Counts summed across CellType (predicted) and batch")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/pseudoBulk_clusterDroplet.tsv.gz", "Counts summed across Cluster (Droplet) and batch")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/pseudoBulk_clusterWell.tsv.gz", "Counts summed across Cluster (Well) and batch")))),
                      fluidRow(column(width = 8, offset = 1, h2("Metadata"))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/metadata_filter.tsv.gz", "Cell Metadata")))),
                      fluidRow(column(width = 8, offset = 1, h2("Counts"))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/counts.Rdata", "Kallisto counts, R sparse matrix (~1.9 GB)")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/cpm.Rdata", "Kallisto counts, cpm scaled, R sparse matrix (~1.9 GB)")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2020_10_19/counts_unfiltered.Rdata", "Kallisto counts, no filtering, R sparse matrix (~1.9 GB)")))),
                      br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),
                      fluidRow(includeHTML("www/footer.html"))
             ),

             navbarMenu('Info', # Info ------
                        tabPanel('Overview', # Overview ------
                                 fluidPage(
                                   fluidRow(column(width = 8, offset = 1, h1('plae v0.50'))),
                                   br(),
                                   fluidRow(column(width = 8, offset = 1, h2(HTML("<b>PL</b>atform for <b>A</b>nalysis of sc<b>E</b>iad")))),
                                   fluidRow(column(width = 8, offset = 1,
                                                   img(src = 'plae_slide.png', width = "600px", alt='plae pronounced play logo. eye ball with arms running to slide made of retina cells')
                                   )),
                                   br(),
                                   fluidRow(column(width = 8, offset = 1, h1('What is scEiaD?'))),
                                   fluidRow(column(width = 8, offset = 1, HTML("<b>s</b>ingle <b>c</b>ell <b>E</b>ye <b>i</b>n <b>a</b> <b>D</b>isk"),)),
                                   br(),
                                   fluidRow(column(width = 8, offset = 1, 'The light-sensitive portion of the eye is the retina. The retina itself is not a monolithic tissue - there are over 10 major cell types. The cones and rods which convert light into signal are supported by a wide variety of neural cell types with distinct roles in interpretting and transmitting the visual signal to the brain. Behind the retina is the RPE and vasculature, which supports the high energetic needs of the rods and cones. scEiaD is a meta-atlas that compiles 1.2 million single-cell back of the eye transcriptomes across 28 studies, 18 publications, and 3 species. Deep metadata mining, rigorous quality control analysis, differential gene expression testing, and deep learning based batch effect correction in a unified bioinformatic framework allow the universe of retina single cell expression information to be analyzed in one location.')),
                                   br(),
                                   fluidRow(column(width = 8, offset = 1, h1('tldr'))),
                                   fluidRow(column(width = 8, offset = 1, 'You can look up gene expression by retina cell type across loads of different studies, three organisms, and multiple developmental stages.')),
                                   br(),
                                   fluidRow(column(width = 8, offset = 1, h1('Data Sources'))),

                                   fluidRow(column(width = 8, offset = 1, includeHTML("www/table_01.html"))),
                                   fluidRow(column(width = 8, offset = 1, h1('scEiaD Curated Published Cell Type Labels'))),
                                   fluidRow(column(width = 5, offset = 1, includeHTML("www/table_02.html"))),
                                   fluidRow(column(width = 8, offset = 1, 'Labelled cell types from published papers were pulled, where possible, from a combination of the Sequence Read Archive (SRA), lab web sites, and personal correspondence, then adjusted to be consistent (e.g. MG to Muller Glia) between all studies.')),
                                   br(),br(),
                                   fluidRow(column(width = 8, offset = 1, h1('scEiaD Machine Learned Cell Type Labels'))),
                                   fluidRow(column(width = 5, offset = 1, includeHTML("www/table_02x.html"))),
                                   fluidRow(column(width = 8, offset = 1, 'The labels above were used to create a machine learning modeled which was used to relabel all* cells in the scEiaD ((above a confidence threshold).'))),
                                 br(),br(),
                                 fluidRow(includeHTML("www/footer.html"))),
                        tabPanel('Contact', # Contact ------
                                 fluidPage(
                                   fluidRow(column(width = 8, offset = 1, h1('Contact'))),
                                   fluidRow(column(width = 8, offset = 1, "If you have questions about scEiaD dataset or the plae application, please contact ", tags$a(href="https://www.nei.nih.gov/research/research-labs-and-branches/ophthalmic-genetics-and-visual-function-branch/bioinformatics-group", "David McGaughey, Ph.D"), " .")),
                                   br(),
                                   fluidRow(column(width = 8, offset = 1, "Otherwise the National Eye Institute's Office of Science Communications, Public Liaison and Education responds directly to requests for information on eye diseases and vision research in English and Spanish. We cannot provide personalized medical advice to individuals about their condition or treatment.")),
                                   br(),
                                   fluidRow(column(width = 8, offset = 1, tags$a(href="mailto:2020@nei.nih.gov",
                                                                                 "2020@nei.nih.gov"))),
                                   fluidRow(column(width = 8, offset = 1, "Phone: 301-496-5248 — English and Spanish")),
                                   fluidRow(column(width = 8, offset = 1, "Mail: National Eye Institute")),
                                   fluidRow(column(width = 8, offset = 1, "Information Office")),
                                   fluidRow(column(width = 8, offset = 1, "31 Center Drive MSC 2510")),
                                   fluidRow(column(width = 8, offset = 1, "Bethesda, MD 20892-2510")),
                                   linebreaks(40),
                                   fluidRow(includeHTML("www/footer.html")))),
                        tabPanel('Change Log', # Change Log ------
                                 fluidRow(column(width = 8, offset = 1, h1('Change log'))),
                                 fluidRow(column(width = 8, offset = 1, '0.50 (2020-12-31): Goodbye 2020! Major update to the scVI-based UMAP projection which improves data quality. Removed non-tissue samples (e.g. organoid/cell lines). They will be added back later once I figure out a logical/simple way to do it. Fixed major bug in QC filtering which failed to remove high mitochondrial count (likely apoptosing cells) cells. Dot plot tweaked to improve relative dot sizes. Cowan et al. dataset added.')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.43 (2020-11-09): Downloadable diff results added to "Data." The diff results reactive data table now has a "Download all ..." button which replaces the "CSV" button that only downloaded the viewable data (100 max).')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.42 (2020-10-16): Alt text added to each button, tweaked UMAP-Tables layout again. Slide logo added. Site went public at the version on 2020-11-02!')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.41 (2020-10-06): Download buttons added for each plot.')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.40 (2020-10-05): UI and text labels tweaked in UMAP-Tables to improve tab selection order. Dot plot given a bar plot to show category size. Error handling improved when user fails to provide a category value to filter on. Data Table help button added. Help buttons moved to bottom of page with consistent visual - tabbing order. Colors of UI elements tweaked to improve contrast.')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.39 (2020-09-18): Fixed calculation error in dotplot where expression not scaled by number of cells in grouping variable.')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.38 (2020-09-02): Contact section and footer added for compliance.')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.37 (2020-08-24): Help pop section populated with text. Put white halo back around text in UMAP - Meta section. Loading circles added to plots. Row names added to tables. Diff testing filtered to only return results with FDR < 0.05 and abs(logFC) > 0.5.')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.36 (2020-08-17): Data download section added. Change log moved to separate section. CSS tweaked to show links in blue. First overview table updated to improve contrast. UMAP plots axis fixed.')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.35 (2020-08-14): Moved Overview tables to html for improved rendering and switched over color-blind friendly palette. Temporarily removed Temporal Plotting section. Improved filtering for Facet Plot. DotPlot plotting fixed and improved. Back-end server.R code moved into separate functions. Colors fixed so they stay consistent when filtering/subsetting the plots. Site now starts from scratch in under 5 seconds with improved fst-based data loading and pre-calculating more operations.')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.34 (2020-08-03): Fixed issue with TabulaMuris labels not appearing. Scanned app with koa11y for 508 compliance - changed headers from h2 to h1 to comply. ')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.33 (2020-07-30): Exp plot can now take space or comma separated Genes as input. User can selected number of columns in Exp Plot. Diff Table formatting improved with rounding and PB_Test can be selected as a drop down now in the data table search.')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.32 (2020-07-29): In situ Projection viz added courtesy of Zachary Batz! It\'s a simulated cross section of the retina with each cell type colored by intensity of scRNA expression! Move table draw button under filtering in UMAP - Tables. Sort Diff Exp results by FDR. Filtering on numeric column now returns slider UI. Remove super dangerous ability to create faceted plots on numeric values.')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.31 (2020-07-24): Re-created scEiaD with better internal (Hufnagel) transwell RPE labelling (there are roughly two groups - mature RPE with high TTR expression and less (?) mature RPE with lower TTR), removal of the SRP166660 study as it was *all* non-normal (injured retina) (confirmed with correspondence with Dr. Poche), removed the pan RGC CellType labelling for the SRP212151 as I see post-hoc that there are LOADS of non-RGC cells. Did the same for SRP186407, which has substantial non-microglia. Generally, FACS != 100% celltype purity. Added differential testing against all Tabula Muris cell types. Removing clusters/cells with high doublet scores. Added cell cycle phase (G1/G2M/S) assignment. More study level metadata.')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.30 (2020-07-20): Huge update. Hundreds of thousands of cells added. The Tabula Muris project data (pan mouse) has been added to faciliate non-eye comparison. Filtering options added to most of the plotting views to allow for quick slicing into this huge dataset. Differential expression testing totally reworked - now uses "pseudoBulk" approach to better utilize the large number of studies we have.')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.23 (2020-06-16): Remove low N cell type from diff expression tables, tweak Overview with spacing alterations and updated text.')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.22 (2020-06-15): Added expression plot by user selected groups plot view. Fixed bug in mean cpm expression calculation for Viz -> UMAP - Table gene tables')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.21 (2020-06-15): Added subcluster diff testing tables, temporal gene expression by celltype plot section.')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.20 (2020-06-06): New 2D UMAP projection that includes the full Yu - Clark Human scRNA dataset. Added tables to "Overview" section showing data stats. Added "filtering" functionality to UMAP plot section.')),
                                 linebreaks(2),
                                 fluidRow(includeHTML("www/footer.html"))
                        ))
             #tags$head(tags$style(".helpAlign{float:right;}"))
  )
)

