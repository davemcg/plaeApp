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
                                                                                   choices=c(1,2,3,6,10),
                                                                                   selected = 1, multiple=FALSE))),
                                                    shinyWidgets::setSliderColor(c("#3269FF"), c(1)),
                                                    fluidRow(column(5,
                                                                    sliderInput("gene_scatter_slider", label = strong("Filter Gene Expression (log2(counts + 1)): "), min =1,
                                                                                max = 10, value = c(2, 10))
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
                                                                                   choices=c(1,2,3,6,10),
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
                                                    # selectizeInput('meta_filter_on', strong('Select: '),
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
                                                             choices = c('Mean log2(Counts + 1)', '% of Cells Detected')))
                                        ),
                                        fluidRow(
                                          column(3,
                                                 (selectizeInput('exp_plot_facet', strong('Facet on: '),
                                                                 choices = c('CellType','cluster','CellType_predict'),
                                                                 selected = 'CellType_predict',
                                                                 multiple = FALSE))),
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
                                          column(3, selectizeInput('exp_filter_on', strong('Select: '),
                                                                   choices = NULL, multiple = TRUE)),
                                          column(3, numericInput('exp_filter_min_cell_number', strong('Minimum # Cells in Group: '),
                                                                 value = 50))
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
                                            column(5, selectizeInput('insitu_filter_on', label = strong('Select: '),
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
                                                                 selectizeInput('facet_filter_on', strong('Select: '),
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
                                          column(4, selectizeInput('dotplot_filter_on', strong('Select: '),
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
             # diff testing  tables ------------
             navbarMenu('Diff Testing',
                        tabPanel('Markers',
                                 fluidPage(column(8,
                                                  fluidRow(tags$h1('Gene Differential Expression Tests')),
                                                  fluidRow('Tests run with ', tags$a(href="https://bioconductor.org/packages/release/bioc/html/scran.html", "scran"), ' findMarkers tool.
                                                  The findMarkers tool uses the ', tags$a(href="https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test", "wilcox"), ' test to generate the
                                                  AUC scores and p values (the wilcox test is a bit more conservative than the ', tags$a(href="https://en.wikipedia.org/wiki/Student%27s_t-test", "t-test"), ' and more robust when outliers are present).
                                                           The log Fold Change (logFC) scores are generated by findMarkers t-test method.'),
                                                  br(),
                                                  fluidRow('What is "AUC"? Area Under the Curve. This reports the power of the gene of interest to distinguish between the base group
                                                           (e.g. Rods) versus the comparison group (e.g. Cones). An AUC of 1 mean thats the marker can perfectly (100%) distinguish cells
                                                           between the two groups with the marker. AUC of 0 (or missing) means that the gene has no power.'),
                                                  br(),
                                                  br(),
                                                  fluidRow(
                                                    selectInput('search_by', strong('Search by: '),
                                                                choices = c('Gene',
                                                                            "CellType (Predict)",
                                                                            "CellType",
                                                                            "Cluster"),
                                                                selected = 'Gene')
                                                  )),
                                           column(8,
                                                  fluidRow(
                                                    conditionalPanel("input.search_by == 'Gene' || input.search_by == 'Haystack'",
                                                                     selectizeInput('diff_gene', strong('Genes: '),
                                                                                    choices =  NULL,
                                                                                    multiple = FALSE)),
                                                    conditionalPanel("input.search_by != 'Gene' & input.search_by != 'Haystack'",
                                                                     selectizeInput('diff_base', strong('Base: '),
                                                                                    choices =  NULL,
                                                                                    multiple = FALSE)),
                                                    conditionalPanel("input.search_by != 'Gene' & input.search_by != 'Haystack'",
                                                                     selectizeInput('diff_against', strong('Against: '),
                                                                                    choices =  NULL,
                                                                                    multiple = FALSE))
                                                  )),
                                           column(12,
                                                  fluidRow(
                                                    fluidRow(
                                                      column(width = 6, div(DT::DTOutput('make_diff_table'), style='font-size:75%')),
                                                      column(width = 6, div(DT::DTOutput('make_diff_table_auc'), style='font-size:75%'))),
                                                    downloadButton("diff_table_download","Download all results as csv"),
                                                    br(), br())),
                                           fluidRow(column(12,
                                                           actionButton("diff_testing_help", "Page Pop Up Info"),
                                                           actionButton("data_table_help3", "Data Table Pop Up Info")))),
                                 linebreaks(10),
                                 fluidRow(includeHTML("www/footer.html"))),
                        # haystack tables ------------
                        tabPanel("singleCellHaystack",
                                 fluidPage(column(8, offset = 0,
                                                  fluidRow(tags$h1('Haystack')),
                                                  fluidRow(tags$a(href="https://www.nature.com/articles/s41467-020-17900-3", "singleCellHaystack"), ' is a
                                                           cluster or cell type independent method for identifying differentially expressed or "interesting" genes'),
                                                  br(),
                                                  fluidRow('Very briefly, it uses the ', tags$a(href="https://en.wikipedia.org/wiki/Kullback–Leibler_divergence",
                                                                                                HTML(paste0("D",tags$sub("KL")))),
                                                           ' divergence across the scVI multidimensional space to find non-randomly expressed genes. The table is ordered
                                                           by the log10(p value) calculated (lower is a lower p value). A higher D_KL score means that the genes is less randomly
                                                           expressed. T counts sums the number of counts (higher is expressed in more cells).'),
                                                  br(),
                                                  fluidRow('The CellType(s) and Cluster columns are the "top" genes which are most differentially expressed in the
                                                           comparison. The idea is to provide a quick way to see what CellType(s) or Cluster are driving the
                                                           singleCellHaystack identified gene.')),
                                           br(),
                                           column(10, offset = 0,
                                                  fluidRow(column(width = 10, div(DT::DTOutput('make_haystack_table'), style='font-size:75%'))))),
                                 linebreaks(10),
                                 fluidRow(includeHTML("www/footer.html"))
                        )
             ),

             tabPanel('Data', # Data ---------
                      fluidRow(column(width = 8, offset = 1, h1("Data"))),
                      br(),
                      fluidRow(column(width = 8, offset = 1, "The codebase for the creation of the scEiaD dataset is on ", tags$a(href="https://www.github.com/davemcg/scEiaD", "github"))),
                      br(),
                      fluidRow(column(width = 8, offset = 1, "If size not given, it is less than 1 GB")),
                      fluidRow(column(width = 8, offset = 1, h2("Run plae Locally"))),
                      fluidRow(column(width = 8, offset = 1, 'If you have 500GB (!) of free hard drive space, you can run plae on your own computer. ', tags$a(href="https://www.github.com/davemcg/plaeApp", "Installation instructions are available in our Github repository"), ' (this is the codebase for the app you are using now).')),
                      fluidRow(column(width = 8, offset = 1, h2("Seurat Objects"))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2021_10_08/scEiaD_all_seurat_v3.Rdata", "All (~25 GB)")))),
                      fluidRow(column(width = 8, offset = 1, h2("AnnData Objects"))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2021_10_08/scEiaD_all_anndata.h5ad", "All (~8 GB)")))),
                      fluidRow(column(width = 8, offset = 1, h2("Diff Testing Results"))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2021_10_08/wilcox_diff_results.tsv.gz", "All Diff Results (~1.1GB)")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2021_10_08/wilcox_diff_resultsCellType.tsv.gz", "CellType (Predict) against CellType testing")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2021_10_08/wilcox_diff_resultsCellType(Predict).tsv.gz", "CellType against CellType (Predict)")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2021_10_08/wilcox_diff_resultsCluster.tsv.gz", "Cluster against Cluster")))),
                      fluidRow(column(width = 8, offset = 1, h2("Metadata"))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2021_10_08/metadata_filter.tsv.gz", "Cell Metadata")))),
                      fluidRow(column(width = 8, offset = 1, h2("Counts"))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2021_10_08/counts.Rdata", "Kallisto counts, R sparse matrix (12 GB)")))),
                      fluidRow(column(width = 8, offset = 1, tags$li(tags$a(href="http://hpc.nih.gov/~mcgaugheyd/scEiaD/2021_10_08/counts_unfiltered.Rdata", "Kallisto counts, no filtering, R sparse matrix (15 GB)")))),
                      br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),
                      fluidRow(includeHTML("www/footer.html"))
             ),

             navbarMenu('Info', # Info ------
                        tabPanel('Overview', # Overview ------
                                 fluidPage(
                                   fluidRow(column(width = 8, offset = 1, h1('plae v0.80'))),
                                   br(),
                                   fluidRow(column(width = 8, offset = 1, h2(HTML("<b>PL</b>atform for <b>A</b>nalysis of sc<b>E</b>iad")))),
                                   fluidRow(column(width = 8, offset = 1,
                                                   img(src = 'plae_slide.png', width = "600px", alt='plae pronounced play logo. eye ball with arms running to slide made of retina cells')
                                   )),
                                   br(),
                                   fluidRow(column(width = 8, offset = 1, h1('What is scEiaD?'))),
                                   fluidRow(column(width = 8, offset = 1, HTML("<b>s</b>ingle <b>c</b>ell <b>E</b>ye <b>i</b>n <b>a</b> <b>D</b>isk"),)),
                                   br(),
                                   fluidRow(column(width = 8, offset = 1, 'The light-sensitive portion of the eye is the retina. The retina itself is not a monolithic tissue - there are over 10 major cell types. The cones and rods which convert light into signal are supported by a wide variety of neural cell types with distinct roles in interpretting and transmitting the visual signal to the brain. Behind the retina is the RPE and vasculature, which supports the high energetic needs of the rods and cones. scEiaD is a meta-atlas that compiles 1.1 million single-cell eye and body tissue transcriptomes across 42 studies, 31 publications, and 4 species. Deep metadata mining, rigorous quality control analysis, differential gene expression testing, and deep learning based batch effect correction in a unified bioinformatic framework allow the universe of retina single cell expression information to be analyzed in one location.')),
                                   br(),
                                   fluidRow(column(width = 8, offset = 1, h1('tldr'))),
                                   fluidRow(column(width = 8, offset = 1, 'You can look up gene expression by retina cell type across loads of different studies, four organisms, and multiple developmental stages.')),
                                   br(),
                                   fluidRow(column(width = 8, offset = 1, h1('Preprint of the data creation and benchmarking now on', tags$a(href="https://www.biorxiv.org/content/10.1101/2021.03.26.437190v1", "bioRxiv!")))),
                                   br(),
                                   fluidRow(column(width = 8, offset = 1, h1('Licensing'))),
                                   fluidRow(column(width = 8, offset = 1, 'This work is released under the CC0 license')),
                                   br(),
                                   fluidRow(column(width = 8, offset = 1, h1('Data Sources'))),

                                   fluidRow(column(width = 8, offset = 1, includeHTML("www/table_01.html"))),
                                   fluidRow(column(width = 8, offset = 1, h1('scEiaD Curated Published Cell Type Labels'))),
                                   fluidRow(column(width = 5, offset = 1, includeHTML("www/table_02.html"))),
                                   fluidRow(column(width = 8, offset = 1, 'Labelled cell types from published papers were pulled, where possible, from a combination of the Sequence Read Archive (SRA), lab web sites, and personal correspondence, then adjusted to be consistent (e.g. MG to Muller Glia) between all studies.')),
                                   br(),br(),
                                   fluidRow(column(width = 8, offset = 1, h1('scEiaD Machine Learned Cell Type Labels'))),
                                   fluidRow(column(width = 5, offset = 1, includeHTML("www/table_03.html"))),
                                   fluidRow(column(width = 8, offset = 1, 'The labels above were used to create a machine learning modeled which was used to relabel all* cells in the scEiaD (*above a confidence threshold of 0.5).'))),
                                 br(),br(),
                                 fluidRow(includeHTML("www/footer.html"))),
                        tabPanel('Analysis and Extension', # Analyses ------
                                 fluidPage(
                                   fluidRow(column(width = 8, offset = 1, (tags$h1('Using and Extending plae and the scEiaD')))),
                                   br(),
                                   fluidRow(column(width = 8, offset = 1, tags$b('All Links are External'))),
                                   br(),
                                   fluidRow(column(width = 8, offset = 1, includeHTML("www/analysis_toc.html"))))),
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
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.80 (2021-10-08): MASSIVE update. Chicken data added. Counts cleaned up with ',
                                 tags$a(href="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1950-6", "DecontX"),
                                 ' to remove (mostly) Rod gene contamination (e.g. Rho *was* everywhere).',  tags$a(href="https://www.nature.com/articles/s41467-020-17900-3", "singleCellHaystack"),
                                 ' table added to Diff Testing. Updated cell filtering with a higher minimum gene count cutoff to improve overall quality. Fixed bug in gene selection where human
                                 genes that mapped to multiple mouse genes were accidently removed. Improved scran-based differential gene expression testing with better parameters and added logFC calculations
                                 to improve interpretability. Fixed bug in dotplot plot where filtered data had the incorrect denominator values.')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.74 (2021-08-12): Remove broken link, fix bug in Expression Plot that was making average expression far too low.')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.73 (2021-06-10): Allow for UMAP plot to show all values (filter now starts at >=1)')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.72 (2021-04-29): Added more content and a table organzation to the Info -> Analysis... section')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.71 (2021-04-14): scEiaD preprint on bioRxiv! Added a filter in the "exp plot" section to remove data points with user-selected (default 50) minimum cells. I was finding that data points (e.g. cell type - study) with low N would often have "outlier" results. Added a new section to the web page - Analyses (under the "Info" tab)!')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.70 (2021-03-22): New scEiaD built with corrected fastq file sets (potential bug in 10x bamtofastq tool resulted ina few datasets getting scrambled barcodes). Removed a macaque dataset (SRR7733526) with odd behavior (clustering in 2D UMAP space largely alone). Tweaked the UMAP 2D gene view with a darker "background" cell color scheme to reduce "over-emphasis" on cells with low expression of a gene. CPM replaced with counts as some odd behaviour was detected in some genes in the UMAP view where there was high "background" expression. Counts have more consistent behavior. Removed hard filter that tossed cells with >2500 detected genes.')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.60 (2021-02-08): New scEiaD built with more studies. Removed several retinal organoid datasets that had snuck in. Added a filter option for the diff searching to search, for example, one cluster directly against another cluster.')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1)),
                                 fluidRow(column(width = 8, offset = 1, '0.52 (2021-01-13): Adding "missing" genes (we had only retained genes which were expressed in all three species, which naturally led to many genes (some important, like OPN1MW) to be dropped. That has been fixed. Tweak "Expression Plot" dot size to prevent crazy tiny point sizes. ')),
                                 br(),
                                 fluidRow(column(width = 8, offset = 1, '0.51 (2021-01-06): Hello 2021! ', tags$a(href="https://adamgayoso.com", "Adam Gayoso"), ' kindly pointed out that I was using scVI in a non-optimal manner, so I updated the scVI modeling to match their recommend "scArches" parameters. This (fortunately for my sanity) only subtly changes the downstream result. The more significant change is that we have totally changed the diff testing section to use the scran findMarkers test instead of our complicated and compute expensive pseudo-Bulk testing which continually gave odd results.')),
                                 br(),
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

