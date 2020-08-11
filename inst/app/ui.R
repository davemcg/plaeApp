print('UI Start')
print(Sys.time())

library(shiny)
library(Cairo)
library(ggplot2)
library(scattermore)
library(pals)
library(shinythemes)
library(cowplot)

# header color
## orig: #2c3e50
## new:  #6439db
# header select color
## orig: #1a242f
## new: #400a91
# text color
## orig: #2c3e50
## new: #0f0f0f


# button and slider styling
# #3399ff

shinyUI(
  navbarPage('Ocular plae',id = 'nav',
             theme = 'flatly_mod.css',
             selected = 'Overview',
             navbarMenu('Viz', # UMAP ----------
                        tabPanel('UMAP - Tables',
                                 fluidPage(
                                   fluidRow(column(6, radioButtons(
                                     inputId = "gene_and_meta_scatter_tech",
                                     label = "Technology",
                                     choices = c("Droplet", "Well"),
                                     selected = "Droplet",
                                     inline = TRUE,
                                   ))),
                                   fluidRow(
                                     # Gene Scatter  ---------------
                                     column(6,
                                            plotOutput('gene_scatter_plot',
                                                       dblclick = "gene_scatter_plot_dblclick",
                                                       brush = brushOpts(
                                                         id = "gene_scatter_plot_brush",
                                                         resetOnNew = TRUE)),
                                            fluidRow(column(5,
                                                            selectizeInput('Gene', strong('Gene: '),
                                                                           choices=NULL, multiple=FALSE)),
                                                     column(5,
                                                            selectizeInput('pt_size_gene', strong('Point Size: '),
                                                                           choices=c(1,3,5,10),
                                                                           selected = 1, multiple=FALSE))),
                                            shinyWidgets::setSliderColor(c("#3399ff"), c(1)),
                                            fluidRow(column(5,
                                                            sliderInput("gene_scatter_slider", label = strong("Expression Range: "), min = 1,
                                                                        max = 15, value = c(1, 15))
                                            )),
                                            fluidRow(column(5,
                                                            selectizeInput('gene_filter_cat', strong('Filter Category: '),
                                                                           choices = NULL, selected = NULL, multiple = TRUE)),
                                                     column(5,
                                                            uiOutput('gene_filter_on_dynamicUI')))),
                                     # Meta Plot ------
                                     column(6,
                                            plotOutput('meta_plot',
                                                       dblclick = "meta_plot_dblclick",
                                                       brush = brushOpts(
                                                         id = "meta_plot_brush",
                                                         resetOnNew = TRUE)),
                                            fluidRow(column(5,
                                                            selectizeInput('meta_column', strong('Color: '),
                                                                           choices= NULL, selected = 'CellType_predict')),
                                                     column(5,
                                                            selectizeInput('pt_size_meta', strong('Point Size: '),
                                                                           choices=c(1,3,5),
                                                                           selected = 1, multiple=FALSE))),
                                            fluidRow(column(5,
                                                            selectInput("label_toggle", label = strong("Label: "),
                                                                        choices = list("None" = 0,
                                                                                       "CellType (published)" = 1,
                                                                                       "CellType (predict)" = 2,
                                                                                       "Cluster" = 3,
                                                                                       "Tabula Muris" = 4), multiple = FALSE,
                                                                        selected = 2)),
                                                     column(2,
                                                            radioButtons('meta_column_transform',
                                                                         label = 'Numeric Transform', inline = FALSE,
                                                                         choices = list("None" = "None", "log2" = "log2")))
                                            ),
                                            fluidRow(column(5,
                                                            selectizeInput('meta_filter_cat', strong('Filter Category: '),
                                                                           choices = NULL, selected = NULL, multiple = TRUE)),
                                                     column(5,
                                                            uiOutput('meta_filter_on_dynamicUI')))
                                            # selectizeInput('meta_filter_on', strong('Filter on: '),
                                            #                choices = NULL, selected = NULL, multiple = TRUE)))
                                     )
                                   ),
                                   fluidRow(
                                     column(6,
                                            fluidRow(
                                              column(12, actionButton('BUTTON_draw_scatter',' Draw Scatter Plot', icon = icon("arrow-up"),
                                                                      style='background-color: #3399ff; color: #ffffff'))),
                                            br(),
                                            selectizeInput('grouping_features', strong('Gene Table Grouping(s)'),
                                                           choices = NULL,
                                                           multiple = TRUE),
                                            fluidRow(
                                              column(12,  actionButton('BUTTON_make_gene_table',' Make Gene Table', icon = icon("arrow-down"),
                                                                       style='background-color: #3399ff; color: #ffffff'))
                                            ),
                                            br(),
                                            div(DT::dataTableOutput('gene_cluster_stats'), style='font-size:75%')),
                                     column(6,
                                            fluidRow(
                                              column(12,
                                                     actionButton('BUTTON_draw_meta',' Draw Meta Plot', icon = icon("arrow-up"),
                                                                  style='background-color: #3399ff; color: #ffffff'))),
                                            br(),
                                            selectizeInput('meta_groupings', strong('Metadata Table Groupings '),
                                                           choices = NULL,
                                                           multiple = TRUE),
                                            fluidRow(
                                              column(12,
                                                     actionButton('BUTTON_make_meta_table',' Make Meta Table', icon = icon("arrow-down"),
                                                                  style='background-color: #3399ff; color: #ffffff'))
                                            ),
                                            br(),
                                            div(DT::dataTableOutput('metadata_stats'), style='font-size:75%'))
                                   )
                                 )
                        ),
                        # exp_plots ------
                        tabPanel('Expression Plot',
                                 column(12,
                                        fluidRow(
                                          column(3,
                                                 (selectizeInput('exp_plot_genes', strong('Gene(s): '),
                                                                 choices = NULL, multiple = TRUE))),
                                          column(3,
                                                 selectizeInput('exp_plot_height', strong('Plot Height: '),
                                                                choices = seq(200, 4000, by = 100),
                                                                selected = 400, multiple = FALSE)),
                                          column(3,
                                                 selectInput('exp_plot_ylab', strong('Value: '),
                                                             choices = c('Mean CPM', '% of Cells Detected')))),
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
                                          #column(3, uiOutput('exp_filter_on_dynamicUI')),
                                          column(3, selectizeInput('exp_filter_on', strong('Filter On: '),
                                                                   choices = NULL, multiple = TRUE)),
                                        ),
                                        fluidRow(column(10, actionButton('BUTTON_draw_exp_plot','Draw Plot', icon = icon("arrow-down"),
                                                                         style='background-color: #3399ff; color: #ffffff'))),
                                        br(),
                                        fluidRow(column(10, plotOutput('exp_plot'))))),

                        # in situ ---------
                        tabPanel('In Situ Projection',
                                 fluidPage(
                                   column(8,
                                          fluidRow(
                                            column(4, selectizeInput('insitu_Gene', strong('Genes: '),
                                                                     choices=NULL, multiple=FALSE)),
                                            column(4, selectizeInput('insitu_height', strong('Plot Height: '),
                                                                     choices = seq(500, 1000, by = 100), selected = 700))
                                          ),
                                          fluidRow(
                                            column(4, selectizeInput('insitu_filter_cat', strong('Filter category: '),
                                                                     choices=NULL, multiple=FALSE)),
                                            column(4, selectizeInput('insitu_filter_on', strong('Filter on: '),
                                                                     choices=NULL, multiple=TRUE)),
                                          ),
                                          fluidRow(
                                            column(4,actionButton('BUTTON_draw_insitu','Draw In Situ Projection!', icon = icon("arrow-down"),
                                                                  style='background-color: #3399ff; color: #ffffff')),
                                            column(4,radioButtons('RADIO_show_insitu_table', "Show data table?", choices = c("Yes"="yes", "No"="no"), selected = "yes", inline=TRUE))
                                          ),

                                          fluidRow(
                                            plotOutput('insitu_img', height = "auto")
                                          ),
                                          conditionalPanel(condition = "input.RADIO_show_insitu_table == 'yes'",
                                                           hr(),
                                                           fluidRow(
                                                             div(DT::dataTableOutput('insitu_gene_stats'), style='font-size:75%'))
                                          )))),
                        tabPanel('Facet UMAP', # Facet UMAP ---------
                                 column(10,
                                        fluidRow(
                                          column(10,
                                                 fluidRow(column(5,
                                                                 selectizeInput('facet', strong('Facet On: '),
                                                                                choices=NULL, multiple=FALSE)),
                                                          column(5,
                                                                 selectizeInput('facet_color', strong('Color On: '),
                                                                                choices=NULL, multiple=FALSE)),
                                                          column(5,
                                                                 selectizeInput('pt_size_facet', strong('Point Size: '),
                                                                                choices=c(1,3,5,10),
                                                                                selected = 1, multiple=FALSE)),
                                                          column(5,
                                                                 selectizeInput('facet_height', strong('Plot Height: '),
                                                                                choices = c(100,200,300,400,600, 800),
                                                                                selected = 400, multiple = FALSE))),
                                                 fluidRow(column(10, actionButton('BUTTON_draw_filter','Draw Plot', icon = icon("arrow-down"),
                                                                                  style='background-color: #3399ff; color: #ffffff'))),
                                                 br(),
                                                 plotOutput('facet_plot'))
                                        )

                                 )),
                        # temporal plot -----
                        tabPanel('Temporal Gene x Cell Type',
                                 column(10,
                                        fluidRow(
                                          column(10,
                                                 fluidRow(column(3, selectizeInput('temporal_gene', strong('Gene(s): '),
                                                                                   choices = NULL, multiple = TRUE)),
                                                          column(3, selectInput('temporal_group', strong('Split on: '),
                                                                                choices = c('CellType', 'CellType (predict)'))),
                                                          column(3, selectInput('temporal_y_val', strong('Value: '),
                                                                                choices = c('Mean CPM', 'Ratio Detected'))),
                                                          column(3, selectInput('temporal_plot_height',strong('Plot Height: '),
                                                                                selected = 1000,
                                                                                choices = seq(400, 2000, 200)) )))),
                                        fluidRow(column(5,
                                                        actionButton('BUTTON_draw_temporal','Draw Plot', icon = icon("arrow-down"),
                                                                     style='background-color: #3399ff; color: #ffffff'))),
                                        br(), br(),
                                        fluidRow(column(10, plotOutput('temporal_plot')))
                                 )),
                        tabPanel('Dotplot', # Dotplot ---------
                                 column(8,
                                        fluidRow(
                                          column(4, selectizeInput('dotplot_Gene', strong('Genes: '),
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
                                                     style='background-color: #3399ff; color: #ffffff'),
                                        br(), br(),
                                        plotOutput('dotplot')))
             ),
             # # diff testing  tables ------------
             tabPanel('Diff Testing',
                      fluidPage(column(8,
                                       fluidRow(
                                         selectInput('search_by', strong('Search by: '),
                                                     choices = c('Gene',
                                                                 "CellType (Predict) against Remaining",
                                                                 "CellType against Remaining",
                                                                 "Cluster (Droplet) against Remaining",
                                                                 "Cluster (Well) against Remaining",
                                                                 "Organism against Organism within CellType",
                                                                 "Organism against Organism within CellType (Predict)",
                                                                 "Organism against Organism within Cluster (Droplet)",
                                                                 "Organism against Organism within Cluster (Well)",
                                                                 "Pairwise CellType (Predict) against CellType (Predict)",
                                                                 "Pairwise CellType against CellType",
                                                                 "Pairwise Cluster against Cluster",
                                                                 "Pairwise Cluster against Cluster (Well)"),
                                                     selected = 'Gene')
                                       )),
                                column(8,
                                       fluidRow(
                                         conditionalPanel("input.search_by == 'Gene'",
                                                          selectizeInput('diff_gene', strong('Genes: '),
                                                                         choices =  NULL,
                                                                         multiple = TRUE)),
                                         conditionalPanel("input.search_by != 'Gene'",
                                                          selectizeInput('diff_term', strong('Term: '),
                                                                         choices =  NULL,
                                                                         multiple = TRUE))
                                       )),
                                column(8,
                                       div(DT::dataTableOutput('make_diff_table'), style='font-size:75%')))),
             tabPanel('Overview', # Overview ------
                      fluidPage(
                        fluidRow(column(width = 8, offset = 1, h1('plae v0.35'))),
                        br(),
                        fluidRow(column(width = 8, offset = 1, h1('Overview'))),
                        fluidRow(column(width = 8, offset = 1, 'The light-sensitive portion of the mammalian eye is the retina. The retina itself is not a monolithic tissue - there are over 10 major cell types. The cones and rods which convert light into signal are supported by a wide variety of neural cell types with distinct roles in interpretting and transmitting the visual signal to the brain. Behind the retina is the RPE and vasculature, which supports the high energetic needs of the rods and cones. plae is a meta-analysis project over 1.2 million single-cell transcriptomes across 28 studies, 18 publications, and 3 species encompassing the back of the eye. Deep metadata minining, rigorous quality control analysis, differential gene expression testing, and deep learning based batch effect correction in a unified bioinformatic framework allow the universe of retina single cell expression information to be analyzed in one location.')),
                        fluidRow(column(width = 8, offset = 1, h1('Data Sources'))),
                        #fluidRow(column(width = 8, offset = 1, formattableOutput("formattable01"))),
                        fluidRow(column(width = 8, offset = 1, img(src='01_table.png', width="700px"))),
                        fluidRow(column(width = 8, offset = 1, h1('Extracted Cell Type Labels'))),
                        #fluidRow(column(width = 6, offset = 1, formattableOutput("formattable02"))),
                        fluidRow(column(width = 8, offset = 1, img(src='02_table.png', width="450px"))),
                        fluidRow(column(width = 8, offset = 1, 'Labelled cell types from published papers were pulled, where possible, from a combination of the Sequence Read Archive (SRA), lab web sites, and personal correspondence, then adjusted to be consistent (e.g. MG to Muller Glia) between all studies.'))),
                      br(),
                      fluidRow(column(width = 8, offset = 1, h1('Change log'))),
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
                      br(), br(), br()
             ))
)

