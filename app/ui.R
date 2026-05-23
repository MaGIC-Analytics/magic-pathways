library(shiny)
require(shinyjs)
library(shinythemes)
require(shinycssloaders)
library(shinyWidgets)

library(DT)
library(tidyverse)
library(data.table)
library(colourpicker)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(DOSE)
library(AnnotationDbi)
library(msigdbr)
## UpSetR removed — enrichplot::upsetplot has signature issues

tagList(
    tags$head(
        includeHTML(("www/GA.html")),
        tags$style(type = 'text/css','.navbar-brand{display:none;}'),
        tags$style(HTML("
            .control-group-panel {
                border: 1px solid #ddd;
                border-radius: 6px;
                padding: 10px 12px;
                margin-bottom: 10px;
                background-color: #f9f9f9;
            }
            .control-group-title {
                font-weight: bold;
                font-size: 14px;
                color: #0F344C;
                margin-bottom: 8px;
            }
            #show_help_float {
                position: fixed;
                bottom: 28px;
                right: 28px;
                z-index: 9999;
                border-radius: 50%;
                width: 46px;
                height: 46px;
                font-size: 20px;
                padding: 0;
                box-shadow: 0 3px 8px rgba(0,0,0,0.25);
            }
        "))
    ),
    ## Global always-visible help button (fixed bottom-right)
    actionButton("show_help_float", label=NULL,
        icon=icon("circle-question"),
        title="Help & documentation",
        class="btn btn-info"
    ),
    fluidPage(theme = shinytheme('yeti'),
            windowTitle = "MaGIC Gene Set Enrichment Tool",
            useShinyjs(),
            titlePanel(
                fluidRow(
                column(2, tags$a(href='http://www.bioinformagic.io/', tags$img(height=75, src="MaGIC_Icon_0f344c.svg")), align='center'),
                column(10, fluidRow(
                    column(10, h1(strong('MaGIC Gene Set Enrichment Tool'), align='center', style="color:#0F344C;"))
                ))
                ),
                windowTitle = "MaGIC Gene Set Enrichment Tool"),
                tags$style(type='text/css', '.navbar{font-size:20px;}'),
                tags$style(type='text/css', '.nav-tabs{padding-bottom:20px;}'),
                tags$style(type='text/css', '.navbar-default{background-color:#0F344C;}'),
                tags$style(type='text/css', HTML('.navbar { background-color: #0F344C;}
                          .tab-panel{ background-color: #0F344C;}
                          .navbar-default .navbar-nav > .active > a,
                           .navbar-default .navbar-nav > .active > a:focus,
                           .navbar-default .navbar-nav > .active > a:hover {
                                color: white;
                                background-color: #008cba;
                            }')
                          ),
                tags$head(tags$style(".modal-dialog{ width:1300px}")),

        navbarPage(title="", id='NAVTABS',

        ## Intro Page
##########################################################################################################################################################
            tabPanel('Introduction',
                fluidRow(
                    column(2),
                    column(8,
                        column(12, align="center", style="margin-bottom:25px;",
                            h3(markdown("Welcome to the Gene Set Enrichment Tool by the
                            [Molecular and Genomics Informatics Core (MaGIC)](http://www.bioinformagic.io)."))),
                        hr(),
                        h4("How to Use This Tool", style="color:#0F344C;"),
                        tags$ol(
                            tags$li(strong("Navigate to the Data Input tab."),
                                " Upload your DESeq2 results table, or click 'Load Demo Data' to explore with a built-in synthetic dataset."),
                            tags$li(strong("Configure your analysis."),
                                " Select your organism, gene ID type, and analysis mode (ORA or GSEA)."),
                            tags$li(strong("Submit your data."),
                                " The Enrichment Analysis tab will become visible once data is successfully loaded."),
                            tags$li(strong("Select a pathway database and run the analysis."),
                                " Choose from GO, KEGG, Reactome, or MSigDB Hallmark, adjust parameters, and click 'Run Enrichment Analysis'."),
                            tags$li(strong("Explore and download results."),
                                " View dot plots, bar plots, enrichment maps, gene-concept networks, and more. Download results as CSV or plots as PDF/PNG.")
                        ),
                        hr(),
                        h4("Over-Representation Analysis (ORA)", style="color:#0F344C;"),
                        p("ORA tests whether a pre-defined set of genes (e.g., your significant DEGs)
                           overlaps with known gene sets more than expected by chance. It uses a
                           hypergeometric test (Fisher's exact test). You supply a ", strong("gene list"),
                           " (typically filtered by padj and log2FC cutoffs) and a ", strong("background"),
                           " (all measured genes)."),
                        hr(),
                        h4("Gene Set Enrichment Analysis (GSEA)", style="color:#0F344C;"),
                        p("GSEA takes a ", strong("ranked list of all genes"), " (not just significant ones)
                           and tests whether genes in a given set tend to cluster toward the top or
                           bottom of the ranking. It is more sensitive than ORA because it uses the
                           full distribution of changes, not a hard cutoff."),
                        hr(),
                        h4("Required Input Data", style="color:#0F344C;"),
                        div(class="control-group-panel",
                            p("Upload a DESeq2-style differential expression results table (CSV or TSV) with the following columns:"),
                            tags$ul(
                                tags$li(strong("Gene identifier column:"), " Gene symbols (e.g., TP53), Ensembl IDs (e.g., ENSG00000141510), or Entrez IDs (e.g., 7157)"),
                                tags$li(strong("log2FoldChange:"), " Log2 fold change between conditions"),
                                tags$li(strong("pvalue:"), " Raw p-value from the statistical test"),
                                tags$li(strong("padj:"), " Adjusted p-value (Benjamini-Hochberg)"),
                                tags$li(strong("baseMean:"), " (Optional) Mean normalized count across all samples")
                            ),
                            tags$pre("gene_symbol, baseMean,   log2FoldChange, pvalue,    padj\nTP53,        1245.32,   2.15,           1.3e-08,   5.6e-07\nBRCA1,       876.41,   -1.82,           3.1e-06,   9.2e-05")
                        ),
                        hr()
                    ),
                    column(2)
                )
            ),


        ## Data Input Page
##########################################################################################################################################################
            tabPanel('Data Input',
                fluidRow(
                    column(3,
                        wellPanel(
                            h2('Input Data', align='center'),
                            hr(),
                            materialSwitch("DemoData", label="Upload custom data", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.DemoData",
                                fileInput('de_file', 'Upload DESeq2 Results Table (CSV/TSV)',
                                    accept=c('text/csv', 'text/comma-separated-values, text/plain', '.csv',
                                             'text/tsv', 'text/tab-separated-values, text/plain', '.tsv'),
                                    multiple=FALSE
                                ),
                                uiOutput('column_selectors'),
                                hr(),
                                selectInput("organism", "Organism:",
                                    choices=c("Homo sapiens"="hsapiens", "Mus musculus"="mmusculus"),
                                    selected="hsapiens"),
                                selectInput("gene_id_type", "Gene ID Type:",
                                    choices=c("Gene Symbol"="SYMBOL", "Ensembl ID"="ENSEMBL", "Entrez ID"="ENTREZID"),
                                    selected="SYMBOL"),
                                selectInput("analysis_mode", "Analysis Mode:",
                                    choices=c("Over-Representation Analysis (ORA)"="ORA",
                                              "Gene Set Enrichment Analysis (GSEA)"="GSEA"),
                                    selected="ORA"),
                                hr(),
                                actionButton('submit', "Submit Data", class='btn btn-info btn-block')
                            ),
                            conditionalPanel("input.DemoData==false",
                                p("Use pre-loaded synthetic DESeq2 results to explore the tool's features."),
                                p(em("Demo dataset: ~5000 human genes with realistic fold-change and p-value distributions (Homo sapiens, SYMBOL IDs).")),
                                hr(),
                                selectInput("demo_mode", "Analysis Mode:",
                                    choices=c("Over-Representation Analysis (ORA)"="ORA",
                                              "Gene Set Enrichment Analysis (GSEA)"="GSEA"),
                                    selected="ORA"),
                                actionButton('demo_submit', "Load Demo Data", class='btn btn-success btn-block')
                            )
                        )
                    ),
                    column(9,
                        tabsetPanel(id='InputTables',
                            tabPanel(title='DESeq2 Results', hr(),
                                withSpinner(type=6, color='#5bc0de',
                                    dataTableOutput('input_table')
                                )
                            ),
                            tabPanel(title='Gene Summary', hr(),
                                withSpinner(type=6, color='#5bc0de',
                                    uiOutput('gene_summary_panel')
                                )
                            )
                        )
                    )
                )
            ),


        ## Enrichment Analysis Page (hidden until data submitted)
##########################################################################################################################################################
            tabPanel('Enrichment Analysis',
                fluidRow(
                    column(3,
                        wellPanel(

                            ## ── Database Selection ──
                            h5(strong("Database Selection"), style="color:#0F344C; margin-top:4px;"),
                            hr(),
                            selectInput("enrichment_db", "Pathway Database:",
                                choices=c(
                                    "GO Biological Process"="GO_BP",
                                    "GO Molecular Function"="GO_MF",
                                    "GO Cellular Component"="GO_CC",
                                    "KEGG Pathways"="KEGG",
                                    "Reactome Pathways"="Reactome",
                                    "MSigDB Hallmark"="Hallmark"
                                ),
                                selected="GO_BP"
                            ),

                            ## ── ORA-specific controls ──
                            conditionalPanel("output.current_mode == 'ORA'",
                                hr(),
                                h5(strong("ORA Filters"), style="color:#0F344C;"),
                                materialSwitch("ora_use_padj", label="Use padj (vs pvalue)",
                                    value=TRUE, right=TRUE, status='info'),
                                conditionalPanel("input.ora_use_padj",
                                    sliderInput("ora_padj_cutoff", "padj cutoff:",
                                        min=0.001, max=0.1, step=0.001, value=0.05)
                                ),
                                conditionalPanel("input.ora_use_padj == false",
                                    sliderInput("ora_pval_cutoff", "pvalue cutoff:",
                                        min=0.001, max=0.1, step=0.001, value=0.05)
                                ),
                                sliderInput("ora_lfc_cutoff", "abs(log2FC) cutoff:",
                                    min=0, max=5, step=0.25, value=1),
                                radioButtons("ora_direction", "Direction:", inline=TRUE,
                                    choices=c("Both"="both", "Up"="up", "Down"="down"),
                                    selected="both"),
                                uiOutput('ora_gene_count_display')
                            ),

                            ## ── GSEA-specific controls ──
                            conditionalPanel("output.current_mode == 'GSEA'",
                                hr(),
                                h5(strong("GSEA Parameters"), style="color:#0F344C;"),
                                selectInput("gsea_metric", "Ranking Metric:",
                                    choices=c(
                                        "log2FoldChange"="log2FC",
                                        "-log10(pvalue) * sign(log2FC)"="signed_pval",
                                        "stat (Wald statistic)"="stat"
                                    ),
                                    selected="log2FC"
                                ),
                                sliderInput("gsea_minGS", "Min gene set size:",
                                    min=5, max=100, step=5, value=15),
                                sliderInput("gsea_maxGS", "Max gene set size:",
                                    min=100, max=2000, step=50, value=500),
                                numericInput("gsea_nperm", "Permutations:",
                                    value=10000, min=1000, max=100000, step=1000)
                            ),

                            ## ── Run button ──
                            hr(),
                            actionButton("run_enrichment", "Run Enrichment Analysis",
                                class="btn btn-info btn-block", icon=icon("play")),
                            hr(),

                            ## ── Display Options ──
                            materialSwitch("show_display", label="Display Options", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.show_display",
                                hr(),
                                sliderInput("top_n", "Top pathways to display:",
                                    min=5, max=50, step=5, value=20),
                                materialSwitch("show_descriptions", label="Show pathway descriptions (vs IDs)",
                                    value=TRUE, right=TRUE, status='info'),
                                h6(strong("Text Sizes"), style="color:#0F344C; margin-top:10px;"),
                                sliderInput("text_size_axis_title", "Axis titles:",
                                    min=6, max=24, step=1, value=12),
                                sliderInput("text_size_x_labels", "X-axis labels:",
                                    min=6, max=20, step=1, value=10),
                                sliderInput("text_size_y_labels", "Y-axis / pathway labels:",
                                    min=6, max=20, step=1, value=10),
                                sliderInput("text_size_legend", "Legend labels:",
                                    min=6, max=20, step=1, value=10),
                                sliderInput("text_size_legend_title", "Legend title:",
                                    min=6, max=20, step=1, value=11),
                                hr(),
                                materialSwitch("wrap_names", label="Wrap long pathway names",
                                    value=TRUE, right=TRUE, status='info'),
                                conditionalPanel("input.wrap_names",
                                    sliderInput("wrap_width", "Wrap width (characters):",
                                        min=15, max=80, step=5, value=40)
                                ),
                                sliderInput("row_height", "Row height (px per pathway):",
                                    min=15, max=80, step=5, value=35)
                            ),

                            ## ── Color Options ──
                            materialSwitch("show_colors", label="Color Options", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.show_colors",
                                hr(),
                                selectInput("color_scale", "Color Scale:",
                                    choices=c("viridis"="viridis", "magma"="magma",
                                              "RdBu (diverging)"="RdBu", "Custom High/Low"="custom"),
                                    selected="viridis"
                                ),
                                conditionalPanel("input.color_scale == 'custom'",
                                    column(6, colourInput("color_high", "High", "#d73027")),
                                    column(6, colourInput("color_low", "Low", "#4575b4")),
                                    div(style="clear:both;")
                                ),
                                selectInput("color_by", "Color by:",
                                    choices=c("p.adjust"="p.adjust", "pvalue"="pvalue"),
                                    selected="p.adjust"
                                )
                            ),

                            ## ── Category Filter ──
                            materialSwitch("show_filter", label="Category Filter", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.show_filter",
                                hr(),
                                uiOutput("category_filter_ui")
                            ),

                            ## ── Resize ──
                            materialSwitch("show_resize", label="Resize Plot", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.show_resize",
                                hr(),
                                sliderInput("plot_height", "Plot height (px):",
                                    min=200, max=2000, step=50, value=700),
                                sliderInput("plot_width", "Plot width (px):",
                                    min=200, max=2000, step=50, value=900)
                            )

                        )# end wellPanel sidebar
                    ),
                    column(9,
                        tabsetPanel(id='EnrichmentTabs',
                            tabPanel(title='Dot Plot', hr(),
                                fluidRow(style="margin: 0 8px 4px 0;",
                                    column(12, align="right",
                                        actionButton("show_code_modal", label=NULL,
                                            icon=icon("file-code"),
                                            title="View R code to reproduce this analysis",
                                            class="btn btn-default btn-sm",
                                            style="border-radius:6px; font-size:16px; padding:4px 8px;"
                                        )
                                    )
                                ),
                                hr(),
                                plotOutput("dotplot_out", height='auto'),
                                uiOutput("dotplot_download_ui")
                            ),
                            tabPanel(title='Bar Plot', hr(),
                                plotOutput("barplot_out", height='auto'),
                                uiOutput("barplot_download_ui")
                            ),
                            tabPanel(title='GSEA Plot', hr(),
                                uiOutput("gsea_pathway_selector"),
                                plotOutput("gseaplot_out", height='auto'),
                                uiOutput("gseaplot_download_ui")
                            ),
                            tabPanel(title='Enrichment Map', hr(),
                                plotOutput("emapplot_out", height='auto'),
                                uiOutput("emapplot_download_ui")
                            ),
                            tabPanel(title='Gene-Concept Network', hr(),
                                plotOutput("cnetplot_out", height='auto'),
                                uiOutput("cnetplot_download_ui")
                            ),
                            tabPanel(title='Results Table', hr(),
                                withSpinner(type=6, color='#5bc0de',
                                    dataTableOutput('enrichment_results_table')
                                ),
                                div(style="margin-top:30px; text-align:center; padding-bottom:50px;",
                                    downloadButton('download_results_csv', 'Download Results (CSV)')
                                )
                            )
                        )
                    )
                )
            ),


        ## Footer
##########################################################################################################################################################
            tags$footer(
                wellPanel(
                    fluidRow(
                        column(4, align='center',
                        tags$a(href="https://github.com/MaGIC-Analytics/magic-enrichment", icon("github", "fa-3x")),
                        tags$h4('GitHub to submit issues/requests')
                        ),
                        column(4, align='center',
                        tags$a(href="http://www.bioinformagic.io/", icon("magic", "fa-3x")),
                        tags$h4('MaGIC Home Page')
                        ),
                        column(4, align='center',
                        tags$a(href="https://github.com/MaGIC-Analytics", icon("address-card", "fa-3x")),
                        tags$h4("Developer's Page")
                        )
                    ),
                    fluidRow(
                        column(12, align='center',
                            HTML('<a href="https://www.youtube.com/watch?v=dQw4w9WgXcQ">
                            <p>&copy;
                                <script language="javascript" type="text/javascript">
                                var today = new Date()
                                var year = today.getFullYear()
                                document.write(year)
                                </script>
                            </p>
                            </a>
                            ')
                        )
                    )
                )
            )
        )# Ends navbarPage
    )# Ends fluidPage
)# Ends tagList
