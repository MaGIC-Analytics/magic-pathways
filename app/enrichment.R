# ═══════════════════════════════════════════════════════════════════════════
# enrichment.R — Core analysis engine, plot renderers, downloads, modals
# ═══════════════════════════════════════════════════════════════════════════

# ─── State Tracking ─────────────────────────────────────────────────────
LastRunMode <- reactiveVal(NULL)

PlotStatus <- reactive({
    res <- tryCatch(EnrichmentResult(), error = function(e) NULL)
    if (is.null(res)) return("none")
    if (!is.null(LastRunMode()) && LastRunMode() != AnalysisMode()) return("stale")
    "ready"
})

# ─── Helpers ────────────────────────────────────────────────────────────
get_text_sizes <- function(input) {
    list(
        axis_title   = input$text_size_axis_title   %||% 12,
        x_labels     = input$text_size_x_labels     %||% 10,
        y_labels     = input$text_size_y_labels     %||% 10,
        legend       = input$text_size_legend       %||% 10,
        legend_title = input$text_size_legend_title %||% 11
    )
}

get_color_scale <- function(scale_name, high, low, aesthetic = "color") {
    if (aesthetic == "fill") {
        switch(scale_name,
            viridis = scale_fill_viridis_c(option = "viridis", direction = -1),
            magma   = scale_fill_viridis_c(option = "magma",   direction = -1),
            RdBu    = scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu"))),
            custom  = scale_fill_gradient(low = low, high = high),
            NULL)
    } else {
        switch(scale_name,
            viridis = scale_color_viridis_c(option = "viridis", direction = -1),
            magma   = scale_color_viridis_c(option = "magma",   direction = -1),
            RdBu    = scale_color_gradientn(colors = rev(brewer.pal(11, "RdBu"))),
            custom  = scale_color_gradient(low = low, high = high),
            NULL)
    }
}

apply_plot_theme <- function(p, sizes, wrap, wrap_width, color_obj = NULL) {
    if (isTRUE(wrap)) {
        p <- p + scale_y_discrete(labels = function(x) str_wrap(x, width = wrap_width))
    }
    if (!is.null(color_obj)) {
        p <- tryCatch(p + color_obj, error = function(e) p)
    }
    p + theme(
        axis.title     = element_text(size = sizes$axis_title),
        axis.text.x    = element_text(size = sizes$x_labels),
        axis.text.y    = element_text(size = sizes$y_labels),
        legend.text    = element_text(size = sizes$legend),
        legend.title   = element_text(size = sizes$legend_title)
    )
}

draw_placeholder <- function(msg, col = "#999999") {
    par(mar = c(0, 0, 0, 0))
    plot.new()
    text(0.5, 0.5, msg, cex = 1.3, col = col, font = 3)
}

save_plot_to_file <- function(p, file, fmt, h, w) {
    if (fmt == "png")       png(file, height = h, width = w, res = 150)
    else if (fmt == "jpeg") jpeg(file, height = h, width = w, res = 150)
    else if (fmt == "tiff") tiff(file, height = h, width = w, res = 150)
    else if (fmt == "pdf")  pdf(file, height = h / 96, width = w / 96)
    else if (fmt == "svg")  svg(file, height = h / 96, width = w / 96)
    else if (fmt == "eps")  { setEPS(); postscript(file, height = h / 96, width = w / 96) }
    tryCatch(print(p), error = function(e) NULL)
    dev.off()
}


# ═══════════════════════════════════════════════════════════════════════════
# Core Enrichment Analysis
# ═══════════════════════════════════════════════════════════════════════════

EnrichmentResult <- eventReactive(input$run_enrichment, {
    df     <- ProcessedData()
    req(df)
    mode   <- AnalysisMode()
    db     <- input$enrichment_db
    orgdb  <- OrgDb()
    kegg   <- OrgKegg()
    react  <- OrgReactome()
    species_name <- OrgSpeciesName()

    withProgress(message = 'Running enrichment analysis...', value = 0.3, {

        if (mode == "ORA") {
            use_padj  <- isTRUE(input$ora_use_padj)
            sig_cut   <- if (use_padj) (input$ora_padj_cutoff %||% 0.05) else (input$ora_pval_cutoff %||% 0.05)
            lfc_cut   <- input$ora_lfc_cutoff %||% 1
            direction <- input$ora_direction %||% "both"

            if (use_padj) {
                sig <- df[!is.na(df$padj) & df$padj < sig_cut & abs(df$log2FC) >= lfc_cut, ]
            } else {
                sig <- df[!is.na(df$pvalue) & df$pvalue < sig_cut & abs(df$log2FC) >= lfc_cut, ]
            }
            if (direction == "up")   sig <- sig[sig$log2FC > 0, ]
            if (direction == "down") sig <- sig[sig$log2FC < 0, ]

            gene_list <- unique(sig$ENTREZID)
            universe  <- unique(df$ENTREZID)

            shiny::validate(need(length(gene_list) >= 5,
                "Fewer than 5 genes pass the current filters. Relax your cutoffs."))

            incProgress(0.3, detail = "Running ORA...")

            result <- tryCatch({
                if (db == "GO_BP") {
                    enrichGO(gene = gene_list, universe = universe, OrgDb = orgdb,
                             keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                             pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
                } else if (db == "GO_MF") {
                    enrichGO(gene = gene_list, universe = universe, OrgDb = orgdb,
                             keyType = "ENTREZID", ont = "MF", pAdjustMethod = "BH",
                             pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
                } else if (db == "GO_CC") {
                    enrichGO(gene = gene_list, universe = universe, OrgDb = orgdb,
                             keyType = "ENTREZID", ont = "CC", pAdjustMethod = "BH",
                             pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
                } else if (db == "KEGG") {
                    res <- enrichKEGG(gene = gene_list, universe = universe,
                                      organism = kegg, keyType = "ncbi-geneid",
                                      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
                    tryCatch(setReadable(res, OrgDb = orgdb, keyType = "ENTREZID"),
                             error = function(e) res)
                } else if (db == "Reactome") {
                    enrichPathway(gene = gene_list, universe = universe,
                                  organism = react, pAdjustMethod = "BH",
                                  pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
                } else if (db == "Hallmark") {
                    h_df <- msigdbr(species = species_name, category = "H")
                    term2gene <- h_df[, c("gs_name", "entrez_gene")]
                    enricher(gene = gene_list, universe = universe, TERM2GENE = term2gene,
                             pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
                }
            }, error = function(e) {
                showNotification(paste("Enrichment error:", e$message), type = 'error', duration = NULL)
                NULL
            })

        } else {
            # ── GSEA: build ranked gene vector ──
            metric <- input$gsea_metric %||% "log2FC"

            if (metric == "log2FC") {
                df$rank_value <- df$log2FC
            } else if (metric == "signed_pval") {
                df$rank_value <- -log10(pmax(df$pvalue, 1e-300)) * sign(df$log2FC)
            } else if (metric == "stat") {
                shiny::validate(need("stat" %in% colnames(df),
                    "The 'stat' column is not present in your data. Choose a different metric."))
                df$rank_value <- df$stat
            }

            df <- df[order(df$rank_value, decreasing = TRUE), ]
            df <- df[!duplicated(df$ENTREZID), ]
            gene_vec <- setNames(df$rank_value, df$ENTREZID)

            minGS <- input$gsea_minGS %||% 15
            maxGS <- input$gsea_maxGS %||% 500
            nperm <- input$gsea_nperm %||% 10000

            incProgress(0.3, detail = "Running GSEA...")

            result <- tryCatch({
                if (db == "GO_BP") {
                    gseGO(geneList = gene_vec, OrgDb = orgdb, keyType = "ENTREZID",
                          ont = "BP", minGSSize = minGS, maxGSSize = maxGS,
                          pvalueCutoff = 0.05, pAdjustMethod = "BH",
                          nPermSimple = nperm, verbose = FALSE)
                } else if (db == "GO_MF") {
                    gseGO(geneList = gene_vec, OrgDb = orgdb, keyType = "ENTREZID",
                          ont = "MF", minGSSize = minGS, maxGSSize = maxGS,
                          pvalueCutoff = 0.05, pAdjustMethod = "BH",
                          nPermSimple = nperm, verbose = FALSE)
                } else if (db == "GO_CC") {
                    gseGO(geneList = gene_vec, OrgDb = orgdb, keyType = "ENTREZID",
                          ont = "CC", minGSSize = minGS, maxGSSize = maxGS,
                          pvalueCutoff = 0.05, pAdjustMethod = "BH",
                          nPermSimple = nperm, verbose = FALSE)
                } else if (db == "KEGG") {
                    res <- gseKEGG(geneList = gene_vec, organism = kegg,
                                   keyType = "ncbi-geneid",
                                   minGSSize = minGS, maxGSSize = maxGS,
                                   pvalueCutoff = 0.05, pAdjustMethod = "BH",
                                   nPermSimple = nperm, verbose = FALSE)
                    tryCatch(setReadable(res, OrgDb = orgdb, keyType = "ENTREZID"),
                             error = function(e) res)
                } else if (db == "Reactome") {
                    gsePathway(geneList = gene_vec, organism = react,
                               minGSSize = minGS, maxGSSize = maxGS,
                               pvalueCutoff = 0.05, pAdjustMethod = "BH",
                               nPermSimple = nperm, verbose = FALSE)
                } else if (db == "Hallmark") {
                    h_df <- msigdbr(species = species_name, category = "H")
                    term2gene <- h_df[, c("gs_name", "entrez_gene")]
                    GSEA(geneList = gene_vec, TERM2GENE = term2gene,
                         minGSSize = minGS, maxGSSize = maxGS,
                         pvalueCutoff = 0.05, pAdjustMethod = "BH",
                         nPermSimple = nperm, verbose = FALSE)
                }
            }, error = function(e) {
                showNotification(paste("GSEA error:", e$message), type = 'error', duration = NULL)
                NULL
            })
        }

        incProgress(0.4, detail = "Done!")
    })

    shiny::validate(need(!is.null(result) && nrow(as.data.frame(result)) > 0,
        "No significant enrichment results found. Try different parameters or a different database."))

    LastRunMode(mode)
    result
})


# ─── Category Filter UI ────────────────────────────────────────────────
output$category_filter_ui <- renderUI({
    res <- tryCatch(EnrichmentResult(), error = function(e) NULL)
    if (is.null(res)) {
        return(p("Run an enrichment analysis first.", style = "color:#999;"))
    }
    df_res <- as.data.frame(res)
    categories <- unique(df_res$Description)
    selectizeInput("category_filter", "Filter pathways:", choices = categories,
        selected = NULL, multiple = TRUE,
        options = list(placeholder = "Type to filter...", maxOptions = 200))
})

# ─── Filtered Result ───────────────────────────────────────────────────
FilteredResult <- reactive({
    res <- EnrichmentResult()
    req(res)
    df_res <- as.data.frame(res)

    if (!is.null(input$category_filter) && length(input$category_filter) > 0) {
        keep_idx <- which(df_res$Description %in% input$category_filter)
        if (length(keep_idx) > 0) res <- res[keep_idx, ]
    }

    if (!isTRUE(input$show_descriptions) && input$enrichment_db == "Hallmark") {
        df_mod <- as.data.frame(res)
        df_mod$Description <- gsub("HALLMARK_", "", df_mod$Description)
        df_mod$Description <- gsub("_", " ", df_mod$Description)
        res@result <- df_mod
    }

    res
})


# ═══════════════════════════════════════════════════════════════════════════
# Plot Reactives
# ═══════════════════════════════════════════════════════════════════════════

dotplot_plotter <- reactive({
    res <- FilteredResult()
    req(res)
    n         <- input$top_n %||% 20
    color_var <- input$color_by %||% "p.adjust"
    sizes     <- get_text_sizes(input)
    wrap_w    <- input$wrap_width %||% 40
    cs <- get_color_scale(input$color_scale %||% "viridis",
                          input$color_high %||% "#d73027",
                          input$color_low  %||% "#4575b4",
                          aesthetic = "color")

    p <- dotplot(res, showCategory = n, color = color_var)
    suppressWarnings(apply_plot_theme(p, sizes, input$wrap_names, wrap_w, cs))
})

barplot_plotter <- reactive({
    res <- FilteredResult()
    req(res)
    n         <- input$top_n %||% 20
    color_var <- input$color_by %||% "p.adjust"
    sizes     <- get_text_sizes(input)
    wrap_w    <- input$wrap_width %||% 40
    cs <- get_color_scale(input$color_scale %||% "viridis",
                          input$color_high %||% "#d73027",
                          input$color_low  %||% "#4575b4",
                          aesthetic = "fill")

    p <- tryCatch(
        barplot(res, showCategory = n, color = color_var),
        error = function(e) {
            df_res <- as.data.frame(res)
            df_res <- head(df_res[order(df_res$p.adjust), ], n)
            df_res$Description <- factor(df_res$Description, levels = rev(df_res$Description))
            x_col <- if ("NES" %in% colnames(df_res)) "NES"
                     else if ("Count" %in% colnames(df_res)) "Count"
                     else "setSize"
            ggplot(df_res, aes(x = .data[[x_col]], y = Description, fill = .data[[color_var]])) +
                geom_col() + theme_minimal() +
                labs(x = x_col, y = NULL, fill = color_var)
        }
    )
    suppressWarnings(apply_plot_theme(p, sizes, input$wrap_names, wrap_w, cs))
})


emapplot_plotter <- reactive({
    res <- FilteredResult()
    req(res)
    n     <- input$top_n %||% 20
    sizes <- get_text_sizes(input)

    res_pw <- tryCatch(pairwise_termsim(res), error = function(e) {
        showNotification("Could not compute term similarity.", type = 'warning')
        NULL
    })
    req(res_pw)

    p <- emapplot(res_pw, showCategory = n, color = input$color_by %||% "p.adjust",
                  layout = "nicely")
    p + theme(
        text           = element_text(size = sizes$legend),
        legend.text    = element_text(size = sizes$legend),
        legend.title   = element_text(size = sizes$legend_title)
    )
})

cnetplot_plotter <- reactive({
    res <- FilteredResult()
    req(res)
    n     <- min(input$top_n %||% 20, 10)
    sizes <- get_text_sizes(input)

    if (AnalysisMode() == "GSEA") {
        df <- ProcessedData()
        fc <- setNames(df$log2FC, df$ENTREZID)
        p <- tryCatch(
            cnetplot(res, showCategory = n, foldChange = fc, node_label = "all"),
            error = function(e) cnetplot(res, showCategory = n, node_label = "all")
        )
    } else {
        p <- cnetplot(res, showCategory = n, node_label = "all")
    }
    p + theme(
        text           = element_text(size = sizes$legend),
        legend.text    = element_text(size = sizes$legend),
        legend.title   = element_text(size = sizes$legend_title)
    )
})

# ─── GSEA Running Score Plot ────────────────────────────────────────────
output$gsea_pathway_selector <- renderUI({
    if (PlotStatus() != "ready" || AnalysisMode() != "GSEA") {
        return(div(style = "text-align:center; color:#999; padding:10px;",
            p(em("Run a GSEA analysis first, then select a pathway to view its running score plot."))
        ))
    }
    res <- tryCatch(FilteredResult(), error = function(e) NULL)
    if (is.null(res)) return(NULL)
    df_res <- as.data.frame(res)
    choices <- setNames(df_res$ID, df_res$Description)
    tagList(
        selectizeInput("gsea_pathway_id", "Select a pathway:",
            choices = choices, selected = choices[1], width = "100%"),
        hr()
    )
})

gseaplot_plotter <- reactive({
    req(AnalysisMode() == "GSEA")
    req(PlotStatus() == "ready")
    res <- FilteredResult()
    req(res)
    pathway_id <- input$gsea_pathway_id
    req(pathway_id)

    df_res <- as.data.frame(res)
    req(pathway_id %in% df_res$ID)

    idx   <- which(df_res$ID == pathway_id)
    title <- df_res$Description[idx]

    gseaplot2(res, geneSetID = pathway_id, title = title,
              color = "green4", pvalue_table = TRUE)
})


# ═══════════════════════════════════════════════════════════════════════════
# Plot Renderers (with placeholder logic)
# ═══════════════════════════════════════════════════════════════════════════

output$dotplot_out <- renderPlot({
    status <- PlotStatus()
    if (status == "none") return(draw_placeholder("Click 'Run Enrichment Analysis'\nto generate plots"))
    if (status == "stale") return(draw_placeholder("Analysis mode changed.\nRerun Enrichment Analysis to update.", "#b8860b"))
    tryCatch(dotplot_plotter(), error = function(e)
        draw_placeholder(paste("Dot plot error:\n", e$message), "#cc0000"))
}, height = function() max(400, (input$top_n %||% 20) * (input$row_height %||% 35)),
   width  = function() input$plot_width %||% 900)

output$barplot_out <- renderPlot({
    status <- PlotStatus()
    if (status == "none") return(draw_placeholder("Click 'Run Enrichment Analysis'\nto generate plots"))
    if (status == "stale") return(draw_placeholder("Analysis mode changed.\nRerun Enrichment Analysis to update.", "#b8860b"))
    tryCatch(barplot_plotter(), error = function(e)
        draw_placeholder(paste("Bar plot error:\n", e$message), "#cc0000"))
}, height = function() max(400, (input$top_n %||% 20) * (input$row_height %||% 35)),
   width  = function() input$plot_width %||% 900)

output$emapplot_out <- renderPlot({
    status <- PlotStatus()
    if (status == "none") return(draw_placeholder("Click 'Run Enrichment Analysis'\nto generate plots"))
    if (status == "stale") return(draw_placeholder("Analysis mode changed.\nRerun Enrichment Analysis to update.", "#b8860b"))
    tryCatch(emapplot_plotter(), error = function(e)
        draw_placeholder(paste("Enrichment map error:\n", e$message), "#cc0000"))
}, height = function() input$plot_height %||% 700,
   width  = function() input$plot_width  %||% 900)

output$cnetplot_out <- renderPlot({
    status <- PlotStatus()
    if (status == "none") return(draw_placeholder("Click 'Run Enrichment Analysis'\nto generate plots"))
    if (status == "stale") return(draw_placeholder("Analysis mode changed.\nRerun Enrichment Analysis to update.", "#b8860b"))
    tryCatch(cnetplot_plotter(), error = function(e)
        draw_placeholder(paste("Gene-concept network error:\n", e$message), "#cc0000"))
}, height = function() input$plot_height %||% 700,
   width  = function() input$plot_width  %||% 900)

output$gseaplot_out <- renderPlot({
    status <- PlotStatus()
    if (status == "none") return(draw_placeholder("Click 'Run Enrichment Analysis'\nto generate plots"))
    if (status == "stale") return(draw_placeholder("Analysis mode changed.\nRerun Enrichment Analysis to update.", "#b8860b"))
    if (AnalysisMode() != "GSEA") return(draw_placeholder("GSEA plots are only available\nin GSEA mode."))
    if (is.null(input$gsea_pathway_id) || input$gsea_pathway_id == "")
        return(draw_placeholder("Select a pathway from the dropdown\nabove to view its running score plot."))
    tryCatch(gseaplot_plotter(), error = function(e)
        draw_placeholder(paste("GSEA plot error:\n", e$message), "#cc0000"))
}, height = function() input$plot_height %||% 700,
   width  = function() input$plot_width  %||% 900)


# ═══════════════════════════════════════════════════════════════════════════
# Download UIs and Handlers (per plot)
# ═══════════════════════════════════════════════════════════════════════════

plot_names <- c("dotplot", "barplot", "emapplot", "cnetplot", "gseaplot")
row_based  <- c("dotplot", "barplot")

for (pn in plot_names) {
    local({
        my_name  <- pn
        is_row   <- my_name %in% row_based
        my_plotter <- switch(my_name,
            dotplot   = dotplot_plotter,
            barplot   = barplot_plotter,
            emapplot  = emapplot_plotter,
            cnetplot  = cnetplot_plotter,
            gseaplot  = gseaplot_plotter)

        output[[paste0(my_name, "_download_ui")]] <- renderUI({
            if (PlotStatus() != "ready") return(NULL)
            if (my_name == "gseaplot" && AnalysisMode() != "GSEA") return(NULL)
            div(style = "margin-top:15px; text-align:center; padding-bottom:20px;",
                div(style = "width:180px; margin:0 auto;",
                    selectInput(paste0(my_name, "_dl_format"), "Format:",
                        choices = c('png','pdf','svg','tiff','jpeg','eps'), width = "100%"),
                    downloadButton(paste0('download_', my_name), 'Download Plot',
                        style = "width:100%;")
                )
            )
        })

        output[[paste0("download_", my_name)]] <- downloadHandler(
            filename = function() {
                fmt <- input[[paste0(my_name, "_dl_format")]] %||% "png"
                paste0("enrichment_", my_name, ".", fmt)
            },
            content = function(file) {
                p <- tryCatch(my_plotter(), error = function(e) NULL)
                req(p)
                fmt <- input[[paste0(my_name, "_dl_format")]] %||% "png"
                h <- if (is_row) max(400, (input$top_n %||% 20) * (input$row_height %||% 35))
                     else input$plot_height %||% 700
                w <- input$plot_width %||% 900
                save_plot_to_file(p, file, fmt, h, w)
            }
        )
    })
}


# ═══════════════════════════════════════════════════════════════════════════
# Results Table
# ═══════════════════════════════════════════════════════════════════════════

output$enrichment_results_table <- DT::renderDataTable({
    res <- FilteredResult()
    req(res)
    df <- as.data.frame(res)
    DT::datatable(df, style = 'bootstrap',
        options = list(pageLength = 20, scrollX = TRUE),
        filter = 'top')
})

output$download_results_csv <- downloadHandler(
    filename = function() paste0("enrichment_results_", input$enrichment_db, ".csv"),
    content = function(file) {
        res <- FilteredResult()
        req(res)
        write.csv(as.data.frame(res), file, row.names = FALSE)
    }
)


# ═══════════════════════════════════════════════════════════════════════════
# Help Modal
# ═══════════════════════════════════════════════════════════════════════════

show_enrichment_help_ui <- function() {
    showModal(modalDialog(
        title     = tagList(icon("circle-question"), " Gene Set Enrichment Tool Help"),
        size      = "l",
        easyClose = TRUE,
        footer    = modalButton("Close"),
        tabsetPanel(
            tabPanel("Overview",
                br(),
                h4("Over-Representation Analysis (ORA)"),
                p("Tests whether a pre-defined list of significant genes overlaps with known
                   gene sets more than expected by chance (hypergeometric test). You define
                   significance cutoffs (padj, log2FC) to select the gene list."),
                h4("Gene Set Enrichment Analysis (GSEA)"),
                p("Takes a ranked list of ALL genes and tests whether members of a gene set
                   tend to appear toward the top or bottom of the ranking. More sensitive than
                   ORA because it uses the entire distribution, not a hard cutoff."),
                h4("When to Use Each"),
                tags$ul(
                    tags$li(strong("ORA:"), " When you have clear significant DEGs and want to know which pathways they belong to."),
                    tags$li(strong("GSEA:"), " When changes are subtle/distributed, or when you want to avoid arbitrary cutoffs.")
                )
            ),
            tabPanel("Input Data",
                br(),
                h4("DESeq2 Results Format"),
                p("Upload a CSV or TSV with columns for gene identifier, log2FoldChange, pvalue, and padj.
                   The tool will map your gene IDs to Entrez IDs for the selected organism."),
                h4("Gene ID Types"),
                tags$ul(
                    tags$li(strong("SYMBOL:"), " Standard gene symbols (e.g., TP53, BRCA1)"),
                    tags$li(strong("ENSEMBL:"), " Ensembl gene IDs (e.g., ENSG00000141510)"),
                    tags$li(strong("ENTREZID:"), " NCBI Entrez gene IDs (e.g., 7157)")
                ),
                h4("Supported Organisms"),
                tags$ul(
                    tags$li("Homo sapiens (human)"),
                    tags$li("Mus musculus (mouse)")
                )
            ),
            tabPanel("Databases",
                br(),
                h4("GO (Gene Ontology)"),
                p("Three sub-ontologies: Biological Process (BP), Molecular Function (MF), Cellular Component (CC)."),
                h4("KEGG"),
                p("Kyoto Encyclopedia of Genes and Genomes pathway database."),
                h4("Reactome"),
                p("Peer-reviewed, curated pathway database."),
                h4("MSigDB Hallmark"),
                p("50 well-defined gene sets representing specific biological states or processes
                   from the Molecular Signatures Database.")
            ),
            tabPanel("Plot Types",
                br(),
                h4("Dot Plot"),
                p("Shows enriched terms with dot size = gene count and color = significance."),
                h4("Bar Plot"),
                p("Bar chart of enriched terms colored by significance."),
                h4("GSEA Plot (GSEA only)"),
                p("Classic running enrichment score plot for a single pathway. Shows the cumulative
                   enrichment score walking down the ranked gene list, with the gene hit positions
                   and the ranked list metric."),
                h4("Enrichment Map"),
                p("Network of enriched terms connected by shared genes (requires pairwise term similarity)."),
                h4("Gene-Concept Network"),
                p("Links enriched terms to their member genes.")
            )
        )
    ))
}

observeEvent(input$show_help_float, { show_enrichment_help_ui() })


# ═══════════════════════════════════════════════════════════════════════════
# Reproducible Code Modal
# ═══════════════════════════════════════════════════════════════════════════

observeEvent(input$show_code_modal, {
    res <- isolate(tryCatch(EnrichmentResult(), error = function(e) NULL))
    if (is.null(res)) {
        code <- "# Run an enrichment analysis first, then click here to get reproducible code."
    } else {
        mode <- isolate(AnalysisMode())
        db   <- isolate(input$enrichment_db)
        orgdb_name <- isolate(if (OrgKegg() == "hsa") "org.Hs.eg.db" else "org.Mm.eg.db")
        kegg_org   <- isolate(OrgKegg())
        react_org  <- isolate(OrgReactome())
        sp_name    <- isolate(OrgSpeciesName())

        lib_lines <- "library(clusterProfiler)\nlibrary(enrichplot)\nlibrary(DOSE)\n"
        if (db == "Reactome") lib_lines <- paste0(lib_lines, "library(ReactomePA)\n")
        if (db == "Hallmark") lib_lines <- paste0(lib_lines, "library(msigdbr)\n")
        lib_lines <- paste0(lib_lines, sprintf("library(%s)\n", orgdb_name))

        if (mode == "ORA") {
            use_padj  <- isolate(isTRUE(input$ora_use_padj))
            sig_cut   <- isolate(if (use_padj) (input$ora_padj_cutoff %||% 0.05) else (input$ora_pval_cutoff %||% 0.05))
            lfc_cut   <- isolate(input$ora_lfc_cutoff %||% 1)
            direction <- isolate(input$ora_direction %||% "both")
            sig_col   <- if (use_padj) "padj" else "pvalue"

            filter_line <- sprintf("sig <- df[df$%s < %s & abs(df$log2FoldChange) >= %s, ]",
                sig_col, sig_cut, lfc_cut)
            if (direction == "up") filter_line <- paste0(filter_line, "\nsig <- sig[sig$log2FoldChange > 0, ]")
            if (direction == "down") filter_line <- paste0(filter_line, "\nsig <- sig[sig$log2FoldChange < 0, ]")

            analysis_call <- if (db == "GO_BP" || db == "GO_MF" || db == "GO_CC") {
                ont <- gsub("GO_", "", db)
                sprintf('result <- enrichGO(\n    gene     = gene_list,\n    universe = universe,\n    OrgDb    = %s,\n    keyType  = "ENTREZID",\n    ont      = "%s",\n    pAdjustMethod = "BH",\n    pvalueCutoff  = 0.05,\n    readable = TRUE\n)', orgdb_name, ont)
            } else if (db == "KEGG") {
                sprintf('result <- enrichKEGG(\n    gene     = gene_list,\n    universe = universe,\n    organism = "%s",\n    keyType  = "ncbi-geneid",\n    pAdjustMethod = "BH",\n    pvalueCutoff  = 0.05\n)\nresult <- setReadable(result, OrgDb = %s, keyType = "ENTREZID")', kegg_org, orgdb_name)
            } else if (db == "Reactome") {
                sprintf('result <- enrichPathway(\n    gene     = gene_list,\n    universe = universe,\n    organism = "%s",\n    pAdjustMethod = "BH",\n    pvalueCutoff  = 0.05,\n    readable = TRUE\n)', react_org)
            } else {
                sprintf('h_df <- msigdbr(species = "%s", category = "H")\nterm2gene <- h_df[, c("gs_name", "entrez_gene")]\nresult <- enricher(\n    gene     = gene_list,\n    universe = universe,\n    TERM2GENE = term2gene,\n    pAdjustMethod = "BH",\n    pvalueCutoff  = 0.05\n)', sp_name)
            }

            code <- paste0(
                lib_lines,
                "\n# ── Load your data ──\n",
                "# df <- read.csv('your_deseq2_results.csv')\n\n",
                "# ── Filter significant genes ──\n",
                filter_line, "\n",
                "gene_list <- unique(sig$ENTREZID)\n",
                "universe  <- unique(df$ENTREZID)\n\n",
                "# ── Run ORA ──\n",
                analysis_call, "\n\n",
                "# ── Visualize ──\n",
                sprintf("dotplot(result, showCategory = %d)\n", isolate(input$top_n %||% 20))
            )
        } else {
            metric <- isolate(input$gsea_metric %||% "log2FC")
            minGS  <- isolate(input$gsea_minGS %||% 15)
            maxGS  <- isolate(input$gsea_maxGS %||% 500)
            nperm  <- isolate(input$gsea_nperm %||% 10000)

            rank_line <- if (metric == "log2FC") {
                "df$rank_value <- df$log2FoldChange"
            } else if (metric == "signed_pval") {
                "df$rank_value <- -log10(pmax(df$pvalue, 1e-300)) * sign(df$log2FoldChange)"
            } else {
                "df$rank_value <- df$stat"
            }

            analysis_call <- if (db == "GO_BP" || db == "GO_MF" || db == "GO_CC") {
                ont <- gsub("GO_", "", db)
                sprintf('result <- gseGO(\n    geneList = gene_vec,\n    OrgDb    = %s,\n    keyType  = "ENTREZID",\n    ont      = "%s",\n    minGSSize = %d,\n    maxGSSize = %d,\n    pvalueCutoff = 0.05,\n    nPermSimple  = %d,\n    verbose  = FALSE\n)', orgdb_name, ont, minGS, maxGS, nperm)
            } else if (db == "KEGG") {
                sprintf('result <- gseKEGG(\n    geneList = gene_vec,\n    organism = "%s",\n    keyType  = "ncbi-geneid",\n    minGSSize = %d,\n    maxGSSize = %d,\n    pvalueCutoff = 0.05,\n    nPermSimple  = %d,\n    verbose  = FALSE\n)\nresult <- setReadable(result, OrgDb = %s, keyType = "ENTREZID")', kegg_org, minGS, maxGS, nperm, orgdb_name)
            } else if (db == "Reactome") {
                sprintf('result <- gsePathway(\n    geneList = gene_vec,\n    organism = "%s",\n    minGSSize = %d,\n    maxGSSize = %d,\n    pvalueCutoff = 0.05,\n    nPermSimple  = %d,\n    verbose  = FALSE\n)', react_org, minGS, maxGS, nperm)
            } else {
                sprintf('h_df <- msigdbr(species = "%s", category = "H")\nterm2gene <- h_df[, c("gs_name", "entrez_gene")]\nresult <- GSEA(\n    geneList  = gene_vec,\n    TERM2GENE = term2gene,\n    minGSSize = %d,\n    maxGSSize = %d,\n    pvalueCutoff = 0.05,\n    nPermSimple  = %d,\n    verbose   = FALSE\n)', sp_name, minGS, maxGS, nperm)
            }

            code <- paste0(
                lib_lines,
                "\n# ── Load your data ──\n",
                "# df <- read.csv('your_deseq2_results.csv')\n\n",
                "# ── Build ranked gene vector ──\n",
                rank_line, "\n",
                "df <- df[order(df$rank_value, decreasing = TRUE), ]\n",
                "df <- df[!duplicated(df$ENTREZID), ]\n",
                "gene_vec <- setNames(df$rank_value, df$ENTREZID)\n\n",
                "# ── Run GSEA ──\n",
                analysis_call, "\n\n",
                "# ── Visualize ──\n",
                sprintf("dotplot(result, showCategory = %d)\n", isolate(input$top_n %||% 20)),
                "gseaplot2(result, geneSetID = 1, title = as.data.frame(result)$Description[1])\n"
            )
        }
    }

    showModal(modalDialog(
        title     = tagList(icon("file-code"), " Reproducible R Code"),
        size      = "l",
        easyClose = TRUE,
        footer    = modalButton("Close"),
        p("Copy this code to reproduce your current analysis in an offline R session.",
          style = "color:#555; margin-bottom:12px;"),
        tags$pre(
            style = paste(
                "background:#1e1e1e; color:#d4d4d4; border-radius:6px;",
                "padding:16px; font-size:12px; max-height:520px; overflow-y:auto;",
                "white-space:pre; font-family:'Courier New', monospace;"
            ),
            code
        )
    ))
})
