# ─── Tab Visibility Management ─────────────────────────────────────────────

# Hide Enrichment Analysis tab on initial load
observe({
    hideTab(inputId="NAVTABS", target="Enrichment Analysis")
})

# Show tab after custom data submit
observeEvent(input$submit, {
    dat <- ProcessedData()
    if (!is.null(dat) && nrow(dat) > 0) {
        showTab(inputId="NAVTABS", target="Enrichment Analysis")
        updateTabsetPanel(session, inputId="NAVTABS", selected="Enrichment Analysis")
        shinyjs::delay(300, shinyjs::runjs("$(window).trigger('resize');"))
    }
})

# Show tab after demo data submit
observeEvent(input$demo_submit, {
    showTab(inputId="NAVTABS", target="Enrichment Analysis")
    updateTabsetPanel(session, inputId="NAVTABS", selected="Enrichment Analysis")
    shinyjs::delay(300, shinyjs::runjs("$(window).trigger('resize');"))
})

# ─── GSEA-only tabs: GSEA Plot ───────────────────────────────────────────
observe({
    mode <- AnalysisMode()
    if (mode == "ORA") {
        hideTab(inputId="EnrichmentTabs", target="GSEA Plot")
    } else {
        showTab(inputId="EnrichmentTabs", target="GSEA Plot")
    }
})

# ─── NES color option: only valid in GSEA mode ──────────────────────────
observe({
    mode <- AnalysisMode()
    if (mode == "ORA") {
        updateSelectInput(session, "color_by",
            choices=c("p.adjust"="p.adjust", "pvalue"="pvalue"))
    } else {
        updateSelectInput(session, "color_by",
            choices=c("p.adjust"="p.adjust", "pvalue"="pvalue", "NES"="NES"))
    }
})

# ─── Notify when organism/gene ID type changes ──────────────────────────
observeEvent(c(input$organism, input$gene_id_type), {
    showNotification("Organism or gene ID type changed. Re-run enrichment to update results.",
        type='warning', duration=6)
}, ignoreInit=TRUE)
