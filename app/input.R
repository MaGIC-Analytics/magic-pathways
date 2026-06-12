# ─── Null-coalescing operator ─────────────���────────────────────────────────
`%||%` <- function(a, b) {
    if (is.null(a)) return(b)
    if (length(a) == 0) return(b)
    if (length(a) == 1 && is.character(a) && !nzchar(a)) return(b)
    a
}

# ─── Auto-detect delimiter ────────────────────────────────────────────────
read_delim_auto <- function(path) {
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c("tsv", "txt")) {
        fread(path, sep="\t")
    } else {
        fread(path, sep=",")
    }
}

# ── msigdbr version-compatibility shim ─────────────────────────────────────────
# msigdbr >= 10 renamed category/subcategory -> collection/subcollection, split
# KEGG into CP:KEGG_LEGACY/CP:KEGG_MEDICUS, and renamed entrez_gene -> ncbi_gene.
# Take the legacy argument names, try the new API first, fall back to the old one,
# and guarantee the gs_name / gene_symbol / entrez_gene columns callers expect.
msigdbr_compat <- function(species, category, subcategory = NULL) {
    sub_new <- if (identical(subcategory, "CP:KEGG")) "CP:KEGG_LEGACY" else subcategory
    df <- tryCatch(
        do.call(msigdbr, c(list(species = species, collection = category),
                           if (!is.null(sub_new)) list(subcollection = sub_new))),
        error = function(e)
            do.call(msigdbr, c(list(species = species, category = category),
                               if (!is.null(subcategory)) list(subcategory = subcategory)))
    )
    if (!"entrez_gene" %in% colnames(df) && "ncbi_gene" %in% colnames(df)) {
        df$entrez_gene <- df$ncbi_gene
    }
    df
}

# ─── Cached demo data (generated once per session) ───────────────────────
DemoDataCache <- reactiveVal(NULL)

generate_demo_data <- function() {
    set.seed(42)

    # ── Hardcoded pathway gene sets (known biology, will enrich) ──────────
    # These are real human gene symbols belonging to well-characterized pathways.
    # Upregulated pathway genes (cell cycle, DNA replication, p53 signaling)
    up_genes <- c(
        # Cell cycle (GO:0007049 / KEGG hsa04110)
        "CDK1", "CDK2", "CDK4", "CDK6", "CCNA2", "CCNB1", "CCNB2", "CCND1",
        "CCNE1", "CCNE2", "CDC6", "CDC20", "CDC25A", "CDC25C", "CDC45",
        "BUB1", "BUB1B", "MAD2L1", "AURKA", "AURKB", "PLK1", "PLK4",
        "E2F1", "E2F2", "RB1", "CHEK1", "CHEK2", "WEE1",
        # DNA replication (GO:0006260)
        "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "ORC1", "ORC6",
        "PCNA", "POLA1", "POLE", "POLE2", "RFC3", "RFC4", "RFC5",
        "RPA1", "RPA2", "PRIM1", "PRIM2", "FEN1", "LIG1",
        # p53 signaling (KEGG hsa04115)
        "TP53", "MDM2", "CDKN1A", "GADD45A", "GADD45B", "BAX", "BBC3",
        "PMAIP1", "SESN1", "SESN2", "RRM2B", "DDB2", "CCNG1",
        # Mitotic spindle / kinetochore
        "KIF11", "KIF23", "CENPA", "CENPF", "CENPE", "NDC80", "NUF2",
        "TPX2", "BIRC5", "TOP2A", "MKI67", "FOXM1", "MYBL2"
    )

    # Downregulated pathway genes (immune/inflammatory response, TNFa/NF-kB)
    down_genes <- c(
        # TNFa signaling via NF-kB (Hallmark / KEGG hsa04668)
        "TNF", "TNFAIP3", "NFKB1", "NFKB2", "RELA", "RELB", "NFKBIA",
        "NFKBIB", "TRAF1", "TRAF2", "TRAF5", "BIRC2", "BIRC3",
        "ICAM1", "VCAM1", "SELE", "CCL2", "CCL5", "CXCL1", "CXCL2",
        "CXCL8", "CXCL10", "IL6", "IL1B", "IL1A", "CSF2",
        "PTGS2", "SOD2", "IRF1", "STAT5A",
        # Inflammatory response
        "TLR2", "TLR4", "MYD88", "IRAK1", "IRAK4", "CD14", "LY96",
        "NLRP3", "CASP1", "IL18", "IL1R1", "IL1RAP", "RIPK2",
        # Complement / innate immunity
        "C3", "C1QA", "C1QB", "C1QC", "CFB", "CFD", "SERPING1",
        "CD55", "CD59", "CR1"
    )

    # ── Get background genes ─────────────────────────────────────────────
    all_symbols <- keys(org.Hs.eg.db, keytype="SYMBOL")
    pathway_genes <- unique(c(up_genes, down_genes))

    # Ensure pathway genes actually exist in the database
    up_genes   <- up_genes[up_genes %in% all_symbols]
    down_genes <- down_genes[down_genes %in% all_symbols]

    # Sample background genes (non-pathway)
    bg_pool <- setdiff(all_symbols, pathway_genes)
    n_bg <- 5000 - length(up_genes) - length(down_genes)
    bg_genes <- sample(bg_pool, min(n_bg, length(bg_pool)))

    all_gene_symbols <- c(up_genes, down_genes, bg_genes)
    n <- length(all_gene_symbols)
    n_up <- length(up_genes)
    n_down <- length(down_genes)
    n_bg_actual <- length(bg_genes)

    # ── Generate expression values ───────────────────────────────────────
    baseMean <- rlnorm(n, meanlog=6, sdlog=1.5)

    # Upregulated genes: strong positive fold change
    lfc_up   <- abs(rnorm(n_up, mean=2.5, sd=0.8))
    # Downregulated genes: strong negative fold change
    lfc_down <- -abs(rnorm(n_down, mean=2.0, sd=0.7))
    # Background genes: small random fold changes (mostly noise)
    lfc_bg   <- rnorm(n_bg_actual, mean=0, sd=0.3)

    log2FoldChange <- c(lfc_up, lfc_down, lfc_bg)

    # Generate p-values correlated with fold change magnitude
    # Pathway genes get very significant p-values; background gets high p-values
    pval_up   <- 10^(-runif(n_up, 4, 12))
    pval_down <- 10^(-runif(n_down, 3, 10))
    pval_bg   <- runif(n_bg_actual, 0.01, 1)

    pvalue <- c(pval_up, pval_down, pval_bg)
    padj   <- p.adjust(pvalue, method="BH")

    # Wald statistic and SE
    lfcSE <- abs(log2FoldChange / qnorm(pmax(1 - pvalue/2, 0.5001)))
    lfcSE[lfcSE < 0.01 | is.na(lfcSE) | is.infinite(lfcSE)] <- 0.1
    stat <- log2FoldChange / lfcSE

    data.table(
        gene_symbol    = all_gene_symbols,
        baseMean       = round(baseMean, 2),
        log2FoldChange = round(log2FoldChange, 4),
        lfcSE          = round(lfcSE, 4),
        stat           = round(stat, 4),
        pvalue         = signif(pvalue, 4),
        padj           = signif(padj, 4)
    )
}

# ─── Column Selector UI (custom upload) ─────────��────────────────────────
output$column_selectors <- renderUI({
    req(input$de_file)
    dat <- tryCatch(read_delim_auto(input$de_file$datapath), error=function(e) NULL)
    req(dat)
    cols <- colnames(dat)

    guess <- function(pattern, fallback_idx) {
        matched <- grep(pattern, cols, ignore.case=TRUE, value=TRUE)
        if (length(matched) > 0) matched[1] else cols[min(fallback_idx, length(cols))]
    }

    tagList(
        hr(),
        h4("Map Columns", style="color:#0F344C;"),
        selectInput("gene_col", "Gene ID column:", choices=cols,
            selected=cols[1]),
        selectInput("lfc_col", "log2FoldChange column:", choices=cols,
            selected=guess("log2F|logFC|lfc", 2)),
        selectInput("pval_col", "pvalue column:", choices=cols,
            selected=guess("^pval", 3)),
        selectInput("padj_col", "padj column:", choices=cols,
            selected=guess("padj|fdr|adj", 4)),
        selectInput("basemean_col", "baseMean column (optional):", choices=c("(none)", cols),
            selected=guess("basem|mean", 5)),
        hr()
    )
})

# ─── Input Reactive ────────���────────────────────────────────────���────────
InputReactive <- reactive({
    if (input$DemoData == TRUE) {
        shiny::validate(need(!is.null(input$de_file), "Please upload a DESeq2 results file."))
        tryCatch(
            read_delim_auto(input$de_file$datapath),
            error=function(e) {
                showNotification(paste("File parse error:", e$message), type='error', duration=NULL)
                NULL
            }
        )
    } else {
        cached <- DemoDataCache()
        if (is.null(cached)) {
            cached <- generate_demo_data()
            DemoDataCache(cached)
        }
        cached
    }
})

# ─── Preview Table ────���─────────────────────────────────────���────────────
output$input_table <- DT::renderDataTable({
    dat <- InputReactive()
    req(dat)
    DT::datatable(dat, style='bootstrap', options=list(pageLength=15, scrollX=TRUE))
})

# ─── Analysis Mode Reactive ─────────���────────────────────────────────────
AnalysisMode <- reactive({
    if (input$DemoData) input$analysis_mode else input$demo_mode
})

output$current_mode <- reactive({ AnalysisMode() })
outputOptions(output, "current_mode", suspendWhenHidden=FALSE)

# ─── Organism DB Reactive ────────────────────────────────────────────────
OrgDb <- reactive({
    org <- if (input$DemoData) input$organism else "hsapiens"
    if (org == "hsapiens") org.Hs.eg.db else org.Mm.eg.db
})

OrgKegg <- reactive({
    org <- if (input$DemoData) input$organism else "hsapiens"
    if (org == "hsapiens") "hsa" else "mmu"
})

OrgReactome <- reactive({
    org <- if (input$DemoData) input$organism else "hsapiens"
    if (org == "hsapiens") "human" else "mouse"
})

OrgSpeciesName <- reactive({
    org <- if (input$DemoData) input$organism else "hsapiens"
    if (org == "hsapiens") "Homo sapiens" else "Mus musculus"
})

# ─── Gene ID Type ───────────���────────────────────────────────────────────
GeneIdType <- reactive({
    if (input$DemoData) input$gene_id_type else "SYMBOL"
})

# ─── Processed DESeq2 Data ────────────��─────────────────────���───────────
ProcessedData <- reactive({
    dat <- InputReactive()
    req(dat)

    if (input$DemoData) {
        gene_col <- input$gene_col
        lfc_col  <- input$lfc_col
        pval_col <- input$pval_col
        padj_col <- input$padj_col
        bm_col   <- input$basemean_col
    } else {
        gene_col <- "gene_symbol"
        lfc_col  <- "log2FoldChange"
        pval_col <- "pvalue"
        padj_col <- "padj"
        bm_col   <- "baseMean"
    }

    df <- data.frame(
        gene       = as.character(dat[[gene_col]]),
        log2FC     = as.numeric(dat[[lfc_col]]),
        pvalue     = as.numeric(dat[[pval_col]]),
        padj       = as.numeric(dat[[padj_col]]),
        stringsAsFactors = FALSE
    )
    if (bm_col != "(none)" && bm_col %in% colnames(dat)) {
        df$baseMean <- as.numeric(dat[[bm_col]])
    }

    # Include stat column if available
    if ("stat" %in% colnames(dat)) {
        df$stat <- as.numeric(dat[["stat"]])
    }

    df <- df[!is.na(df$log2FC) & !is.na(df$pvalue), ]

    # Convert gene IDs to ENTREZID
    id_type <- GeneIdType()
    orgdb   <- OrgDb()

    if (id_type != "ENTREZID") {
        withProgress(message='Mapping gene identifiers...', value=0.5, {
            mapping <- tryCatch(
                AnnotationDbi::select(orgdb,
                    keys=unique(df$gene),
                    keytype=id_type,
                    columns="ENTREZID"),
                error=function(e) {
                    showNotification(paste("Gene ID mapping error:", e$message), type='error', duration=NULL)
                    NULL
                }
            )
        })
        if (!is.null(mapping)) {
            mapping <- mapping[!is.na(mapping$ENTREZID) & !duplicated(mapping[[id_type]]), ]
            df <- merge(df, mapping, by.x="gene", by.y=id_type, all.x=FALSE)
        }
    } else {
        df$ENTREZID <- df$gene
    }

    df <- df[!is.na(df$ENTREZID) & df$ENTREZID != "", ]
    df <- df[!duplicated(df$ENTREZID), ]
    df
})

# ─── Gene Summary Panel ──────────────────────────────────────────────────
output$gene_summary_panel <- renderUI({
    df <- ProcessedData()
    req(df)
    n_total <- nrow(df)
    n_up    <- sum(df$log2FC > 0 & !is.na(df$padj) & df$padj < 0.05, na.rm=TRUE)
    n_down  <- sum(df$log2FC < 0 & !is.na(df$padj) & df$padj < 0.05, na.rm=TRUE)
    n_ns    <- n_total - n_up - n_down

    tagList(
        h4("Gene Summary", style="color:#0F344C;"),
        div(class="control-group-panel",
            p(strong("Total genes with Entrez mapping:"), n_total),
            p(strong("Significant UP (padj < 0.05, LFC > 0):"), n_up, style="color:#d73027;"),
            p(strong("Significant DOWN (padj < 0.05, LFC < 0):"), n_down, style="color:#4575b4;"),
            p(strong("Not significant:"), n_ns, style="color:#999;")
        )
    )
})

# ─── ORA Gene Count Display ──────────────────────────────────────────────
output$ora_gene_count_display <- renderUI({
    df <- ProcessedData()
    req(df)
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

    div(class="control-group-panel",
        p(strong(nrow(sig)), " genes selected from ", strong(nrow(df)), " total",
          style="color:#0F344C; font-size:13px;")
    )
})
