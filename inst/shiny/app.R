# autoqmsp point-and-click app.
# Run with autoqmsp::run_app(), or deploy this folder to shinyapps.io /
# Posit Connect Cloud (both have free tiers).

library(shiny)
library(autoqmsp)

options(shiny.maxRequestSize = 200 * 1024^2)

default_reference <- "ACTB|B.?ACTIN|BETA.?ACTIN"
default_ntc <- "NTC|dH2O|dH20|water|blank|^NK"
default_positive <- "H460|A549|HT29|positive|^PC\\b"

css <- "
  body { background:#f7f7f5; }
  .card { background:#fff; border-radius:8px; padding:16px 20px; margin-bottom:16px;
          box-shadow:0 1px 2px rgba(0,0,0,.06); }
  .card h4 { margin-top:0; }
  .hint { color:#656d76; font-size:13px; }
  .upload-box .form-group { margin-bottom:0; }
  .upload-box .btn-file { font-size:16px; padding:10px 18px; }
  .step-next { margin-top:8px; font-weight:600; }
  .navbar-brand { font-weight:700; }
  iframe.report { width:100%; border:0; min-height:600px; background:#fff; }
  .code-box pre { background:#f3f3f0; }
"

ui <- navbarPage(
  title = "autoqmsp",
  id = "tabs",
  header = tags$head(tags$style(HTML(css))),

  tabPanel(
    "1. Upload", value = "upload",
    fluidRow(column(
      8, offset = 2,
      div(class = "card upload-box",
          h3("Upload your qPCR run files"),
          p(class = "hint",
            "Select one or more QuantStudio / Applied Biosystems .eds files ",
            "(the files saved after the run). Files stay on this server only ",
            "for this session."),
          fileInput("files", NULL, multiple = TRUE, accept = ".eds",
                    buttonLabel = "Choose .eds files...",
                    placeholder = "No files selected", width = "100%")),
      uiOutput("upload_summary")
    ))
  ),

  tabPanel(
    "2. Settings", value = "settings",
    fluidRow(
      column(6,
        div(class = "card",
            h4("When is a gene methylated?"),
            p("A gene is ", strong("methylated"), " when its ", strong("Ct is 40 or less"),
              ". The Ct is read where the amplification curve crosses ",
              strong(HTML("&Delta;Rn = 10,000")), " (for every gene, including the reference gene)."),
            p(class = "hint",
              "These are the lab's fixed settings. A gene is 'not determined' when ",
              "the sample's reference gene has no Ct of 40 or less, the gene's ",
              "water control amplified (for a methylated result), its positive ",
              "control did not amplify (for an unmethylated result), or the ",
              "replicates disagree.")),
        div(class = "card",
            h4("Cancer risk decision"),
            numericInput("min_genes",
                         "A sample is 'Potential cancer' when at least this many genes are methylated",
                         value = 1, min = 1, step = 1),
            uiOutput("panel_input"))
      ),
      column(6,
        div(class = "card",
            h4("Controls"),
            uiOutput("control_inputs")),
        div(class = "card",
            checkboxInput("advanced", strong("Show advanced settings")),
            conditionalPanel(
              "input.advanced",
              numericInput("ref_ct_warn",
                           "Reference gene: warn 'low DNA input' if Ct is above",
                           35, min = 20, max = 50, step = 0.5)
            ))
      )
    ),
    fluidRow(column(12, p(class = "step-next",
                          "Results update automatically. Open the Report tab.")))
  ),

  tabPanel(
    "3. Report", value = "report",
    div(class = "card",
        fluidRow(
          column(8, uiOutput("report_status")),
          column(4, align = "right",
                 downloadButton("download_report", "Download report (HTML)"),
                 downloadButton("download_excel", "Download Excel"))
        )),
    uiOutput("report_frame"),
    div(class = "card code-box",
        checkboxInput("show_code", "Show the R code for this analysis"),
        conditionalPanel("input.show_code", verbatimTextOutput("r_code")))
  ),

  tabPanel(
    "4. Details", value = "details",
    tabsetPanel(
      tabPanel("Heatmap",
               radioButtons("heat_value", NULL, inline = TRUE,
                            choices = c("Call" = "call", "Ct" = "ct")),
               plotOutput("heatmap", height = "650px")),
      tabPanel("All results", tableOutput("results")),
      tabPanel("Amplification curves",
               fluidRow(
                 column(4, selectInput("curve_run", "Run", NULL)),
                 column(8, selectInput("curve_target", "Genes (empty = all)",
                                       NULL, multiple = TRUE))),
               plotOutput("curves", height = "650px")),
      tabPanel("Plate layout",
               selectInput("plate_run", "Run", NULL),
               plotOutput("plate", height = "550px"))
    )
  )
)

server <- function(input, output, session) {
  runs <- reactive({
    req(input$files)
    read_eds_files(input$files$datapath,
                   names = tools::file_path_sans_ext(input$files$name))
  })

  genes <- reactive(sort(unique(runs()$wells$target)))
  samples <- reactive(sort(unique(runs()$wells$sample)))
  reference <- reactive({
    if (is.null(input$reference) || input$reference == "(none)") NULL
    else input$reference
  })
  panel_genes <- reactive(setdiff(genes(), reference()))

  # ---- upload ----------------------------------------------------------------
  output$upload_summary <- renderUI({
    if (is.null(input$files)) return(NULL)
    x <- tryCatch(runs(), error = function(e) e)
    if (inherits(x, "error")) {
      return(div(class = "card", style = "border-left:4px solid #c0392b",
                 strong("Could not read the files: "), conditionMessage(x)))
    }
    rows <- lapply(unique(x$wells$run), function(r) {
      w <- x$wells[x$wells$run == r, ]
      tags$tr(tags$td(r), tags$td(length(unique(w$sample))),
              tags$td(paste(unique(w$target), collapse = ", ")),
              tags$td(nrow(w)))
    })
    div(class = "card",
        h4(icon("check"), " ", length(rows), " run(s) loaded"),
        tags$table(class = "table table-condensed",
                   tags$thead(tags$tr(tags$th("Run"), tags$th("Samples"),
                                      tags$th("Genes"), tags$th("Wells"))),
                   tags$tbody(rows)),
        actionButton("go_settings", "Next: check the settings",
                     class = "btn-primary"))
  })
  observeEvent(input$go_settings, updateNavbarPage(session, "tabs", "settings"))

  # ---- settings ---------------------------------------------------------------
  output$control_inputs <- renderUI({
    if (is.null(input$files)) return(p(class = "hint", "Upload files first."))
    g <- genes()
    s <- samples()
    ref_default <- g[grepl(default_reference, g, ignore.case = TRUE)][1]
    tagList(
      selectInput("reference", "Reference gene (e.g. beta-actin)",
                  c("(none)", g),
                  selected = if (is.na(ref_default)) "(none)" else ref_default),
      selectizeInput("ntc_samples", "No-template (water) controls", s,
                     multiple = TRUE,
                     selected = s[grepl(default_ntc, s, ignore.case = TRUE,
                                        perl = TRUE)]),
      selectizeInput("pos_samples", "Positive (methylated) controls", s,
                     multiple = TRUE,
                     selected = s[grepl(default_positive, s, ignore.case = TRUE,
                                        perl = TRUE)]),
      p(class = "hint", "Controls were detected from the sample names. ",
        "Add or remove samples if needed.")
    )
  })

  output$panel_input <- renderUI({
    if (is.null(input$files)) return(NULL)
    checkboxGroupInput("panel", "Genes that count for the decision",
                       panel_genes(), selected = panel_genes(), inline = TRUE)
  })

  exact <- function(x) {
    if (!length(x)) return("(?!)")   # matches nothing
    paste0("^(", paste(gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x),
                       collapse = "|"), ")$")
  }

  settings <- reactive({
    req(input$files, input$reference)
    list(
      reference = if (is.null(reference())) NULL else exact(reference()),
      ntc = exact(input$ntc_samples),
      positive = exact(input$pos_samples),
      ct_cutoff = 40, threshold = 10000, ref_ct_max = 40,
      ref_ct_warn = input$ref_ct_warn,
      min_methylated_genes = input$min_genes,
      panel = input$panel
    )
  })

  result <- reactive({
    s <- settings()
    do.call(analyze_qmsp, c(list(runs()), s))
  })

  # ---- report -----------------------------------------------------------------
  output$report_status <- renderUI({
    if (is.null(input$files)) {
      return(p("Upload .eds files in the ", strong("Upload"), " tab first."))
    }
    r <- result()$report
    v <- table(r$verdict)
    p(strong(nrow(r), " samples: "),
      paste(names(v), v, sep = " ", collapse = " · "))
  })

  output$report_frame <- renderUI({
    req(input$files)
    tags$iframe(
      class = "report",
      srcdoc = autoqmsp:::report_html(result()),
      onload = "this.style.height = (this.contentWindow.document.body.scrollHeight + 40) + 'px';"
    )
  })

  output$download_report <- downloadHandler(
    filename = function() paste0("qmsp_report_", Sys.Date(), ".html"),
    content = function(file) write_report(result(), file)
  )
  output$download_excel <- downloadHandler(
    filename = function() paste0("qmsp_results_", Sys.Date(), ".xlsx"),
    content = function(file) export_results(result(), file)
  )

  output$r_code <- renderText({
    req(input$files)
    s <- settings()
    arg <- function(v) paste(deparse(v, width.cutoff = 500), collapse = "")
    args <- vapply(names(s), function(n) sprintf("  %s = %s", n, arg(s[[n]])),
                   character(1))
    paste0(
      "library(autoqmsp)\n\n",
      "runs <- read_eds_files(", arg(input$files$name), ")\n\n",
      "res <- analyze_qmsp(\n  runs,\n", paste(args, collapse = ",\n"), "\n)\n\n",
      "res$report                                # verdict per sample\n",
      "write_report(res, \"qmsp_report.html\")\n",
      "export_results(res, \"qmsp_results.xlsx\")\n"
    )
  })

  # ---- details ----------------------------------------------------------------
  observeEvent(runs(), {
    rn <- unique(runs()$wells$run)
    updateSelectInput(session, "curve_run", choices = rn)
    updateSelectInput(session, "plate_run", choices = rn)
  })
  observeEvent(input$curve_run, {
    w <- runs()$wells
    updateSelectInput(session, "curve_target",
                      choices = unique(w$target[w$run == input$curve_run]))
  })

  output$heatmap <- renderPlot(plot_methylation(result(), input$heat_value))
  output$results <- renderTable({
    r <- result()$results
    r$call <- as.character(r$call)
    r[, c("run", "sample", "role", "target", "ct", "ref_ct", "call", "notes")]
  }, digits = 2)
  output$curves <- renderPlot({
    req(input$curve_run)
    plot_amplification(result(), run = input$curve_run,
                       targets = input$curve_target)
  })
  output$plate <- renderPlot({
    req(input$plate_run)
    plot_plate(result(), run = input$plate_run, fill = "ct")
  })
}

shinyApp(ui, server)
