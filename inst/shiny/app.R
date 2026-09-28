# autoqmsp point-and-click app.
# Run with autoqmsp::run_app(), or deploy this folder to shinyapps.io /
# Posit Connect Cloud (both have free tiers).

library(shiny)
library(autoqmsp)

options(shiny.maxRequestSize = 100 * 1024^2)

ui <- fluidPage(
  titlePanel("autoqmsp: qMSP analysis from QuantStudio .eds files"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      fileInput("files", "1. Upload .eds file(s)", multiple = TRUE,
                accept = ".eds"),
      h5("2. Settings"),
      textInput("reference", "Reference gene (regex, empty = none)",
                "ACTB|B.?ACTIN|BETA.?ACTIN"),
      textInput("ntc", "No-template control names (regex)",
                "NTC|dH2O|dH20|water|blank|^NK"),
      textInput("positive", "Positive control names (regex)",
                "H460|A549|HT29|positive|^PC\\b"),
      numericInput("ct_cutoff", "Gene Ct cutoff (methylated if <=)", 40,
                   min = 20, max = 50, step = 0.5),
      numericInput("ref_ct_max", "Reference Ct max (invalid if >)", 40,
                   min = 20, max = 50, step = 0.5),
      numericInput("ref_ct_warn", "Reference Ct warning (low input if >)", 35,
                   min = 20, max = 50, step = 0.5),
      numericInput("min_cq_conf", "Minimum Cq confidence", 0.5,
                   min = 0, max = 1, step = 0.05),
      numericInput("min_plateau",
                   "Minimum curve height (fraction of positive control)", 0.2,
                   min = 0, max = 1, step = 0.05),
      hr(),
      downloadButton("download", "3. Download Excel")
    ),
    mainPanel(
      width = 9,
      uiOutput("status"),
      tabsetPanel(
        tabPanel("Summary",
                 plotOutput("heatmap", height = "600px"),
                 tableOutput("summary")),
        tabPanel("Results", tableOutput("results")),
        tabPanel("Controls", tableOutput("controls")),
        tabPanel("Wells to review", tableOutput("review")),
        tabPanel("Curves",
                 fluidRow(
                   column(6, selectInput("curve_run", "Run", NULL)),
                   column(6, selectInput("curve_target", "Gene", NULL,
                                         multiple = TRUE))
                 ),
                 plotOutput("curves", height = "650px")),
        tabPanel("Plate",
                 selectInput("plate_run", "Run", NULL),
                 plotOutput("plate", height = "550px"))
      )
    )
  )
)

server <- function(input, output, session) {
  runs <- reactive({
    req(input$files)
    read_eds_files(input$files$datapath,
                   names = tools::file_path_sans_ext(input$files$name))
  })

  result <- reactive({
    ref <- if (nzchar(input$reference)) input$reference else NULL
    analyze_qmsp(runs(), reference = ref, ntc = input$ntc,
                 positive = input$positive, ct_cutoff = input$ct_cutoff,
                 ref_ct_max = input$ref_ct_max, ref_ct_warn = input$ref_ct_warn,
                 min_cq_conf = input$min_cq_conf,
                 min_plateau = input$min_plateau)
  })

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

  output$status <- renderUI({
    if (is.null(input$files)) {
      return(p("Upload one or more .eds files to start."))
    }
    r <- result()
    bad <- r$controls[r$controls$ntc_status %in% c("Fail", "Review") |
                        r$controls$positive_status == "Fail", ]
    tagList(
      if (nrow(bad)) div(
        style = "background:#fff3cd;padding:8px;border-radius:4px;",
        strong("Control problems: "),
        paste0(bad$run, " / ", bad$target, " (NTC ", bad$ntc_status,
               ", positive ", bad$positive_status, ")", collapse = "; ")
      ),
      p(sprintf("%d wells need manual review.", sum(r$wells$result == "Review")))
    )
  })

  output$heatmap <- renderPlot(plot_methylation(result()))
  output$summary <- renderTable(results_wide(result(), "call"))
  output$results <- renderTable({
    r <- result()$results
    r$call <- as.character(r$call)
    r[, c("run", "sample", "role", "target", "ct", "ref_ct", "ref_status",
          "delta_ct", "ratio", "pmr", "call", "notes")]
  }, digits = 2)
  output$controls <- renderTable(result()$controls, digits = 2)
  output$review <- renderTable({
    w <- result()$wells
    w[w$result == "Review",
      c("run", "well", "sample", "target", "ct", "cq_conf", "flags")]
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

  output$download <- downloadHandler(
    filename = function() paste0("qmsp_results_", Sys.Date(), ".xlsx"),
    content = function(file) export_results(result(), file)
  )
}

shinyApp(ui, server)
