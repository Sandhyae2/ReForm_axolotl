library(shiny)
library(bslib)

ui <- page_sidebar(
  title = div(
  style = "display:flex;align-items:center;gap:10px;",
  tags$img(
    src = "axolotl.png",
    style = "width:48px;height:48px;object-fit:contain;border-radius:50%;"
  ),
  div(
    span(
      "ReForm",
      style = "font-weight:850;letter-spacing:1px;color:white;font-size:25px;"
    ),
    br(),
    span(
      "The science of becoming whole again.",
      style = "font-size:12px;color:#b9e5dc;"
    )
  )
),

  theme = bs_theme(
    version = 5,
    bg = "#153d35", ###"#f5f8f5",
    fg = "#183d35",
    primary = "#287c6b",
    secondary = "#b9e5dc",
    base_font = font_google("DM Sans"),
    heading_font = font_google("DM Sans")
  ),

  sidebar = sidebar(
    width = 263,
    bg = "#153d35",
    fg = "#ffffff",

    div(
      style = "padding:10px 4px 20px;",
      p(
        "AXOLOTL LIMB REGENERATION",
        style = paste(
          "color:#8fd9c7;font-size:10px;",
          "letter-spacing:2px;font-weight:800;"
        )
      ),
      p(
        "Explore the molecular story behind regeneration.",
        style = "color:#d0e8df;font-size:13px;"
      )
    ),

    navset_pill_list(
      id = "page",
      widths = c(12, 12),

      nav_panel(
        title = "Home",
        value = "home",
        icon = icon("house")
      ),
      nav_panel(
        title = "Differential Expression",
        value = "de",
        icon = icon("chart-bar")
      ),
      nav_panel(
        title = "Gene Expression",
        value = "genes",
        icon = icon("dna")
      ),
      nav_panel(
        title = "Heatmap & PCA",
        value = "heatmap",
        icon = icon("th")
      ),
      nav_panel(
        title = "Functional Enrichment",
        value = "enrichment",
        icon = icon("project-diagram")
      ),
      nav_panel(
        title = "DEG Explorer",
        value = "deg",
        icon = icon("search")
      ),
      nav_panel(
        title = "About the Project",
        value = "about",
        icon = icon("book-open")
      )
    ),

    br(),
    p(
      "ReForm | RNA-seq exploration",
      style = "font-size:11px;color:#b9e5dc;"
    )
  ),

  tags$head(
    tags$style(HTML("
      body {
        background: #f5f8f5;
      }
	
	body {
        background: #f5f8f5;
      }

      /* Dark teal background only behind the sidebar-width header area */
      @media (min-width: 768px) {
        .bslib-page-sidebar > .navbar > .container-fluid {
          background: transparent;
        }
      }

      .hero {
        background: linear-gradient(120deg, #153d35, #287c6b);
        color: white;
        padding: 36px;
        border-radius: 22px;
        margin-bottom: 28px;
      }

      .hero {
        background: linear-gradient(120deg, #153d35, #287c6b);
        color: white;
        padding: 36px;
        border-radius: 22px;
        margin-bottom: 28px;
      }

      .hero h1 {
        font-size: clamp(30px, 4vw, 48px);
        font-weight: 850;
        letter-spacing: -1.5px;
      }

      .hero p {
        color: #d0e8df;
        line-height: 1.8;
        max-width: 720px;
      }

      .eyebrow {
        color: #8fd9c7;
        font-size: 11px;
        font-weight: 800;
        letter-spacing: 2px;
      }

      .section-title {
        font-weight: 800;
        margin-top: 34px;
        margin-bottom: 8px;
        color: #183d35;
      }

      .info-card, .stage-card, .candidate, .workflow-card {
        background: white;
        border: 1px solid #e1ebe5;
        border-radius: 16px;
        padding: 20px;
        height: 100%;
        transition: transform .2s ease, box-shadow .2s ease;
      }

      .info-card:hover, .stage-card:hover,
      .candidate:hover, .workflow-card:hover {
        transform: translateY(-4px);
        box-shadow: 0 10px 25px rgba(24, 61, 53, .10);
      }

      .card-symbol {
        font-size: 26px;
        color: #287c6b;
        margin-bottom: 12px;
      }

      .info-card h4, .stage-card h4,
      .workflow-card h4, .candidate strong {
        font-weight: 800;
        color: #183d35;
      }

      .info-card p, .stage-card p,
      .workflow-card p, .candidate p {
        color: #71837c;
        font-size: 13px;
        line-height: 1.7;
      }

      .stage-number {
        display: inline-flex;
        align-items: center;
        justify-content: center;
        width: 34px;
        height: 34px;
        border-radius: 50%;
        background: #e2f4ed;
        color: #287c6b;
        font-weight: 800;
        margin-bottom: 12px;
      }

      .focus-panel {
        background: #e8f3ed;
        border-left: 5px solid #287c6b;
        border-radius: 14px;
        padding: 25px;
        margin-top: 12px;
      }

      .focus-panel h4 {
        font-weight: 800;
      }

      .focus-panel p {
        color: #526b61;
        line-height: 1.8;
      }

      .footer-band {
        background: #153d35;
        color: #d0e8df;
        padding: 18px 22px;
        border-radius: 14px;
        margin-top: 28px;
        font-size: 13px;
      }

      .small-label {
        font-size: 10px;
        letter-spacing: 1.5px;
        font-weight: 800;
        color: #287c6b;
      }

      @media (prefers-reduced-motion: reduce) {
        .info-card, .stage-card,
        .candidate, .workflow-card {
          transition: none;
        }
      }
    "))
  ),

  uiOutput("home_content")
)

server <- function(input, output, session) {

  output$home_content <- renderUI({

    if (input$page != "home") {
      return(
        div(
          class = "focus-panel",
          h3("This section is coming next"),
          p(
            "The homepage is ready. The data tables, ",
            "differential-expression results and plots ",
            "will be connected in the next stage."
          )
        )
      )
    }

    tagList(

      # HERO
      div(
        class = "hero",
        p(
          class = "eyebrow",
          "AXOLOTL LIMB REGENERATION .RNA-seq"
        ),
        h1("From amputation to regeneration"),
        br(),
        p(
          paste(
            "The axolotl(Ambystoma mexicanum) can regenerate lost limbs,
		 making it a powerful model to study the molecular and cellular
             mechanisms of regeneration."
          )
        ),
        actionButton(
          "explore_data",
          "Explore the data",
          icon = icon("arrow-right"),
          class = "btn-light"
        )
      ),

      p(
        "Navigate the data and analysis sections.",
        style = "color:#71837c;"
      ),

      layout_columns(
        div(
          class = "workflow-card",
          div(class = "card-symbol", icon("table")),
          h4("Metadata"),
          p("Explore sample information and experimental groups."),
          actionButton(
            "go_metadata",
            "Explore metadata",
            class = "btn-outline-success btn-sm"
          )
        ),

        div(
          class = "workflow-card",
          div(class = "card-symbol", icon("database")),
          h4("Count matrix"),
          p("Inspect the gene-expression count data."),
          actionButton(
            "go_counts",
            "View counts",
            class = "btn-outline-success btn-sm"
          )
        ),

        div(
          class = "workflow-card",
          div(class = "card-symbol", icon("chart-bar")),
          h4("DESeq2 results"),
          p("Explore differential-expression comparisons."),
          actionButton(
            "go_de",
            "Explore results",
            class = "btn-outline-success btn-sm"
          )
        ),

        div(
          class = "workflow-card",
          div(class = "card-symbol", icon("chart-line")),
          h4("Visualisations"),
          p("Explore expression patterns and analysis plots."),
          actionButton(
            "go_plots",
            "Explore plots",
            class = "btn-outline-success btn-sm"
          )
        ),

        col_widths = c(6, 6, 6, 6)
      ),

      # CANDIDATE GENES
      h3(
        class = "section-title",
        "Six candidate genes"
      ),
      p(
        "Six genes selected for further exploration in the regeneration project.",
        style = "color:#71837c;"
      ),

      layout_columns(
        div(
          class = "candidate",
          strong("FLRT3"),
          p("Cell interactions and adhesion")
        ),

        div(
          class = "candidate",
          strong("MDK"),
          p("Growth signalling and cell responses")
        ),

        div(
          class = "candidate",
          strong("MMP11"),
          p("Extracellular matrix remodelling")
        ),

        div(
          class = "candidate",
          strong("SALL4"),
          p("Developmental regulation")
        ),

        div(
          class = "candidate",
          strong("TGFB2"),
          p("TGF-β signalling")
        ),

        div(
          class = "candidate",
          strong("VTN"),
          p("Cell adhesion and extracellular matrix")
        ),

        col_widths = c(4, 4, 4, 4, 4, 4)
      ),

      # PROJECT FOCUS
      h3(
        class = "section-title",
        "The focus of this project"
      ),

      div(
        class = "focus-panel",
        h4("What changes during regeneration?"),
        p(
          paste(
            "This project investigates gene-expression",
            "changes across the 3 and 5 days",
            "post-amputation stages, with the aim of",
            "identifying molecular patterns and candidate",
            "genes associated with regeneration."
          )
        ),
        p(
          paste(
            "The next sections will let you explore",
            "the metadata, count matrix, differential-expression",
            "results and visualisations interactively."
          )
        )
      ),

      div(
        class = "footer-band",
        icon("leaf"),
        " Explore the data. Investigate the genes. Understand regeneration."
      )
    )
  })

  # NAVIGATION BUTTONS
  observeEvent(input$explore_data, {
    updateNavsetPillList(
      session, "page", selected = "de"
    )
  })

  observeEvent(input$go_metadata, {
    updateNavsetPillList(
      session, "page", selected = "about"
    )
  })

  observeEvent(input$go_counts, {
    updateNavsetPillList(
      session, "page", selected = "genes"
    )
  })

  observeEvent(input$go_de, {
    updateNavsetPillList(
      session, "page", selected = "de"
    )
  })

  observeEvent(input$go_plots, {
    updateNavsetPillList(
      session, "page", selected = "heatmap"
    )
  })
}

shinyApp(ui = ui, server = server)
