# Required packages ####

if (!require("shiny", quietly = TRUE)) {
  install.packages("shiny")
}

if (!require("bslib", quietly = TRUE)) {
  install.packages("bslib")
}

if (!require("DT", quietly = TRUE)) {
  install.packages("DT")
}

if (!require("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

if (!require("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}

if (!require("plotly", quietly = TRUE)) {
  install.packages("plotly")
}

if (!require("rsconnect", quietly = TRUE)) {
  install.packages("rsconnect")
}
if (!require("readr", quietly = TRUE)) {
  install.packages("readr")
}
if (!require("zip", quietly = TRUE)) {
  install.packages("zip")
}

library(readr)
library(shiny)
library(bslib)
library(DT)
library(dplyr)
library(tidyr)
library(plotly)
library(rsconnect)


# Load & sanitize weights ####
pesi_path <- ("pesi_2.csv")
raw_pesi <- tryCatch({
  readr::read_delim(
    pesi_path,
    delim = NULL,        # AUTO-DETECT separator
    locale = locale(decimal_mark = "."), 
    show_col_types = FALSE,
    progress = FALSE
  ) |> as.data.frame()
}, error = function(e) {
  stop("Cannot read pesi file: ", e$message)
})
names(raw_pesi) <- gsub("﻿", "", names(raw_pesi))
names(raw_pesi) <- trimws(names(raw_pesi))

pesi <- raw_pesi
pesi[] <- lapply(pesi, function(col) {
  if (is.character(col)) {
    tmp <- gsub(",", ".", col, fixed = TRUE)
    num <- suppressWarnings(as.numeric(tmp))
    if (!all(is.na(num))) return(num) else return(col)
  } else {
    return(col)
  }
})

expected_inds <- c("GS","CS","DW","RN","SP","DB","RW","SLC","FD","IF")
expected_cols <- c(expected_inds, "pagg_c")

pesi_cols <- names(pesi)
missing_cols <- setdiff(expected_cols, pesi_cols)
extra_cols <- setdiff(pesi_cols, expected_cols)
if (length(missing_cols) > 0) {
  for (mc in missing_cols) pesi[[mc]] <- 0
  message("Weights file missing columns: ", paste(missing_cols, collapse = ", "),
          ". These columns were added with 0 values temporarily. Update pesi_2.csv to fix permanently.")
}
if (length(extra_cols) > 0) {
  message("Weights file contains extra columns (ignored by calc): ", paste(extra_cols, collapse = ", "))
}
pesi <- pesi[, intersect(c(expected_cols, names(pesi)), names(pesi)), drop = FALSE]

# Helper to fetch weight (returns numeric or 0)

get_weight_or_zero <- function(pesi_df, row_idx, col_name) {
  if (!(col_name %in% names(pesi_df))) return(0)
  val <- pesi_df[[col_name]][row_idx]
  if (is.list(val)) { val <- val[[1]] }
  num <- suppressWarnings(as.numeric(val))
  if (is.na(num)) return(0)
  return(num)
}


# Normalization functions ####

normalize_positive <- function(x, na.rm = TRUE) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  max_x <- max(x, na.rm = na.rm)
  min_x <- min(x, na.rm = na.rm)
  if (is.na(max_x) || is.na(min_x) || max_x == min_x) return(rep(0, length(x)))
  (x - min_x) / (max_x - min_x)
}
normalize_negative <- function(x, na.rm = TRUE) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  max_x <- max(x, na.rm = na.rm)
  min_x <- min(x, na.rm = na.rm)
  if (is.na(max_x) || is.na(min_x) || max_x == min_x) return(rep(0, length(x)))
  (x - max_x) / (min_x - max_x)
}

# calculation of smartness indeces function

calc_smartness <- function(df, pesi_df) {
  positive_indicators <- c("GS", "CS", "DW", "RN", "SP", "DB", "RW")
  negative_indicators <- c("SLC", "FD", "IF")
  all_possible <- c(positive_indicators, negative_indicators)
  all_indicators <- intersect(colnames(df), all_possible)
  if (length(all_indicators) == 0) {
    stop("No recognized indicator columns found in uploaded CSV.")
  }
  for (ind in all_indicators) {
    if (ind %in% positive_indicators) {
      df[[paste0(ind, "_norm")]] <- normalize_positive(df[[ind]])
    } else {
      df[[paste0(ind, "_norm")]] <- normalize_negative(df[[ind]])
    }
  }
  results <- lapply(seq_len(nrow(df)), function(i) {
    row <- df[i, , drop = FALSE]
    norm_cols <- grep("_norm$", names(row), value = TRUE)
    base_names <- gsub("_norm$", "", norm_cols)
    vals <- as.vector(unlist(row[norm_cols]))
    valid_inds <- base_names[!is.na(vals)]
    if (length(valid_inds) == 0) {
      return(data.frame(csf = NA_real_, csf_mit = NA_real_, csf_adp = NA_real_, csf_sd = NA_real_))
    }
    mit <- 0; adp <- 0; sd <- 0
    for (ind in valid_inds) {
      norm_val <- as.numeric(row[[paste0(ind, "_norm")]])
      w_mit <- get_weight_or_zero(pesi_df, 1, ind)
      w_adp <- get_weight_or_zero(pesi_df, 2, ind)
      w_sd  <- get_weight_or_zero(pesi_df, 3, ind)
      mit <- mit + norm_val * w_mit
      adp <- adp + norm_val * w_adp
      sd  <- sd  + norm_val * w_sd
    }
    # Weights of the three dimensions
    pagg_mit <- get_weight_or_zero(pesi_df, 1, "pagg_c")
    pagg_adp <- get_weight_or_zero(pesi_df, 2, "pagg_c")
    pagg_sd  <- get_weight_or_zero(pesi_df, 3, "pagg_c")
    
    # Sum of weights (used for normalization)
    tot_weight <- pagg_mit + pagg_adp + pagg_sd
    
    # Weighted scores
    csf_mit <- mit * pagg_mit
    csf_adp <- adp * pagg_adp
    csf_sd  <- sd * pagg_sd
    
    csf_total <- (csf_mit + csf_adp + csf_sd)/tot_weight
    
    data.frame(
      csf = round(as.numeric(csf_total), 2),
      csf_mit = round(as.numeric(mit), 2),
      csf_adp = round(as.numeric(adp), 2),
      csf_sd = round(as.numeric(sd), 2)
    )
  })
  do.call(rbind, results)
}

# UI ####
ui <- fluidPage(
  theme = bs_theme(
    bootswatch = "flatly",
    primary = "#003D39", 
    secondary = "#6FBC85",
    base_font = font_google("Open Sans"),
    heading_font = font_google("Merriweather")
  ),
  # theme style 
  tags$head(
    tags$style(HTML("
    :root{
      --deep-green: #003D39;
      --bright-green: #6FBC85;
      --fauna-teal: #00A19A;
      --stream-blue: #1488CA;
      --sun-yellow: #FFED99;
      --white: #FFFFFF;
    }

    body {
      background: linear-gradient(180deg, var(--deep-green), #002f2d) !important;
      color: var(--white) !important;
    }

    /* Title: very explicit selectors + high specificity */
    .app-title,
    .title-panel h1, 
    .title-panel .shiny-title,
    .shiny-title,
    h1.title,
    .navbar-brand,
    .bslib-navbar-brand {
      color: var(--white) !important;
      font-weight: 700 !important;
      text-shadow: 0 1px 2px rgba(0,0,0,0.45) !important;
    }

    /* Tabs */
    .nav-tabs .nav-link {
      color: var(--fauna-teal) !important;
      font-weight: 600;
    }
    .nav-tabs .nav-link.active {
      background-color: rgba(255,255,255,0.12) !important;
      color: var(--white) !important;
    }

    /* Left control panel */
    .left-panel {
      background-color: var(--bright-green) !important;
      color: var(--deep-green) !important;
      padding: 16px;
      border-radius: 10px;
      box-shadow: 0 3px 10px rgba(0,0,0,0.18);
    }
    .left-panel h4, .left-panel strong { color: var(--deep-green) !important; }

    /* Dataset card */
    .dataset-card {
      background-color: var(--fauna-teal) !important;
      color: var(--white) !important;
      padding: 12px;
      border-radius: 8px;
      margin-bottom: 12px;
    }

    /* Smartness card */
    .smartness-card {
      background-color: var(--stream-blue) !important;
      color: var(--white) !important;
      padding: 10px;
      border-radius: 8px;
      margin-top: 8px;
    }

    /* Ternary description */
    .ternary-desc {
      background-color: var(--white) !important;
      color: var(--deep-green) !important;
      border-radius: 8px;
    }

    .inner-white {
      background-color: var(--white) !important;
      color: var(--deep-green) !important;
      padding: 8px;
      border-radius: 6px;
    }

    /* Buttons */
    #reset_all { background-color: var(--sun-yellow) !important; color: var(--deep-green) !important; border: none !important; }
    #download_all { background-color: var(--deep-green) !important; color: var(--white) !important; border: none !important; }

    /* Plot layout */
    .ternary-row { display: flex; gap: 12px; align-items: stretch; }
    .ternary-plot-panel { flex: 1 1 60%; }
    .ternary-desc-panel { flex: 1 1 40%; }
    .ternary-plot-panel .inner-white, .ternary-desc-panel .ternary-desc { height: 600px; overflow: hidden; }
    .ternary-desc-panel .ternary-desc { overflow-y: auto; padding: 12px; }

    /* Ensure Plotly background remains transparent */
    .plotly, .plot-container, .plotly .main-svg { background-color: transparent !important; }

    /* Custom title row with logos */
    .title-with-logos {
      display: flex;
      align-items: center;
      justify-content: center;
      gap: 20px;
      margin-bottom: 10px;
    }
    .title-with-logos img {
      height: 60px;
      border-radius: 8px;
    }
    /* Help text color (match main text color) */
    .help-block, .help-text, .form-text {
      color: var(--deep-green) !important;
      font-size: 13px !important;
    }

/* tabella */
.weights-table-wrapper table {
  min-width: 700px;
  background-color: transparent !important;
  color: var(--white) !important;
  border-collapse: collapse;
}

/* celle */
.weights-table-wrapper td,
.weights-table-wrapper th {
  background-color: transparent !important;
  border: 1px solid rgba(255,255,255,0.15) !important;
  font-size: 13px;
  white-space: nowrap;
}

/* header */
.weights-table-wrapper th {
  font-weight: 600;
}

/* riga criteri (prima riga rbind) */
.weights-table-wrapper tr:first-child {
  font-style: italic;
  background-color: rgba(255,255,255,0.05);
}

.weights-table-wrapper {
  width: 100%;
  max-width: 100%;
  overflow-x: auto;
  display: block;
}
}"))
  ),
  
  titlePanel(
    div(
      style = "display:flex; align-items:center; justify-content:space-between; gap:20px; width:100%;",
      # Left spacer to help centering the title
      div(style = "flex: 1 1 20%;"),
      # Center title
      div(style = "flex: 0 1 60%; display:flex; justify-content:center; align-items:center;",
          h1("Climate-Smart Forestry Assessment Tool (CSF-AT App)", class = "app-title", style = "margin:0; text-align:center;")
      ),
      # Right: column with two logos side by side and AEDIT below, more right-aligned
      div(style = "flex: 0 0 20%; display:flex; flex-direction:column; align-items:flex-end; justify-content:center; gap:8px; padding-left:8px;",
          # top row: Unimol + Forwards (side by side, right aligned)
          div(style = "display:flex; align-items:center; gap:10px; justify-content:flex-end; width:100%;",
              img(src = "https://www3.unimol.it/assets/images/unimol/images/header/unimol_on.svg",
                  alt = "Unimol logo",
                  style = "max-height:80px; height:auto; width:auto; display:block;"),
              img(src = "https://forwards-project.eu/wp-content/themes/forwards_theme/images/FWD-LOGO-RGB_MAIN_COLOR.svg",
                  alt = "Forwards logo",
                  style = "max-height:80px; height:auto; width:auto; display:block;")
          ),
          # bottom row: AEDIT (smaller, right-aligned under the two)
          div(style = "display:flex; justify-content:flex-end; width:100%;",
              img(src = "https://www.aedit.it/wp-content/themes/aedit-theme/img/logo-aedit.png",
                  alt = "AEDIT logo",
                  style = "max-height:40px; height:auto; width:auto; display:block;")
          )
      )
    )
  ),
  
  tabsetPanel(
    # App description ####
    tabPanel("ℹ️ App Description",
             
             tags$style(HTML("
    .info-card {
      background-color: rgba(255,255,255,0.05);
      border-radius: 10px;
      padding: 18px;
      height: 240px; /* altezza uguale per tutte le card */
      display: flex;
      flex-direction: column;
      justify-content: flex-start;
      color: var(--white);
    }
    
    /* top row: due card con altezza uguale */
    .info-card.top-card { height: 320px; }
    /* bottom row: due card con altezza uguale (diversa dalla top) */
    .info-card.bottom-card { height: 380px; }
    
    .info-card h4 {
      color: var(--sun-yellow);
      margin-top: 0;
      margin-bottom: 10px;
    }
    .info-card ul, .info-card ol {
      margin-top: 4px;
      padding-left: 18px;
    }
    
    /* assicurare che su schermi piccoli le card siano leggibili */
    @media (max-width: 767px) {
      .info-card.top-card, .info-card.bottom-card { height: auto; }
    }
    
  ")),
             
             fluidRow(
               # Row 1
               column(6,
                      tags$div(class = "info-card top-card",
                               tags$h4("What is Climate-Smart Forestry (CSF)?"),
                               tags$p("Climate-Smart Forestry (CSF) is an integrated approach to forest management that balances three complementary objectives:"),
                               tags$ul(
                                 tags$li(tags$b("Mitigation"), " – enhancing carbon stocks and reducing greenhouse gas emissions."),
                                 tags$li(tags$b("Adaptation"), " – improving forest resilience to climate change impacts."),
                                 tags$li(tags$b("Social dimension"), " – sustaining production and social benefits.")
                               )
                      )
               ),
               column(6,
                      tags$div(class = "info-card top-card",
                               tags$h4("Purpose of this app"),
                               tags$p("CSF-AT translates the CSF concept into a practical and reproducible assessment tool. The CSF-AT is an interactive and user-friendly application that allows users to assess the dregree of 'climate smartness' of forest stands. It provides a structured assessment of the three main pillars of CSF.")
                      )
               )
             ),
             
             fluidRow(style = "margin-top:12px;",
                      # Row 2
                      column(6,
                             tags$div(class = "info-card bottom-card",
                                      tags$h4("Method & data"),
                                      tags$p("The application computes a Composite Climate-Smart index (ICSF) from ecological and structural indicators through a hierarchical and analytical method, allowing consistent comparisons across stands, management scenarios, and time (Alfieri et al., 2024). The method consists of:"),
                                      tags$ol(
                                        tags$li(tags$b("Normalization"), " – scaling indicators to a 0–1 range."),
                                        tags$li(tags$b("Weighting"), " – applying expert-derived weights to indicators and criteria."),
                                        tags$li(tags$b("Aggregation"), " – using the Analytic Hierarchy Process (AHP) to produce ICSF scores for Mitigation, Adaptation, and Social dimensions.")
                                      )
                             )
                      ),
                      column(6,
                             tags$div(class = "info-card bottom-card",
                                      tags$h4("How to use & Main outputs"),
                                      tags$ol(
                                        tags$li(tags$b("Upload data"), " – Import one or more forest datasets in .csv format. Each file should represent one forest stand."),
                                        tags$li(tags$b("Calculate smartness values"), " – After uploading the dataset, click Calculate Smartness. The application computes: Plot-level outputs (displayed in a table) and Stand-level output (displayed in the 'Overall Smartness Values' section)"),
                                        tags$li(tags$b("Visualize results"), " – In the Graphic Comparison panel, stand-level values are shown using an interactive ternary plot, where each point represents the balance among Mitigation, Adaptation, and Social dimension."),
                                        tags$li(tags$b("Export results"), " – All outputs (plot-level and stand-level) can be downloaded as CSV files, and graphs can be saved as images.")
                                      )
                             )
                      )
             ),
             
             # Footer immediately after the cards, with reduced margins
             fluidRow(
               column(12,
                      tags$div(
                        style = "
          margin-top:6px;    /* piccolo spazio sopra il footer */
          padding:8px 12px;
          color:var(--white);
          font-size:12px;
        ",
                        tags$hr(style = "border:none; border-top:1px solid rgba(255,255,255,0.08); margin-bottom:8px;"),
                        tags$p(
                          tags$strong("References:"),
                          tags$br(),
                          " Alfieri, D., Tognetti, R., & Santopuoli, G. (2024). Exploring climate-smart forestry in Mediterranean forests through an innovative composite climate-smart index. ",
                          tags$em("Journal of Environmental Management, 368"), ", 122002.",
                          tags$br(),
                          "Bowditch, E., Santopuoli, G., Binder, F., Del Rio, M., La Porta, N., Kluvankova, T., Lesinski, J., Motta, R., Pach, M., Panzacchi, P., & others. (2020). What is Climate-Smart Forestry? A definition from a multinational collaborative process focused on mountain regions of Europe.",
                          tags$em("Ecosystem Services, 43,"), ", 101113",
                          tags$br(),
                          "Nabuurs, G.-J., Verkerk, P. J., Schelhaas, M.-J., González Olabarria, J. R., Trasobares, A., & Cienciala, E. (2018). Climate-Smart Forestry: mitigation impacts in three European regions (Vol. 6). EFI Helsinki, Finland",
                          tags$br(),
                          "Weatherall, A., Nabuurs, G., Velikova, V., Santopuoli, G., Neroj, B., Bowditch, E., ... & Tognetti, R. (2022). Defining climate-smart forestry. In Climate-Smart Forestry in Mountain Regions (Vol. 40, pp. 35-58). Springer.",
                          style = "margin-bottom:10px;"
                        ),
                        tags$p(
                          tags$strong("Developers and Maintainers:"), " Ricerca forestale - University of Molise - Diana Alfieri, Concetta Lisella, Serena Antonucci and Giovanni Santopuoli",
                          style = "margin-top:0;"
                        )
                      )
               )
             )
    ),
    # Smartness Calculation ####
    tabPanel("⚙️ Smartness Calculation",
             fluidRow(
               column(3,
                      tags$div(class = "left-panel",
                               numericInput("n_datasets", "Number of datasets", 2, min = 1),
                               helpText("Upload your datasets in CSV format (semicolon separator)."),
                               
                               br(),
                               downloadButton("download_all", "Download all", class = "btn"),
                               actionButton("reset_all", "Reset App", class = "btn mt-2"),
                               
                               br(), br(),
                               
                               tags$div(style = "font-size: 13px;",
                                        HTML(
                                          "<strong>CSV format:</strong><br> Please provide one CSV file per forest stand. All values should be numeric where present. Minimum: 2 plots per dataset. <br>
                                   <br>
                                   <strong>Required columns (include exact headers):</strong>
                                   <ul style='margin-top:6px; margin-bottom:6px;'>
                                   <li><em>id_plot</em> — Identification number of the sample plot</li>
                                   <li><em>year</em> — Year of survey <em>(if available)</em> </li>
                                   <li><em>CS</em> — Carbon Stock*</li>
                                   <li><em>GS</em> — Growing Stock*</li>
                                   <li><em>SP</em> — Tree Species Composition*</li>
                                   <li><em>DB</em> — Diameter Distribution*</li>
                                   <li><em>SLC</em> — Slenderness Coefficient*</li>
                                   <li><em>DW</em> — Deadwood</li>
                                   <li><em>RW</em> — Roundwood</li>
                                   <li><em>RN</em> — Regeneration</li>
                                   <li><em>FD</em> — Forest Damage</li>
                                   <li><em>IF</em> — Increment-Harvest</li>
                                   </ul>
                                   <p></p>
                                   <strong>*</strong> indicates the <strong>mandatory indicators</strong> for computing the smartness values.
                                   <br>
                                   If an indicator is not available or could not be measured, please still include the corresponding column header and fill its values with <strong><em>NA</em>.</strong><br>
                                   An <strong><em>example dataset</em></strong> is available in the <strong><em>GitHub repository</em></strong> <em>(https://github.com/ForestryLab/CSF-AT/tree/main/tables/example_datasets)</em> to guide users in preparing input data with the correct structure required by the application.<br><br>
                                   <strong>Annual evaluation:</strong> If your dataset contains a <em>year</em> column, the app can calculate CSF values for each year separately. <br><br>
                                   <strong>Table description:</strong> After clicking <code>Calculate Smartness</code>, a table is generated showing the uploaded indicators together with four additional columns reporting the climate-smartness scores computed for each stand:<br><br>
                                   <ul>
                                   <li><code>csf</code> – Composite Climate-Smart Index (overall smartness)</li>
                                   <li><code>csf_mit</code> – Mitigation score</li>
                                   <li><code>csf_adp</code> – Adaptation score</li>
                                   <li><code>csf_sd</code> – Social Dimension score</li>
                                   </ul>
                                   <br>
                                   <code>Select year to visualize:</code> After uploading and calculating smartness, this selector allows you to choose which year to display in the summary boxes. Each 'overall smartness value' will be the mean of all plots measured in that year. If only one plot is available for a year, its values will be used directly. You can also select 'All years' to see the overall average values across all available years.<br><br>
                                   ")
                               ),
                               tags$hr(),
                               
                               tags$div(style = "font-size: 13px;",
                                        tags$strong("Indicator weights used in the model"),
                                        br(),
                                        helpText("These weights are derived from expert judgments provided by scientific experts in the forestry sector at European level. Details on the composition and geographic distribution of experts involved in the weighting process are reported in the reference article on the development of the method (Alfieri et al., 2024)."),
                                        tags$div(class = "weights-scroll",
                                                 div(class = "weights-table-wrapper",
                                                     DTOutput("weights_table")
                                                 )
                                        )
                               )
                      )
               ),
               column(9,
                      uiOutput("dynamic_ui")
               )
             )
    ),
    # Graphic comparison ####
    tabPanel("📊 Graphic Comparison",
             br(), br(),
             tags$div(class = "ternary-row",
                      tags$div(class = "ternary-plot-panel",
                               div(class = "inner-white", style = "padding:6px; height:600px;",
                                   checkboxInput("ternary_show_years", "Show all years as separate points if available", value = TRUE),
                                   plotlyOutput("ternary_plot", height = "480px")
                               )
                      ),
                      tags$div(class = "ternary-desc-panel",
                               tags$div(class = "ternary-desc",
                                        HTML(
                                          "<strong>How to read the ternary plot</strong><br><br>
     This plot shows the relative composition of the three Climate-Smart Forestry dimensions (Mitigation, Adaptation, Social Dimension) for each forest stand.<br><br>

     <p><strong>What each point represents</strong>: each point represents a forest stand in a ternary space defined by the relative contributions of Mitigation, Adaptation, and Social Dimension.
     The coordinates are normalized to their sum and therefore reflect proportional shares rather than absolute CSF values. Point positions depend only on these relative contributions, 
     meaning that absolute CSF levels do not influence the location of points in the plot. When hovering over a point, the absolute CSF values and their standard deviations are displayed for reference.</p>
     
     <p>If your dataset contains a <em>year</em> column and you enable <code>Show all years as separate points</code>, the app will display one point per stand-year (representing the mean of plots measured in that year).</p>
    
     <strong>Interacting</strong>: Move the cursor over a point to see details, use Plotly controls to zoom or pan, click legend items to show or hide series, and use the Plotly toolbar to download the chart (Save as PNG).<br><br>

    <div style='margin-top:10px;'>
    This image shows an example of how to read the ternary plot and how the values should be interpreted along the three axes (Mitigation, Adaptation, Social Dimension).
    </div>"),
                                        tags$img(
                                          src = "ternary_guide.png",
                                          style = "width:100%; margin-top:10px; border-radius:8px;"
                                        )
                               )
                      )
             ),
             br(), br(),
             
             fluidRow(
               column(12,
                      div(class = "inner-white",
                          
                          h4("Smartness Distribution"),
                          
                          uiOutput("boxplot_desc"),
                          
                          uiOutput("box_dataset_ui"),
                          
                          plotlyOutput("boxplot_smartness", height = "450px")
                      )
               )
             ),
             
             br(), br(),
             fluidRow(
               column(12,
                      div(class = "inner-white",
                          
                          h4("Temporal Trend of Smartness"),
                          
                          uiOutput("trend_desc"),
                          
                          plotlyOutput("trend_plot", height = "400px")
                      )
               )
             )
    )
  )
)


# Server

server <- function(input, output, session) {
  
  data_list <- reactiveValues()
  filenames <- reactiveValues()
  meta <- reactiveValues()
  
  # Weights table ####
  
  output$weights_table <- renderDT({
    req(pesi)
    
    indicators <- setdiff(names(pesi), "pagg_c")
    
    indicators_df <- data.frame(
      Element = indicators,
      Mitigation = signif(sapply(indicators, function(ind)
        get_weight_or_zero(pesi, 1, ind)), 4),
      Adaptation = signif(sapply(indicators, function(ind)
        get_weight_or_zero(pesi, 2, ind)), 4),
      Social_Dimension = signif(sapply(indicators, function(ind)
        get_weight_or_zero(pesi, 3, ind)), 4)
    )
    
    criteria_df <- data.frame(
      Element = "CRITERIA",
      Mitigation = signif(get_weight_or_zero(pesi, 1, "pagg_c"), 4),
      Adaptation = signif(get_weight_or_zero(pesi, 2, "pagg_c"), 4),
      Social_Dimension = signif(get_weight_or_zero(pesi, 3, "pagg_c"), 4)
    )
    
    final_df <- rbind(criteria_df, indicators_df)
    datatable(
      final_df,
      options = list(
        paging = FALSE,
        searching = FALSE,
        dom = "t",
        autoWidth = TRUE,
        scrollX = TRUE,
        columnDefs = list(
          list(className = 'dt-center', targets = "_all"),
          list(width = 'auto', targets = "_all")
        )
      ),
      rownames = FALSE
    )
  })
  
  
  
  observeEvent(input$reset_all, { session$reload() })
  
  detect_col <- function(df, patterns) {
    nm_clean <- toupper(gsub("[^A-Z0-9]", "", names(df)))
    for (p in patterns) {
      p_clean <- toupper(gsub("[^A-Z0-9]", "", p))
      idx <- which(nm_clean == p_clean)
      if (length(idx) == 1) return(names(df)[idx])
    }
    for (p in patterns) {
      idx <- which(grepl(p, nm_clean, ignore.case = TRUE))
      if (length(idx) == 1) return(names(df)[idx])
    }
    return(NA_character_)
  }
  
  # dynamic UI: fileInput, action button, table, year select####
  output$dynamic_ui <- renderUI({
    req(input$n_datasets)
    lapply(seq_len(input$n_datasets), function(i) {
      ns <- paste0("ds", i)
      card(
        title = paste("Dataset", i),
        fileInput(ns, paste("Upload dataset", i), accept = ".csv"),
        actionButton(paste0("run_", ns), paste("Calculate Smartness")),
        DTOutput(paste0("table_", ns)),
        uiOutput(paste0("year_ui_", ns)),
        card("Overall Smartness Values",
             value_box("Smartness (mean ± sd)", textOutput(paste0("mean_csf_", ns)), theme = "success"),
             layout_columns(
               value_box("Mitigation (mean ± sd)", textOutput(paste0("mean_mit_", ns)), theme = "cyan"),
               value_box("Adaptation (mean ± sd)", textOutput(paste0("mean_adp_", ns)), theme = "blue"),
               value_box("Social Dimension (mean ± sd)", textOutput(paste0("mean_sd_", ns)), theme = "purple")
             )
        )
      )
    })
  })
  
  # handler for "Calculate Smartness" ####
  observe({
    req(input$n_datasets)
    for (i in seq_len(input$n_datasets)) {
      local({
        ii <- i
        ns <- paste0("ds", ii)
        run_id <- paste0("run_", ns)
        file_input_id <- ns
        
        observeEvent(input[[run_id]], {
          file_input <- input[[file_input_id]]
          req(file_input)
          df <- tryCatch({
            readr::read_delim(
              file_input$datapath,
              delim = NULL,  
              locale = locale(decimal_mark = ".", grouping_mark = ""),
              show_col_types = FALSE,
              progress = FALSE
            ) |> as.data.frame()
          }, error = function(e) {
            showModal(modalDialog(title = "File error", paste("Could not read the uploaded file:", e$message), easyClose = TRUE))
            return(NULL)
          })
          if (is.null(df)) return()
          
          names(df) <- gsub("﻿", "", names(df)); names(df) <- trimws(names(df))
          expected_cols <- c("CS","GS","DW","RN","SP","DB","RW","FD","IF","SLC")
          cm <- gsub("\\s+", "", toupper(names(df)))
          for (ec in expected_cols) {
            if (!(ec %in% names(df))) {
              idx <- which(cm == ec)
              if (length(idx) == 1) names(df)[idx] <- ec
            }
          }
          
          year_name <- detect_col(df, c("year", "anno"))
          if (!is.na(year_name) && !(year_name == "year")) names(df)[names(df) == year_name] <- "year"
          idplot_name <- detect_col(df, c("id_plot", "idplot", "plot_id", "plot", "id"))
          if (!is.na(idplot_name) && !(idplot_name == "id_plot")) names(df)[names(df) == idplot_name] <- "id_plot"
          
          missing_final <- expected_cols[!(expected_cols %in% names(df))]
          if (length(missing_final) > 0) {
            showModal(modalDialog(title = "Missing columns", paste0("The uploaded CSV must contain these columns: ", paste(expected_cols, collapse = ", "), ". Missing: ", paste(missing_final, collapse = ", ")), easyClose = TRUE))
            return()
          }
          
          df <- df %>% mutate(across(all_of(expected_cols), ~ suppressWarnings(as.numeric(gsub(",", ".", as.character(.))))))
          if ("year" %in% names(df)) df <- df %>% mutate(year = as.character(year))
          
          df <- df %>% filter(rowSums(is.na(select(., all_of(expected_cols)))) < length(expected_cols))
          if (nrow(df) < 2) {
            showModal(modalDialog(title = "Insufficient data", "Each dataset must contain at least 2 valid plot rows.", easyClose = TRUE))
            return()
          }
          
          res <- tryCatch({ calc_smartness(df, pesi) }, error = function(e) {
            showModal(modalDialog(title = "Calculation error", paste("Error during smartness calculation:", e$message), easyClose = TRUE))
            return(NULL)
          })
          if (is.null(res)) return()
          res[] <- lapply(res, function(x) suppressWarnings(as.numeric(as.character(x))))
          full <- cbind(df, res)
          
          # NOTS: removed creation of 'year_num' column so it won't appear in tables or downloads
          
          data_list[[ns]] <- full
          filenames[[ns]] <- tools::file_path_sans_ext(basename(file_input$name))
          meta[[paste0(ns, "_has_year")]] <- ("year" %in% names(full))
          meta[[paste0(ns, "_id_plot")]] <- ifelse("id_plot" %in% names(full), "id_plot", NA_character_)
          
          # populate the year select UI for this dataset (renderUI)
          output[[paste0("year_ui_", ns)]] <- renderUI({
            if (meta[[paste0(ns, "_has_year")]]) {
              years <- sort(unique(full$year))
              selectInput(paste0("year_sel_", ns),
                          label = "Select year to visualize:",
                          choices = c("All years", years),
                          selected = "All years")
            } else {
              tags$div(style = "margin-bottom:6px; color: #666;", "No 'year' column detected in this dataset.")
            }
          })
          
          # render datatable with scrollbars
          output[[paste0("table_", ns)]] <- renderDT({
            req(data_list[[ns]])
            datatable(
              data_list[[ns]],
              options = list(
                scrollX = TRUE,
                scrollY = "260px",
                paging = TRUE,
                pageLength = 10
              ),
              rownames = FALSE,
              class = 'stripe hover'
            )
          })
          
        }) # end observeEvent run_
      }) # end local
    } # end for
  }) # end observe
  
  # reactive outputs when select year: mean boxes
  observe({
    req(input$n_datasets)
    for (i in seq_len(input$n_datasets)) {
      local({
        ii <- i
        ns <- paste0("ds", ii)
        year_input_id <- paste0("year_sel_", ns)
        
        output[[paste0("mean_csf_", ns)]] <- renderText({
          df <- data_list[[ns]]; req(df)
          if (!is.null(meta[[paste0(ns, "_has_year")]]) && meta[[paste0(ns, "_has_year")]] &&
              !is.null(input[[year_input_id]]) && input[[year_input_id]] != "All years") {
            df_f <- df %>% filter(as.character(year) == as.character(input[[year_input_id]]))
          } else df_f <- df
          m <- mean(df_f$csf, na.rm = TRUE)
          s <- sd(df_f$csf, na.rm = TRUE)
          
          if (is.nan(m)) {
            "NA"
          } else {
            paste0(round(m, 2), " ± ", round(s, 2))
          }
        })
        
        output[[paste0("mean_mit_", ns)]] <- renderText({
          df <- data_list[[ns]]; req(df)
          if (!is.null(meta[[paste0(ns, "_has_year")]]) && meta[[paste0(ns, "_has_year")]] &&
              !is.null(input[[year_input_id]]) && input[[year_input_id]] != "All years") {
            df_f <- df %>% filter(as.character(year) == as.character(input[[year_input_id]]))
          } else df_f <- df
          m <- mean(df_f$csf_mit, na.rm = TRUE)
          s <- sd(df_f$csf_mit, na.rm = TRUE)
          
          if (is.nan(m)) "NA" else paste0(round(m, 2), " ± ", round(s, 2))
        })
        
        output[[paste0("mean_adp_", ns)]] <- renderText({
          df <- data_list[[ns]]; req(df)
          if (!is.null(meta[[paste0(ns, "_has_year")]]) && meta[[paste0(ns, "_has_year")]] &&
              !is.null(input[[year_input_id]]) && input[[year_input_id]] != "All years") {
            df_f <- df %>% filter(as.character(year) == as.character(input[[year_input_id]]))
          } else df_f <- df
          m <- mean(df_f$csf_adp, na.rm = TRUE)
          s <- sd(df_f$csf_adp, na.rm = TRUE)
          
          if (is.nan(m)) "NA" else paste0(round(m, 2), " ± ", round(s, 2))
        })
        
        output[[paste0("mean_sd_", ns)]] <- renderText({
          df <- data_list[[ns]]; req(df)
          if (!is.null(meta[[paste0(ns, "_has_year")]]) && meta[[paste0(ns, "_has_year")]] &&
              !is.null(input[[year_input_id]]) && input[[year_input_id]] != "All years") {
            df_f <- df %>% filter(as.character(year) == as.character(input[[year_input_id]]))
          } else df_f <- df
          m <- mean(df_f$csf_sd, na.rm = TRUE)
          s <- sd(df_f$csf_sd, na.rm = TRUE)
          
          if (is.nan(m)) "NA" else paste0(round(m, 2), " ± ", round(s, 2))
        })
        
      })
    }
  })
  
  # ternary plot ####
  output$ternary_plot <- renderPlotly({
    req(input$n_datasets)
    show_per_year <- isTRUE(input$ternary_show_years)
    
    symbols_vec <- c('circle','square','diamond','cross','triangle-up','triangle-down',
                     'triangle-down-open','star','x','hourglass')
    
    rows_list <- lapply(seq_len(input$n_datasets), function(i) {
      ns <- paste0("ds", i)
      df <- data_list[[ns]]
      if (is.null(df)) return(NULL)
      
      if (show_per_year && !is.null(meta[[paste0(ns, "_has_year")]]) && meta[[paste0(ns, "_has_year")]]) {
        years <- sort(unique(df$year))
        per_year <- lapply(years, function(yy) {
          df_y <- df %>% filter(as.character(year) == as.character(yy))
          nplots <- nrow(df_y)
          mit <- mean(df_y$csf_mit, na.rm = TRUE)
          adp <- mean(df_y$csf_adp, na.rm = TRUE)
          sd  <- mean(df_y$csf_sd, na.rm = TRUE)
          tot <- mit + adp + sd
          if (is.na(tot) || tot == 0) tot <- 1
          data.frame(
            Dataset = if (!is.null(filenames[[ns]])) filenames[[ns]] else paste0("Dataset ", i),
            Year = as.character(yy),
            Nplots = nplots,
            Mitigation = (mit / tot) * 100,
            Adaptation = (adp / tot) * 100,
            Socio_economic = (sd / tot) * 100,
            csf = round(mean(df_y$csf, na.rm = TRUE), 3),
            csf_sd = round(sd(df_y$csf, na.rm = TRUE), 3),
            
            mit_raw = round(mit, 2),
            mit_sd = round(sd(df_y$csf_mit, na.rm = TRUE), 2),
            
            adp_raw = round(adp, 2),
            adp_sd = round(sd(df_y$csf_adp, na.rm = TRUE), 2),
            
            sd_raw = round(sd, 2),
            sd_sd = round(sd(df_y$csf_sd, na.rm = TRUE), 2),
            Symbol = symbols_vec[(i-1) %% length(symbols_vec) + 1],
            stringsAsFactors = FALSE
          )
        }) %>% bind_rows()
        return(per_year)
      } else {
        year_input_id <- paste0("year_sel_", ns)
        df_f <- df
        if (!is.null(meta[[paste0(ns, "_has_year")]]) && meta[[paste0(ns, "_has_year")]] &&
            !is.null(input[[year_input_id]]) && input[[year_input_id]] != "All years") {
          df_f <- df %>% filter(as.character(year) == as.character(input[[year_input_id]]))
        }
        mit <- mean(df_f$csf_mit, na.rm = TRUE)
        adp <- mean(df_f$csf_adp, na.rm = TRUE)
        sd  <- mean(df_f$csf_sd, na.rm = TRUE)
        tot <- mit + adp + sd
        if (is.na(tot) || tot == 0) tot <- 1
        data.frame(
          Dataset = if (!is.null(filenames[[ns]])) filenames[[ns]] else paste0("Dataset ", i),
          Year = if (!is.null(meta[[paste0(ns, "_has_year")]]) && meta[[paste0(ns, "_has_year")]]) as.character(input[[year_input_id]]) else "All",
          Nplots = nrow(df_f),
          Mitigation = (mit / tot) * 100,
          Adaptation = (adp / tot) * 100,
          Socio_economic = (sd / tot) * 100,
          
          csf = round(mean(df_f$csf, na.rm = TRUE), 3),
          csf_sd = round(sd(df_f$csf, na.rm = TRUE), 3),
          
          mit_raw = round(mit, 2),
          mit_sd = round(sd(df_f$csf_mit, na.rm = TRUE), 2),
          
          adp_raw = round(adp, 2),
          adp_sd = round(sd(df_f$csf_adp, na.rm = TRUE), 2),
          
          sd_raw = round(sd, 2),
          sd_sd = round(sd(df_f$csf_sd, na.rm = TRUE), 2),
          
          Symbol = symbols_vec[(i-1) %% length(symbols_vec) + 1],
          stringsAsFactors = FALSE
        )
      }
    })
    
    ternary_df <- bind_rows(rows_list)
    req(nrow(ternary_df) > 0)
    
    hover_txt <- ~paste(
      "Dataset:", Dataset,
      "<br>Year:", Year,
      "<br>N plots:", Nplots,
      
      "<br><br>Mitigation:",
      sprintf("%.2f ± %.2f", mit_raw, mit_sd),
      
      "<br>Adaptation:",
      sprintf("%.2f ± %.2f", adp_raw, adp_sd),
      
      "<br>Social Dimension:",
      sprintf("%.2f ± %.2f", sd_raw, sd_sd),
      
      "<br>CSF:",
      sprintf("%.2f ± %.2f", csf, csf_sd)
    )
    
    if (show_per_year) {
      p <- plot_ly(
        data = ternary_df,
        type = 'scatterternary',
        mode = 'markers+text',
        a = ~Mitigation, b = ~Adaptation, c = ~Socio_economic,
        color = ~Year,
        symbol = ~Dataset,
        text = ~Year,
        textposition = 'top center',
        hoverinfo = 'text',
        hovertext = hover_txt,
        marker = list(size = 9)
      )
    } else {
      p <- plot_ly(
        data = ternary_df,
        type = 'scatterternary',
        mode = 'markers',
        a = ~Mitigation, b = ~Adaptation, c = ~Socio_economic,
        color = ~Dataset,
        symbol = ~Dataset,
        text = hover_txt,
        hoverinfo = 'text',
        marker = list(size = 12)
      )
    }
    
    p %>% layout(
      ternary = list(
        sum = 100,
        aaxis = list(title = "Mitigation", min = 0, max = 100, ticksuffix = "%"),
        baxis = list(title = "Adaptation", min = 0, max = 100, ticksuffix = "%"),
        caxis = list(title = "Social Dimension", min = 0, max = 100, ticksuffix = "%")
      ),
      margin = list(t = 80, b = 40, l = 40, r = 40),
      paper_bgcolor = 'rgba(0,0,0,0)',
      plot_bgcolor = 'rgba(0,0,0,0)'
    ) %>% config(displayModeBar = TRUE, displaylogo = FALSE,
                 toImageButtonOptions = list(format = "png", filename = "ternary_plot", height = 900, width = 1200, scale = 2))
  })
  
  # boxplot ####
  
  output$box_dataset_ui <- renderUI({
    req(input$n_datasets)
    
    has_year_any <- any(sapply(seq_len(input$n_datasets), function(i) {
      ns <- paste0("ds", i)
      !is.null(meta[[paste0(ns, "_has_year")]]) &&
        meta[[paste0(ns, "_has_year")]]
    }))
    
    # If no dataset has year → no selector
    if (!has_year_any) {
      return(NULL)
    }
    
    # Otherwise show selectInput
    choices <- lapply(seq_len(input$n_datasets), function(i) {
      ns <- paste0("ds", i)
      if (!is.null(filenames[[ns]])) filenames[[ns]] else paste0("Dataset ", i)
    })
    
    names(choices) <- choices
    
    selectInput("box_dataset", "Select dataset:", choices = choices)
  })
  
  output$boxplot_smartness <- renderPlotly({
    req(input$n_datasets)
    
    # Check if at least one dataset has years
    has_year_any <- any(sapply(seq_len(input$n_datasets), function(i) {
      ns <- paste0("ds", i)
      !is.null(meta[[paste0(ns, "_has_year")]]) &&
        meta[[paste0(ns, "_has_year")]]
    }))
    
    # 1° CASE: No YEAR → comparison between datasets
    
    if (!has_year_any) {
      
      df_all <- lapply(seq_len(input$n_datasets), function(i) {
        ns <- paste0("ds", i)
        df <- data_list[[ns]]
        if (is.null(df)) return(NULL)
        
        data.frame(
          Dataset = if (!is.null(filenames[[ns]])) filenames[[ns]] else paste0("Dataset ", i),
          csf = df$csf
        )
      }) %>% bind_rows()
      
      req(nrow(df_all) > 0)
      
      return(
        plot_ly(
          data = df_all,
          x = ~Dataset,
          y = ~csf,
          type = "box",
          boxpoints = "all",
          jitter = 0.3,
          pointpos = 0
        ) %>%
          layout(
            xaxis = list(title = "Dataset"),
            yaxis = list(title = "Smartness (CSF)"),
            paper_bgcolor = 'rgba(0,0,0,0)',
            plot_bgcolor = 'rgba(0,0,0,0)'
          )
      )
    }
    
    # 2° CASE: THERE IS YEAR → dropdown
    req(input$box_dataset)
    
    selected_name <- input$box_dataset
    
    ns_selected <- NULL
    for (i in seq_len(input$n_datasets)) {
      ns <- paste0("ds", i)
      nm <- if (!is.null(filenames[[ns]])) filenames[[ns]] else paste0("Dataset ", i)
      if (nm == selected_name) {
        ns_selected <- ns
        break
      }
    }
    
    req(ns_selected)
    df <- data_list[[ns_selected]]
    req(df)
    
    if ("year" %in% names(df)) {
      
      plot_ly(
        data = df,
        x = ~as.factor(year),
        y = ~csf,
        type = "box",
        boxpoints = "all",
        jitter = 0.3,
        pointpos = 0
      ) %>%
        layout(
          xaxis = list(title = "Year"),
          yaxis = list(title = "Smartness (CSF)"),
          paper_bgcolor = 'rgba(0,0,0,0)',
          plot_bgcolor = 'rgba(0,0,0,0)'
        )
      
    } else {
      
      plot_ly(
        data = df,
        y = ~csf,
        type = "box",
        boxpoints = "all",
        jitter = 0.3,
        pointpos = 0
      ) %>%
        layout(
          xaxis = list(title = "Plots"),
          yaxis = list(title = "Smartness (CSF)"),
          paper_bgcolor = 'rgba(0,0,0,0)',
          plot_bgcolor = 'rgba(0,0,0,0)'
        )
    }
  })
  
  output$boxplot_desc <- renderUI({
    req(input$n_datasets)
    
    has_year_any <- any(sapply(seq_len(input$n_datasets), function(i) {
      ns <- paste0("ds", i)
      !is.null(meta[[paste0(ns, "_has_year")]]) &&
        meta[[paste0(ns, "_has_year")]]
    }))
    
    HTML(paste0("
  <div style='margin-bottom:10px; font-size:13px;'>
  
  <strong>How to read this boxplot</strong><br><br>
  
  This plot shows the <strong>distribution of Climate-Smartness (CSF) values</strong> at plot level for each forest stand.<br><br>
  
  Each box summarizes the variability of smartness values:
  <ul>
    <li><strong>Median</strong> – central value</li>
    <li><strong>Box</strong> – interquartile range (variability)</li>
    <li><strong>Points</strong> – individual plot values</li>
  </ul>
  ",
                
                if (has_year_any) {
                  "
    <br>
    If datasets include a <em>year</em> column, the boxplot shows the distribution of values across <strong>years</strong> for the selected forest stand.<br>
    The x-axis represents years, and you can use the selector above to choose which dataset (forest stand) to visualize.
    "
                } else {
                  "
    <br>
    If no datasets include a <em>year</em> column, the boxplot compares the distributions across <strong>different forest stands</strong>.<br>
    "
                },
                
                "</div>
  "))
  })
  
  
  # trend plot ####
  
  output$trend_plot <- renderPlotly({
    req(input$n_datasets)
    
    has_year_any <- any(sapply(seq_len(input$n_datasets), function(i) {
      ns <- paste0("ds", i)
      !is.null(meta[[paste0(ns, "_has_year")]]) &&
        meta[[paste0(ns, "_has_year")]]
    }))
    
    # if there are no years → no plot
    if (!has_year_any) return(NULL)
    
    rows_list <- lapply(seq_len(input$n_datasets), function(i) {
      ns <- paste0("ds", i)
      df <- data_list[[ns]]
      if (is.null(df)) return(NULL)
      if (!"year" %in% names(df)) return(NULL)
      
      dataset_name <- if (!is.null(filenames[[ns]])) filenames[[ns]] else paste0("Dataset ", i)
      
      df %>%
        group_by(year) %>%
        summarise(
          mean_csf = mean(csf, na.rm = TRUE),
          sd_csf = sd(csf, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(Dataset = dataset_name)
    })
    
    trend_df <- bind_rows(rows_list)
    req(nrow(trend_df) > 0)
    
    plot_ly(
      data = trend_df,
      x = ~year,
      y = ~mean_csf,
      color = ~Dataset,
      text = ~Dataset,
      type = "scatter",
      mode = "lines+markers",
      hovertemplate = paste(
        "Dataset: %{text}<br>",
        "Year: %{x}<br>",
        "Mean CSF: %{y:.2f}<extra></extra>"
      )
    ) %>%
      layout(
        xaxis = list(title = "Year"),
        yaxis = list(title = "Mean Smartness (CSF)"),
        legend = list(orientation = "v", x = 1.02, y = 1, xanchor = "left"),
        paper_bgcolor = 'rgba(0,0,0,0)',
        plot_bgcolor = 'rgba(0,0,0,0)'
      )
  })
  
  output$trend_desc <- renderUI({
    
    has_year_any <- any(sapply(seq_len(input$n_datasets), function(i) {
      ns <- paste0("ds", i)
      !is.null(meta[[paste0(ns, "_has_year")]]) &&
        meta[[paste0(ns, "_has_year")]]
    }))
    
    if (!has_year_any) {
      return(HTML("
      <div style='margin-bottom:10px; font-size:13px;'>
      <em>No temporal trend available: none of the uploaded datasets contains a <strong>year</strong> column.</em>
      </div>
    "))
    }
    
    HTML("
  <div style='margin-bottom:10px; font-size:13px;'>
  
  <strong>How to read this plot</strong><br><br>
  
  This graph shows the <strong>temporal evolution</strong> of CSF values.<br><br>
  
  Each line represents a <strong>forest stand (dataset)</strong>, and each point corresponds to the 
  <strong>mean CSF value</strong> of all plots measured in a given year.<br><br>
  
  This visualization helps to:
  <ul>
    <li>Identify increasing or decreasing trends in smartness</li>
    <li>Compare trajectories among different stands</li>
    <li>Understand how smartness evolves over time</li>
  </ul>
  
  </div>
  ")
  })
  
  # download all ####
  output$download_all <- downloadHandler(
    
    filename = function() {
      paste0("smartness_results_", Sys.Date(), ".zip")
    },
    
    content = function(zipfile) {
      
      tmp_dir <- tempfile("csf_export_")
      dir.create(tmp_dir)
      
      file_list <- c()
      
      for (i in seq_len(input$n_datasets)) {
        
        ns <- paste0("ds", i)
        
        df <- data_list[[ns]]
        
        if (is.null(df))
          next
        
        name_prefix <- if (!is.null(filenames[[ns]])) {
          filenames[[ns]]
        } else {
          paste0("dataset_", i)
        }
        
        # --------------------------
        # Plot-level results
        # --------------------------
        
        export_df <- df %>%
          select(-matches("_norm$"))
        
        plot_file <- file.path(
          tmp_dir,
          paste0(name_prefix, "_plot_level.csv")
        )
        
        write.csv(
          export_df,
          plot_file,
          row.names = FALSE
        )
        
        file_list <- c(file_list, plot_file)
        
        # --------------------------
        # Stand-level results
        # --------------------------
        
        if (isTRUE(meta[[paste0(ns, "_has_year")]])) {
          
          means_df <- df %>%
            group_by(year) %>%
            summarise(
              Nplots = n(),
              
              Smartness = mean(csf, na.rm = TRUE),
              Smartness_sd = sd(csf, na.rm = TRUE),
              
              Mitigation = mean(csf_mit, na.rm = TRUE),
              Mitigation_sd = sd(csf_mit, na.rm = TRUE),
              
              Adaptation = mean(csf_adp, na.rm = TRUE),
              Adaptation_sd = sd(csf_adp, na.rm = TRUE),
              
              Socio_economic = mean(csf_sd, na.rm = TRUE),
              Socio_economic_sd = sd(csf_sd, na.rm = TRUE),
              
              .groups = "drop"
            )
          
        } else {
          
          means_df <- data.frame(
            year = NA_character_,
            Nplots = nrow(df),
            
            Smartness = mean(df$csf, na.rm = TRUE),
            Smartness_sd = sd(df$csf, na.rm = TRUE),
            
            Mitigation = mean(df$csf_mit, na.rm = TRUE),
            Mitigation_sd = sd(df$csf_mit, na.rm = TRUE),
            
            Adaptation = mean(df$csf_adp, na.rm = TRUE),
            Adaptation_sd = sd(df$csf_adp, na.rm = TRUE),
            
            Socio_economic = mean(df$csf_sd, na.rm = TRUE),
            Socio_economic_sd = sd(df$csf_sd, na.rm = TRUE),
            
            stringsAsFactors = FALSE
          )
        }
        
        stand_file <- file.path(
          tmp_dir,
          paste0(name_prefix, "_stand_level.csv")
        )
        
        write.csv(
          means_df,
          stand_file,
          row.names = FALSE
        )
        
        file_list <- c(file_list, stand_file)
      }
      
      # No file available
      if (length(file_list) == 0) {
        stop("No results available for download.")
      }
      
      # ZIP Creation
      zip::zipr(
        zipfile = zipfile,
        files = file_list
      )
    }
  )
}

# Launch the app ####
shinyApp(ui, server)
