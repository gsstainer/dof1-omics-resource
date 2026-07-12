
################################################################################
#
# DO-F1 Microbiome & Atherosclerosis Interactive Data Resource
#
# Author: Myungsuk Kim
# Date: June 20, 2025
#
# Description:
# This Shiny application, presented as an interactive dashboard, serves as a
# supplement for the manuscript "Dissecting Gene-Microbiome Interactions...".
# It allows users to dynamically explore microbiota composition, diversity,
# heritability, and expanded QTL mapping results from the DO-F1 mouse study.
#
################################################################################


# --- 1. LOAD LIBRARIES ---
# Ensure these packages are installed:
# install.packages(c("shiny", "shinydashboard", "plotly", "ggplot2", "dplyr", "tidyr", "ggpubr", "ape", "DT", "BiocManager"))
# BiocManager::install("qtl2")
options(repos = BiocManager::repositories()) #to resolve the renv_snapshot_validate_report error

library(shiny)
library(shinydashboard) # For the dashboard layout
library(plotly)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(qtl2)
library(phyloseq)
library(vegan)
library(qiime2R)
library(ape)
library(DT)


# --- 2. LOAD AND PREPARE DATA ---
# IMPORTANT: Place all required data files in a subdirectory named 'data/'

# --- Load Metadata and Pre-processed Tables ---
# Using suppressWarnings to hide minor parsing messages on load
metadata <- read.csv("data/metadata.csv")
lmm_results <- read.csv("data/LMM_Results_Microbiome_vs_Traits.csv")
h2_f1_data <- read.csv("data/h2_DOProF1.csv")
h2_dof1_data <- read.csv("data/h2_DOF1.csv")
load("data/pso_rare.RData")

# --- Load QTL-related R Data Objects ---
cross2 <- readRDS("data/DOF1POST_retype_cross2_mb.rds")
apr <- readRDS("data/DOF1_retype_probs.rds")
kloco <- readRDS("data/DOF1_retype_kLOCO.rds")
covar <- readRDS("data/DOF1_covar.rds")

# Load scan results
scan_female <- readRDS("data/DOF1POST_fullOut_female_mb_origin_f9.rds")
scan_male <- readRDS("data/DOF1POST_fullOut_male_mb_origin_f9.rds")
scan_additive <- readRDS("data/DOF1POST_fullOut_mb_origin_f9.rds")

# --- NOTE: The creation of 'qtl_choices' has been moved into the server function ---
# This makes the app more robust during deployment.

# Define mouse subsets
all_mice <- rownames(cross2$covar)
female_mice <- rownames(cross2$covar)[cross2$covar$sex == "F"]
male_mice <- rownames(cross2$covar)[cross2$covar$sex == "M"]

# Fix tree structure to prevent warnings
tree <- phy_tree(pso_rare)
tree <- ape::multi2di(tree)
phy_tree(pso_rare) <- tree


# --- 3. UI DEFINITION ---
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "DO-F1 Microbiome"),
  
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("About", tabName = "about", icon = icon("info-circle")),
      menuItem("Cohort Analysis", tabName = "cohort", icon = icon("users")),
      menuItem("Heritability", tabName = "heritability", icon = icon("dna")),
      menuItem("Microbe-Trait LMM", tabName = "lmm", icon = icon("chart-line")),
      menuItem("Microbiome QTL", tabName = "qtl", icon = icon("search-location"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # About Tab
      tabItem(tabName = "about",
              fluidRow(
                valueBoxOutput("mice_count", width=3),
                valueBoxOutput("female_count", width=3),
                valueBoxOutput("male_count", width=3),
                valueBoxOutput("genus_count", width=3)
              ),
              fluidRow(
                box(
                  title = "About this Resource", status = "primary", solidHeader = TRUE, width = 12,
                  h3("Welcome to the DO-F1 Microbiome & Atherosclerosis Data Resource"),
                  p("This application serves as an interactive supplement to the manuscript: ", 
                    strong("Dissecting Gene-Microbiome Interactions in Atherosclerosis using a Systems Genetics Approach in a Large-Scale Mouse Population.")),
                  p("Here, you can dynamically explore the key findings of our study, including microbiota composition, diversity, heritability, linear mixed-model associations, and quantitative trait locus (QTL) mapping results from the DO-F1 mouse study."),
                  h4("Data & Code Availability:"),
                  p("Raw sequencing data are available at SRA (PRJNA686143). Processed data and scripts are available at our GitHub repository [LINK].")
                )
              ),
              
              fluidRow(
                box(
                  title = "Contact Information", status = "primary", solidHeader = TRUE, width = 12,
                  p("For questions about the study or data, please contact the Principal Investigator:"),
                  p(strong("Brian Bennett, Ph.D."), br(), "USDA Western Human Nutrition Research Center", br(), "Email: brian.bennett@usda.gov", tags$a(href="mailto:brian.bennett@usda.gov")),
                  hr(),
                  p("For questions regarding this Shiny application, please contact the developer:"),
                  p(strong("Myungsuk Kim"), br(), "Email: g-sstainer@kist.re.kr", tags$a(href="mailto:g-sstainer@kist.re.kr"))
                )
              )
      ), 
      # Cohort Differences Tab
      tabItem(tabName = "cohort",
              fluidRow(
                box(title = "Controls", status = "warning", solidHeader = TRUE, width = 12,
                    column(6, selectInput("alpha_metric_cohort", "Choose Alpha Diversity Metric:",
                                          choices = c("Shannon" = "shannon", "Faith's PD" = "faith_pd", "Observed ASVs" = "observed_otus"))),
                    column(6, selectInput("beta_metric_cohort", "Choose Beta Diversity Metric:",
                                          choices = c("Weighted UniFrac" = "wunifrac", "Unweighted UniFrac" = "uunifrac", "Bray-Curtis" = "bray")))
                )
              ),
              fluidRow(
                box(title = "Alpha Diversity", status = "primary", width = 6, plotlyOutput("alphaPlot_cohort")),
                box(title = "Beta Diversity", status = "primary", width = 6, plotOutput("betaPlot_cohort"))
              )
      ),
      
      # Heritability Tab
      tabItem(tabName = "heritability",
              fluidRow(
                box(title = "Controls", status = "warning", solidHeader = TRUE, width = 12,
                    radioButtons("h2_cohort", "Select Cohort:", choices = c("DO-F1 Progeny", "F1 Founder Crosses"), inline = TRUE)
                )
              ),
              fluidRow(
                box(title = "Heritability Estimates (h²)", status = "primary", width = 12, plotlyOutput("heritabilityPlot", height = "800px"))
              )
      ),
      
      # Microbe-Trait LMM Tab
      tabItem(tabName = "lmm",
              fluidRow(
                box(title = "Controls", status = "warning", solidHeader = TRUE, width = 12,
                    column(6, selectInput("lmm_trait", "Select Cardiometabolic Trait:", choices = c("Aortic lesion area", "Plasma total cholesterol"))),
                    column(6, selectInput("lmm_genus", "Select Microbial Genus:", choices = NULL))
                )
              ),
              fluidRow(
                box(title = "Association Plot", status = "primary", width = 7, plotOutput("lmmPlot")),
                box(title = "LMM Statistics", status = "primary", width = 5, DT::dataTableOutput("lmmTable"))
              )
      ),
      
      # Microbiome QTL Tab
      tabItem(tabName = "qtl",
              fluidRow(
                box(title = "QTL Controls", status = "warning", solidHeader = TRUE, width = 12,
                    # Choices are set to NULL, to be populated by the server.
                    column(6, selectInput("qtl_trait", "Select Microbial Genus:", choices = NULL)),
                    column(6, radioButtons("qtl_model", "Select QTL Model:", choices = c("Female", "Male", "Additive"), inline = TRUE))
                )
              ),
              fluidRow(
                box(title = "Genome-wide LOD Scan", status = "primary", solidHeader = TRUE, width = 12, 
                    plotOutput("qtlScanPlot", height="350px")),
                box(title = "Allele Effects at Peak LOD", status = "primary", solidHeader = TRUE, width = 12, 
                    plotOutput("qtlAllelePlot", height="350px"))
              )
      )
    )
  )
)


# --- 4. SERVER LOGIC ---
server <- function(input, output, session) { # Added 'session' for updateSelectInput
  
  # --- MOVED FROM GLOBAL SCOPE ---
  # This section is now inside the server function. This ensures that these
  # data-dependent objects are created only when a user session starts,
  # making the app more robust, especially when deployed.
  scanned_qtl_traits <- colnames(scan_female)
  known_metadata_fields <- c(
    "SampleID", "ID", # Common ID names
    "Sex", "ngen", "Cage_ID", "Dam_ID", # Covariates used in plots
    "shannon", "faith_pd", "observed_otus", # Alpha diversity metrics
    "lesionArea", "folChol", "liverTC", "folGLC", "folTG" # Other phenotypes
  )
  genera_in_metadata <- colnames(metadata)[!colnames(metadata) %in% known_metadata_fields]
  qtl_choices <- intersect(scanned_qtl_traits, genera_in_metadata)
  # --- END OF MOVED BLOCK ---
  
  # --- About Page Logic ---
  output$mice_count <- renderValueBox({
    valueBox(nrow(metadata), "Total Mice", icon = icon("mouse"), color = "purple")
  })
  output$female_count <- renderValueBox({
    valueBox(length(female_mice), "Female Mice", icon = icon("venus"), color = "orange")
  })
  output$male_count <- renderValueBox({
    valueBox(length(male_mice), "Male Mice", icon = icon("mars"), color = "blue")
  })
  output$genus_count <- renderValueBox({
    # Use the 'qtl_choices' object which is now defined within the server function
    valueBox(length(qtl_choices), "Microbial Genera with QTLs", icon = icon("bacteria"), color = "green")
  })
  
  # --- Cohort Differences Logic ---
  output$alphaPlot_cohort <- renderPlotly({
    req(input$alpha_metric_cohort)
    color_map <- c("Wave 1" = "#FFDC00", "Wave 2" = "#888888")
    y_max <- max(metadata[[input$alpha_metric_cohort]], na.rm = TRUE)
    
    p <- ggplot(metadata, aes_string(x = "ngen", y = input$alpha_metric_cohort, fill = "ngen")) +
      geom_boxplot(notch = TRUE, outlier.alpha = 0.5) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      scale_fill_manual(values = color_map, name = "Cohort") +
      stat_compare_means(method = "wilcox.test", label.y = y_max * 1.05) +
      labs(title = paste("Alpha Diversity:", input$alpha_metric_cohort), x = "Cohort", y = "Value") +
      theme_bw() + theme(legend.position = "none")
    ggplotly(p)
  })
  
  beta_results <- reactive({
    req(input$beta_metric_cohort, pso_rare)
    dist_matrix <- phyloseq::distance(pso_rare, method = input$beta_metric_cohort)
    sample_df <- data.frame(sample_data(pso_rare))
    perm_res <- adonis2(dist_matrix ~ ngen, data = sample_df, permutations = 999)
    pcoa <- ordinate(pso_rare, "PCoA", distance = dist_matrix)
    list(pcoa = pcoa, r2 = perm_res$R2[1], pval = perm_res$`Pr(>F)`[1])
  })
  
  output$betaPlot_cohort <- renderPlot({
    res <- beta_results()
    perm_label <- paste0("PERMANOVA\nR² = ", round(res$r2, 3), ", p = ", format.pval(res$pval, digits = 3))
    plot_ordination(pso_rare, res$pcoa, color = "ngen") +
      geom_point(size = 3, alpha = 0.7) +
      stat_ellipse(type = "norm", level = 0.95) +
      scale_color_manual(values = c("Wave 1" = "#FFDC00", "Wave 2" = "#888888"), name = "Cohort") +
      annotate("text", x = Inf, y = -Inf, label = perm_label, hjust = 1.1, vjust = -0.5, size = 4) +
      theme_bw() + labs(title = paste("Beta Diversity:", input$beta_metric_cohort))
  })
  
  # --- Heritability Logic ---
  output$heritabilityPlot <- renderPlotly({
    df <- if (input$h2_cohort == "F1 Founder Crosses") h2_f1_data else h2_dof1_data
    colnames(df)[1:2] <- c("Trait", "Heritability")
    df <- df %>% filter(!is.na(Heritability)) %>% arrange(Heritability)
    
    p <- ggplot(df, aes(x = reorder(Trait, Heritability), y = Heritability, fill = Heritability)) +
      geom_col() + coord_flip() +
      scale_fill_gradient(low = "lightblue", high = "navy") +
      labs(x = NULL, y = "Heritability (%)") +
      theme_bw() + theme(legend.position = "none", axis.text.y = element_text(face="italic", size=8))
    ggplotly(p, tooltip = c("x", "y"))
  })
  
  # --- Microbe-Trait LMM Logic ---
  
  selected_trait_prog <- reactive({
    if (input$lmm_trait == "Aortic lesion area") "lesionArea" else "folChol"
  })
  
  observe({
    req(selected_trait_prog())
    prog_trait <- selected_trait_prog()
    
    available_genera <- lmm_results %>%
      filter(Trait == prog_trait) %>%
      pull(Microbial_Feature) %>%
      unique() %>%
      sort()
    
    updateSelectInput(session, "lmm_genus",
                      choices = available_genera,
                      selected = available_genera[1])
  })
  
  output$lmmPlot <- renderPlot({
    req(input$lmm_genus)
    prog_trait <- selected_trait_prog()
    
    df_lmm <- metadata %>%
      mutate(across(all_of(c(input$lmm_genus, prog_trait)), ~ log1p(as.numeric(.x))))
    
    ggplot(df_lmm, aes_string(x = paste0("`", input$lmm_genus, "`"), y = prog_trait, color = "Sex")) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE, aes(fill = Sex)) +
      scale_color_manual(values = c("Female" = "red", "Male" = "blue")) +
      scale_fill_manual(values = c("Female" = "pink", "Male" = "lightblue")) +
      stat_cor(method = "spearman", aes(label = ..r.label..), label.y.npc = 0.95) +
      labs(x = paste("log(", input$lmm_genus, "Abundance + 1)"), y = paste("log(", input$lmm_trait, "+ 1)")) +
      theme_bw() + theme(legend.position = "bottom")
  })
  
  output$lmmTable <- DT::renderDataTable({
    req(input$lmm_genus)
    prog_trait <- selected_trait_prog()
    
    lmm_results %>%
      filter(Trait == prog_trait, Microbial_Feature == input$lmm_genus) %>%
      select(Estimate, CI_Low, CI_High, P_value, R2_Marginal, Heritability_h2) %>%
      mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
      DT::datatable(options = list(dom = 't', pageLength = 5, searching = FALSE), rownames = FALSE)
  })
  
  
  # --- Microbiome QTL Logic ---
  
  # Update the QTL trait dropdown once the server starts.
  # This happens only once.
  observeEvent(TRUE, {
    updateSelectInput(session, "qtl_trait", choices = qtl_choices, selected = qtl_choices[1])
  })
  
  reactive_scan <- reactive({
    switch(input$qtl_model,
           "Female" = scan_female,
           "Male" = scan_male,
           "Additive" = scan_additive)
  })
  
  output$qtlScanPlot <- renderPlot({
    req(input$qtl_trait, reactive_scan())
    scan_obj <- reactive_scan()
    plot_scan1(scan_obj, cross2$pmap, lodcolumn = input$qtl_trait, 
               main = paste("QTL Scan for", input$qtl_trait, "| Model:", input$qtl_model))
    abline(h = 7.25, col = "red", lty = 2) # Example significance threshold
  })
  
  output$qtlAllelePlot <- renderPlot({
    req(input$qtl_trait, input$qtl_model, reactive_scan())
    
    scan_obj <- reactive_scan()
    peak_chr <- find_peaks(scan_obj, cross2$pmap, threshold=4, prob=0.95) %>%
      filter(lodcolumn == input$qtl_trait) %>%
      arrange(desc(lod)) %>%
      pull(chr) %>%
      first()
    
    if (is.na(peak_chr) || length(peak_chr) == 0) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      text(1, 1, "No significant QTL peak found for this trait/model.", cex = 1.2)
      return()
    }
    
    mice_subset <- switch(input$qtl_model,
                          "Female" = female_mice,
                          "Male" = male_mice,
                          "Additive" = all_mice)
    
    add_covar <- if(input$qtl_model == "Additive") {
      covar[, c("ngen", "sex"), drop = FALSE]
    } else {
      covar[, "ngen", drop = FALSE]
    }
    
    common_mice <- Reduce(intersect, list(mice_subset, rownames(add_covar),
                                          dimnames(apr[[peak_chr]])[[1]],
                                          rownames(kloco[[peak_chr]])))
    
    pheno_subset <- cross2$pheno[common_mice, input$qtl_trait, drop = FALSE]
    complete_cases <- !is.na(pheno_subset[, 1])
    
    if(sum(complete_cases) < 10) {
      plot(1, type="n", axes=F, xlab="", ylab=""); text(1,1,"Not enough data for allele effects plot."); return()
    }
    final_mice <- common_mice[complete_cases]
    
    blup <- scan1blup(apr[final_mice, peak_chr],
                      cross2$pheno[final_mice, input$qtl_trait, drop = FALSE],
                      kloco[[peak_chr]][final_mice, final_mice],
                      addcovar = add_covar[final_mice, , drop = FALSE])
    
    plot_coefCC(blup, cross2$pmap[peak_chr], 
                main = paste("Allele Effects for", input$qtl_trait, "on Chr", peak_chr))
  })
}


# --- 5. RUN THE APP ---
shinyApp(ui = ui, server = server)

