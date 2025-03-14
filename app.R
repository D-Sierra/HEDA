# HLA Eplet-based Delisting Assistant (HEDA) isa program designed to 
# facilitate the delisting process of prohibited HLA alleles in highly 
# sensitized patients using eplet incompatibility as the base criterion.
# Copyright (C) 2024  Daniel Álvarez-Sierra and David San Segundo
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

library(shiny)
library(shinydashboard)
library(ggplot2)
library(tidyverse)
library(hlapro)
library(readxl)
library(magrittr)
library(gridExtra)
library(DT)
library(readr)
library(htmltools)
library(rmarkdown)
library(knitr)
library(kableExtra)
library(shinyjs)

inputs_confirmation <- c("Yes", "No", "All")
inputs_exposition <- c("High", "Intermediate", "Low", "Very Low", "All")

#Loading of all alleles that can be evaluated using Luminex for DL2/3 candidate selection
luminex_alleles <- readRDS(file.path(getwd(), "data", "luminex_alleles.rds"))


#App dashbord
ui <- dashboardPage(
  
  dashboardHeader(title = "HEDA"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data input", tabName = "Data_input", icon = icon("droplet"), startExpanded = T,
               menuSubItem('Input', tabName = "Input", icon = NULL),
               menuSubItem('Summary', tabName = "Summary", icon = NULL)),
      
      menuItem("Delisting 1", tabName = "DL1", icon = icon("list-check")),
      
      menuItem("Delisting 2/3", tabName = "DL23", icon = icon("list-check")),
      
      menuItem("Tutorial", tabName = "help", icon = icon("question-circle")),
      
      menuItem("Contact", tabName = "Contact", icon = icon("envelope"))
    )
  ),
  
  dashboardBody(
    useShinyjs(),
    tabItems(
      tabItem(tabName = "Input",
              fluidRow(
                box(
                    title = "HLA typing",
                    color = "black",
                    width = 12,
                    tags$p("HLA typing information for the recipient and previous donors must be provided in high resolution. Example: A*26:01, A*26:01, C*02:02,	C*07:02, B*07:02, B*27:05, DRB3*02:02, DRB3*02:02, DRB1*11:01, DRB1*14:01, DQB1*03:01, DQB1*05:03"),
                    textInput(inputId = "id_receptor", label = "Receptor ID", width = "200px", placeholder = "73497839"),
                    textInput(inputId = "receptorHLA", label = "Receptor HLa typing", placeholder = "A*26:01, A*01:01, C*02:02,	C*07:02, B*07:02, B*27:05, DRB3*02:02, DRB3*02:02, DRB1*11:01, DRB1*14:01, DQB1*03:01, DQB1*05:03"),
                    span(textOutput(outputId = "QC_receptor"),  style= "color:red"),
                    textInput(inputId = "donor1HLA", label = "Previous donors HLA typings", placeholder = "A*01:01,	A*26:01, C*05:01,	C*08:02, B*14:02,	B*44:02, DRB3*03:01,	DRB4*01:01, DRB1*04:01,	DRB1*14:01, DQB1*03:01,	DQB1*06:09"),
                    span(textOutput(outputId = "QC_donor1"),  style= "color:red"),
                    textInput(inputId = "donor2HLA", label = NULL, placeholder = "A*01:01, A*26:01, C*05:01,	C*08:02, B*14:02,	B*27:01, DRB1*04:01,	DRB1*11:01, DQB1*02:01,	DQB1*06:09"),
                    span(textOutput(outputId = "QC_donor2"),  style= "color:red"),
                    textInput(inputId = "donor3HLA", label = NULL, value = NA),
                    span(textOutput(outputId = "QC_donor3"),  style= "color:red"),
                    textInput(inputId = "donor4HLA", label = NULL, value = NA),
                    span(textOutput(outputId = "QC_donor4"),  style= "color:red"),
                    textInput(inputId = "donor5HLA", label = NULL, value = NA),
                    span(textOutput(outputId = "QC_donor5"),  style= "color:red"),
                )),
              fluidRow(  
                box(
                    title = "Delisting 1 (DL1) file",
                    color = "black",
                    width = 4,
                    tags$p("Anti-HLA antibody follow-up data allows for automatic analysis of candidate specificities to be delisted at the DL1 stage."),
                    fileInput(inputId = "fileDL1", label = "Anti-HLA antibody data", accept = c(".txt", ".csv")),
                    span(textOutput(outputId = "QC_alleles"),  style= "color:red; font-size:18px"),
                    span(textOutput(outputId = "QC_dates"),  style= "color:red; font-size:18px")
                  ),
                box(
                  title = "Complement-fixing donor specific antibodies (DSA)",
                  color = "black",
                  width = 8,
                  tags$p("Manual selection of complement-fixing DSA."),
                  tabBox(
                    width = NULL,
                    tabPanel("Locus A", checkboxGroupInput(inputId = "LocusA_compl_fix", label = "Locus A",choices = luminex_alleles[[1]], inline = TRUE)),
                    tabPanel("Locus B", checkboxGroupInput(inputId = "LocusB_compl_fix", label = "Locus B",choices = luminex_alleles[[2]], inline = TRUE)),
                    tabPanel("Locus C", checkboxGroupInput(inputId = "LocusC_compl_fix", label = "Locus C",choices = luminex_alleles[[3]], inline = TRUE)),
                    tabPanel("Locus DRB1", checkboxGroupInput(inputId = "LocusDRB1_compl_fix", label = "Locus DRB1",choices = luminex_alleles[[4]], inline = TRUE)),
                    tabPanel("Locus DRB345", checkboxGroupInput(inputId = "LocusDRB345_compl_fix", label = "Locus DRB345",choices = luminex_alleles[[5]], inline = TRUE)),
                    tabPanel("Locus DQB1", checkboxGroupInput(inputId = "LocusDQB1_compl_fix", label = "Locus DQB1",choices = luminex_alleles[[6]], inline = TRUE)),
                    tabPanel("Locus DQA1", checkboxGroupInput(inputId = "LocusDQA1_compl_fix", label = "Locus DQA1",choices = luminex_alleles[[7]], inline = TRUE)),
                    tabPanel("Locus DPB1", checkboxGroupInput(inputId = "LocusDPB1_compl_fix", label = "Locus DPB1",choices = luminex_alleles[[8]], inline = TRUE)),
                    tabPanel("Locus DPA1", checkboxGroupInput(inputId = "LocusDPA1_compl_fix", label = "Locus DPA1",choices = luminex_alleles[[9]], inline = TRUE))
                  )
                ))),
              
      tabItem(tabName = "Summary",
              fluidRow(
                box(
                  title = "Imputation of DQA1 locus",
                  status = "success",
                  solidHeader = TRUE,
                  width = 8,
                  DTOutput(outputId = "typing_table")
                ),
                column(
                  width = 4,
                  box(
                    title = "HLA missmatch with previous donors",
                    status = "warning",
                    solidHeader = TRUE,
                    width = 12,
                    tags$p("Unique alleles present in previous donors but not in the recipient, which will be considered to identify high-risk eplets."),
                    textOutput(outputId = "allele_MM")),
                  box(
                    title = "Complement-fixing DSA",
                    status = "danger",
                    solidHeader = TRUE,
                    width = 12,
                    textOutput(outputId = "alleles_compl")))),
              fluidRow(
                box(
                  title = "Eplet assigment",
                  status = "success",
                  solidHeader = TRUE,
                  width = 12,
                  tags$p("List of unique eplets present in the recipient and previous donors."),
                  tags$p("These listings are calculated according to the user-defined eplet filtering parameters in the Delisting 2/3 section."),
                  DTOutput(outputId = "eplet_table")))
              ),
      
      tabItem(tabName = "DL1",
              fluidRow(
                box(
                  width = 12,
                  title = "Configuración",
                column(
                  width = 5,
                  fluidRow(numericInput(inputId = "MFI_positive", label = "MFI positive cut-off", value = 5000, min = 0, max = 50000, step = 500)), 
                  fluidRow(tags$p("The selected value establishes the minimum MFI cut-off point to determine that the patient's serum has been positive for a particular HLA allele at least at some point in time in the assessed time series 
                                  (e.g. for the default value of 5000, all HLA alleles with MFIs always below this value will be considered negative, and will not be considered for delisting).")),
                  fluidRow(numericInput(inputId = "MFI_negative", label = "MFI negative cut-off", value = 5000, min = 0, max = 50000, step = 500)),
                  fluidRow(tags$p("The selected value determines the MFI cut-off point that the algorithm will take into account to identify HLA alleles that, after being positive at some point in time, have become negative and are considered for delisting. 
                                  (e.g. for the default value of 5000, only alleles that have been positive at some point in time and in recent months are below this value will be considered negative)."))),
                column(width = 1),
                column(
                  width = 5,
                  fluidRow(numericInput(inputId = "timeRange", label = "Time range", value = 24, min = 6, max = 48, step = 3)),
                  fluidRow(tags$p("Indicates the period of time during which a previously positive allele has to show a negative value (below the negative cut-off) to be considered a candidate for delisting.")),
                  fluidRow(radioButtons(inputId = "ignoreDP", label = "Ignore DP locus:", choices = c("Yes", "No"))),
                  fluidRow(tags$p("Option to ignore or take into account the DP locus for the assessment of DL1."))
                ))
                ),
              fluidRow(
                box(
                  title = "DELISTING - LEVEL 1",
                  width = 12,
                  status = "success",
                  solidHeader = TRUE,
                  fluidRow(
                    column(
                      width = 6,
                      DTOutput(outputId = "DL1_df"),
                      span(textOutput(outputId = "warningDL1"), style= "color:red; font-size:18px")
                    ),
                    column(
                      width = 6,
                      plotOutput(outputId = "DL1_plot")
                    ))
                )),
              fluidRow(  
                column(
                  width = 12,
                  align="center",
                  tags$p(style = "font-size: 16px",
                         tags$b("Download a PDF report with the results shown on this page.")),
                  column(
                    align = "right",
                    width = 6,
                    textInput(inputId = "responsableDL1", label = NULL, placeholder = "Name of the person responsible for the analysis")),
                  column(
                    align = "left",
                    width = 6,
                    downloadButton(outputId = "downloadDL1",label =  "Download DL1 results"))),
                  column(
                    width = 12,
                    align="center",
                    textOutput(outputId = "error_button_DL1"))
                )
      ),
      
      tabItem(tabName = "DL23",
              fluidRow(
                box(
                  title = "Configuration",
                  color = "black",
                  width = 12,
                  column(
                    tags$p("Selection of candidate specificities to be delisted."),
                    width = 8,
                    tabBox(
                      width = NULL,
                      tabPanel("Locus A", checkboxGroupInput(inputId = "LocusA_DL23", label = "Locus A",choices = luminex_alleles[[1]], inline = TRUE)),
                      tabPanel("Locus B", checkboxGroupInput(inputId = "LocusB_DL23", label = "Locus B",choices = luminex_alleles[[2]], inline = TRUE)),
                      tabPanel("Locus C", checkboxGroupInput(inputId = "LocusC_DL23", label = "Locus C",choices = luminex_alleles[[3]], inline = TRUE)),
                      tabPanel("Locus DRB1", checkboxGroupInput(inputId = "LocusDRB1_DL23", label = "Locus DRB1",choices = luminex_alleles[[4]], inline = TRUE)),
                      tabPanel("Locus DRB345", checkboxGroupInput(inputId = "LocusDRB345_DL23", label = "Locus DRB345",choices = luminex_alleles[[5]], inline = TRUE)),
                      tabPanel("Locus DQB1", checkboxGroupInput(inputId = "LocusDQB1_DL23", label = "Locus DQB1",choices = luminex_alleles[[6]], inline = TRUE)),
                      tabPanel("Locus DQA1", checkboxGroupInput(inputId = "LocusDQA1_DL23", label = "Locus DQA1",choices = luminex_alleles[[7]], inline = TRUE)),
                      tabPanel("Locus DPB1", checkboxGroupInput(inputId = "LocusDPB1_DL23", label = "Locus DPB1",choices = luminex_alleles[[8]], inline = TRUE)),
                      tabPanel("Locus DPA1", checkboxGroupInput(inputId = "LocusDPA1_DL23", label = "Locus DPA1",choices = luminex_alleles[[9]], inline = TRUE))
                  )),
                  column(
                    color = "black",
                    width = 4,
                    height = "400px",
                    radioButtons(inputId = "ab_verification", "Exclusive use of antibody-verified eplets:", inputs_confirmation),
                    radioButtons(inputId = "exposition", "Subset of eplets analysed according to their exposure:", inputs_exposition, selected = "All")
                  )
              )),
              fluidRow(
                box(
                  title = "Prohibited Eplets",
                  status = "warning",
                  width = 12,
                  solidHeader = TRUE,
                  tags$p("Unique eplets present in previous donors but not in the recipient that meet the selected criteria."),
                  textOutput(outputId = "eplet_MM"))
              ),
              fluidRow(
                box(
                  title = "DELISTING - LEVELS 2 AND 3",
                  width = 12,
                  status = "success",
                  solidHeader = TRUE,
                  DTOutput(outputId = "selectedDL23")
                )
              ),
              fluidRow(  
                column(
                  width = 12,
                  align="center",
                  tags$p(style = "font-size: 16px",
                         tags$b("Download a PDF report with the results shown on this page.")),
                  column(
                    align = "right",
                    width = 6,
                    textInput(inputId = "responsableDL23", label = NULL, placeholder = "Name of the person responsible for the analysis")),
                  column(
                    align = "left",
                    width = 6,
                    downloadButton(outputId = "downloadDL23",label =  "Download DL2/3 results"))
              ),
                column(
                  width = 12,
                  align="center",
                  textOutput(outputId = "error_button_DL23")
              ))
      ),
      
      tabItem(tabName = "help",
              tutorial()),
      
      tabItem(tabName = "Contact",
              style = "font-size: 16px",
              align = "center",
              tags$div("Daniel Álvarez de la Sierra", 
                       tags$br(),
                       "dalvarez@idival.org",
                       tags$br(),
                       "Immunology Group",
                       tags$br(),
                       "University Hospital Marqués de Valdecilla (HUMV)",
                       tags$br(),
                       "Santander (Spain)")))
  )
)



server <- function(input, output) {
  
  ##################################################################################################################
  ##########################################      QUALITY CONTROL       ############################################
  ##################################################################################################################
  
  #QC for alleles in input HLA typings 
  output$QC_receptor <- renderText({
    req(input$receptorHLA)
    check_alleles(input$receptorHLA)
  })
  output$QC_donor1 <- renderText({
    req(input$donor1HLA)
    check_alleles(input$donor1HLA)
  })
  output$QC_donor2 <- renderText({
    req(input$donor2HLA)
    check_alleles(input$donor2HLA)
  })
  output$QC_donor3 <- renderText({
    req(input$donor3HLA)
    check_alleles(input$donor3HLA)
  })
  output$QC_donor4 <- renderText({
    req(input$donor4HLA)
    check_alleles(input$donor4HLA)
  })
  output$QC_donor5 <- renderText({
    req(input$donor5HLA)
    check_alleles(input$donor5HLA)
  })

  #QC of the allele names contained in the file for DL1
  output$QC_alleles <- renderText({
    req(input$fileDL1)
    df_luminex <-  read.delim(input$fileDL1$datapath, sep = "\t", header = FALSE, na.strings	= "")
    DL1_QC_alleles(df_luminex)
  })
 
  #Quality check of the dates contained in the column names of the file for DL1
  output$QC_dates <- renderText({
    req(input$fileDL1)
    df_luminex <-  read.delim(input$fileDL1$datapath, sep = "\t", header = FALSE, na.strings	= "")
    DL1_QC_dates(df_luminex)
  })
  
  
  ##################################################################################################################
  ##########################################      DATA INPUT       #################################################
  ##################################################################################################################  
    
  #Determination of prohibited alleles according to previous donors
  allele_MM <- reactive({HLA_MM(receptor_HLA = input$receptorHLA,
                                donor1_HLA = input$donor1HLA,
                                donor2_HLA = input$donor2HLA,
                                donor3_HLA = input$donor3HLA,
                                donor4_HLA = input$donor4HLA,
                                donor5_HLA = input$donor5HLA)[[1]]})
  output$allele_MM <- renderText(allele_MM())
  
  #Reacitve vector that contains complement-fixing DSA
  selected_alleles_compl <- reactive({
    alleles <- c(input$LocusA_compl_fix, input$LocusB_compl_fix, input$LocusC_compl_fix, input$LocusDRB1_compl_fix, input$LocusDRB345_compl_fix, 
                 input$LocusDQB1_compl_fix, input$LocusDQA1_compl_fix, input$LocusDPB1_compl_fix, input$LocusDPA1_compl_fix)
    alleles <- alleles[!is.na(alleles)]
    alleles <- paste(alleles, collapse = ", ")})
  output$alleles_compl <- renderText(selected_alleles_compl())

  #Creation of a table with all the typings to check that they have been loaded correctly.
  output$typing_table <- renderDT(HLA_dataframe(receptor_HLA = input$receptorHLA,
                                                donor1_HLA = input$donor1HLA,
                                                donor2_HLA = input$donor2HLA,
                                                donor3_HLA = input$donor3HLA,
                                                donor4_HLA = input$donor4HLA,
                                                donor5_HLA = input$donor5HLA), options = list(scrollX = TRUE))
  
  #Reactive database that selects the eplets to be used according to the selected configuration
  eplet_database <- reactive({
    #Loading of eplets from HLA eplet registry
    #df_eplets <- hlapro::load_eplet_registry()                    #Code for loading by scraping from HLA Eplet Registry via hlapro
    df_eplets <- readRDS(file.path(getwd(), "data", "eplets.rds")) #Code for loading from local database
    req(input$ab_verification)
    req(input$exposition)
    #Filtering of eplet table according to the arguments ‘verified’ and ‘exposition’.
    verified = input$ab_verification
    exposition = input$exposition
    if(verified == "All") {verified <-  c("Yes", "No")}
    if(exposition == "All") {exposition <-  c("High", "Intermediate", "Low", "Very Low")}
    df_filtered <- df_eplets[df_eplets$confirmation %in% verified & df_eplets$exposition %in% exposition, ]
  })
  
  #Reactive table storing the unique eplets for the recipient and each of the donors
  eplet_table <- reactive({sample_eplet_assignment(hres_df = HLA_dataframe(receptor_HLA = input$receptorHLA,
                                                                           donor1_HLA = input$donor1HLA,
                                                                           donor2_HLA = input$donor2HLA,
                                                                           donor3_HLA = input$donor3HLA,
                                                                           donor4_HLA = input$donor4HLA,
                                                                           donor5_HLA = input$donor5HLA),
                                                   df_eplets = eplet_database())[, c("Sample","HLA_all_eplets")]})
  #Reactive element that collects the prohibited eplets for the recipient based on previous donors.
  eplet_missmatch <- reactive(eplet_MM(eplet_table()))
  
  #Output for single eplet vectors and for eplet missmatch for the receiver
  output$eplet_table <- renderDT({eplet_table()})
  output$eplet_MM <- renderText({eplet_missmatch()})
  

  ##################################################################################################################
  ##########################################      DELISTING       ##################################################
  ##################################################################################################################
  
  #Filtering of the file that contains luminex MFIs and saving as a reactive table
  DL1_file <- reactive({
    req(input$fileDL1)
    
    DL1_df <- read.delim(input$fileDL1$datapath, sep = "\t", header = FALSE, na.strings	= "")
  
    DL1_filtered <- DL1_filter(df_luminex = DL1_df,
                              allele_MM = allele_MM(),
                              allele_compl = selected_alleles_compl(),
                              cut_off_positive = input$MFI_positive,
                              cut_off_negative = input$MFI_negative, 
                              time_range = input$timeRange,
                              ignore_DP = input$ignoreDP)
  })
  #Reactive object for conditional colour formatting of DL1_file()
  DL1_file_format <- reactive(datatable(DL1_file(), options = list(scrollX = TRUE))%>% 
                                formatStyle(c('comp_fixation', 'allele_MM'), backgroundColor = styleEqual("Yes",'#f06565')))
                              
  #Creation of a new reactive table for the graphical representation of DL1 candidate specificities.
  DL1_file_plot <- reactive({
    req(DL1_file())
    DL1_df <- DL1_file()
    DL1_long <- DL1_df[,-(1:2)] %>%
      rownames_to_column("Alelos") %>%
      pivot_longer(-Alelos, names_to = "Fecha", values_to = "Valor") %>%
      arrange(as.Date(Fecha))
  })
  
  #Output table with candidate specificities and their MFIs over time
  output$DL1_df <- renderDT({tryCatch({
      DL1_file_format()
      }, error = function(e){
        datatable(data.frame())
    })
  })
  output$warningDL1 <- renderText(warningDL1_f(DL1_file())) #Error generation when nrow(DL1_file_format()) == 0
  
  #Plot showing the evolution of DL1 candidate specificities in MFIs
  DL1_plot_export <- reactive({
    ggplot( DL1_file_plot(), aes(x = as.Date(Fecha), y = Valor, color = Alelos, group = Alelos)) +
      geom_line(size = 1.5) +
      geom_point(shape = 1, size = 5) +
      labs(x = "Fecha",
           y = "MFI",
           color = NULL) +
      theme_bw() +
      theme(text = element_text(size = 12),
            axis.text = element_text(size = 12),
            plot.title = element_text(size = 14)) +
      geom_hline(yintercept = input$MFI_negative, linetype = "dashed", color = "gray", size = 1) })
  
  output$DL1_plot <- renderPlot(DL1_plot_export())
  
  #Reactive object that collects the selected candidate alleles for DL23 and converts it into a table
  selected_alleles_DL23 <- reactive({
    alleles <- c(input$LocusA_DL23, input$LocusB_DL23, input$LocusC_DL23, input$LocusDRB1_DL23, input$LocusDRB345_DL23, 
                 input$LocusDQB1_DL23, input$LocusDQA1_DL23, input$LocusDPB1_DL23, input$LocusDPA1_DL23)
    alleles <- alleles[!is.na(alleles)] #Remove possible NA
    #Empty vector to store the eplets and loop to populate it with the eplets associated to each allele according to the filter criteria
    eplets <- c()
    for (i in 1:length(alleles)){
      eplets[i] <- paste(hlapro::lookup_eplets(eplet_database(), alleles[i])[[1]], collapse = " ")
    }
    #Vector containing for each allele the list of eplets
    missmatch_eplet <- c()
    for (i in 1:length(eplets)){
      missmatch_eplet[i] <- forbidden_eplets(eplets[i], eplets =  eplet_missmatch())
    }
    #Vector indicating whether the allele is a candidate for DL2 or DL3
    dl_level <- c()
    for (i in 1:length(eplets)){
      dl_level[i] <- ifelse(missmatch_eplet[i] == "No eplet missmatch", "DL2", "DL3" )
    }
    #Dataframe for DL2 DL3
    datatable(data.frame(Alleles = alleles, Eplets = eplets, MM = missmatch_eplet, Delisting = dl_level)) %>% 
      formatStyle('Delisting', backgroundColor = styleEqual(c("DL2", "DL3"), c('#e4eb8a', '#f27024')))
  })
  
  output$alleles_DL23_table <- renderDataTable({
    selected_alleles_DL23()
  })
  
  #Same table as above but creates a df object instead of datatable so that it can be exported in the PDF report
  selected_alleles_DL23_df <- reactive({
    alleles <- c(input$LocusA_DL23, input$LocusB_DL23, input$LocusC_DL23, input$LocusDRB1_DL23, input$LocusDRB345_DL23, 
                 input$LocusDQB1_DL23, input$LocusDQA1_DL23, input$LocusDPB1_DL23, input$LocusDPA1_DL23)
    alleles <- alleles[!is.na(alleles)] #Remove possible NA
    
    eplets <- c()
    for (i in 1:length(alleles)){
      eplets[i] <- paste(hlapro::lookup_eplets(eplet_database(), alleles[i])[[1]], collapse = " ")
    }
    
    missmatch_eplet <- c()
    for (i in 1:length(eplets)){
      missmatch_eplet[i] <- forbidden_eplets(eplets[i], eplets = eplet_missmatch())
    }
    
    dl_level <- c()
    for (i in 1:length(eplets)){
      dl_level[i] <- ifelse(missmatch_eplet[i] == "No eplet missmatch", "DL2", "DL3")
    }
    
    data.frame(Alleles = alleles, Eplet_missmatch = missmatch_eplet, Delisting = dl_level)
  })
  
  #Output for eplet table
  output$selectedDL23 <- renderDT({tryCatch({
        selected_alleles_DL23()
      }, error = function(e){
        datatable(data.frame())
    })
  })
 
  ##################################################################################################################
  ##################################       GENERACION DE INFORMES       ############################################
  ##################################################################################################################
  
  #Generation of the report in PDF format for DL1
  output$downloadDL1 <- downloadHandler(
    filename = function() {
      ID_file_name <- input$id_receptor
      paste(ID_file_name, " DL1 ", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      #Generation of temporary file and copying of the template to the temporary directory
      tempReport <- file.path(tempdir(), "DL1_report.Rmd")
      file.copy("www/DL1_report.Rmd", tempReport, overwrite = TRUE)
      
      #Parameters that can be used as variables in the markdown file
      DL1_table <- DL1_file()
      params <- list(
        max = input$MFI_positive,
        min = input$MFI_negative,
        time = input$timeRange,
        ID = input$id_receptor,
        receptorHLA = input$receptorHLA,
        DL1_table = DL1_table,
        DL1_plot = DL1_plot_export(),
        responsable = input$responsableDL1
      )
      
      #PDF rendering
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv()))
    }
  )
  
  #Conditional so that the download button of the DL1 report is only enabled if the receiver ID, the name of the person responsible for the analysis and the DL1 file are provided
  observe({
    if (input$id_receptor == "" || input$responsableDL1 == "" || is.null(input$fileDL1)) {
      shinyjs::disable("downloadDL1")
      output$error_button_DL1 <- renderText("It is necessary to enter the recipient's ID, attach the file with the recipient's anti-HLA antibody historic data and the name of the person responsible for the analysis in order to generate the corresponding report.")
    } else {
      shinyjs::enable("downloadDL1")
      output$error_button_DL1 <- renderText("")
    }
  })
  
  
  #Generation of the report in PDF format for DL23
  output$downloadDL23 <- downloadHandler(
    filename = function() {
      ID_file_name <- input$id_receptor
      paste(ID_file_name, " DL2-3 ", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      #Generation of temporary file and copying of the template to the temporary directory
      tempReport <- file.path(tempdir(), "DL23_report.Rmd")
      file.copy("www/DL23_report.Rmd", tempReport, overwrite = TRUE)
      
      #Parameters that can be used as variables in the markdown file
      DL23_table <- selected_alleles_DL23_df()
      verified <- input$ab_verification
      exposition <- input$exposition
      params <- list(
        responsable = input$responsableDL23,
        ID = input$id_receptor,
        verified = verified,
        exposition = input$exposition,
        eplet_MM = eplet_missmatch(),
        DL23_table = DL23_table
      )
      
      #PDF rendering
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv()))
    }
  )

  #Conditional so that the download button of the DL1 report is only enabled if the receiver ID, the name of the person responsible for the analysis and the DL1 file are provided
  observe({
    alleles <- c(input$LocusA_DL23, input$LocusB_DL23, input$LocusC_DL23, input$LocusDRB1_DL23, input$LocusDRB345_DL23, 
                 input$LocusDQB1_DL23, input$LocusDQA1_DL23, input$LocusDPB1_DL23, input$LocusDPA1_DL23)
    if (input$id_receptor == "" || input$responsableDL23 == "" || length(alleles) == 0) {
      shinyjs::disable("downloadDL23")
      output$error_button_DL23 <- renderText("It is necessary to enter the recipient's ID, the name of the person responsible for the analysis and at least one candidate specificity in order to generate the corresponding report.")
    } else {
      shinyjs::enable("downloadDL23")
      output$error_button_DL23 <- renderText("")
    }
  })
}


shinyApp(ui = ui, server = server)
