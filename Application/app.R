#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

source("functions.R")

library(shiny)
library(shinydashboard)
library(colourpicker)
library(plotly)
library(DT)
library(shinyjs)
library(shinycssloaders)
library(shinyWidgets)
library(gplots)
library(RColorBrewer)


# Header title-----------------------------------------------------------------------------
header <- dashboardHeader(title = "Sequence Analysis")#,
                          #titleWidth = 450)



# Dashboard Sidebar-------------------------------------------------------------------------
sidebar <- dashboardSidebar(
    sidebarMenu(
        id = "tabs",
        menuItem("About", tabName = "about", icon = icon("question")),
        menuItem("Samples", tabName = "samples", icon = icon("vial")),
        menuItem("Counts",tabName = "counts", icon = icon("cubes")),
        menuItem("Differential Expression", tabName = "de", icon = icon("chart-line")),
        menuItem("Visualization of Gene Expression", tabName = "expression", icon = icon("dna"))
        # menuItem()
    )
)



# Main Body of data---------------------------------------------------------
body<- dashboardBody(
    tabItems(
        tabItem(tabName = "about",
                h1("mRNA-seq Data Analysis of patients with Huntington's Disease Post-mortem"),
                h3(strong("Credits")),
                p("This project was developed by Aravind Sundaravadivelu"),
                p("This data is obtained from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810"),
                h3(strong("Summary")),
                p("Huntington’s Disease (HD) is a devastating neurodegenerative disorder that is caused by an expanded CAG trinucleotide repeat in the Huntingtin (HTT) gene. Transcriptional dysregulation in the human HD brain has been documented but is incompletely understood. Here we present a genome-wide analysis of mRNA expression in human prefrontal cortex from 20 HD and 49 neuropathologically normal controls using next generation high-throughput sequencing. Surprisingly, 19% (5,480) of the 28,087 confidently detected genes are differentially expressed (FDR<0.05) and are predominantly up-regulated. A novel hypothesis-free geneset enrichment method that dissects large gene lists into functionally and transcriptionally related groups discovers that the differentially expressed genes are enriched for immune response, neuroinflammation, and developmental genes. Markers for all major brain cell types are observed, suggesting that HD invokes a systemic response in the brain area studied. Unexpectedly, the most strongly differentially expressed genes are a homeotic gene set (represented by Hox and other homeobox genes), that are almost exclusively expressed in HD, a profile not widely implicated in HD pathogenesis. The significance of transcriptional changes of developmental processes in the HD brain is poorly understood and warrants further investigation. The role of inflammation and the significance of non-neuronal involvement in HD pathogenesis suggest anti-inflammatory therapeutics may offer important opportunities in treating HD."),
                h3(strong("Citation")),
                p(
                tags$ol(
                    tags$li("Labadorf A, Hoss AG, Lagomarsino V, Latourelle JC et al. RNA Sequence Analysis of Human Huntington Disease Brain Reveals an Extensive Increase in Inflammatory and Developmental Gene Expression. PLoS One 2015;10(12):e0143563. PMID: 26636579"),
                    tags$li("Labadorf A, Choi SH, Myers RH. Evidence for a Pan-Neurodegenerative Disease Response in Huntington's and Parkinson's Disease Expression Profiles. Front Mol Neurosci 2017;10:430. PMID: 29375298"),
                    tags$li("Agus F, Crespo D, Myers RH, Labadorf A. The caudate nucleus undergoes dramatic and unique transcriptional changes in human prodromal Huntington's disease brain. BMC Med Genomics 2019 Oct 16;12(1):137. PMID: 31619230")

                )
                )

        ),
## Samples Tab--------
        tabItem(tabName = "samples",
                h1("Sample Analysis"),
                h3(strong("Information about page:")),
                tags$p("This page takes as the input the text file containing the sample information or the GEOID of the online repository and displays:"),
                tags$ol(
                    tags$li("Whole table"),
                    tags$li("Summarized Table"),
                    tags$li("Plots of necessary values")
                ),
                fluidRow(
                    box(title = "Input", collapsible = TRUE,
                        width = 12,
                        fileInput("samplefile","Choose Sample info file to read",accept = ".txt"),
                        textInput(inputId = "geoid",label = "Enter the GEO ID of the dataset",placeholder = "E.g. GSE64810"),
                        actionButton(inputId = "samplesubmit",label = "Submit")
                        ),
                    conditionalPanel("input.samplesubmit",
                        box(title = "Data explorer", collapsible = TRUE, width = 12,
                            withSpinner(dataTableOutput("table"))
                        ),
                        box(title = "Data Summary", collapsible = TRUE, width = 12,
                            textOutput("dimensions"),
                            withSpinner(dataTableOutput("summarytable"))
                        )
                    )
                    ),
                # Checks for the output of summarytable is present and then generated this tab.
                conditionalPanel("output.summarytable",
                                 fluidRow(
                                     box(title = "Generate Plots", collapsible = TRUE, width = 4,
                                         radioButtons(inputId = "plot_choice", "Plot Type",c("Histogram","Density","Violin"), inline = TRUE),
                                         splitLayout(
                                             cellWidths = c("50%", "50%"),
                                             uiOutput("column_selector"),
                                             radioButtons(inputId = "groupby", label = "Select a column to group by", choices = c("diagnosis","vonsattel grade"))
                                         ),
                                         actionButton(inputId = "plotsubmit",label = "Plot")
                                     ),
                                     box(title = "Plot", collapsible = TRUE, width = 8,
                                         plotlyOutput("sampleplot"))
                                     )
                                 )

        ),
## Counts Tab----------------
        tabItem(tabName = 'counts',
                h1("Sample Counts Analysis"),
                h3(strong("Information about page:")),
                tags$p("This page explores the count data in various ways and shows its information"),
                fluidRow(
                    box(title = "Input", collapsible = TRUE, width = 8,
                        fileInput("countsfile","Choose counts file to read",accept = ".txt"),
                        sliderInput(inputId = "varianceslider",label = "Filter out genes containing variance percentile less than:",min = 0 ,max = 100,value = 50, step = 0.1),
                        sliderInput(inputId = "zeroslider",label = "Filter out genes with zero counts greater than:",min = 0 ,max = 100,value = 50, step = 1),
                        actionButton(inputId = "countsubmit",label = "Submit")
                        ),
                    conditionalPanel("input.countsubmit > 0",
                                     box(title = "Summary of input", collapsible = TRUE, width = 4,
                                         uiOutput("countsummary")
                                         )
                                     ),
                    box(title = "Counts table", collapsible = TRUE, width = 12,
                        dataTableOutput("filteredcounts")
                        )
                    ),
                # conditionalPanel()
                fluidRow(
                    tabBox(
                        title = "Filtered Counts Analysis", id = "fcanal",side = "right", width = 12,height = 900,

                        tabPanel(title = "Scatter Plot",
                                 p(strong("This tab displays the scatter plot showing the genes that have been filtered out")),
                                 plotOutput("plot_variance"),
                                 plotOutput("plotzeros")
                                 ),
                        tabPanel(title = "Heatmap",
                                 fluidRow(
                                          column(width = 2, radioButtons("log_trans","Enabling log-transforming counts for visualization:",choices = c("TRUE", "FALSE")),
                                                 actionButton(inputId = "heatmapsubmit", label = "Submit")),
                                          column(width = 10,
                                                 conditionalPanel("input.heatmapsubmit",
                                                                  withSpinner(plotOutput("heatmapplot"))
                                                                  )
                                                 )

                                          )
                                 ),
                        tabPanel(title = "PCA plots",
                                 fluidRow(
                                     column(width = 4, uiOutput("pca_options")),
                                     column(width = 8, plotOutput("pcaplot"))
                                 )
                                 # splitLayout(
                                 #     cellWidths = c("30%", "70%"),
                                 #     uiOutput("pca_options"),
                                 #     plotOutput("pcaplot")
                                     # )
                                 )
                    )
                )

                ),
## Differential Expression Tab---------------
        tabItem(tabName = "de",
                # shinyAppDir("DESEQ")
                h1("DESeq Data Analysis"),
                h3(strong("Information about page:")),
                tags$p("This page takes in DESeq data and filters and displays the data in suitable ways"),
                fluidRow(
                    box(title = "Input",width = 3,collapsible = TRUE,
                        fileInput("deseqfile","Choose File to read",accept = ".txt"),
                        HTML("<p>A volcano plot can be generated with log2 fold-change on the x-axis and p-adjusted on the y-axis.</p>"),
                        radioButtons("x_axis", "Choose the column for X-axis", choices=c("baseMean","HD.mean","Control.mean","log2FoldChange","lfcSE","stat","pvalue","padj"),selected ="log2FoldChange" ),
                        radioButtons("y_axis", "Choose the column for Y-axis", choices=c("baseMean","HD.mean","Control.mean","log2FoldChange","lfcSE","stat","pvalue","padj"),selected ="padj" ),
                        colourpicker::colourInput(inputId = "base", label = "Base Point Color:", value = "black"),
                        colourpicker::colourInput(inputId = "highlight", label = "Highlight Point Color:", value = "#F56767"),
                        sliderInput(inputId = "slider",label = "Select the magnitude of the p adjusted coloring:",min = -100 ,max = 0,value = -10, step = 1),
                        div(actionButton(inputId = "desubmit", label = "Plot and Filter", style = "color: white; background-color: #0072B2;"), style = "text-align: center;")
                        ),
                    conditionalPanel("input.desubmit",
                                     tabBox(title = "Results of DESeq",width = 9,side = "right",selected = "Filtered Table",
                                            tabPanel("Volcano Plot",withSpinner(plotlyOutput("devolcano"))),
                                            tabPanel("Filtered Table",dataTableOutput("detable"))
                                            )
                                     )
                    )
                ),


## Individual Gene exploration tab----------------------------------
        tabItem(tabName = "expression",
                h1("Visualization of Individual Gene Expression(s)"),
                h3(strong("Information about page:")),
                tags$p("This page takes input from the previous pages and helps you to visualise the gene expression.Use this page only after you have used the samples and count page"),
                conditionalPanel("output.summarytable && output.filteredcounts ",
                                 fluidRow(box(title = "Input",width = 4,
                                              selectInput("genedrop", "Choose a variable to search by:", choices = c("age of death","age of onset","cag",
                                                                                                                 "diagnosis","duration","h-v cortical score",
                                                                                                                 "h-v striatal score","mrna-seq reads","pmi",
                                                                                                                 "rin","tissue","vonsattel grade")),
                                              # selectizeInput("gene_search", "Enter a gene:",
                                              #                choices = ensembl, multiple = FALSE,
                                              #                options = list(create = FALSE,
                                              #                               maxOptions = 10,
                                              #                               placeholder = "Search for a gene..."),
                                              #                selected = NULL),
                                              withSpinner(uiOutput("gene_search_ui")),

                                              actionButton("expressionsubmit",label = "Submit")
                                              ),
                                          conditionalPanel("input.expressionsubmit",
                                                           box(title = "Gene Expression Analysis",width = 8,
                                                               plotlyOutput("g_expression")
                                                               )
                                                           )
                                          )
                                 )
                )


        # tabItem() add comma
    )
)




ui <- dashboardPage(
    header,
    sidebar,
    body)






#Server-----------
server <- function(input, output,session) {
    source("functions.R")
    options(shiny.maxRequestSize = 50*1024^2)
    ## Page 1: Sample Data Exploration----------------------
    observeEvent(input$samplesubmit, {
        # Call the function that processes the input
        # This is a global variable
        if(!is.null(input$samplefile$datapath)){
            sampletable <<- GEO_data(filepath = file.path(input$samplefile$datapath))
        } else if(!is.null(input$geoid)){
            sampletable <<- GEO_data(geo_id = input$geoid)
        }
        else(showNotification("Please provide a file path or GEO ID", type = "warning"))

        summary <- summaryfunction(sampletable)
        n_rows <- nrow(sampletable)
        n_cols <- ncol(sampletable)
        dim_sentence <- paste("The table has", n_rows, "rows and", n_cols, "columns.")

        output$column_selector <- renderUI({
            radioButtons("column", "Select a column to plot",
                         choices = colnames(sampletable)[-(1:3)],
                         selected = colnames(sampletable)[4])
        })

        # Switch to the target tab and display the output
        # updateTabItems(session, "samples", selected = "output")
        output$table <- renderDataTable({
            sampletable
            }, width = "100%", options = list(scrollX = TRUE))

        output$summarytable <- renderDataTable({
            summary
        }, width = "100%", options = list(scrollX = TRUE))

        output$dimensions <- renderText({
            dim_sentence
        })
    })

    observeEvent(input$plotsubmit,{
        column<-input$column
        groupby <- input$groupby
        plotchoice <- input$plot_choice
        # shinyjs::toggle(id = "tab2")

        output$sampleplot <- renderPlotly({
            plot <- sampleplot(data = sampletable, column = column, group = groupby,plottype = plotchoice)
        })

    })

    ## Page 2 Counts data exploration-----------------
    observeEvent(input$countsubmit,{
        #Counts data
        countstable <<- countsread(file.path(input$countsfile$datapath))

        c_dim <- dim(countstable)
        #Filtered counts data
        filtered_data <<- filter_genes(data = countstable, var_percentile = input$varianceslider,zero_filter_value = input$zeroslider)

        f_dim <- dim(filtered_data)


        output$filteredcounts <- renderDataTable({
            filtered_data
        },width = "100%", options = list(scrollX = TRUE))

        # output$countsummary <- renderUI({
        #     tags$div(
        #         HTML(paste0("1) The original table had ", c_dim[1], " genes & ", c_dim[2], " samples <br>")),
        #         HTML(paste0("2) The filtered table had ", f_dim[1], " genes & ", f_dim[2], " samples <br>")),
        #         HTML(paste0("3) The number of genes filtered is ", (c_dim[1] - f_dim[1]), "<br>")),
        #         HTML(paste0("4) Percentage of genes that passed: ", (f_dim[1] / c_dim[1]) * 100))
        #     )
        # })

        output$countsummary <- renderUI({

            # Create a data frame with summary information
            summary_table <- data.frame(
                Metric = c("Original genes", "Original samples", "Filtered genes", "Filtered samples", "Number of genes filtered", "Percentage of genes passed"),
                Value = c(c_dim[1], c_dim[2], f_dim[1], f_dim[2], (c_dim[1] - f_dim[1]), round((f_dim[1] / c_dim[1]) * 100, 2))
            )

            # Convert the data frame to a table using the DT::datatable() function
            datatable(summary_table, rownames = FALSE)

        })




        output$pca_options <- renderUI({
            pca_object <<- performPCA(filtered_data)
            my_list <- pca_object$components
            new <- list(
                selectInput("PC_A", "Select 1st PC to plot:", choices = my_list),
                selectInput("PC_B", "Select 2nd PC to plot:", choices = my_list),
                actionButton("PCsubmit","Submit")
            )

        })

        output$plot_variance <- renderPlot({
            med_vs_var(counts_tib = countstable,perc_var = input$varianceslider)
        })


        output$plotzeros <- renderPlot({
            med_vs_nz(#counts_tib = , nz_genes = input$zeroslider
                      dataf = countstable, slider_var = input$varianceslider, slider_num = input$zeroslider)
        })

        observeEvent(input$PCsubmit,{
            output$pcaplot <- renderPlot({
                pcaplot(pca_obj = pca_object,pca = input$PC_A, pcb = input$PC_B, loadings = FALSE )
            })
        })



        observeEvent(input$heatmapsubmit,{
            output$heatmapplot <- renderPlot({
                plot_heatmap(filtered_data,input$log_trans)},width = 1000, height = 600)

        })



    })

    ## Page 3 DESEQ exploration-----------------
    observeEvent(input$desubmit,{
        mydata <- load_data(path = input$deseqfile$datapath)

        output$devolcano <- renderPlotly({
            volcano_plot(dataf = mydata,x_name = input$x_axis,y_name = input$y_axis,slider = input$slider,color1 = input$base, color2 = input$highlight)
            })

        output$detable <- renderDataTable({
            draw_table(dataf = mydata,slider = input$slider)
            },width = "100%", options = list(scrollX = TRUE))
    })

    ## Page 4 Individual Gene exploration-----------------
    update_gene_search_choices <- reactive({
        selectizeInput(inputId = "gene_search", label = "Enter a gene:",
                       choices = ensembl, multiple = FALSE,
                       options = list(create = FALSE, maxOptions = 10,
                                      placeholder = "Search for a gene..."),
                       selected = NULL)
    })

    output$gene_search_ui <- renderUI({
        update_gene_search_choices()
    })

    observeEvent(input$expressionsubmit,{
        output$g_expression <- renderPlotly({
            geneexpression(sampledata = sampletable,countsdata = countstable, gene = input$gene_search,variable = input$genedrop)
        })
    })

    # Adding a reactive expression to check the conditions
    # if(output$table == TRUE && output$filteredcounts == TRUE){
    #     panel_visible <- reactive({
    #         output$table && output$filteredcounts
    #     })
    # }
    #
    # observeEvent(panel_visible(), {
    #     if (panel_visible()) {
    #         ensembl <<- countstable[[1]]
    #         output$gene_search <- renderUI({
    #             selectizeInput("gene_search", "Enter a gene:",
    #                            choices = ensembl, multiple = FALSE,
    #                            options = list(create = FALSE,
    #                                           maxOptions = 10,
    #                                           placeholder = "Search for a gene..."))
    #         })
    #     }
    # })



    # WORKS!!!
    # observeEvent(input$samplesubmit, {
    #     updateTabItems(session, "sidebar", "expression")
    #     output$expression_content <- renderUI({
    #         fluidRow(
    #             column(12,
    #                    checkboxInput("checkbox", "Check me!"))
    #         )
    #     })
    # })
    #UI here!!!
    # fluidRow(
    #     box(title = "Input Check",collapsible = TRUE,
    #         # uiOutput("checkbox1_ui")
    #         uiOutput("expression_content")
    #
    #     )
    # )




}

# Run the application
shinyApp(ui = ui, server = server)


## Notes
