## Author: Aravind Sundaravadivelu
## saravind@bu.edu
## BU BF591
## Assignment 7

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(magrittr)
library(tidyverse)
library(ggplot2)
library(colourpicker) # you might need to install this.
library(glue)


# Define UI for application that draws a histogram
ui <- fluidPage(
  title = "DESeq Data Analyser",
  theme ="bootstrap.css",
  titlePanel("BF591 Assignment 7"),
  sidebarLayout(
    sidebarPanel(
      HTML("<p>Load differential expression results</p>"),
      fileInput("file","Choose File to read",accept = ".csv"),
      HTML("<p>A volcano plot can be generated with log2 fold-change on the x-axis and p-adjusted on the y-axis.</p>"),
      radioButtons("x_axis", "Choose the column for X-axis", choices=c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"),selected ="log2FoldChange" ),
      radioButtons("y_axis", "Choose the column for Y-axis", choices=c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"),selected ="padj" ),
      colourpicker::colourInput(inputId = "base", label = "Base Point Color:", value = "black"),
      colourpicker::colourInput(inputId = "highlight", label = "Highlight Point Color:", value = "#F56767"),
      sliderInput(inputId = "slider",label = "Select the magnitude of the p adjusted coloring:",min = -300 ,max = 0,value = -150, step = 1),
      div(actionButton(inputId = "Plot", label = "Plot and Filter", style = "color: white; background-color: #0072B2;"), style = "text-align: center;")
    ),
    mainPanel(
      tabsetPanel(tabPanel("Volcano Plot",plotOutput("volcano")), tabPanel("Filtered Table",tableOutput("table")))
    )
  )
)
# tableOutput()

# Define server logic required to draw a histogram
server <- function(input, output, session) {

    #' load_Data
    #'
    #' @details Okay this one is a little weird but bear with me here. This is
    #' still a "function", but it will take no arguments. The `reactive({})` bit
    #' says "if any of my inputs (as in, input$...) are changed, run me again".
    #' This is useful when a user clicks a new button or loads a new file. In
    #' our case, look for the uploaded file's datapath argument and load it with
    #' read.csv. Return this data frame in the normal return() style.
    load_data <- reactive({
      req(input$file)
      file <- read_csv(input$file$datapath) %>%
        dplyr::rename("gene" = "...1")
      return(file)
    })

    #' Volcano plot
    #'
    #' @param dataf The loaded data frame.
    #' @param x_name The column name to plot on the x-axis
    #' @param y_name The column name to plot on the y-axis
    #' @param slider A negative integer value representing the magnitude of
    #' p-adjusted values to color. Most of our data will be between -1 and -300.
    #' @param color1 One of the colors for the points.
    #' @param color2 The other colors for the points. Hexadecimal strings: "#CDC4B5"
    #'
    #' @return A ggplot object of a volcano plot
    #' @details I bet you're tired of these plots by now. Me too, don't worry.
    #' This is _just_ a normal function. No reactivity, no bells, no whistles.
    #' Write a normal volcano plot using geom_point, and integrate all the above
    #' values into it as shown in the example app. The testing script will treat
    #' this as a normal function.
    #'
    #' !!sym() may be required to access column names in ggplot aes().
    #'
    #' @examples volcano_plot(df, "log2fc", "padj", -100, "blue", "taupe")
    volcano_plot <-
        function(dataf, x_name, y_name, slider, color1, color2) {
          dataf <- dataf %>% drop_na()
          plot <- ggplot(data = dataf, aes(x = !!sym(x_name), y = -log10(!!sym(y_name))))+
            geom_point(aes(color = if_else(!!sym(y_name) < 1*10^slider, TRUE,FALSE)))+
            labs(title="DESeq Results", x = glue('{x_name}'), y = glue('-log10({y_name})')) +
            scale_color_manual(name = glue("{x_name} < 1*10^{slider}"),
                               values = c("FALSE" = color1,
                                          "TRUE" = color2),
                               labels = c("FALSE", "TRUE"))+
            theme_classic()+
            theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")

            return(plot)
        }



    #' Draw and filter table
    #'
    #' @param dataf Data frame loaded by load_data()
    #' @param slider Negative number, typically from the slider input.
    #'
    #' @return Data frame filtered to p-adjusted values that are less than
    #' 1 * 10^slider, columns for p-value and p-adjusted value have more digits
    #' displayed.
    #' @details Same as above, this function is a standard R function. Tests will
    #' evaluate it normally. Not only does this function filter the data frame to
    #' rows that are above the slider magnitude, it should also change the format
    #' of the p-value columns to display more digits. This is so that it looks
    #' better when displayed on the web page. I would suggest the function
    #' `formatC()`
    #'
    #' @examples draw_table(deseq_df, -210)
    #'    X  baseMean     log2FC     lfcSE      stat       pvalue         padj
    #'gene1 11690.780   9.852926 0.2644650  37.25607 8.45125e-304 1.54472e-299
    #'gene2  3550.435  -6.183714 0.1792708 -34.49369 9.97262e-261 9.11398e-257
    draw_table <- function(dataf, slider) {
        filtered_table <- dataf %>% dplyr::filter(padj <= 1*10^slider) %>%
          drop_na()
        # filtered_table$padj <- sprintf("%.2f", filtered_table$padj)
        filtered_table$padj <- format(filtered_table$padj, digits = 6)
        filtered_table$pvalue <- format(filtered_table$pvalue, digits = 6)
        return(filtered_table)
    }

    #' These outputs aren't really functions, so they don't get a full skeleton,
    #' but use the renderPlot() and renderTabel() functions to return() a plot
    #' or table object, and those will be displayed in your application.
    # observeEvent(input$Plot,{
    #
    #     isolate({mydata <- load_data()})})

    output$volcano <- renderPlot({
      input$Plot
      mydata <- load_data()
      isolate({volcano_plot(dataf = mydata,x_name = input$x_axis,y_name = input$y_axis,slider = input$slider,color1 = input$base, color2 = input$highlight)})})

    output$table <- renderTable({
      input$Plot
      mydata <- load_data()
      isolate({draw_table(dataf = mydata,slider = input$slider)})})
    # renderTable()

    # {output$volcano <- NULL}
    # {output$table <- NULL}

    # observeEvent(input$Plot,{
    #
    #   isolate({mydata <- load_data()})
    #   if (is.null(data)) return (NULL)
    #     output$volcano <- renderPlot({
    #       isolate({volcano_plot(dataf = mydata,x_name = input$xaxis,y_name = input$yaxis,slider = input$slider,color1 = input$color1, color2 = input$color2)})
    #         })
    #     output$table <- renderDataTable({
    #       isolate({draw_table(dataf = mydata,slider = input$slider)})
    #         })
    # })


}

# Run the application
shinyApp(ui = ui, server = server)
