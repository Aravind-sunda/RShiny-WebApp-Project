#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Post Mortem Huntington Data Analysis"),

    tabsetPanel(
        tabPanel("File Input",
                 tabsetPanel(
                     tabPanel("Tab 1A", "Content for Tab 1A"),
                     tabPanel("Tab 1B", "Content for Tab 1B")
                 )
        ),
        tabPanel("Tab 2", "Content for Tab 2"),
        tabPanel("Tab 3",
                 tabsetPanel(
                     tabPanel("Tab 3A", "Content for Tab 3A"),
                     tabPanel("Tab 3B",
                              tabsetPanel(
                                  tabPanel("Tab 3B1", "Content for Tab 3B1"),
                                  tabPanel("Tab 3B2", "Content for Tab 3B2")
                              )
                     )
                 )
        )
    )

    # Sidebar with a slider input for number of bins
    # sidebarLayout(
    #     sidebarPanel(
    #         sliderInput("bins",
    #                     "Number of bins:",
    #                     min = 1,
    #                     max = 50,
    #                     value = 30)
    #     ),

        # Show a plot of the generated distribution
        # mainPanel(
        #    plotOutput("distPlot")
        # )

    # )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    source("functions.R")

}

# Run the application
shinyApp(ui = ui, server = server)
