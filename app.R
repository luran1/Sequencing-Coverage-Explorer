#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#Loading in the required libraries 
library(shiny)
library(tidyverse)
library(DT)

# Define UI for application that draws a histogram and displayed the Coverage maximums.
ui <- fluidPage(

    # Application title
    titlePanel("Max Coverage Distribution"),
    
    
    # Sidebar with a file input button and slider input for number of bins. 
    sidebarLayout(
        sidebarPanel(
            
            #file upload button
            fileInput("upload",
                    "Position file upload",
                    buttonLabel = "upload...",
                    multiple = FALSE,
                    accept = ".tsv"),
            
            #slider input for number of bins.
            sliderInput("bins",
                        "Number of bins:",
                        min = 1,
                        max = 50,
                        value = 30)

        ),
        
        
        # Show a plot of the generated distribution and the Coverage table.
        mainPanel(
          tabsetPanel(
            type = "tabs",
            #Histogram plot
            tabPanel("Plot",plotOutput("distPlot")),
            #Max Coverage table
            tabPanel("Data", dataTableOutput("content")),
            tabPanel("Summary",dataTableOutput("summary"))
          )
           
        )
        
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  # Generate the Coverage table from positions file.     
    output$content <- renderDataTable({
      
      #input$upload will be NULL initially.
      req(input$upload)
       
      #transform the raw data to display only the Reference Name and the Coverage maximum.
      positions1 <- read_tsv(input$upload$datapath)%>%
        select("Ref Name", Coverage)%>%
        group_by(`Ref Name`)%>%
        mutate(Coverage_Max = max(Coverage))%>%
        slice_max(order_by = Coverage, n=1)%>%
        ungroup()%>%
        select(`Ref Name`, Coverage_Max)%>%
        distinct()
    })
  
  # Generate the distribution graph (histogram)    
    output$distPlot <- renderPlot({
      
      #input$upload will be NULL initially.
      req(input$upload)  
      
      #transform the raw data to display only the Reference Name and the Coverage maximum.
      positions1 <- read_tsv(input$upload$datapath)%>%
          select("Ref Name", Coverage)%>%
          group_by(`Ref Name`)%>%
          mutate(Coverage_Max = max(Coverage))%>%
          slice_max(order_by = Coverage, n=1)%>%
          ungroup()%>%
          select(`Ref Name`, Coverage_Max)%>%
          distinct()
      
      # take the Coverage Maximum list and determine bin values
      x <- positions1[,2]
      bins <- seq(min(x), max(x), length.out = input$bins + 1)

      # draw the histogram with the specified number of bins
      hist(as.numeric(unlist(x)), breaks = bins, labels = TRUE, col = 'darkgray', border = 'white',
           xlab = 'Coverage (reads)',
           main = 'Histogram of Coverage Maximums of Amplicons of Interest')
    })

  # Generate the Sequence Summary Table : Average Coverage,Average Q-score, STD for each, total amplicons in Lib
    output$summary <- renderDataTable({
      #input$upload will be NULL initially.
      req(input$upload) 
      summary_stats <- read_tsv(input$upload$datapath)%>%
        group_by(`Ref Name`)%>%
        mutate(Coverage_MAX = max(Coverage))%>%
        mutate(Coverage_AVG = mean(Coverage))%>%
        mutate("Q-Score_AVG" = mean(`Avg. Quality`))%>%
        mutate(Coverage_STD = sd(Coverage))%>%
        mutate("Q-Score_STD" = sd(`Avg. Quality`))%>%
        select(`Ref Name`,
               Coverage_MAX,
               Coverage_AVG,
               Coverage_STD,
               `Q-Score_AVG`,
               `Q-Score_STD`)%>%
        distinct()%>%
        mutate_if(is.numeric,round,digits = 2)
    })    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
