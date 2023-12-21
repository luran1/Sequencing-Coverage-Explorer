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
    titlePanel("Sequence Coverage Explorer"),
    
    
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
            conditionalPanel(condition = "input.conditionedPanels == 'Plots'",
                             sliderInput("bins",
                                         "Number of bins:",
                                         min = 1,
                                         max = 50,
                                         value = 30)),

            #Select Amplicon Sequence to review.
            conditionalPanel(condition = "input.conditionedPanels == 'Amplicon Viewer'",
                             uiOutput('Amplicon'),
                             sliderInput("Percent", 
                                         "minimum Base Accuracy",
                                         min = 0.1, value = 0.9, max = 1))
        ),    

        
        
        # Show a plot of the generated distribution and the Coverage table.
        mainPanel(
          tabsetPanel(
            type = "tabs",
            #Histogram plot
            tabPanel("Plots",plotOutput("distPlot"),
                     plotOutput("boxplot")),
            #Max Coverage table
            tabPanel("Amplicon Viewer", dataTableOutput("content")),
            tabPanel("Summary",dataTableOutput("summary")),
            id = "conditionedPanels"
          )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

# Reactive select menu for Amplicon Viewer
    output$Amplicon = renderUI({
      #input$upload will be NULL initially.
      req(input$upload)
      
      summary_stats <- read_tsv(input$upload$datapath)%>%
        select(`Ref Name`)
      
      selectInput('Amplicon_Name', 'Amplicon', summary_stats[,1])
    })
    
# Generate the selected Amplicon by base position with Accuracy.     
    output$content <- renderDataTable({
      
      #input$upload and input$Amplicon_Name will be NULL initially. 
      req(input$upload)
      req(input$Amplicon_Name)
      
      
      #Display only information from selected Amplicon.
      positions1 <- read_tsv(input$upload$datapath)%>%
        filter(`Ref Name` == input$Amplicon_Name)%>%
        select(!c(`Ref Name`,"Base 1", "Base 2", "Base 3", "Base 4", "Base 2 %", "Base 3 %", "Base 4 %") )%>%
        rename("Base Accuracy (%)" = "Base 1 %")%>%
        relocate("Base Accuracy (%)", .after = "Base")%>%
        relocate("Position", .before = "Base")%>%
        filter("Base Accuracy (%)" > input$Percent)%>%
        column_to_rownames(var = "Position")
      
      # Change Base Accuracy to % value
          
      #transpose the table to show each base as a column value.
      positions1 <- as.data.frame(t(positions1))
      
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

# Generate the distribution graph (box and whisker plot)
    output$boxplot <- renderPlot({
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
        distinct()%>%
        mutate("library"="library")%>%
        select(Coverage_Max,`library`)
      
      # draw box and whisker plot 
      p1 <- ggplot(positions1, aes(x=`library`, y=Coverage_Max)) +
        geom_boxplot() +
        stat_summary(fun.y = mean, geom="point", shape=10, size=2) +
        ggtitle("Distribution of Amplicon Coverage for Library")
      p1
      
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
