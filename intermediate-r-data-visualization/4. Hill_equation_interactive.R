library(shiny)
library(ggplot2)

#Get a range of substrate concentration values. 
S <- 0:10000*0.01

#Half-maximal concentration constant.
K <- 50

# Define User Interface for app ----
ui <- pageWithSidebar(
  
  # App title ----
  headerPanel("Hill equation"),
  
  # Sidebar panel for inputs ----
  sidebarPanel(
    
    sliderInput(inputId = "Hill_coef", 
                label = "Hill coefficient:",
                min = 1, max = 50,
                value = 1)
    
  ),
  
  # Main panel for displaying outputs ----
  mainPanel(
    
    plotOutput(outputId = "Hill_figure")
    
  )
)

# Define server logic to plot various variables against mpg ----
server <- function(input, output) {
  
  
  output$Hill_figure <- renderPlot({
    i <- input$Hill_coef
    y_vals <- (S^i)/((S^i)+(K^i)) #For formula refer Wikipedia page for Hill coefficient.
    this_dat <- data.frame(Substrate = S, Rate = y_vals)
    
    ggplot(this_dat, aes(x = Substrate, y = Rate)) + 
      geom_line() +
      coord_cartesian(xlim = c(0, 100), ylim = c(0, 1))+
      theme(panel.border = element_rect(size = c(1,1,1,1), color = "black", 
                                        fill = NA))+
      ggtitle(paste("Hill coefficient =", i))
  })
  
}

shinyApp(ui, server)


