#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(leaflet)
library(igraph)
library(rgdal)
library(maptools)
library(sp)
library(rgeos)
library(viridis)
file="Essai"
# source('Modele non parallelise©_MA_map.R', encoding = 'UTF-8')
# source("carteinteractive.R")
load("habitats_Dep64.RData")
habitats <- spTransform(habitats, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
InfPerCells<-read.table(paste(file,"-FDayOutInfPerCell.txt",sep=''),sep=' ')
InfPerCells[InfPerCells$V1==1,]->InfPerCells
ui <- fluidPage(
    # leafletOutput("map1"),
    # plotOutput("distPlot"),
    # 
    # p(),
    hr(), 
    # tmp<<- as.Date(sample(df$date,1)),
    
    # Sidebar with a slider input
    # sidebarPanel(
    fluidRow(
      sliderInput(inputId = "days", label = "days:", min = 730, max = 1850, step = 1, value = 1, animate = TRUE)
       ),
    # selectInput("OutputSelection", "Choose a output type:",c("Map","Prevalence","Within Herd Dynamics"),width=200)),
    mainPanel(
        tabsetPanel(

            tabPanel(title = "MapI", leafletOutput("mapI"))

      )
    )
    #   leafletOutput("map1"),width = 9)
)

server <- function(input, output, session) {
  
    # 
    # output$dayOut <- renderUI({
    #     sliderInput(inputId = "days", label = "days:", min = 730, max = 1850, step = 30, value = 1, animate = TRUE)
    #     # sliderInput("days", "days:",1,60,1)
    # })
    # 
  mapI <- leaflet(habitats) %>%
    addTiles() %>% 
    addPolygons(color = 'black', weight = 1, smoothFactor = 0.5,
                opacity = 1.0, fillColor = "green", fillOpacity = 0.5)
  
  output$mapI <- renderLeaflet(mapI)
  leafletProxy('mapI')
    # prevalence<<-NULL
     observe({
      day<<-input$days
      


         habitats$InfStatus <- 0
         habitats$alpha <- 0
         tmp <- subset.data.frame(InfPerCells,InfPerCells$V2 == day)
         habitats$InfStatus[match(tmp[,3],habitats$ID)] <- 1
         habitats$alpha[match(tmp[,3],habitats$ID)] <- tmp[,4]/tmp[,5]
         col <- rev(heat.colors(10))[c(1,10)]
         habitats$col <- col[habitats$InfStatus+1]
         # habitats$col[habitats$alpha>0] <- adjustcolor(habitats$col[habitats$alpha>0],
         #                                            habitats$alpha[habitats$alpha>0])
       
         pal <- colorBin("plasma", bins = seq(0,1,0.1))
         leafletProxy("mapI") %>% addPolygons(data = habitats, color="black",
                                              fillColor = ~pal(round(habitats$alpha,1)),
                                              popup = paste("Infected WB:", habitats@data[,'alpha']), fillOpacity =1)
         })
}


shinyApp(ui = ui, server = server)
