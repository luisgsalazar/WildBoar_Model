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

load("Inputs/habitats_Dep64.RData")
#load("Inputs/habitats_EA9.RData")
habitats <- spTransform(habitats, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))


for (Scen in 1:1) {
  
  load(paste0("./Outputs/Scenario_", Scen, "/", Scen, ".RData"))
  do.call(rbind, res) -> results
  
  res_newinf <- results[ ,c("NewInfAnimals", "NewInfCarcass", "NewInfGroups", "InfectedPerCells")]
  apply(res_newinf, 2, function(x) do.call(rbind, x)) -> z
  InfPerCells <- z$InfectedPerCells
  InfPerCells <- apply(InfPerCells, 2, as.numeric)
  colnames(InfPerCells) <- c("Iter", "gTime", "Cell_ID", "Inf_WB", "Pop_Cell")
  InfPerCells <- as.data.frame(InfPerCells)
}

ui <- fluidPage(
  
  titlePanel('Wild boar Model Results', windowTitle = "WB Model Results"),    
  
  hr(), 
  
 sidebarLayout(
   
  # Sidebar with a slider input
  sidebarPanel(
    
  radioButtons(inputId = 'scenario', label = 'Department 64 / France-Belgium Border', 
                 choices = c('Department 64', 'France-Belgium Border')),
  
  conditionalPanel(condition = "inpute.scenario = 'Department 64'",
                   radioButtons(inputId = 'season', label = "Season", choices = c("ASF introduction Winter", "ASF introduction Summer"))
  
  ),
                   conditionalPanel(
              
              condition = "input.season == 'ASF introduction Winter'",
              radioButtons(inputId = 'hunting_Dep', label = "Hunting", choices = c("No Hunting", "Hunting Season"))
              
            ),
            
            conditionalPanel(
              
              condition = "input.season == 'ASF introduction Summer'",
              radioButtons(inputId = 'hunting_FB', label = "Hunting", choices = c("No Hunting", "Hunting Season"))
              
            ),
    
  # fluidRow(
  #     sliderInput(inputId = "days", label = "days:", min = 730, max = 1850, step = 1, value = 1, animate = TRUE)
  #      ),
  
  uiOutput(outputId = 'DayOut')
   
 ), #Sidebar Panel
    mainPanel(
      
        tabsetPanel(

            tabPanel(title = "Infected WB Map", leafletOutput(outputId = "mapI"))

      )
    )
  )#closes Sidebar layout (inside sidebar panel (where inputs are) + main panel)
 )

server <- function(input, output, session) {
  
  output$DayOut <- renderUI({
    
    sliderInput(inputId = "days", label = "days:", min = 730, max = 1850, step = 1, value = 1, animate = TRUE)
    
  })
  
  mapI <- leaflet(habitats) %>%
    addTiles() %>% 
    addPolygons(color = 'black', weight = 1, smoothFactor = 0.1,
                opacity = 1.0, fillColor = "green", fillOpacity = 0.5)
  
  output$mapI <- renderLeaflet(mapI)
  leafletProxy('mapI')

     observe({
       
      day <<- input$days

         habitats$InfStatus <- 0
         habitats$alpha     <- 0
         tmp <- subset.data.frame(InfPerCells, InfPerCells$gTime == day)
         habitats$InfStatus[match(tmp[ ,3], habitats$ID)] <- 1
         habitats$alpha[match(tmp[ ,"Cell_ID"], habitats$ID)] <- (tmp[ ,"Inf_WB"]/tmp[ ,"Pop_Cell"]*100)
         col <- rev(heat.colors(10))[c(1, 10)]
         habitats$col <- col[habitats$InfStatus + 1]
         # habitats$col[habitats$alpha>0] <- adjustcolor(habitats$col[habitats$alpha>0],
         #                                            habitats$alpha[habitats$alpha>0])
       
         pal <- colorBin("plasma", bins = seq(0, 100, 0.10), alpha = 0.5)
         leafletProxy("mapI") %>% addPolygons(data = habitats, color = "black",
                                              fillColor = ~pal(round(habitats$alpha, 1)),
                                              popup = paste("Infected WB:", habitats@data[ ,'alpha']), fillOpacity = 2)
         })
}


shinyApp(ui = ui, server = server)
