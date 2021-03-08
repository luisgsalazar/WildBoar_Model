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
    
    # Sidebar with a slider input
    # sidebarPanel(
    fluidRow(
      sliderInput(inputId = "days", label = "days:", min = 730, max = 1850, step = 1, value = 1, animate = TRUE)
       ),
   
    mainPanel(
        tabsetPanel(

            tabPanel(title = "MapI", leafletOutput("mapI"))

      )
   )
)

server <- function(input, output, session) {
  
  mapI <- leaflet(habitats) %>%
    addTiles() %>% 
    addPolygons(color = 'black', weight = 1, smoothFactor = 0.1,
                opacity = 1.0, fillColor = "green", fillOpacity = 0.5)
  
  output$mapI <- renderLeaflet(mapI)
  leafletProxy('mapI')

     observe({
       
      day <<- input$days

         habitats$InfStatus <- 0
         habitats$alpha <- 0
         tmp <- subset.data.frame(InfPerCells, InfPerCells$gTime == day)
         habitats$InfStatus[match(tmp[ ,3], habitats$ID)] <- 1
         habitats$alpha[match(tmp[ ,3], habitats$ID)] <- tmp[ ,4]/tmp[ ,5]
         col <- rev(heat.colors(10))[c(1, 10)]
         habitats$col <- col[habitats$InfStatus + 1]
         # habitats$col[habitats$alpha>0] <- adjustcolor(habitats$col[habitats$alpha>0],
         #                                            habitats$alpha[habitats$alpha>0])
       
         pal <- colorBin("plasma", bins = seq(0, 1, 0.1), alpha = 0.3)
         leafletProxy("mapI") %>% addPolygons(data = habitats, color = "black",
                                              fillColor = ~pal(round(habitats$alpha, 1)),
                                              popup = paste("Infected WB:", habitats@data[ ,'alpha']), fillOpacity = 2)
         })
}


shinyApp(ui = ui, server = server)
