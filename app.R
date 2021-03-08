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

#Department 64 Data
load("Inputs/habitats_Dep64.RData")
Dep64_habitats <- spTransform(habitats, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#France-Belgium border data
load("Inputs/habitats_EA9.RData")
FrBel_Habitats <- spTransform(habitats, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

remove(habitats)

#Scenario Results

Scenario <- 4

for (Scen in Scenario) {
  
  load(paste0("./Outputs/Scenario_", Scen, "/", Scen, ".RData"))
  do.call(rbind, res) -> results
  
  res_newinf <- results[ ,c("NewInfAnimals", "NewInfCarcass", "NewInfGroups", "InfectedPerCells")]
apply(res_newinf, 2, function(x) do.call(rbind, x)) -> z
InfPerCells <- z$InfectedPerCells
InfPerCells <- apply(InfPerCells, 2, as.numeric)
colnames(InfPerCells) <- c("Iter", "gTime", "Cell_ID", "Inf_WB", "Pop_Cell")

}

InfPerCells4 <- InfPerCells

#InfPerCells[InfPerCells$V1 == 1, ] -> InfPerCells
ui <- fluidPage(
  
  titlePanel('Wild boar model results', windowTitle = "WB Model Results"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      radioButtons(inputId = 'scenario', label = 'Departmen 64 / Fra-Bel Border', 
                   choices = c('Department 64', 'Fra-Bel Border')),
      
      # Only show this panel if Department 64
      conditionalPanel(
        
        condition = "input.infection == 'Department 64'",
        radioButtons(inputId = 'hunting_Dep', label = "Hunting", choices = c("No Hunting", "Hunting Season"))
        
      ),
      
      conditionalPanel(
        
        condition = "input.infection == 'Fra-Bel Border'",
        radioButtons(inputId = 'hunting_FB', label = "Hunting", choices = c("No Hunting", "Hunting Season"))
        
      ),
      
    fluidRow(
      
      sliderInput(inputId = "days", label = "days:", min = 730, max = 1850, step = 1, value = 1, animate = TRUE)
       
      ),
   
    mainPanel(
      
        tabsetPanel(

            tabPanel(title = "MapI", leafletOutput("mapI"))

        ) #closes tabsetpanel
        
      ) #closes main panel
    
    ) #closes sidebar panel
    
  ) #closes sidebar layout (inside sidebar panel and main panel)
  
) #closes fluidPage

server <- function(input, output, session) {
  
  mapI <- leaflet(habitats) %>%
    addTiles() %>% 
    addPolygons(color = 'black', weight = 1, smoothFactor = 0.5,
                opacity = 1.0, fillColor = "green", fillOpacity = 0.5)
  
  output$mapI <- renderLeaflet(mapI)
  leafletProxy('mapI')
     observe({
      day<<-input$days
      
         habitats$InfStatus <- 0
         habitats$alpha <- 0
         tmp <- subset.data.frame(InfPerCells, InfPerCells[ ,"gTime"] == day)
         habitats$InfStatus[match(tmp[ ,3], habitats$ID)] <- 1
         habitats$alpha[match(tmp[ ,3], habitats$ID)] <- tmp[ ,4]/tmp[ ,5]
         col <- rev(heat.colors(10))[c(1, 10)]
         habitats$col <- col[habitats$InfStatus + 1]
       
         pal <- colorBin("plasma", bins = seq(0, 1, 0.1))
         leafletProxy("mapI") %>% addPolygons(data = habitats, color = "black",
                                              fillColor = ~pal(round(habitats$alpha, 1)),
                                              popup = paste("Infected WB:", habitats@data[ ,'alpha']), fillOpacity = 1)
         })
}


shinyApp(ui = ui, server = server)
