# LIBRARY ####
library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyWidgets)
library(DT)
library(readr)
library(leaflet)
library(mapview)
library(plyr)
library(dplyr)
library(ggplot2)
library(sp)
library(stringr)
library(data.table)
library(rlang)
library(sf)
library(RColorBrewer)
#

#Note: If map fails to load need to install older version of mapview
#library(devtools)
#install_version('mapview',version='2.7.8')
#install_version('rlang',version='0.4.9')

# Get Data ----------------------------------------------------------------

MyTab <- fread("data_with_effect_sizes.csv")
#MyTab <- fread("shiny/ArcticHerbivoreDiversitySystematicReview/data_with_effect_sizes.csv")#Turn off in app


MyTab$coordinates_E<-as.numeric(MyTab$coordinates_E)
MyTab$coordinates_N<-as.numeric(MyTab$coordinates_N)

#Remove points with NA in coordinates
MyTab<-MyTab[!is.na(MyTab$coordinates_E),]
#Remove points with NA in yi
MyTab<-MyTab[!is.na(MyTab$yi_smd),]

#Truncate for plotting
MyTab$yi_smd[MyTab$yi_smd>6]<-6
MyTab$yi_smd[MyTab$yi_smd< -6]<-(-6)

#Rounding and truncating text for popups
MyTab$EffectSize<-round(MyTab$yi_smd,1)
MyTab$FirstAuthor<-sapply(strsplit(MyTab$author_list," "), `[`, 1)


#List variables to play with

varA <- c("year",
          "change.long_f",
          "measured_response_variable_new",
          "country",
          "study_design",
          "diversity_contrast"
)


# UI ----------------------------------------------------------------------

ui <- dashboardPage(
  
   dashboardHeader(title = "Arctic Herbivore Diversity Systematic Review",
                  titleWidth = 400),
  
  dashboardSidebar(disable = T),
 
  dashboardBody(
    
    # Include custom CSS to add padding and auto-opacity to filtering-box
    tags$head(
      includeCSS("styles.css")
    ),
    
    
    
    # The Map -----------------------------------------------------------------
    
    shinydashboard::box(width = NULL,
                        leafletOutput('theMap')),
    
    # Abs.Panel ---------------------------------------------------------------
    
    absolutePanel(id = "controls", 
                  class = "panel panel-default", 
                  fixed = F,
                  draggable = TRUE, 
                  top = 80, left = 70, right = 'auto', bottom = "auto",
                  width = 300, height = "auto",
                  
                  # Changing the colour of the slider from blue to orange
                  tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: orange} .irs-from, .irs-to, .irs-single {background: orange }")),   
                  
                  
                  # . Year and Colour ---------------------------------------------------------
                  
                  
                  fluidRow(
                    column(width = 7,
                           br(),
                           switchInput(
                             inputId = "filteron",
                             label = "Filters", 
                             value = TRUE,
                             labelWidth = "80px")
                    )#,
                    #column(width=5,
                    #       pickerInput(inputId = "colour", 
                    #                   label = "Colour by:", 
                    #                   choices =varA,
                    #                   selected = "country",
                    #                   options = list(size=10)))
                    ),
                  
                  
                  
                  
                  # . Filters -----------------------------------------------------------------
                  
                  sliderInput("year", 
                              label = h5("Publication year"), 
                              min = 1940, #min(MyTab$year), 
                              max = 2021, #min(MyTab$year), 
                              value = c(1940, 2021),
                              step=1, sep="", ticks=F, width = "100%",
                              animate=F
                  ),
                  
                  # fluidRow(
                  column(width = 6,
                         # .. Country ---------------------------------------------------------------
                         
                         pickerInput(
                           inputId = "country",
                           choices = unique(MyTab$country),
                           selected = unique(MyTab$country),
                           multiple = TRUE,
                           width = "100%",
                           options = list(
                             `actions-box` = TRUE, `selected-text-format`= "static", title = "Country",size=6)),
                         
                         
                           
                         
                         # .. study design ----------------------------------------------------------
                         
                         pickerInput(inputId = 'study_design',
                                     choices = unique(MyTab$study_design),
                                     selected = unique(MyTab$study_design),
                                     multiple = TRUE,
                                     width = "100%",
                                     options = list(
                                       `actions-box` = TRUE, 
                                       `selected-text-format`= "static",
                                       title = "Study Design")
                         ),
                         
                         
                         
                         # .. diversity contrast ----------------------------------------------------------
                         
                         pickerInput(
                           inputId = 'diversity_contrast',
                           choices = unique(MyTab$diversity_contrast),
                           selected = unique(MyTab$diversity_contrast),
                           multiple = TRUE,
                           width = "100%",
                           options = list(
                             `actions-box` = TRUE,
                             `selected-text-format`= "static",
                             title = "Diversty contrast")
                         ),

                  
                  
                  # .. measured response vairable ----------------------------------------------------------
                  
                  pickerInput(
                    inputId = 'measured_response_variable_new',
                    choices = unique(MyTab$measured_response_variable_new),
                    selected = unique(MyTab$measured_response_variable_new),
                    multiple = TRUE,
                    width = "100%",
                    options = list(
                      `actions-box` = TRUE, 
                      `selected-text-format`= "static",
                      title = "Response Variable")
                  ),
                       
                         
                         # Remaining datapoints ----------------------------------------------------
                         
                         verbatimTextOutput('remaining')
                         
                         
                  ) # column
                  
    ), # abs.panel1
    
    
    # Description text and Active Filters ####
    h5("This is an interactive user interface of the systematic review of herbivore diversity impacts in the Arctic tundra.",tags$br(),
       "Above, the spatial distribution of studies is shown. These can be filtered  by a range of variables.",tags$br(),
       "The full data set can be downloaded below.",tags$br(),
       "The systematic review is published as Barbero-Palacios et al. 2024 Environmental Evidence", tags$a(href="https://www.youtube.com/watch?v=dQw4w9WgXcQ","here"), " [PLACEHOLDER]" ,tags$br(),
       "The systematic review protocol is published as Barrio et al. 2022 Environmental Evidence and can be downloaded", tags$a(href="https://environmentalevidencejournal.biomedcentral.com/articles/10.1186/s13750-022-00257-z", "here"),tags$br(),
       "Scripts and data are also available on the project ",tags$a(href="https://github.com/JamesDMSpeed/ArcticHerbivoreDiversitySystematicRev","GitHub repository"),tags$br(),
       "For further information contact", tags$a(href="mailto:James.Speed@ntnu.no", "James Speed"), "or",tags$a(href="mailto:isabel@lbhi.is",'Isabel Barrio'),tags$br(),
       "",tags$br(),
      ),
    uiOutput('activeFilters'), br(),
    
    
    
 
    # Download ----------------------------------------------------------------
    
    
    
    h4("Press the download button to download a speadsheet copy of the raw data:"),
    downloadButton("downloadData", "DownloadData")#,
    #h4("Press the download button to download the description of variables and coding:"),
    #downloadButton("downloadAppendix", "DownloadAppendix")
    
    
  ) # body
  
  
  
) # page





# SERVER ------------------------------------------------------------------



server <- function(input, output, session){ 
  
  # LOOKUP TABLE -------------------------------------------
  
  # oneFeatureTab <- reactive({
  #   Value <- t(MyTab[MyTab$article_ID==input$feature,])
  #   Variable <- rownames(Value)
  #   temp <- as.data.frame(cbind(Variable, Value))
  #   colnames(temp) <- c("Variable", "Value")
  #   temp
  # })
  # 
  # output$oneFeature <- renderTable(
  #   oneFeatureTab(), rownames = FALSE
  # )
  
  
  # FILTERED DATASET -------------------------------------
  
  datR <- reactive({
       MyTab[
     # MyTab$incl == TRUE &
        dplyr::between(MyTab$year, input$year[1], input$year[2]) &
        MyTab$country %in% input$country &
        MyTab$study_design %in% input$study_design &
        MyTab$diversity_contrast %in% input$diversity_contrast &
        MyTab$measured_response_variable_new %in% input$measured_response_variable_new #&
       # MyTab$change.long_f %in% input$change.long_f #&
      #  MyTab$experimental_design %in% input$expdesign
      ,]
    
  })
  
  
  # Remaining datapoints ----------------------------------------------------
  output$remaining <- renderText(nrow(datR()))
  
  
  
  # THE MAP ####
  output$theMap <- renderLeaflet({
    
    ifelse(input$filteron == TRUE, dat <- datR(), dat <- MyTab)
 
    dat2<-st_as_sf(dat,coords=c("coordinates_E","coordinates_N"),crs=CRS("+proj=longlat +datum=WGS84"))
    dat2 <- st_jitter(dat2, factor = 0.00001)
    
   pal <-  mapviewPalette("mapviewSpectralColors")
    
    m <- mapview::mapview(dat2["yi_smd"],
                          layer.name = "Effect size",
                          map.types = c("Esri.WorldImagery"),
                          cex = 5,#+abs(dat2$yi_smd),
                          alpha.regions = 0.8,
                          #na.color=grey(0.8,alpha=0.8),
                          #zcol=input$color,
                          popup = leafpop::popupTable(dat2, 
                                                      row.numbers = F, feature.id = F,
                                                      zcol = c("EffectSize",
                                                               "article_ID",
                                                               "FirstAuthor",
                                                               "year",
                                                               "journal",
                                                               "study_design"
                                                               )),
                          col.regions=pal(100),
                          #color=pal(100),
                          legend=T)
    
    m@map
    #addLegend(m@map,position='right',pal=colpoints,values=input$colour)
  })
  
  
  #Responsive Table ####
  output$responsiveTable <-
    renderDataTable({
      ifelse(input$filteron == TRUE, dat <- datR(), dat <- MyTab)
      DT::datatable(dat,
                    options = list(
                      scrollX = TRUE,
                      columnDefs = list(list(
                        targets = "_all",
                        render = JS(
                          "function(data, type, row, meta) {",
                          "return type === 'display' && data != null && data.length > 30 ?",
                          "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
                          "}")
                      ))),
                    class = "display")
    })


  # # Counts ####
  # 
  # 
  # 
  # output$univ <- renderPlot({
  #   ifelse(input$filteron == TRUE, dat <- datR(), dat <- MyTab)
  #   ggplot(data = dat, aes_string(x=input$uni))+
  #     geom_bar()+
  #     geom_bar(data = MyTab, alpha =0.3)+
  #     coord_flip()+
  #     theme_bw()+
  #     theme(text = element_text(size = 20))
  # })
  # 
  
  # pairwise ####
  # output$space <- renderPlot({
  #   
  #   
  #   ifelse(input$filteron == TRUE, dat <- datR(), dat <- MyTab)
  #   
  #   
  #   myG <- ggplot(data = dat, 
  #                 aes_string(x = input$var1, y = input$var2, size = input$var3))+    
  #     theme_bw()
  #   # remove size legend if size = NULL
  #   if(input$var3=="3"){
  #     myG <- myG + guides(size = FALSE)
  #   }
  #   
  #   if(input$var4!="NULL"){
  #     myG <- myG+    geom_point(aes_string(fill=input$var4),
  #                               shape=21,  stroke = 1, alpha=0.5, colour="black")+
  #       guides(fill = guide_legend(override.aes = list(size = 5)))+
  #       scale_color_viridis_d(option = "B", 
  #                             direction = -1,
  #                             aesthetics = c("fill"),
  #                             begin=0.4)} else{
  #                               myG <- myG+    geom_point(shape=21,  fill="blue", 
  #                                                         stroke = 1, alpha=0.5, colour="black")
  #                             }
  #   # add background points is var 1 and 2 are not GPW, Footprint ...                         
  #   #  if(!c(input$var1, input$var2) %in% EEvars2 ){
  #   myG <- myG +geom_point(data = range, aes_string(x = input$var1, y = input$var2),
  #                          alpha=input$alpha, size=2)  #}
  #   
  #   
  #   myG
  # })
  # 
  # TRENDS ####
  # output$trends <- renderPlot({
  #   ifelse(input$filteron == TRUE, dat <- datR(), dat <- MyTab)
  #   ggplot(data = dat, aes_string(x=input$cont))+
  #     geom_histogram()+
  #     geom_histogram(data = MyTab, alpha = 0.3)+
  #     theme_bw()+
  #     theme(text = element_text(size = 20))
  # })
  # 
  
  
  
  
   # DOWNLOAD ####
   output$downloadData <- downloadHandler(
     filename = function() {"data_with_effect_sizes.csv"},
     content = function(file) {
     write.csv( MyTab, file,  row.names = FALSE)}
  )
  # 
  # output$downloadAppendix<-downloadHandler(
  #   filename="Appendix 3.xlsx",
  #   content= function(file) {
  #     file.copy("www/Appendix 3_Coding template, including descriptions of coded variables.xlsx", file)}
  # )
  # 
  # output$downloadProtocol<-downloadHandler(
  #   filename='Soininen_et_al._2016_Protocol.pdf',
  #   content=function(file){
  #     file.copy('www/Soininen_2018_Protocol.pdf')
  #   }
  # )
  
 
 }
shinyApp(ui = ui, server = server)

