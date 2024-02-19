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
       "Above, the spatial distribution of evidence points is shown. These can be filtered  by a range of variables.",tags$br(),
       "The full data set can be downloaded below.",tags$br(),
       "The systematic review is published as Barbero-Palacios et al. 2024 Environmental Evidence", tags$a(href="https://www.youtube.com/watch?v=dQw4w9WgXcQ","here"), " [PLACEHOLDER]" ,tags$br(),
       "The systematic review protocol is published as Barrio et al. 2022 Environmental Evidence and can be downloaded", tags$a(href="https://environmentalevidencejournal.biomedcentral.com/articles/10.1186/s13750-022-00257-z", "here"),tags$br(),
       "For further information contact", tags$a(href="mailto:James.Speed@ntnu.no", "James Speed"), "or",tags$a(href="mailto:isabel@lbhi.is",'Isabel Barrio'),tags$br(),
       "",tags$br(),
      ),
    uiOutput('activeFilters'), br(),
    
    
    
    # # Lower Tabbox ------------------------------------------------------------
    # 
    # 
    # tabBox(width = NULL, id = 'additionals',
    #        
    #        
    #        # # . count cases -----------------------------------------------------------
    #        # 
    #        # 
    #        # tabPanel('Count cases',
    #        #          h5("This figure is responsive to the filters applied above. When filters are applied, the unfiltered data is shown as light gray bars.", tags$br(),
    #        #             "Download Appendix below for full description of variables and coding"),
    #        #          selectInput(inputId = "uni", 
    #        #                      label = "Select variable",
    #        #                      choices = varA,
    #        #                      selected = c("country"),
    #        #                      selectize = TRUE),
    #        #          plotOutput('univ')),  
    #        # 
    #        # 
    #        # # . trends ----------------------------------------------------------------
    #        # 
    #        # 
    #        # tabPanel('Univariate continuous variables',
    #        #          h5("This figure is responsive to the filters applied above", tags$br(),
    #        #             "Download Appendix below for full description of variables and coding"),
    #        #          
    #        #          pickerInput(
    #        #            inputId = "cont",
    #        #            label = "Select variable",
    #        #            choices = c("year",
    #        #                        "coordinates_N",
    #        #                        "coordinates_E",
    #        #                        "first_year_of_study",
    #        #                        "last_year_of_study",
    #        #                        EEvars)),
    #        #          plotOutput('trends')),
    #        # 
    #        # 
    #        # 
    #        # # . pairwise plots --------------------------------------------------------
    #        # 
    #        # 
    #        # tabPanel('Pairwise Plots',
    #        #          fluidRow(h5("The coloured circles are the evidence points 
    #        #        identified in the systematic review. The grey 
    #        #        background points show all possible values for each 
    #        #        variable inside the study region", tags$br(),
    #        #                      "Download Appendix below for full description of variables and coding")),
    #        #          fluidRow(
    #        #            column(width = 2,
    #        #                   pickerInput(
    #        #                     inputId = "var1",
    #        #                     label = "X variable", 
    #        #                     choices = EEvars,
    #        #                     selected = "Annual_Mean_Temperature",
    #        #                     options = list(size=10)
    #        #                   )),
    #        #            column(width = 2,
    #        #                   pickerInput(
    #        #                     inputId = "var2",
    #        #                     label = "Y variable", 
    #        #                     choices = EEvars,#[!EEvars %in% EEvars2],  # avoid weird error by removing 4 variables that don't have background points
    #        #                     selected = "Annual_Precipitation",
    #        #                     options = list(size=10)
    #        #                   )),
    #        #            column(width = 2,
    #        #                   pickerInput(
    #        #                     inputId = "var3",
    #        #                     label = "Size variable", 
    #        #                     choices = c("NULL" = 3, EEvars),
    #        #                     selected = "ArcticHerbivore_Species.richness",
    #        #                     options = list(size=10)
    #        #                   )),
    #        #            column(width=2,
    #        #                   pickerInput(
    #        #                     inputId = "var4",
    #        #                     label = "Colour",
    #        #                     options = list(size=10),
    #        #                     choices = c("NULL", varA) #"extent_of_spatial_scale", "study_design", "experimental_design")
    #        #                   )),
    #        #            column(width=3,
    #        #                   sliderInput('alpha', "Transparency of background points",
    #        #                               min=0, max = 1, step = 0.1, value = 0.2)
    #        #                   #switchInput(
    #        #                   #  inputId = "var5",
    #        #                   #  label = "Background points", 
    #        #                   #  value = TRUE,
    #        #                   #  labelWidth = "80px")
    #        #            )),
    #        #          fluidRow(     
    #        #            plotOutput('space')),
    #        #          fluidRow(
    #        #            h5("For information on the variables and coding please see Appendix 3" ,tags$br(),
    #        #               "For information about the climatic variables, 
    #        #   go to: https://www.worldclim.org/data/bioclim.html"))),
    #        # 
    #        # 
    #        # 
    #        # . lookup table ----------------------------------------------------------
    #        
    #        tabPanel('Lookup table',
    #                 selectInput('feature', 'Type the evidence point ID of the record you want to show. Get the unique ID by first clicking on the point on the map.', 
    #                             choices = MyTab$evidence_point_ID, multiple=FALSE, selectize=TRUE),
    #                 tableOutput('oneFeature')         
    #        ),
    #        
    #        
    #        # . responsive table ------------------------------------------------------
    #        
    #        
    #        tabPanel('Responsive table',
    #                 h5("This figure is responsive to the filters applied above"),
    #                 
    #                 DTOutput('responsiveTable'))
    #        
    # ),  # End panel box
    # 
    # 
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
  
#   # Active filters --------------------------------------------------------------
#   output$years <- renderPrint({
#     t <- input$year
#     #t <- strsplit(t, "\\s+")
#     #t <- noquote(unlist(t))
#     
#     t
#   })
#   
#   output$countries <- renderPrint({
#     t <- input$country
#     if(length(t) == length(unique(MyTab$country))) {t <- "All"} else{
#       #t <- gsub(pattern="\\[",replacement = "",t)
#       #t <- stringr::str_replace(t, " \\[.*\\]", "")
#       t <- strsplit(t, "\\s+")
#       t <- noquote(unlist(t))
#     }
#     t
#   })
#   
#   
#   output$sp <- renderPrint({
#     t <- input$species
#     if(length(t) == length(species4)) {t <- "All"} else{
#       t <- strsplit(t, "\\s+")
#       t <- noquote(unlist(t))
#     }
#     
#     t
#   })
#   
#   output$la <- renderPrint({
#     t <- input$language
#     if(length(t) == length(unique(MyTab$language))) {t <- "All"} else{
#       t <- strsplit(t, "\\s+")
#       t <- noquote(unlist(t))
#     }
#     t
#   })
#   
#   output$la <- renderPrint({
#     t <- input$language
#     if(length(t) == length(unique(MyTab$language))) {t <- "All"} else{
#       t <- strsplit(t, "\\s+")
#       t <- noquote(unlist(t))
#     }
#     t
#   })
#   
#   output$he <- renderPrint({
#     t <- input$herbivore
#     if(length(t) == length(unique(MyTab$herbivore_type))) {t <- "All"} else{
#       t <- strsplit(t, "\\s+")
#       t <- noquote(unlist(t))
#     }
#     t
#   })
#   
#   output$ty <- renderPrint({
#     t <- input$studydesign
#     if(length(t) == length(unique(MyTab$study_design))) {t <- "All"} else{
#       t <- strsplit(t, "\\s+")
#       t <- noquote(unlist(t))
#     }
#     t
#   })
#   output$me <- renderPrint({
#     t <- input$studymethod
#     if(length(t) == length(unique(MyTab$study_method))) {t <- "All"} else{
#       t <- strsplit(t, "\\s+")
#       t <- noquote(unlist(t))
#     }
#     t
#   })
#   output$ex <- renderPrint({
#     t <- input$expdesign
#     if(length(t) == length(unique(MyTab$experimental_design))) {t <- "All"} else{
#       t <- strsplit(t, "\\s+")
#       t <- noquote(unlist(t))
#     }
#     t
#   })
#   
#   output$activeFilters <- renderUI({
#     dropdownButton(
#       h3("Active Filters"),
#       h4("Years: "), textOutput('years'),
#       h4("Countries: "), textOutput('countries'),
#       h4("Species: "), textOutput('sp'),
#       h4("Language: "), textOutput('la'),
#       h4("Herbivore types: "), textOutput('he'),
#       h4("Study design types: "), textOutput('ty'),
#       h4("Study method: "), textOutput('me'),
#       h4("Experimental design: "), textOutput('ex'),
#       circle = TRUE, status = "info",
#       icon = icon("list"), width = "300px",
#       tooltip = tooltipOptions(title = "Click to see which filters are active")
#     )
#   })
#   
 }
shinyApp(ui = ui, server = server)

