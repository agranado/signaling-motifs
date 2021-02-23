
# Goal: To make an interactive subplot that updates when hovering on a different subplot
# This app reads file in the transcriptome/ folder

library(shiny)
library(plotly)

library(shiny)
library(tibble)
library(pheatmap)


ui <- fluidPage(
  #theme = shinytheme("cerulean"),
  #sliderbar for heatmap
  titlePanel("Global map of signaling motifs") ,

  # conditionalPanel(condition = "$('li.active a').first().html()==='Heatmap'",
  #                  h2 ("Click the button to produce heatmap"),
  #                  actionButton('getHmap', 'get heatmap')),


  mainPanel (
      tabsetPanel (
            tabPanel('Heatmap',
                     fluidRow(column(8, offset = 1,
                                     h2("Heatmap"),
                                     plotlyOutput("theheatmap", width = "1200px", height = "400px"),
                                    )
                            ), # fluidRow

                     fluidRow(column(8, offset = 1,

                                      h2('UMAP'),
                                      plotlyOutput("plot", width = "900px", height = "600px"),
                                    )


                             ) #fluidRow
                    ),# tabPanel
            tabPanel("UMAP", # Second tab panel

                        sidebarPanel(
                              tags$h3("Select motifs to show:"),
                              #checkboxGroupInput('dlist', 'Display List:', displayList, selected = displayList[1:2])
                              uiOutput('checkbox') #dynamically created in the server section
                              ), # sidebarPanel
                        mainPanel(
                              h1("Signaling motifs distribution"),

                              h4("Global UMAP coordinates"),
                              #verbatimTextOutput("txtout"),

                              # conditionalPanel(condition = "input.dlist.indexOf('1') > -1",
                              #                 p("List 1 selected")
                              #               ), #conditionalPanel

                              # UMAP
                              plotlyOutput('motif_umap', width = "900px", height = "600px"),
                              # Barplot
                              plotlyOutput('barplot', width = "500px", height = "200px"),

                              # Table with cell type information
                              dataTableOutput('table')
                        ) # mainPanel
                    )
      ) #tabsetPanel
  ) # mainPanel

)

server <- function(input, output){

  # load the data on start
  # in principle we could load the data as reactive or with user input i.e. differect pathways ?
  raw_counts = read.csv('app/filtered_counts.csv', header = T, row.names = 1)
  ann_counts = read.csv('app/annotated_counts.csv', header = T, row.names = 1)

  ann_counts$global_cluster <- as.character(ann_counts$global_cluster) # for compatibility
  row.names(ann_counts) <- ann_counts$global_cluster
  # 3D umap labeled by motif (includes all cells)
  motifs <- read.csv(file= 'app/global_transcriptome_motifLabeled.csv', sep =",", header= T)

  # Merge with motif labels:
    # Note: this could be done before exporting the .csv
  motifs_fil <- motifs[motifs$motif_label>0, ] # only cell types with motif annotation
  motifs_fil$global_cluster <- as.character(motifs_fil$global_cluster)
  # only profiles with label from recursive clsutering
  ann_counts %>% dplyr::filter(global_cluster %in% motifs_fil$global_cluster) -> ann_counts
  raw_counts %>% dplyr::filter(global_cluster %in% motifs_fil$global_cluster) -> raw_counts


  x_mat<-as.matrix(raw_counts[-1])
  row.names(x_mat) <- raw_counts$global_cluster # for heatmap


  ann_counts %>% left_join(motifs_fil %>%
            dplyr::select(global_cluster, motif_label), by='global_cluster') -> motifs_ann
  motifs_ann$motif_label <- as.character(motifs_ann$motif_label) # for colors to show up
  row.names(motifs_ann)<- motifs_ann$global_cluster

  # Barplot
  tidy_pathway = gather(raw_counts, "gene", "Expr", -c(global_cluster))

  # Colors
  colors$motif_label = makeQualitativePal(length(motifs_ann$motif_label %>% unique() ))
  names(colors$motif_label) <- motifs_ann$motif_label %>% unique %>% sort

  ###############
  #### FUNCTIONS
  data <- reactive({
    mouse_event <- event_data('plotly_click', source = 'main_umap', priority= 'event')
    if(length(mouse_event)){
      ann_counts %>% dplyr::filter(global_cluster == mouse_event$customdata) -> this_profile
      return(this_profile)
    }
  })

  plotdata <- eventReactive(input$getHmap, {
    #d = data()
    data <- as.matrix(data()[-1] )
    row.names(data) <- data()$global_cluster

    #dataplot[is.na(dataplot)] <- 0
    data
  })

  ### MAKE PLOTS

  ## UI
  # this will render automatically
  output$checkbox <- renderUI({
    labels = motifs$motif_label %>% unique()
    choice = labels[labels>0] %>% sort()
    checkboxGroupInput("checkbox","Select motif", choices = choice, selected = choice[1:10])

  }) #output$checkbox will render a checkbox based on the dataset

  output$theheatmap = renderPlotly({
    # 8. Final heatmap and save dendrogram
   # pheatmap(t(x_mat), # with reactive it would be t(plotdata() )
   #                  clustering_distance_cols = dist(x_mat %>% as.matrix),
   #                  annotation_col = motifs_ann %>% dplyr::select( Tissue,dataset, motif_label),
   #                  clustering_method = 'ward.D2',
   #                  annotation_colors = colors,
   #                  show_colnames = F, cluster_rows = F, fontsize = 12,
   #                  cutree_cols = 10  , col = blues_pal(100))


      heatmaply(
        t(x_mat),
        seriate = 'OLO',
        col_side_colors = motifs_ann %>% dplyr::select(motif_label),
        Rowv =F,
        colors  = blues_pal(100)
      )

  })


  # PAGE 1
  # UMAP
  output$plot <- renderPlotly({

    # Select meta columns from main data.frame
    meta_master = motifs %>% dplyr::select(cell_id, Tissue, cell_ontology_class, age, dataset)

    # this function is in the workspace
    # ideally we need to load it in the server
    colors<-makeColorsAll(meta_master, list() )

    fig <- plot_ly(motifs, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
          marker = list(size = 3), color = ~Cell_class,
          colors = colors$Cell_class,
          text=~paste("Tissue:",Tissue,"<br>Age:",age,"<br>dataset:",
                        dataset,"<br>Cell type:", cell_ontology_class,
                        "<br>Motif", motif_label),
          source = 'tissue_umap', customdata = ~global_cluster)

      fig %>% add_markers()
  })

  # SECOND PAGE of the app
  output$txtout2 <-renderText({
                        paste(input$dlist," " ,sep = "")
  })

  output$motif_umap <- renderPlotly({

    local_colors = colors
    # Select meta columns from main data.frame
    meta_master = motifs %>% dplyr::select(cell_id, Tissue, cell_ontology_class, age, dataset)

    motifs$motif_id = motifs$motif_label
    # from the inputs in the check box, color the ones selected
    motifs[which(!motifs$motif_label %in% input$checkbox), ]$motif_id = 0
    motifs$motif_id <- as.character(motifs$motif_id)
    # this function is in the workspace
    # ideally we need to load it in the server
    #colors<-makeColorsAll(meta_master, list() )
    # Make a new color for the 'unselected' class
    local_colors$motif_label  = c('0' = "#BEBEBE", local_colors$motif_label)
    this_colors <- local_colors$motif_label[motifs$motif_id %>% unique() %>% sort() %>% as.character() ]

    fig <- plot_ly(motifs, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
          marker = list(size = 3), color = ~motif_id,
          colors = this_colors,
          text=~paste("Tissue:",Tissue," Age:",age,"<br>dataset:",dataset,
                      "<br>Cell type:", cell_ontology_class,
                      "<br>Motif:", motif_label),
          source = 'main_umap', customdata = ~global_cluster)

      fig %>% add_markers()
  })


  output$barplot <- renderPlotly({
    # generating graph dynamically on hover:
    # 1. Detect mouse event to know where the mouse pointer is
    mouse_event <- event_data('plotly_click', source = 'main_umap', priority= 'event')

    print(mouse_event)
    # mouse event contains the x, y values of plot1
    # which in thise case are x = years, y = J.D
    # the x attribute of mouse event in the third column:
    #year <- mouse_event[3]
    # Filter the global cluster selected by the mouse event
    if(length(mouse_event)){
      print(motifs[motifs$global_cluster == mouse_event$customdata, ])
      # based on this event, create a subset of the dataset:
      #monthly_subset <- monthly[monthly$Year == year$x,]
      tidy_pathway %>% dplyr::filter(global_cluster == mouse_event$customdata) -> this_profile
      ann_counts %>%  dplyr::filter(global_cluster == mouse_event$customdata) -> this_meta
      fig<- plot_ly(this_profile, x = ~gene, y = ~Expr, type   = 'bar')
      fig %>% layout(title = paste(this_meta$age, this_meta$Tissue,":",this_meta$cell_ontology_class, sep = " "),
                yaxis = list(title = 'Norm Expression',
                          range=c(0,max(this_profile$Expr))
                        ),
                xaxis = list(
                          title ='Gene'
                        )
                    )

    }# mouse_event

  }) # barplot

  output$table <- renderDataTable({
      DT::datatable(data())
  })

}
# Run the app
shinyApp(ui, server)
