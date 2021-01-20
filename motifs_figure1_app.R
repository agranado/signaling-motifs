# Goal: To make an interactive subplot that updates when hovering on a different subplot

library(shiny)
library(plotly)


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
                               plotOutput("theheatmap", width = "1200px", height = "400px"),
                              )
                      ),

              fluidRow(column(8, offset = 1,

                h2('UMAP'),
                plotlyOutput("plot", width = "900px", height = "600px"),
              )


             )
    )
  )
)

)

server <- function(input, output){

  # load the data on start
  # in principle we could load the data as reactive or with user input i.e. differect pathways ?
  raw_counts = read.csv('app/filtered_counts.csv', header = T, row.names = 1)
  ann_counts = read.csv('app/annotated_counts.csv', header = T, row.names = 1)


  # 3D umap labeled by motif (includes all cells)
  motifs <- read.csv(file= 'app/global_transcriptome_motifLabeled.csv', sep =",", header= T)

  # we do it outside the reactive functions so it is executed only once
  x <- rownames_to_column(raw_counts, 'global_cluster')
  ann_counts <- rownames_to_column(ann_counts, 'global_cluster') # goes to first column
  x_mat<-as.matrix(x[-1])
  row.names(x_mat) <- x$global_cluster

  # Merge with motif labels:
    # Note: this could be done before exporting the .csv
  motifs_fil <- motifs[motifs$motif_label>0, ] # only cell types with motif annotation
  motifs_fil$global_cluster <- as.character(motifs_fil$global_cluster)
  ann_counts %>% left_join(motifs_fil %>% dplyr::select(global_cluster, motif_label), by='global_cluster') -> motifs_ann
  motifs_ann$motif_label <- as.character(motifs_ann$motif_label) # for colors to show up


  data <- reactive({
    #x <- head(mtcars)
    #x <- rownames_to_column(x, "Name")
    # puts the new column on the first column
    x <- rownames_to_column(raw_counts,'global_cluster')
    #x <- x %>% dplyr::select(-c(Tissue,age, dataset, cell_ontology_class) )
    return(x)
  })

  plotdata <- eventReactive(input$getHmap, {
    #d = data()
    data <- as.matrix(data()[-1] )
    row.names(data) <- data()$global_cluster

    #dataplot[is.na(dataplot)] <- 0
    data
  })


  output$theheatmap = renderPlot({
    # 8. Final heatmap and save dendrogram
   pheatmap(t(x_mat), # with reactive it would be t(plotdata() )
                    clustering_distance_cols = 'euclidean',
                    annotation_col = motifs_ann %>% dplyr::select( Tissue,dataset, motif_label),
                    clustering_method = 'ward.D2',
                    annotation_colors = colors,
                    show_colnames = F, cluster_rows = F, fontsize = 12,
                    cutree_cols = 10  , col = blues_pal(100))
  })


  # we can in principle detect a mouse_event here
  output$plot <- renderPlotly({


    # Select meta columns from main data.frame
    meta_master = motifs %>% dplyr::select(cell_id, Tissue, cell_ontology_class, age, dataset)

    # this function is in the workspace
    # ideally we need to load it in the server
    colors<-makeColorsAll(meta_master, list() )

    fig <- plot_ly(motifs, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
          marker = list(size = 3), color = ~Tissue,
          colors = colors$Tissue,
          text=~paste("Tissue:",Tissue,"<br>Age:",age,"<br>dataset:",dataset,"<br>Cell type:", cell_ontology_class),
                      shoverinfo = 'text',
          source = 'main_umap', customdata = ~global_cluster)

      fig %>% add_markers()
  })

}
# Run the app
shinyApp(ui, server)
