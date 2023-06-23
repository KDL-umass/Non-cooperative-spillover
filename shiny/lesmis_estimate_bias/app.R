library(shiny)

source('../lesmis.R')
source('../../experiments/adversary-experiment.R')


test.single.config <- function(idx, configs, trials=1, all=TRUE) { 
  cat("Running", idx, "\n")
  print(configs[idx,])
  
  results <- data.frame(index=numeric(), size.of.dom=logical(), method=character(), 
                        pt.uncovered=numeric(), adversary.influence=numeric(), ATE.true=numeric(), 
                        variable=numeric(), value=numeric(), pt.covered=numeric(), n=numeric(), 
                        graph.type=character(), power=numeric(), degree=numeric(), p=numeric(), 
                        mu=numeric(), ncoms=numeric(), maxc=numeric(), minc=numeric(), 
                        lambda_0=numeric(), lambda_1=numeric(), lambda_2=numeric(), stringsAsFactors=FALSE)
  
  graph.params <- build.graph.params(configs, idx)
  adversary.params <- list()
  adversary.params$model <- reduction.adv.model
  adversary.params$all <- all
  outcome.params <- build.outcome.params(configs[idx,"lambda_0"], configs[idx,"lambda_1"], configs[idx,"lambda_2"], configs[idx,"sd.noise"])
  clustering <- "infomap"
  
  for(i in 1:trials) {
    graph.params$ind <- i
    
    cat("trial", i, "\n")
    bias.behavior.ATE <- adversary.experiment(graph.params, clustering, adversary.params, outcome.params)
    bias.behavior.ATE$adversary.influence <- as.numeric(bias.behavior.ATE$adversary.influence)
    bias.behavior.ATE$gui.beta <- as.numeric(bias.behavior.ATE$gui.beta)
    bias.behavior.ATE$gui.gamma <- as.numeric(bias.behavior.ATE$gui.gamma)
    
    bias.behavior.ATE <- add.graph.params(bias.behavior.ATE, graph.params)
    bias.behavior.ATE <- add.outcome.params(bias.behavior.ATE, outcome.params)
    bias.behavior.ATE$graph.id <- configs[idx,"graph.no"]
    bias.behavior.ATE$adv.bias <- bias.behavior.ATE$nonadv.ATE - bias.behavior.ATE$ATE.adv.gui
    
    results <- rbind(results, bias.behavior.ATE)
    #write.csv(results, paste0("results/adversary-results-", graph.params$graph.type, "-", idx, ".csv"))
  }
}


# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Bias estimates in Les Mis!"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "n_ncps",
                  label = "Number of NCPs:",
                  min = 0,
                  max = graph.properties$n,
                  value = 11)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotOutput(outputId = "exposureGraphPlot")
      
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  # Plot of Les Mis dialog graph ----
  # NCPs marked in orange
  # unexposed nodes marked in light blue 
  # NCP exposure 
  # with requested number of NCPs
  
  # This expression that generates a graph plot is wrapped in a call
  # to renderPlot to indicate that:
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$n_ncps) change
  # 2. Its output type is a plot
  
  output$exposureGraphPlot <- renderPlot({
    
    #adversaries.deg <- unlist(list(determine.adversaries(graph.properties, ncp.params)))
    #adversaries <- matrix(0,1,graph.properties$n)
    #adversaries[which(adversaries.deg==1)] <- 1
    
    #random NCP selection
    n_adversaries <- input$n_ncps
    random.adversaries <- sample(1:graph.properties$n, n_adversaries, replace=FALSE)
    adversaries <- matrix(0,1,graph.properties$n)
    adversaries[random.adversaries] <- 1
    
    treatment <- treatment.assignment(graph.properties$g, clusters)
    treatment.assignments <- treatment[clusters]
    
    uncovered.vertices <- 1 - adversaries %*% graph.properties$adj - adversaries
    pt.uncovered <- sum(uncovered.vertices == 1)/graph.properties$n
    #prepare.for.plots(g, adversaries, ncp.params$ncp.exposure, treatment.assignments, labels=TRUE)
    exposure.params <- exposure.probs(ncp.params, graph.properties, treatment.assignments, adversaries)
    print(exposure.params)
    #V(g)$color <- adversaries.deg
    
    bdr <- rep("black", length(V(g)))
    if(length(treatment.assignments) > 2) bdr <- ifelse(treatment.assignments, "green", "black")
    
    g$palette <- grey.colors(100)
    # ncp.params$ncp.exposure.neighbors
    V(g)$color <- 100 - exposure.params$ncp.exposure.neighbors * 100 # grey.colors runs dark to light
    V(g)$color[which(exposure.params$ncp.exposure.neighbors == 0)] <- "lightblue"
    V(g)$color[which(adversaries > 0)] <- "orange"
    
    
    lo <- layout_with_kk(g) # create a layout
    lo <- norm_coords(lo, ymin=-1, ymax=1, xmin=-1, xmax=1)
    #plot(g, layout=lo*2, vertex.color=V(g)$color, vertex.frame.color=bdr)
    plot(g, layout=lo*2, vertex.color=V(g)$color)
    
  })
  
}

shinyApp(ui = ui, server = server)