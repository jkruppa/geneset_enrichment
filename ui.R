library(shiny)


shinyUI(fluidPage(
  titlePanel("Geneset Enrichment Analysis"),
  sidebarLayout(
      ## sidebar panel
      sidebarPanel(
          ## Controls
          h3("Controls"),
          sliderInput("degs", 
                      "Number of differentially expression genes (of 10.000)", 
                      min = 1,
                      max = 1000,
                      value = 50),
          sliderInput("pathgenes", 
                      "Number of genes related to pathway", 
                      min = 2,
                      max = 500,
                      value = 10),
          sliderInput("perpath", 
                      "Percentage of pathway genes that are significant", 
                      min = 0,
                      max = 100,
                      value = 20),
          helpText("ADD HELP TEXT"),
          helpText("ADD HELP TEXT"),
          hr(),
          ## Dependencies
          h3("Dependencies"),
          p("The following R packages are needed for the shiny app."),
          code('install.packages("shiny")'),
          br(),
          code('install.packages("mvtnorm")'), 
          br(),
          code('install.packages("scatterplot3d")'), 
          hr()
          ## References
          ## h3("References"),
          ## p("Cui X and Churchill GA (2003). 'Statistical tests for differential expression in cDNA microarray experiments'. Genome Biology, 4(4), 210.")
      ),
      ## main panel
      ## Show a plot of the generated distribution
      mainPanel(
          plotOutput("plots"),
          plotOutput("plots2")
      ),
  )
))
