# Chapter 8
# Bayesian Variable Selection shiny app
library(tidyverse)
library(reshape2)

# Load all results
load("data/gamma90/bvs-res-gam90")  # res.gam90
load("data/gamma75/bvs-res-gam75")  # res.gam75
load("data/gamma50/bvs-res-gam50")  # res.gam50
load("data/gamma25/bvs-res-gam25")  # res.gam25
# load("data/gamma10/bvs-res-gam10")  # res.gam10

## ---- bvs.sims.res ----
count_false_choices <- function(x) {
  x.0to2 <- mean(x <= 2)
  x.6 <- mean(x >= 6)
  x.3to5 <- 1 - x.0to2 - x.6
  ress <- c(x.0to2, x.3to5, x.6)
  names(ress) <- c("0-2", "3-5", ">5")
  ress
}

se_false_choices <- function(x) {
  n <- length(x)
  x.0to2 <- (sd(x <= 2) * (n - 1) / n) / sqrt(n)
  x.6 <- (sd(x >= 6) * (n - 1) / n) / sqrt(n)
  x.3to5 <- (x > 2) & (x < 6)
  x.3to5 <- (sd(x.3to5) * (n - 1) / n) / sqrt(n)
  ress <- c(x.0to2, x.3to5, x.6)
  names(ress) <- c("0-2", "3-5", ">5")
  ress
}

bvs_res <- function(x = res, type = mean) {
  res <- lapply(x, function(y) t(sapply(y, function(z) apply(z, 2, type))))

  if (ncol(res[[1]]) > 3) {
    res <- lapply(res, function(y) {
      y <- y[, 7:9];
      colnames(y) <- c("0-2", "3-5", ">5");
      y
    })
  }

  res
}

create_bvs_plot_df <- function(res, stage = 1) {
  rbind(
    data.frame(res[[1]][[stage]], SNR = "90%"),
    data.frame(res[[2]][[stage]], SNR = "75%"),
    data.frame(res[[3]][[stage]], SNR = "50%"),
    data.frame(res[[4]][[stage]], SNR = "25%"),
    data.frame(res[[5]][[stage]], SNR = "10%")
  )
}

theme_set(theme_bw())

plot_res <- function(res, stage = 1, type = "false") {
  if (type == "False choices") type <- "false"
  if (type == "False inclusions") type <- "false.inc"
  if (type == "False exclusions") type <- "false.exc"

  if (type == "false") {
    p <- ggplot(create_bvs_plot_df(res, stage),
                aes(false, ..density.., col = SNR)) +
      ggtitle(paste0("I-prior (false choices) stage = ", stage))
  } else if (type == "false.inc") {
    p <- ggplot(create_bvs_plot_df(res, stage),
                aes(false.inc, ..density.., col = SNR)) +
      ggtitle(paste0("I-prior (false inclusions) stage = ", stage))
  } else if (type == "false.exc") {
    p <- ggplot(create_bvs_plot_df(res, stage),
                aes(false.exc, ..density.., col = SNR)) +
      ggtitle(paste0("I-prior (false exclusions) stage = ", stage))
  }

  p +
    geom_freqpoly(binwidth = 1, size = 0.8) +
    scale_y_continuous(limits = c(0, 0.85)) +
    scale_x_continuous(limits = c(0, 90)) +
    labs(x = "Number of false choices", y = "Proportion")
}


# Define UI for application that draws a histogram
ui <- fluidPage(

  # Application title
  titlePanel("I-prior BVS results"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      radioButtons("type",
                  "Plot to show",
                  c("False choices", "False inclusions", "False exclusions"))
    ),

    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("prior = 0.90",
          verticalLayout(
            plotOutput("gam90.1"),
            plotOutput("gam90.2")
          )
        ),
        tabPanel("prior = 0.75",
                 verticalLayout(
                   plotOutput("gam75.1"),
                   plotOutput("gam75.2")
                 )
        ),
        tabPanel("prior = 0.50",
                 verticalLayout(
                   plotOutput("gam50.1"),
                   plotOutput("gam50.2")
                 )
        ),
        tabPanel("prior = 0.25",
                 verticalLayout(
                   plotOutput("gam25.1"),
                   plotOutput("gam25.2")
                 )
        ),
        tabPanel("prior = 0.10",
                 verticalLayout(
                   plotOutput("gam10.1"),
                   plotOutput("gam10.2")
                 )
        )
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$gam90.1 <- renderPlot({
    plot_res(res.gam90, stage = 1, type = input$type)
  })

  output$gam90.2 <- renderPlot({
    plot_res(res.gam90, stage = 2, type = input$type)
  })

  output$gam75.1 <- renderPlot({
    plot_res(res.gam75, stage = 1, type = input$type)
  })

  output$gam75.2 <- renderPlot({
    plot_res(res.gam75, stage = 2, type = input$type)
  })

  output$gam50.1 <- renderPlot({
    plot_res(res.gam50, stage = 1, type = input$type)
  })

  output$gam50.2 <- renderPlot({
    plot_res(res.gam50, stage = 2, type = input$type)
  })

  output$gam25.1 <- renderPlot({
    plot_res(res.gam25, stage = 1, type = input$type)
  })

  output$gam25.2 <- renderPlot({
    plot_res(res.gam25, stage = 2, type = input$type)
  })

  output$gam10.1 <- renderPlot({
    plot_res(res.gam10, stage = 1, type = input$type)
  })

  output$gam10.2 <- renderPlot({
    plot_res(res.gam10, stage = 2, type = input$type)
  })

}

# Run the application
shinyApp(ui = ui, server = server)
