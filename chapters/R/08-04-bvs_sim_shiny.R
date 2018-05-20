# Chapter 6
# Bayesian Variable Selection shiny app (analysis)
library(tidyverse)
library(reshape2)
chapter.no <- "06"

# Load all results
load("data/bvs-res-gam90")  # res.gam90
load("data/bvs-res-gam75")  # res.gam75
load("data/bvs-res-gam50")  # res.gam50
load("data/bvs-res-gam25")  # res.gam25
load("data/bvs-res-gam10")  # res.gam10

# Re-arrange data
res.snr90 <- res.gam90
res.snr90[[1]] <- res.gam90[[1]]  # SNR = 0.90, gamma = 0.90
res.snr90[[2]] <- res.gam75[[1]]  # SNR = 0.90, gamma = 0.75
res.snr90[[3]] <- res.gam50[[1]]  # SNR = 0.90, gamma = 0.50
res.snr90[[4]] <- res.gam25[[1]]  # SNR = 0.90, gamma = 0.25
res.snr90[[5]] <- res.gam10[[1]]  # SNR = 0.90, gamma = 0.10

res.snr75 <- res.gam75
res.snr75[[1]] <- res.gam90[[2]]  # SNR = 0.75, gamma = 0.90
res.snr75[[2]] <- res.gam75[[2]]  # SNR = 0.75, gamma = 0.75
res.snr75[[3]] <- res.gam50[[2]]  # SNR = 0.75, gamma = 0.50
res.snr75[[4]] <- res.gam25[[2]]  # SNR = 0.75, gamma = 0.25
res.snr75[[5]] <- res.gam10[[2]]  # SNR = 0.75, gamma = 0.10

res.snr50 <- res.gam50
res.snr50[[1]] <- res.gam90[[3]]  # SNR = 0.50, gamma = 0.90
res.snr50[[2]] <- res.gam75[[3]]  # SNR = 0.50, gamma = 0.75
res.snr50[[3]] <- res.gam50[[3]]  # SNR = 0.50, gamma = 0.50
res.snr50[[4]] <- res.gam25[[3]]  # SNR = 0.50, gamma = 0.25
res.snr50[[5]] <- res.gam10[[3]]  # SNR = 0.50, gamma = 0.10

res.snr25 <- res.gam25
res.snr25[[1]] <- res.gam90[[4]]  # SNR = 0.25, gamma = 0.90
res.snr25[[2]] <- res.gam75[[4]]  # SNR = 0.25, gamma = 0.75
res.snr25[[3]] <- res.gam50[[4]]  # SNR = 0.25, gamma = 0.50
res.snr25[[4]] <- res.gam25[[4]]  # SNR = 0.25, gamma = 0.25
res.snr25[[5]] <- res.gam10[[4]]  # SNR = 0.25, gamma = 0.10

res.snr10 <- res.gam10
res.snr10[[1]] <- res.gam90[[5]]  # SNR = 0.10, gamma = 0.90
res.snr10[[2]] <- res.gam75[[5]]  # SNR = 0.10, gamma = 0.75
res.snr10[[3]] <- res.gam50[[4]]  # SNR = 0.10, gamma = 0.50
res.snr10[[4]] <- res.gam25[[5]]  # SNR = 0.10, gamma = 0.25
res.snr10[[5]] <- res.gam10[[5]]  # SNR = 0.10, gamma = 0.10

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
    data.frame(res[[1]][[stage]], gamma = "0.90"),
    data.frame(res[[2]][[stage]], gamma = "0.75"),
    data.frame(res[[3]][[stage]], gamma = "0.50"),
    data.frame(res[[4]][[stage]], gamma = "0.25"),
    data.frame(res[[5]][[stage]], gamma = "0.10")
  )
}

plot_res2 <- function(stage = 1, type = "false") {

  tmp <- function(res) {
    res <- matrix(rapply(res, function(x) apply(x, 2, mean, na.rm = TRUE)),
                  ncol = 4, byrow = TRUE)
    if (stage == 1) res <- res[c(1, 3, 5, 7, 9), ]
    if (stage == 2) res <- res[c(2, 4, 6, 8, 10), ]
    data.frame(res, gamma = c(0.90, 0.75, 0.50, 0.25, 0.10))
  }
  df <- data.frame(rbind(
    cbind(tmp(res.snr90), SNR = "90%"),
    cbind(tmp(res.snr75), SNR = "75%"),
    cbind(tmp(res.snr50), SNR = "50%"),
    cbind(tmp(res.snr25), SNR = "25%"),
    cbind(tmp(res.snr10), SNR = "10%")
  ))
  df <- df[, c(2, 1, 3, 4, 5, 6)]
  colnames(df)[1:4] <- c("false.exc", "false.inc", "false", "brier")
  df$SNR <- factor(df$SNR, levels = c("90%", "75%", "50%", "25%", "10%"))
  df <- reshape2::melt(df, id = c("gamma", "SNR"))
  levels(df$variable)[1:2] <- c("False exclusion", "False inclusion")

  ggplot(subset(df, df$var == "False inclusion" | df$var == "False exclusion"),
         aes(x = gamma, group = SNR)) +
    geom_line(aes(y = value,  col = SNR)) +
    facet_grid(. ~ variable) +
    labs(y = "Count", x = expression(paste("Hyperprior setting for ", pi[j]))) +
    theme_bw()
}
plot_res2(2)

ggsave("figure/06-sens_analysis.pdf", plot_res2(2), "pdf",
       width = 7, height = 3.2)
move_fig_to_chapter()


plot_res <- function(res, stage = 1, type = "false") {
  if (type == "False choices") type <- "false"
  if (type == "False inclusions") type <- "false.inc"
  if (type == "False exclusions") type <- "false.exc"

  if (type == "false") {
    p <- ggplot(create_bvs_plot_df(res, stage),
                aes(false, ..density.., col = gamma)) +
      ggtitle(paste0("I-prior (false choices) stage = ", stage))
  } else if (type == "false.inc") {
    p <- ggplot(create_bvs_plot_df(res, stage),
                aes(false.inc, ..density.., col = gamma)) +
      ggtitle(paste0("I-prior (false inclusions) stage = ", stage))
  } else if (type == "false.exc") {
    p <- ggplot(create_bvs_plot_df(res, stage),
                aes(false.exc, ..density.., col = gamma)) +
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
        tabPanel("SNR = 0.90",
          verticalLayout(
            plotOutput("gam90.1"),
            plotOutput("gam90.2")
          )
        ),
        tabPanel("SNR = 0.75",
                 verticalLayout(
                   plotOutput("gam75.1"),
                   plotOutput("gam75.2")
                 )
        ),
        tabPanel("SNR = 0.50",
                 verticalLayout(
                   plotOutput("gam50.1"),
                   plotOutput("gam50.2")
                 )
        ),
        tabPanel("SNR = 0.25",
                 verticalLayout(
                   plotOutput("gam25.1"),
                   plotOutput("gam25.2")
                 )
        ),
        tabPanel("SNR = 0.10",
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
    plot_res(res.snr90, stage = 1, type = input$type)
  })

  output$gam90.2 <- renderPlot({
    plot_res(res.snr90, stage = 2, type = input$type)
  })

  output$gam75.1 <- renderPlot({
    plot_res(res.snr75, stage = 1, type = input$type)
  })

  output$gam75.2 <- renderPlot({
    plot_res(res.snr75, stage = 2, type = input$type)
  })

  output$gam50.1 <- renderPlot({
    plot_res(res.snr50, stage = 1, type = input$type)
  })

  output$gam50.2 <- renderPlot({
    plot_res(res.snr50, stage = 2, type = input$type)
  })

  output$gam25.1 <- renderPlot({
    plot_res(res.snr25, stage = 1, type = input$type)
  })

  output$gam25.2 <- renderPlot({
    plot_res(res.snr25, stage = 2, type = input$type)
  })

  output$gam10.1 <- renderPlot({
    plot_res(res.snr10, stage = 1, type = input$type)
  })

  output$gam10.2 <- renderPlot({
    plot_res(res.snr10, stage = 2, type = input$type)
  })

}

# Run the application
shinyApp(ui = ui, server = server)
