library(shiny)
library(bslib)
library(ggplot2)
library(purrr)
library(viridis)

# ============================================================================
# QURS (Quantile Unit Rayleigh-Sin) DISTRIBUTION FUNCTIONS
# ============================================================================

# Probability Density Function
dRURS <- function(y, mu, tau, log = FALSE) {
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(tau <= 0) | any(tau >= 1)) stop(paste("tau must be between 0 and 1", "\n", ""))
  if (any(y <= 0) | any(y >= 1)) stop(paste("y must be between 0 and 1", "\n", ""))
  alpha <- log((2 / pi) * asin(tau / (2 - tau))) / log(mu)
  fy <- (pi * alpha * y^(alpha - 1) * cos((pi * y^alpha) / 2)) /
    (sin((pi * y^alpha) / 2) + 1)^2
  if (log == FALSE) pdf <- fy else pdf <- log(fy)
  return(pdf)
}

# Cumulative Distribution Function
pRURS <- function(q, mu, tau, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(tau <= 0) | any(tau >= 1)) stop(paste("tau must be between 0 and 1", "\n", ""))
  if (any(q <= 0) | any(q >= 1)) stop(paste("q must be between 0 and 1", "\n", ""))
  alpha <- log((2 / pi) * asin(tau / (2 - tau))) / log(mu)
  Fy <- (2 * sin((pi / 2) * q^alpha)) / (1 + sin((pi / 2) * q^alpha))
  if (lower.tail == TRUE) cdf <- Fy else cdf <- 1 - Fy
  if (log.p == FALSE) cdf <- cdf else cdf <- log(cdf)
  return(cdf)
}

# Quantile Function
qRURS <- function(p, mu, tau) {
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(tau <= 0) | any(tau >= 1)) stop(paste("tau must be between 0 and 1", "\n", ""))
  if (any(p <= 0) | any(p >= 1)) stop(paste("p must be between 0 and 1", "\n", ""))
  alpha <- log((2 / pi) * asin(tau / (2 - tau))) / log(mu)
  Q <- (2 / pi * asin(p / (2 - p)))^(1 / alpha)
  return(Q)
}

# Random Generation
rRURS <- function(n, mu, tau) {
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(tau <= 0) | any(tau >= 1)) stop(paste("tau must be between 0 and 1", "\n", ""))
  u <- runif(n)
  alpha <- log((2 / pi) * asin(tau / (2 - tau))) / log(mu)
  y <- (2 / pi * asin(u / (2 - u)))^(1 / alpha)
  return(y)
}

# Survival Function
sRURS <- function(q, mu, tau) {
  1 - pRURS(q, mu, tau)
}

# Hazard Function
hRURS <- function(y, mu, tau) {
  dRURS(y, mu, tau) / sRURS(y, mu, tau)
}

# Bowley Skewness (mu and tau only)
bowley <- function(mu, tau) {
  Q1 <- qRURS(p = 1 / 4, mu = mu, tau = tau)
  Q2 <- qRURS(p = 2 / 4, mu = mu, tau = tau)
  Q3 <- qRURS(p = 3 / 4, mu = mu, tau = tau)
  k1 <- Q3 - 2 * Q2 + Q1
  k2 <- Q3 - Q1
  bowley <- k1 / k2
  return(bowley)
}

# Moors Kurtosis (mu and tau only)
moors <- function(mu, tau) {
  Q78 <- qRURS(p = 7 / 8, mu = mu, tau = tau)
  Q58 <- qRURS(p = 5 / 8, mu = mu, tau = tau)
  Q38 <- qRURS(p = 3 / 8, mu = mu, tau = tau)
  Q18 <- qRURS(p = 1 / 8, mu = mu, tau = tau)
  Q34 <- qRURS(p = 3 / 4, mu = mu, tau = tau)
  Q14 <- qRURS(p = 1 / 4, mu = mu, tau = tau)
  k1 <- Q78 - Q58 + Q38 - Q18
  k2 <- Q34 - Q14
  moors <- k1 / k2
  return(moors)
}

# ============================================================================
# USER INTERFACE
# ============================================================================

ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "journal"),
  withMathJax(),
  titlePanel("QURS - Quantile Unit Rational Sine distribution"),
  h5("Skewness, Kurtosis and Distribution"),
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(
        condition = "input.tab_ativa == 'Skewness and Kurtosis'",
        h4("Skewness and Kurtosis"),
        helpText("$$\\mu_{\\tau}$$"),
        sliderInput("mu_max", NULL, min = 0.01, max = 0.99, value = 0.50, step = 0.01),
        helpText("$$\\tau$$"),
        sliderInput("tau_max", NULL, min = 0.01, max = 0.99, value = 0.50, step = 0.01)
      ),
      conditionalPanel(
        condition = "input.tab_ativa == 'Distribution'",
        h4("Distribution"),
        helpText("$$\\mu_\\tau$$"),
        sliderInput("mu_d", NULL, min = 0.01, max = 0.99, value = 0.50, step = 0.01),
        helpText("$$\\tau$$"),
        sliderInput("tau_d", NULL, min = 0.01, max = 0.99, value = 0.50, step = 0.01),
        helpText("$$y$$"),
        sliderInput("y_max", NULL, min = 0.01, max = 0.99, value = 0.99, step = 0.01)
      )
    ),
    mainPanel(
      tabsetPanel(
        id = "tab_ativa",
        tabPanel(
          "Skewness and Kurtosis",
          plotOutput("bowley"),
          br(),
          plotOutput("moors")
        ),
        tabPanel(
          "Distribution",
          fluidRow(
            column(6, plotOutput("dist_dens")),
            column(6, plotOutput("dist_cdf")),
            column(6, plotOutput("dist_surv")),
            column(6, plotOutput("dist_haz"))
          )
        )
      )
    )
  )
)

# ============================================================================
# SERVER
# ============================================================================

server <- function(input, output) {
  
  # Data for Bowley - mu on x-axis, tau as grouping
  db <- reactive({
    req(input$mu_max, input$tau_max)
    mu_seq <- seq(0.01, input$mu_max, length.out = 100)
    tau_seq <- seq(0.01, input$tau_max, length.out = 8)
    
    map_df(tau_seq, function(t) {
      bowley_vals <- sapply(mu_seq, function(m) {
        tryCatch(bowley(mu = m, tau = t),
                 error = function(e) NA_real_)
      })
      data.frame(
        mu = mu_seq,
        v = bowley_vals,
        tau = t
      )
    })
  })
  
  # Data for Moors - mu on x-axis, tau as grouping
  dm <- reactive({
    req(input$mu_max, input$tau_max)
    mu_seq <- seq(0.01, input$mu_max, length.out = 100)
    tau_seq <- seq(0.01, input$tau_max, length.out = 8)
    
    map_df(tau_seq, function(t) {
      moors_vals <- sapply(mu_seq, function(m) {
        tryCatch(moors(mu = m, tau = t),
                 error = function(e) NA_real_)
      })
      data.frame(
        mu = mu_seq,
        v = moors_vals,
        tau = t
      )
    })
  })
  
  # Data for Distribution plots
  dd <- reactive({
    req(input$mu_d, input$tau_d, input$y_max)
    y_start <- min(0.001, input$y_max / 2)
    y_seq <- seq(y_start, input$y_max, length.out = 200)
    
    dens <- tryCatch(dRURS(y_seq, input$mu_d, input$tau_d), error = function(e) NA_real_)
    cdf  <- tryCatch(pRURS(y_seq, input$mu_d, input$tau_d), error = function(e) NA_real_)
    surv <- tryCatch(sRURS(y_seq, input$mu_d, input$tau_d), error = function(e) NA_real_)
    haz  <- tryCatch(hRURS(y_seq, input$mu_d, input$tau_d), error = function(e) NA_real_)
    
    df <- data.frame(
      y = y_seq,
      dens = as.numeric(dens),
      cdf  = as.numeric(cdf),
      surv = as.numeric(surv),
      haz  = as.numeric(haz)
    )
    df[is.finite(df$dens) & is.finite(df$cdf) &
         is.finite(df$surv) & is.finite(df$haz), , drop = FALSE]
  })
  
  # Helper for Bowley/Moors plots
  plot_gg <- function(df, ylab, tau_max) {
    req(df, nrow(df) > 0)
    tau_seq <- seq(0.01, tau_max, length.out = 8)
    fin <- df[is.finite(df$v), ]
    req(nrow(fin) > 0)
    
    ggplot(fin, aes(x = mu, y = v, colour = as.factor(tau))) +
      geom_line(size = 0.8) +
      scale_color_manual(
        values = viridis::mako(length(tau_seq), alpha = 0.8, begin = 0, end = 0.9),
        name = expression(tau),
        labels = sprintf("%.2f", tau_seq)
      ) +
      labs(y = ylab, x = expression(mu[tau]), title = "") +
      theme_bw(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box = "vertical",
        legend.position = "right"
      ) +
      guides(colour = guide_legend(
        direction = "vertical", title.position = "top",
        override.aes = list(size = 1)
      )) +
      coord_cartesian(
        ylim = c(min(fin$v, na.rm = TRUE), max(fin$v, na.rm = TRUE)),
        xlim = c(0, input$mu_max)
      )
  }
  
  output$bowley <- renderPlot(plot_gg(db(), "Bowley Skewness", input$tau_max))
  output$moors  <- renderPlot(plot_gg(dm(), "Moors Kurtosis", input$tau_max))
  
  # Distribution plots
  output$dist_dens <- renderPlot({
    df <- dd()
    req(nrow(df) > 0)
    ggplot(df, aes(x = y, y = dens)) +
      geom_line(color = "#440154", size = 1) +
      labs(title = "Probability Density Function", x = "y", y = "f(y)") +
      theme_minimal()
  })
  
  output$dist_cdf <- renderPlot({
    df <- dd()
    req(nrow(df) > 0)
    ggplot(df, aes(x = y, y = cdf)) +
      geom_line(color = "#3B528B", size = 1) +
      labs(title = "Cumulative Distribution Function", x = "y", y = "F(y)") +
      theme_minimal()
  })
  
  output$dist_surv <- renderPlot({
    df <- dd()
    req(nrow(df) > 0)
    ggplot(df, aes(x = y, y = surv)) +
      geom_line(color = "#21918C", size = 1) +
      labs(title = "Survival Function", x = "y", y = "S(y)") +
      theme_minimal()
  })
  
  output$dist_haz <- renderPlot({
    df <- dd()
    req(nrow(df) > 0)
    ggplot(df, aes(x = y, y = haz)) +
      geom_line(color = "brown", size = 1) +
      labs(title = "Hazard Function", x = "y", y = "h(y)") +
      theme_minimal()
  })
}

shinyApp(ui = ui, server = server)