library(shiny)
library(tidyverse)
library(lubridate)

# ---- Load and combine all yearly CSV files ----

csv_files <- list.files("data", pattern = "^WIZ_NZ_[0-9]{4}\\.csv$", full.names = TRUE)

seismic_all <- list()

for (f in csv_files) {
  d <- read_csv(f, show_col_types = FALSE)
  
  # Build the timestamp in UTC first, then convert to NZ time.
  # (Passing tz = "Pacific/Auckland" directly into as.POSIXct() here would
  # shift every timestamp by ~12-13 hours, since the "1970-01-01" origin
  # would itself be treated as NZ time instead of UTC.)
  utc_time <- as.POSIXct(d$unix_timestamp, origin = "1970-01-01", tz = "UTC")
  d$datetime_nz <- with_tz(utc_time, tzone = "Pacific/Auckland")
  
  seismic_all[[f]] <- d
}

seismic_all <- bind_rows(seismic_all)

# Derive year, month, and quarter all from the same NZ-local timestamp.
# (Tagging "year" from the filename instead would cause a mismatch: NZDT is
# UTC+13, so the last UTC readings of Dec 31 land on Jan 1 NZ time of the
# *next* year. Those rows would then carry the old file's year but a "Jan"
# month, which drags the Jan/Q1 facet's x-axis out to cover the whole year.)
seismic_all$year <- year(seismic_all$datetime_nz)
seismic_all$month <- month(seismic_all$datetime_nz, label = TRUE)
seismic_all$quarter <- paste0("Q", quarter(seismic_all$datetime_nz))

# ---- UI ----

ui <- fluidPage(
  titlePanel("WIZ Seismic Amplitude Time Series"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("year", "Select year:", choices = sort(unique(seismic_all$year))),
      
      selectInput("variable", "Select variable:",
                  choices = c("Displacement average" = "displacement_avg_m",
                              "RSAM average" = "rsam_avg")),
      
      selectInput("resolution", "Temporal view:",
                  choices = c("Daily", "Monthly", "Quarterly")),
      
      selectInput("display_type", "Display type:",
                  choices = c("Separate panels", "Drop-down period")),
      
      uiOutput("period_ui")
    ),
    
    mainPanel(
      plotOutput("time_plot", height = "700px"),
      tableOutput("summary_table")
    )
  )
)

# ---- Server ----

server <- function(input, output) {
  
  # Data for the selected year
  year_data <- reactive({
    seismic_all %>% filter(year == input$year)
  })
  
  # Drop-down for picking one specific month/quarter (only shown when needed)
  output$period_ui <- renderUI({
    if (input$display_type != "Drop-down period") return(NULL)
    
    if (input$resolution == "Monthly") {
      selectInput("period", "Select month:", choices = month.abb)
    } else if (input$resolution == "Quarterly") {
      selectInput("period", "Select quarter:", choices = c("Q1", "Q2", "Q3", "Q4"))
    }
  })
  
  # Data actually used for the plot and table
  plot_data <- reactive({
    df <- year_data()
    
    if (input$display_type == "Drop-down period") {
      if (input$resolution == "Monthly")   df <- df %>% filter(month == input$period)
      if (input$resolution == "Quarterly") df <- df %>% filter(quarter == input$period)
    }
    
    df
  })
  
  output$time_plot <- renderPlot({
    df <- plot_data()
    
    p <- ggplot(df, aes(x = datetime_nz, y = .data[[input$variable]])) +
      geom_line() +
      labs(
        title = paste(input$resolution, "view of", input$variable, "-", input$year),
        x = "Date",
        y = input$variable
      ) +
      theme_minimal(base_size = 14)
    
    if (input$display_type == "Separate panels") {
      if (input$resolution == "Monthly")   p <- p + facet_wrap(~ month, scales = "free_x")
      if (input$resolution == "Quarterly") p <- p + facet_wrap(~ quarter, scales = "free_x")
    }
    
    p
  })
  
  output$summary_table <- renderTable({
    df <- plot_data()
    req(nrow(df) > 0)
    
    x <- df[[input$variable]]
    x <- x[!is.na(x)]
    req(length(x) > 0)
    
    data.frame(
      Statistic = c("Minimum", "Maximum", "Average"),
      Value = formatC(c(min(x), max(x), mean(x)), format = "e", digits = 4)
    )
  })
}

shinyApp(ui, server)