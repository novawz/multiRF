library(shiny)
library(ggplot2)
library(cluster)

server <- function(input, output) {

  data0 <- eventReactive(input$example, {
    tcga_brca
  })

  output$datasummary <- renderTable({
    data <- data0()
    dim_df <- data %>%
      purrr::map(~dim(.)) %>%
      Reduce(rbind,.) %>%
      t() %>%
      data.frame()
    data_name <- names(data)
    rownames(dim_df) <- c("Number of Samples", "Number of Features")
    colnames(dim_df) <- data_name

    dim_df
  }, rownames = T, bordered = T, hover = T, striped = T)

  output$dataplot1 <- renderPlot({
    data <- data0()
    all_samples <- Reduce(union, data %>% purrr::map(~rownames(.)))

    mat <- matrix("Available", ncol = length(all_samples), nrow = length(data))
    rownames(mat) <- names(data)
    sapply(1:length(data), function(i) {
      mat[i, which(!all_samples %in% rownames(data[[i]]))] <<- "Missing"
    })
    plotDecoding <- reshape2::melt(mat)
    ggplot(data.frame(plotDecoding), aes(x = Var2, y = Var1)) +
      theme_bw() +
      theme(legend.position = "right",
              legend.title = element_blank()) +
      geom_tile(aes(fill = value), color = "white") +
      xlab("Samples") +
      ylab("Cohorts") +
      ggsci::scale_fill_d3()
  }, width = 800, height = 300, execOnResize = T)


  mod <- eventReactive(input$run, {
    data <- data0()
    if("Scale" %in% input$para) {
      scaled <- T
    } else {
      scaled <- F
    }

    if("Normalized weights" %in% input$para) {
      normalized <- T
    } else {
      normalized <- F
    }

    if("Direct" %in% input$para) {
      direct <- T
    } else {
      direct <- F
    }
    mrf3_init(data, ntree = input$ntree,
              scale = scaled,
              normalized = normalized,
              direct = direct,
              yprob = input$yprob)
  })
  output$initialModel <- renderText({
    model <- mod()
    paste0("Finish running model! Optimal connection is: ", paste0(model$connection, collapse = ", "), ". Please click on the Results tab." )
  })

  # Render scatter plot for Variable Selection section - Model Fitting
  output$scatterPlot1 <- renderPlot({
    model <- mod()
    plot_weights(
      model$weights
    )
  })

  output$newSelect1 <- renderUI({
    if (input$method == "test") {
      return(sliderInput("slider", "Level:", 0.01, 0.1, 0.05, width = "80%"))
    }
    if(input$method == "mixture") {
      return(selectInput("m1", "Model 1", c("normal", "truncn"), width = "50%"))
    }
  })

  output$newSelect2 <- renderUI({
    if(input$method == "mixture") {
      selectInput("m2", "Model 2", c("normal", "gamma"), width = "50%")
    }
  })

  modvs <- eventReactive(input$run2, {
    data <- data0()
    mod0 <- mod()
    if("Scale" %in% input$para) {
      scaled <- T
    } else {
      scaled <- F
    }
    mrf3_vs(mod0, data,
            method = input$method,
            c1 = input$m1,
            c2 = input$m2,
            level = input$slider,
            scale = scaled)
  })



  # Render scatter plot for Variable Selection section
  output$table2 <- renderTable({
    model <- modvs()
    data.frame(
      "Data" = names(model$weights),
      "Number of Features Selected" = unlist(model$weights %>% purrr::map(~length(.[. != 0])))
    )
  }, rownames = F, bordered = T, hover = T, striped = T, align = "c", spacing = "l")
  output$scatterPlot2 <- renderPlot({
    model <- modvs()
    plot_weights(model$weights, top = NULL)
  })

  output$dataplot2 <- renderPlot({
    model <- modvs()
    g <- plyr::llply(
      names(model$dat.list),
      .fun = function(d) {
        df <- model$dat.list[[d]]
        p <- plot_heatmap(df, breaks = seq(-min(c(-min(df[df < 0]),max(df[df > 0]))),
                                          min(c(-min(df[df < 0]),max(df[df > 0]))),
                                          length.out = 100),
                          main = d)
        p[[4]]
      }
    )

    gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs= g,ncol=3))
  })

  # Render clustering plot for Clustering Analysis section
  output$clusterPlot <- renderPlot({

  })
}
