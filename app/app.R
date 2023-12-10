library(shiny)
library(DT)
library(tidyverse)
library(dplyr)
library(matrixStats)
library(gplots)
library(ggplot2)
library(gridExtra)
library(data.table)
library(beeswarm)
library(DESeq2)
library(igraph)
library(EnhancedVolcano)

choices_list <- c("csv", "tsv", "text")

ui <- fluidPage(
  tags$head(
    tags$style(
      HTML("
        .tab-content {
          height: calc(100vh - 50px); /* Adjust 50px based on your layout */
          overflow-y: auto; /* Add scrollbar if needed */
        }
      ")
    )
  ),
  navbarPage(
    "Project for BF591",
    tabPanel(
      "Information",
      uiOutput("info_text")
    ),
    tabPanel(
      "Samples",
      actionButton("samples", "Load Samples"),
      conditionalPanel(
        condition = "input.samples > 0",
        selectInput("samples_format", "Select file format", choices = choices_list),
        fileInput("samples_data", "Upload counts data file to create the sample data and view its summary."),
        tabsetPanel(
          tabPanel("Summary", DTOutput("sample_summary_table")),
          tabPanel("Sample Information", DTOutput("sample_info_table"),
                   downloadButton("downloadSampleInfo", "Download Sample MetaData")),
          tabPanel("Plots", plotOutput("sample_plots"))
        )
      )
    ),
    tabPanel(
      "Counts",
      actionButton("counts", "Load Counts"),
      conditionalPanel(
        condition = "input.counts > 0",
        selectInput("counts_format", "Select file format", choices = choices_list),
        fileInput("counts_data", "Upload counts data file"),
        sliderInput("variance_slider", "Select genes with at least X percentile of variance:", min = 0, max = 100, value = 50),
        sliderInput("non_zero_slider", "Select genes with at least X samples that are non-zero:", min = 0, max = 100, value = 50, step = 4),
        tabsetPanel(
          tabPanel("Summary of Filtering", dataTableOutput("counts_summary_table")),
          tabPanel("Scatter Plots", plotOutput("scatter_plot_counts")),
          tabPanel("HeatMap", plotOutput("heat_map_counts")),
          tabPanel("PCA",
                   textInput("selected_pcs", "Enter Principal Components (comma-separated):", value = "2"),
                   plotOutput("pca_scatter_counts"),
                   DTOutput("pca_table"))
        )
      )
    ),
    tabPanel(
      "Differential Expression",
      actionButton("DE", "Start DE Analysis"),
      conditionalPanel(
        condition = "input.DE > 0",
        selectInput("counts_de_format", "Select counts file format", choices = choices_list),
        fileInput("counts_data_de", "Upload counts data file"),
        selectInput("col_format", "Select metadata file format", choices = choices_list),
        fileInput("col_data", "Upload metadata file"),
        selectInput("condition1", "Select the first condition", choices = NULL),
        selectInput("condition2", "Select the second condition", choices = NULL),
        tabsetPanel(
          tabPanel("DESeq2 Results Table", DTOutput("de_table")),
          tabPanel("Volcano Plot", plotOutput("volcano_plot"))
        )
      )
    ),
    tabPanel(
      "Correlation Network Analysis",
      actionButton("CNA", "Start Correlation Network Analysis"),
      conditionalPanel(
        condition = "input.CNA > 0",
        selectInput("counts_cna_format", "Select counts file format", choices = choices_list),
        fileInput("counts_data_cna", "Upload counts data file"),
        textInput("selected_genes", "Enter gene IDs (comma-separated or space separated):", value = "Gnai3"),
        sliderInput("correlation_slider", "Select the minimum correlation for drawing an edge:", min = 0, max = 1, value = 0.5),
        tabsetPanel(
          tabPanel("Clustered HeatMap", plotOutput("clustered_heatmap")),
          tabPanel("Correlation Network", plotOutput("correlation_network")),
          tabPanel("Metrics", DTOutput("information_table"))
        )
      )
    )
  )
)


server <- function(input, output, session) {
  
  load_samples <- eventReactive(input$samples, {
    if (is.null(input$samples_data))
      return(NULL)
    else if (input$samples_format == "csv") {
      df <- read.csv(input$samples_data$datapath)
      return(df)
    } else {
      df <- read.delim(input$samples_data$datapath)
      return(df)
    }
  })
  
  get_samples_meta_data <- eventReactive(input$samples,{
    if (is.null(input$samples_data))
      return(NULL)
    else if (input$samples_format == "csv") {
      df <- read.csv(input$samples_data$datapath)
      columns_to_remove <- c("Gene", "GeneID", "Coordinates")
      df <- df[, !tolower(colnames(df)) %in% tolower(columns_to_remove), drop = FALSE]
      meta_data <- meta_info_from_labels(colnames(df))
      return(meta_data)
    } else {
      df <- read.delim(input$samples_data$datapath)
      columns_to_remove <- c("Gene", "GeneID", "Coordinates")
      df <- df[, !tolower(colnames(df)) %in% tolower(columns_to_remove), drop = FALSE]
      meta_data <- meta_info_from_labels(colnames(df))
      return(meta_data)
    }
  })
  
  
  get_col_data <- eventReactive(input$DE, {
    if (is.null(input$col_data))
      return(NULL)
    else {
      if (input$col_format == "csv") {
        df <- read.csv(input$col_data$datapath)
      } else {
        df <- read.delim(input$col_data$datapath)
      }
      return(df)
    }
  })
  
  
  sample_summary <- reactive({
    if (is.null(get_samples_meta_data()))
      return(NULL)
    
    sample_data <- get_samples_meta_data()
    
    summary_df <- tibble(
      Column_Name = colnames(sample_data),
      Type = sapply(sample_data, class),
      Distinct_Values = sapply(sample_data, function(x) paste(unique(x)))
    )
    
    return(summary_df)
  })
  
  
  load_counts <- eventReactive(input$counts, {
    if (is.null(input$counts_data))
      return(NULL)
    else if (input$counts_format == "csv") {
      df <- read.csv(input$counts_data$datapath)
      return(df)
    } else {
      df <- read.delim(input$counts_data$datapath)
      return(df)
    }
  })
  
  load_counts_de <- eventReactive(input$DE, {
    if (is.null(input$counts_data_de))
      return(NULL)
    else if (input$counts_de_format == "csv") {
      df <- read.csv(input$counts_data_de$datapath)
      return(df)
    } else {
      df <- read.delim(input$counts_data_de$datapath)
      return(df)
    }
  })
  
  observe({
    # Update the max value of non_zero_slider based on the actual sample size
    if (!is.null(load_counts())) {
      sample_size <- ncol(load_counts()[,-c(1:3)])
      updateSliderInput(session, "non_zero_slider", max = sample_size, value = 0, step = 4)
    }
  })
  
 condition_choices <- reactive({
    coldata <- get_col_data()
    if (!is.null(coldata)) {
      unique(coldata$timepoint)
    } else {
      character(0)
    }
  })
  
  # Update the choices whenever the reactive expression changes
  observe({
    updateSelectInput(session, "condition1", choices = condition_choices())
    updateSelectInput(session, "condition2", choices = condition_choices())
  })
  
  load_counts_cna <- eventReactive(input$CNA, {
    if (is.null(input$counts_data_cna))
      return(NULL)
    else if (input$counts_cna_format == "csv") {
      df <- read.csv(input$counts_data_cna$datapath)
      return(df)
    } else {
      df <- read.delim(input$counts_data_cna$datapath)
      return(df)
    }
  })
  
  # Reactive function to filter genes based on variance and non-zero sample counts
  filtered_genes <- reactive({
    if (is.null(load_counts()))
      return(NULL)
    
    counts_data <- load_counts()
    
    non_zero_threshold <- input$non_zero_slider
    variance_threshold <- input$variance_slider
    
    # Compute variances for each gene
    variances <- apply(counts_data[, -c(1:3)], 1, var, na.rm = TRUE)
    
    # Create a data frame with Gene and variance
    gene_variances <- data.frame(Gene = counts_data$Gene, Variance = variances)
    
    # Filter genes based on variance
    genes_filtered_by_variance <- filter(gene_variances, Variance > quantile(Variance, variance_threshold / 100))
    
    # Filter genes based on non-zero sample counts
    non_zero_counts <- rowSums(counts_data[, -c(1:3)] > 0, na.rm = TRUE)
    non_zero_genes <- filter(gene_variances, non_zero_counts >= non_zero_threshold)
    
    # Combine both filters
    final_filtered_genes <- intersect(genes_filtered_by_variance$Gene, non_zero_genes$Gene)
    
    # Return the filtered genes
    return(counts_data[counts_data$Gene %in% final_filtered_genes, ])
  })
  
  
  
  create_diagnostics <- reactive({
    if (is.null(load_counts()))
      return(NULL)
    
    counts_data <- load_counts()
    filtered_data <- filtered_genes()
    log_transform <- function(x) log2(x + 1)
    counts_data$median_count <- rowMedians(as.matrix(counts_data[,-c(1:3)]), na.rm=TRUE)
    counts_data$log_variance <- log_transform(apply(counts_data[,-c(1:3)], 1, var, na.rm = TRUE))
    
    # Check if filtered_data has the "Gene" column
    counts_data$color <- ifelse(counts_data$Gene %in% filtered_data$Gene, "Passing Filter", "Filtered Out")
    
    # Calculate num_zeros
    counts_data$num_zeros <- rowSums(counts_data[,-c(1:3)] == 0, na.rm = TRUE)
    
    
    # Set a more visually appealing color palette
    my_palette <- c("Filtered Out" = "#CCCCCC", "Passing Filter" = "#3498db")
    
    p1 <- ggplot(counts_data, aes(x = median_count, y = log_variance, color = color)) +
      geom_point(size = 3, alpha = 0.7, shape = 16) +
      scale_color_manual(values = my_palette) +
      labs(title = "Median Count vs Variance", x = "Median Count", y = "Log Variance") +
      theme_minimal() +
      theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
    
    p2 <- ggplot(counts_data, aes(x = median_count, y = num_zeros, color = color)) +
      geom_point(size = 3, alpha = 0.7, shape = 16) +
      scale_color_manual(values = my_palette) +
      labs(title = "Median Count vs Number of Zeros", x = "Median Count", y = "Number of Zeros") +
      theme_minimal() +
      theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
    
    return(list(p1, p2))
  })
  
  
  create_heatmap <- reactive({
    if (is.null(load_counts()))
      return(NULL)
    counts_data <- load_counts()
    filtered_data <- filtered_genes()
    genes_to_include <- filtered_data$Gene
    
    # Assuming your counts_data has a column for gene names (adjust if needed)
    counts_data_filtered <- counts_data[counts_data$Gene %in% genes_to_include, ]
    
    # Assuming "gene_column" is the name of the column in counts_data with gene names
    log_counts <- log2(counts_data_filtered[, -c(1:3)] + 1)
    
    heatmap_output <- heatmap.2(as.matrix(log_counts),
                                col = viridis::viridis(75),
                                scale = "row",
                                cluster.cols = FALSE,
                                main = "Clustered Heatmap of Selected Genes",
                                labRow = counts_data_filtered$GeneID)
    
    return(heatmap_output)
})
  
  # Inside the create_pca_plot reactive function
  create_pca <- reactive({
    if (is.null(filtered_genes()))
      return(NULL)
    
    filtered_data <- filtered_genes()
    
    # Assuming "Gene" is the first column; adjust if needed
    filtered_data<-  filtered_data[, -c(1:3)]
    
    transposed_data = transpose(data.frame(filtered_data))
    
    
    # Perform PCA directly using prcomp
    pca_result <- prcomp(transposed_data, scale. = TRUE)
    
    # Extract the top PCs
    selected_pcs <- as.numeric(strsplit(input$selected_pcs, ",")[[1]])
    selected_pcs <- selected_pcs[selected_pcs <= ncol(pca_result$x)]
    top_pcs <- pca_result$x[, selected_pcs]
    
    # Add % variance explained to the column names
    var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 2)
    label <- paste0("PC", selected_pcs, " (", var_explained[selected_pcs], "%)")
    
    # Create scatter plot
    beeswarm_plot <- beeswarm(
      as.matrix(top_pcs),
      pch = 16,
      col = viridis::viridis(length(selected_pcs), option = "D"),
      main = "Beeswarm Plot of Selected PCs"
    )
    
    table_data <- data.frame(
      "PC" = paste0("PC", selected_pcs),
      "Percent Variance Explained" = var_explained[selected_pcs]
    )
    
    table_output <- datatable(table_data, caption = "Percentage of Variance Explained by Selected PCs")
    
    return(list(top_pcs = top_pcs, plot = beeswarm_plot,table = table_output))
  })
  
  
  run_de <- reactive({
    
    input_condition1 <- shiny::isolate(input$condition1)
    input_condition2 <- shiny::isolate(input$condition2)
    
    if (is.null(load_counts_de()) || is.null(get_col_data()) || is.null(input$condition1) || is.null(input$condition2)
      || input$condition1 == input$condition2){
      return(NULL)}
    
    og_counts <- load_counts_de()
    
    
    columns_to_remove <- c("Gene", "GeneID", "Coordinates")
    counts <- og_counts[, !tolower(colnames(og_counts)) %in% tolower(columns_to_remove), drop = FALSE]
    

    
    colData <- get_col_data()
    
    colnames(colData)[colnames(colData) == "timepoint"] <- "condition"
    
    colData$condition = factor(colData$condition)
    
    
    condition1 <- input$condition1
    condition2 <- input$condition2
    
    # Filter colData and counts to include only samples from the specified conditions
    coldata_subset <- colData[colData$condition %in% c(condition1, condition2), ]
    counts_subset <- counts[, colnames(counts) %in% coldata_subset$sample_name]
    
    
    
   
    
    # Create DESeqDataSet
    dds <- DESeqDataSetFromMatrix(countData = counts_subset,
                                  colData = coldata_subset,
                                  design = ~ condition)
    
    # Filtering based on count threshold (adjust as needed)
    dds <- dds[rowSums(counts(dds)) >= 10, ]
    
    gene_info <- og_counts$gene
    
    # Run DESeq2 analysis
    dds <- DESeq(dds)
    
    res <- results(dds)
    
    # Merge results and gene_info by row names
    res_with_gene <- merge(data.frame(res), gene_info, by = "row.names")
    
    colnames(res_with_gene)[ncol(res_with_gene)] <- "Gene"
    
    # Move the "GeneName" column to the first position
    res_with_gene <- res_with_gene[, c("Gene", setdiff(colnames(res_with_gene), "Gene"))]
    
    # Remove the redundant row names column
    res_with_gene$Row.names <- NULL
    
    title <- paste("Conditions Tested: ", input_condition1, " and ", input_condition2)
    
    # Add the title to the data frame
    res_with_gene$title <- title
    
    return(res_with_gene)
  })
  
  debounced_run_de <- debounce(run_de, 500)
  
  create_volcano_plot <- reactive({
    # Check if DE results are available
    if (is.null(debounced_run_de())) {
      return(NULL)
    }
    
    de_results <- run_de()
    
    # Create a volcano plot
    volcano_plot <- EnhancedVolcano(
      de_results,
      lab = de_results$Gene,
      x = 'log2FoldChange',
      y = 'padj',
      xlim = c(-5, 5),
      title = 'Volcano Plot of Differential Expression',
      subtitle = 'Adjusted p-value vs. Log2 Fold Change',
      pCutoff = 0.05,
      FCcutoff = 2,
      pointSize = 3,
      labSize = 2.5,
      legendLabels = c('Not Significant', 'Significant Up', 'Significant Down'),
      legendPosition = 'right'
    )
    
    return(volcano_plot)
  })
  
  # Reactive function to generate clustered heatmap for selected genes
  create_clustered_heatmap <- reactive({
    # Check if counts are available
    if (is.null(load_counts_cna())) {
      return(NULL)
    }
    
    counts_data <- load_counts_cna()
    
    selected_genes <- tools::toTitleCase(strsplit(input$selected_genes, "[,]")[[1]])
    
    not_found_genes <- setdiff(selected_genes, counts_data$GeneID)
    if (length(not_found_genes) > 0) {
      stop(paste("Genes '", paste(not_found_genes, collapse = "', '"), "' not found in the counts data.", sep = ""))
    }
    
    counts_data_selected <- counts_data[counts_data$GeneID %in% selected_genes, ]
    
    # Assuming "gene_column" is the name of the column in counts_data with gene names
    log_counts <- log2(counts_data_selected[, -c(1:3)] + 1)
    
    # Create clustered heatmap
    heatmap.2(as.matrix(log_counts),
            col = viridis::viridis(75),
            scale = "row",
            cluster.cols = FALSE,
            main = "Clustered Heatmap of Selected Genes",
            labRow = counts_data_selected$GeneID)
  })
  
  
  create_correlation_network <- reactive({
    if (is.null(load_counts_cna())) {
      return(NULL)
    }
    counts_data <- load_counts_cna()
    
    # Assuming your counts_data has a column for gene names (adjust if needed)
    selected_genes <- tools::toTitleCase(strsplit(input$selected_genes, ",")[[1]])
    
    not_found_genes <- setdiff(selected_genes, counts_data$GeneID)
    if (length(not_found_genes) > 0) {
      stop(paste("Genes '", paste(not_found_genes, collapse = "', '"), "' not found in the counts data.", sep = ""))
    }
    
    counts_data_selected <- counts_data[counts_data$GeneID %in% selected_genes, ]
    gene_ids_selected <- counts_data_selected$GeneID
    threshold <- input$correlation_slider
    
    correlation_matrix <- cor(t(counts_data_selected[, -c(1:3)]))
    correlation_matrix[abs(correlation_matrix) < threshold] <- 0
    correlation_matrix[is.na(correlation_matrix)] <- 0
    
    graph <- graph.adjacency(correlation_matrix, weighted = TRUE ,mode = "directed")
    
    V(graph)$name <- gene_ids_selected
    
    constant_color <- "orange"
    V(graph)$color <- constant_color
    
    crnp <- plot(
      graph,
      layout = layout_with_fr(graph),
      vertex.label = V(graph)$name,
      vertex.color = V(graph)$color,
      main = "Correlation Network"
    )
  return(list(p_obj=p,graph_obj=graph))
})
  
  create_graph_information_table <- reactive({
    graph <- create_correlation_network()$graph_obj
    
    degree <- degree(graph)
    closeness <- closeness(graph)
    betweenness <- betweenness(graph)
    
    information_df <- data.frame(
      Gene = V(graph)$name,
      Degree = degree,
      ClosenessCentrality = closeness,
      BetweennessCentrality = betweenness
    )
    return(information_df)
  })
  
  output$information_table <- renderDataTable({
    datatable(create_graph_information_table(), caption = "Metrics for Each Gene in the Graph Network")
  })
  
  output$correlation_network <- renderPlot({
    create_correlation_network()$p_obj
  })
  
  # Render the clustered heatmap
  output$clustered_heatmap <- renderPlot({
    create_clustered_heatmap()
  })
  
  
  output$counts_summary_table <- renderDataTable({
    counts_data <- load_counts()
    filtered_genes_data <- filtered_genes()
    
    if (is.null(filtered_genes_data)) {
      return(NULL)
    }
    
    total_genes <- nrow(counts_data)
    genes_passing_filter <- nrow(filtered_genes_data)
    genes_not_passing_filter <- total_genes - genes_passing_filter
    
    summary_df <- data.frame(
      "Variables" = c(
        "Number of Samples",
        "Total Number of Genes",
        "Number of Genes Passing Filter",
        "Percent of Genes Passing Filter",
        "Number of Genes Not Passing Filter",
        "Percent of Genes Not Passing Filter"
      ),
      "Values" = c(
        ncol(counts_data)-3,
        total_genes,
        genes_passing_filter,
        round(100 * (genes_passing_filter / total_genes), 2),
        genes_not_passing_filter,
        round(100 * (genes_not_passing_filter / total_genes), 2)
      )
    )
    
    datatable(summary_df, caption = "Summary Table")
  })
  

output$sample_summary_table <- renderDT({
  datatable(sample_summary(), caption = "Sample Summary Table")
})

output$scatter_plot_counts <- renderPlot({
    scatter_plots <- create_diagnostics()
    grid.arrange(grobs = scatter_plots, ncol = 2)
})

output$heat_map_counts <- renderPlot({
  create_heatmap()
})

output$pca_scatter_counts <- renderPlot({
  create_pca()$plot
})

output$pca_table <- renderDT({
  create_pca()$table
})

output$de_table<- renderDT({
  debounced_run_de()
})

output$volcano_plot <- renderPlot({
  create_volcano_plot()
})
# Display sample information in a table
output$sample_info_table <- renderDT({
  data_to_show <- if (!is.null(get_samples_meta_data())) {
    get_samples_meta_data()
  } else if (!is.null(load_samples())) {
    load_samples()
  } else {
    NULL
  }
})

output$downloadSampleInfo <- downloadHandler(
  filename = function() {
    paste("sample_info_table", ".csv", sep="")
  },
  content = function(file) {
    if (!is.null(get_samples_meta_data())) {
      write.csv(get_samples_meta_data(), file, row.names = FALSE)
    } else if (!is.null(load_samples())) {
      write.csv(load_samples(), file, row.names = FALSE)
    }
  },
  contentType = 'text/csv; charset=utf-8'
)


output$sample_plots <- renderPlot({
  if (!is.null(get_samples_meta_data())) {
    sample_data <- get_samples_meta_data()
    counts <- table(sample_data$timepoint, sample_data$replicate)
    barplot(counts, 
            main = "Number of Samples by Time Point and Replicate",
            xlab = "Time Point",
            ylab = "Number of Samples",
            col = c("skyblue", "lightgreen", "lightcoral"), 
            legend = FALSE,
            beside = TRUE,
            names.arg = sample_data$timepoint)
    replicate_levels <- unique(sample_data$replicate)
    legend("topright", legend = replicate_levels, fill = c("skyblue", "lightgreen", "lightcoral"),title = "Replicates")
  }
})

output$info_text <- renderUI({
  HTML(paste(
    "<div><h3>About the App</h3>",
    "<p>This Shiny app is designed for Project BF591 and aims to assist with the analysis of Bulk RNA Sequencing Data.</p>",
    
    
    "<h3>Data Source</h3>",
    "<p>The analysis in this app is based on data from the paper:</p>",
    "<p><b>Paper Title:</b> Transcriptional Reversion of Cardiac Myocyte Fate During Mammalian Cardiac Regeneration</p>",
    "<p><b>Authors:</b> [O’Meara et al.,2015]</p>",
    "<p><b>Journal:</b> [Circulation Research]</p>",
    "<p><b>The paper can be accessed at :</b> <a href='https://www.ahajournals.org/doi/10.1161/CIRCRESAHA.116.304269' target='_blank'>[Link to Paper]</a></p>",
    "<p><b>The dataset with GEO accession number as GSE64403 is available at:</b> <a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64403' target='_blank'>[Link to Dataset]</a></p>",
    
    "<h3>Information Tab</h3>",
    "<p>This tab provides information about the app and a brief overview of each of the main tabs.</p></div>",
    
    "<h3>Tab 1: Samples</h3>",
    "<p>This tab allows you to load and explore information about your samples, including a summary, individual details, and plots.
         It creates this sample metadata from the sample names present in the FPKM data.</p>",
    
    "<h3>Tab 2: Counts</h3>",
    "<p>Load and filter counts data, explore summary statistics, scatter plots, heatmaps, and PCA analysis.</p>",
    
    "<h3>Tab 3: Differential Expression</h3>",
    "<p>Perform DE analysis using DESeq2 with options to upload counts data and metadata. View results in a table and a volcano plot.</p>",
    
    "<h3>Tab 4: Correlation Network Analysis</h3>",
    "<p>Explore correlation networks for selected genes, including a clustered heatmap and network metrics.</p>"
  ))
})



}



timepoint_from_sample <- function(x) {
  return(substring(x, 2, regexpr("_", x) - 1))
}

sample_replicate <- function(x) {
  return(as.integer(substring(x, regexpr("_", x) + 1)))
}

meta_info_from_labels <- function(sample_names) {
  replicate <- sapply(sample_names, FUN = sample_replicate)
  timepoint <- sapply(sample_names, FUN = timepoint_from_sample)
  final <- tibble(
    sample_name = sample_names,
    timepoint = timepoint,
    replicate = replicate
  )
  return(final)
}



options(shiny.maxRequestSize = 50 * 1024^2)
shinyApp(ui, server)
