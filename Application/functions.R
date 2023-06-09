#Libraries
library(data.table)
library(tidyverse)
library(GEOquery)
library(plotly)
library(shiny)
library(glue)
library(PCAtools)

ensembl <<- read_rds(file = "../ensembl_ids.rds")

# TAB 1 - Sample Exploration
## Function 1- Reading the data
GEO_data <- function(filepath = NULL, geo_id = NULL) {
    if (!is.null(geo_id)) {
        GEO_obj <- getGEO(GEO = geo_id)
        GEO_obj <- GEO_obj[[1]]
        phenodata <- pData(GEO_obj)
        phenodata <- filtering(phenodata)
        return(phenodata)
    }
    else if (!is.null(filepath)) {
        GEO_obj <- GEOquery::parseGEO(fname = filepath)
        phenodata <- pData(GEO_obj)
        phenodata <- filtering(phenodata)
        return(phenodata)
    }
    else {
        stop("Either geo_id or filepath must be provided and valid")
    }
}

## Function 2 Filtering the table
filtering <- function(data){
    subset <- data %>% select(title, geo_accession,type,ends_with(":ch1"))
    colnames(subset) <- gsub(":ch1", "", colnames(subset))
    cols_to_change <- c("age of death", "age of onset", "cag", "duration","h-v cortical score","h-v striatal score","pmi","rin","mrna-seq reads")

    # Loop through the selected columns and change their types to numeric
    for (col in cols_to_change) {
        subset[[col]] <- as.numeric(subset[[col]])
    }
    return(subset)
}

## Function 3 Summary Table generation
summaryfunction <- function(data){
    new_df <- data.frame(col_name = character(),
                         type = character(),
                         mean_or_uniquevalues = character(),
                         sd = character(),
                         stringsAsFactors = FALSE)

    for (col_name in colnames(data)) {
        # Check if the column is numeric
        if (typeof(data[[col_name]]) == "double") {
            col_mean <- mean(data[[col_name]],na.rm = TRUE)
            col_sd <- sd(data[[col_name]],na.rm = TRUE)
            # Add the output to the new data frame as a row
            new_row <- data.frame(column = col_name, type = "double" ,mean_or_uniquevalues = col_mean, sd = col_sd)
            new_df <- rbind(new_df, new_row)
        }
        if (typeof(data[[col_name]]) == "character"){
            unique <-unique(data[[col_name]])
            new_row <- data.frame(column = col_name, type = "character" ,mean_or_uniquevalues = paste(unique, collapse = ", "), sd = "null")
            new_df <- rbind(new_df, new_row)
        }
    }
    return(new_df)
}

## Plotting the graphs
sampleplot <- function(data = NULL,column = "age of death",group = "diagnosis",plottype = "Histogram"){
    title <- paste(plottype,"plot of",column,"grouped by",group)

    if(plottype == "Density"){
        plot <- ggplot(data = data)+
            geom_density(mapping = aes(x = data[[column]],fill = data[[group]]), alpha = 0.5)+
            labs(title = title,x = column,y = "Counts", fill = group)+
            theme_linedraw()
    }
    if(plottype =="Histogram"){
        plot <- ggplot(data = data)+
            geom_histogram(mapping = aes(x = data[[column]],fill = data[[group]]), alpha = 0.5)+
            labs(title = title,x = column,y = "Counts", fill = group)+
            theme_linedraw()
    }

    if(plottype == "Violin"){
        plot <- ggplot(data = data)+
            geom_violin(mapping = aes(x = data[[group]], y = data[[column]], fill = data[[group]]), alpha = 0.5)+
            labs(title = title,x = group,y = column, fill = group)+
            theme_linedraw()
    }

    # if(plottype %in% c("Density","Histogram")){
    #     plot <- ggplotly(plot)
    # }

    return(plot)
}



# TAB 2 - Counts Exploration
## Function 1
countsread <- function(filepath){
    data <- fread(filepath)%>%
        dplyr::rename("genes"="V1")
    return(data)
}
## Function 2 - Filter
filter_genes <- function(data, zero_filter_value, var_percentile) {
    # Calculate gene variances
    gene_variances <- apply(data[, -1], MARGIN = 1, var)

    # Calculate gene variance percentiles
    gene_percentiles <- rank(gene_variances) / length(gene_variances)

    # Filter genes by zero counts and variance percentiles
    zero_counts <- rowSums(data[, -1] == 0) # Assuming first column contains gene names
    filtered_data <- data[zero_counts <= zero_filter_value & gene_percentiles >= var_percentile/100, ]

    return(filtered_data)
}

## Function 3--> Diagnostic plot
## a median vs variance
med_vs_var <- function(counts_tib, perc_var){

        #make a plot tibble - calculate variance and median
        plot_tib <- counts_tib%>%
            mutate(Median = apply(counts_tib[,-1], MARGIN = 1, FUN = median),
                   Variance = apply(counts_tib[,-1], MARGIN = 1, FUN = var))
        perc_val <- quantile(plot_tib$Variance, probs = perc_var/100)   #calculate percentile
        plot_tib <- plot_tib %>% mutate(threshold = case_when(Variance >= perc_val ~ "TRUE", TRUE ~ "FALSE")) #sort genes by percentile
        #plot scatter plot
        cols <- c("FALSE" = "red", "TRUE" = "black")
        scatter <- ggplot(plot_tib, aes(Median, Variance))+
            geom_point(aes(color=threshold), alpha=0.75)+
            scale_color_manual(values = cols)+
            labs(title = 'Plot of Median vs Variance.', subtitle = "Genes filtered out are in red. X and Y axes are log-scaled.")+
            scale_y_log10()+
            scale_x_log10()+
            theme_linedraw()+
            theme(legend.position = 'bottom')
        return(scatter)
}

## b
# med_vs_nz <- function(counts_tib, nz_genes){
#     tot_samp <- ncol(counts_tib)-1  #store original number of samples
#     #make a plot tibble
#     plot_tib <- counts_tib %>%
#         mutate(Median = na_if(apply(counts_tib[,-1], MARGIN = 1, FUN = median), 0))  #calculate median, convert 0 to NA
#     plot_tib$no_zeros <- rowSums(is.na(plot_tib))  #make new col, with counts.
#     plot_tib <- plot_tib %>% mutate(thresh = case_when(no_zeros <= nz_genes ~ "TRUE", TRUE ~ "FALSE")) #sort genes based on number of zeroes
#     #plot scatter plot
#     cols <- c("FALSE" = "red", "TRUE" = "black")
#     scatter <- ggplot(plot_tib, aes(Median, no_zeros))+
#         geom_point(aes(color=thresh), alpha=0.75)+
#         scale_color_manual(values = cols)+
#         scale_x_log10()+
#         labs(title = 'Plot of Median vs Number of Non-Zero genes', subtitle = "Genes filtered out are in red. X-axis is log scaled.")+
#         theme_bw()+
#         ylab('Number of samples with zero count')+
#         theme(legend.position = 'bottom')
#     return(scatter)
# }

med_vs_nz<- function(dataf, slider_var, slider_num){
    stats_df <- dataf
    stats_df$variance <- apply(dataf[,-1], 1, var)
    stats_df <- arrange(stats_df, by=variance)
    stats_df$rank <- rank(stats_df$variance)
    stats_df$nonzero <- rowSums(dataf[,-1]!=0)
    stats_df$median <- apply(dataf[,-1], 1, median)
    filtered_status <- stats_df %>%
        mutate(filter_status=ifelse(rank >= nrow(stats_df)*(slider_var/100) & nonzero >= slider_num, 'pass filter', 'filtered out'))
    # p1 <- ggplot(filtered_status, aes(x=variance, y=median, color=filter_status)) +
    #     geom_point() +
    #     scale_y_log10() +
    #     scale_x_log10() +
    #     labs(x='Percentile of Variance', y='Median Number of Counts', color='Filter Status', main='Median Counts vs. Percentile Variance') +
    #     theme_dark() +
    #     scale_color_brewer(palette="Greens")

    p2 <- ggplot(filtered_status, aes(x= nonzero , y=median, color=filter_status)) +
        geom_point() +
        scale_y_log10() +
        labs(x='Number of Non-Zero Samples', y='Median Number of Counts', color='Filter Status', main='Median Counts vs. Number of Non-Zero Samples') +
        theme_linedraw() +
        scale_color_brewer(palette="Set1")+
        theme(legend.position = 'bottom')
    # grid.arrange(p1, p2, nrow = 1)
    return(p2)
}




## Function 4 --> Heatmap
plot_heatmap <- function(filtered_df,log_trans) {
    # remove the "gene" column and transpose the dataframe to
    # display genes on the y-axis and the samples on the x-axis of the heatmap
    heatmap_df <- as.matrix(filtered_df[, -1])

    # transform counts to log10 values if log_trans is TRUE
    if(log_trans) {
        heatmap_df <- log10(heatmap_df + 1)
    } else {
        heatmap_df <- heatmap_df
    }
    # create heatmap with color bar
    return(heatmap.2(heatmap_df,
              col=brewer.pal(9, "YlOrRd"),
              main = "Clustered Heatmap of Filtered Counts"))
}

## Function 5 --> PCA
performPCA <- function(fdata){
    object <- PCAtools::pca(fdata[,-1])
    return(object)
}
# use object to get the options

## Function 6 --> plot PCA
pcaplot <- function(pca_obj, pca, pcb, loadings = FALSE ){
    plot <- biplot(pca_obj, x = pca, y=pcb, showLoadings = loadings)
    return(plot)
}


# Tab 3 --> DESeq Exploration ----

load_data <- function(path){
    file <- read_tsv(file = path) %>%
        dplyr::rename("genes" = "...1")
    return(file)
}

draw_table <- function(dataf, slider) {
    filtered_table <- dataf %>% dplyr::filter(padj <= 1*10^slider) %>%
        drop_na()
    # filtered_table$padj <- sprintf("%.2f", filtered_table$padj)
    filtered_table$padj <- format(filtered_table$padj, digits = 6)
    filtered_table$pvalue <- format(filtered_table$pvalue, digits = 6)
    return(filtered_table)
}

volcano_plot <-
    function(dataf, x_name, y_name, slider, color1, color2) {
        dataf <- dataf %>% drop_na()
        plot <- ggplot(data = dataf, aes(x = !!sym(x_name), y = -log10(!!sym(y_name))))+
            geom_point(aes(color = if_else(!!sym(y_name) < 1*10^slider, TRUE,FALSE)))+
            labs(title="DESeq Results", x = glue('{x_name}'), y = glue('-log10({y_name})')) +
            scale_color_manual(name = glue("{x_name} < 1*10^{slider}"),
                               values = c("FALSE" = color1,
                                          "TRUE" = color2),
                               labels = c("FALSE", "TRUE"))+
            theme_linedraw()+
            theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")

        return(plot)
    }

# Tab 4 --> Expression Exploration ----

geneexpression <- function(sampledata = NULL, countsdata = NULL, gene = NULL, variable = NULL) {
    v_gene <- countsdata %>% filter(genes == gene) %>% select(-1) %>% t() %>% tibble(count = as.numeric(.)) %>% select(-1) %>% pull(count)
    v_variable <- sampledata %>% pull(variable)
    disease_status <- sampledata %>% pull(diagnosis)
    sample <- sampledata %>% pull(title)
    data_table <- tibble(counts = v_gene, variable = v_variable, status = disease_status, sample = sample)

    p <- plot_ly(data_table, x = ~variable, y = ~counts, color = ~status, text = ~sample) %>%
        add_markers() %>%
        layout(xaxis = list(title = variable), yaxis = list(title = "Counts"), showlegend = TRUE)

    return(p)
}


