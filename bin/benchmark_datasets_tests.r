setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
# Sample data
set.seed(123)
options(repos = list(CRAN="http://cran.rstudio.com/"))

source("utils.r")

# Load necessary libraries
library("pacman")
p_load("cluster" ,"factoextra", "clValid", "clusterCrit", "entropy", "FCPS", "clusterSim", "fpc", "dplyr", "crayon", "gridExtra")

# Function to compute Gap statistic
gap_statistic <- function(data, max_k, B = 100) {
  # Compute the total within-cluster sum of squares for k-means
  wss <- numeric(max_k)
  for (k in 1:max_k) {
    kmeans_result <- kmeans(data, centers = k, nstart = 25)
    wss[k] <- kmeans_result$tot.withinss
  }

  # Generate reference datasets
  reference_wss <- numeric(max_k)
  for (b in 1:B) {
    # Create a random uniform dataset
    random_data <- matrix(runif(nrow(data) * ncol(data)), ncol = ncol(data))
    for (k in 1:max_k) {
      kmeans_result <- kmeans(random_data, centers = k, nstart = 25)
      reference_wss[k] <- reference_wss[k] + kmeans_result$tot.withinss
    }
  }
  reference_wss <- reference_wss / B  # Average over B simulations

  # Calculate Gap statistic
  gap <- log(reference_wss) - log(wss)
  return(gap[max_k])
}




# applyKMeansAndCalculateMetrics
applyKMeansAndCalculateMetrics <- function(data, dataset_name, topK) {

    startK <- 2
    endK <- 10

    cat("Dataset considered: ", dataset_name, "\n", sep="")

    silhouette_array <- c()
    chi_array <- c()
    rev_dbi_array <- c()
    dunn_array <- c()
    rev_entropy_array <- c()
    gap_array <- c()

    for (k in seq(startK, endK)) {

        NUMBER_OF_CLUSTERS <- k

        # Apply k-means clustering
        kmeans_result <- kmeans(data, centers=NUMBER_OF_CLUSTERS)

        # Compute Silhouette index
        silhouette_score <- silhouette(kmeans_result$cluster, dist(data))
        mean_silhouette <- mean(silhouette_score[, ncol(silhouette_score)])

        # Compute Calinski-Harabasz Index (CHI)
        chi_index <- cluster.stats(dist(data), kmeans_result$cluster)$ch

        # Compute Davies-Bouldin Index (DBI)
        db_index <- clusterSim::index.DB(data, kmeans_result$cluster)$"DB"
        reverse_dbi <- 1/db_index

        # Compute Dunn index
        dunn_index <- dunn(dist(data), kmeans_result$cluster)

        # Compute Shannon entropy
        # Here, we calculate the proportion of points in each cluster
        proportions <- table(kmeans_result$cluster) / length(kmeans_result$cluster)
        shannon_entropy <- -sum(proportions * log(proportions))
        reverse_entropy <- 1/shannon_entropy

        # gap_stat <- gap_statistic(data, NUMBER_OF_CLUSTERS)
        gap_stat <- clusGap(data, FUNcluster = kmeans, nstart = 20, K.max = 2, B = 60)
        gap <- gap_stat$Tab[gap_stat$Tab %>% nrow(),3]

        topKstring <- ""
        if(k == topK) topKstring <- "*****"

        if(k == 2) cat("(topK = ", topK,  ") Sil \t CHI \t   rDBI    rSE \t    Dunn   Gap \n",  sep="")
        if(k == topK) cat(green("(k = ", NUMBER_OF_CLUSTERS, ") ", dec_five(mean_silhouette), " ", dec_five(chi_index), " ", dec_five(reverse_dbi), " ", dec_five(reverse_entropy), " ", dec_five(dunn_index), " ", dec_five(gap), " ", topKstring, "\n", sep=""))
        else cat("(k = ", NUMBER_OF_CLUSTERS, ") ", dec_five(mean_silhouette), " ", dec_five(chi_index), " ", dec_five(reverse_dbi), " ", dec_five(reverse_entropy), " ", dec_five(dunn_index), " ", dec_five(gap), " ", topKstring, "\n", sep="")

        silhouette_array <- c(silhouette_array, mean_silhouette)
        chi_array <- c(chi_array, chi_index)
        rev_dbi_array <- c(rev_dbi_array, reverse_dbi)
        dunn_array <- c(dunn_array, dunn_index)
        rev_entropy_array <- c(rev_entropy_array, reverse_entropy)
        gap_array <- c(gap_array, gap)

#         # Output the results
#         cat("(k = ", NUMBER_OF_CLUSTERS, ") Silhouette coefficient: \t", mean_silhouette, "\n", sep="")
#         cat("(k = ", NUMBER_OF_CLUSTERS, ") Calinski-Harabasz index: \t", chi_index, "\n", sep="")
#         cat("(k = ", NUMBER_OF_CLUSTERS, ") reverse Davies-Bouldin index: \t", (1/db_index), "\n", sep="")
#         cat("(k = ", NUMBER_OF_CLUSTERS, ") reverse Shannon entropy: \t", (1/shannon_entropy), "\n", sep="")
#         cat("(k = ", NUMBER_OF_CLUSTERS, ") Dunn index: \t\t\t", dunn_index, "\n", sep="")
#         cat("(k = ", NUMBER_OF_CLUSTERS, ") Gap statistic: \t\t\t", gap, "\n\n", sep="")

    }

    output_dataframe <- data.frame(k=seq(startK, endK), silhouette_array, chi_array, rev_dbi_array, dunn_array, rev_entropy_array, gap_array)

    return(output_dataframe)

}

# plotResultsScatterplot
plotResultsScatterplot <- function(this_dataframe, dataset_name, thisTopK, saveSinglePlot, saveArrangedPlot) {

    xlim_start <- 0
    xlim_end <- nrow(this_dataframe) + 1
    dist_label <-  0.9
    this_linewidth <- 0.5
    this_vjust <- 1.5
    this_vertical_label_y <- 0.2

    # Silhouette
    max_silhouette_y <- max(this_dataframe$silhouette_array)
    max_silhouette_x <- this_dataframe$k[which.max(this_dataframe$silhouette_array)]

    silhouette_plot <-  ggplot(data = this_dataframe, aes(x = k, y = silhouette_array)) +
    geom_point() + geom_line() +
    geom_vline(xintercept = max_silhouette_x, color = "blue", linetype = "longdash", linewidth = this_linewidth-0.2) +  # Add vertical line
    geom_vline(xintercept = thisTopK, color = "green", linetype = "dashed", linewidth = this_linewidth) +  # Add vertical line
    annotate("text", x = (thisTopK+dist_label), y = this_vertical_label_y, label = "\"correct\"", color = "green")  +  # Add label
    annotate("text", x = (max_silhouette_x+0.4), y = this_vertical_label_y+0.2, label = "top", color = "blue")  +  # Add label
    labs(title = paste0(dataset_name, " dataset"), x = "k (number of clusters)", y = "Silhouette score") +
    theme_minimal() + ylim(0,1)   +  scale_x_continuous(breaks = seq(xlim_start, xlim_end, by = 1)) + # Set x-axis breaks to integer values


    if(saveSinglePlot) {

          fileName <- paste0("../results/", dataset_name, "_silhouette.pdf")
          ggsave(silhouette_plot, file=fileName)
          cat("saved file ", fileName, "\n", sep="")
    }


    # Calinski-Harabas index
    max_chi_y <- max(this_dataframe$chi_array)
    max_chi_x <- this_dataframe$k[which.max(this_dataframe$chi_array)]

    CHI_plot <-  ggplot(data = this_dataframe, aes(x = k, y = chi_array)) +
    geom_point() + geom_line() +
    geom_vline(xintercept = max_chi_x, color = "blue", linetype = "longdash", linewidth = this_linewidth-0.2) +  # Add vertical line
    geom_vline(xintercept = thisTopK, color = "green", linetype = "dashed", linewidth = this_linewidth) +  # Add vertical line
    labs(title = paste0(dataset_name, " dataset"), x = "k (number of clusters)", y = "Calinski-Harabasz index") +
    theme_minimal()  +  scale_x_continuous(breaks = seq(xlim_start, xlim_end, by = 1)) + # Set x-axis breaks to integer values

    if(saveSinglePlot) {

          fileName <- paste0("../results/", dataset_name, "_CHI_.pdf")
          ggsave(CHI_plot, file=fileName)
          cat("saved file ", fileName, "\n", sep="")
    }

    # Davies-Bouldin index
    max_dbi_y <- max(this_dataframe$rev_dbi_array)
    max_dbi_x <- this_dataframe$k[which.max(this_dataframe$rev_dbi_array)]

    DBI_plot <-  ggplot(data = this_dataframe, aes(x = k, y = rev_dbi_array)) +
    geom_point() + geom_line() +
    geom_vline(xintercept = max_dbi_x, color = "blue", linetype = "longdash", linewidth = this_linewidth-0.2) +  # Add vertical line
    geom_vline(xintercept = thisTopK, color = "green", linetype = "dashed", linewidth = this_linewidth) +  # Add vertical line
    labs(title = paste0(dataset_name, " dataset"), x = "k (number of clusters)", y = "reverse Davies-Bouldin index") +
    theme_minimal() +  scale_x_continuous(breaks = seq(xlim_start, xlim_end, by = 1)) + # Set x-axis breaks to integer values

    if(saveSinglePlot) {

          fileName <- paste0("../results/", dataset_name, "_DBI_.pdf")
          ggsave(DBI_plot, file=fileName)
          cat("saved file ", fileName, "\n", sep="")
    }

    # Dunn index
    max_dunn_y <- max(this_dataframe$dunn_array)
    max_dunn_x <- this_dataframe$k[which.max(this_dataframe$dunn_array)]

    Dunn_plot <-  ggplot(data = this_dataframe, aes(x = k, y = dunn_array)) +
    geom_point() + geom_line() +
    geom_vline(xintercept = max_dunn_x, color = "blue", linetype = "longdash", linewidth = this_linewidth-0.2) +  # Add vertical line
    geom_vline(xintercept = thisTopK, color = "green", linetype = "dashed", linewidth = this_linewidth) +  # Add vertical line
    labs(title = paste0(dataset_name, " dataset"), x = "k (number of clusters)", y = "Dunn index") +
    theme_minimal() +  scale_x_continuous(breaks = seq(xlim_start, xlim_end, by = 1)) + # Set x-axis breaks to integer values

    if(saveSinglePlot) {

          fileName <- paste0("../results/", dataset_name, "_Dunn_.pdf")
          ggsave(Dunn_plot, file=fileName)
          cat("saved file ", fileName, "\n", sep="")
    }


    # Shannon entropy
    max_entropy_y <- max(this_dataframe$rev_entropy_array)
    max_entropy_x <- this_dataframe$k[which.max(this_dataframe$rev_entropy_array)]

    entropy_plot <-  ggplot(data = this_dataframe, aes(x = k, y = rev_entropy_array)) +
    geom_point() + geom_line() +
    geom_vline(xintercept = max_entropy_x, color = "blue", linetype = "longdash", linewidth = this_linewidth-0.2) +  # Add vertical line
    geom_vline(xintercept = thisTopK, color = "green", linetype = "dashed", linewidth = this_linewidth) +  # Add vertical line
    labs(title = paste0(dataset_name, " dataset"), x = "k (number of clusters)", y = "reverse Shannon entroy") +
    theme_minimal() +  scale_x_continuous(breaks = seq(xlim_start, xlim_end, by = 1)) + # Set x-axis breaks to integer values

    if(saveSinglePlot) {

          fileName <- paste0("../results/", dataset_name, "_entropy_.pdf")
          ggsave(entropy_plot, file=fileName)
          cat("saved file ", fileName, "\n", sep="")
    }

    # Gap
    max_gap_y <- max(this_dataframe$gap_array)
    max_gap_x <- this_dataframe$k[which.max(this_dataframe$gap_array)]

    gap_plot <-  ggplot(data = this_dataframe, aes(x = k, y = gap_array)) +
    geom_point() + geom_line() +
    geom_vline(xintercept = max_gap_x, color = "blue", linetype = "longdash", linewidth = this_linewidth-0.2) +  # Add vertical line
    geom_vline(xintercept = thisTopK, color = "green", linetype = "dashed", linewidth = this_linewidth) +  # Add vertical line
    labs(title = paste0(dataset_name, " dataset"), x = "k (number of clusters)", y = "Gap statistic") +
    theme_minimal() +  scale_x_continuous(breaks = seq(xlim_start, xlim_end, by = 1)) + # Set x-axis breaks to integer values

    if(saveSinglePlot) {

          fileName <- paste0("../results/", dataset_name, "_gap_.pdf")
          ggsave(gap_plot, file=fileName)
          cat("saved file ", fileName, "\n", sep="")
    }


    arranged_plot <- grid.arrange(silhouette_plot, CHI_plot, DBI_plot, Dunn_plot, entropy_plot, gap_plot, ncol = 2)
    if(saveArrangedPlot) {

          fileName <- paste0("../results/", dataset_name, "_arranged_plot.pdf")
          ggsave(arranged_plot, file=fileName)
          cat("saved file ", fileName, "\n", sep="")
    }


}

thisSavePlot <- FALSE
thisSaveArrangedPlot <- TRUE

cat(" = = = =  = = = =  = = = =  = = = = \n")
Ktop <- 2
this_dataset_name <- "TwoDiamonds"
cat("dataset: ", this_dataset_name, ": the \"correct\" number of clusters here is 2\n")
all_results <- applyKMeansAndCalculateMetrics(TwoDiamonds$"Data", this_dataset_name, Ktop)
plotResultsScatterplot(all_results, this_dataset_name, Ktop,  thisSavePlot, thisSaveArrangedPlot)

cat(" = = = =  = = = =  = = = =  = = = = \n")
Ktop <- 2
this_dataset_name <- "WingNut"
cat("dataset: ", this_dataset_name, ": the \"correct\" number of clusters here is 2\n")
all_results <- applyKMeansAndCalculateMetrics(WingNut$"Data", this_dataset_name, Ktop)
plotResultsScatterplot(all_results, this_dataset_name, Ktop,  thisSavePlot, thisSaveArrangedPlot)

cat(" = = = =  = = = =  = = = =  = = = = \n")
Ktop <- 4
this_dataset_name <- "Tetra"
cat("dataset: ", this_dataset_name, ": the \"correct\" number of clusters here is 2\n")
all_results <- applyKMeansAndCalculateMetrics(Tetra$"Data", this_dataset_name, Ktop)
plotResultsScatterplot(all_results, this_dataset_name, Ktop,  thisSavePlot, thisSaveArrangedPlot)

cat(" = = = =  = = = =  = = = =  = = = = \n")
Ktop <- 3
this_dataset_name <- "Lsun3D"
cat("dataset: ", this_dataset_name, ": the \"correct\" number of clusters here is 2\n")
all_results <- applyKMeansAndCalculateMetrics(Lsun3D$"Data", this_dataset_name, Ktop)
plotResultsScatterplot(all_results, this_dataset_name, Ktop,  thisSavePlot, thisSaveArrangedPlot)

cat(" = = = =  = = = =  = = = =  = = = = \n")
Ktop <- 7
this_dataset_name <- "Hepta"
cat("dataset: ", this_dataset_name, ": the \"correct\" number of clusters here is 2\n")
all_results <- applyKMeansAndCalculateMetrics(Hepta$"Data", this_dataset_name, Ktop)
plotResultsScatterplot(all_results, this_dataset_name, Ktop,  thisSavePlot, thisSaveArrangedPlot)

