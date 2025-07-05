setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
# Sample data
set.seed(123)
options(repos = list(CRAN="http://cran.rstudio.com/"))

source("utils.r")

# Load necessary libraries
library("pacman")
p_load("cluster" ,"factoextra", "clValid", "clusterCrit", "entropy", "FCPS", "clusterSim", "fpc", "dplyr", "crayon")

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

        # Compute Dunn index
        dunn_index <- dunn(dist(data), kmeans_result$cluster)

        # Compute Shannon entropy
        # Here, we calculate the proportion of points in each cluster
        proportions <- table(kmeans_result$cluster) / length(kmeans_result$cluster)
        shannon_entropy <- -sum(proportions * log(proportions))

        # gap_stat <- gap_statistic(data, NUMBER_OF_CLUSTERS)
        gap_stat <- clusGap(data, FUNcluster = kmeans, nstart = 20, K.max = 2, B = 60)
        gap <- gap_stat$Tab[gap_stat$Tab %>% nrow(),3]

        topKstring <- ""
        if(k == topK) topKstring <- "*****"

        if(k == 2) cat("(topK = ", topK,  ") Sil \t CHI \t   rDBI    rSE \t    Dunn   Gap \n",  sep="")
        if(k == topK) cat(green("(k = ", NUMBER_OF_CLUSTERS, ") ", dec_five(mean_silhouette), " ", dec_five(chi_index), " ", dec_five(1/db_index), " ", dec_five(1/shannon_entropy), " ", dec_five(dunn_index), " ", dec_five(gap), " ", topKstring, "\n", sep=""))
        else cat("(k = ", NUMBER_OF_CLUSTERS, ") ", dec_five(mean_silhouette), " ", dec_five(chi_index), " ", dec_five(1/db_index), " ", dec_five(1/shannon_entropy), " ", dec_five(dunn_index), " ", dec_five(gap), " ", topKstring, "\n", sep="")

#         # Output the results
#         cat("(k = ", NUMBER_OF_CLUSTERS, ") Silhouette coefficient: \t", mean_silhouette, "\n", sep="")
#         cat("(k = ", NUMBER_OF_CLUSTERS, ") Calinski-Harabasz index: \t", chi_index, "\n", sep="")
#         cat("(k = ", NUMBER_OF_CLUSTERS, ") reverse Davies-Bouldin index: \t", (1/db_index), "\n", sep="")
#         cat("(k = ", NUMBER_OF_CLUSTERS, ") reverse Shannon entropy: \t", (1/shannon_entropy), "\n", sep="")
#         cat("(k = ", NUMBER_OF_CLUSTERS, ") Dunn index: \t\t\t", dunn_index, "\n", sep="")
#         cat("(k = ", NUMBER_OF_CLUSTERS, ") Gap statistic: \t\t\t", gap, "\n\n", sep="")

    }

}

cat(" = = = =  = = = =  = = = =  = = = = \n")
Ktop <- 2
cat("dataset: TwoDiamonds: the \"correct\" number of clusters here is 2\n")
applyKMeansAndCalculateMetrics(TwoDiamonds$"Data", "TwoDiamonds", Ktop)

cat(" = = = =  = = = =  = = = =  = = = = \n")
Ktop <- 2
cat("dataset: WingNut: the \"correct\" number of clusters here is 2\n")
applyKMeansAndCalculateMetrics(WingNut$"Data", "WingNut", Ktop)

cat(" = = = =  = = = =  = = = =  = = = = \n")
Ktop <- 4
cat("dataset: Tetra: the \"correct\" number of clusters here is 4\n")
applyKMeansAndCalculateMetrics(Tetra$"Data", "Tetra", Ktop)

cat(" = = = =  = = = =  = = = =  = = = = \n")
Ktop <- 3
cat("dataset: Lsun3D: the \"correct\" number of clusters here is 3\n")
applyKMeansAndCalculateMetrics(Lsun3D$"Data", "Lsun3D", Ktop)

cat(" = = = =  = = = =  = = = =  = = = = \n")
Ktop <- 7
cat("dataset: Hepta: the \"correct\" number of clusters here is 7\n")
applyKMeansAndCalculateMetrics(Hepta$"Data", "Hepta", Ktop)


