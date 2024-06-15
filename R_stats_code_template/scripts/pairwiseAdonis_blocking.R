perform_pairwise_adonis_blocking <- function(dist_object, data, grouping_var, blocking_var) {
  # Get unique levels of the blocking factor
  blocking_levels <- unique(data[[blocking_var]])
  
  # Initialize a list to store the pairwise comparison results
  pairwise_results <- list()
  
  # Loop through each level of the blocking factor and perform pairwise comparisons
  for (level in blocking_levels) {
    # Obtain the indices corresponding to the current level of the blocking factor
    indices <- which(data[[blocking_var]] == level)
    
    # Obtain the names of observations for the current level
    obs_names <- rownames(data)[indices]
    
    # Subset the data for the current level
    subset_data <- data[indices, ]
    
    # Convert dist object to a square matrix
    square_matrix <- as.matrix(dist_object)
    
    # Subset square matrix
    subset_matrix <- square_matrix[obs_names, obs_names]
    
    # Convert matrix back to dist
    subset_dist_object <- as.dist(subset_matrix)
    
    # Construct the formula for pairwise comparison
    formula <- paste("subset_dist_object ~", grouping_var)
    
    # Perform pairwise comparison for the subset
    pairwise_result <- pairwise.adonis2(
      as.formula(formula),  # subset dist and specify the grouping variable
      data = subset_data,  
      nperm = 9999,
      p.adjust.methods = "BH" # Benjamini-Hochberg correction
    )
    
    # Store the pairwise result
    pairwise_results[[as.character(level)]] <- pairwise_result
  }
  
  # Return the pairwise comparison results
  return(pairwise_results)
}

