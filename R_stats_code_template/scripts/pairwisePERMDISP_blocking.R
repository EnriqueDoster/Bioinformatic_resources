perform_pairwise_PERMDISP_blocking <- function(dist_object, data, grouping_var, blocking_var) {
  # Get unique levels of the blocking factor
  blocking_levels <- unique(data[[blocking_var]])
  
  # Initialize a list to store the pairwise comparison results
  permdisp_results <- list()
  
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
    
    # Perform betadisper for the subset
    disper_result <- betadisper(
      subset_dist_object, subset_data[[grouping_var]]  # subset dist and specify the grouping variable
    )
    
    #Permute p-values for betadisper
    permdisp_result <- permutest(
      disper_result,
      permutations = 9999,
      pairwise=T
    )
    
    # Store the pairwise result
    permdisp_results[[as.character(level)]] <- permdisp_result
  }
  
  # Return the pairwise comparison results
  return(permdisp_results)
}

