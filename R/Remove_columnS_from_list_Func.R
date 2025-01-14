remove_columns_from_list <- function(complex_list, ...) {
  # Capture column names as a character vector
  columns_to_remove <- c(...)
  
  lapply(complex_list, function(element) {
    if (is.data.frame(element)) {
      # Check if any of the specified columns exist and remove them
      cols_to_keep <- setdiff(colnames(element), columns_to_remove)
      element <- element[, cols_to_keep, drop = FALSE]
      return(element)
    } else if (is.list(element)) {
      # Recursively process sub-lists
      return(remove_columns_from_list(element, ...))
    } else {
      # Return the element as-is if not a data frame or list
      return(element)
    }
  })
}

# Example usage:
# remove_columns_from_list(complex_list, "sim.1", "extra")
