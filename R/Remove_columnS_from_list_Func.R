remove_columnS_from_list <- function(complex_list, columns_to_remove = c("col1", "col2")) {
  lapply(complex_list, function(element) {
    if (is.data.frame(element)) {
      # Check if any of the specified columns exist and remove them
      cols_to_keep <- setdiff(colnames(element), columns_to_remove)
      element <- element[, cols_to_keep, drop = FALSE]
      return(element)
    } else if (is.list(element)) {
      # Recursively process sub-lists
      return(remove_columns_from_list(element, columns_to_remove))
    } else {
      # Return the element as-is if not a data frame or list
      return(element)
    }
  })
}

## Example usage:
#columns_to_remove <- c("sim.1", "extra")
#cleaned_list <- remove_columns_from_list(complex_list, columns_to_remove)