# Function to recursively clean the list
remove_column_from_list <- function(complex_list, column_to_remove = "sim.1") {
  lapply(complex_list, function(element) {
    if (is.data.frame(element)) {
      # Check if the column exists and remove it
      if (column_to_remove %in% colnames(element)) {
        element <- element[, !colnames(element) %in% column_to_remove, drop = FALSE]
      }
      return(element)
    } else if (is.list(element)) {
      # Recursively process sub-lists
      return(remove_column_from_list(element, column_to_remove))
    } else {
      # Return the element as-is if not a data frame or list
      return(element)
    }
  })
}

## Example usage
#cleaned_list <- remove_column_from_list(complex_list, column_to_remove = "sim.1")
