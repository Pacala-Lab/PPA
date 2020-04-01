## Utility Functions
## anything that doesn't really belong somewhere else, or is generally useful

## function to extract a single column from a simulation data matrix
extract_column <- function(colname, data) {

  data_matrix <- data[[colname]]

  if(!is.matrix(data_matrix)) {
    print("Error: column matrix must be of class matrix")
    return()
  }

  column <- as.numeric(data_matrix)
  return(column)

}


## utility function to calculate the number of species of a given crown_class at a given timepoint
spp_populations <- function(data, spp, i, crown_class) {
  if(missing(crown_class)) {
    out <- sum(data[["spp"]][,i] == spp)
  }
 else {
    out <- sum(data[["crown_class"]][, i] == crown_class & data[["spp"]][, i] == spp)
  }
  return(out)
}
