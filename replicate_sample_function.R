# Function to subset random 4 biological replicates of a group
library(stringr)

replicate_sample <- function(samples, group, size = 4) {
  gp_index <- grep(paste0("^", group), samples)
  sample_index <- sample(gp_index, size)
  return(sample_index)
}

replicate_sample_all <- function(samples, size = 4) {
  gps <- unique(str_remove(samples, "[0-9]+$"))
  full_index <- lapply(gps, FUN = function(x) replicate_sample(samples, x))
  full_index <- unlist(full_index)
  return(full_index)
}


#### If time, come back to playing around with more general function that can deal with different input sample sizes.
#If replace - FALSE, return only the same samples as original, or if replace = TRUE, use replacement
replicate_sample_2 <- function(samples, group, size = 4, replace = FALSE) {
  gp_index <- grep(paste0("^", group), samples)
  if (length(gp_index) < size) {
    print(paste0("Warning, group ", group, " has less than ", size, "replicates."))
  }
  sample_index <- sample(gp_index, size)
  return(sample_index)
}

