

#' Model Light Competition
#'
#' This function simulates a stand of trees using the perfect plasticity approximation for light competition.
#'
#'@param allometry an object of class: allometry created using the create_allometry method
#'@param spp a vector of species names to simulate
#'@param max_t the last timestep in the simulation
#'@param timestep the length of each timestep
#'@param dnot the initial diameter for each tree
#'@param plot.area the total area of the plot to be simulated
#'@param n_start a named list of the starting abundances. Names must correspond to the spp names specified in spp.
#'@param full.allom TRUE/FALSE specifying whether you would like the model output to include height and structural biomass
#'@param full.data TRUE/FALSE specifying whether you would like the model output to include data from each timestep or just the last timestep
#'
#'@return a list object of class "model_light" which includes simulated growth and allometric data


## simulator for light competition model.
model_light <- function(allometry = create_allometry.model_light(n.spp = 1),
                           spp = 1,
                           max_t,
                           timestep,
                           dnot = rnorm(n = n_start, mean = 0.0004, sd = 0.00005),
                           plot.area,
                           n_start,
                           full.allom = TRUE,
                           full.data = TRUE) {

  ## stop on error
  stopifnot(ncol(allometry) >= length(spp),
            is.numeric(max_t),
            is.numeric(timestep),
            is.numeric(dnot),
            is.numeric(plot.area),
            timestep <= max_t,
            is.logical(full.allom),
            !is.null(names(n_start)))

  ## ensure spp name is character not numeric
  spp <- as.character(spp)

  ## initialize data list
  data <- list(spp = matrix(1, nrow = sum(unlist(n_start)), ncol = max_t/timestep + 1),
               diameter = matrix(1, nrow = sum(unlist(n_start)), ncol = max_t/timestep+1),
               crown_class = matrix(1, nrow = sum(unlist(n_start)), ncol = max_t/timestep + 1),
               diameter_growth = matrix(1, nrow = sum(unlist(n_start)), ncol = max_t/timestep + 1))

  data[["diameter"]][,1] <- dnot
  data[["spp"]][,1] <- c(rep(spp[1], times = n_start[spp[1]]), rep(spp[2], times = n_start[spp[2]]))

  A <- calc.A(allometry = allometry, n.crown.class = 2, spp = spp)

  ## loop over days in simulation
  for (i in seq(1, max_t, by = timestep)) {

    ## ------------------------- KILL PLANTS ----------------------------- ##

    ## pre-calculate # of each spp in each crown class
    spp_pops <- list("1" = sapply(FUN = spp_populations, spp, i = i, crown_class = 1, data = data, simplify = "vector"),
                     "2" = sapply(FUN = spp_populations, spp, i = i, crown_class = 1, data = data, simplify = "vector"))

    ## create lists of death outcomes based on each spp's individual death probability for canopy and understor
    mu_c <- mapply(rbinom, spp_pops[["1"]], matrix(allometry["mu_c", 2:length(spp)], nrow = 1), size = 1)

    if(length(spp_pops[["2"]]) > 0) {
      mu_u <- mapply(rbinom, spp_pops[["2"]], matrix(allometry["mu_u", 2:length(spp)], nrow = 1), size = 1)
      mu <- list(mu_c)
    }
    else {
        mu <- list(mu_c, mu_u) ## combine into one list
    }
    ## reorder data to easily index and align death vector and data matrices

    data <- order_simulation(data = data, by = c("crown_class", "spp"), decreasing = TRUE, i = i)

    ## kill trees by removing all trees with mu == 0 from data. In future could include tracking of dead trees, or perhaps mortality rates?
    ## Though not particularly interesting with "density independent" death process
    data[["spp"]] <- data[["spp"]][unlist(mu) == 0, ]
    data[["diameter"]] <- data[["diameter"]][unlist(mu) == 0, ]
    data[["crown_class"]] <- data[["crown_class"]][unlist(mu) == 0, ]
    data[["diameter_growth"]] <- data[["diameter_growth"]][unlist(mu) == 0, ]

    ## ------------------------- ASSIGN NEW SPP ID ------------------------##

    data[["spp"]][, i + 1] <- data[["spp"]][, i]

    ## ------------------------- GROW PLANTS ------------------------------ ##

    ## grow plants
    ## calculate diameter_growth for overstory and understory
    ## currently goint to put this into a loop. There should be a more efficient way to do this using apply and/or lookup table, but my head hurts
    for (j in 1:length(spp)) {

      if (sum(data[["crown_class"]][, i] == 1 & data[["spp"]][, i] == spp[j]) > 0) {

        data[["diameter_growth"]][data[["crown_class"]][, i] == 1 & data[["spp"]][, i] == spp[j], i] <- sapply(data[["diameter"]][data[["crown_class"]][, i] == 1 & data[["spp"]][, i] == spp[j], i],
                                                                                                              calc.dd, allometry = allometry,
                                                                                                              l = allometry["l_c", spp[j]],
                                                                                                              r <- allometry["r_c", spp[j]],
                                                                                                              spp = spp[j], A = A[spp[j], "A_c"],
                                                                                                              simplify = "vector")

      }

      if (sum(data[["crown_class"]][, i] == 2 & data[["spp"]][, i] == spp[j]) > 0) {

        if (class(data[["diameter_growth"]]) !=  "matrix") print(c(i, j))

        data[["diameter_growth"]][data[["crown_class"]][, i] == 2 & data[["spp"]][, i] == spp[j], i] <- sapply(data[["diameter"]][data[["crown_class"]][, i] == 2 & data[["spp"]][, i] == spp[j], i],
                                                                                                              calc.dd, allometry = allometry,
                                                                                                              l = allometry["l_u", spp[j]],
                                                                                                              r <- allometry["r_u", spp[j]],
                                                                                                              spp = spp[j], A = A[spp[j], "A_u"],
                                                                                                              simplify = "vector")
      
      }

    }

    ## this is just to improve interpretability for summary functions, likely there is a cleaner way to do this?
    if (i == max_t) {

      data[["diameter_growth"]][, i+1] <- NA

    }

    ## calculate diameter for next timestep
    data[["diameter"]][, i+1] <- data[["diameter"]][, i] + data[["diameter_growth"]][, i]

    ## ------------------- CANOPY CLASS REASSIGNMENT ----------------- ##

    CA <- sum(mapply(FUN = calc.CA,
                     alpha_w = allometry["alpha_w", spp],
                     gamma = allometry["gamma", spp],
                     spp = spp,
                     MoreArgs = list(diam_data = data[["diameter"]],
                                     spp_data = data[["spp"]],
                                     i = i),
                     SIMPLIFY = "vector"))

    ## check if plants need to added to the understory
    if (CA <= plot.area) {

      data[["crown_class"]][, i+1] <- 1

    }

    else {
      ## order data from largest to smallest
      data <- order_simulation(data = data, by = c("diameter"), decreasing = TRUE, i = i)

      ## create vector of cumulative CA sums, used to determine which trees to add/remove from canopy
      ca.sum <- vector(length = nrow(data[["diameter"]]), mode = "numeric")

      for (j in 1:nrow(data[["diameter"]])) {

        if (j == 1) {

          ca.sum[j] <- allometry["alpha_w", data[["spp"]][j, i+1]] * data[["diameter"]][j, i+1]^allometry["gamma", data[["spp"]][j, i+1]]

        }

        else {

          ca.sum[j] <- sum((allometry["alpha_w", data[["spp"]][j, i+1]] * data[["diameter"]][j, i+1]^allometry["gamma", data[["spp"]][j, i+1]]), ca.sum[j - 1])

        }

        ## if gap in canopy trees, add understory trees to canopy and vis versa
        data[["crown_class"]][ca.sum <= plot.area, i+1] <- 1
        data[["crown_class"]][ca.sum > plot.area, i+1] <- 2
      }
    }

  }

  ## assign informative values to list items
  data[["max_t"]] <- max_t
  data[["timestep"]] <- timestep
  data[["plot_area"]] <- plot.area

  ## assign class of return object to class simulation
  class(data) <- "model_light"

  ## if user specifies full allometry, calculate and add to the data
  if (full.allom) {

    data <- calculate_allometry.model_light(data, allometry = allometry, spp = spp)

  }
  
  ## if user specifies that full data should be reporeted in single matrix, initialize and populate matrix
  if (full.data) {

    colnames <- names(data)[!(names(data) %in% c("avgs", "max_t", "timestep", "plot_area"))]
    columns <- lapply(colnames, extract_column, data = data)
    full_data <- matrix(unlist(columns), ncol = length(colnames))
    time_column <- rep(seq(1, (max_t+timestep), by = timestep), each = nrow(data[["diameter"]]))
    ind_column <- rep(seq(1, nrow(data[["diameter"]]), by = 1), times = (max_t/timestep)+1)
    full_data <- cbind(time_column, ind_column, full_data)
    colnames(full_data) <- c("timestep", "individual", colnames)
    data[["full_data"]] <- full_data
    
  }

  return(data)
}
