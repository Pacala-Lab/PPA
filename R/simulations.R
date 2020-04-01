
## code for the class: simulation

## this is the primary simulation function. Should we rename class/function to PPA?

simulate_stand <- function(allometry = create_allometry(n.spp = 1),
                           spp = 1,
                           max_t,
                           timestep,
                           dnot = rnorm(n = n_start, mean = 0.0004, sd = 0.00005),
                           plot.area,
                           n_start,
                           full.allom = TRUE,
                           avgs = TRUE,
                           full.data = TRUE) {

  ## very script-heavy method. Will change to stopifnot eventually
  if(ncol(allometry) < length(spp)) {
    print("Error: allometry object must include values for each spp in simulation")
    return()


  if(!is.numeric(max_t)) {
    print("Error: max_t must be numeric")
    return()
  }

  if(!is.numeric(timestep)) {
    print("Error: timestep must be numeric")
    return()
  }
  ##stopornot/stopifnot
  if(!is.numeric(dnot)) {
    print("Error: dnot must be numeric")
    return()
  }

  if(!is.numeric(plot.area)) {
    print("Error: plot.area must be numeric")
    return()
  }

  if(timestep > max_t) {
    print("Error: timestep cannot be larger than max_t")
    return()
  }

  if(!is.logical(full.allom)) {
    print("Error: full.allom must be TRUE/FALSE")
    return()
  }

  if(length(spp) > 1 && n_start <= 1) {
    print("Error: n_start must be same length as spp.")
    return()
  }

  if(is.null(names(n_start))) {
    print("Error: n_start must be a named vector, names must match species names")
    return()
  }

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

    pre-calculate # of each spp in each crown class
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
    ## reorder data to easily index align death vector and data matrices


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
  class(data) <- "simulation"

  ## if user specifies full allometry, calculate and add to the data
  if (full.allom) {
    data <- calculate_allometry(data, spp = spp)
  }

  ## if user specifies average values are wanted, calculate and add as list element
  if(avgs) {
    data[["avgs"]] <- simulation_averages(data)
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



#######----------------- METHODS AND CLASS SPECIFIC UTILITY FUNCTIONS ---------------------########
###################################################################################################

## define is. method for simulation class
is.simulation <- function(x) inherits(x, "simulation")

## function to order the simulation data object
order.simulation <- function(data, by, decreasing = TRUE, i) {
  ## currently an inelegant solution with if and, could be better?
  if (length(by) == 1) {
    data[["spp"]] <- data[["spp"]][order(data[[by]][,i], decreasing = decreasing), ]
    data[["diameter"]] <- data[["diameter"]][order(data[[by]][,i], decreasing = decreasing), ]
    data[["diameter_growth"]] <- data[["diameter_growth"]][order(data[[by]][,i], decreasing = decreasing), ]
    data[["crown_class"]] <- data[["crown_class"]][order(data[[by]][,i], decreasing = decreasing), ]
  }
  else if (length(by) == 2) {
    data[["spp"]] <- data[["spp"]][order(data[[by[1]]][,i], data[[by[2]]][,i], decreasing = decreasing), ]
    data[["diameter"]] <- data[["diameter"]][order(data[[by[1]]][,i], data[[by[2]]][,i], decreasing = decreasing), ]
    data[["diameter_growth"]] <- data[["diameter_growth"]][order(data[[by[1]]][,i], data[[by[2]]][,i], decreasing = decreasing), ]
    data[["crown_class"]] <- data[["crown_class"]][order(data[[by[1]]][,i], data[[by[2]]][,i], decreasing = decreasing), ]
  }
  else {
    print("Error: function currently only supports two by arguments")
    return()
  }

  return(data)
}

## not an S3 method. This gets wrapped into the plot.simulation method
plot_simulation_result <- function(data, y.val) {
  ggplot(aes_string(x="timestep", y.val, color = "spp", linetype = "crown_class"), data = data) +
    geom_point(size = 0.5) +
    geom_smooth(method = gam, formula = y~s(x, bs="cs", sp=3), se = FALSE) +
    xlab("Time Step") +
    ylab(y.val) +
    scale_color_discrete(name = "species") +
    scale_linetype_discrete(name = "crown class") +
    theme(title = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", size = 1))
}


## plot method for simulation class
plot.simulation <- function(simulation, avgs = FALSE, y.values) {
  if (avgs) ind <- "avgs"
  else ind <- "full_data"

  if (avgs & is.null(simulation[["avgs"]])) {
    print("Error: avgs must be defined in simulation")
    return()
  }

  if (ind == "full_data" & is.null(simulation[["full_data"]])) {
    print("Error: full data must be defined in simulation")
    return()
  }

  data <- as.data.frame(simulation[[ind]])
  data[, "crown_class"] <- factor(data[, "crown_class"])
  data[, "spp"] <- factor(data[, "spp"])

  if (missing(y.values)) {
    y.values <- as.list(colnames(data)[!(colnames(data) %in% c("timestep", "individual","crown_class", "spp"))])
  }

  plot.list <- lapply(y.values, plot_simulation_result, data = data)
  eval.list <- character()

  for (i in 1:length(y.values)) {
    if (i == length(y.values)) eval.list <- paste0(eval.list, "plot.list[[", i, "]]")
    else eval.list <- paste0(eval.list, "plot.list[[", i, "]], ")
  }

  eval(parse(text = paste0("grid.arrange(", eval.list, ", nrow = 2, top = textGrob('Stand Simulation Output'))")))

}


## this function is broken/needs improvement in general
## want it to print essential/useful imformation from simulation output.
summary.simulation <- function(simulation) {

  if(!full.data) {
    print("Error: simulation must have full.data element for summary method")
    return()
  }

  ## generate table of means and sds
  varnames <- c("diameter", "diameter_growth", "height", "crown_area", "struct_biomass")[c("diameter", "diameter_growth", "height", "crown_area", "struct_biomass") %in% names(simulation)]
  table_text <- matrix(0, nrow = length(varnames), ncol = 2)
  colnames(table_text) <- c("mean", "std. deviation")
  rownames(table_text) <- varnames
  for(i in 1:length(varnames)) {
    table_text[i,2] <- mean(simulation[[varname]][simulation[["full.data"]][,"timestep"]==max_t+1,varname])
    table_text[i,3] <- sd(simulation[[varname]][simulation[["full.data"]][,"timestep"]==max_t+1,varname])
  }

  cat("Simulation:", "\n",
      "Max time: ", as.character(simulation[["max_t"]]), "\n",
      "Time Step: ", as.character(simulation[["timestep"]]), "\n",
      "Iterations: ", as.character(simulation[["max_t"]]/simulation[["timestep"]]), "\n",
     "Plot Area:", as.character(simulation[["plot.area"]]), "\n",
      "\n",
      "Means and standard deviations for variables at last time step", "\n")
  return(x)
}

## this function calculates average values for parameters of interest at each time step.
## its of limited utility and I am thinking about deleting it.

simulation_averages <- function(simulation) {

  ## check that input is of the class simulation
  if(!is.simulation(simulation)) {
    print("Error: input must be class simulations")
    return()
  }

  ## create list of potential variable names, initialize matrix to populate with averages
  varnames <- c("diameter", "diameter.growth", "height", "crown_area", "struct_biomass")
  cnames <- c("timestep","crown.class", varnames[varnames %in% names(simulation)])
  avgs <- matrix(1, nrow = 2*((simulation[["max_t"]]/simulation[["timestep"]])+1), ncol = length(cnames))
  colnames(avgs) <- cnames

  avgs[,2] <- c(rep(1, times = (simulation[["max_t"]]/simulation[["timestep"]])+1), rep(2, times = (simulation[["max_t"]]/simulation[["timestep"]])+1))
  avgs[,1] <- rep(seq(1, simulation[["max_t"]]+1, by = simulation[["timestep"]]), times = 2)
  for (i in 3:length(cnames)) {
    for(j in 1:(nrow(avgs)/2)) {
      ## no diameter growth in last timestep
      if(cnames[i]=="diameter.growth" & j == nrow(avgs)/2) {
        avgs[j,i] <- NA
        avgs[j+(nrow(avgs)/2),i] <- NA
      }
      else {
        avgs[j,i] <- mean(simulation[[cnames[i]]][simulation[["crown.class"]][,j]==1,j])

        if (any(simulation[["crown.class"]][,j]==2)) {
          avgs[j+(nrow(avgs)/2),i] <- mean(simulation[[cnames[i]]][simulation[["crown.class"]][,j]==2,j])
        }
        else avgs[j+(nrow(avgs)/2),i] <- NA
      }
    }
  }
  return(avgs)
}
