## Utility Functions
## anything that doesn't really belong somewhere else, or is generally useful

#' EXTRACT_COLUMN
#'
#' function to extract a single column from a simulation data matrix
#'
#'@param data the simulation dataset that contains the column you want to extract
#'@param colname the name of the column to extract from data
#'
#'@return returns the column specified

extract_column <- function(colname, data) {

  data_matrix <- data[[colname]]

  if(!is.matrix(data_matrix)) {
    print("Error: column matrix must be of class matrix")
    return()
  }

  column <- as.numeric(data_matrix)
  return(column)

}

#' SPP_POPULATIONS
#'
#' function to calculate the number of species of a given crown_class at a given timepoint
#'
#'@param data the simulation dataset
#'@param spp the name of the species whose population you want to extract
#'@param i the timestep at which you want to extract the populatioon
#'@param crown_class which crown class's population do you wish to extract
#'
#'@return returns the population for the species, timestep and crown class you specified
spp_populations <- function(data, spp, i, crown_class) {
  if(missing(crown_class)) {
    out <- sum(data[["spp"]][,i] == spp)
  }
 else {
    out <- sum(data[["crown_class"]][, i] == crown_class & data[["spp"]][, i] == spp)
  }
  return(out)
}

#'CALC.CA
#'
#'A function to calculate the crown area of a plant given its allometric constants and diameter. Works specifically within
#'the context of stand simulators
#'
#'@param alpha_w allometric constant governing the relationship between diameter and crown area
#'@param gamma allometric exponent governing the relationship between diameter and crown area
#'@param spp the name of the species for which you are calculating CA
#'@param diam_data diameter data for the species, specifically from within a model_ call
#'@param spp_data spp data from model_ call
#'@param i current timestep
#'
#'@return returns a numeric value, the crown area for the specied plant or plants

## need to separate diameter and spp data in order to work with mapply funciton - this is annoying, would rather be able to past data as arg, but mapply tries to interpret it as multiple args
calc.CA <- function(alpha_w, gamma, spp, diam_data, spp_data, i) {

  out <- sum(alpha_w * (diam_data[spp_data[, i + 1] == spp, i + 1]^gamma))
  return(out)

}


#'CALC.A
#'
#' this function calculates max photosynthetic rates given allometric constants, crown classes and the spp info
#'
#'@param allometry an object of class "allometry"
#'@param n.crown.class The number of crown classes to calculate A for, currently only 2 supported
#'@param spp a vector of species names, similar to the input for model_ calls
#'
#'@return returns the maximum photosynthetic rate

calc.A <- function(allometry, n.crown.class = 2, spp) {

  spp <- as.character(spp)
  out <- data.frame(matrix(1, nrow = length(spp), ncol = n.crown.class+1))
  colnames(out) <- c("spp","A_c", "A_u")[1:(n.crown.class+1)]
  out[, "spp"] <- spp

  ## calculate maximum light level for understory trees
  L_u <- allometry["L_0", spp] * (exp(-(allometry["k", spp]) * allometry["l_c", spp]))

  out[, "A_c"] <- as.numeric((allometry["V", spp] / allometry["k", spp]) * (1 + log((allometry["a_f", spp] * allometry["L_0", spp])/allometry["V", spp]) - ((allometry["a_f", spp] * allometry["L_0", spp]) / allometry["V", spp]) * exp(-allometry["k", spp] * allometry["l_c", spp])))

  if(n.crown.class == 2) {
    out[, "A_u"] <- as.numeric((allometry["a_f", spp] * L_u) / allometry["k", spp] * (1 - exp(-allometry["k", spp] * allometry["l_u", spp])))
  }
  if(n.crown.class > 2) {
    print("error: only 2 crown classes currently supported.")
    return()
  }
  return(out)
}

#'CALC.DD
#'function to calculate diameter growth
#'
#'@param diameter the current diameter of the plant
#'@param allometry an object of class "allometry"
#'@param spp the spp name
#'@param l the leaf area of the plant
#'@param r the root area of the plant'
#'
#'@return returns a diameter increment

calc.dd <- function(diameter, allometry, spp, l, r, A) {
  dd <- 1/((allometry["alpha_s", spp] * (allometry["gamma", spp]+1) * (1+allometry["c_bg", spp]) * (1/allometry["alpha_w", spp])) + ((allometry["gamma", spp]/diameter) * (l * allometry["c_lb", spp] + r * allometry["c_rb", spp]))) * (A - (l*allometry["c_l", spp]) - (r*allometry["c_r", spp]))
  return(dd)
}


#'ORDER_SIMULATION
#'
#' function to order the simulation data object
#' currently written like a method, could be moved to a methods file if relevant for multiple classes
#'
#'@param data dataset to reorder
#'@param by dimensional value to reorder by
#'@param decreasing if TRUE, ordered from large to small
#'@param i current timestep
#'
#'@return returns a reordered set of data matrices

order_simulation <- function(data, by, decreasing = TRUE, i) {
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


#'SIMULATION_AVERAGES
#'
#'this function calculates average values for parameters of interest at each time step.
#'
#'@param simulation an object that is a result of a model_ call
#'
#'@return returns the average dimensional values of the plants at each timestep in the simulation

## function is of limited utility and we might consider getting rid of it
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
