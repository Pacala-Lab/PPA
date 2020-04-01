## This file contains functions related to allometric scaling of the plants


## eventually this function should have methods for each allometry type supported by the package:
## i.e. create_allometry.clonal, create_allometry.annual, create_allometry.forest (idk how specific we want to get? spp?)

create_allometry <- function(n.spp, custom, stochastic = TRUE) {
  ## allometry will be a dataframe of allometric variables, with rows for vars and columns for spp
  ## this function also assigns values for non-allometric, but relevant variables. Not sure if should be in separate file

  ## for now, almost all numbers will be same for all species, some will be drawn from random dist.
  out <- as.data.frame(matrix(1, nrow = 25, ncol = n.spp))
  colnames(out) <- c(as.character(seq(1, n.spp, by = 1)))
  rownames(out) <- c("l_c", "l_u", "r_c", "r_u", "mu_c", "mu_u", "r", "H", "alpha_s", "alpha_w", "gamma", "c_lb",
               "c_rb", "c_bg", "tau_l", "tau_r", "p_r", "p_l", "p_sw", "c_l", "c_r", "V", "k", "L_0", "a_f")

  ## function accepts a list of named vectors for each species with custom allometric constants
  if(!missing(custom)) {
    list_classes <- lapply(custom, class)
    list_names <- lapply(custom, names)
    if(class(custom) != "list") {
      print("Error: custom variables must be a list of named numeric vectors")
      return()
    }
    else if(any(list_classes != "numeric")) {
      print("Error: custom variables must be a list of named numeric vectors")
      return()
    }
    else if(any(is.null(list_names))) {
      print("Error: custom variables must be a list of named numeric vectors")
      return()
    }
  }

  for (spp in as.character(1:n.spp)) {

    #### allometric variables
    out["l_c", spp] <- 4               ## leaf area index of tree canopy tree
    out["l_u", spp] <- 1               ## leaf are index of understory tree
    out["r_c", spp] <- 3
    out["r_u", spp] <- 2
    out["mu_c", spp] <- 0.005
    out["mu_u", spp] <- 0.015
    out["r", spp] <- 6.5             ## fine-root surface area per unit crown area
    out["H", spp] <- 3.6             ## allometric constant for height
    out["alpha_s", spp] <- 0.0815        ## allometric constant for sapwood
    out["alpha_w", spp] <- 0.20          ## allometric constant for crown area
    out["gamma", spp] <- 1.5         ## allometric exponent

    #### carbon accumulation and allocation variables
    out["c_lb", spp] <- 0.07656      ## cost of building leaf biomass per LAI
    out["c_rb", spp] <- 0.02933      ## cost of builing root biomass per LAI
    out["c_bg", spp] <- 0.2          ## building cost of structural biomass
    out["tau_l", spp] <- 1           ## average lifetime of a unit carbon in leaf
    out["tau_r", spp] <- 2           ## average lifetime of a unit carbon in roots
    out["p_r", spp] <- 0.02933       ## respiration rate of roots per unit surface area
    out["p_l", spp] <- 0.0638        ## respiration rate of leaves per unit surface area
    out["p_sw", spp] <- 0.0466       ## respiration rate of sapwood
    out["c_l", spp] <- out["c_lb", spp] / out["tau_l", spp] + out["p_l", spp] + out["p_r", spp] ## cost of building and maintaining leaf biomass per unit LAI
    out["c_r", spp] <- out["c_rb", spp] / out["tau_r", spp] + out["p_r", spp]       ## cost of building and maintaining root biomass per unit RAI
    out["V", spp] <- 0.6             ## maximum rate of carbon fixation - kgC/m^2/day
    out["k", spp] <- 0.33            ## light extinction coefficient
    out["L_0", spp] <- 1200          ## light level above highest canopy
    out["a_f", spp] <- 0.001         ## conversion rate from photons to carbohydrates
  }

  ## reassign custom variables. Perhaps there is a way to do his that doesnt require overwriting?
  if(!missing(custom)) {
    for (i in 1:length(custom)) {
      for (j in 1:length(custom[[i]])) {
        out[var = names(custom[[i]])[j], i] <- custom[[i]][j]
      }
    }
  }


  ## crude implementation of random variation for now
  if (stochastic) {
    for (spp in as.character(1:n.spp)) {
      out[, spp] <- rnorm(nrow(out), mean = out[, spp], sd = 0.1*out[, spp])

    }
  }

  return(out)
}


########-------------------------METHODS AND CLASS RELATED UTILITY FUNCTIONS -------------------#######
#######################################################################################################

is.allometry <- function(x) inherits(x, "allometry")

## function to calculate full allometry from object of class simulation
calculate_allometry <- function(simulation, spp) {

  if(!is.simulation(simulation)) {
    print("Error: input must be of class simulation")
    return()
  }
  
  simulation[["height"]] <- matrix(1, nrow = nrow(simulation[["diameter"]]), ncol = ncol(simulation[["diameter"]]))
  simulation[["crown_area"]] <- matrix(1, nrow = nrow(simulation[["diameter"]]), ncol = ncol(simulation[["diameter"]]))
  simulation[["struct_biomass"]] <- matrix(1, nrow = nrow(simulation[["diameter"]]), ncol = ncol(simulation[["diameter"]]))

  for (s in spp) {

    simulation[["height"]][simulation[["spp"]] == s] <- allometry["H", s] * (simulation[["diameter"]][simulation[["spp"]] == s] ^ (allometry["gamma", s] - 1))
    simulation[["crown_area"]][simulation[["spp"]] == s] <- allometry["alpha_w", s] * (simulation[["diameter"]][simulation[["spp"]] == s] ^ allometry["gamma", s])
    simulation[["struct_biomass"]][simulation[["spp"]] == s] <- allometry["alpha_s", s] * (simulation[["diameter"]][simulation[["spp"]] == s] ^ (allometry["gamma", s] + 1))

  }

  return(simulation)
}


## need to separate diameter and spp data in order to work with mapply funciton - this is annoying, would rather be able to past data as arg, but mapply tries to interpret it as multiple args
calc.CA <- function(alpha_w, gamma, spp, diam_data, spp_data, i) {

  out <- sum(alpha_w * (diam_data[spp_data[, i + 1] == spp, i + 1]^gamma))
  return(out)

}

## this function calculates max photosynthetic rates given allometric constants, crown classes and the spp info
calc.A <- function(allometry, n.crown.class = 2, spp) {

  spp <- as.character(spp)
  out <- data.frame(matrix(1, nrow = length(spp), ncol = n.crown.class+1))
  colnames(out) <- c("spp","A_c", "A_u")[1:(n.crown.class+1)]
  out[, "spp"] <- spp

  ## calculate maximum light level for understory trees
  L_u <- allometry["L_0", spp] * (exp(-(allometry["k", spp]) * allometry["l_c", spp]))

  out[, "A_c"] <- as.numeric((allometry["V", spp] / allometry["k", spp]) * (1 + log((allometry["a_f", spp] * allometry["L_0", spp])/allometry["V", spp]) - ((allometry["a_f", spp] * allometry["L_0", spp]) / allometry["V", spp]) * exp(-allometry["k", spp] * allometry["l_c", spp])))

  if(n.crown.class == 2) {
    out[, "A_u"] <- as.numeric((allometry["a_f", spp] * L_u) / allometry["V", spp] * exp(-allometry["k", spp] * allometry["l_u", spp]))
  }
  if(n.crown.class > 2) {
    print("error: only 2 crown classes currently supported.")
    return()
  }
  return(out)
}

## function to calculate diameter growth
calc.dd <- function(diameter, allometry, spp, l, r, A) {
  dd <- 1/((allometry["alpha_s", spp] * (allometry["gamma", spp]+1) * (1+allometry["c_bg", spp]) * (1/allometry["alpha_w", spp])) + ((allometry["gamma", spp]/diameter) * (l * allometry["c_lb", spp] + r * allometry["c_rb", spp]))) * (A - (l*allometry["c_l", spp]) - (r*allometry["c_r", spp]))
  return(dd)
}
