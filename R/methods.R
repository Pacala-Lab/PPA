
## Methods file

## this file contains functions which define methods for the S3 classes defined in this package

## --------------------- IS. METHODS --------------------------

## methods which define class inheritance, each class requires one

is.model_light <- function(x) inherits(x, "model_light")

is.allometry <- function(x) inherits(x, "allometry")

## --------------------- SUMMARY. METHODS --------------------------

## methods which print summary output for S3 class objects:

#'SUMMARY.model_light
#'
#'summary method for model_light object
#'
#'@param simulation a simulation generate by model_light()
#'
#'@return a summary of relevant data from the input simulation

summary.model_light <- function(simulation) {

  ## generate table of means and sds

  varnames <- c("diameter", "diameter_growth", "height", "crown_area", "struct_biomass")[c("diameter", "diameter_growth", "height", "crown_area", "struct_biomass") %in% names(simulation)]
  table_text <- matrix(0, nrow = length(varnames), ncol = 2)
  colnames(table_text) <- c("mean", "std. deviation")
  rownames(table_text) <- varnames
  for(i in 1:length(varnames)) {
    if (varnames[i] == "diameter_growth") t <- simulation$max_t
    else t <- simulation$max_t+1

    table_text[i,1] <- mean(simulation[[varnames[i]]][, t], na.rm = TRUE)
    table_text[i,2] <- sd(simulation[[varnames[i]]][, t], na.rm = TRUE)
  }

  cat("Simulation:", "\n",
      "Max time: ", as.character(simulation[["max_t"]]), "\n",
      "Time Step: ", as.character(simulation[["timestep"]]), "\n",
      "Iterations: ", as.character(simulation[["max_t"]]/simulation[["timestep"]]), "\n",
     "Plot Area:", as.character(simulation[["plot_area"]]), "\n",
      "\n",
     "Means and standard deviations for variables at last time step:", "\n")
  print(table_text)
}



summary.model_clonal <- function(simulation) {



}


## --------------------- PLOT. METHODS --------------------------

## methods which plot S3 class objects:


#'PLOT_ONERESULT
#'
#'This function plots the time trajectory of a single dimensional quantity from a simulation.
#'
#'@param data full.data object from simulation
#'@param y.val The value to plot on the y axis
#'
#'@return returns a plot of y.val vs time

## not an S3 method. This gets wrapped into the plot.simulation method
plot_oneresult <- function(data, y.val) {
  p <- ggplot(aes_string(x="timestep", y.val, color = "individual", linetype = "spp"), data = data) +
    geom_line(size = 0.5) +
    xlab("Time Step") +
    ylab(y.val) +
    scale_color_discrete(guide = FALSE) +
    scale_linetype_discrete(name = "species") +
    theme_bw()

  legend <- gtable_filter(ggplotGrob(p), "guide-box")

  p <- p + theme(legend.position = "none")

  return(list(p, legend))
}

#'plot.model_light
#'
#'plot method for model_light objects
#'
#'@param simulation an object of class "model_light"
#'@param y.values which dimensional qualities would you like to plot
#'
#'@return returns a grid of time trajectory plots
#'

## plot method for simulation class
plot.model_light <- function(simulation, y.values) {
  ind <- "full_data"

  if (ind == "full_data" & is.null(simulation[["full_data"]])) {
    print("Error: full data must be defined in simulation")
    return()
  }

  data <- as.data.frame(simulation[[ind]])
  data[, "crown_class"] <- factor(data[, "crown_class"])
  data[, "spp"] <- factor(data[, "spp"])
  data[, "individual"] <- factor(data[, "individual"])

  if (missing(y.values)) {
    y.values <- as.list(colnames(data)[!(colnames(data) %in% c("timestep", "individual","crown_class", "spp", "diameter_growth"))])
  }

  plot.list <- lapply(y.values, plot_oneresult, data = data)
  eval.list <- character()

  legend <- plot.list[[1]][[2]]

  for (i in 1:length(y.values)) {
    if (i == length(y.values)) eval.list <- paste0(eval.list, "plot.list[[", i, "]][[1]]")
    else eval.list <- paste0(eval.list, "plot.list[[", i, "]][[1]], ")
  }

  eval(parse(text = paste0("grid.arrange(arrangeGrob(", eval.list, ", nrow = 2, top = textGrob('Stand Simulation Output')), legend, widths=unit.c(unit(1, 'npc') - legend$width, legend$width), nrow = 1)")))

}



## plot methods for model_clonal

plot.model_clonal <- function() {



}


## --------------------- INDEXING. METHODS --------------------------

## define indexing method for allometry class.
"[.allometry" <- function(x, i, j) {

  output <- x$allometry[i, j]
  return(output)
  
}


## --------------------- ALLOMETRY. METHODS --------------------------

## allometry methods for the various simulation class types:


## eventually this function should have methods for each allometry type supported by the package:
## i.e. create_allometry.clonal, create_allometry.annual, create_allometry.forest (idk how specific we want to get? spp?)

## this should eventually have different methods for the different model classes
## (i.e. create_allometry.forest, create_allometry.clonal, create_allometry.hydro)

## some of the parameters in the resulting output are not, strictly speaking, allometric parameters.
## Perhaps this should be called "defaults", as Marco suggested, or those non-allometric parameters could
## be specified in a seperate defaults method.


#'CREATE_ALLOMTRY.MODEL_LIGHT
#'
#'creates an allometry object to be used in the model_light simulator
#'
#'@param n.spp the number of spp to create allometry for
#'@param custon a named light supplying custom values for any or all allometric parameters
#'@param stochastic logical, if true inject random variation into the parameters
#'
#'@return returns an object of class "allometry" to be used in model_light simulator


create_allometry.model_light <- function(n.spp, custom, stochastic = TRUE) {
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
  ## needs major refinement
  if (stochastic) {
    for (spp in as.character(1:n.spp)) {
      out[, spp] <- rnorm(nrow(out), mean = out[, spp], sd = 0.05*out[, spp])

    }
  }
  outlist <- list(allometry = out)
  class(outlist) <- "allometry"

  return(outlist)
}

## placeholder for illustration purposes"
create_allometry.clonal <- function() {


}

#'CALCULATE_ALLOMETRY.model_light
#'
#'function to calculate full allometry from object of class simulation
#'
#'@param simulation a simulation generated by model_light
#'@param allometry an object of class "allometry"
#'@param spp a named list of species
#'
#'@return returns datasets for height, crown_area and biomass

calculate_allometry.model_light <- function(simulation, allometry, spp) {

  if(!is.model_light(simulation)) {
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
