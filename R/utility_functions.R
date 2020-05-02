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

#'CALC.A_CLONAL
#'
#' this function calculates max photosynthetic rates given allometric constants, crown classes and the spp info for model_clonal
#'
#'@param allometry an object of class "allometry"
#'@param n.crown.class The number of crown classes to calculate A for, currently only 2 supported
#'@param spp a vector of species names, similar to the input for model_ calls
#'
#'@return returns the maximum photosynthetic rate

calc.A_clonal <- function(allometry, n.crown.class = 2, spp) {
  
  spp <- as.character(spp)
  out <- data.frame(matrix(1, nrow = length(spp), ncol = n.crown.class+1))
  colnames(out) <- c("spp","A_c", "A_u")[1:(n.crown.class+1)]
  out[, "spp"] <- spp
  
  ## calculate maximum light level for understory trees
  L_u <- allometry["L_0", spp] * (exp(-(allometry["k", spp]) * allometry["l_c", spp]))
  
  out[, "A_c"] <- as.numeric((allometry["V", spp] / allometry["k", spp]) * (1 - exp(-allometry["k", spp] * allometry["l_c", spp])))
  
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

#'CALC.DD_NONCLONAL
#'function to calculate diameter growth rate, i.e., dD/dt for a non-clonal plant in model_clonal
#'
#'@param diameter the current diameter of the plant
#'@param allometry an object of class "allometry" 
#'@param spp the spp name
#'@param understory logical, if true calculate the diameter growth rate for an understory plant; if false calculate the diameter growth rate for a canopy tree.
#'
#'@return returns a diameter growth rate at a particular diameter size

calc.dd_nonclonal <- function(diameter, allometry, spp, understory=TRUE) {
  spp <- as.character(spp)
  CarbonGain<-calc.A_clonal(allometry = allometry,spp=spp)
  A<-CarbonGain[spp==spp,"A_c"]
  dDdt<-(allometry["phi_l",spp]*diameter^allometry["theta_l",spp]*
           (A-(allometry["epsilon_l",spp]+allometry["r_l",spp])*
              allometry["l_c",spp]-(allometry["epsilon_r",spp]+allometry["r_r",spp])*
              allometry["delta_rl",spp]*allometry["l_c",spp]/allometry["omega",spp]-allometry["f", spp]))/
    (allometry["theta_s",spp]*allometry["phi_s",spp]*diameter^(allometry["theta_s",spp]-1)+(allometry["sigma",spp]+
                                                                                              allometry["delta_rl",spp]*allometry["l_c",spp]*allometry["omega",spp])*allometry["phi_l",spp]*allometry["theta_l",spp]*diameter^(-1))
  
  if(understory){
    A<-CarbonGain[spp==spp,"A_u"]
    dDdt=(allometry["phi_l",spp]*diameter^allometry["theta_l",spp]*
            (A-(allometry["epsilon_l",spp]+allometry["r_l",spp])*
               allometry["l_u",spp]-(allometry["epsilon_r",spp]+allometry["r_r",spp])*
               allometry["delta_rl",spp]*allometry["l_u",spp]/allometry["omega",spp]))/
      (allometry["theta_s",spp]*allometry["phi_s",spp]*diameter^(allometry["theta_s",spp]-1)+(allometry["sigma",spp]+
                                                                                                allometry["delta_rl",spp]*allometry["l_u",spp]*allometry["omega",spp])*allometry["phi_l",spp]*allometry["theta_l",spp]*diameter^(-1))
  }
  ifelse(dDdt<0,print("Error: diameter growth cannot be negative"),return(dDdt))
}

####This function below really needs to be faster. Currently I am using a loop because I don't know a better way...

#'EQUIL_NONCLONAL
#'function to pre-caculate the equilibrium diameter size for a non-clonal plant in the understory over time in model_clonal.
#'
#'@param dnot the initial diameter size of a seedling to begin with, in meters
#'@param allometry an object of class "allometry" 
#'@param spp the spp name
#'@param max_t the total number of time steps of the simulation
#'@param timestep the length of each timestep, in years
#'
#'@return returns the equilibrium values of diameter and the number of timesteps to reach the canopy at population equilibrium

equil_nonclonal<-function(dnot,allometry,spp,max_t,timestep){
  out <- as.data.frame(matrix(1, nrow =max_t, ncol = 2))
  colnames(out) <- c("diameter","LRS_test")
  rownames(out) <- c(as.character(seq(1,max_t, by = 1)))
  
  ##intergrand is the equation for the number of seeds a tree reproduces during its time in the canopy
  out[1,"diameter"]<-dnot
  for (i in 2:max_t){
    out[i,"diameter"]<-out[i-1,"diameter"]+calc.dd_nonclonal(diameter=out[i-1,"diameter"],allometry=allometry,spp = spp)*timestep
    intergrand<-function(t){
      value<-exp(-(allometry["mu_c",spp])*t)*out[i,"diameter"]^allometry["theta_l",spp]*allometry["F",spp]*allometry["g",spp]*allometry["phi_l",spp]
      return(value)
    }
    out[i,"LRS_test"]<-exp(-(allometry["mu_u",spp])*timestep*i)*integrate(intergrand,lower = 0,upper=Inf)$value
    if(out[i,"LRS_test"]>1){
    }else{
      break
    }
  }
  ##dstar is the diameter size when a tree reaches canopy in its population at equilibrium, and tstar is the corresponding time it needs
  tstar<-i-1
  dstar<-out[tstar,"diameter"]
  understory_diameter<-as.matrix(out[1:tstar,1])
  colnames(understory_diameter)<-"diameter"
  outlist<-list("understory_diameter"=understory_diameter,"tstar"=tstar,"dstar"=dstar)
  return(outlist)
}


#'GET.GROWTH_NONCLONAL
#'function to calculate the diameter growth rates dD/dt in the understory and in the canopy based on equilibrium diameter size for a non-clonal plant in model_clonal
#'
#'@param dnot the initial diameter size of a seedling to begin with, in meters
#'@param allometry an object of class "allometry" 
#'@param spp the spp name
#'@param max_t the total number of time steps of the simulation
#'@param timestep the length of each timestep, in years
#'@param full.allom logical, if true returns allometric scaling results along with diameter sizes
#'
#'@return returns a data frame containing diameter, bolewood mass, crown area, height based on defined allometry.

get.growth_nonclonal<-function(dnot,allometry,spp,max_t,timestep,full.allom=TRUE){
  out <- as.data.frame(matrix(2, nrow =max_t, ncol = 2))
  colnames(out) <- c("diameter","crown_class") ##class 2 is understory, 1 is canopy
  rownames(out) <- c(as.character(seq(1,max_t, by = 1)))
  
  ##use equil_nonclonal to generate dstar and understory diameter sizes 
  u.growth.list<-equil_nonclonal(dnot=dnot,allometry=allometry,spp = spp,max_t=max_t,timestep = timestep)
  dstar<-u.growth.list[["dstar"]]
  tstar<-u.growth.list[["tstar"]]
  out[1:tstar,"diameter"]<-u.growth.list[["understory_diameter"]]
  out[((tstar+1):max_t),"crown_class"]<-1
  for(i in (tstar+1):max_t){
    out[i,"diameter"]<-out[i-1,"diameter"]+calc.dd_nonclonal(diameter=out[i-1,"diameter"],allometry=allometry, spp=spp, understory = FALSE)*timestep
  }
  if(full.allom){
    out[,"time"]<-NA
    out[(tstar+1):max_t,"time"]<-seq(1,max_t-tstar)*timestep
    out[,"height"]<-allometry["phi_z",spp]*out[,"diameter"]^allometry["theta_z",spp]
    out[,"crown_area"]<-allometry["phi_l",spp]*out[,"diameter"]^allometry["theta_l",spp]
    out[,"bolewood_mass"]<-allometry["phi_s",spp]*out[,"diameter"]^allometry["theta_s",spp]
    out[1:tstar,"step_fecundity"]<-0
    out[(tstar+1):max_t,"step_fecundity"]<- exp(-allometry["mu_u",spp]*tstar*timestep)*exp(-allometry["mu_c",spp]*out[(tstar+1):max_t,"time"])*allometry["F",spp]*allometry["g",spp]*out[(tstar+1):max_t,"crown_area"]*timestep
    out[,"cumsum"]<-cumsum(out[,"step_fecundity"])
    }
  return(out)
}

####eventually, the function below should be able to calculate dD/dt for each stem for n-stemmed plants. I'm still working on the math, and so far 
####no.stem=2 for a simple root suckering tree. There can be more variants of this function in the future, but now a clonal tree starts building
####its second trunk when the first trunk reaches size D*. At present, cloning only happens when the primary stem is in full light.

#'CALC.DD_CLONAL_NURSING
#'function to calculate diameter growth rate dD/dt for a clonal plant's nth trunk/ramet before disconnection happens.
#'
#'@param diameter.prim the diameter size of the primary stem 
#'@param diameter.ramet a vector of the diameter sizes of the secondary, third... nth stem
#'@param no.stem total number of stems/ramets in the plant's clone
#'@param allometry an object of class "allometry" 
#'@param spp the spp name
#'here the plant's status depends on whether the first ramet/stem reaches dstar. 
#'
#'@return returns diameter growth rates for each stem at a particular time point

calc.dd_clonal_nursing<- function(diameter.prim,diameter.ramet,no.stem=2,allometry,spp){
  spp <- as.character(spp)
  CarbonGain<-calc.A_clonal(allometry = allometry,spp=spp)
  A<-CarbonGain[spp==spp,"A_c"]
  ##canopy diameter growth rate for primary stem
  dDdt<-(allometry["phi_l",spp]*diameter.prim^allometry["theta_l",spp]*
           (A-(allometry["epsilon_l",spp]+allometry["r_l",spp])*
              allometry["l_c",spp]-(allometry["epsilon_r",spp]+allometry["r_r",spp])*
              allometry["delta_rl",spp]*allometry["l_c",spp]/allometry["omega",spp]-allometry["f", spp]-allometry["v", spp]))/
    (allometry["theta_s",spp]*allometry["phi_s",spp]*diameter.prim^(allometry["theta_s",spp]-1)+(allometry["sigma",spp]+
                                                                                                   allometry["delta_rl",spp]*allometry["l_c",spp]*allometry["omega",spp])*allometry["phi_l",spp]*allometry["theta_l",spp]*diameter.prim^(-1))
  ##nursed growth rate for secondary stem 
  A_ramet<-CarbonGain[spp==spp,"A_u"]
  dDdt_ramet<-(allometry["v",spp]*allometry["phi_l",spp]*diameter.prim^allometry["theta_l",spp]+allometry["phi_l",spp]*diameter.ramet^allometry["theta_l",spp]*
                 (A_ramet-(allometry["epsilon_l",spp]+allometry["r_l",spp])*
                    allometry["l_u",spp]-(allometry["epsilon_r",spp]+allometry["r_r",spp])*
                    allometry["delta_rl",spp]*allometry["l_u",spp]/allometry["omega",spp]))/
    (allometry["theta_s",spp]*allometry["phi_s",spp]*diameter.ramet^(allometry["theta_s",spp]-1)+(allometry["sigma",spp]+
                                                                                                    allometry["delta_rl",spp]*allometry["l_u",spp]*allometry["omega",spp])*allometry["phi_l",spp]*allometry["theta_l",spp]*diameter.ramet^(-1))
  if(dDdt<=0) {
    print("Error: diameter growth cannot be negative")
  }
  if (dDdt_ramet<=0){
    print("Error: diameter growth cannot be negative")
  }
  vector<-c(dDdt,dDdt_ramet)
  names(vector)<-c("prim","second")
  return(vector)
}

#'GET.GROWTH_CLONAL
#'function to calculate the diameter growth rates dD/dt in the understory and in the canopy based on the threshold diameter size in model_clonal. When a ramet reaches the threshold size, it disconnects from the primary stem.
#'
#'@param dnot the initial diameter size of the primary stem
#'@param dnot_ramet the initial size of a ramet's diameter 
#'@param allometry an object of class "allometry" 
#'@param spp the spp name
#'@param max_t the total number of time steps of the simulation
#'@param timestep the length of each timestep, in years
#'@param custom a vector supplying custom values for dstar and tstar--the size and time to start cloning as well as end nursing
#'@param full.allom logical, if true returns allometric scaling results along with diameter sizes
#'
#'@return returns a data frame containing diameter, bolewood mass, crown area, height based on defined allometry for the clone.

get.growth_clonal<-function(dnot,dnot_ramet, allometry,spp,max_t,timestep,custom,full.allom=TRUE){
  out <- as.data.frame(matrix(1, nrow =max_t, ncol = 2))
  colnames(out) <- c("diameter.prim","diameter.ramet") 
  rownames(out) <- c(as.character(seq(1,max_t, by = 1)))
  
  ##use equil_nonclonal or customized values to generate dstar and understory diameter sizes 
  if(!missing(custom)) {
    dstar<-custom[1]
    tstar<-custom[2]
    if(any(custom != "numeric")) {
      print("Error: custom variables must be a vector of named numeric vectors")
    }
  }
  u.growth.list<-equil_nonclonal(dnot=dnot,allometry=allometry,spp = spp,max_t=max_t,timestep = timestep)
  dstar<-u.growth.list[["dstar"]]
  tstar<-u.growth.list[["tstar"]]
  out[1:tstar,"diameter.prim"]<-u.growth.list[["understory_diameter"]]
  out[1:tstar,"diameter.ramet"]<-0
  out[tstar,"diameter.ramet"]<-dnot_ramet 
  
  ##calculate diameter sizes for the primary and secondary stems during the connected nursing period until the ramet reaches dclonal 
  for(i in (tstar+1):max_t){
    out[i,]<-out[i-1,]+calc.dd_clonal_nursing(diameter.prim=out[i-1,"diameter.prim"],diameter.ramet=out[i-1,"diameter.ramet"],allometry=allometry, spp=spp)*timestep
    if(out[i,"diameter.ramet"]>dstar){
      break
    }
  }
  
  ##calculate diameter sizes for the primary and secondary stems during the disconnected period in the canopy
  for(t in (i+1):max_t){
    out[t,]<-out[t-1,]+calc.dd_nonclonal(diameter=out[t-1,],allometry=allometry, spp=spp, understory = FALSE)*timestep
  }
  
  ##assign statuses: U is before cloning, N is nursing, C is both stems in canopy
  out[,"time"]<-NA ##time is time since canopy
  out[(tstar+1):max_t,"time"]<-seq(1,max_t-tstar)*timestep
  out[1:tstar,"status"]<-"U"
  out[(tstar+1):i,"status"]<-"N"
  out[(i+1):max_t,"status"]<-"C"
  out[,"step_fecundity.prim"]<-0
  out[,"step_fecundity.ramet"]<-0
  
  ##get tbreak
  tbreak<-nrow(out[out[,"status"]!="C",])
  
  if(full.allom){
    out[,"height.prim"]<-allometry["phi_z",spp]*out[,"diameter.prim"]^allometry["theta_z",spp]
    out[,"height.ramet"]<-allometry["phi_z",spp]*out[,"diameter.ramet"]^allometry["theta_z",spp]
    out[,"crown_area.prim"]<-allometry["phi_l",spp]*out[,"diameter.prim"]^allometry["theta_l",spp]
    out[,"crown_area.ramet"]<-allometry["phi_l",spp]*out[,"diameter.ramet"]^allometry["theta_l",spp]
    out[,"bolewood_mass.total"]<-allometry["phi_s",spp]*out[,"diameter.prim"]^allometry["theta_s",spp]+allometry["phi_s",spp]*out[,"diameter.ramet"]^allometry["theta_s",spp]
    out[out[,"status"]!="U","step_fecundity.prim"]<- exp(-allometry["mu_u",spp]*tstar*timestep)*exp(-allometry["mu_c",spp]*out[out[,"status"]!="U","time"])*allometry["F",spp]*allometry["g",spp]*out[out[,"status"]!="U","crown_area.prim"]*timestep
    out[out[,"status"]=="C","step_fecundity.ramet"]<- exp(-allometry["mu_u",spp]*tstar*timestep)*exp(-allometry["mu_cl",spp]*(tbreak-tstar)*timestep)*exp(-allometry["mu_c",spp]*out[out[,"status"]=="C","time"])*allometry["F",spp]*allometry["g",spp]*out[out[,"status"]=="C","crown_area.ramet"]*timestep
    out[,"cumsum.prim"]<-cumsum(out[,"step_fecundity.prim"])
    out[,"cumsum.ramet"]<-cumsum(out[,"step_fecundity.ramet"])
    }
  return(out)
}

##eventually this function will return a vector with all the disconnection times for a multi-stemmed plant if the primary stem did not die in nursing.
##But now it only supports the 2-stem root suckering tree.
#'TBREAK_CLONAL
#'function to calculate the number of timesteps when the disconnection between stems happen.
#'
#'@data ClonalLH a data frame of a clonal plant's life history of three statuses, U,N, and C. Can be generated using get.growth_clonal
#'
#'@return returns the developmentally programmed breaking time.
tbreak_clonal<-function(ClonalLH){
  max(as.integer(rownames(ClonalLH[which(ClonalLH$status=="N"),])))
}

##Below are a group of functions to calculate fecundity of the clonal tree during each of its life stage, understory/nursing/canopy depending on
##the time when the primary stem dies and the time when the secondary stem dies. They are by no means efficient, but I couldn't figure out
## a better way to do this besides looking through tables. Please let me know if you have any idea on making this cumbersome calculation of
##LRS easier. Also these functions might belong to methods.R by forming a new class? 

#'NURSING.PRIM
#'
#' function to calculate the fecundity of a clonal tree during its primary stem nursing period
#' 
#'@data ClonalLH a data frame of a clonal plant's life history of three statuses, U,N, and C. Can be generated using get.growth_clonal
#'@param tbreak a time point, in number of timesteps when the primary stem disconnects from the secondary trunk and stops nursing it with carbohydrates. It can be the time when the primary stem 
#'dies during nursing, the time when the secondary stem dies during nursing, YBD.
#'
#'@return returns the fecundity of a clonal tree during the nursing stage

Nursing.prim<-function(ClonalLH,tbreak){
  return(ClonalLH$cumsum.prim[tbreak])
}

#'SWITCHING.PRIM
#'
#' function to calculate the fecundity of a clonal tree's primary stem during the "switching" stage. If a ramet dies before reaching threshold size dstar, the primary stem stops nursing immediately and switches back to normal reproduction.
#' 
#'@data NonclonalLH a data frame of a nonclonal plant's life history
#'@data ClonalLH a data frame of a clonal plant's life history 
#'@param tbreak a time point, in number of timesteps, when the primary stem disconnects from the secondary trunk and stops nursing it with carbohydrates
#'@param YMD year mom dies---a time point, in number of timesteps, when the primary stem dies
#'
#'@return returns the fecundity of a clonal tree's primary stem during the switching stage

Switching.prim<-function(NonclonalLH,ClonalLH,tbreak,YMD){
 #find the diameter size of primary stem at tbreak in the clonal plant's life history data to switch a clonal plant back to non-clonal reproduction in the canopy
 breakt<-which(abs(NonclonalLH[,"diameter"]-ClonalLH[tbreak,"diameter.prim"])==min(abs(NonclonalLH[,"diameter"]-ClonalLH[tbreak,"diameter.prim"])))
  ifelse(YMD-breakt<2,return(0),return(NonclonalLH$cumsum[YMD-(tbreak-breakt)])) 
}

#'TBREAK_SWITCHING.RAMET
#'
#' function to calculate the expected disconnection time, tbreak, when a ramet is left alone to grow in the understory without support
#' 
#'@data NonclonalLH a data frame of a nonclonal plant's life history
#'@data ClonalLH a data frame of a clonal plant's life history 
#'@param YBD year baby dies---a time point, in number of timesteps, when the secondary stem dies
#'
#'@return returns the expected tbreak, actual time needed to reach the canopy
tbreak_switching.ramet<-function(NonclonalLH,ClonalLH,YMD,tstar=tstar){
  breakt<-which(abs(NonclonalLH[,"diameter"]-ClonalLH[YMD,"diameter.ramet"])==min(abs(NonclonalLH[,"diameter"]-ClonalLH[YMD,"diameter.ramet"])))
  return((tstar-breakt)+YMD)
}

#'SWITCHING.RAMET
#'
#' function to calculate the fecundity of a clonal tree's secondary stem during the "switching" stage. If a primary stem dies during nursing while the ramet survives into the canopy, its fecundity can be calculated.
#' 
#'@data NonclonalLH a data frame of a nonclonal plant's life history
#'@data ClonalLH a data frame of a clonal plant's life history 
#'@param tbreak a time point, in number of timesteps, when the primary stem disconnects from the secondary trunk and stops nursing it with carbohydrates
#'@param YBD year baby dies---a time point, in number of timesteps, when the secondary stem dies
#'
#'@return returns the fecundity of a clonal tree's ramet during the switching stage

Switching.ramet<-function(NonclonalLH,ClonalLH,YMD,YBD,tstar=tstar){
 #find the diameter size of ramet at break in the clonal plant's life history data to switch a clonal plant back to non-clonal understory
  breaktime<-tbreak_switching.ramet(NonclonalLH = NonclonalLH,ClonalLH = ClonalLH,YMD=YMD,tstar=tstar)
  ifelse(YBD-(tstar-breaktime)>nrow(NonclonalLH),return(NonclonalLH$cumsum[max_t]),return(NonclonalLH$cumsum[YBD-(tstar-breaktime)]))
}

#'SUCCESS.RAMET
#'
#' function to calculate the fecundity of a clonal tree's secondary stem if it successfully get nursed into an independent canopy tree. 
#' 
#'@data ClonalLH a data frame of a clonal plant's life history 
#'@param YBD year baby dies---a time point, in number of timesteps, when the secondary stem dies
#'
#'@return returns the fecundity of a clonal tree's ramet 

Success.ramet<-function(ClonalLH,YBD){
  return(ClonalLH$cumsum.ramet[YBD])
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
