##devtools::load_all("/Users/aiyuzheng/PPA/")
## This file contain code for the function model_clonal.

##model_clonal is not an explicit simulator at present. It assumes infinite space. It's for simulating the LRS for a clonal tree when the mother/sucker are killed randomly.
##but it outputs growth rates, allometries, and lifetime reproductive sucess. I hope we can eventually incorporate clonal allometries into model_light
##or model_water, but I am still working on ranking clonal crowns. So how far we want to go will depend on collective decision. I will put in the codes for invasion analysis and visualization
##after my general exam :)

#' This function simulates a clonal tree invading a population of non-clonal trees at equilibrium or vice versa using the perfect plasticity approximation for light competition.
#'
#'@param allometry an object of class: clonal or non-clonal allometries created using the create_allometry method 
#'@param max_t the last timestep in the simulation
#'@param timestep the length of each timestep
#'@param N the initial population size of the clonal trees' understory seedlings
#'@param dnot the initial diameter for each tree seedling
#'@param dnot_ramet the initial size of a ramet's diameter 
#'@param nstem the number of stems a clonal plant develops 
#'@param full.lifehistory TRUE/FALSE specifying whether you would like the model output to include life histories of a clonal tree and a non-clonal tree
#'@param recip TRUE/FALSE specifying the direction of invasion. TRUE indicates the clonal tree being the resident, and FALSE indicates the clonal tree being the invader.
#'
#'@return a table containing life time reproductive success for each simulated tree or the final average LRS 


## simulator for light competition model.
model_clonal <- function(allometry = create_allometry.model_clonal(n.spp = 1,stochastic = FALSE),
                        spp,
                        max_t,
                        timestep,
                        N,
                        dnot,
                        dnot_ramet,
                        nstem=2,
                        full.lifehistory = TRUE,
                        recip = FALSE) {
  ## stop on error
  stopifnot(is.numeric(max_t),
            is.numeric(timestep),
            is.numeric(dnot),
            is.numeric(dnot_ramet),
            is.numeric(nstem),
            timestep <= max_t,
            is.logical(full.lifehistory),
            is.logical(recip))
  
  ## ensure spp name is character not numeric
  spp <- as.character(spp)
  
  ##generate life history data for a non-clonal tree and a clonal tree
  
  NonclonalLH<-get.growth_nonclonal(dnot=dnot,allometry=allometry,spp=spp,max_t=max_t,timestep=timestep)

  ClonalLH<-get.growth_clonal(dnot=dnot,dnot_ramet=dnot_ramet,allometry = allometry,spp=spp,max_t=max_t,timestep=timestep)
  
  ##get tstar
  tstar<-equil_nonclonal(dnot=dnot,allometry=allometry,spp = spp,max_t=max_t,timestep = timestep)[["tstar"]]
  ##get tbreak
  tbreak<-tbreak_clonal(ClonalLH=ClonalLH)
  
  ##N trees to start with
  
  ####kill some trees in the understory 
  N_understory<-rexp(N, rate = allometry["mu_u",spp])
  N_canopy<-round(N_understory[N_understory>(tstar*timestep)],digits=1)
  
  ####assign years of death to a clonal tree's primary stem and secondary stem in the canopy
  Death<-matrix(ncol = 4,nrow = length(N_canopy))
  colnames(Death)<-c("YMD","YBD","theme","LRS")
  rownames(Death)<-as.character(seq(1,length(N_canopy)),by=1)
  Death[,"YMD"]<-tstar+round(rexp(length(N_canopy), rate = allometry["mu_c",spp]),digits = 1)/timestep
  Death[,"YBD"]<-tstar+round(rexp(length(N_canopy), rate = allometry["mu_cl",spp]),digits=1)/timestep
 # Death[Death[,"YBD"]>tbreak*timestep]<-tbreak+round(rexp(length(Death[Death[,"YBD"]>tbreak*timestep]), rate = allometry["mu_c",spp]),digits = 1)/timestep
  Death[Death > max_t] <- max_t ##make sure trees don't live longer than max_t
  Death[,"theme"]<-NA
  Death[,"LRS"]<-NA
  
  ##based on data frame death, sum up the fecundity for each tree 
  for (i in 1:nrow(Death)){
    if(Death[i,"YMD"]<=tbreak){
      if(Death[i,"YBD"]<=tbreak){
        if(Death[i,"YMD"]>Death[i,"YBD"]){
          Death[i,"theme"]<-1
          Death[i,"LRS"]<-Nursing.prim(ClonalLH = ClonalLH,tbreak=Death[i,"YBD"])+Switching.prim(NonclonalLH = NonclonalLH, ClonalLH=ClonalLH, tbreak = Death[i,"YBD"], YMD=Death[i,"YMD"])
        }else{
          Death[i,"theme"]<-2
          Death[i,"LRS"]<-Nursing.prim(ClonalLH = ClonalLH,tbreak=Death[i,"YMD"])
        }
      }else{
        if(tbreak_switching.ramet(NonclonalLH = NonclonalLH,ClonalLH = ClonalLH, YMD=Death[i,"YMD"],tstar=tstar)>=Death[i,"YBD"]){
          Death[i,"theme"]<-2
          Death[i,"LRS"]<-Nursing.prim(ClonalLH = ClonalLH,tbreak=Death[i,"YMD"])
        }else{
          Death[i,"theme"]<-4
          Death[i,"LRS"]<-Nursing.prim(ClonalLH = ClonalLH,tbreak=Death[i,"YMD"])+Switching.ramet(NonclonalLH = NonclonalLH,ClonalLH = ClonalLH,YMD = Death[i,"YMD"],YBD=Death[i,"YBD"],tstar=tstar)
        }
      }
    }else{
      if(Death[i,"YBD"]<=tbreak){
        Death[i,"theme"]<-3 
        Death[i,"LRS"]<-Nursing.prim(ClonalLH = ClonalLH,tbreak=Death[i,"YBD"])+Switching.prim(NonclonalLH = NonclonalLH, ClonalLH= ClonalLH,tbreak = Death[i,"YBD"], YMD=Death[i,"YMD"])
      }else{
        Death[i,"theme"]<-5 
        Death[i,"LRS"]<-Nursing.prim(ClonalLH = ClonalLH,tbreak=Death[i,"YMD"])+Success.ramet(ClonalLH = ClonalLH,YBD=Death[i,"YBD"])
      }
    }
   }
  
  ## if user specifies that full data should be reported in single matrix, initialize and populate matrix
  if(full.lifehistory){
    outlist<-list("clonal_lifehistory"=ClonalLH,"non-clonal_lifehistory"=NonclonalLH,"simulated_lifehistory"=Death)
    class(outlist) <- "model_clonal"
    return(outlist)
  }else{
    output <- Death
    class(output) <- "model_clonal"
    return(output) 
  }
}




