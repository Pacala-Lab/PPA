


## This script is a sandbox for writing out demo's for team members to easily
## examine functionality. As an example I will write out a basic script to demonstrate
## how the light competition model I wrote works.

## First, load the package using:

devtools::load_all(##path/to/your/package/PPA)


## Then create an allometry object, which will contain all the parameters needed to simulate
## a stand with two species.
## The only argument you need to supply is the number of species. The
## stochastic argument, which is TRUE by default, automatically injects
## some variation into the allometries (right now in a very crude way).
## There is an additional parameter, custom, that can be used to provide
## custom allometries. Eventually, it could be nice to build in some set
## species allometries
test_allometry <- create_allometry.model_light(n.spp = 2, stochastic = TRUE)

## print test allometry to look inside. Note that its just a named matrix containing
## the relevant parameters.
test_allometry


## Next, we can run a simulation using the function: model_light.
test_simulation <- model_light(allometry = test_allometry[["allometry"]], ## supply the allometry object we just created
                               spp = c(1, 2), ## spp must be a vector of spp names that matches the names in the allometry object
                               max_t = 100, ## last timestep in simulation
                               timestep = 1, ## length of a timestep
                               plot.area = 1000, ## how big is the simulation area, in square meters
                               n_start = list("1" = 25, "2" = 25), ## how many individuals do we start with
                               full.data = TRUE, ## do we want the model output to include all the steps on the way?
                               full.allom = TRUE) ## dow we want the model output to include all of the allometric components? (i.e. biomass, height etc)

## now we can take a peak at the model_light output. As you can see it is a list of multiple datasets
## This is it in its most verbose form, as we specified both full.data and full.allom
str(test_simulation)


## now we can use the summary method to check out what we got:
summary(test_simulation)

## currently the summary.model_light method prints some basic attributes about the simulation as well as the average values
## for some of the dimensional attributes of the individuals. I don't necessarily think this is the best stuff to display, but
## I think deciding that should be a collaborative effort.

## finally we can plot the results of the our simulation using plot()
plot(test_simulation)

## This plot output plots the trajectories of the individual trees through time. Again, the plot output
## should probably be revised collaboratively
