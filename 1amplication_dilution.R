rm(list=ls())
source("lab_paths.R")
setwd(local.path)
ncores <- 3

setwd("cascades-meadows/analysis/parasites")
source("src/misc.R")
source("src/writeResultsTable.R")
source("src/makeMultiLevelData.R")
source("src/runParasiteModels.R")
source("src/standardize_weights.R")
source("src/runPlotFreqModelDiagnostics.R")

## site or lat as the geographic variable
site.or.lat <- "lat"

## all of the variables that are explanatory variables and thus need
## to be centered
vars_yearsr <- c("VegAbundance",
                 "VegDiversity",
                 "SiteBeeDiversity",
                 "SiteBeeAbundance",
                 "SiteBombusAbundance",
                 "SiteHBAbundance"
                 )
vars_yearsrsp <- "rare.degree"

vars_site <- c("Proximity_150.200ha", "Proximity_750.1000ha", "Area")

variables.to.log <- c("rare.degree","SiteBeeAbundance", "Area", "VegAbundance",
                      "Proximity_150.200ha", "Proximity_750.1000ha")

## some zeros in data
variables.to.log.1 <- c("SiteHBAbundance", "SiteBombusAbundance")

## loads specimen data
source("src/init.R")


## because only Bombus and apis models converge, set the rest of
## the trait data to NA so that the variables scale properly
screened.bombus <- unique(spec.net$GenusSpecies[spec.net$Apidae == 1 &
                                                spec.net$Genus == "Bombus"])
screened.bombus <- screened.bombus[!is.na(screened.bombus)]

spec.net$rare.degree[!spec.net$GenusSpecies %in%
                     c("Apis mellifera", screened.bombus)] <- NA

dim(spec.net)
## raw, non standardized data for plotting
spec.orig <- prepDataSEM(spec.net, variables.to.log = NULL , variables.to.log.1,
                         standardize=FALSE)

## Make SEM weights and standardize data.
spec.net <- prepDataSEM(spec.net, variables.to.log = NULL, 
                        variables.to.log.1 = variables.to.log.1,
                        vars_yearsr = vars_yearsr, vars_sp = NULL,
                        vars_yearsrsp = vars_yearsrsp,
                        vars_site=vars_site)

spec.net$Site <- as.character(spec.net$Site)
## otherwise levels with no data are not properly dropped using subset
spec.net$Year <- as.character(spec.net$Year)
spec.net$GenusSpecies <- as.character(spec.net$GenusSpecies)

## bombus only data
spec.bombus <- makeGenusSubset(spec.net, "Bombus")
## apis only data
spec.apis <- makeGenusSubset(spec.net, "Apis")
## can repeat for other genera but have deleted since the models do
## not converge

## define all the formulas for the different parts of the models
source("src/plant_poll_models.R")

save(spec.net, spec.orig, file="saved/spec_weights.Rdata")

## Load phylogeny
load("../../data/community_phylogeny.Rdata")
## Species that are not in the phylogeny are not used. brms is not
## allowing an incomplete phylogeny, to avoid the error we changed the
## species not present to one that is in the phylogeny.  We chose a
## species for which we did not do parasite screening and should not
## influence results.
not_in_phylo <- unique(spec.net$GenusSpecies[!spec.net$GenusSpecies
                                             %in%
                                             phylo$tip.label])

## only bombus model is muti species given the others do not converge
spec.bombus$GenusSpecies[spec.bombus$GenusSpecies %in%
                         not_in_phylo]<- "Andrena topazana"



## **********************************************************
## community model, check assumptions first before adding parasites
## **********************************************************

run_plot_freq_model_diagnostics(remove_subset_formula(formula.flower.div),
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="students", site.lat=site.or.lat)

run_plot_freq_model_diagnostics(remove_subset_formula(formula.flower.abund),
                                 this_data=spec.net[spec.net$Weights == 1,],
                                 this_family="students",site.lat=site.or.lat)

run_plot_freq_model_diagnostics(remove_subset_formula(formula.bee.abund),
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="students", site.lat=site.or.lat)

run_plot_freq_model_diagnostics(remove_subset_formula(formula.bombus.abund),
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="students", site.lat=site.or.lat)

run_plot_freq_model_diagnostics(remove_subset_formula(formula.HB.abund),
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="students", site.lat=site.or.lat)

run_plot_freq_model_diagnostics(remove_subset_formula(formula.bee.div),
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="gaussian",
                                site.lat=site.or.lat)

## **********************************************************
## Parasite models set up
## **********************************************************
## Bombus predictor variables
## phylogeny must be last in all xvar sets

xvars.fd <-  c("VegDiversity",
               "YearPar",
               "SRDoy",
               "Area",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.fa <-  c("VegAbundance",
               "YearPar",
               "SRDoy",
               "Area",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.bd <-  c("SiteBeeDiversity",
               "YearPar",
               "SRDoy",
               "Area",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.ba <-  c("SiteBombusAbundance",
               "SRDoy", "YearPar",
               "Area",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")


xvars.ha <-  c("SiteHBAbundance",
               "YearPar",
               "SRDoy",
               "Area",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.d <-  c("rare.degree",
              "YearPar",
              "SRDoy",
              "Area",
              "(1|Site)",
              "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.p150 <-  c("Proximity_150.200ha",
              "(1|Site)",
              "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.p750 <-  c("Proximity_750.1000ha",
                 "(1|Site)",
                 "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.a <-  c("Area",
              "(1|Site)",
              "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.sr <-  c("SRDoy",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.y <-  c("Year",
              "(1|Site)",
              "(1|gr(GenusSpecies, cov = phylo_matrix))")


## **********************************************************
## Parasite presence
## **********************************************************
## Bombus
## **********************************************************

bombus.fd <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.fd,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="floral_div")

bombus.fa <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.fa,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="floral_abund")

bombus.bd <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.bd,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="bee_div")


bombus.ba <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.ba,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="bombus_abund")

bombus.ha <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.ha,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="hb_abund")

bombus.d <- runCombinedParasiteModels(spec.data= spec.bombus,
                                      species.group="bombus",
                                      xvars=xvars.d,
                                      ncores=ncores,
                                      data2= list(phylo_matrix=phylo_matrix),
                                      xvar.name="degree")

bombus.l <- runCombinedParasiteModels(spec.data= spec.bombus,
                                      species.group="bombus",
                                      xvars=xvars.l,
                                      ncores=ncores,
                                      data2= list(phylo_matrix=phylo_matrix),
                                      xvar.name="lat")

bombus.a <- runCombinedParasiteModels(spec.data= spec.bombus,
                                      species.group="bombus",
                                      xvars=xvars.a,
                                      ncores=ncores,
                                      data2= list(phylo_matrix=phylo_matrix),
                                      xvar.name="area")

bombus.sr <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.sr,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="SRDOY")

bombus.y <- runCombinedParasiteModels(spec.data= spec.bombus,
                                      species.group="bombus",
                                      xvars=xvars.y,
                                      ncores=ncores,
                                      data2= list(phylo_matrix=phylo_matrix),
                                      xvar.name="year")

bombus.p150 <- runCombinedParasiteModels(spec.data= spec.bombus,
                                      species.group="bombus",
                                      xvars=xvars.y,
                                      ncores=ncores,
                                      data2= list(phylo_matrix=phylo_matrix),
                                      xvar.name="proximity150")

bombus.p750 <- runCombinedParasiteModels(spec.data= spec.bombus,
                                      species.group="bombus",
                                      xvars=xvars.y,
                                      ncores=ncores,
                                      data2= list(phylo_matrix=phylo_matrix),
                                      xvar.name="proximity750")
