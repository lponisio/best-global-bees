## **********************************************************
## Load libraries
## **********************************************************

## this project is so fun! woah!


rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path
setwd(local.path)
setwd("pnw_survey/analyses/microbes")

## Toggle on or off to run individual genus models
run.bombus = TRUE

library(picante)
library(plyr)
library(bayesplot)
library(pscl)
library(brms)
library(performance)
library(R2admb)
#library(shinystan)
## If not already installed
# install.packages("glmmADMB", 
#                  repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#                          getOption("repos")),
#                  type="source")


## **********************************************************
## Standardize, center, and transform data
## **********************************************************

## all of the variables that are explanatory variables and thus need
## to be centered

variables.to.log <- c("ContExtractedSeverity",
                      #"HighSeverity", 
                      "HighSeverityArea",
                      "ForageDist_km",
                      "BombusRichness"
                      #"FlowerRareRichness"
)

variables.to.log.1 <- c("BombusAbundance",
                        "FlowerRareRichness",
                        "VegAbundance",
                        "BombusDiversity",
                        "BeeAbundance",
                        "BeeDiversity",
                        "VegDiversity")

vars_yearsr <- c("VisitedFloralDiversity",
                 "FlowerRareRichness",
                 "BeeAbundance",
                 "BeeDiversity",
                 "VegAbundance",
                 "BombusAbundance",
                 "BombusRichness",
                 "VegDiversity"
)

vars_yearsrsp <- "rare.degree"

## not ready yet
#vars_sp <- "MeanITD"

# Continuous fire variables
vars_site <- c("HighSeverityArea",
               #"HighSeverity",
               #"FireSeverity",
               #"CatExtractedSeverity",
               "ContExtractedSeverity"
               )

## **********************************************************
## Load Data
## **********************************************************

load("../../data/spec_RBCL_16s.Rdata")

ncores <- 1

## **********************************************************
## Source files
## **********************************************************

source("src/misc_microbe.R")
#source("src/misc.R")
source('src/makeMultiLevelData.R')
source("src/standardize_weights_microbes.R") #TODO fix title
source("src/init_microbe.R")
source("src/writeResultsTable.R")
#source("src/runPlotFreqModelDiagnostics.R")


ncores <- 1

spec.net <- spec.net %>%
  filter(FireName == "Holiday") %>%
  mutate(FireSeverity = ifelse(FireSeverity == "", "Unburned", FireSeverity)) %>%
  mutate(BinnedSeverity = ifelse(CatExtractedSeverity == 1 | CatExtractedSeverity == 2, 0, 1))

## Binning severity by unburned/low and mod/high

## **********************************************************
## Flower diversity and abundance
## **********************************************************

flower.div.vars <- c("HighSeverityArea*BinnedSeverity",
                      "DoyStart",
                      "I(DoyStart^2)",
                      "Year",
                      "(1|Stand)"
                      )

flower.div.x <- paste(flower.div.vars, collapse="+")
flower.div.y <- "VegDiversity | subset(Weights)"
formula.flower.div <- as.formula(paste(flower.div.y, "~",
                                        flower.div.x))

flower.abund.vars <- c("HighSeverityArea*BinnedSeverity",
                     "DoyStart",
                     "I(DoyStart^2)",
                     "Year",
                     "(1|Stand)"
)

flower.abund.x <- paste(flower.abund.vars, collapse="+")
flower.abund.y <- "VegAbundance | subset(Weights)"
formula.flower.abund <- as.formula(paste(flower.abund.y, "~",
                                       flower.abund.x))

## **********************************************************
## Bee abundance
## **********************************************************

bombus.abund.vars <- c("VegAbundance",
                       "TempCStart",
                       "I(TempCStart^2)",
                       "HighSeverityArea*BinnedSeverity",
                       "Year",
                       "(1|Stand)")

bombus.abund.x <- paste(bombus.abund.vars, collapse="+")
bombus.abund.y <- "BombusAbundance | subset(Weights)"
formula.bombus.abund <- as.formula(paste(bombus.abund.y, "~",
                                         bombus.abund.x))

tot.bee.abund.vars <- c("VegAbundance",
  "TempCStart",
  "I(TempCStart^2)",
  "HighSeverityArea*BinnedSeverity",
  "Year",
  "(1|Stand)")

tot.bee.abund.x <- paste(tot.bee.abund.vars, collapse="+")
tot.bee.abund.y <- "BeeAbundance | subset(Weights)"
formula.tot.bee.abund <- as.formula(paste(tot.bee.abund.y, "~",tot.bee.abund.x))

## **********************************************************
## Bee diversity
## **********************************************************

## bombus richness
bombus.div.vars <- c("VegDiversity",
                      "TempCStart",
                      "I(TempCStart^2)",
                     "HighSeverityArea*BinnedSeverity",
                      "Year",
                      "(1|Stand)")

bombus.div.x <- paste(bombus.div.vars, collapse="+")
bombus.div.y <- "BombusDiversity | subset(Weights)"
formula.bombus.div <- as.formula(paste(bombus.div.y, "~",
                                        bombus.div.x))



tot.bee.div.vars <-  c("VegDiversity",
                       "TempCStart",
                       "I(TempCStart^2)",
                       "HighSeverityArea*BinnedSeverity",
                       "Year",
                       "(1|Stand)")

tot.bee.div.x <- paste(tot.bee.div.vars, collapse="+")
tot.bee.div.y <- "BeeDiversity | subset(Weights)"
formula.tot.bee.div <- as.formula(paste(tot.bee.div.y, "~",tot.bee.div.x))

## **********************************************************
## convert formulas to brms format
## **********************************************************
bf.fdiv <- bf(formula.flower.div)
bf.fabund <- bf(formula.flower.abund)
bf.tot.babund <- bf(formula.tot.bee.abund)
bf.tot.bdiv <- bf(formula.tot.bee.div)
bf.bom.div <- bf(formula.bombus.div)
bf.bom.abund <- bf(formula.bombus.abund)

## **********************************************************
## change NAs to 0 to prevent brms dropping missing data
## **********************************************************

spec.net[is.na(spec.net)] <- 0

spec.bombus <- spec.net[spec.net$Genus == "Bombus",]

spec.bombus$CatExtractedSeverity <- as.factor(spec.bombus$CatExtractedSeverity)

## **********************************************************
## Bombus microbe PD model
## **********************************************************

## all microbes
microbe.bombus.vars <- c("BeeAbundance",
                            "BeeDiversity",
                            "VegAbundance",
                            "VegDiversity",
                         "HighSeverityArea*BinnedSeverity",
                            "rare.degree", 
                            "(1|Stand)",
                            "(1|gr(GenusSpecies, cov = phylo_matrix))")



microbe.bombus.x <- paste(microbe.bombus.vars, collapse="+")
microbe.bombus.y <- "PD | subset(WeightsMicrobe)"
formula.microbe.bombus <- as.formula(paste(microbe.bombus.y, "~",
                                              microbe.bombus.x))

bf.microbe.bombus <- bf(formula.microbe.bombus)



## Obligate PD model
ob.microbe.bombus.vars <- c(#"BeeAbundance",
                            #"BeeDiversity",
                            "BombusDiversity",
                            "BombusAbundance",
                            "VegAbundance",
                            "VegDiversity",
                            "HighSeverityArea*BinnedSeverity",
                            "rare.degree", 
                            "(1|Stand)",
                            "(1|gr(GenusSpecies, cov = phylo_matrix))")


ob.microbe.bombus.x <- paste(ob.microbe.bombus.vars, collapse="+")
## distribution looks better unlogged
ob.microbe.bombus.y <- "PD.obligate | subset(WeightsObligateBombus)"
formula.ob.microbe.bombus <- as.formula(paste(ob.microbe.bombus.y, "~",
                                           ob.microbe.bombus.x))

bf.ob.microbe.bombus <- bf(formula.ob.microbe.bombus)

## Transient PD model
tr.microbe.bombus.vars <- c(#"BeeAbundance",
                            #"BeeDiversity",
                            "BombusDiversity",
                            "BombusAbundance",
                            "VegAbundance",
                            "VegDiversity",
                            "HighSeverityArea*BinnedSeverity",
                            "rare.degree", 
                            "(1|Stand)",
                            "(1|gr(GenusSpecies, cov = phylo_matrix))")


tr.microbe.bombus.x <- paste(tr.microbe.bombus.vars, collapse="+")
tr.microbe.bombus.y <- "PD.transient.log | subset(WeightsTransientBombus)"
formula.tr.microbe.bombus <- as.formula(paste(tr.microbe.bombus.y, "~",
                                              tr.microbe.bombus.x))

bf.tr.microbe.bombus <- bf(formula.tr.microbe.bombus)

## Facultative model
# non.ob.microbe.bombus.vars <- c("BeeAbundance",
#                             "BeeDiversity",
#                             "Lat", 
#                             "MeanFloralDiversity",
#                             "MeanITD",
#                             "rare.degree",
#                             "(1|Site)",
#                             "(1|gr(GenusSpecies, cov = phylo_matrix))") 
# 
# non.ob.microbe.bombus.x <- paste(non.ob.microbe.bombus.vars, collapse="+")
# non.ob.microbe.bombus.y <- "PD.transient.log | subset(WeightsTransientBombus)"
# formula.non.ob.microbe.bombus <- as.formula(paste(non.ob.microbe.bombus.y, "~",
#                                               non.ob.microbe.bombus.x))
# 
# bf.non.ob.microbe.bombus.student <- bf(formula.non.ob.microbe.bombus, family=student())
# 
# Chain 1:   Error evaluating the log probability at the initial value.
# Chain 1: Exception: normal_id_glm_lpdf: Vector of dependent variables[42] is -inf, but must be finite! (in 'string', line 53, column 4 to column 55)
# [1] "Error : Initialization failed."
# [1] "error occurred during calling the sampler; sampling not done"
# Error in h(simpleError(msg, call)) : 
#   error in evaluating the argument 'x' in selecting a method for function 'as.data.frame': non-numeric argument to mathematical function
# Called from: h(simpleError(msg, call))

## combined model
bform.bombus <- bf.ob.microbe.bombus + bf.tr.microbe.bombus +
  bf.fdiv +
  bf.fabund +
  bf.bom.abund +
  bf.bom.div +
  set_rescor(FALSE)


if(run.bombus){
  fit.microbe.bombus <- brm(bform.bombus, spec.bombus,
                          cores=ncores,
                          iter = 10000,
                          chains =1,
                          thin=1,
                          init=0,
                          open_progress = FALSE,
                          control = list(adapt_delta = 0.99,
                                         max_treedepth = 15
                                         ),
                          save_pars = save_pars(all = TRUE),
                          data2 = list(phylo_matrix=phylo_matrix))

  write.ms.table(fit.microbe.bombus, "bombus_microbe_strongweak_full_binned")
  r2loo.bombus <- loo_R2(fit.microbe.bombus)
  r2.bombus <- rstantools::bayes_R2(fit.microbe.bombus)
  save(fit.microbe.bombus, spec.net, r2.bombus, r2loo.bombus,
       file="saved/fullMicrobeBombusStrongWeakFit_full_binned.Rdata")
}

## aug == using bee abund/div (only pareto warning)
## aug2 == using bombus abund/div (only pareto warning too -- keep this one)


## initial model results (obligate PD negative relationship to bee richness, transient pd positive relationship to bee abundance)

run.bombus.all = TRUE
if(run.bombus.all){
  fit.microbe.bombus <- brm(bf.microbe.bombus, spec.bombus,
                            cores=ncores,
                            iter = 10000,
                            chains =1,
                            thin=1,
                            init=0,
                            open_progress = FALSE,
                            control = list(adapt_delta = 0.99,
                                           max_treedepth = 15
                            ),
                            save_pars = save_pars(all = TRUE),
                            data2 = list(phylo_matrix=phylo_matrix))
  
  write.ms.table(fit.microbe.bombus, "bombus_microbe_allPD")
  r2loo.bombus <- loo_R2(fit.microbe.bombus)
  r2.bombus <- rstantools::bayes_R2(fit.microbe.bombus)
  save(fit.microbe.bombus, spec.net, r2.bombus, r2loo.bombus,
       file="saved/fullMicrobeBombusAllPDFit.Rdata")
}

# prelim run warnings: 1: There were 4 divergent transitions after warmup.

## TODO next steps -- run other model layers combined model
## also make models for evenness
## model diagnostics -- also examine model fams used for sky islands modeling
