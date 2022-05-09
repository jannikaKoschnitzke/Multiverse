library(tidyverse)
library(multiverse) 
library(mverse)
library(heplots) # for calculating effect sizes



# This code is based on data and code provided by the multiverse R-package:
# data: https://rdrr.io/github/MUCollective/multiverse/man/durante.html
# code: https://mucollective.github.io/multiverse/articles/visualising-multiverse.html

rm(list = ls()) 
data("durante") #read in data

## Prepare dataset ----------------------------------------------------------------

for (i in c("Abortion", "StemCell", "Marijuana", "RichTax", "StLiving", "Profit")){
  col_num = which(names(durante)== i)
  durante[,col_num] = abs(7 - durante[,col_num]) + 1
}

rawDat <- durante %>%
  mutate(
    RelComp = round((Rel1 + Rel2 + Rel3)/3, 2),
    FiscConsComp = (FreeMarket + PrivSocialSec + RichTax + StLiving + Profit)/5,
    SocConsComp = (Marriage + RestrictAbortion + Abortion + StemCell + Marijuana)/5
  )

# format dates
rawDat$StartDateNext <- as.Date(rawDat$StartDateNext, format = "%m/%d/%y")
rawDat$DateTesting <- as.Date(rawDat$DateTesting, format = "%m/%d/%y")
rawDat$StartDateofLastPeriod <- as.Date(rawDat$StartDateofLastPeriod, format = "%m/%d/%y")
rawDat$StartDateofPeriodBeforeLast <- as.Date(rawDat$StartDateofPeriodBeforeLast, format = "%m/%d/%y")



## Create Multiverse ---------------------------------------------------------------------

multiverse = create_multiverse(rawDat)

inside(multiverse, {
  data <- rawDat %>%
    mutate( ComputedCycleLength = StartDateofLastPeriod - StartDateofPeriodBeforeLast 
    ) %>%
    # 1. Exclusion of women based on cycle length (ECL)
    filter( branch(CycleLength, 
                   "cycle_length1" ~ TRUE,
                   "cycle_length2" ~ ComputedCycleLength > 25 & ComputedCycleLength < 35,
                   "cycle_length3" ~ ReportedCycleLength > 25 & ReportedCycleLength < 35
    )) %>%
    # 2. Next menstrual onset (NMO)
    # Depending on CycleLength, some combination of options where excluded due to inconsistency
    mutate( Menstruation = branch(menstrual_calculation, 
                                  "NextMenstrualOnset1" %when% (CycleLength != "cycle_length3") ~ StartDateofLastPeriod + ComputedCycleLength,
                                  "NextMenstrualOnset2" %when% (CycleLength != "cycle_length2") ~ StartDateofLastPeriod + ReportedCycleLength,
                                  "NextMenstrualOnset3" ~ StartDateNext
    )) %>%
    # 3. Assessment of relationship status (R) (single vs relationship)
    mutate(RelationshipStat = branch( relationship_status, 
                                      "Relationship1" ~ factor(ifelse(Relationship==1 | Relationship==2, 'Single', 'Relationship')),
                                      "Relationship2" ~ factor(ifelse(Relationship==1, 'Single', 'Relationship')),
                                      "Relationship3" ~ factor(ifelse(Relationship==1, 'Single', ifelse(Relationship==3 | Relationship==4, 'Relationship', NA))) )
    ) %>%
    mutate(
      CycleDay = 28 - (Menstruation - DateTesting),
      CycleDay = ifelse(CycleDay > 1 & CycleDay < 28, CycleDay, ifelse(CycleDay < 1, 1, 28))
    ) %>%
    # 4. Exclusion of women based on certainty ratings of start dates 
    filter( branch(Certainty,
                   "Certainty1" ~ TRUE,
                   "Certainty2" ~ Sure1 > 6 | Sure2 > 6
    )) %>%
    # 5. Assessment of fertility (F)—high vs low
    mutate( Fertility = branch( fertile,
                                "Fertility1" ~ factor( ifelse(CycleDay >= 7 & CycleDay <= 14, "high", ifelse(CycleDay >= 17 & CycleDay <= 25, "low", NA)) ),
                                "Fertility2" ~ factor( ifelse(CycleDay >= 6 & CycleDay <= 14, "high", ifelse(CycleDay >= 17 & CycleDay <= 27, "low", NA)) ),
                                "Fertility3" ~ factor( ifelse(CycleDay >= 9 & CycleDay <= 17, "high", ifelse(CycleDay >= 18 & CycleDay <= 25, "low", NA)) ),
                                "Fertility4" ~ factor( ifelse(CycleDay >= 8 & CycleDay <= 14, "high", "low") ),
                                "Fertility5" ~ factor( ifelse(CycleDay >= 8 & CycleDay <= 17, "high", "low") )
    ))
  
})


## Define Models -----------------------------------------------------------------

inside(multiverse, {
  
  ## Define models 
  # here: Anova
  mod_FCC <- lm(FiscConsComp~Fertility*RelationshipStat, data = data) 
  mod_SCC <- lm(SocConsComp~Fertility*RelationshipStat, data = data) 
  mod_RC  <- lm(RelComp~Fertility*RelationshipStat, data = data) 
  # here: logistic Regression 
  mod_D  <- glm(Donate~Fertility*RelationshipStat, data = data, family = binomial(link = "logit") )
  mod_V  <- glm(Vote~Fertility*RelationshipStat, data = data, family = binomial(link = "logit") )
  
  
  ## Calculating Effect Sizes (Eta2)
  eff_FCC <-  c(NA, etasq(mod_FCC)$`Partial eta^2`[-length(etasq(mod_FCC)$`Partial eta^2`)])
  eff_SCC <-  c(NA, etasq(mod_SCC)$`Partial eta^2`[-length(etasq(mod_SCC)$`Partial eta^2`)])
  eff_RC  <-  c(NA, etasq(mod_RC)$`Partial eta^2`[-length(etasq(mod_RC)$`Partial eta^2`)])
  eff_D   <-  c(NA, etasq(mod_D)$`Partial eta^2`[-length(etasq(mod_D)$`Partial eta^2`)])
  eff_V   <-  c(NA, etasq(mod_V)$`Partial eta^2`[-length(etasq(mod_V)$`Partial eta^2`)])
  
  ## save summary statistics in result object + effect sizes
  res_FCC <- mod_FCC  %>% 
    broom::tidy( conf.int = TRUE )
  res_FCC <- tibble(res_FCC, "effect_size" = eff_FCC)
  
  res_SCC <- mod_SCC  %>% 
    broom::tidy( conf.int = TRUE )
  res_SCC <- tibble(res_SCC, "effect_size" = eff_SCC)
  
  res_RC  <- mod_RC   %>% 
    broom::tidy( conf.int = TRUE )
  res_RC <- tibble(res_RC, "effect_size" = eff_RC)
  
  res_D   <- mod_D    %>% 
    broom::tidy( conf.int = TRUE )
  res_D <- tibble(res_D, "effect_size" = eff_D)
  
  res_V   <- mod_V    %>% 
    broom::tidy( conf.int = TRUE )
  res_V <- tibble(res_V, "effect_size" = eff_V)
})

multiverse = mverse::execute_multiverse(multiverse) 


#save to rda file
save(multiverse, file = "MV_computation/MultiVerse.rda")
