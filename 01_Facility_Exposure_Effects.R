## Maggie Lind
## Generate the matched dataset for the close contact analysis
  ## matching Cell and Dorm residents based on exposure (known close contact) on facility, room type, vaccination status 

library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(hrbrthemes)
library(ggridges)
library(ggsci)
library(survival)
library(forestplot)
library(AF)
library(survminer)
library(multcomp)
library(splines)

## Load data 
survdata <- read.csv(LOAD_DATA.csv)

## subset to the variant of interest (Delta or Omicron in manuscript) (Delta used as exampled here)
survdata <- subset(survdata, variant == "Delta")

## format the survival data 
survdata$exposure <- as.factor(survdata$exposure)
survdata$stfactor <- paste(survdata$Housing.Unit.Code, survdata$LOC_FACILITY_CODE)
survdata$clusterid <- paste(survdata$Housing.Unit.Code, survdata$contact_date)
  
## Immunity Event Analysis: INFECTIOUSNESS ########################################################
## VACCINATION #### 
  ## run model 
    res.cox <- coxph(Surv(time, outcome) ~ exposure + ns(Age, df = 4) + vax_status + priorinf +
                       Race + ns(ymd(start_time), df = 4) + roomsize + ns(unitsize, df = 4) + 
                       strata(LOC_FACILITY_CODE) + cluster(clusterid), data =  survdata, method = "breslow")
    summary(res.cox)
    r <- summary(res.cox)
  
    ## KM plots data
      sfit <- survfit(Surv(time, outcome) ~ exposure , data =  survdata)
    
    ## View results
      ggsurvplot(sfit, data = survdata, risk.table = TRUE, conf.int = TRUE, palette = c("#2C628D", "#DB9F37", "#982932"),
               legend.title = "", legend.labs = c("Unexposed", "Within Cell Contact", "Within Unit Contact"),
               ylim = c(0.75,1), xlim = c(0, 14), break.x.by = 2, tables.theme = theme_cleantable(), risk.table.title = "Risk Table (Absolute Number)", 
           title = "Kaplan Meier Curve Comparing Survival Between Residents with and without \nClose Contact Exposure (Survival defined as absence of SARS-CoV-2 infection)", xlab = "Time since exposure (days)", ylab = "Survival") 
    
  ## SAVE THE HR 
    t <- as.data.frame(summary(res.cox)$coefficients)
    t$names <- rownames(t)

    ## Format the saved results
      rt <- data.frame(exposure = t$names[1:2], fit = t$coef[1:2], se = t$`robust se`[1:2], point = t$`exp(coef)`[1:2], lower = r$conf.int[1:2,3], 
                       upper = r$conf.int[1:2,4]) %>% mutate(seed = iteration + seed, analysis = "primary", group = "none", pvalue = t$`Pr(>|z|)`[1:2])
      
      rt$number_exposed_0 <- nrow(survdata[which(survdata$exposure == 0),])
      rt$number_exposed_1 <- nrow(survdata[which(survdata$exposure == 1),])
      rt$number_exposed_2 <- nrow(survdata[which(survdata$exposure == 2),])
      rt$number_exposed_0_cases <- nrow(survdata[which(survdata$exposure == 0 & survdata$outcome == 1),])
      rt$number_exposed_1_cases <- nrow(survdata[which(survdata$exposure == 1 & survdata$outcome == 1),])
      rt$number_exposed_2_cases <- nrow(survdata[which(survdata$exposure == 2 & survdata$outcome == 1),])
      rt$mean_time_0 <- mean(survdata$time[which(survdata$exposure == 0)])
      rt$mean_time_1 <- mean(survdata$time[which(survdata$exposure == 1)])
      rt$mean_time_2 <- mean(survdata$time[which(survdata$exposure == 2)])
      rt$sd_time_0 <- sd(survdata$time[which(survdata$exposure == 0)])
      rt$sd_time_1 <- sd(survdata$time[which(survdata$exposure == 1)])
      rt$sd_time_2 <- sd(survdata$time[which(survdata$exposure == 2)])
      rt$exposure <- gsub("as.factor[(]exposure[)]", "", rt$exposure)
    
  ## get the scholfelds
    test.ph <- cox.zph(res.cox)$table

## save results 
  write.csv(rt, filename) ## save results
  write.csv(test.ph, filename) ## save residuals
  