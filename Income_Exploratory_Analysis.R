##############################################################################
## Title: UKBB Income Exploratory Analysis
## Author: Regina Manansala
## Date Created: 08-October-2018
## Date Modified: 24-October-2018
##############################################################################

library(data.table)
library(dplyr)
library(tidyr)
library(MASS)

## Import UKBB Income data
source('ukb23544.r')
ukbb_inc <- bd

## Import Main UKBB data
# source('ukb21755.r')

## Match variables between income and main data
# inc_names <- as.data.frame(names(ukbb_inc))
# names(inc_names)[1] <- "Var"
# bd_names <- as.data.frame(names(bd))
# names(bd_names)[1] <- "Var"
# dup_var_names <- inner_join(inc_names, bd_names, by="Var")
# dup_var_names <- as.list(dup_var_names)

## Rename Variables (Income Data)
names(ukbb_inc)[names(ukbb_inc) == "f.eid"] <- "ID"
names(ukbb_inc)[names(ukbb_inc) == "f.738.0.0"] <- "income"
names(ukbb_inc)[names(ukbb_inc) == "f.31.0.0"] <- "sex" 
names(ukbb_inc)[names(ukbb_inc) == "f.54.0.0"] <- "location"
names(ukbb_inc)[names(ukbb_inc) == "f.34.0.0"] <- "birth_yr"
names(ukbb_inc)[names(ukbb_inc) == "f.52.0.0"] <- "birth_mo"
names(ukbb_inc)[names(ukbb_inc) == "f.6142.0.0"] <- "empl_status"
names(ukbb_inc)[names(ukbb_inc) == "f.21001.0.0"] <- "BMI"

## Create disease variables
ukbb_inc$diabetes <- ifelse(is.na(ukbb_inc$f.2443.0.0) |ukbb_inc$f.2443.0.0 == "No" | ukbb_inc$f.2443.0.0 == "Prefer not to answer" | ukbb_inc$f.2443.0.0 == "Do not know", 0, 1) #diabetes
ukbb_inc$asthma <- ifelse(ukbb_inc$f.3786.0.0 == -3 | ukbb_inc$f.3786.0.0 == -1 | is.na(ukbb_inc$f.3786.0.0), 0, 1)
ukbb_inc$mi_flag <- ifelse(is.na(ukbb_inc$f.42000.0.0), 0, 1) #myocardial infarction
ukbb_inc$cancer_dx_dt <- ifelse(is.na(ukbb_inc$f.40005.0.0), 0, 1) #cancer diagnosis date
ukbb_inc$cancer_sr <- ifelse(is.na(ukbb_inc$f.20001.0.0), 0, 1) #self-reported cancer

## Recategorize race/ethnicity variable
ukbb_inc$race[ukbb_inc$f.21000.0.0 == "British"] <- "British"
ukbb_inc$race[ukbb_inc$f.21000.0.0 == "Irish"] <- "Irish"
ukbb_inc$race[ukbb_inc$f.21000.0.0 == "Any other white background" | ukbb_inc$f.21000.0.0 == "White"] <- "Other White"
ukbb_inc$race[ukbb_inc$f.21000.0.0 == "Chinese"] <- "Chinese"
ukbb_inc$race[ukbb_inc$f.21000.0.0 == "Caribbean" | ukbb_inc$f.21000.0.0 == "African" | ukbb_inc$f.21000.0.0 == "Any other Black background" | ukbb_inc$f.21000.0.0 == "Black or Black British"] <- "Black"
ukbb_inc$race[ukbb_inc$f.21000.0.0 == "Mixed" | ukbb_inc$f.21000.0.0 == "White and Black Caribbean" | ukbb_inc$f.21000.0.0 == "White and Black African" | ukbb_inc$f.21000.0.0 == "White and Asian" | ukbb_inc$f.21000.0.0 == "Any other mixed background" | ukbb_inc$f.21000.0.0 == "Other ethnic group" | ukbb_inc$f.21000.0.0 == "Do not know" | ukbb_inc$f.21000.0.0 == "Prefer not to answer" | ukbb_inc$f.21000.0.0 == "Any other Asian background" | ukbb_inc$f.21000.0.0 == "Asian or Asian British"] <- "Other"
ukbb_inc$race[ukbb_inc$f.21000.0.0 == "Indian" | ukbb_inc$f.21000.0.0 == "Pakistani" | ukbb_inc$f.21000.0.0 == "Bangladeshi"] <- "South Asian"

## Drop unnecessary variables
ukbb_inc_sub <- dplyr::select(ukbb_inc, -starts_with("f."))

## Remove NA income
inc_analysis <- subset(ukbb_inc_sub, !is.na(ukbb_inc_sub$income) & ukbb_inc_sub$income != "Prefer not to answer" & ukbb_inc_sub$income != "Do not know")
inc_analysis <- inc_analysis[inc_analysis$empl_status != "Retired" & inc_analysis$empl_status != "Full or part-time student" & inc_analysis$empl_status != "Doing unpaid or voluntary work"
 & inc_analysis$empl_status != "None of the above" & inc_analysis$empl_status != "Prefer not to answer", ]
inc_analysis$income <- factor(inc_analysis$income)
inc_analysis$empl_status <- factor(inc_analysis$empl_status)

## Tabulate variables
table(inc_analysis$income)
table(inc_analysis$empl_status)

table(inc_analysis$income, inc_analysis$race, exclude=NULL)
prop.table(table(inc_analysis$income, inc_analysis$race, exclude=NULL))
table(inc_analysis$income, inc_analysis$sex, exclude=NULL)
prop.table(table(inc_analysis$income, inc_analysis$sex, exclude=NULL))
table(inc_analysis$income, inc_analysis$diabetes, exclude=NULL)
prop.table(table(inc_analysis$income, inc_analysis$diabetes, exclude=NULL))
table(inc_analysis$income, inc_analysis$asthma, exclude=NULL)
prop.table(table(inc_analysis$income, inc_analysis$asthma, exclude=NULL))
table(inc_analysis$income, inc_analysis$mi_flag, exclude=NULL)
prop.table(table(inc_analysis$income, inc_analysis$mi_flag, exclude=NULL))
table(inc_analysis$income, inc_analysis$cancer_dx_dt, exclude=NULL)
prop.table(table(inc_analysis$income, inc_analysis$cancer_dx_dt, exclude=NULL))
table(inc_analysis$income, inc_analysis$cancer_sr, exclude=NULL)
prop.table(table(inc_analysis$income, inc_analysis$cancer_sr, exclude=NULL))

table(inc_analysis$income, inc_analysis$empl_status, exclude=NULL)
prop.table(table(inc_analysis$income, inc_analysis$empl_status, exclude=NULL))

table(inc_analysis$income, inc_analysis$location, exclude=NULL)
prop.table(table(inc_analysis$income, inc_analysis$location, exclude=NULL))

table(inc_analysis$income, inc_analysis$birth_mo, exclude=NULL)
prop.table(table(inc_analysis$income, inc_analysis$birth_mo, exclude=NULL))


## Subsets for Continuous Variables
# inc_bmi <- inc_analysis[, c("income", "BMI")]
# inc_bmi %>% group_by(income) %>% summarize(mean=mean(BMI, na.rm = TRUE), sd=sd(BMI, na.rm = TRUE),min=min(BMI, na.rm = TRUE), max=max(BMI, na.rm = TRUE), n=n())
# write.table(inc_bmi, "inc_bmi.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)
# 
# inc_birthyr <- inc_analysis[, c("income", "birth_yr", "empl_status")]
# write.table(inc_birthyr, "inc_birthyr.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)
# 
# inc_birthmo <- inc_analysis[, c("income", "birth_mo")]
# write.table(inc_birthmo, "inc_birthmo.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)
# 
# inc_loc <- inc_analysis[, c("income", "location")]
# write.table(inc_loc, "inc_loc.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)

## Plots
# 
# inc_bmi$income = factor(inc_bmi$income, levels=c('Less than 18,000','18,000 to 30,999','31,000 to 51,999','52,000 to 100,000', 'Greater than 100,000'))
# inc_birthyr$income = factor(inc_birthyr$income, levels=c('Less than 18,000','18,000 to 30,999','31,000 to 51,999','52,000 to 100,000', 'Greater than 100,000'))
# inc_birthmo$income = factor(inc_birthmo$income, levels=c('Less than 18,000','18,000 to 30,999','31,000 to 51,999','52,000 to 100,000', 'Greater than 100,000'))
# inc_loc$income = factor(inc_loc$income, levels=c('Less than 18,000','18,000 to 30,999','31,000 to 51,999','52,000 to 100,000', 'Greater than 100,000'))
# 
# table(inc_birthyr$income)
# head(inc_birthyr)
# 
# foo <- as.data.frame(table(inc_birthyr$income, inc_birthyr$birth_yr))
# names(foo) <- c("income", "birth_yr", "num")
# foo$income <- factor(foo$income, levels=c("Less than 18,000","18,000 to 30,999","31,000 to 51,999","52,000 to 100,000","Greater than 100,000"))
# 
# 
# ggplot(inc_bmi, aes(x=BMI))+geom_histogram(binwidth = 1, colour="black", fill="white")+
#   scale_x_continuous(breaks=unique(inc_loc$location))+
#   theme(axis.text.x = element_text(face = "italic", size = 8, angle = 90, vjust = 0.5))+
#   facet_wrap(~inc_bmi$income)
# 
# ggplot(foo, aes(x=birth_yr, y=num, color=income)) + geom_line(aes(group = income)) +theme(axis.text.x = element_text(face = "italic", size = 8, angle = 90, vjust = 0.5))
# 
# ggplot(inc_birthyr, aes(x=birth_yr))+geom_bar(aes(y = (..count..)/sum(..count..)))+facet_grid(~inc_birthyr$income) 
# + scale_y_continuous(labels = scales::percent)
# 
# ggplot(inc_birthyr, aes(x=birth_yr))+geom_bar()+facet_grid(~inc_birthyr$income)
# 
# ggplot(inc_birthmo, aes(x=birth_mo))+geom_bar( )+
#   scale_x_discrete(limits = month.name)+
#   theme(axis.text.x = element_text(face = "italic", size = 8, angle = 90, vjust = 0.5))+
#   facet_grid(~inc_birthmo$income)
# 
# ggplot(inc_loc, aes(x=location))+geom_bar( )+
#   scale_x_continuous(breaks=unique(inc_loc$location))+
#   theme(axis.text.x = element_text(face = "italic", size = 8, angle = 90, vjust = 0.5))+
#   facet_grid(~inc_loc$income)
  
## Regression
inc_analysis$inc_bin[inc_analysis$income == 'Less than 18,000' | inc_analysis$income == '18,000 to 30,999'] <- 1 # "Low"
inc_analysis$inc_bin[inc_analysis$income == '52,000 to 100,000' | inc_analysis$income == 'Greater than 100,000'] <- 0 # "High"
inc_analysis$inc_bin[inc_analysis$income == '31,000 to 51,999'] <- NA
table(inc_analysis$inc_bin, exclude=NULL)
inc_analysis$inc_bin <- as.numeric(inc_analysis$inc_bin)

inc_analysis$race_fct <- as.factor(inc_analysis$race)

inc_analysis$employment <- as.character(inc_analysis$empl_status)
inc_analysis$employment <- as.factor(inc_analysis$employment)

inc_analysis$sexfct <- as.character(inc_analysis$sex)
inc_analysis$sexfct <- as.factor(inc_analysis$sexfct)

a <- glm(inc_bin ~ sexfct + birth_yr + birth_mo + location + employment + BMI + diabetes + asthma + mi_flag + cancer_dx_dt + relevel(race_fct, ref = "British"), data=inc_analysis, family = binomial())
summary(a)
exp(cbind(coef(a), confint(a)))

b <- glm(inc_bin ~ sexfct + birth_yr + location + employment + BMI + diabetes + asthma + mi_flag + cancer_dx_dt + relevel(race_fct, ref = "British"), data=inc_analysis, family = binomial())
summary(b)
exp(cbind(coef(b), confint(b)))

summary(glm(inc_bin ~ sexfct, data=inc_analysis, family = binomial()))
summary(glm(inc_bin ~ relevel(race_fct, ref="British"), data=inc_analysis, family = binomial()))
summary(glm(inc_bin ~ birth_yr, data=inc_analysis, family = binomial()))
summary(glm(inc_bin ~ birth_mo, data=inc_analysis, family = binomial()))
summary(glm(inc_bin ~ location, data=inc_analysis, family = binomial()))
summary(glm(inc_bin ~ employment, data=inc_analysis, family = binomial()))
summary(glm(inc_bin ~ BMI, data=inc_analysis, family = binomial()))
summary(glm(inc_bin ~ diabetes, data=inc_analysis, family = binomial()))
summary(glm(inc_bin ~ asthma, data=inc_analysis, family = binomial()))
summary(glm(inc_bin ~ mi_flag, data=inc_analysis, family = binomial()))
summary(glm(inc_bin ~ cancer_dx_dt, data=inc_analysis, family = binomial()))

summary(glm(sexfct ~ employment, data=inc_analysis, family = binomial()))
summary(glm(diabetes ~ employment, data=inc_analysis, family = binomial()))
summary(glm(asthma ~ employment, data=inc_analysis, family = binomial()))
summary(glm(mi_flag ~ employment, data=inc_analysis, family = binomial()))
summary(glm(cancer_dx_dt ~ employment, data=inc_analysis, family = binomial()))

# Model Selection (Stepwise)
# selectmod <- step(bmi.mod)
# summary(selectmod)

#Chi-Sq & ANOVA
fit <- aov(BMI ~ income, data=inc_analysis)
summary(fit)

# write.table(inc_analysis, "inc_analysis.txt", sep="\t", row.names=F, col.names=T, quote=FALSE)
# 
# yr_inc_dist <- as.data.frame(t(table(inc_analysis$income, inc_analysis$birth_yr)))
# foo <- spread(yr_inc_dist, Var2, Freq)