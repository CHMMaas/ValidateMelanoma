# erase memory
rm(list = ls(all.names = TRUE))

# for reproducibility
set.seed(100)

# package
library(dplyr)
library(naniar)
# remotes::install_github("CHMMaas/PredictionTools")
library(PredictionTools)
library(survival)
library(survminer)
library(dcurves)
library(grid)
library(ggpubr)

# load data file
# in Excel: changed DOB 1974 instead of 2974, and 1960 instead of 2960.
orig.val.dat <- read.csv(file="Z:/Project Melanoom/PaperValidation/Data/Roger_ext_validation - July 2024.csv",
                         stringsAsFactors = TRUE)

###
### DATA INSPECTION
###
Rec.df <- orig.val.dat %>% filter(Recurrence=="Yes")
tab.Rec <- table(Rec.df$Status, Rec.df$Cause_of_death)
no.Rec.df <- orig.val.dat %>% filter(Recurrence=="No")
miss.df <- orig.val.dat %>% filter(Recurrence=="")
outcome.df <- data.frame(Recurrence=rep(c("Yes", "No", ""), each=32),
                         Cause_of_death=rep(row.names(table(Rec.df$Cause_of_death, Rec.df$Status)), each=8),
                         Status=rep(row.names(tab.Rec), 12),
                         N=c(tab.Rec,
                             table(no.Rec.df$Status, no.Rec.df$Cause_of_death),
                             table(miss.df$Status, miss.df$Cause_of_death)))
outcome.df <- outcome.df %>% mutate(CO=case_when(Recurrence=="Yes" |
                                                   Status=="Alive with disease" |
                                                   Status=="Deceased with disease" |
                                                   Status=="Deceased without disease" |
                                                   Status=="Deceased with unknown disease status" |
                                                   Cause_of_death=="Melanoma" |
                                                   Cause_of_death=="Other case (not melanoma)" ~ 1,
                                                 TRUE ~ 0))

###
### CREATE VARIABLES
###
# diagnose year
orig.val.dat$diag.year <- as.numeric(format(as.POSIXct(orig.val.dat$diagnose_date, form="%Y-%m-%d"), "%Y"))

# Age at SN procedure
orig.val.dat$Age <- as.numeric(floor((as.POSIXct(orig.val.dat$SNdate, form="%Y-%m-%d") -
                                        as.POSIXct(orig.val.dat$DOB, form="%Y-%m-%d"))/(24*60*60)/365.25))

# Rec.composite, dead, MSM, and merge Histology categories that are small
orig.val.dat <- orig.val.dat %>%
  mutate(Rec.composite = case_when(Recurrence=="Yes" |
                                     Status=="Alive with disease" |
                                     Status=="Deceased with disease" |
                                     Status=="Deceased without disease" |
                                     Status=="Deceased with unknown disease status" |
                                     Cause_of_death=="Melanoma" |
                                     Cause_of_death=="Other case (not melanoma)" ~ 1,
                                   TRUE ~ 0),
         dead = case_when(Status=="Deceased with disease" |
                            Status=="Deceased without disease" |
                            Status=="Deceased with unknown disease status" |
                            Cause_of_death=="Melanoma" |
                            Cause_of_death=="Other case (not melanoma)" ~ 1,
                          TRUE ~ 0),
         MSM = case_when(Cause_of_death=="Melanoma" ~ 1, 
                         TRUE ~ 0),
         Histology = case_when(
           Histology=="Melanoma of childhood" |
             Histology=="Persistent melanoma" |
             Histology=="Melanoma arising from blue" |
             Histology=="Melanoma arising in giant" |
             Histology=="Neurotropic" ~ "Other (please specify)",
           TRUE ~ Histology))

# difference diagnose date and SN date
orig.val.dat$diag_SN_wk <- as.numeric((as.POSIXct(orig.val.dat$SNdate, form="%Y-%m-%d") -
                                         as.POSIXct(orig.val.dat$diagnose_date, form="%Y-%m-%d"))/(24*60*60)/52)

# logRec.composite.time
orig.val.dat$Recurrence.time <- as.numeric((as.POSIXct(orig.val.dat$Date_1st_Rec, form="%Y-%m-%d") -
                                              as.POSIXct(orig.val.dat$SNdate, form="%Y-%m-%d"))/(24*60*60)/365.25)
orig.val.dat$FU.time <- as.numeric((as.POSIXct(orig.val.dat$Last_FU, form="%Y-%m-%d") -
                                      as.POSIXct(orig.val.dat$SNdate, form="%Y-%m-%d"))/(24*60*60)/365.25)
orig.val.dat$logRecurrence.time <- log(orig.val.dat$Recurrence.time)
orig.val.dat$logFU.time <- log(orig.val.dat$FU.time)

# add time until composite outcome and country
orig.val.dat <- orig.val.dat %>% mutate(
  logRec.composite.time = case_when(Recurrence=="Yes" ~ logRecurrence.time,
                                    TRUE ~ logFU.time),
  Country = case_when(Site=="Brown U/RI Hospital/LifeSpan Oncology" ~ "USA",
                      Site=="California Pacific Med Ctr/UCSF" ~ "USA",
                      Site=="Carolinas Medical Center" ~ "USA",
                      Site=="City of Hope" ~ "USA",
                      Site=="Istituto Nazionale Tumori" ~ "Italy", 
                      Site=="Loma Linda University" ~ "USA",
                      Site=="Moffitt Cancer Center" ~ "USA",
                      Site=="Oregon Health & Science University" ~ "USA",
                      Site=="Rush University" ~ "USA",
                      Site=="Sahlgrenska University Hospital" ~ "Sweden",
                      Site=="Southern Arizona VA/University of Arizona" ~ "USA",
                      Site=="St. Joseph Hospital" ~ "USA",
                      Site=="Tel Aviv Sourasky Medical Center" ~ "Israel",
                      Site=="University Medical Center Groningen" ~ "Netherlands",
                      Site=="University of Miami/Sylvester CCC (Jackson Memorial Hospital)" ~ "USA",
                      Site=="University of South Alabama" ~ "USA",
                      Site=="Ventura County Medical Center" ~ "USA",
                      startsWith(as.character(DatabaseNR), "MCA") ~ "USA", # those who have unknown country have ID starting with MCA which stands for Mayo Clinic Arizona
                      TRUE ~ "Unknown"))

###
### EXCLUSION
###
melanoma <- orig.val.dat %>%
  filter(Other_cancer_yn!="Yes")
incl.disease <- melanoma %>%
  filter(melanoma$diag.year>=1991 & # immunotherapy is excluded, so we can include these people, KM doesn't change much too
           melanoma$diag.year<=2018 &
           Microsatellitosis!="Yes" &
           Intra_lymphatic_metastasis!="Yes" &
           n_staging!="N1b" &
           n_staging!="N1c" &
           n_staging!="N2b" &
           n_staging!="N2c" &
           n_staging!="N3b" &
           n_staging!="N3c" &
           m_staging=="M0") 
incl.treat <- incl.disease %>%
  filter(neoadj!="Yes" &
           Immuno_therapy!="Yes" &
           Immuno_therapy_date=="" &
           Radiation_therapy!="Yes" &
           Radiation_therapy_date=="" &
           Chemo_therapy!="Yes" &
           Chemo_therapy_date=="" &
           Target_therapy!="Yes" &
           Target_therapy_date=="")
sel.val.dat <- incl.treat %>%
  filter(SNdate!="" &
           Last_FU!="" &
           !(incl.treat$Recurrence=="" & incl.treat$Cause_of_death=="" & incl.treat$Status=="") &
           !(incl.treat$Recurrence=="" & incl.treat$Cause_of_death=="" & incl.treat$Status=="Lost to follow-up") &
           !(incl.treat$Recurrence=="" & incl.treat$Cause_of_death=="Unknown" & incl.treat$Status=="") &
           !(incl.treat$Recurrence=="" & incl.treat$Cause_of_death=="Unknown" & incl.treat$Status=="Lost to follow-up") &
           FU.time > 0 &
           !(Recurrence=="Yes" & Date_1st_Rec=="") &
           !(Recurrence=="Yes" & Recurrence.time <= 0) &
           diag_SN_wk>-14 &
           DOB!="" &
           Age > 0)
cat(" Number of patients in original data set:", nrow(orig.val.dat), "\n",
    "   Exclude number of patients with other cancer:",
    sum(orig.val.dat$Other_cancer_yn=="Yes"), "\n",
    "Number of patients after excluding other cancers:", nrow(melanoma), "\n",
    "   Exclude number of patients diagnosed before 1991:",
    sum(melanoma$diag.year<1991, na.rm=TRUE), "\n",
    "   Exclude number of patients diagnosed after 2018:",
    sum(melanoma$diag.year>2018, na.rm=TRUE), "\n",
    "   Exclude number of patients with microsatellitosis, intra lymphatic metastasis, N-stage N1b, 2b, 2c, 3b, 3c:",
    sum(melanoma$Microsatellitosis=="Yes")+ 
    sum(melanoma$Intra_lymphatic_metastasis=="Yes")+ 
    sum(melanoma$n_staging=="N1b")+ 
    sum(melanoma$n_staging=="N2b")+ 
    sum(melanoma$n_staging=="N2c")+ 
    sum(melanoma$n_staging=="N3b")+ 
    sum(melanoma$n_staging=="N3c"), "\n",
    "   Exclude number of patients with M-stage M:",
    sum(melanoma$m_staging=="M1")+
      sum(melanoma$m_staging=="M1a")+
      sum(melanoma$m_staging=="M1b"), "\n",
    "   Exclude number of patients with unknown M-stage:",
    sum(melanoma$m_staging=="Mx" | melanoma$m_staging=="6" | melanoma$m_staging=="7" | melanoma$m_staging==""), "\n",
    "Number of patients after excluding other diseases:", nrow(incl.disease), "\n",
    "   Exclude number of patients with neo adjuvant therapy or immunotherapy:",
    sum(incl.disease$Immuno_therapy=="Yes" | incl.disease$Immuno_therapy_date!="") +
      sum(incl.disease$neoadj=="Yes"), "\n",
    "   Exclude number of patients with radiation therapy:",
    sum(incl.disease$Radiation_therapy=="Yes" | incl.disease$Radiation_therapy_date!=""), "\n",
    "   Exclude number of patients with chemotherapy:",
    sum(incl.disease$Chemo_therapy=="Yes" | incl.disease$Chemo_therapy_date!=""), "\n",
    "   Exclude number of patients with targeted therapy:",
    sum(incl.disease$Target_therapy=="Yes" | incl.disease$Target_therapy_date!=""), "\n",
    "Number of patients after excluding other treatments:", nrow(incl.treat), "\n",
    "   Exclude number of patients with unknown SN date and unknown last follow-up date:",
    sum(incl.treat$SNdate=="")+
      sum(incl.treat$Last_FU==""), "\n",
    "   Exclude number of patients with unknown Recurrence, COD and Status:",
    sum((incl.treat$Recurrence=="" & incl.treat$Cause_of_death=="" & incl.treat$Status=="") |
          (incl.treat$Recurrence=="" & incl.treat$Cause_of_death=="" & incl.treat$Status=="Lost to follow-up") |
          (incl.treat$Recurrence=="" & incl.treat$Cause_of_death=="Unknown" & incl.treat$Status=="") |
          (incl.treat$Recurrence=="" & incl.treat$Cause_of_death=="Unknown" & incl.treat$Status=="Lost to follow-up")), "\n",
    "   Exclude number of patients with date of follow-up before or on the same day as SN date, \n unknown date of recurrence, date of recurrence before or on the same day as SN date:",
    sum(incl.treat$FU.time<=0, na.rm=TRUE)+
      sum(incl.treat$Recurrence=="Yes" & incl.treat$Date_1st_Rec=="", na.rm=TRUE)+
      sum(incl.treat$Recurrence=="Yes" & incl.treat$Recurrence.time<=0, na.rm=TRUE), "\n",
    "   Exclude number of patients with diagnose date 14 weeks before SN date:",
    sum(incl.treat$diag_SN_wk<-14), "\n",
    "   Exclude number of patients with SN date before DOB or at same date:",
    sum(incl.treat$Age<=0, na.rm=TRUE), "\n",
    "   Exclude number of patients with unknown DOB:",
    sum(incl.treat$DOB==""), "\n",
    "Number of patients after excluding registration issues:", nrow(sel.val.dat), "\n")

# save sites sample size
openxlsx::write.xlsx(data.frame(table(sel.val.dat$Site)),
                     file="Z:/Project Melanoom/PaperValidation/Results/Sites.xlsx")

###
### DATA VISUALIZATION
###
# define outcomes
S.Rec <- survival::Surv(exp(sel.val.dat$logRec.composite.time),
                        sel.val.dat$Rec.composite)
S.MSM <- survival::Surv(exp(sel.val.dat$logFU.time),
                        sel.val.dat$MSM)
all.med.FU <- stats::quantile(prodlim::prodlim(prodlim::Hist(exp(logFU.time),
                                                             dead)~1,
                                               data=sel.val.dat,
                                               reverse=TRUE))$quantile

# median follow-up time
med.FU <- stats::quantile(prodlim::prodlim(prodlim::Hist(exp(logFU.time),
                                                         dead)~Country,
                                           data=sel.val.dat,
                                           reverse=TRUE))
df.med.FU <- data.frame(Country=med.FU$Country,
                        q=med.FU$q,
                        quantile=med.FU$quantile) |> filter(q==0.25 | q==0.5 | q==0.75)
result.med.FU <- c()
for (country.name in unique(df.med.FU$Country)){
  med.FU.country <- as.numeric(df.med.FU %>% filter(Country==country.name & q == 0.5) %>% select(quantile))
  lower.FU.country <- as.numeric(df.med.FU %>% filter(Country==country.name & q == 0.75) %>% select(quantile))
  upper.FU.country <- as.numeric(df.med.FU %>% filter(Country==country.name & q == 0.25) %>% select(quantile))
  result.med.FU <- c(result.med.FU, paste0(sprintf("%.1f", med.FU.country), ", IQR:",
                                           sprintf("%.1f", lower.FU.country), "-",
                                           sprintf("%.1f", upper.FU.country)))
}
openxlsx::write.xlsx(data.frame(country=c("All", unique(df.med.FU$Country)),
                                FU=c(paste0(sprintf("%.1f", all.med.FU[3]), ", IQR:",
                                            sprintf("%.1f", all.med.FU[4]), "-",
                                            sprintf("%.1f", all.med.FU[2])),
                                     result.med.FU)),
                     file="Z:/Project Melanoom/PaperValidation/Results/med.FU.xlsx")

# limit to 5 years
horizon <- 5
S.Rec.5 <- S.Rec
S.Rec.5[S.Rec[,1]>horizon, 1] <- horizon
S.Rec.5[S.Rec[,1]>horizon, 2] <- 0
S.MSM.5 <- S.MSM
S.MSM.5[S.MSM[,1]>horizon, 1] <- horizon
S.MSM.5[S.MSM[,1]>horizon, 2] <- 0

# make Table with number of events and KM estimates
est.CI.n <- function(KM=NULL, horizon=5){
  t <- max(KM$time[KM$time<=horizon])
  est <- 100*KM$surv[KM$time==t]
  est.LCI <- 100*KM$lower[KM$time==t]
  est.UCI <- 100*KM$upper[KM$time==t]
  n <- sum(KM$n.event[KM$time<=horizon])
  return(list(t=t, est=est, est.LCI=est.LCI, est.UCI=est.UCI, n=n))
}

# sample sizes
KM.est <- c()
country.names <- c("USA", "Sweden", "Israel")
for (population in c("Total", country.names)){
  KM.est <- rbind(KM.est, c(population, rep("", 4)))
  for (patients in c("all", "pos", "neg")){
    if (patients=="pos" & population=="Total"){
      sel.patients <- sel.val.dat$SNstatus=="Positive"
    } else if (patients=="neg" & population=="Total"){
      sel.patients <- sel.val.dat$SNstatus=="Negative"
    } else if (patients=="all" & population=="Total"){
      sel.patients <- rep(TRUE, nrow(sel.val.dat))
    } else if (patients=="pos" & population!="Total"){
      sel.patients <- sel.val.dat$SNstatus=="Positive" & sel.val.dat$Country==population
    } else if (patients=="neg" & population!="Total"){
      sel.patients <- sel.val.dat$SNstatus=="Negative" & sel.val.dat$Country==population
    } else if (patients=="all" & population!="Total"){
      sel.patients <- sel.val.dat$Country==population
    }
    
    n <- sum(sel.patients)
    
    row <- c(paste0(ifelse(patients=="all", "All",
                           ifelse(patients=="pos", "Positive SN", "Negative SN")),
                    " (N=", n, ")"))
    for (outcome in c("Rec", "MSM")){
      S <- eval(parse(text=paste0("S.", outcome, ".5")))
      out <- est.CI.n(survival::survfit(S~1,
                                        subset=sel.patients,
                                        data=sel.val.dat),
                      horizon=horizon)
      row <- c(row, paste0(out$n, " (",
                           round(out$n/n*100, 0), "%)"),
               paste0(sprintf("%.1f", out$est), "% (",
                      sprintf("%.1f", out$est.LCI), "-",
                      sprintf("%.1f", out$est.UCI), ")"))
      
    }
    KM.est <- rbind(KM.est, row)
  }
}
openxlsx::write.xlsx(data.frame(KM.est),
                     file="Z:/Project Melanoom/PaperValidation/Results/KM.est.xlsx")

# save KM plots
surv.Rec <- survival::survfit(S.Rec.5~Country,
                              subset=Country!="Netherlands",
                              data=sel.val.dat)
ggplot2::ggsave(file="Z:/Project Melanoom/PaperValidation/Results/RFS.KM.png",
                plot=ggpubr::ggarrange(survminer::ggsurvplot(surv.Rec,
                                                             ylab="Recurrence-free survival probability",
                                                             conf.int=TRUE,
                                                             legend.title="",
                                                             risk.table=FALSE)$plot,
                                       survminer::ggsurvplot(surv.Rec,
                                                             risk.table=TRUE)$table,
                                       nrow=2, heights=c(4, 1), align="v"),
                width=8, height=8, dpi=300)
surv.MSM <- survival::survfit(S.MSM.5~Country,
                              subset=Country!="Netherlands",
                              data=sel.val.dat)
ggplot2::ggsave(file="Z:/Project Melanoom/PaperValidation/Results/MSS.KM.png",
                plot=ggpubr::ggarrange(survminer::ggsurvplot(surv.MSM,
                                                             ylab="Melanoma-specific survival probability",
                                                             conf.int=TRUE,
                                                             legend.title="",
                                                             risk.table=FALSE)$plot,
                                       survminer::ggsurvplot(surv.MSM,
                                                             risk.table=TRUE)$table,
                                       nrow=2, heights=c(4, 1), align="v"),
                width=8, height=8, dpi=300)

###
### VARIABLE CODING
###
# select variables
to.be.imp.data <- sel.val.dat %>% select(SNstatus, Sex, Ulceration, Loc_CAT, Age,
                                         Histology, Breslow, Rdamcrit, Dewar, Mitosis,
                                         Rec.composite, logRec.composite.time,
                                         MSM, logFU.time, Country)
# not in data: multiple.fields, Tot_SNs_neg, Tot_SNs_pos

# code missingness
to.be.imp.data <- to.be.imp.data %>%
  mutate(SNstatus = case_when(SNstatus=="" | SNstatus=="Unknown" ~ NA,
                              TRUE ~ SNstatus),
         Sex = case_when(Sex=="" | Sex=="Unknown" ~ NA,
                         TRUE ~ Sex),
         Loc_CAT = case_when(Loc_CAT=="Upper extremity" ~ "arm",
                             Loc_CAT=="Lower extremity" ~ "leg",
                             Loc_CAT=="Trunk" ~ "trunk",
                             Loc_CAT=="Head/neck" ~ "head & neck",
                             TRUE ~ NA),
         Ulceration = case_when(Ulceration=="" | Ulceration=="Unknown" ~ NA,
                                TRUE ~ Ulceration),
         Histology = case_when(Histology=="" | Histology=="Unknown" ~ NA,
                               TRUE ~ Histology),
         Breslow = case_when(Breslow==999 | Breslow==997 | Breslow==301 | Breslow==99 | Breslow==0 ~ NA,
                             TRUE ~ Breslow),
         Rdamcrit = case_when(Rdamcrit=="<0.1mm" ~ 0.05,
                              Rdamcrit==">=0.1 mm to <=1.0 mm" ~ 0.55,
                              Rdamcrit==">1.0 mm" ~ 3.93,
                              Rdamcrit=="Unknown" & SNstatus=="Negative" ~ 1,
                              Rdamcrit=="" & SNstatus=="Negative" ~ 1,
                              Rdamcrit=="Unknown" & SNstatus=="Positive" ~ NA, # impute
                              Rdamcrit=="" & SNstatus=="Positive" ~ NA,        # impute
                              Rdamcrit=="Unknown" & is.na(SNstatus) ~ NA,      # impute
                              Rdamcrit=="" & is.na(SNstatus) ~ NA,             # impute
                              TRUE ~ 1),
         Dewar = case_when((Dewar=="" | Dewar=="Unknown") & SNstatus=="Positive" ~ NA,
                           (Dewar=="" | Dewar=="Unknown") & is.na(SNstatus) ~ NA,    
                           (Dewar=="" | Dewar=="Unknown") & SNstatus=="Negative" ~ NA,
                           TRUE ~ Dewar),
         Mitosis = case_when(Mitosis=="" | Mitosis=="Unknown" ~ NA,
                             TRUE ~ Mitosis))

# drop levels
to.be.imp.data <- droplevels(to.be.imp.data)

# missingness
colSums(is.na(to.be.imp.data))

# make factors
to.be.imp.data$Country <- as.factor(to.be.imp.data$Country)
to.be.imp.data$Loc_CAT <- as.factor(to.be.imp.data$Loc_CAT)
to.be.imp.data$Histology <- as.factor(to.be.imp.data$Histology)
to.be.imp.data$Dewar <- as.factor(to.be.imp.data$Dewar)

###
### DESCRIPTIVES
###
save(to.be.imp.data, file="Z:/Project Melanoom/PaperValidation/Data/to.be.imp.val.data.Rdata")

# remove Dewar due to very small groups
to.be.imp.data <- to.be.imp.data |> select(-Dewar)

###
### SINGLE IMPUTATION
###
m <- 5
new.imputation <- FALSE
if (new.imputation){
  mice.data <- mice::mice(to.be.imp.data, m=m, maxit=5, print=FALSE, seed=500)
  print(mice.data$loggedEvents)
  imp.data <- mice::complete(mice.data, m)
  save(mice.data, imp.data, file="Z:/Project Melanoom/PaperValidation/Data/imp.val.data.Rdata")
} else{
  load("Z:/Project Melanoom/PaperValidation/Data/imp.val.data.Rdata")
}

###
### MAKE PREDICTIONS
###
predict.RFS.MSS <- function(input){
  p.RFS <- c()
  p.MSS <- c()
  for (patient in 1:nrow(input)){
    # get input from patient i
    input.patient <- input[patient,]
    input.numeric <- c(as.numeric(input.patient$SNstatus=="Positive"),
                       as.numeric(input.patient$Age),
                       as.numeric(input.patient$Ulceration=="Yes"),
                       as.numeric(input.patient$Loc_CAT=="leg"),
                       as.numeric(input.patient$Loc_CAT=="trunk"),
                       as.numeric(input.patient$Loc_CAT=="head & neck"),
                       as.numeric(log(input.patient$Breslow) - 0.71779603),
                       as.numeric(input.patient$SNstatus=="Positive")*as.numeric(log(input.patient$Rdamcrit)),
                       as.numeric(input.patient$SNstatus=="Positive")*as.numeric(log(input.patient$Breslow) - 0.71779603))
    
    # calculate linear predictor for patient i
    lp <- -1.282591+
      c(1.08621071, 0.01366932, 0.47862830, 0.12834775, 0.25329321, 0.63572751,
        0.83026128, 0.21189902, -0.43670541) %*% input.numeric
    
    # calculate 5-year RFS and 5-YEAR mss
    p.RFS <- c(p.RFS, exp(-0.22846279*exp(lp)))
    p.MSS <- c(p.MSS, exp(-0.090758801*exp(1.1384077*lp)))
  }
  
  return(list(p.RFS=p.RFS, p.MSS=p.MSS))
}

# change name to p.Rec and p.MSM
pred.Rec <- c()
pred.MSM <- c()
for (i in 1:m){
  pred <- predict.RFS.MSS(mice::complete(mice.data, i))
  pred.Rec <- cbind(pred.Rec, 1-pred$p.RFS)
  pred.MSM <- cbind(pred.MSM, 1-pred$p.MSS)
}

# calculate AJCC8 
cindex <- c()
cindex.se <- c()
for (i in 1:m){
  AJCC.i <- data.frame(Breslow=mice::complete(mice.data, i)$Breslow,
                       Ulceration=mice::complete(mice.data, i)$Ulceration) |>
    mutate(AJCC8 = case_when(
      Breslow<=0.8 & Ulceration=="No" ~ 1.1,
      Breslow<=0.8 & Ulceration=="Yes" ~ 1.2,
      Breslow>0.8 & Breslow<=1 ~ 1.2,
      Breslow>1 & Breslow<=2 & Ulceration=="No" ~ 2.1,
      Breslow>1 & Breslow<=2 & Ulceration=="Yes" ~ 2.2,
      Breslow>2 & Breslow<=4 & Ulceration=="No" ~ 3.1,
      Breslow>2 & Breslow<=4 & Ulceration=="Yes" ~ 3.2,
      Breslow>4 & Ulceration=="No" ~ 4.1,
      Breslow>4 & Ulceration=="Yes" ~ 4.2,
    ))
  
  # C-index of AJCC8
  rc <- Hmisc::rcorr.cens(-AJCC.i$AJCC8, S.Rec.5)
  cindex <- c(cindex, rc["C Index"])
  cindex.se <- c(cindex.se, rc["S.D."]/2)
}
cindex.mi<-Rubin.combine(cindex, cindex.se)
cat("C-index AJCC8:", round(cindex.mi$est, 2),
    "95% CI:", round(cindex.mi$est-1.96*cindex.mi$se, 2),
    "-", round(cindex.mi$est+1.96*cindex.mi$se, 2), "\n")

# save calibration figures
for (site in c("Total", country.names)){
  for (status in c("All", "Positive", "Negative")[1]){
    if (site=="Total" & status=="All"){
      sel.patients <- rep(TRUE, nrow(to.be.imp.data))
    } else if (site=="Total" & status!="All"){
      sel.patients <- imp.data$SNstatus==status
    } else if (site!="Total" & status=="All"){
      sel.patients <- to.be.imp.data$Country==site
    } else if (site!="Total" & status!="All"){
      sel.patients <- to.be.imp.data$Country==site & imp.data$SNstatus==status
    }
    cat("Make calibration plot for", site, "patients, i.e.,", sum(sel.patients), "\n")
    png(file=paste0("Z:/Project Melanoom/PaperValidation/Results/Calibration plots/",
                    "cal.RFS.", site, ".", status, ".png"),
        width=16, height=16, units="cm", res=300)
    PredictionTools::val.surv.mi(p=pred.Rec[sel.patients, ],
                                 y=S.Rec.5[sel.patients, ], 
                                 time=horizon,
                                 show.metrics=c(rep(TRUE, 5), TRUE),
                                 main=ifelse(site=="Total" & status=="All", "All patients",
                                             ifelse(site=="Total" & status=="Positive", "Patients with positive SN",
                                                    ifelse(site=="Total" & status=="Negative", "Patients with negative SN", site))),
                                 CI.metrics=TRUE, n.sim=100)
    dev.off()
    png(file=paste0("Z:/Project Melanoom/PaperValidation/Results/Calibration plots/cal.MSS.",
                    site, ".", status, ".png"),
        width=16, height=16, units="cm", res=300)
    PredictionTools::val.surv.mi(p=pred.MSM[sel.patients, ],
                                 y=S.MSM.5[sel.patients, ], 
                                 time=horizon,
                                 show.metrics=c(rep(TRUE, 5), TRUE), 
                                 main=ifelse(site=="Total" & status=="All", "All patients",
                                             ifelse(site=="Total" & status=="Positive", "Patients with positive SN",
                                                    ifelse(site=="Total" & status=="Negative", "Patients with negative SN", site))),
                                 CI.metrics=TRUE, n.sim=100)
    dev.off()
    
    if (status=="All"){
      for (outcome in c("RFS", "MSS")){
        cat("Make decision curve for", outcome, "in", site, "patients, i.e.,",
            sum(sel.patients), "\n")
        S <- eval(parse(text=paste0("S.", ifelse(outcome=="RFS", "Rec.5", "MSM.5"))))
        pred <- eval(parse(text=paste0("pred.", ifelse(outcome=="RFS", "Rec", "MSM"))))
        
        # dca plot
        dca.plot <- dcurves::dca(y ~ Model,
                                 time=horizon,
                                 data=data.frame(y=S[sel.patients, ],
                                                 Model=rowMeans(pred[sel.patients, ]))) %>%
          plot(smooth = TRUE) 
        
        # histogram
        hist.pred <- ggplot(data=data.frame(Predictions=rowMeans(pred.Rec)),
                            aes(x=Predictions)) +
          geom_histogram(fill="#619bff", color="grey", alpha=0.6, binwidth=0.01) +
          ylab("Count") +
          theme_minimal() +
          theme(panel.background = element_rect(fill="white", color=NA),
                plot.background = element_rect(fill="white", color=NA))
        
        # save
        ggsave(file=paste0("Z:/Project Melanoom/PaperValidation/Results/dca/dca.",
                           outcome, ".", site, ".png"),
               plot=ggarrange(dca.plot, hist.pred,
                              nrow=2, ncol=1, heights=c(2, 1), align="v"),
               width=16, height=16, units="cm", dpi=300)
        
        # save net benefit at certain treatment thresholds
        openxlsx::write.xlsx(data.frame(dca.plot$data |>
                                          filter(threshold==0.44 |
                                                   threshold==0.45 |
                                                   threshold==0.5 |
                                                   threshold==0.55 |
                                                   threshold==0.6) |>
                                          arrange(threshold) |>
                                          select(label, threshold, net_benefit) |>
                                          mutate(threshold=threshold*100) |>
                                          mutate_if(is.numeric, round, 2)),
                             file=paste0("Z:/Project Melanoom/PaperValidation/Results/dca/dca.",
                                         outcome, ".", site, ".xlsx"))
      }
    }
  }
}

###
### Manual calculation of DCA
###
dca.out <- dcurves::dca(y ~ Model,
                        time=horizon,
                        data=data.frame(y=S.Rec,
                                        Model=rowMeans(pred.Rec))) %>%
  plot(smooth = TRUE)
NB.table <- c()
for (t in c(0.44, 0.45, 0.5, 0.55, 0.6)){
  # deduct TP, FP, TN, FN from output
  POS_rate <- as.numeric(dca.out$data |> filter(threshold==t & variable=="Model") |> select(pos_rate))
  TP_rate <- as.numeric(dca.out$data |> filter(threshold==t & variable=="Model") |> select(tp_rate))
  FP_rate <- as.numeric(dca.out$data |> filter(threshold==t & variable=="Model") |> select(fp_rate))
  N <- nrow(S.Rec.5)
  TP <- TP_rate*N
  FP <- FP_rate*N
  TN <- (POS_rate-TP_rate)*N
  FN <- N*(1-(FP_rate+POS_rate))
  
  # calculate metrics
  sensitivity <- sprintf("%.1f", TP/(TP+FN)*100)
  specificity <- sprintf("%.1f", TN/(TN+FP)*100)
  PPV <- sprintf("%.1f", TP/(TP+FP)*100)
  NPV <- sprintf("%.1f", TN/(TN+FN)*100)
  NB_model <- sprintf("%.2f", as.numeric(dca.out$data |> filter(threshold==t & variable=="Model") |> select(net_benefit)))
  NB_all <- sprintf("%.2f", as.numeric(dca.out$data |> filter(threshold==t & variable=="all") |> select(net_benefit)))
  NB_none <- sprintf("%.2f", as.numeric(dca.out$data |> filter(threshold==t & variable=="none") |> select(net_benefit)))
  
  # save
  NB.table <- rbind(NB.table, 
                    c(t*100, sensitivity, specificity, 
                      sprintf("%.0f", TP), sprintf("%.0f", FP), 
                      sprintf("%.0f", TN), sprintf("%.0f", FN), 
                      PPV, NPV, 
                      NB_model, NB_all, NB_none))
}
NB.table <- data.frame(NB.table)
colnames(NB.table) <- c("Probability threshold (%)", 
                        "Sensitivity (%)", "Specificity (%)", "TP", "FP", "TN", "FN", 
                        "PPV (%)", "NPV (%)", "Net benefit of the model", 
                        "Net benefit treat all", "Net benefit treat none")
openxlsx::write.xlsx(NB.table,
                     file="Z:/Project Melanoom/PaperValidation/Results/NB.table.xlsx")

###
### Combine calibration plots for countries
###
for (status in c("All", "Positive", "Negative")[1]){
  for (outcome in c("RFS", "MSS")){
    png(file=paste0("Z:/Project Melanoom/PaperValidation/Results/",
                    ifelse(status!="All", "Calibration plots/", ""),
                    "cal.combined.",
                    outcome, ".", status, ".png"),
        width=16, height=16, units="cm", res=300)
    par(mar=rep(0, 4))
    layout(matrix(1:4, ncol=2, byrow=TRUE))
    for (site in c("Total", country.names)){
      plot(NA, xlim=0:1, ylim=0:1, xaxt="n", yaxt="n", bty="n")
      img <- png::readPNG(paste0("Z:/Project Melanoom/PaperValidation/Results/Calibration plots/cal.",
                                 outcome, ".", site, ".", status, ".png"))
      rasterImage(img, 0, 0, 1, 1)
    }
    grDevices::dev.off()
  }
}

###
### Combine all, positive, and negative SN calibration plots
###
for (outcome in c("RFS", "MSS")){
  for (patients in c("All", "Positive", "Negative")[1]){
    img <- png::readPNG(paste0("Z:/Project Melanoom/PaperValidation/Results/Calibration plots/cal.",
                               outcome, ".Total.", patients, ".png"))
    img_grob <- grid::rasterGrob(img, interpolate=TRUE)
    assign(paste0("img.", patients), img_grob)
  }
  ggsave(file=paste0("Z:/Project Melanoom/PaperValidation/Results/cal.",
                     outcome, ".png"),
         plot=ggpubr::ggarrange(img.All,
                                ggpubr::ggarrange(img.Positive, img.Negative,
                                                  ncol=2),
                                ncol=1, nrow=2,
                                heights=c(1, 1), widths=c(1, 1)),
         width=10, height=10, dpi=300)
}