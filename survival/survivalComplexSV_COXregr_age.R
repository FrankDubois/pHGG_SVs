library("survival")
library("survminer")
library(data.table)
library("splines")
library("lattice")

metaDsurvSig <- fread('/Users/fdubois/Dropbox/pHGG_Data/data/20210107survival/20210107metaD_survSig_TP53.txt')
metaDsurvSig$Overall_Survival_Months <- as.numeric(metaDsurvSig$Overall_Survival_Months)
metaDsurvSig[,status := ifelse(test = Surviva_Status ==  'DECEASED', yes = 2, no = 1)]
metaDsurvSig$ComplexSVsign <- as.numeric(metaDsurvSig$ComplexSVsign)
metaD_surv <- metaDsurvSig[!is.na(Overall_Survival_Months)]

table(metaDsurvSig$Location)
metaD_survMid <- metaD_surv[Location == 'Midline']
table(metaD_survMid$Histone_group)
metaD_survMid <- metaD_survMid[!is.na(ComplexSVsign)]
fit <- coxph(Surv(Overall_Survival_Months, status) ~  Histone_group + Age + TP53snv  + ComplexSVsign, data = metaD_survMid)
summary(fit)
anova(fit)
ggforest(fit, data = metaD_survMid, fontsize = 1.5) 

check_PH <- cox.zph(fit, transform = "km", terms = TRUE)
check_PH
ggcoxzph(check_PH)


##Effect Plot####
ND <- with(metaD_survMid, expand.grid(
  ComplexSVsign = seq(0, 1, length.out = 25),
  Histone_group = levels(as.factor(Histone_group)),
  TP53snv = levels(as.factor(TP53snv))
))
head(ND)
prs <- predict(fit, newdata = ND, type = "lp", se.fit = TRUE)
ND$pred <- prs[[1]]
ND$se <- prs[[2]]
ND$lo <- ND$pred - 1.96 * ND$se
ND$up <- ND$pred + 1.96 * ND$se
ND$TP53lb <- factor(ND$TP53snv, labels = paste("TP53snv =", sort(unique(ND$TP53snv))))
xyplot(pred + lo + up ~ ComplexSVsign | TP53lb * Histone_group, data = ND, 
       type = "l", col = "black", lwd = 2, lty = c(1, 2, 2),
       abline = list(h = 0, lty = 2, lwd = 2, col = "red"),
       xlab = "ComplexSVsign", ylab = "log Hazard Ratio")

#####

fit_noSigna <- coxph(Surv(Overall_Survival_Months, status) ~  Histone_group  + TP53snv  , data = metaD_survMid)
summary(fit_noSigna)
anova(fit_noSigna)
anova(fit_noSigna, fit)

# TP53 hetloss#########
colnames(metaD_surv)[47] = 'TP53hetloss'
metaD_survMid <- metaD_surv[Location == 'Midline']
fitTP53allni <- coxph(Surv(Overall_Survival_Months, status) ~  Histone_group + Age + TP53snv + TP53hetloss + ComplexSVsign , data = metaD_survMid)
summary(fitTP53allni)
anova(fitTP53allni)
ggforest(fitTP53allni, data = metaD_survMid, fontsize = 1.5) 


# Age #########
metaD_surv$Age <- as.numeric(metaD_surv$Age)
metaD_surv[Age < 1, AgeGroup := '<1']
metaD_surv[Age >= 1 & Age < 3 , AgeGroup := '1-3']
metaD_surv[Age >= 3 & Age < 10, AgeGroup := '3-10']
metaD_surv[Age > 10, AgeGroup := '>10']

metaD_survMid <- metaD_surv[Location == 'Midline']
fitTP53allniAgeCat <- coxph(Surv(Overall_Survival_Months, status) ~  Histone_group + AgeGroup + TP53snv + TP53hetloss + ComplexSVsign , data = metaD_survMid)
summary(fitTP53allniAgeCat)
anova(fitTP53allniAgeCat)
ggforest(fitTP53allniAgeCat, data = metaD_survMid, fontsize = 1.5) 

checkPH_TP53allniAgeCat <- cox.zph(fitTP53allniAgeCat, transform = "km", terms = TRUE)
checkPH_TP53allniAgeCat
ggcoxzph(checkPH_TP53allniAgeCat)

fitAgeCat <- coxph(Surv(Overall_Survival_Months, status) ~  Histone_group + AgeGroup + TP53snv  + ComplexSVsign , data = metaD_survMid)
summary(fitAgeCat)
anova(fitAgeCat)
ggforest(fitAgeCat, data = metaD_survMid, fontsize = 1.5) 

checkPH_TP53allniAgeCat <- cox.zph(fitTP53allniAgeCat, transform = "km", terms = TRUE)
checkPH_TP53allniAgeCat
ggcoxzph(checkPH_TP53allniAgeCat)
# fit <- coxph(Surv(Overall_Survival_Months, status) ~  Histone_group + ComplexSVsign + TP53snvORdel  , data = metaD_survMid)
# summary(fit)
# anova(fit)

# fit <- coxph(Surv(Overall_Survival_Months, status) ~ Location + Histone_group + ComplexSVsign + TP53snv + TP53hetloss  , data = metaD_surv)
# summary(fit)
# anova(fit)

