wants <- c('psych', 'dplyr','readxl','psy','olsrr','MVN','parameters')
has <- wants %in% rownames(installed.packages())
if (any(!has)) install.packages(wants[!has])
options("mc.cores"=parallel::detectCores())

knitr::opts_chunk$set(echo = TRUE)
defaultW <- getOption("warn")
options(warn = -1)

# Prerequisites

## Loading libs

library(readr)
library(dplyr)
library(psych)
library(lavaan)
library(olsrr)
library(semTools)
library(mirt)

## Loading data

raw_data <- read.csv("data-dfs.csv")
dfs <- as.data.frame(sapply(select(raw_data, starts_with("Q")), as.integer))

## Checking assumptions 

### Evaluating normality

(mvn_mod <- MVN::mvn(dfs))

### Performing Kaiser, Meyer, Olkin (KMO) Measure of Sampling Adequacy 

(kmo_mod <- KMO(dfs)) 

### Summarizing descriptives and KMOs

t(round(cbind(mvn_mod$Descriptives[,c("Mean","Std.Dev","Median","Skew","Kurtosis")],"MSA"=kmo_mod$MSAi),2))


### Performing Bartlett's test of sphericity

(parameters::check_sphericity(dfs)) # *bartlett.test()* performs test of homogeneity of variances

# Construct validity

## Structural validity

### CFA in the multicorrelated model

Modelo1 <- '
CSB =~ Q1 + Q10 + Q19 + Q28
MAA =~ Q2 + Q11 + Q20 + Q29
CG  =~ Q3 + Q12 + Q21 + Q30
UF  =~ Q4 + Q13 + Q22 + Q31
CTH =~ Q5 + Q14 + Q23 + Q32
SC  =~ Q6 + Q15 + Q24 + Q33
LSC =~ Q7 + Q16 + Q25 + Q34
TT  =~ Q8 + Q17 + Q26 + Q35
AE  =~ Q9 + Q18 + Q27 + Q36

CSB ~~ MAA
CSB ~~ CG
CSB ~~ UF
CSB ~~ CTH
CSB ~~ SC
CSB ~~ LSC
CSB ~~ TT
CSB ~~ AE

MAA ~~ CG
MAA ~~ UF
MAA ~~ CTH
MAA ~~ SC
MAA ~~ LSC
MAA ~~ TT
MAA ~~ AE

CG ~~ UF
CG ~~ CTH
CG ~~ SC
CG ~~ LSC
CG ~~ TT
CG ~~ AE

UF ~~ CTH
UF ~~ SC
UF ~~ LSC
UF ~~ TT
UF ~~ AE

CTH ~~ SC
CTH ~~ LSC
CTH ~~ TT
CTH ~~ AE

SC ~~ LSC
SC ~~ TT
SC ~~ AE

LSC ~~ TT
LSC ~~ AE

TT ~~ AE
'

fit1 <-cfa(Modelo1, data=dfs,estimator="WLSMV", std.lv=TRUE)
standardizedSolution(fit1)

Modelo1a <- '
CSB =~ Q1 + Q10 + Q19 + Q28
MAA =~ Q2 + Q11 + Q20 + Q29
CG  =~ Q3 + Q12 + Q21 + Q30
UF  =~ Q4 + Q13 + Q22 + Q31
CTH =~ Q5 + Q14 + Q23 + Q32
SC  =~ Q6 + Q15 + Q24 + Q33
LSC =~ Q7 + Q16 + Q25 + Q34
TT  =~ Q8 + Q17 + Q26 + Q35
AE  =~ Q9 + Q18 + Q27 + Q36

CSB ~~ 1*MAA
CSB ~~ 1*CG
CSB ~~ UF
CSB ~~ CTH
CSB ~~ SC
CSB ~~ 1*LSC
CSB ~~ 1*TT
CSB ~~ 1*AE

MAA ~~ 1*CG
MAA ~~ UF
MAA ~~ CTH
MAA ~~ SC
MAA ~~ LSC
MAA ~~ 1*TT
MAA ~~ 1*AE

CG ~~ 1*UF
CG ~~ CTH
CG ~~ SC
CG ~~ LSC
CG ~~ TT
CG ~~ AE

UF ~~ CTH
UF ~~ SC
UF ~~ LSC
UF ~~ TT
UF ~~ AE

CTH ~~ 1*SC
CTH ~~ 1*LSC
CTH ~~ 1*TT
CTH ~~ AE

SC ~~ 1*LSC
SC ~~ 1*TT
SC ~~ AE

LSC ~~ 1*TT
LSC ~~ 1*AE

TT ~~ 1*AE
'

fit1a <-cfa(Modelo1a, data=dfs,estimator="WLSMV", std.lv=TRUE)
standardizedSolution(fit1a)

### CFA in the 2nd-order model

Modelo2 <- '
CSB =~ Q1 + Q10 + Q19 + Q28
MAA =~ Q2 + Q11 + Q20 + Q29
CG  =~ Q3 + Q12 + Q21 + Q30
UF  =~ Q4 + Q13 + Q22 + Q31
CTH =~ Q5 + Q14 + Q23 + Q32
SC  =~ Q6 + Q15 + Q24 + Q33
LSC =~ Q7 + Q16 + Q25 + Q34
TT  =~ Q8 + Q17 + Q26 + Q35
AE  =~ Q9 + Q18 + Q27 + Q36
F1  =~ CSB + MAA + CG + UF + CTH + SC + LSC + AE

CSB ~~ 0*MAA
CSB ~~ 0*CG
CSB ~~ 0*UF
CSB ~~ 0*CTH
CSB ~~ 0*SC
CSB ~~ 0*LSC
CSB ~~ 0*TT
CSB ~~ 0*AE

MAA ~~ 0*CG
MAA ~~ 0*UF
MAA ~~ 0*CTH
MAA ~~ 0*SC
MAA ~~ 0*LSC
MAA ~~ 0*TT
MAA ~~ 0*AE

CG ~~ 0*UF
CG ~~ 0*CTH
CG ~~ 0*SC
CG ~~ 0*LSC
CG ~~ 0*TT
CG ~~ 0*AE

UF ~~ 0*CTH
UF ~~ 0*SC
UF ~~ 0*LSC
UF ~~ 0*TT
UF ~~ 0*AE

CTH ~~ 0*SC
CTH ~~ 0*LSC
CTH ~~ 0*TT
CTH ~~ 0*AE

SC ~~ 0*LSC
SC ~~ 0*TT
SC ~~ 0*AE

LSC ~~ 0*TT
LSC ~~ 0*AE

TT ~~ 0*AE
' 

fit2 <-cfa(Modelo2, data=dfs,estimator="WLSMV", std.lv=TRUE)
standardizedSolution(fit2)

### CFA in the short version of DFS-BR (DFS-short)

Modelo3 <- 'DFS =~ Q19 + Q29 + Q12 + Q22 + Q32 + Q6 + Q7 + Q17 + Q36'

fit3 <-cfa(Modelo3, data=dfs,estimator="WLSMV", std.lv=TRUE)
standardizedSolution(fit3)

### CFA in an alternative short version of DFS-BR (DFS-short2)

Modelo4 <- 'DFS =~ Q10 + Q20 + Q21 + Q4 + Q32 + Q24 + Q7 + Q17 + Q36'

fit4 <-cfa(Modelo4, data=dfs,estimator="WLSMV", std.lv=TRUE)
standardizedSolution(fit4)

sl <- standardizedSolution(fit4)
sl <- sl$est.std[sl$op == "=~"]
(describe(sl))

### Summarize the fit indexes and draw the diagrams of CFA models

df_fit <- do.call(rbind, lapply(c(fit1a,fit2,fit3,fit4), FUN = function(fit) {
  dat <- as.list(round(fitMeasures(fit, c("chisq","df","rni","cfi","tli","srmr"
                                          ,"rmsea","rmsea.ci.lower","rmsea.ci.upper")), 3))
  rbind(c(dat[c("chisq","df")],"chisq/df"=round(dat$chisq/dat$df,3)
          , dat[c("rni","cfi","tli","srmr","rmsea")]
          , "rmsea.ci" = paste0("[",dat$rmsea.ci.lower,"; ",dat$rmsea.ci.upper,"]")))
}))
rownames(df_fit) <- c("multicorr","2nd-order","unidim1", "unidim2")
(df_fit)

anova(fit1a,fit2)


## Internal consistency (reliability)

reliability(fit1a, return.total = T)
reliabilityL2(fit2, 'F1')

reliability(fit4, return.total = T)

Fatores <- list(
  'CSB'   = c("Q1","Q10","Q19","Q28")
  , 'MAA' = c("Q2","Q11","Q20","Q29")
  , 'CG'  = c("Q3","Q12","Q21","Q30")
  , 'UF'  = c("Q4","Q13","Q22","Q31")
  , 'CTH' = c("Q5","Q14","Q23","Q32")
  , 'SC'  = c("Q6","Q15","Q24","Q33")
  , 'LSC' = c("Q7","Q16","Q25","Q34")
  , 'TT'  = c("Q8","Q17","Q26","Q35")
  , 'AE'  = c("Q9","Q18","Q27","Q36")
)

compReliability <- function(fit, lfactors = c(), return.total = F) {
  toReturn <- sapply(lfactors, FUN = function(x) {
    sl <- standardizedSolution(fit)
    sl <- sl$est.std[sl$op == "=~" & sl$rhs %in% x]
    names(sl) <- x
    
    re <- 1 - sl^2
    sum(sl)^2 / (sum(sl)^2 + sum(re))
  })
  if (return.total) {
    sl <- standardizedSolution(fit)
    sl <- sl$est.std[sl$op == "=~"]
    re <- 1 - sl^2
    toReturn <- c(toReturn, total=sum(sl)^2 / (sum(sl)^2 + sum(re)))
  }
  toReturn
}

(compReliability(fit1a, Fatores, return.total = T))
(compReliability(fit4, c(), return.total = T))


df.rel <- rbind(reliability(fit1a, return.total = T),
                "CR"=compReliability(fit1a, Fatores, return.total = T))
(round(df.rel, 3))

## Convergent validity and discriminant validity

convergentDiscriminantValidity <- function(fit, lvn, factors, dat) {
  library(olsrr)
  library(semTools)
  
  CR <- compReliability(fit, factors)
  
  AVE <- reliability(fit)[c("avevar"),]
  
  for (f1 in names(AVE)) dat[[f1]] <- rowSums(dat[,factors[[f1]]])
  dat[['F']] <- rowSums(dat[,names(AVE)])
  
  mdl <- lm(as.formula(paste0('F ~', paste0(names(factors), collapse = '+'))), data = dat)
  VIF <- ols_vif_tol(mdl)$VIF
  VIF.i <- sapply(names(factors), FUN = function(f1) {
    mdl <- lm(as.formula(paste0(f1,' ~ ',paste0(factors[[f1]], collapse = '+'))), data = dat)
    max(ols_vif_tol(mdl)$VIF)
  })
  
  corr.df <- as.table(inspect(fit, "cor.lv"))
  corr.df[upper.tri(corr.df)] <- NA
  htmt.df <- as.table(semTools::htmt(lvn, dfs))
  htmt.df[lower.tri(htmt.df)] <- NA
  
  df <- corr.df
  df[upper.tri(corr.df)] <- htmt.df[upper.tri(df)]
  for (cname in names(AVE)) df[cname,cname] <-  sqrt(AVE[[cname]])
  
  as.data.frame(cbind(CR,AVE, VIF, VIF.i, df))
}

convergentDiscriminantValidity(fit1a, Modelo1a, Fatores, dat = dfs)

# Item quality analysis based on IRT

## IRT analysis in the long version of the DFS-BR

dados_tri<-mirt(dfs, 1, itemtype='graded')

### Estimate parameters

(params <- coef(dados_tri, simplify=TRUE, IRTpars=TRUE))

### Graphs of test information curve

plot(dados_tri, type='infotrace')
plot(dados_tri, type='info')
plot(dados_tri, type='infoSE')

## IRT analysis in the short version of the DFS-BR

dados_tri2 <-mirt(dfs[,c("Q10","Q20","Q21","Q4","Q32","Q24","Q7","Q17","Q36")], 1, itemtype='graded')

### Estimate parameters

(params2 <- coef(dados_tri2, simplify=TRUE, IRTpars=TRUE))

### Graphs of test information curve

plot(dados_tri2,type='infotrace')
plot(dados_tri2,type='info')
plot(dados_tri2,type='infoSE')
