
# factors affecting endophytic fungal diversity and communities
###################################
###################################

## Pakages we need for the analysis

library(mvabund)
library(vegan)
library(boral)
library(MASS)

##################################
##################################

## Data input

EndoAbun= read.csv (file="matrix_endo.csv", 
                   header = T, row.names = 1)

EndoMetaData= read.csv (file="metadata_end.csv", header = T, row.names = 1)

# are the rownames matching?
rownames(EndoAbun) == rownames(EndoMetaData)

# setting the explaining predictors as factors

Tissue <- factor(EndoMetaData$SUBSTRATE, levels = c("Leaf", "Branch"))
Time <- factor(EndoMetaData$TIME, levels = c("1","2","3"))
Locality<- factor(EndoMetaData$LOCALITY, levels = c("BISOTON", "HASAN ABAD", "KHOSRO ABAD",
                                                   "KEREND","SORKHE DIZE"))
Temperature<-factor(EndoMetaData$TEMPRATURE, levels = c("25 Degree","4 Degree"))
IR= EndoMetaData$IR ### isolation rate
## Data distribution
hist(EndoMetaData$IR)
hist(log(EndoMetaData$IR))
boxplot(EndoMetaData$IR~EndoMetaData$SUBSTRATE)
boxplot(EndoMetaData$IR~EndoMetaData$TIME)
boxplot(EndoMetaData$IR~EndoMetaData$LOCALITY)
boxplot(EndoMetaData$IR~EndoMetaData$TEMPRATURE)

#########################
###### 1- Isolation rate
#########################
# zero-inflated model for success
# 
#   ZerSuc.m1= zeroinfl(EndoMetaData$SUCCESS ~ Substrate, data = EndoMetaData, dist = "poisson")
# AIC(ZerSuc.m1)
# 
# anova(ZerSuc.m2)
# 
# ZerSuc.m2= zeroinfl(EndoMetaData$SUCCESS ~ Time + Locality + Temperature + Substrate, data = EndoMetaData, dist = "negbin")
# AIC(ZerSuc.m2)

## not working ## I don't enderstand this

### GLM modles selection for IR

Ir.m= glm(IR ~ Tissue*Time + Tissue*Locality+ Tissue*Temperature,
           data = EndoMetaData)
plot(Ir.m)
AIC(Ir.m)
anova(Ir.m, test = "Chisq")
summary(Ir.m)
## don't know what kind of model or family is good for this kind of data non-integer 
#with lots of ZEROs
########################################
# ### 2-Richness and Hill Diversities
########################################
# Remove zero observations for diversity and richness calculation
NotZero = EndoMetaData$SUCCESS > 0 #filter for zero-observation samples
EndoAbunZero = EndoAbun[NotZero,]
EndoMetaZero = EndoMetaData[NotZero,]
##aranging the factors with the new datasets
LocalityZ<- factor(EndoMetaZero$LOCALITY,levels = c("BISOTON", "HASAN ABAD", "KHOSRO ABAD",
                               "KEREND","SORKHE DIZE"))
TissueZ <- factor(EndoMetaZero$SUBSTRATE, levels = c("Leaf", "Branch"))
TimeZ <- factor(EndoMetaZero$TIME, levels = c("1","2","3"))
TemperatureZ<-factor(EndoMetaZero$TEMPRATURE, levels = c("25 Degree","4 Degree"))
#EndoRichness = specnumber(EndoAbunZero)
#hist(EndoRichness)

EndoHill = renyi(EndoAbunZero, scale=c(0,1,2), hill=T)
Endohill.1 = EndoHill$"0"#this is richness
Endohill.2 = EndoHill$"1"#antilogarithm of the Shannon representing the abundant species
Endohill.3 = EndoHill$"2"#inverse Simpson representing the very abundant species
hist(Endohill.1)
hist(Endohill.2)
hist(Endohill.3)

### MODELS
### I want to see if the time,loclaity and temperature have differensial effects on
# different tissues diversity

############### First hill
EHill1.m= glm(Endohill.1~ TissueZ*LocalityZ+TissueZ*TimeZ+TissueZ*TemperatureZ 
               , data = EndoMetaZero, family = poisson(link = "log"))
AIC(EHill1.m)
par(mfrow= c(2,2))
plot(EHill1.m)
dev.off()

EHill1.m.anova= anova(EHill1.m, test = "Chisq")
EHill1.m.summary= summary(EHill1.m)# why this doesn't show the significance?

# See if glm.nb is a better fit
# EHill1.m.nb= glm.nb(Endohill.1~ TissueZ*LocalityZ+TissueZ*TimeZ+TissueZ*TemperatureZ 
#               , data = EndoMetaZero, link = "log")
# AIC(EHill1.m.nb)
# par(mfrow= c(2,2))
# plot(EHill1.m.nb)
# dev.off()
## These models are not that different but due to the warnings that I get from glm.nb
## I choose the poisson glm model

############### second Hill
## use the same model orders here

EHill2.m= glm(Endohill.2~ TissueZ*LocalityZ+TissueZ*TimeZ+TissueZ*TemperatureZ 
              ,data = EndoMetaZero, family =Gamma(link = "log"))
AIC(EHill2.m)
par(mfrow= c(2,2))
plot(EHill2.m)
dev.off()

EHill2.m.anova= anova(EHill2.m,test = "F")
EHill2.m.summary= summary(EHill2.m)
## try the glm.nb
# Ehill2.m.nb= glm.nb(Endohill.2~ TissueZ*LocalityZ+TissueZ*TimeZ+TissueZ*TemperatureZ 
#                     ,data = EndoMetaZero)
# warnings()

############# Third Hill
EHill3.m= glm(Endohill.3~ TissueZ*LocalityZ+TissueZ*TimeZ+TissueZ*TemperatureZ 
              ,data = EndoMetaZero, family =Gamma(link = "log"))
AIC(EHill3.m)
par(mfrow= c(2,2))
plot(EHill3.m)
dev.off()

EHill3.m.anova= anova(EHill3.m,test = "F")
EHill3.m.summary= summary(EHill3.m)

##############################
#### 3- Community composition
##############################

### Define CORE OTUs

## Summarize OTU observation
TotalCount = apply(EndoAbun,2,sum)

## The average observation of OTUs
MeanCount=apply(EndoAbun,2,function(vec) mean(vec[vec>0]))

## In how many samples is an OTU present?
TotalPresent = apply(EndoAbun,2,function(vec) sum(vec>0))

## The highest number of an OTU in a sample
MaximumCount=apply(EndoAbun,2,max)

## Plotting observation against abundance
plot(TotalPresent, MaximumCount, xlab="OTU Observation",
     ylab="OTU Maximum Abundance", pch=20)

plot(TotalPresent, log(MaximumCount), xlab="OTU Observation",
     ylab="log(OTU Maximum Abundance)", pch=20)

## Create a smoothed trendline
gam1 = gam(log(MaximumCount)~s(TotalPresent))

plot(gam1, residuals=T, shade=T, rug=F, cex=2.6,
     xlab="OTU Observation", ylab="log Mean Abundance") # , xaxp=c(0,150,15)

## keep core OTUs
OTUobserv = TotalPresent > 10
EndoCorAbun = EndoAbun[,OTUobserv]

### name of the Core OTUs
COREOTUS=colnames(EndoCorAbun)

#### Remove the Zero samples from the Core OTU abbundnace object and metadata
IsolSucc = apply(EndoCorAbun,1, sum)
NotZero2= IsolSucc>0
ECorAbunZero = EndoCorAbun[NotZero2,]
ECorMetaZero = EndoMetaData[NotZero2,]
row.names(ECorAbunZero)==row.names(ECorMetaZero)

TissueC<-factor(ECorMetaZero$SUBSTRATE, levels = c("Leaf", "Branch"))
LocalityC<- factor(ECorMetaZero$LOCALITY,levels = c("BISOTON","HASAN ABAD","KHOSRO ABAD",
                                                    "KEREND","SORKHE DIZE"))
TimeC<- factor(ECorMetaZero$TIME, levels = c("1", "2","3"))
TemperatureC<- factor(ECorMetaZero$TEMPRATURE, levels = c("25 Degree","4 Degree"))

### Multispecies Model 

ECrMvabund= mvabund(ECorAbunZero)
plot(ECrMvabund)

EndoMV.m= manyglm (ECrMvabund~TissueC*LocalityC+TissueC*TimeC+TissueC*TemperatureC,
                   data = ECorMetaZero, family = "negative.binomial", show.residuals=T)

plot.manyglm(EndoMV.m)
EndoMV.m.anova= anova.manyglm(EndoMV.m,nBoot=100, test="LR", p.uni="adjusted", 
                              resamp="montecarlo")
EndoMV.m.sum= summary.manyglm(EndoMV.m, nBoot=100, test="LR",p.uni="adjusted", 
                              resamp="montecarlo")

##############################
##### Tissue specifity test
##############################


