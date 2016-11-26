
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

### GLM model for IR

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
OTUobserv = TotalPresent > 5
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

### Multispecies Model for Core OTUs
ECrMvabund= mvabund(ECorAbunZero)
plot(ECrMvabund)

EndoMV.m= manyglm (ECrMvabund~TissueC*LocalityC+TissueC*TimeC+TissueC*TemperatureC,
                   data = ECorMetaZero, family = "negative.binomial", show.residuals=T)

plot.manyglm(EndoMV.m)
EndoMV.m.anova= anova.manyglm(EndoMV.m,nBoot=100, test="LR", p.uni="adjusted", 
                              resamp="montecarlo")
EndoMV.m.sum= summary.manyglm(EndoMV.m, nBoot=100, test="LR",p.uni="adjusted", 
                              resamp="montecarlo")
## Which OTUs are significaantly affected
EnAnova <- as.data.frame(EndoMV.m.anova$uni.p)
otuTissueEf<-colnames(EnAnova)[EnAnova["TissueC",]<= 0.05]
otuLocEf<- colnames(EnAnova)[EnAnova["LocalityC",]<= 0.05]
otuTimEf<- colnames(EnAnova)[EnAnova["TimeC",]<= 0.05]
otuTempEf<- colnames(EnAnova)[EnAnova["TemperatureC",]<= 0.05]
### try to visualize these effects

## Tissue effects
OTUtissu<- c("Byssochlamys.spectabilis.","Gnomoniaceae.sp..66","Microsphaeriopsis.olivacea",
              "Penicillium.sp..A21","Preussia.sp..A31")
TissuABUN<- ECorAbunZero[OTUtissu]##Keeping only tissue affected OTUs
# get the mean valuse for each OTU in each tisse
Tissuemean <- aggregate(. ~ ECorMetaZero$SUBSTRATE, TissuABUN, mean)
#CReat a data frame of the mean valuse
TissuMeanfram<- as.data.frame(Tissuemean,optional=TRUE)
attr(TissuMeanfram, "row.names")<- c("Branch", "Leaf")

### creat a matrix of mean observation of OTUs affected by tissue for ploting
Tissudata<- data.matrix (TissuMeanfram[2:6],rownames.force = NA )

pdf(file = "Effect of Tissue on OTU observation.pdf", paper = "a4", width = 7, height = 4)
barplot(Tissudata, legend.text =TRUE, beside = TRUE,ylab= "mean observation per sample",
        names.arg= c("B. spectabilis", "Gnomoniaceae sp.", "M. olivacea","Penicillium sp.",
                     "Preussia sp."), axes= TRUE,ylim= c(0,1),  cex.names = 0.8,
              args.legend = list(x = "topright",bty= "n"), border = "Black" )
dev.off()

### Temprature effects                                                                                       
OTUtemp<- c ("Alternaria.sp..A25","Alternaria.sp..A9","Aspergillus.sp..A20",
            "Aureobasidium.sp..A17","Byssochlamys.spectabilis.","Microsphaeriopsis.olivacea",
            "Preussia.sp..A31")
TempABUN<- ECorAbunZero[OTUtemp]##Keeping only Temp affected OTUs
# get the mean valuse for each OTU in each temp
Tempmean <- aggregate(. ~ ECorMetaZero$TEMPRATURE, TempABUN, mean)
#CReat a data frame of the mean valuse
TempMeanfram<- as.data.frame(Tempmean,optional=TRUE)
attr(TempMeanfram, "row.names")<- c("25 Degree", "4 Degree")

### creat a matrix of mean observation of OTUs affected by temp for ploting
Tempdata<- data.matrix (TempMeanfram[2:8],rownames.force = NA )
pdf(file = "Effect of Temprature on OTU observation.pdf", paper = "a4", width = 7, height = 4)
barplot(Tempdata,legend.text =TRUE, beside = TRUE,ylab= "mean observation per sample" ,
        names.arg= c ("A25","A9","A20",
                      "A17","B.spec","M.oliv",
                      "A31"), axes= TRUE,ylim= c(0,1.6), cex.names = 0.8,
        args.legend = list(x = "topleft",bty= "n"), border = "Black", 
        width = 0.5)
dev.off()

## Locality effects
plot(ECorAbunZero$Byssochlamys.spectabilis.~ ECorMetaZero$LOCALITY)

# time effects 
plot(ECorAbunZero$Alternaria.sp..A25~ ECorMetaZero$TIME)

################################
#### 5- NMDS and similarities
################################
ENdoNMDS<-metaMDS(ECorAbunZero, distance = "bray", k= 2, trymax = 20)
## Plot NMDS for localities
dev.off()
plot(ENdoNMDS$points, xlab="dimension 1", ylab="dimension 2")
ordiplot(ENdoNMDS, type = "n", display = "sites",xlab="Dimension 1", ylab="Dimension 2"
         )
points(ENdoNMDS$points, pch=20, col= "black", cex=0.5)

with(ECorMetaZero,ordiellipse(ENdoNMDS,ECorMetaZero$LOCALITY,cex=.5, 
                       draw="polygon", col="blue",
                       alpha=100,kind="se",conf=0.95, 
                       show.groups=(c("HASAN ABAD"))))
with(ECorMetaZero,ordiellipse(ENdoNMDS, ECorMetaZero$LOCALITY,cex=.5, 
                 draw="polygon", col=c("black"),
                   alpha=100,kind="se",conf=0.95, 
                   show.groups=(c("KHOSRO ABAD"))))
with(ECorMetaZero,ordiellipse(ENdoNMDS,ECorMetaZero$LOCALITY,cex=.5, 
                draw="polygon", col=c("green"),
                  alpha=100,kind="se",conf=0.95, 
                show.groups=(c("KEREND"))))
with(ECorMetaZero,ordiellipse(ENdoNMDS, ECorMetaZero$LOCALITY,cex=.5, 
            draw="polygon", col=c("yellow"),
                      alpha=100,kind="se",conf=0.95, 
                   show.groups=(c("SORKHE DIZE"))))
with(ECorMetaZero,ordiellipse(ENdoNMDS, ECorMetaZero$LOCALITY,cex=1.5,
                              draw="polygon", col= "darkred",
                              alpha=100, kind="se", conf=0.95, 
                              show.groups=(c("BISOTON"))))#red
mylegend = legend("topright", c("BISOTON","HASAN ABAD","KHOSRO ABAD",
                                 "KEREND","SORKHE DIZE"), cex=0.8,
                  fill=c("red","blue","black","green","yellow"), border="white", bty="n")

## I can't see the differences in this plot
## try 3D NMDS
NMDS.3<-metaMDS(ECorAbunZero, distance = "bray", k= 3, trymax = 20)

### Plot nmds1 &2
dev.off()
#pdf(file = "3D NMDS for localities.pdf", paper = "a4", width = 7, height = 4)
par(mfrow= c(1,3))
NMDS1.2=ordiplot(NMDS.3,choices=c(1,2), type = "n", display = "sites",xlab="Dimension 1", ylab="Dimension 2"
         ,ylim = c(-2,2), xlim = c(-3,3))
points(NMDS.3$points[,1],NMDS.3$points[,2], pch=20, col= "black", cex= 0.3)
with(ECorMetaZero,ordiellipse(NMDS1.2,ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col="blue",
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("HASAN ABAD"))))
with(ECorMetaZero,ordiellipse(NMDS1.2, ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col=c("black"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("KHOSRO ABAD"))))
with(ECorMetaZero,ordiellipse(NMDS1.2,ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col=c("green"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("KEREND"))))
with(ECorMetaZero,ordiellipse(NMDS1.2, ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col=c("yellow"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("SORKHE DIZE"))))
with(ECorMetaZero,ordiellipse(NMDS1.2, ECorMetaZero$LOCALITY,cex=1.5,
                              draw="polygon", col= "darkred",
                              alpha=100, kind="se", conf=0.95, 
                              show.groups=(c("BISOTON"))))#red
## plot nmds2&3
NMDS2.3=ordiplot(NMDS.3,choices=c(2,3), type = "n", display = "sites",xlab="Dimension 2", 
                 ylab="Dimension 3"
                 ,ylim = c(-1.5,1.5), xlim = c(-2,2))
points(NMDS.3$points[,2],NMDS.3$points[,3], pch=20, col= "gray", cex=0.3)
with(ECorMetaZero,ordiellipse(NMDS2.3,ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col="blue",
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("HASAN ABAD"))))
with(ECorMetaZero,ordiellipse(NMDS2.3, ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col=c("black"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("KHOSRO ABAD"))))
with(ECorMetaZero,ordiellipse(NMDS2.3,ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col=c("green"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("KEREND"))))
with(ECorMetaZero,ordiellipse(NMDS2.3, ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col=c("yellow"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("SORKHE DIZE"))))
with(ECorMetaZero,ordiellipse(NMDS2.3, ECorMetaZero$LOCALITY,cex=1.5,
                              draw="polygon", col= "darkred",
                              alpha=100, kind="se", conf=0.95, 
                              show.groups=(c("BISOTON"))))#red
### plot nmds 1&3
NMDS1.3=ordiplot(NMDS.3,choices=c(1,3), type = "n", display = "sites",xlab="Dimension 1", 
                 ylab="Dimension 3"
                 ,ylim = c(-1.5,1.5), xlim = c(-2,2))
points(NMDS.3$points[,1],NMDS.3$points[,3], pch=20, col= "gray",cex=0.3)
with(ECorMetaZero,ordiellipse(NMDS1.3,ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col="blue",
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("HASAN ABAD"))))
with(ECorMetaZero,ordiellipse(NMDS1.3, ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col=c("black"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("KHOSRO ABAD"))))
with(ECorMetaZero,ordiellipse(NMDS1.3,ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col=c("green"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("KEREND"))))
with(ECorMetaZero,ordiellipse(NMDS1.3, ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col=c("yellow"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("SORKHE DIZE"))))
with(ECorMetaZero,ordiellipse(NMDS1.3, ECorMetaZero$LOCALITY,cex=1.5,
                              draw="polygon", col= "darkred",
                              alpha=100, kind="se", conf=0.95, 
                              show.groups=(c("BISOTON"))))#red
mylegend = legend("topright", c("BISOTON","HASAN ABAD","KHOSRO ABAD",
                                "KEREND","SORKHE DIZE"), cex=0.5,
                  fill=c("red","blue","black","green","yellow"), 
                  border="white", bty="n")
dev.off()
####### NMDS for time of sampling

dev.off()
plot(ENdoNMDS$points, xlab="dimension 1", ylab="dimension 2")
ordiplot(ENdoNMDS, type = "n", display = "sites",xlab="Dimension 1", ylab="Dimension 2", 
         ylim = c(-2,2), xlim = c(-3,3))
points(ENdoNMDS$points, pch=20, col= as.numeric(ECorMetaZero$TIME))
#ordispider(ENdoNMDS,ECorMetaZero$TIME, col=c("grey"))

with(ECorMetaZero,ordiellipse(ENdoNMDS,ECorMetaZero$TIME,cex=.5, 
                              draw="polygon", col="black",
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("1"))))
with(ECorMetaZero,ordiellipse(ENdoNMDS, ECorMetaZero$TIME,cex=.5, 
                              draw="polygon", col=c("red"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("2"))))
with(ECorMetaZero,ordiellipse(ENdoNMDS,ECorMetaZero$TIME,cex=.5, 
                              draw="polygon", col=c("green"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("3"))))
legend("topright", c("May 2015","Jun 2015","July 2015"), 
                  fill=c("black","red","green"), 
                  border="white", bty="n")


############################
#### 6- community analysis at family level
############################
## Reporting which factors are affecting the communities at family level
# Data input
### many OTUs were identified at famiy level so I chose this level to do this analysis
### the OTU abundances for each family was merged together
### out of 59 OTUs we were not able to identify 15 OTU and were grouped according to
#the assignment level

FmilyAbun= read.csv (file="matrix_family.csv", 
                    header = T, row.names = 1)

FamilyMVABUND= mvabund(FmilyAbun)
plot(FamilyMVABUND)

FamilyMV.m= manyglm(FamilyMVABUND ~ Locality+ Time+ Temperature+ Tissue, data = EndoMetaData,
                    family = "negative.binomial", show.residuals=T)
FamilyMV.Anova= anova.manyglm(FamilyMV.m,nBoot=100, test="LR", p.uni="adjusted", 
                              resamp="montecarlo" )
FamilyMV.summ= summary.manyglm(FamilyMV.m,nBoot=100, test="LR", p.uni="adjusted", 
                              resamp="montecarlo")
### which families are significantly affected?

FamilyAnova <- as.data.frame(FamilyMV.Anova$uni.p)
FmilyTissue<-colnames(FamilyAnova)[FamilyAnova["Tissue",]<= 0.05]
FmilyLoc<-colnames(FamilyAnova)[FamilyAnova["Locality",]<= 0.05]
FmilyTim<-colnames(FamilyAnova)[FamilyAnova["Time",]<= 0.05]
FmilyTemp<-colnames(FamilyAnova)[FamilyAnova["Temperature",]<= 0.05]








