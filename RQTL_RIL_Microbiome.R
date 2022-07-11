#This is an example of the workflow for QTL analysis in R/qtl for shoot and root traits measured in the P.hallii RIL mapping population



#loading rQTL library
library(qtl)


#read in data 
cross<-read.cross("csv", dir="",file="TableS2.csv", genotypes=c("FF","HH"))
#Seperating markers that are in the same location
cross<- jittermap(cross)
#Convert a cross to type "riself" (RIL by selfing)
class(cross)[1]<-"riself"
#Imputation of genotype data
cross<-calc.genoprob(cross, step=2,
                     error.prob=0.0001,
                     map.function="kosambi")


# select the trait (Example on Specific Root Length in Corpus Inoculated microbial treatment))

phe<-"SRL_CI"
plot(cross, pheno.col=phe)



#Scantwo analysis to calaculate penalties for stepwise QTL

SRL_CI_1000perms<-scantwo(cross, pheno.col=phe, method="hk", model="normal", n.perm=1000)

save(SRL_CI_1000perms, file="SRL_CI_1000perms.RData")


print(pens<-calc.penalties(SRL_CI_1000perms, alpha=0.1))
#Stepwise QTL using penalties from scantwo permutations

SRL_CI_stepout<-stepwiseqtl(cross, pheno.col=phe, 
                           method="hk", model="normal", penalties=pens, max.qtl=7,
                           refine.locations=TRUE, keeptrace=TRUE, keeplodprofile=TRUE)

plotLodProfile(SRL_CI_stepout)

title (main="SRL, Scantwo, stepwiseqtl")
plot(SRL_CI_stepout)
plotModel(SRL_CI_stepout)

summary(fitqtl(cross, qtl=SRL_CI_stepout, formula=formula(SRL_CI_stepout), pheno.col=phe,
               method="hk", model="normal", get.ests=T, dropone=T))


#Visualizing the effects for epistatic interactions

effectplot(cross, pheno.col=phe, mname1="03@4.0",mname2= "03@51.9" )


#Calculate confidence interval for QTLs
library(qtlTools)


stats<-qtlStats(cross, pheno.col = phe, mod = SRL_CI_stepout)

print(stats)

# Visualizing QTL on the MAP

cis2<-stats[c("phenotype","LOD","chr","pos","lowposition","highposition")]

cis<-cis2[complete.cases(cis2),c("phenotype","LOD","chr",
                                 "pos","lowposition","highposition")]
rainbowcols <- rainbow(1)
with(cis, segmentsOnMap(cross, phe = phenotype, chr = chr, l = lowposition, h = highposition,col=rainbowcols,
                        peaklod = LOD, peakcM = pos,  showPeaks = TRUE,
                        chrBuffer = c(0.1,0.2) , tick.width=0.05,lwd=3,
                        leg.inset=.55, legendCex=0.6, legendPosition="topleft"))





#Loop for Visualizing QTL on the MAP (include all traits)

#loading qtlTools library
library(qtlTools)

# Selecting the list of traits (here named as X1, X2, X3, X4 ...)
phes <- c("X1",
          "X2",
          "X3",
          "X4")

#loading previously saved scantwo permutation files for each trait

perm.files = c("X1perms.RData", 
               "X2perms.RData",
               "X3perms.RData",
               "X4perms.RData")

for(i in perm.files) load(i)

# add all perm objects in same order as "phes"

perms = list(  X1perms,
               X2perms,
               X3perms,
               X4perms)


#Loop for Stepwise QTL using penalties from scantwo permutations (include all previously selected traits)

stats.list = lapply(1:length(phes), function(x){
  print(phes[x])
  print(pens<-calc.penalties(perms[[x]], alpha=0.1))
  mod<-stepwiseqtl(cross, pheno.col=phes[x], 
                   method="hk", model="normal", penalties=pens, max.qtl=7,
                   additive.only=FALSE, refine.locations=TRUE,
                   keeptrace=TRUE, keeplodprofile=TRUE)
  stats<-qtlStats(cross, pheno.col = phes[x], mod = mod)
  return(stats)
})

stats.df<-do.call(rbind, stats.list)

print(stats.df)

cis2<-stats.df[c("phenotype","LOD","chr","pos","lowposition","highposition")]

cis<-cis2[complete.cases(cis2),c("phenotype","LOD","chr",
                                 "pos","lowposition","highposition")]

print(cis)


#Visualize the positions of QTL on a genetic map.
with(cis, segmentsOnMap(cross, phe = phenotype, chr = chr, l = lowposition, h = highposition,
                        peaklod = LOD, peakcM = pos,  showPeaks = TRUE,
                        chrBuffer = c(0.1,0.2) , tick.width=0.05,lwd=3,
                        leg.inset=.55, legendCex=0.6, legendPosition="topleft"))



