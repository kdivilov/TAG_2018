library(qtl)

#Single phenotype QTL analysis for RH family breeding values (obtained from linear mixed models)

mycross = read.cross(format=c("csv"),
                     file="RHallphenoandmap.csv",
                     na.strings=("-"),
                     genotypes=NULL,
                     estimate.map = FALSE)

mycross = jittermap(mycross)
#Calculate conditional genotype probabilities. This is needed in order to use method="hk".
mycross = calc.genoprob(mycross, step=1, error.prob=0.001)

names(mycross$pheno)

#1000 permutations to obtain LOD threshold levels
permutations = scanone(mycross,pheno.col=2:4, method="hk", n.perm=1000,n.cluster=3)
#LOD thresholds/penalties from permutations. Note that the penalties will be slightly different for different sets of 1000 permutations.
summary(permutations, alpha=0.05)


####### model creation - sporulation (computer vision)
#Search for QTL using stepwise regression
step1 = stepwiseqtl(mycross,pheno.col=2,method="hk",additive.only = T,penalties = 3.65)
#Fit a model with significant QTL to obtain statistics for the QTL
model = makeqtl(mycross, chr=c(14,30,33), pos=c(54.045,66,70),what="prob")
fqtl = fitqtl(mycross, pheno.col=2, qtl=model, formula=y~Q1+Q2+Q3, method="hk",get.ests = T)
#QTL statistics
summary(fqtl, pvalues=FALSE)

#Approximate Bayesian credible intervals
#Use the find.marker() function to find the closest marker if only a genetic map position is given.
bayesint(step1,qtl.index=1,expandtomarkers=T)
bayesint(step1,qtl.index=2,expandtomarkers=T)
bayesint(step1,qtl.index=3,expandtomarkers=T)


####### model creation - sporulation (manual)

step2 = stepwiseqtl(mycross,pheno.col=3,method="hk",additive.only = T,penalties = 3.59)

model = makeqtl(mycross, chr=c(14,30,33,37), pos=c(59.797,45.08,66.459,47.791),what="prob")
fqtl = fitqtl(mycross, pheno.col=3, qtl=model, formula=y~Q1+Q2+Q3+Q4, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step2,qtl.index=1,expandtomarkers=T)
bayesint(step2,qtl.index=2,expandtomarkers=T)
bayesint(step2,qtl.index=3,expandtomarkers=T)
bayesint(step2,qtl.index=4,expandtomarkers=T)

####### model creation - hypersensitive response (HR)

step3 = stepwiseqtl(mycross,pheno.col=4,method="hk",additive.only = T,penalties = 3.6)

model = makeqtl(mycross, chr=c(11,27,30), pos=c(17.084,42,60.234),what="prob")
fqtl = fitqtl(mycross, pheno.col=4, qtl=model, formula=y~Q1+Q2+Q3, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step3,qtl.index=1,expandtomarkers=T)
bayesint(step3,qtl.index=2,expandtomarkers=T)
bayesint(step3,qtl.index=3,expandtomarkers=T)


#Single phenotype QTL analysis for HC family breeding values (obtained from linear mixed models)

mycross = read.cross(format=c("csv"),
                     file="HCallphenoandmap.csv",
                     na.strings=("-"),
                     genotypes=NULL,
                     estimate.map = FALSE)

mycross = jittermap(mycross)
mycross = calc.genoprob(mycross, step=1, error.prob=0.001)

names(mycross$pheno)

permutations = scanone(mycross,pheno.col=2:5, method="hk", n.perm=1000,n.cluster=3)
summary(permutations, alpha=0.05)


####### model creation - sporulation (computer vision)

step1 = stepwiseqtl(mycross,pheno.col=2,method="hk",additive.only = T,penalties = 3.56)

#NO QTL

####### model creation - sporulation (manual)

step2 = stepwiseqtl(mycross,pheno.col=3,method="hk",additive.only = T,penalties = 3.43)

model = makeqtl(mycross, chr=c(5,7,8), pos=c(5.5857,11,66.7886),what="prob")
fqtl = fitqtl(mycross, pheno.col=3, qtl=model, formula=y~Q1+Q2+Q3, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step2,qtl.index=1,expandtomarkers=T)
bayesint(step2,qtl.index=2,expandtomarkers=T)
bayesint(step2,qtl.index=3,expandtomarkers=T)

####### model creation - hypersensitive response (HR)

step3 = stepwiseqtl(mycross,pheno.col=4,method="hk",additive.only = T,penalties = 3.44)

model = makeqtl(mycross, chr=c(5,6,8), pos=c(4.2434,2.5067,73.3909),what="prob")
fqtl = fitqtl(mycross, pheno.col=4, qtl=model, formula=y~Q1+Q2+Q3, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step3,qtl.index=1,expandtomarkers=T)
bayesint(step3,qtl.index=2,expandtomarkers=T)
bayesint(step3,qtl.index=3,expandtomarkers=T)

####### model creation - leaf trichomes

step4 = stepwiseqtl(mycross,pheno.col=5,method="hk",additive.only = T,penalties = 3.42)

model = makeqtl(mycross, chr=c(5,8,15), pos=c(5.5857,54.7672,63),what="prob")
fqtl = fitqtl(mycross, pheno.col=5, qtl=model, formula=y~Q1+Q2+Q3, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step4,qtl.index=1,expandtomarkers=T)
bayesint(step4,qtl.index=2,expandtomarkers=T)
bayesint(step4,qtl.index=3,expandtomarkers=T)
