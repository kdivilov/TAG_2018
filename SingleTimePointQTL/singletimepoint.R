library(qtl)

#RH 2015 sporulation (computer vision)

mycross = read.cross(format=c("csv"),
                     file="RH2015CV.csv",
                     na.strings=("-"),
                     genotypes=NULL,
                     estimate.map = FALSE)

mycross = jittermap(mycross)
mycross = calc.genoprob(mycross, step=1, error.prob=0.001)

names(mycross$pheno)

permutations = scanone(mycross,pheno.col=2:13, method="hk", n.perm=1000,n.cluster=3)
summary(permutations, alpha=0.05)


####### model creation - exp1-3dpi

step1 = stepwiseqtl(mycross,pheno.col=2,method="hk",additive.only = T,penalties = 3.56)

#NO QTL

####### model creation - exp1-4dpi

step2 = stepwiseqtl(mycross,pheno.col=3,method="hk",additive.only = T,penalties = 3.56)

model = makeqtl(mycross, chr=c(12,14,16,27,33), pos=c(19.404,54.045,52,65.597,69.096),what="prob")
fqtl = fitqtl(mycross, pheno.col=3, qtl=model, formula=y~Q1+Q2+Q3+Q4+Q5, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step2,qtl.index=1,expandtomarkers=T)
bayesint(step2,qtl.index=2,expandtomarkers=T)
bayesint(step2,qtl.index=3,expandtomarkers=T)
bayesint(step2,qtl.index=4,expandtomarkers=T)
bayesint(step2,qtl.index=5,expandtomarkers=T)

####### model creation - exp1-5dpi

step3 = stepwiseqtl(mycross,pheno.col=4,method="hk",additive.only = T,penalties = 3.56)

model = makeqtl(mycross, chr=c(7,12,14,16,27,29,33), pos=c(75.690,20.393,49,52,60,23.250,79.967),what="prob")
fqtl = fitqtl(mycross, pheno.col=4, qtl=model, formula=y~Q1+Q2+Q3+Q4+Q5+Q6+Q7, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step3,qtl.index=1,expandtomarkers=T)
bayesint(step3,qtl.index=2,expandtomarkers=T)
bayesint(step3,qtl.index=3,expandtomarkers=T)
bayesint(step3,qtl.index=4,expandtomarkers=T)
bayesint(step3,qtl.index=5,expandtomarkers=T)
bayesint(step3,qtl.index=6,expandtomarkers=T)
bayesint(step3,qtl.index=7,expandtomarkers=T)

####### model creation - exp1-6dpi

step4 = stepwiseqtl(mycross,pheno.col=5,method="hk",additive.only = T,penalties = 3.67)

model = makeqtl(mycross, chr=c(12,16,29,33), pos=c(21.358,53,55,80),what="prob")
fqtl = fitqtl(mycross, pheno.col=5, qtl=model, formula=y~Q1+Q2+Q3+Q4, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step4,qtl.index=1,expandtomarkers=T)
bayesint(step4,qtl.index=2,expandtomarkers=T)
bayesint(step4,qtl.index=3,expandtomarkers=T)
bayesint(step4,qtl.index=4,expandtomarkers=T)

####### model creation - exp2-3dpi

step5 = stepwiseqtl(mycross,pheno.col=6,method="hk",additive.only = T,penalties = 3.66)

#NO QTL

####### model creation - exp2-4dpi

step6 = stepwiseqtl(mycross,pheno.col=7,method="hk",additive.only = T,penalties = 3.58)

#NO QTL

####### model creation - exp2-5dpi

step7 = stepwiseqtl(mycross,pheno.col=8,method="hk",additive.only = T,penalties = 3.55)

model = makeqtl(mycross, chr=c(14,30), pos=c(42.966,47.793),what="prob")
fqtl = fitqtl(mycross, pheno.col=8, qtl=model, formula=y~Q1+Q2, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step7,qtl.index=1,expandtomarkers=T)
bayesint(step7,qtl.index=2,expandtomarkers=T)

####### model creation - exp2-6dpi

step8 = stepwiseqtl(mycross,pheno.col=9,method="hk",additive.only = T,penalties = 3.7)

model = makeqtl(mycross, chr=c(14,16,33), pos=c(37,54.852,72.433),what="prob")
fqtl = fitqtl(mycross, pheno.col=9, qtl=model, formula=y~Q1+Q2+Q3, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step8,qtl.index=1,expandtomarkers=T)
bayesint(step8,qtl.index=2,expandtomarkers=T)
bayesint(step8,qtl.index=3,expandtomarkers=T)

####### model creation - exp3-3dpi

step9 = stepwiseqtl(mycross,pheno.col=10,method="hk",additive.only = T,penalties = 3.56)

#NO QTL

####### model creation - exp3-4dpi

step10 = stepwiseqtl(mycross,pheno.col=11,method="hk",additive.only = T,penalties = 3.45)

#NO QTL

####### model creation - exp3-5dpi

step11 = stepwiseqtl(mycross,pheno.col=12,method="hk",additive.only = T,penalties = 3.61)

model = makeqtl(mycross, chr=c(14,33,37), pos=c(48.824,70,25.043),what="prob")
fqtl = fitqtl(mycross, pheno.col=12, qtl=model, formula=y~Q1+Q2+Q3, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step11,qtl.index=1,expandtomarkers=T)
bayesint(step11,qtl.index=2,expandtomarkers=T)
bayesint(step11,qtl.index=3,expandtomarkers=T)

####### model creation - exp3-6dpi

step12 = stepwiseqtl(mycross,pheno.col=13,method="hk",additive.only = T,penalties = 3.53)

model = makeqtl(mycross, chr=c(14,33,37), pos=c(54.045,70.841,25.043),what="prob")
fqtl = fitqtl(mycross, pheno.col=13, qtl=model, formula=y~Q1+Q2+Q3, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step12,qtl.index=1,expandtomarkers=T)
bayesint(step12,qtl.index=2,expandtomarkers=T)
bayesint(step12,qtl.index=3,expandtomarkers=T)


#RH 2015 sporulation (manual) and HR

mycross = read.cross(format=c("csv"),
                     file="RH2015man.csv",
                     na.strings=("-"),
                     genotypes=NULL,
                     estimate.map = FALSE)

mycross = jittermap(mycross)
mycross = calc.genoprob(mycross, step=1, error.prob=0.001)

names(mycross$pheno)

permutations = scanone(mycross,pheno.col=2:16, method="hk", n.perm=1000,n.cluster=3)
summary(permutations, alpha=0.05)

####### model creation - exp1-3dpi

step1 = stepwiseqtl(mycross,pheno.col=2,method="hk",additive.only = T,penalties = 2.9)

model = makeqtl(mycross, chr=c(12,18,18,27), pos=c(5,19,21,71.355),what="prob")
fqtl = fitqtl(mycross, pheno.col=2, qtl=model, formula=y~Q1+Q2+Q3+Q4, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step1,qtl.index=1,expandtomarkers=T)
bayesint(step1,qtl.index=2,expandtomarkers=T)
bayesint(step1,qtl.index=3,expandtomarkers=T)
bayesint(step1,qtl.index=4,expandtomarkers=T)

####### model creation - exp1-4dpi

step2 = stepwiseqtl(mycross,pheno.col=3,method="hk",additive.only = T,penalties = 3.5)

model = makeqtl(mycross, chr=c(12,14,27,33), pos=c(19,55,65.597,56.954),what="prob")
fqtl = fitqtl(mycross, pheno.col=3, qtl=model, formula=y~Q1+Q2+Q3+Q4, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step2,qtl.index=1,expandtomarkers=T)
bayesint(step2,qtl.index=2,expandtomarkers=T)
bayesint(step2,qtl.index=3,expandtomarkers=T)
bayesint(step2,qtl.index=4,expandtomarkers=T)

####### model creation - exp1-5dpi

step3 = stepwiseqtl(mycross,pheno.col=4,method="hk",additive.only = T,penalties = 3.57)

#NO QTL

####### model creation - exp1-6dpi

step4 = stepwiseqtl(mycross,pheno.col=5,method="hk",additive.only = T,penalties = 3.47)

#NO QTL

####### model creation - exp2-3dpi

step5 = stepwiseqtl(mycross,pheno.col=6,method="hk",additive.only = T,penalties = 3.27)

model = makeqtl(mycross, chr=c(16), pos=c(55),what="prob")
fqtl = fitqtl(mycross, pheno.col=6, qtl=model, formula=y~Q1, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step5,qtl.index=1,expandtomarkers=T)

####### model creation - exp2-4dpi

step6 = stepwiseqtl(mycross,pheno.col=7,method="hk",additive.only = T,penalties = 3.6)

model = makeqtl(mycross, chr=c(33), pos=c(47.712),what="prob")
fqtl = fitqtl(mycross, pheno.col=7, qtl=model, formula=y~Q1, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step6,qtl.index=1,expandtomarkers=T)

####### model creation - exp2-5dpi

step7 = stepwiseqtl(mycross,pheno.col=8,method="hk",additive.only = T,penalties = 3.62)

model = makeqtl(mycross, chr=c(16), pos=c(54.852),what="prob")
fqtl = fitqtl(mycross, pheno.col=8, qtl=model, formula=y~Q1, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step7,qtl.index=1,expandtomarkers=T)

####### model creation - exp2-6dpi

step8 = stepwiseqtl(mycross,pheno.col=9,method="hk",additive.only = T,penalties = 3.53)

model = makeqtl(mycross, chr=c(16), pos=c(54.852),what="prob")
fqtl = fitqtl(mycross, pheno.col=9, qtl=model, formula=y~Q1, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step8,qtl.index=1,expandtomarkers=T)

####### model creation - exp3-3dpi

step9 = stepwiseqtl(mycross,pheno.col=10,method="hk",additive.only = T,penalties = 3.02)

#NO QTL

####### model creation - exp3-4dpi

step10 = stepwiseqtl(mycross,pheno.col=11,method="hk",additive.only = T,penalties = 3.39)

#NO QTL

####### model creation - exp3-5dpi

step11 = stepwiseqtl(mycross,pheno.col=12,method="hk",additive.only = T,penalties = 3.72)

model = makeqtl(mycross, chr=c(14,37), pos=c(47,25.043),what="prob")
fqtl = fitqtl(mycross, pheno.col=12, qtl=model, formula=y~Q1+Q2, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step11,qtl.index=1,expandtomarkers=T)
bayesint(step11,qtl.index=2,expandtomarkers=T)

####### model creation - exp3-6dpi

step12 = stepwiseqtl(mycross,pheno.col=13,method="hk",additive.only = T,penalties = 3.65)

model = makeqtl(mycross, chr=c(14,32,33,37), pos=c(47.858,81.877,69.096,25.043),what="prob")
fqtl = fitqtl(mycross, pheno.col=13, qtl=model, formula=y~Q1+Q2+Q3+Q4, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step12,qtl.index=1,expandtomarkers=T)
bayesint(step12,qtl.index=2,expandtomarkers=T)
bayesint(step12,qtl.index=3,expandtomarkers=T)
bayesint(step12,qtl.index=4,expandtomarkers=T)

####### model creation - exp1-2dpiHR

step13 = stepwiseqtl(mycross,pheno.col=14,method="hk",additive.only = T,penalties = 3.53)

model = makeqtl(mycross, chr=c(8,27,30), pos=c(25.636,42,64.722),what="prob")
fqtl = fitqtl(mycross, pheno.col=14, qtl=model, formula=y~Q1+Q2+Q3, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step13,qtl.index=1,expandtomarkers=T)
bayesint(step13,qtl.index=2,expandtomarkers=T)
bayesint(step13,qtl.index=3,expandtomarkers=T)

####### model creation - exp2-2dpiHR

step14 = stepwiseqtl(mycross,pheno.col=15,method="hk",additive.only = T,penalties = 3.53)

model = makeqtl(mycross, chr=c(11,24,27,30,36), pos=c(17.0839,5.1009,31.3505,60.2343,62),what="prob")
fqtl = fitqtl(mycross, pheno.col=15, qtl=model, formula=y~Q1+Q2+Q3+Q4+Q5, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step14,qtl.index=1,expandtomarkers=T)
bayesint(step14,qtl.index=2,expandtomarkers=T)
bayesint(step14,qtl.index=3,expandtomarkers=T)
bayesint(step14,qtl.index=4,expandtomarkers=T)
bayesint(step14,qtl.index=5,expandtomarkers=T)

####### model creation - exp3-2dpiHR

step15 = stepwiseqtl(mycross,pheno.col=16,method="hk",additive.only = T,penalties = 3.56)

model = makeqtl(mycross, chr=c(27,30), pos=c(41.802,60.234),what="prob")
fqtl = fitqtl(mycross, pheno.col=16, qtl=model, formula=y~Q1+Q2, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step15,qtl.index=1,expandtomarkers=T)
bayesint(step15,qtl.index=2,expandtomarkers=T)


#RH 2016 sporulation (computer vision)

mycross = read.cross(format=c("csv"),
                     file="RH2016CV.csv",
                     na.strings=("-"),
                     genotypes=NULL,
                     estimate.map = FALSE)

mycross = jittermap(mycross)
mycross = calc.genoprob(mycross, step=1, error.prob=0.001)

names(mycross$pheno)

permutations = scanone(mycross,pheno.col=2:9, method="hk", n.perm=1000,n.cluster=3)
summary(permutations, alpha=0.05)

####### model creation - exp1-3dpi

step1 = stepwiseqtl(mycross,pheno.col=2,method="hk",additive.only = T,penalties = 3.48)

model = makeqtl(mycross, chr=c(21), pos=c(5.0617),what="prob")
fqtl = fitqtl(mycross, pheno.col=2, qtl=model, formula=y~Q1, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step1,qtl.index=1,expandtomarkers=T)

####### model creation - exp1-4dpi

step2 = stepwiseqtl(mycross,pheno.col=3,method="hk",additive.only = T,penalties = 3.68)

model = makeqtl(mycross, chr=c(14), pos=c(53.106),what="prob")
fqtl = fitqtl(mycross, pheno.col=3, qtl=model, formula=y~Q1, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step2,qtl.index=1,expandtomarkers=T)

####### model creation - exp1-5dpi

step3 = stepwiseqtl(mycross,pheno.col=4,method="hk",additive.only = T,penalties = 3.59)

model = makeqtl(mycross, chr=c(14), pos=c(53.106),what="prob")
fqtl = fitqtl(mycross, pheno.col=4, qtl=model, formula=y~Q1+Q2+Q3+Q4+Q5+Q6+Q7, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step3,qtl.index=1,expandtomarkers=T)

####### model creation - exp1-6dpi

step4 = stepwiseqtl(mycross,pheno.col=5,method="hk",additive.only = T,penalties = 3.54)

model = makeqtl(mycross, chr=c(14,16), pos=c(53.106,65.366),what="prob")
fqtl = fitqtl(mycross, pheno.col=5, qtl=model, formula=y~Q1+Q2, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step4,qtl.index=1,expandtomarkers=T)
bayesint(step4,qtl.index=2,expandtomarkers=T)

####### model creation - exp2-4dpi

step5 = stepwiseqtl(mycross,pheno.col=6,method="hk",additive.only = T,penalties = 3.63)

#NO QTL

####### model creation - exp2-5dpi

step6 = stepwiseqtl(mycross,pheno.col=7,method="hk",additive.only = T,penalties = 3.52)

#NO QTL

####### model creation - exp2-6dpi

step7 = stepwiseqtl(mycross,pheno.col=8,method="hk",additive.only = T,penalties = 3.51)

model = makeqtl(mycross, chr=c(7), pos=c(28.207),what="prob")
fqtl = fitqtl(mycross, pheno.col=8, qtl=model, formula=y~Q1, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step7,qtl.index=1,expandtomarkers=T)

####### model creation - exp2-7dpi

step8 = stepwiseqtl(mycross,pheno.col=9,method="hk",additive.only = T,penalties = 3.53)

#NO QTL


#RH 2016 sporulation (manual) and HR

mycross = read.cross(format=c("csv"),
                     file="RH2016man.csv",
                     na.strings=("-"),
                     genotypes=NULL,
                     estimate.map = FALSE)

mycross = jittermap(mycross)
mycross = calc.genoprob(mycross, step=1, error.prob=0.001)

names(mycross$pheno)

permutations = scanone(mycross,pheno.col=2:11, method="hk", n.perm=1000,n.cluster=3)
summary(permutations, alpha=0.05)


####### model creation - exp1-3dpi

step1 = stepwiseqtl(mycross,pheno.col=2,method="hk",additive.only = T,penalties = 3.47)

model = makeqtl(mycross, chr=c(37), pos=c(35.361),what="prob")
fqtl = fitqtl(mycross, pheno.col=2, qtl=model, formula=y~Q1, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step1,qtl.index=1,expandtomarkers=T)

####### model creation - exp1-4dpi

step2 = stepwiseqtl(mycross,pheno.col=3,method="hk",additive.only = T,penalties = 3.72)

model = makeqtl(mycross, chr=c(14), pos=c(53.106),what="prob")
fqtl = fitqtl(mycross, pheno.col=3, qtl=model, formula=y~Q1, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step2,qtl.index=1,expandtomarkers=T)

####### model creation - exp1-5dpi

step3 = stepwiseqtl(mycross,pheno.col=4,method="hk",additive.only = T,penalties = 3.59)

#NO QTL

####### model creation - exp1-6dpi

step4 = stepwiseqtl(mycross,pheno.col=5,method="hk",additive.only = T,penalties = 3.47)

#NO QTL

####### model creation - exp2-4dpi

step5 = stepwiseqtl(mycross,pheno.col=6,method="hk",additive.only = T,penalties = 3.35)

#NO QTL

####### model creation - exp2-5dpi

step6 = stepwiseqtl(mycross,pheno.col=7,method="hk",additive.only = T,penalties = 3.64)

#NO QTL

####### model creation - exp2-6dpi

step7 = stepwiseqtl(mycross,pheno.col=8,method="hk",additive.only = T,penalties = 3.68)

model = makeqtl(mycross, chr=c(28), pos=c(64.708),what="prob")
fqtl = fitqtl(mycross, pheno.col=8, qtl=model, formula=y~Q1, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step7,qtl.index=1,expandtomarkers=T)

####### model creation - exp2-7dpi

step8 = stepwiseqtl(mycross,pheno.col=9,method="hk",additive.only = T,penalties = 3.62)

#NO QTL

####### model creation - exp1-2dpiHR

step9 = stepwiseqtl(mycross,pheno.col=10,method="hk",additive.only = T,penalties = 3.48)

#NO QTL

####### model creation - exp2-2dpiHR

step10 = stepwiseqtl(mycross,pheno.col=11,method="hk",additive.only = T,penalties = 3.51)

#NO QTL


#HC 2015 sporulation (computer vision)

mycross = read.cross(format=c("csv"),
                     file="HC2015CV.csv",
                     na.strings=("-"),
                     genotypes=NULL,
                     estimate.map = FALSE)

mycross = jittermap(mycross)
mycross = calc.genoprob(mycross, step=1, error.prob=0.001)

names(mycross$pheno)

full = read.csv("HC2015CV.csv", header = T)
full = as.matrix(full)
full[full == "-"] = NA
full = full[-c(1,2),]
full = apply(full,2,as.numeric)

#exp1-2dpitrichomes
covar1 = full[,10]
#exp2-2dpitrichomes
covar2 = full[,11]


permutations1 = scanone(mycross,pheno.col=2:5,addcovar = covar1,method="hk", n.perm=1000, n.cluster=3)
summary(permutations1, alpha=0.05)
permutations2 = scanone(mycross,pheno.col=6:9,addcovar = covar2, method="hk", n.perm=1000, n.cluster=3)
summary(permutations2, alpha=0.05)


####### model creation - exp1-4dpi

step1 = stepwiseqtl(mycross,pheno.col=2,method="hk",additive.only = T,penalties = 3.48,covar = as.data.frame(covar1))

#NO QTL

####### model creation - exp1-5dpi

step2 = stepwiseqtl(mycross,pheno.col=3,method="hk",additive.only = T,penalties = 3.5,covar = as.data.frame(covar1))

#NO QTL

####### model creation - exp1-6dpi

step3 = stepwiseqtl(mycross,pheno.col=4,method="hk",additive.only = T,penalties = 3.47,covar = as.data.frame(covar1))

#NO QTL

####### model creation - exp1-7dpi

step4 = stepwiseqtl(mycross,pheno.col=5,method="hk",additive.only = T,penalties = 3.45,covar = as.data.frame(covar1))

model = makeqtl(mycross, chr=c(7), pos=c(85.892),what="prob")
fqtl = fitqtl(mycross, pheno.col=5, qtl=model, formula=y~covar1+Q1, method="hk",get.ests = T,covar = as.data.frame(covar1))
summary(fqtl, pvalues=FALSE)

bayesint(step4,qtl.index=1,expandtomarkers=T)

####### model creation - exp2-4dpi

step5 = stepwiseqtl(mycross,pheno.col=6,method="hk",additive.only = T,penalties = 3.3,covar = as.data.frame(covar2))

model = makeqtl(mycross, chr=c(38,38), pos=c(1,6.0856),what="prob")
fqtl = fitqtl(mycross, pheno.col=6, qtl=model, formula=y~covar2+Q1+Q2, method="hk",get.ests = T,covar = as.data.frame(covar2))
summary(fqtl, pvalues=FALSE)

bayesint(step5,qtl.index=1,expandtomarkers=T)
bayesint(step5,qtl.index=2,expandtomarkers=T)

#Neither of the two loci found pass the LOD threshold. Therefore, we do not consider them QTL.

####### model creation - exp2-5dpi

step6 = stepwiseqtl(mycross,pheno.col=7,method="hk",additive.only = T,penalties = 3.46,covar = as.data.frame(covar2))

#NO QTL

####### model creation - exp2-6dpi

step7 = stepwiseqtl(mycross,pheno.col=8,method="hk",additive.only = T,penalties = 3.62,covar = as.data.frame(covar2))

#NO QTL

####### model creation - exp2-7dpi

step8 = stepwiseqtl(mycross,pheno.col=9,method="hk",additive.only = T,penalties = 3.33,covar = as.data.frame(covar2))

model = makeqtl(mycross, chr=c(18), pos=c(35.015),what="prob")
fqtl = fitqtl(mycross, pheno.col=9, qtl=model, formula=y~covar2+Q1, method="hk",get.ests = T,covar = as.data.frame(covar2))
summary(fqtl, pvalues=FALSE)

bayesint(step8,qtl.index=1,expandtomarkers=T)


#HC 2016 cv

mycross = read.cross(format=c("csv"),
                     file="HC2016CV.csv",
                     na.strings=("-"),
                     genotypes=NULL,
                     estimate.map = FALSE)

mycross = jittermap(mycross)
mycross = calc.genoprob(mycross, step=1, error.prob=0.001)

names(mycross$pheno)

full = read.csv("HC2016CV.csv", header = T)
full = as.matrix(full)
full[full == "-"] = NA
full = full[-c(1,2),]
full = apply(full,2,as.numeric)

#exp1-2dpitrichomes
covar1 = full[,10]
#exp2-2dpitrichomes
covar2 = full[,11]

permutations1 = scanone(mycross,pheno.col=2:5,addcovar = covar1,method="hk", n.perm=1000, n.cluster=3)
summary(permutations1, alpha=0.05)
permutations2 = scanone(mycross,pheno.col=6:9,addcovar = covar2, method="hk", n.perm=1000, n.cluster=3)
summary(permutations2, alpha=0.05)


####### model creation - exp1-4dpi

step1 = stepwiseqtl(mycross,pheno.col=2,method="hk",additive.only = T,penalties = 3.1,covar = as.data.frame(covar1))

#NO QTL

####### model creation - exp1-5dpi

step2 = stepwiseqtl(mycross,pheno.col=3,method="hk",additive.only = T,penalties = 3.27,covar = as.data.frame(covar1))

#NO QTL

####### model creation - exp1-6dpi

step3 = stepwiseqtl(mycross,pheno.col=4,method="hk",additive.only = T,penalties = 3.25,covar = as.data.frame(covar1))

#NO QTL

####### model creation - exp1-7dpi

step4 = stepwiseqtl(mycross,pheno.col=5,method="hk",additive.only = T,penalties = 3.18,covar = as.data.frame(covar1))

#NO QTL

####### model creation - exp2-4dpi

step5 = stepwiseqtl(mycross,pheno.col=6,method="hk",additive.only = T,penalties = 3.47,covar = as.data.frame(covar2))

#NO QTL

####### model creation - exp2-5dpi

step6 = stepwiseqtl(mycross,pheno.col=7,method="hk",additive.only = T,penalties = 3.45,covar = as.data.frame(covar2))

#NO QTL

####### model creation - exp2-6dpi

step7 = stepwiseqtl(mycross,pheno.col=8,method="hk",additive.only = T,penalties = 3.54,covar = as.data.frame(covar2))

#NO QTL

####### model creation - exp2-7dpi

step8 = stepwiseqtl(mycross,pheno.col=9,method="hk",additive.only = T,penalties = 3.51,covar = as.data.frame(covar2))

#NO QTL


#HC 2015 manual

mycross = read.cross(format=c("csv"),
                     file="HC2015man.csv",
                     na.strings=("-"),
                     genotypes=NULL,
                     estimate.map = FALSE)

mycross = jittermap(mycross)
mycross = calc.genoprob(mycross, step=1, error.prob=0.001)

names(mycross$pheno)

permutations = scanone(mycross,pheno.col=2:10, method="hk", n.perm=1000,n.cluster=3)
summary(permutations, alpha=0.05)


####### model creation - exp1-4dpi

step1 = stepwiseqtl(mycross,pheno.col=2,method="hk",additive.only = T,penalties = 3.6)

model = makeqtl(mycross, chr=c(5,8), pos=c(5,73.391),what="prob")
fqtl = fitqtl(mycross, pheno.col=2, qtl=model, formula=y~Q1+Q2, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step1,qtl.index=1,expandtomarkers=T)
bayesint(step1,qtl.index=2,expandtomarkers=T)


####### model creation - exp1-5dpi

step2 = stepwiseqtl(mycross,pheno.col=3,method="hk",additive.only = T,penalties = 3.62)

model = makeqtl(mycross, chr=c(5,8), pos=c(5,54.767),what="prob")
fqtl = fitqtl(mycross, pheno.col=3, qtl=model, formula=y~Q1+Q2, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step2,qtl.index=1,expandtomarkers=T)
bayesint(step2,qtl.index=2,expandtomarkers=T)


####### model creation - exp1-6dpi

step3 = stepwiseqtl(mycross,pheno.col=4,method="hk",additive.only = T,penalties = 3.54)

model = makeqtl(mycross, chr=c(7,8,34), pos=c(11,66,39),what="prob")
fqtl = fitqtl(mycross, pheno.col=4, qtl=model, formula=y~Q1+Q2+Q3, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step3,qtl.index=1,expandtomarkers=T)
bayesint(step3,qtl.index=2,expandtomarkers=T)
bayesint(step3,qtl.index=3,expandtomarkers=T)


####### model creation - exp1-7dpi

step4 = stepwiseqtl(mycross,pheno.col=5,method="hk",additive.only = T,penalties = 3.53)

model = makeqtl(mycross, chr=c(5,8), pos=c(5,54.767),what="prob")
fqtl = fitqtl(mycross, pheno.col=5, qtl=model, formula=y~Q1+Q2, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step4,qtl.index=1,expandtomarkers=T)
bayesint(step4,qtl.index=2,expandtomarkers=T)

####### model creation - exp2-4dpi

step5 = stepwiseqtl(mycross,pheno.col=6,method="hk",additive.only = T,penalties = 3.33)

model = makeqtl(mycross, chr=c(7,10), pos=c(57.126,14.946),what="prob")
fqtl = fitqtl(mycross, pheno.col=6, qtl=model, formula=y~Q1+Q2, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step5,qtl.index=1,expandtomarkers=T)
bayesint(step5,qtl.index=2,expandtomarkers=T)

#Neither of the two loci found pass the LOD threshold. Therefore, we do not consider them QTL.

####### model creation - exp2-5dpi

step6 = stepwiseqtl(mycross,pheno.col=7,method="hk",additive.only = T,penalties = 3.5)

model = makeqtl(mycross, chr=c(7), pos=c(9.3215),what="prob")
fqtl = fitqtl(mycross, pheno.col=7, qtl=model, formula=y~Q1, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step6,qtl.index=1,expandtomarkers=T)

####### model creation - exp2-6dpi

step7 = stepwiseqtl(mycross,pheno.col=8,method="hk",additive.only = T,penalties = 3.56)

model = makeqtl(mycross, chr=c(7), pos=c(9.9252),what="prob")
fqtl = fitqtl(mycross, pheno.col=8, qtl=model, formula=y~Q1, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step7,qtl.index=1,expandtomarkers=T)

####### model creation - exp2-7dpi

step8 = stepwiseqtl(mycross,pheno.col=9,method="hk",additive.only = T,penalties = 3.56)

model = makeqtl(mycross, chr=c(7), pos=c(9.9252),what="prob")
fqtl = fitqtl(mycross, pheno.col=9, qtl=model, formula=y~Q1, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step8,qtl.index=1,expandtomarkers=T)

####### model creation - exp2-2dpiHR

step9 = stepwiseqtl(mycross,pheno.col=10,method="hk",additive.only = T,penalties = 3.56)

model = makeqtl(mycross, chr=c(5,6), pos=c(5.5857,9.7241),what="prob")
fqtl = fitqtl(mycross, pheno.col=10, qtl=model, formula=y~Q1+Q2, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step9,qtl.index=1,expandtomarkers=T)
bayesint(step9,qtl.index=2,expandtomarkers=T)


#HC 2016 manual

mycross = read.cross(format=c("csv"),
                     file="HC2016man.csv",
                     na.strings=("-"),
                     genotypes=NULL,
                     estimate.map = FALSE)

mycross = jittermap(mycross)
mycross = calc.genoprob(mycross, step=1, error.prob=0.001)

names(mycross$pheno)

permutations = scanone(mycross,pheno.col=2:11, method="hk", n.perm=1000,n.cluster=3)
summary(permutations, alpha=0.05)


####### model creation - exp1-4dpi

step1 = stepwiseqtl(mycross,pheno.col=2,method="hk",additive.only = T,penalties = 3.24)

model = makeqtl(mycross, chr=c(5,8,8), pos=c(2.9007,48.2423,54.7672),what="prob")
fqtl = fitqtl(mycross, pheno.col=2, qtl=model, formula=y~Q1+Q2+Q3, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step1,qtl.index=1,expandtomarkers=T)
bayesint(step1,qtl.index=2,expandtomarkers=T)
bayesint(step1,qtl.index=3,expandtomarkers=T)


####### model creation - exp1-5dpi

step2 = stepwiseqtl(mycross,pheno.col=3,method="hk",additive.only = T,penalties = 3.44)

model = makeqtl(mycross, chr=c(5), pos=c(3),what="prob")
fqtl = fitqtl(mycross, pheno.col=3, qtl=model, formula=y~Q1, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step2,qtl.index=1,expandtomarkers=T)


####### model creation - exp1-6dpi

step3 = stepwiseqtl(mycross,pheno.col=4,method="hk",additive.only = T,penalties = 3.41)

model = makeqtl(mycross, chr=c(5), pos=c(4.2434),what="prob")
fqtl = fitqtl(mycross, pheno.col=4, qtl=model, formula=y~Q1, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step3,qtl.index=1,expandtomarkers=T)


####### model creation - exp1-7dpi

step4 = stepwiseqtl(mycross,pheno.col=5,method="hk",additive.only = T,penalties = 3.55)

model = makeqtl(mycross, chr=c(5), pos=c(4.2434),what="prob")
fqtl = fitqtl(mycross, pheno.col=5, qtl=model, formula=y~Q1, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step4,qtl.index=1,expandtomarkers=T)

####### model creation - exp2-4dpi

step5 = stepwiseqtl(mycross,pheno.col=6,method="hk",additive.only = T,penalties = 3.24)

model = makeqtl(mycross, chr=c(5), pos=c(3),what="prob")
fqtl = fitqtl(mycross, pheno.col=6, qtl=model, formula=y~Q1, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step5,qtl.index=1,expandtomarkers=T)

####### model creation - exp2-5dpi

step6 = stepwiseqtl(mycross,pheno.col=7,method="hk",additive.only = T,penalties = 3.39)

model = makeqtl(mycross, chr=c(5), pos=c(3),what="prob")
fqtl = fitqtl(mycross, pheno.col=7, qtl=model, formula=y~Q1, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step6,qtl.index=1,expandtomarkers=T)

####### model creation - exp2-6dpi

step7 = stepwiseqtl(mycross,pheno.col=8,method="hk",additive.only = T,penalties = 3.39)

model = makeqtl(mycross, chr=c(5,7), pos=c(3,5.343),what="prob")
fqtl = fitqtl(mycross, pheno.col=8, qtl=model, formula=y~Q1+Q2, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step7,qtl.index=1,expandtomarkers=T)
bayesint(step7,qtl.index=2,expandtomarkers=T)

####### model creation - exp2-7dpi

step8 = stepwiseqtl(mycross,pheno.col=9,method="hk",additive.only = T,penalties = 3.43)

model = makeqtl(mycross, chr=c(5,7,8,35), pos=c(3,4,75.408,19),what="prob")
fqtl = fitqtl(mycross, pheno.col=9, qtl=model, formula=y~Q1+Q2+Q3+Q4, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step8,qtl.index=1,expandtomarkers=T)
bayesint(step8,qtl.index=2,expandtomarkers=T)
bayesint(step8,qtl.index=3,expandtomarkers=T)
bayesint(step8,qtl.index=4,expandtomarkers=T)

####### model creation - exp1-2dpiHR

step9 = stepwiseqtl(mycross,pheno.col=10,method="hk",additive.only = T,penalties = 3.46)

model = makeqtl(mycross, chr=c(5,6), pos=c(4.2434,2.5067),what="prob")
fqtl = fitqtl(mycross, pheno.col=10, qtl=model, formula=y~Q1+Q2, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step9,qtl.index=1,expandtomarkers=T)
bayesint(step9,qtl.index=2,expandtomarkers=T)

####### model creation - exp2-2dpiHR

step10 = stepwiseqtl(mycross,pheno.col=11,method="hk",additive.only = T,penalties = 3.47)

model = makeqtl(mycross, chr=c(5,8), pos=c(4.2434,73.3909),what="prob")
fqtl = fitqtl(mycross, pheno.col=11, qtl=model, formula=y~Q1+Q2, method="hk",get.ests = T)
summary(fqtl, pvalues=FALSE)

bayesint(step10,qtl.index=1,expandtomarkers=T)
bayesint(step10,qtl.index=2,expandtomarkers=T)




