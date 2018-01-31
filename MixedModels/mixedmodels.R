library(EMMREML)

#RH sporulation (computer vision)

pheno = read.csv("RHallsporCVpheno.csv")
pheno = pheno[complete.cases(pheno),]
pheno$Genotype = as.factor(pheno$Genotype)
pheno$Year = as.factor(pheno$Year)
pheno$Experiment = as.factor(pheno$Experiment)

#Create the design matrix for fixed effects
X = model.matrix(~ Year:Experiment+Dpi:Year:Experiment, data = pheno)
#Visualize the design matrix (optional)
#heatmap(X, Rowv = NA, Colv = NA, scale="column")
#X needs to be full rank, i.e., the columns of X cannot be linear combinations of each other
X = X[,-c(2,3,9)]
#Create design matricies for genetic and genetic x year effects
Zg = model.matrix(~ Genotype-1, data = pheno)
Zgy = model.matrix(~ Genotype:Year-1, data = pheno)

E = diag(length(unique(pheno$Year)))
I = diag(ncol(Zg))
#The Kronecker product K of E and I gives the covariance matrix for the genetic x year effects.
K = E %x% I

#K assumes no missing data for genotypes across years, which was not the case for our experiments due to cold injury of vines in 2016.
missingdata = apply(Zgy, 2, function(x) all(x == 0))
K = K[!missingdata, !missingdata]
Zgy = Zgy[,!missingdata]

#Fit the linear mixed model
fit = emmremlMultiKernel(y = pheno$Phenotype, X = X, Zlist = list(Zg,Zgy), Klist = list(I,K))
#Obtain the variance components
fitvarcomp = fit$weights * fit$Vu
#Calculate heritability
heritability = fitvarcomp[1]/(fitvarcomp[1]+(fitvarcomp[2]/2)+(fit$Ve/5))
#Obtain the estimated genetic effects, or breeding values (BVs).
breedingvalues = as.matrix(fit$uhat)
row.names(breedingvalues)[1:ncol(Zg)] = colnames(Zg)

#RH sporulation (manual)

pheno = read.csv("RHallspormanpheno.csv")
pheno = pheno[complete.cases(pheno),]
pheno$Genotype = as.factor(pheno$Genotype)
pheno$Year = as.factor(pheno$Year)
pheno$Experiment = as.factor(pheno$Experiment)


X = model.matrix(~ Year:Experiment+Dpi:Year:Experiment, data = pheno)
#heatmap(X, Rowv = NA, Colv = NA, scale="column")
X = X[,-c(2,3,9)]

Zg = model.matrix(~ Genotype-1, data = pheno)
Zgy = model.matrix(~ Genotype:Year-1, data = pheno)

E = diag(length(unique(pheno$Year)))
I = diag(ncol(Zg))
K = E %x% I

missingdata = apply(Zgy, 2, function(x) all(x == 0))
K = K[!missingdata, !missingdata]
Zgy = Zgy[,!missingdata]

fit = emmremlMultiKernel(y = pheno$Phenotype, X = X, Zlist = list(Zg,Zgy), Klist = list(I,K))

fitvarcomp = fit$weights * fit$Vu

heritability = fitvarcomp[1]/(fitvarcomp[1]+(fitvarcomp[2]/2)+(fit$Ve/5))

breedingvalues = as.matrix(fit$uhat)
row.names(breedingvalues)[1:ncol(Zg)] = colnames(Zg)


#RH HR

pheno = read.csv("RHallHRpheno.csv")
pheno = pheno[complete.cases(pheno),]
pheno$Genotype = as.factor(pheno$Genotype)
pheno$Year = as.factor(pheno$Year)
pheno$Experiment = as.factor(pheno$Experiment)


X = model.matrix(~ Year:Experiment, data = pheno)
#heatmap(X, Rowv = NA, Colv = NA, scale="column")
X = X[,-c(2,3)]

Zg = model.matrix(~ Genotype-1, data = pheno)
Zgy = model.matrix(~ Genotype:Year-1, data = pheno)

E = diag(length(unique(pheno$Year)))
I = diag(ncol(Zg))
K = E %x% I

missingdata = apply(Zgy, 2, function(x) all(x == 0))
K = K[!missingdata, !missingdata]
Zgy = Zgy[,!missingdata]

fit = emmremlMultiKernel(y = pheno$Phenotype, X = X, Zlist = list(Zg,Zgy), Klist = list(I,K))

fitvarcomp = fit$weights * fit$Vu

heritability = fitvarcomp[1]/(fitvarcomp[1]+(fitvarcomp[2]/2)+(fit$Ve/5))

breedingvalues = as.matrix(fit$uhat)
row.names(breedingvalues)[1:ncol(Zg)] = colnames(Zg)


#HC sporulation (computer vision)

pheno = read.csv("HCallsporCVpheno.csv")
pheno = pheno[complete.cases(pheno),]
pheno$Genotype = as.factor(pheno$Genotype)
pheno$Year = as.factor(pheno$Year)
pheno$Experiment = as.factor(pheno$Experiment)


X = model.matrix(~ Year:Experiment+Dpi:Year:Experiment+Trichomes:Year:Experiment, data = pheno)
#heatmap(X, Rowv = NA, Colv = NA, scale="column")
X = X[,-c(2,3,6,9,12,15,18)]

Zg = model.matrix(~ Genotype-1, data = pheno)
Zgy = model.matrix(~ Genotype:Year-1, data = pheno)

E = diag(length(unique(pheno$Year)))
I = diag(ncol(Zg))
K = E %x% I

missingdata = apply(Zgy, 2, function(x) all(x == 0))
K = K[!missingdata, !missingdata]
Zgy = Zgy[,!missingdata]

fit = emmremlMultiKernel(y = pheno$Phenotype, X = X, Zlist = list(Zg,Zgy), Klist = list(I,K))

fitvarcomp = fit$weights * fit$Vu

heritability = fitvarcomp[1]/(fitvarcomp[1]+(fitvarcomp[2]/2)+(fit$Ve/4))

breedingvalues = as.matrix(fit$uhat)
row.names(breedingvalues)[1:ncol(Zg)] = colnames(Zg)

#HC sporulation (manual)

pheno = read.csv("HCallspormanpheno.csv")
pheno = pheno[complete.cases(pheno),]
pheno$Genotype = as.factor(pheno$Genotype)
pheno$Year = as.factor(pheno$Year)
pheno$Experiment = as.factor(pheno$Experiment)


X = model.matrix(~ Year:Experiment+Dpi:Year:Experiment, data = pheno)
#heatmap(X, Rowv = NA, Colv = NA, scale="column")
X = X[,-c(2,3,6,9,12)]

Zg = model.matrix(~ Genotype-1, data = pheno)
Zgy = model.matrix(~ Genotype:Year-1, data = pheno)

E = diag(length(unique(pheno$Year)))
I = diag(ncol(Zg))
K = E %x% I

missingdata = apply(Zgy, 2, function(x) all(x == 0))
K = K[!missingdata, !missingdata]
Zgy = Zgy[,!missingdata]

fit = emmremlMultiKernel(y = pheno$Phenotype, X = X, Zlist = list(Zg,Zgy), Klist = list(I,K))

fitvarcomp = fit$weights * fit$Vu

heritability = fitvarcomp[1]/(fitvarcomp[1]+(fitvarcomp[2]/2)+(fit$Ve/4))

breedingvalues = as.matrix(fit$uhat)
row.names(breedingvalues)[1:ncol(Zg)] = colnames(Zg)

#HC HR

pheno = read.csv("HCallHRpheno.csv")
pheno = pheno[complete.cases(pheno),]
pheno$Genotype = as.factor(pheno$Genotype)
pheno$Year = as.factor(pheno$Year)
pheno$Experiment = as.factor(pheno$Experiment)


X = model.matrix(~ Year:Experiment, data = pheno)
#heatmap(X, Rowv = NA, Colv = NA, scale="column")
X = X[,-c(2,4)]

Zg = model.matrix(~ Genotype-1, data = pheno)
Zgy = model.matrix(~ Genotype:Year-1, data = pheno)

E = diag(length(unique(pheno$Year)))
I = diag(ncol(Zg))
K = E %x% I

missingdata = apply(Zgy, 2, function(x) all(x == 0))
K = K[!missingdata, !missingdata]
Zgy = Zgy[,!missingdata]

fit = emmremlMultiKernel(y = pheno$Phenotype, X = X, Zlist = list(Zg,Zgy), Klist = list(I,K))

fitvarcomp = fit$weights * fit$Vu

heritability = fitvarcomp[1]/(fitvarcomp[1]+(fitvarcomp[2]/2)+(fit$Ve/3))

breedingvalues = as.matrix(fit$uhat)
row.names(breedingvalues)[1:ncol(Zg)] = colnames(Zg)

#HC trichomes

pheno = read.csv("HCalltrichomespheno.csv")
pheno = pheno[complete.cases(pheno),]
pheno$Genotype = as.factor(pheno$Genotype)
pheno$Year = as.factor(pheno$Year)
pheno$Experiment = as.factor(pheno$Experiment)


X = model.matrix(~ Year:Experiment, data = pheno)
#heatmap(X, Rowv = NA, Colv = NA, scale="column")
X = X[,-c(2,3,6)]

Zg = model.matrix(~ Genotype-1, data = pheno)
Zgy = model.matrix(~ Genotype:Year-1, data = pheno)

E = diag(length(unique(pheno$Year)))
I = diag(ncol(Zg))
K = E %x% I

missingdata = apply(Zgy, 2, function(x) all(x == 0))
K = K[!missingdata, !missingdata]
Zgy = Zgy[,!missingdata]

fit = emmremlMultiKernel(y = pheno$Phenotype, X = X, Zlist = list(Zg,Zgy), Klist = list(I,K))

fitvarcomp = fit$weights * fit$Vu

heritability = fitvarcomp[1]/(fitvarcomp[1]+(fitvarcomp[2]/2)+(fit$Ve/4))

breedingvalues = as.matrix(fit$uhat)
row.names(breedingvalues)[1:ncol(Zg)] = colnames(Zg)



