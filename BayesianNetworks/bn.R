#Part of the code here comes from http://www.bnlearn.com/research/genetics14/

library(bnlearn)

#The phenotypes are breeding values (BVs) from the linear mixed models.
phenogeno = read.csv("RHallphenoandmap.csv")
phenogeno = phenogeno[-c(1:2),]
phenogeno[phenogeno=="-"] <- NA
phenogeno = phenogeno[complete.cases(phenogeno),]
phenogeno = apply(phenogeno,2,as.numeric)

pheno = phenogeno[,c(2:4)]
geno = phenogeno[,-c(1:4)]


#Obtain 1000 networks using a random set of 90% of the data. One can also load existing data below.
RHnets = c()

for(i in 1:1000){
  #Take a random sample of 90% of the data
  train <- sample(1:164,148)
  
  spormanBVtrain = pheno[train,2]
  HRBVtrain = pheno[train,3]
  genotrain = geno[train,]
  
  datatrain = as.data.frame(cbind(spormanBVtrain,HRBVtrain,genotrain))
  colnames(datatrain)[1:2] = c("spormanBV","HRBV")

  #Find Markov blankets for each trait
  spormanMB = learn.nbr(datatrain, node = "spormanBV", debug = FALSE, method = "si.hiton.pc", test = "cor", alpha = 0.01)
  HRMB = learn.nbr(datatrain, node = "HRBV", debug = FALSE, method = "si.hiton.pc", test = "cor", alpha = 0.01)
  
  nodes = unique(c("spormanBV","HRBV",spormanMB,HRMB))
  #Restrict manual sporulation from affecting HR and either trait from affecting SNPs
  blacklist = tiers2blacklist(list(nodes[-c(1:2)], "HRBV","spormanBV"))
  #Learn the network structure
  bn = hc(datatrain[,nodes], blacklist = blacklist)
  #Save networks
  RHnets[[i]] = bn
}

#Save networks (optional)
#saveRDS(RHnets, "RHnets.rds")
#Load networks (optional)
#RHnets <- readRDS("RHnets.rds")

#Arcs found in all the models.
arclist = list()
for(i in 1:1000){
  arclist[[i]] = arcs(RHnets[[i]])
}
#Obtain unique nodes from all 1000 networks
nodes = unique(unlist(arclist))
#Calculate arc strength of all arcs
strength = custom.strength(arclist, nodes = nodes)
#Average the networks
averaged = averaged.network(strength,threshold = 0.5)
#Keep only arcs above the strength threshold of 0.5
relevant.nodes = nodes(averaged)[sapply(nodes, degree, object = averaged) > 0]
averaged2 = subgraph(averaged, relevant.nodes)
strength2 = strength[(strength$from %in% relevant.nodes) & (strength$to %in% relevant.nodes), ]
#Visualize the averaged network
averagedgraph = strength.plot(averaged2, strength2, shape = "rectangle", layout = "fdp")


dataall = as.data.frame(cbind(pheno,geno))
colnames(dataall)[2:3] = c("spormanBV","HRBV")
#Obtain individual node distributions. They are linear models (compare the "fitted" object and the linear models below; both also have the effect size information).
fitted = bn.fit(averaged2, dataall[,nodes(averaged2)])

#Calculate the percent variance explained by (Type III SS/Total SS)*100

#Manual sporulation percent variance explained by S14_29543081
lmman = lm(spormanBV~S14_29543081,data=dataall[,nodes(averaged2)])
anova(lmman)
(1.2644/(1.2644+10.3587))*100 #10.88

#HR percent variance explained by S8_11656031
lmHR = lm(HRBV~S11_15397463+S16_22124674+S8_11656031,data=dataall[,nodes(averaged2)])
anova(lmHR)
(0.8588/(1.0238+0.4045+0.8588+4.4264))*100 #12.79
#HR percent variance explained by S11_15397463
lmHR = lm(HRBV~S16_22124674+S8_11656031+S11_15397463,data=dataall[,nodes(averaged2)])
anova(lmHR)
(0.8625/(0.6030+0.8215+0.8625+4.4264))*100 #12.85
#HR percent variance explained by S16_22124674
lmHR = lm(HRBV~S11_15397463+S8_11656031+S16_22124674,data=dataall[,nodes(averaged2)])
anova(lmHR)
(0.3618/(1.0238+0.9014+0.3618+4.4264))*100 #5.39

#Find arcs that appear in at least 5% of the networks. These will be used to find confidence intervals.
averagedconf = averaged.network(strength,threshold = 0.05)
averagedconf$arcs[,1]

#Confidence intervals
#Note that SNPs are divided between two parental maps. Each SNP below should be checked to see on what parental map it is on.
averagedconf$arcs[grepl("S14",averagedconf$arcs[,1]),][which(averagedconf$arcs[grepl("S14",averagedconf$arcs[,1]),][,2]=="spormanBV"),]
averagedconf$arcs[grepl("S8",averagedconf$arcs[,1]),][which(averagedconf$arcs[grepl("S8",averagedconf$arcs[,1]),][,2]=="HRBV"),]
averagedconf$arcs[grepl("S11",averagedconf$arcs[,1]),][which(averagedconf$arcs[grepl("S11",averagedconf$arcs[,1]),][,2]=="HRBV"),]
averagedconf$arcs[grepl("S16",averagedconf$arcs[,1]),][which(averagedconf$arcs[grepl("S16",averagedconf$arcs[,1]),][,2]=="HRBV"),]


#HC family

phenogeno = read.csv("HCallphenoandmap.csv")
phenogeno = phenogeno[-c(1:2),]
phenogeno[phenogeno=="-"] <- NA
phenogeno = phenogeno[complete.cases(phenogeno),]
phenogeno = apply(phenogeno,2,as.numeric)

pheno = phenogeno[,c(3:5)]
geno = phenogeno[,-c(1:5)]


HCnets = c()

for(i in 1:1000){
  train <- sample(1:152,137)
  
  spormanBVtrain = pheno[train,1]
  HRBVtrain = pheno[train,2]
  trichomeBVtrain = pheno[train,3]
  genotrain = geno[train,]
  
  datatrain = as.data.frame(cbind(spormanBVtrain,HRBVtrain,trichomeBVtrain,genotrain))
  colnames(datatrain)[1:3] = c("spormanBV","HRBV","trichomeBV")
  
  spormanMB = learn.nbr(datatrain, node = "spormanBV", debug = FALSE, method = "si.hiton.pc", test = "cor", alpha = 0.01)
  HRMB = learn.nbr(datatrain, node = "HRBV", debug = FALSE, method = "si.hiton.pc", test = "cor", alpha = 0.01)
  trichomeMB = learn.nbr(datatrain, node = "trichomeBV", debug = FALSE, method = "si.hiton.pc", test = "cor", alpha = 0.01)
  
  nodes = unique(c("spormanBV","HRBV","trichomeBV",spormanMB,HRMB,trichomeMB))
  blacklist = tiers2blacklist(list(nodes[-c(1:3)],"trichomeBV","HRBV","spormanBV"))
  bn = hc(datatrain[,nodes], blacklist = blacklist)
  HCnets[[i]] = bn
}

#saveRDS(HCnets, "HCnets.rds")
#HCnets <- readRDS("HCnets.rds")

arclist = list()
for(i in 1:1000){
  arclist[[i]] = arcs(HCnets[[i]])
}
nodes = unique(unlist(arclist))
strength = custom.strength(arclist, nodes = nodes)
averaged = averaged.network(strength,threshold = 0.5)
relevant.nodes = nodes(averaged)[sapply(nodes, degree, object = averaged) > 0]
averaged2 = subgraph(averaged, relevant.nodes)
strength2 = strength[(strength$from %in% relevant.nodes) & (strength$to %in% relevant.nodes), ]
averagedgraph = strength.plot(averaged2, strength2, shape = "rectangle", layout = "fdp")


dataall = as.data.frame(cbind(pheno,geno))
colnames(dataall)[1:3] = c("spormanBV","HRBV","trichomeBV")
fitted = bn.fit(averaged2, dataall[,nodes(averaged2)])

#Manual sporulation percent variance explained by trichomeBV
lmman = lm(spormanBV~HRBV+S7_2610743+trichomeBV,data=dataall[,nodes(averaged2)])
anova(lmman)
(3.4577/(10.5343+4.2885+3.4577+14.5401))*100 #10.54
#Manual sporulation percent variance explained by HRBV
lmman = lm(spormanBV~S7_2610743+trichomeBV+HRBV,data=dataall[,nodes(averaged2)])
anova(lmman)
(1.2276/(6.0958+10.9572+1.2276+14.5401))*100 #3.74
#Manual sporulation percent variance explained by S7_2610743
lmman = lm(spormanBV~HRBV+trichomeBV+S7_2610743,data=dataall[,nodes(averaged2)])
anova(lmman)
(2.5126/(10.5343+5.2336+2.5126+14.5401))*100 #7.66

#HR percent variance explained by trichomeBV
lmHR = lm(HRBV~S6_6640223+trichomeBV,data=dataall[,nodes(averaged2)])
anova(lmHR)
(2.4539/(1.2646+2.4539+4.1224))*100 #31.30
#HR percent variance explained by S6_6640223
lmHR = lm(HRBV~trichomeBV+S6_6640223,data=dataall[,nodes(averaged2)])
anova(lmHR)
(0.5070/(3.2115+0.5070+4.1224))*100 #6.47


#Trichome percent variance explained by S7_2610743
lmLt = lm(trichomeBV~S8_17766225+S15_17663308+S7_2610743,data=dataall[,nodes(averaged2)])
anova(lmLt)
(0.009123/(0.036193+0.013246+0.009123+0.108017))*100 #5.48
#Trichome percent variance explained by S8_17766225
lmLt = lm(trichomeBV~S15_17663308+S7_2610743+S8_17766225,data=dataall[,nodes(averaged2)])
anova(lmLt)
(0.027994/(0.017951+0.012616+0.027994+0.108017))*100 #16.81
#Trichome percent variance explained by S15_17663308
lmLt = lm(trichomeBV~S8_17766225+S7_2610743+S15_17663308,data=dataall[,nodes(averaged2)])
anova(lmLt)
(0.014577/(0.036193+0.007791+0.014577+0.108017))*100 #8.75

averagedconf = averaged.network(strength,threshold = 0.05)
averagedconf$arcs[,1]

averagedconf$arcs[grepl("S6",averagedconf$arcs[,1]),][which(averagedconf$arcs[grepl("S6",averagedconf$arcs[,1]),][,2]=="HRBV"),]
averagedconf$arcs[grepl("S7",averagedconf$arcs[,1]),][which(averagedconf$arcs[grepl("S7",averagedconf$arcs[,1]),][,2]=="spormanBV"),]
averagedconf$arcs[grepl("S7",averagedconf$arcs[,1]),][which(averagedconf$arcs[grepl("S7",averagedconf$arcs[,1]),][,2]=="trichomeBV"),]
averagedconf$arcs[grepl("S8",averagedconf$arcs[,1]),][which(averagedconf$arcs[grepl("S8",averagedconf$arcs[,1]),][,2]=="trichomeBV"),]
averagedconf$arcs[grepl("S15",averagedconf$arcs[,1]),]







