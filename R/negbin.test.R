#!/bin/env
R --vanilla << "EOF" # Pipe all subsequent lines into R.

library(MASS)
library(lattice)
library(multcomp)

set.seed(123)

## Set a fake experiment dataset
features = c("feat1","feat2","feat3","feat4","feat5","feat6","feat7","feat8","feat9","feat10","feat11","feat12")
group1 = c(rnegbin(3, 200,2),rnegbin(3,100,1),rnegbin(3,700,2),rnegbin(3,400,3))
group2 = c(rnegbin(3, 200,0.2),rnegbin(3,100,0.11),rnegbin(3,300,0.7),rnegbin(3,400,0.4))
group3 = c(rnegbin(3, 700,2),rnegbin(3,400,1),rnegbin(3,200,2),rnegbin(3,100,1))
data = data.frame(g1=group1,g2=group2,g3=group3, row.names=features)

temp <- as.factor(rep(rep(c(30,30,30,40,40,40),each=2),3))
treatment <- rep(rep(rep(c("yes", "no"), each=3), 2),3)
feature <- rep(c("feat1", "feat2", "feat3"), each=12)
dtf <- data.frame(feature=feature, exprs=exprs,temp=temp, treatment=treatment)


## Read real dataset
design = "design.txt"
contrast = "contrast.txt"
design <- read.table(design, header=T, row.names=1,  sep=",")

i = 1
Y = as.vector(data[i,])
X = as.matrix(design) # Y = t(X).b
g = data.frame(Y=t(Y),X=t(X))

nbmodel = glm.nb(feat1 ~ X.temp + X.treat, data=g)

for (i in length(data[,1])){ # for each feature
    Y = as.vector(data[i,])
    X = as.matrix(design)
    nbmodel = glm.nb(Y ~ t(X))
}


main <- function(data.file, design, contrast){

  featureCount <- read.table(data.file, header=T, row.names=1, sep="\t")
  design <- read.table(design, header=T, row.names=1,  sep=",")

  contrast <- read.table(contrast.file, header=F, sep="\t")
  contrast.matrix <- data.matrix(contrast)
  model <- compareFactors(featureCount, design, contrast)

  result = data.frame()
  write.table(result, "ouput")

}


compareFactors <- function(featureCount, design, contrast){

  ## check the design matrix

  ## compute betas
  nbmodel <- glm.nb(exprs ~ treatment + temp, data=featureCount)

  ## check the contrast matrix

  ## compute contrasts
  C = glht(nbmodel, linfct=mcp(temp="Tukey"))

  #beta = nbmodel$coef
  #comparisons = transpose(contrast.matrix) %*% beta

  print(featureCount)
  print(summary(nbmodel))
  print(C)

}



print(bwplot(exprs ~ treatment | temp + feature, dtf))

by(dtf, feature, compareFactors, 1 , 1) #Apply a Function to a Data Frame Split by Factors




