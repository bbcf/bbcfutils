library(MASS)
library(lattice)
library(multcomp)

set.seed(123)

## Set a fake experiment dataset
exprs <- c(rnegbin(3, 200,2),rnegbin(3,100,1),rnegbin(3,700,2),rnegbin(3,400,3),
           rnegbin(3, 200,0.2),rnegbin(3,100,0.11),rnegbin(3,300,0.7),rnegbin(3,400,0.4),
           rnegbin(3, 700,2),rnegbin(3,400,1),rnegbin(3,200,2),rnegbin(3,100,1)
           )
temp <- as.factor(rep(rep(c(30,30,30,40,40,40),each=2),3))
treatment <- rep(rep(rep(c("yes", "no"), each=3), 2),3)
feature = rep(c("feat1", "feat2", "feat3"), each=12)
dtf <- data.frame(feature=feature, exprs=exprs,temp=temp, treatment=treatment)







compareFactors <- function(featureCount, design, contrast){

  ## check the design matrix
  
  ## compute betas
  dtf.nbmodel <- glm.nb(exprs ~ treatment + temp, data=featureCount)
  
  ## check the contrast matrix
  
  ## compute contrasts
  C = glht(dtf.nbmodel, linfct=mcp(temp="Tukey"))
  
  print(featureCount)
  print(summary(dtf.nbmodel))
  print(C)
}



print(bwplot(exprs ~ treatment | temp + feature, dtf))

by(dtf, feature, compareFactors, 1 , 1)




