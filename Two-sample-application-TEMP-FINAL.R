require(nlme)
require(multcomp)
require(lattice)
load("Data.RData")
setwd("/YY Go Data")  
source("/functions.R")
source("/functions-simulation-SG2-v2.R")

# SET UP
  Data.YYG <- read.delim("TEMP-YYG.txt")
  getRows_A <- Data.YYG[which(!is.na(Data.YYG$Temp)),]
  getCols_A <- data.frame(Group='A', DPI=getRows_A$DPI, Subject=getRows_A$Stallion, Response=getRows_A$Temp, Group='A')
  getCols_A$Subject <- droplevels(getCols_A$Subject)
  Subjects_A <- as.vector(unique(getCols_A$Subject))

  getRows_B <- Data[which(!is.na(Data$Temp) & Data$DPI<=14),] 
  getCols_B <- data.frame(Group='B', DPI=getRows_B$DPI, Subject=getRows_B$Stallion, Response=getRows_B$Temp, Group='B')
  getCols_B$Subject <- droplevels(getCols_B$Subject)
  Subjects_B <- as.vector(unique(getCols_B$Subject))

  horse <- rbind(getCols_A,getCols_B)
  horse.g <- groupedData(Response~DPI|Subject, data=horse)

# VISUALIZE THE DATA
  hist(horse$Response)
  with(horse, xyplot(Response~DPI|Subject, group=Group, type="b", pch=19, cex=0.6, xlab="Days post-infection", ylab="Temperature (F)"))
  with(horse[horse$Group=='A',], xyplot(Response~DPI|Subject, type="b", pch=19, xlab="Days post-infection", ylab="Temperature (F)",ylim=c(97,106),layout=c(4,2)))
  with(horse[horse$Group=='B',], xyplot(Response~DPI|Subject, type="b", pch=19, xlab="Days post-infection", ylab="Temperature (F)",ylim=c(97,106),layout=c(4,2)))
  # ILLUSTRATION FOR POWERPOINT
    with(horse[horse$Group=='B' & horse$Subject=='L138',], plot(DPI,Response,type="o",pch=19,col="blue",xaxt="n")) 
    abline(h=100.8,col="gray",lty=2); abline(h=103.1,col="gray",lty=2)
    axis(side=1, at=seq(-7,15))

# ASSESS WHICH PARAMETERS DON'T NEED RANDOM EFFECTS
  fit.list <- nlsList(Response ~ xGauss.log(DPI,mu_,Z,B,I)|Group/Subject,
                      data=horse.g,
                      start = c(B=99, I = 4, mu_ = 5, Z=1),
                      control=nls.control(maxiter=500))
  plot(intervals(fit.list),layout=c(4,1))

# FIT MODEL
  fit.full <- nlme(Response ~ xGauss.log(DPI,mu_,Z,B,I), method="ML",
              data=horse.g,
              #groups=~Subject,
              fixed=B+I+mu_+Z~Group,
              random=pdDiag(B+I+mu_~1), 
              #weight=varPower(form=~fitted(.)|Group),
              #correlation=corAR1(),
              #random=B~1|Subject,   # For non-grouped data
              start = c(B=c(99,0), I=c(4,0), mu_=c(5,0), Z=c(1,0)), 
              control=nlmeControl(maxIter=500))
  summary(fit.full)

  fit.reduced <- nlme(Response ~ xGauss.log(DPI,mu_,Z,B,I), method="ML",
                 data=horse,
                 groups=~Subject,
                 fixed=B+I+mu_+Z~1,
                 random=pdDiag(B+I+mu_~1), 
                 start = c(B=c(B=99, I = 4, mu_ = 5, Z=1)), 
                 control=nlmeControl(maxIter=500))
  summary(fit.reduced)

# DIAGNOSTICS
  diagnose.G(fit.full)

# PLOTS
  plot.G2(fit.full, "DPI","Temperature (F)",c(0,42),c(0,1),c(2,4))
  plot.G2.compare(fit.full, "DPI","Temperature (F)",c(-7,14),c(99,104),c(2,4),means=FALSE,legend.names=c('Bucyrus strain in mares','KY-84 strain in stallions'))
  plot(augPred(fit.full), cex=0.6, xlab="Days post-infection", ylab="Temperature (F)")

# TESTS (BE SURE TO CHECK THAT LOG.S PARAMETER IS CORRECTLY SPECIFIED)
  test.G2(fit.reduced, fit.full, log.s=TRUE)

# DURATION
  bootstrap.duration.G2(fit.full, 100, alpha=0.05)

## RM-ANOVA
  #fit.2wayANOVA.lme(horse)
  Data.Combined$DPI.f <- factor(Data.Combined$DPI)
  #fit.lme <- lme(Response~Group*DPI.f, random=~1|Subject, data=horse)
  #confint(glht(fit.lme, linfct=mcp(DPI.f="Tukey")))
  fit.aov <- aov(Response~Group*DPI.f + Error(Subject), Data.Combined)
  summary(fit.aov)

