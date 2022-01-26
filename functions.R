##################
## GIRMM models ##
##################

# Function that returns the value of an extended Gaussian with given parameters and time value (DPI)  

xGauss <- function(DPI,mu,sigma,B,I) {
  return(B+I*exp(-(DPI-mu)^2/(2*sigma^2)))
}

xGauss.log <- function(DPI,mu,sigma,B,I) {
  # exp(sigma) in place of sigma for estimating Z=log(sigma) instead of sigma.
  # This amounts to the same as putting a lognormal distribution on the sigma_j
  # instead of normal, which assures that the sigma_j will be positive.
    return(B+I*exp(-(DPI-mu)^2/(2*exp(sigma)^2)))
}

xGauss2 <- function(DPI,mu,sigma1,sigma2,B,I) {
  return(B+I*exp(-(DPI-mu)^2/(2*ifelse(DPI<=as.numeric(mu),sigma1,sigma2)^2)))
}

xGauss2.log <- function(DPI,mu,sigma1,sigma2,B,I) {
  return(B+I*exp(-(DPI-mu)^2/(2*ifelse(DPI<=as.numeric(mu),exp(sigma1),exp(sigma2))^2)))
}

xGauss2.test1 <- function(DPI,mu,sigma1,sigma2,B,I) {
  return(B+I*exp(-(DPI-mu)^2/(2*((DPI<=as.numeric(mu))*sigma1+(DPI>as.numeric(mu))*sigma2)^2)))
}

xGauss3 <- function(DPI,mu,sigma1,sigma2,B1,B2,I1) {
  v1 <- NULL; v2 <- NULL;
  # Split DPI vector into times less than mu and times greater than mu
    v1 <- DPI[DPI<=as.numeric(mu)]; v2 <- DPI[DPI>as.numeric(mu)]
  # Compute functional values for each partial vector
    f1 <- B1+I1*exp(-(v1-mu)^2/(2*sigma1^2)); f2 <- B2+(B1-B2+I1)*exp(-(v2-mu)^2/(2*sigma2^2))
  # Concatenate partial vectors of functional values and return
    return(c(f1,f2))
}

xGauss4 <- function(DPI,mu,sigma1,sigma2,B,I,onset,recov) {
  onset <- as.numeric(onset); recov <- as.numeric(recov)
  pre.t <- DPI[DPI<onset]; post.t <- DPI[DPI>recov]
  rsp.t <- DPI[DPI >= onset & DPI <= recov]
  pre.f <- pre.t*0; post.f <- post.t*0
  rsp.f <- B+I*exp(-(rsp.t-mu)^2/(2*ifelse(rsp.t<=as.numeric(mu),sigma1,sigma2)^2))
  return(c(pre.f,rsp.f,post.f))
}

xGamma <- function(DPI,alpha,beta,B,I) {
  return(B+I*DPI^(alpha-1)*exp(-beta*DPI))
}

xTricube <- function(u,B,I,t.onset,t.rec) {
  return(B+I*(1-abs(1 + 2*(t.rec-u)/(t.onset-t.rec))^3)^3*(abs(1 + 2*(t.rec-u)/(t.onset-t.rec))<1))
}

xTricube2 <- function(t,p,s1,s2,B,I) {
	return(B + I*(1-abs((t-p)/s1)^3)^3*((p-t>=0) & (t-p+s1>=0)) + I*(1-abs((t-p)/s2)^3)*((t-p>=0) & (p-t+s2>=0)))
}

xTricube.fixed.B <- function(t, t.peak, I, lam.l, lam.r) {
  t.peak.n <- as.numeric(t.peak)
  u <- (t > t.peak.n)*(t-t.peak.n)/lam.r + (t <= t.peak.n)*(t-t.peak.n)/lam.l
  (u > -1) * (u < 1) * I * (1-abs(u^3))^3
}

xTricube.fixed.B.linear.onset <- function(t, t.peak, I, lam.r) {
  t.peak.n <- as.numeric(t.peak)
  t[t<0] <- 0
  u <- (t > t.peak.n)*(t-t.peak.n)/lam.r + (t <= t.peak.n)*0
  y <- (u > -1) * (u < 1) * ifelse(t>t.peak.n, I * (1-abs(u^3))^3, (I/t.peak)*t)
  return(y)
}

xTricube.fixed.B.Symmetric <- function(t, t.peak, I, lam) {
  t.peak.n <- as.numeric(t.peak)
  u <- (t-t.peak.n)/lam
  (u > -1) * (u < 1) * I * (1-abs(u^3))^3
}

# Epanechnikov kernel function
  
  # For a response variable with homogeneous variance, e.g. body temperature
    xEpa2 <- function(t,p,s1,s2,B,I) {
    return(B + I*(1-((t-p)/s1)^2)*((p-t>=0) & (t-p+s1>=0)) + I*(1-((t-p)/s2)^2)*((t-p>=0) & (p-t+s2>=0)))
  }

  # For a response variable with heterogenous variance, specifically a fixed baseline, e.g. viral load
    xEpa.fixed.B <- function(t, t.peak, I, lam.l, lam.r) {
      t.peak.n <- as.numeric(t.peak)
      u <- (t > t.peak.n)*(t-t.peak.n)/lam.r + (t <= t.peak.n)*(t-t.peak.n)/lam.l
      (u > -1) * (u < 1) * I * (1-u^2)
    }

  # Fixed baseline/viral load, symmetric
    xEpa.fixed.B.Symmetric <- function(t, t.peak, I, lam) {
      t.peak.n <- as.numeric(t.peak)
      u <- (t-t.peak.n)/lam
      (u > -1) * (u < 1) * I * (1-u^2)
    }  
    
xRectangular <- function(t,t_o,t_r,B,I) {
  return(B*(t_o-t>=0) + I*(t-t_o>0)*(t_r-t>0) + B*(t-t_r>=0))
}

########################################
# POST-FIT PROCESSING FOR GIRMM MODELS #
########################################

# Takes a fitted one-sample symmetric GIRMM model and returns a vector of estimated 
# means for individual subjects (level=1) or for the group (level=0). Predictions are 
# returned for time points starting with the earliest time point in the data set and at 
# increments according to "incr" parameter up through the latest time point in the data. 
# NOTE: IN THE DATA SET ON WHICH MODEL WAS FIT MUST, THE TIME VARIABLE MUST BE
# CALLED "DPI" AND THE GROUPING VARIABLE MUST BE CALLED "Subject".

  predict.G1 <- function(fm, level="subjects", incr=1, t.start=-999, t.end=-999) {
    
  # Set up
    time.var <- all.vars(getCall(fm))[2]
    data <- getData(fm)
    #start.time <- min(data[,time.var])
    #end.time <- max(data[,time.var])
    start.time <- ifelse(t.start==-999, min(data[,time.var]), t.start)
    end.time <- ifelse(t.end==-999, max(data[,time.var]), t.end)
    DPI.seq <- seq(start.time, end.time, by=incr)
    subjects <- as.character(unique(data$Subject))
      
  # Predictions for individual subjects (level=1)
    if (level=="subjects") {
      DPI.xx <- rep(DPI.seq, length(subjects))
      subjects.xx <- NULL
      for (i in subjects) subjects.xx <- c(subjects.xx, rep(i,length(DPI.seq)))
      df.xx <- data.frame(DPI=DPI.xx, Subject=subjects.xx)
      Pred <- predict(fm, newdata=df.xx)
      result <- cbind(df.xx, Pred)
    }
  
  # Predictions for group (level=0)
    if (level=="group") {
      DPI.xx <- DPI.seq
      df.xx <- data.frame(DPI=DPI.xx)
      Pred <- predict(fm, newdata=df.xx, level=0)
      result <- cbind(df.xx, Pred)
    }
    
    return(result)
    
  }

# NOTE: FOR PREDICT.G2, IN THE DATA SET ON WHICH MODEL WAS FIT MUST, 
# THE TIME VARIABLE MUST BE CALLED "DPI", THE VARIABLE FOR GROUPING OBSERVATIONS
# BY SUBJECT MUST BE CALLED "Subject", AND THE VARIABLE FOR CLASSIFYING SUBJECTS
# INTO THE TWO GROUPS MUST BE CALLED "Group".

  predict.G2 <- function(fm, level="subjects", incr=1) {
  
  # Set up
    time.var <- all.vars(getCall(fm))[2]
    data <- getData(fm)
    start.time <- min(data[,time.var])
    end.time <- max(data[,time.var])
    DPI.seq <- seq(start.time, end.time, by=incr)
  
  # Predictions for individual subjects (level=1)
    if (level=="subjects") {
      df <- unique(data.frame(Group=data$Group, Subject=data$Subject))
      df$Group <- as.character(df$Group); df$Subject <- as.character(df$Subject)    
      DPI.xx <- rep(DPI.seq, length(df$Subject))
      subjects.xx <- NULL; groups.xx <- NULL
      for (i in df$Subject) {
        subjects.xx <- c(subjects.xx, rep(i,length(DPI.seq)))
        groups.xx <- c(groups.xx, rep(df$Group[df$Subject==i],length(DPI.seq)))
      }
      df.xx <- data.frame(DPI=DPI.xx, Group=groups.xx, Subject=subjects.xx)
      Pred <- predict(fm, newdata=df.xx)
    }
  
  # Predictions for group (level=0)
     if (level=="group") {
       df <- unique(data.frame(Group=data$Group))
       df$Group <- as.character(df$Group)
       DPI.xx <- rep(DPI.seq, length(unique(df$Group)))
       groups.xx <- NULL
       for (i in df$Group) groups.xx <- c(groups.xx, rep(i,length(DPI.seq)))
       df.xx <- data.frame(DPI=DPI.xx, Group=groups.xx)
       Pred <- predict(fm, newdata=df.xx, level=0)
    }
    
    result <- cbind(df.xx, Pred)
    return(result)
  
  }

# DIAGNOSTICS FOR ONE-SAMPLE SYMMETRIC GIRMM MODEL

diagnose.G <- function(fm) {

  cat("\n\n***** SUMMARY (fm) *****\n\n")
  print(summary(fm))      # Anything strange in estimates of fixed effects/SE's, correlations, etc?
  cat("\n\n***** INTERVALS(fm) *****\n\n")
  print(intervals(fm))    # Anything strange in confidence intervals, esp. for variance component parameters?
  print(plot(fm))         # Anything strange in residual plots?
  print(plot(fm, form=resid(., type="p")~fitted(.)|Subject, abline=0))
   
}

# PLOTS FOR ONE-SAMPLE SYMMETRIC GIRMM MODEL

plot.G1 <- function(fm, xlab_, ylab_, xlim_, ylim_, panel, baseline=-999, sub.over=TRUE, sub.sep=TRUE, tstart_=-999, tend_=-999) {
  
  if (!"package:lattice" %in% search()) require(lattice)
  pred.sub <- predict.G1(fm, incr=0.1, t.start=tstart_, t.end=tend_)
  pred.grp <- predict.G1(fm, incr=0.1, level="group", t.start=tstart_, t.end=tend_)
  plot.data <- getData(fm)

# Data plus individual curves plus population curve - overlaid single plot
  plot(plot.data$DPI,plot.data$Response,pch=19,cex=0.5,col="gray",xlim=xlim_, ylim=ylim_, xlab=xlab_,ylab=ylab_)
  if (sub.over==TRUE) {
    for (i in unique(pred.sub$Subject)) {
    lines(pred.sub$DPI[pred.sub$Subject==i], pred.sub$Pred[pred.sub$Subject==i],col="gray")
    }
  }
  lines(pred.grp$DPI, pred.grp$Pred, lwd=2, col="black")

# Onset/recovery plot (same as previous but with onset/recoverY added)
  plot(plot.data$DPI,plot.data$Response,pch=19,cex=0.5, col="gray", xlim=xlim_, ylim=ylim_, xlab=xlab_,ylab=ylab_)
  if (sub.over==TRUE) {
    for (i in unique(pred.sub$Subject)) {
    lines(pred.sub$DPI[pred.sub$Subject==i], pred.sub$Pred[pred.sub$Subject==i],col="gray")
    }
  }
  lines(pred.grp$DPI, pred.grp$Pred, lwd=2, col="black")
  ORD <- duration.G1(fm, baseline)
  polygon(c(min(pred.grp$DPI),max(pred.grp$DPI),max(pred.grp$DPI),min(pred.grp$DPI)),c(ORD$B.lower, ORD$B.lower, ORD$B.upper, ORD$B.upper),col=rgb(0.977,0.922,0.828,0.5),border=NA)
  abline(v=ORD$onset,col="deepskyblue3")
  abline(v=ORD$recovery,col="deepskyblue3")

# Data plus individual curves - separated on single plot
  par(mfrow=c(2,4))
  for (i in unique(pred.sub$Subject)) {
    plot(plot.data$DPI[plot.data$Subject==i],plot.data$Response[plot.data$Subject==i],pch=19,cex=0.5,main=i,xlim=xlim_,ylim=ylim_,xlab=xlab_,ylab=ylab_)
    lines(pred.sub$DPI[pred.sub$Subject==i], pred.sub$Pred[pred.sub$Subject==i])
  }
  #with(pred.sub, xyplot(Pred~DPI|Subject, type="l", col="red"))
  
# Data plus individual curves - on seperate plots
  par(mfrow=c(1,1))
  if (sub.sep==TRUE) {
    for (i in unique(pred.sub$Subject)) {
      plot(plot.data$DPI[plot.data$Subject==i],plot.data$Response[plot.data$Subject==i],pch=19,cex=0.5,main=i,xlim=xlim_,ylim=ylim_,xlab=xlab_,ylab=ylab_)
      lines(pred.sub$DPI[pred.sub$Subject==i], pred.sub$Pred[pred.sub$Subject==i])
    }
  }

}

# PLOTS FOR TWO-SAMPLE SYMMETRIC GIRMM MODEL

plot.G2 <- function(fm, xlab_, ylab_, xlim_, ylim_, panel, baseline=-999) {
  
  pred.sub <- predict.G2(fm, incr=0.1)
  pred.grp <- predict.G2(fm, incr=0.1, level="group")
  plot.data <- getData(fm)

# Data plus individual curves plus population curve - overlaid single plot
  plot(plot.data$DPI[plot.data$Group=="A"],plot.data$Response[plot.data$Group=="A"],pch=19,cex=0.5,xlim=xlim_, ylim=ylim_,col="gray", main="Group A", xlab=xlab_,ylab=ylab_)
  for (i in unique(pred.sub$Subject[pred.sub$Group=="A"])) {
    lines(pred.sub$DPI[pred.sub$Subject==i], pred.sub$Pred[pred.sub$Subject==i],col="gray")
  }
  lines(pred.grp$DPI[pred.grp$Group=="A"], pred.grp$Pred[pred.grp$Group=="A"], lwd=2, col="black")
  
  plot(plot.data$DPI[plot.data$Group=="B"],plot.data$Response[plot.data$Group=="B"],pch=19,cex=0.5,xlim=xlim_, ylim=ylim_,col="gray", main="Group B", xlab=xlab_,ylab=ylab_)
  for (i in unique(pred.sub$Subject[pred.sub$Group=="B"])) {
    lines(pred.sub$DPI[pred.sub$Subject==i], pred.sub$Pred[pred.sub$Subject==i],col="gray")
  }
  lines(pred.grp$DPI[pred.grp$Group=="B"], pred.grp$Pred[pred.grp$Group=="B"], lwd=2, col="black")

# Onset/recovery plot (same as previous but with onset/recoverY added)
  ORD <- duration.G2(fm, baseline)

  plot(plot.data$DPI[plot.data$Group=="A"],plot.data$Response[plot.data$Group=="A"],pch=19,cex=0.5,xlim=xlim_, ylim=ylim_,col="gray", main="Group A", xlab=xlab_,ylab=ylab_)
  for (i in unique(pred.sub$Subject[pred.sub$Group=="A"])) {
    lines(pred.sub$DPI[pred.sub$Subject==i], pred.sub$Pred[pred.sub$Subject==i],col="gray")
  }
  lines(pred.grp$DPI[pred.grp$Group=="A"], pred.grp$Pred[pred.grp$Group=="A"], lwd=2, col="black")
  polygon(c(min(pred.grp$DPI),max(pred.grp$DPI),max(pred.grp$DPI),min(pred.grp$DPI)),c(ORD$B.lower.A, ORD$B.lower.A, ORD$B.upper.A, ORD$B.upper.A),col=rgb(0.977,0.922,0.828,0.5),border=NA)
  abline(v=ORD$onset.A,col="deepskyblue3")
  abline(v=ORD$recovery.A,col="deepskyblue3")
  
  plot(plot.data$DPI[plot.data$Group=="B"],plot.data$Response[plot.data$Group=="B"],pch=19,cex=0.5, xlim=xlim_, ylim=ylim_,col="gray", main="Group B", xlab=xlab_, ylab=ylab_)
  for (i in unique(pred.sub$Subject[pred.sub$Group=="B"])) {
    lines(pred.sub$DPI[pred.sub$Subject==i], pred.sub$Pred[pred.sub$Subject==i],col="gray")
  }
  lines(pred.grp$DPI[pred.grp$Group=="B"], pred.grp$Pred[pred.grp$Group=="B"], lwd=2, col="black")
  polygon(c(min(pred.grp$DPI),max(pred.grp$DPI),max(pred.grp$DPI),min(pred.grp$DPI)),c(ORD$B.lower.B, ORD$B.lower.B, ORD$B.upper.B, ORD$B.upper.B),col=rgb(0.977,0.922,0.828,0.5),border=NA)
  abline(v=ORD$onset.B,col="deepskyblue3")
  abline(v=ORD$recovery.B,col="deepskyblue3")

# Data plus individual curves - separated on single plot
  par(mfrow=panel)
  for (i in unique(pred.sub$Subject[pred.sub$Group=="A"])) {
    plot(plot.data$DPI[plot.data$Subject==i],plot.data$Response[plot.data$Subject==i],pch=19,cex=0.5, main=paste("A,",i), xlim=xlim_, ylim=ylim_ , xlab=xlab_ ,ylab=ylab_)
    lines(pred.sub$DPI[pred.sub$Subject==i], pred.sub$Pred[pred.sub$Subject==i])
  } 
  
  for (i in unique(pred.sub$Subject[pred.sub$Group=="B"])) {
    plot(plot.data$DPI[plot.data$Subject==i],plot.data$Response[plot.data$Subject==i],pch=19,cex=0.5, main=paste("B,",i), xlim=xlim_, ylim=ylim_ , xlab=xlab_ ,ylab=ylab_)
    lines(pred.sub$DPI[pred.sub$Subject==i], pred.sub$Pred[pred.sub$Subject==i])
  }

# Data plus individual curves - on seperate plots
  par(mfrow=c(1,1))
  for (i in unique(pred.sub$Subject[pred.sub$Group=="A"])) {
    plot(plot.data$DPI[plot.data$Subject==i],plot.data$Response[plot.data$Subject==i],pch=19,cex=0.5, main=paste("A",i), xlim=xlim_, ylim=ylim_, xlab=xlab_ , ylab=ylab_)
    lines(pred.sub$DPI[pred.sub$Subject==i], pred.sub$Pred[pred.sub$Subject==i])
  }
  
    for (i in unique(pred.sub$Subject[pred.sub$Group=="B"])) {
    plot(plot.data$DPI[plot.data$Subject==i],plot.data$Response[plot.data$Subject==i],pch=19,cex=0.5, main=paste("B",i), xlim=xlim_, ylim=ylim_, xlab=xlab_ , ylab=ylab_)
    lines(pred.sub$DPI[pred.sub$Subject==i], pred.sub$Pred[pred.sub$Subject==i])
  }

}

plot.G2.compare <- function(fm, xlab_, ylab_, xlim_, ylim_, panel, means=TRUE, legend.names, baseline=-999) {
  
  pred.sub <- predict.G2(fm, incr=0.1)
  pred.grp <- predict.G2(fm, incr=0.1, level="group")
  plot.data <- getData(fm)
  
  #plot(plot.data$DPI[plot.data$Group=="A"],plot.data$Response[plot.data$Group=="A"],pch=19,cex=0.5,xlim=xlim_, ylim=ylim_,col="gray", main="Group A", xlab=xlab_,ylab=ylab_)
  #for (i in unique(pred.sub$Subject[pred.sub$Group=="A"])) {
   # lines(pred.sub$DPI[pred.sub$Subject==i], pred.sub$Pred[pred.sub$Subject==i],col="gray")
  #}
  plot(pred.grp$DPI[pred.grp$Group=="A"], pred.grp$Pred[pred.grp$Group=="A"], type="l", lwd=1, lty=1, col="black", pch=19,cex=0.5,xlim=xlim_, ylim=ylim_, xlab=xlab_,ylab=ylab_, xaxt="n",yaxt="n")
  
  #plot(plot.data$DPI[plot.data$Group=="B"],plot.data$Response[plot.data$Group=="B"],pch=19,cex=0.5,xlim=xlim_, ylim=ylim_,col="gray", main="Group B", xlab=xlab_,ylab=ylab_)
  #for (i in unique(pred.sub$Subject[pred.sub$Group=="B"])) {
  #  lines(pred.sub$DPI[pred.sub$Subject==i], pred.sub$Pred[pred.sub$Subject==i],col="gray")
  #}
  lines(pred.grp$DPI[pred.grp$Group=="B"], pred.grp$Pred[pred.grp$Group=="B"], lwd=1, lty=1, col="darkgray")
  legend("topleft", legend.names, lty=c(1,1), col=c("blue","orange"))
    
  if (means==TRUE) {
    t <- unique(plot.data$DPI)
    mean.group.A <- rep(NA,length(t))
    for (i in 1:length(t)) mean.group.A[i] <- mean(plot.data$Response[plot.data$Group=="A" & plot.data$DPI==t[i]])
    mean.group.B <- rep(NA,length(t))
    for (i in 1:length(t)) mean.group.B[i] <- mean(plot.data$Response[plot.data$Group=="B" & plot.data$DPI==t[i]])
    lines(t,mean.group.A,type="l",col="blue"); points(t,mean.group.A,col="blue")
    lines(t,mean.group.B,type="l",col="orange"); points(t,mean.group.B,col="orange")
  }
  
}

# TESTS FOR ONE-SAMPLE SYMMETRIC GIRMM MODEL
# NOTE: Assumes lognormal distribution was specified for dispersion parameter (s)
# If normal distribution was specified instead, set log.s=FALSE

test.G1 <- function(fm.reduced, fm.full, log.s=TRUE) {

  # LRT - is there any difference in fit between the full and reduced models?
    cat("LIKELIHOOD RATIO TEST\n\n")
    print(anova(fm.reduced, fm.full))

  # Conditional F tests
    cat("\nCONDITIONAL F-TEST, TERMS = INTENSITY, TIME-TO-PEAK AND DISPERSION\n\n")
    ifelse(log.s==TRUE, print(anova(fm.full, Terms=c("I","mu_","Z"))), print(anova(fm.full, Terms=c("I","mu_","sigma"))))
    cat("\nCONDITIONAL F-TEST, TERMS = INTENSITY ONLY\n\n")
    print(anova(fm.full, Terms=c("I")))

}

# TESTS FOR TWO-SAMPLE SYMMETRIC GIRMM MODEL
# NOTE: Assumes lognormal distribution was specified for dispersion parameter (s)
# If normal distribution was specified instead, set log.s=FALSE

test.G2 <- function(fm.reduced, fm.full, log.s=TRUE) {

  # LRT - is there any difference in fit between the full and reduced models?
    cat("***** LIKELIHOOD RATIO TEST *****\n\n")
    print(anova(fm.reduced, fm.full))

  # Conditional F tests
    cat("\n***** CONDITIONAL F-TEST, TERMS = INTENSITY, TIME-TO-PEAK AND DISPERSION *****\n\n")
    ifelse(log.s==TRUE, print(anova(fm.full, Terms=c("I.Group","mu_.Group","Z.Group"))), print(anova(fm.full, Terms=c("I.Group","mu_.Group","sigma.Group"))))   # fit.full$fixDF$terms shows the term names
    cat("\n***** CONDITIONAL F-TEST, TERMS = INTENSITY ONLY *****\n\n")
    print(anova(fm.full, Terms=c("I.Group")))

}

# TESTS FOR TWO-SAMPLE ASYMMETRIC GIRMM MODEL
# NOTE: Assumes lognormal distribution was specified for dispersion parameter (s)
# If normal distribution was specified instead, set log.s=FALSE

test.G2.asymm <- function(fm.reduced, fm.full, log.s=TRUE) {
  
  # LRT - is there any difference in fit between the full and reduced models?
  cat("***** LIKELIHOOD RATIO TEST *****\n\n")
  print(anova(fm.reduced, fm.full))
  
  # Conditional F tests
  cat("\n***** CONDITIONAL F-TEST, TERMS = INTENSITY, TIME-TO-PEAK AND DISPERSION *****\n\n")
  ifelse(log.s==TRUE, print(anova(fm.full, Terms=c("I.Group","mu_.Group","logl.Group","logr.Group"))), print(anova(fm.full, Terms=c("I.Group","mu_.Group","l.Group","r.Group"))))   # fit.full$fixDF$terms shows the term names
  cat("\n***** CONDITIONAL F-TEST, TERMS = INTENSITY ONLY *****\n\n")
  print(anova(fm.full, Terms=c("I.Group")))
  cat("\n***** CONDITIONAL F-TEST, TERMS = TIME-TO-PEAK ONLY *****\n\n")
  print(anova(fm.full, Terms=c("mu_.Group")))
  
}


##################################################################
# ONSET, RECOVERY AND DURATION ESTIMATION AND BOOTSTRAPPED SE/CI #
##################################################################

# Takes a fitted one- or two-sample GIRMM model (nlmeObject) and returns a boostrapped sample

  get.bootstrap.sample.G <- function(fm) {
    
    df <- cbind(fm$groups, fm$fitted[,2], fm$resid[,2])
    colnames(df) <- c("Subject","Fitted","Residual")
    resid.sample <- numeric()
    for (i in levels(df$Subject)) {
      v <- df$Residual[df$Subject==i]
      resid.sample <- c(resid.sample, sample(v, size=length(v), replace=TRUE))
    }
    return(df$Fitted + resid.sample)
    
  }
  
# Takes a fitted one-sample GIRMM model and returns a numeric vector with components 
# correponding to estimates of response onset, recovery

  duration.G1 <- function(fm, baseline=-999) {
    
  # Set up
 	  s <- summary(fm)
 	  ifelse(baseline==-999, B.hat <- s$tTable["B",1], B.hat <- baseline)
    I.hat <- s$tTable["I",1]
 	  v <- s$sigma
 	  B.lower <- B.hat - 1.96*v; B.upper <- B.hat + 1.96*v
    response <- predict.G1(fm,level="group",0.1)
     
  # Estimate onset, recovery and duration of response
  	O <- ifelse(I.hat>0, min(response$DPI[which(response$Pred > B.upper)]), 
                             min(response$DPI[which(response$Pred < B.lower)]))
  	R <- ifelse(I.hat>0, max(response$DPI[which(response$Pred > B.upper)]),
  	                            max(response$DPI[which(response$Pred < B.lower)]))
    D <- R-O
    return(list(onset=O, recovery=R, duration=D, B.lower=B.lower, B.upper=B.upper))
     
  }
  
# Takes a fitted two-sample GIRMM model and returns a numeric vector with components 
# correponding to estimates of response onset, recovery
# NOTE: ASSUMES THAT REFERENCE GROUP IS LABELED "A" AND OTHER GROUP IS LABELED "B"

  duration.G2 <- function(fm, baseline=-999) {
    
    response <- predict.G2(fm, level="group", 0.1)
    response.A <- response[response$Group=="A",]
    response.B <- response[response$Group=="B",]    
    
  # Set up - reference group
 	  s <- summary(fm)
 	  ifelse(baseline==-999, B.hat <- s$tTable["B.(Intercept)",1], B.hat <- baseline)
    I.hat <- s$tTable["I.(Intercept)",1]
 	  v <- s$sigma + VarCorr(fm)["B.(Intercept)",1]
 	  B.lower.A <- B.hat - 1.96*v; B.upper.A <- B.hat + 1.96*v
     
  # Estimate onset, recovery and duration of response
  	O.A <- ifelse(I.hat>0, min(response.A$DPI[which(response.A$Pred > B.upper.A)]), 
                             min(response.A$DPI[which(response.A$Pred < B.lower.A)]))
  	R.A <- ifelse(I.hat>0, max(response.A$DPI[which(response.A$Pred > B.upper.A)]),
  	                            max(response.A$DPI[which(response.A$Pred < B.lower.A)]))
    D.A <- R.A-O.A
    
  # Set up - other group
 	  s <- summary(fm)
    ifelse(baseline==-999, B.hat <- s$tTable["B.(Intercept)",1] + s$tTable["B.GroupB",1], B.hat <- baseline)
    I.hat <- s$tTable["I.(Intercept)",1] + s$tTable["I.GroupB",1]
 	  v <- s$sigma + VarCorr(fm)["B.(Intercept)",1]
 	  B.lower.B <- B.hat - 1.96*v; B.upper.B <- B.hat + 1.96*v
     
  # Estimate onset, recovery and duration of response
  	O.B <- ifelse(I.hat>0, min(response.B$DPI[which(response.B$Pred > B.upper.B)]), 
                             min(response.B$DPI[which(response.B$Pred < B.lower.B)]))
  	R.B <- ifelse(I.hat>0, max(response.B$DPI[which(response.B$Pred > B.upper.B)]),
  	                            max(response.B$DPI[which(response.B$Pred < B.lower.B)]))
    D.B <- R.B-O.B
     
    return(list(onset.A=O.A, recovery.A=R.A, duration.A=D.A, 
                onset.B=O.B, recovery.B=R.B, duration.B=D.B,
                B.lower.A=B.lower.A, B.upper.A=B.upper.A,
                B.lower.B=B.lower.B, B.upper.B=B.upper.B))
     
  }

# Takes a fitted one-sample GIRMM model and calculates bootstrapped SE's and (1-alpha)*100% CI's 
# for mean onset, recovery duration of infection response

  bootstrap.duration.G1 <- function(fm, iter, alpha, baseline=-999) {

    require(nlme)
    boots <- data.frame(onset=numeric(), recovery=numeric(), duration=numeric())
    for (i in 1:iter) {
      displayProgressBar(i,iter,10)
      bs <- get.bootstrap.sample.G(fm)
    	dep.var <- all.vars(getCall(fm))[1]
    	bootstrap.df <<- getData(fm)
    	bootstrap.df[dep.var] <- bs
      try({
    	  fm.new <- update(fm, data=bootstrap.df)
        ORD <- duration.G1(fm.new, baseline)
    	  boots <- rbind(boots, c(ORD$onset, ORD$recovery, ORD$duration))
      })
    }
    
    # Clean up
      rm(list=c("bootstrap.df"), pos=".GlobalEnv")
    
    # Calculate bootstrapped SE's and CI's and return in a list along with number of successful fits
      colnames(boots) <- c("O","R","D")
      results <- data.frame(estimate=numeric(), SE=numeric(), lower=numeric(), upper=numeric())
      est <- duration.G1(fm, baseline)
      results <- rbind(results, c(est$onset, sd(boots$O), est$onset-qnorm(1-alpha/2)*sd(boots$O), est$onset+qnorm(1-alpha/2)*sd(boots$O)))
      results <- rbind(results, c(est$recovery, sd(boots$R), est$recovery-qnorm(1-alpha/2)*sd(boots$R), est$recovery+qnorm(1-alpha/2)*sd(boots$R)))
      results <- rbind(results, c(est$duration, sd(boots$D), est$duration-qnorm(1-alpha/2)*sd(boots$D), est$duration+qnorm(1-alpha/2)*sd(boots$D)))
      rownames(results) <- c("Onset","Recovery","Duration")
      colnames(results) <- c("Estimate","SE","Lower","Upper")
      return(list(table=results, fits=nrow(boots)))
  }
  
# Takes a fitted two-sample GIRMM model and calculates bootstrapped SE's and (1-alpha)*100% CI's 
# for mean onset, recovery duration of infection response

  bootstrap.duration.G2 <- function(fm, iter, alpha, baseline=-999) {

    require(nlme)
    boots <- data.frame(onset.A=numeric(), recovery.A=numeric(), duration.A=numeric(),
                          onset.B=numeric(), recovery.B=numeric(), duration.B=numeric())
    for (i in 1:iter) {
      displayProgressBar(i,iter,10)
    	bs <- get.bootstrap.sample.G(fm)
    	dep.var <- all.vars(getCall(fm))[1]
    	bootstrap.df <<- getData(fm)
    	bootstrap.df[dep.var] <- bs
      try({
        fm.new <- update(fm, data=bootstrap.df)
        ORD <- duration.G2(fm.new, baseline)
        boots <- rbind(boots, c(ORD$onset.A, ORD$recovery.A, ORD$duration.A,ORD$onset.B, ORD$recovery.B, ORD$duration.B))
      })
    }
    
    # Clean up
      rm(list=c("bootstrap.df"), pos=".GlobalEnv")
    
    # Calculate bootstrapped SE's and CI's and return in a list along with number of successful fits
      colnames(boots) <- c("O.A","O.B","R.A","R.B","D.A","D.B")
      results.A <- data.frame(estimate=numeric(), SE=numeric(), lower=numeric(), upper=numeric())
      est <- duration.G2(fm, baseline)
      results.A <- rbind(results.A, c(est$onset.A, sd(boots$O.A), est$onset.A-qnorm(1-alpha/2)*sd(boots$O.A), est$onset.A+qnorm(1-alpha/2)*sd(boots$O.A)))
      results.A <- rbind(results.A, c(est$recovery.A, sd(boots$R.A), est$recovery.A-qnorm(1-alpha/2)*sd(boots$R.A), est$recovery.A+qnorm(1-alpha/2)*sd(boots$R.A)))
      results.A <- rbind(results.A, c(est$duration.A, sd(boots$D.A), est$duration.A-qnorm(1-alpha/2)*sd(boots$D.A), est$duration.A+qnorm(1-alpha/2)*sd(boots$D.A)))
      rownames(results.A) <- c("Onset","Recovery","Duration")
      colnames(results.A) <- c("Estimate","SE","Lower","Upper")
      results.B <- data.frame(estimate=numeric(), SE=numeric(), lower=numeric(), upper=numeric())
      results.B <- rbind(results.B, c(est$onset.B, sd(boots$O.B), est$onset.B-qnorm(1-alpha/2)*sd(boots$O.B), est$onset.B+qnorm(1-alpha/2)*sd(boots$O.B)))
      results.B <- rbind(results.B, c(est$recovery.B, sd(boots$R.B), est$recovery.B-qnorm(1-alpha/2)*sd(boots$R.B), est$recovery.B+qnorm(1-alpha/2)*sd(boots$R.B)))
      results.B <- rbind(results.B, c(est$duration.B, sd(boots$D.B), est$duration.B-qnorm(1-alpha/2)*sd(boots$D.B), est$duration.B+qnorm(1-alpha/2)*sd(boots$D.B)))
      rownames(results.B) <- c("Onset","Recovery","Duration")
      colnames(results.B) <- c("Estimate","SE","Lower","Upper")
      return(list(table.A=results.A, table.B=results.B, fits=nrow(boots)))
    
  }

# Takes a fitted one-sample GIRMM model and calculates bootstrapped SE's and (1-alpha)*100% CI's 
# for the mean response at each time point

bootstrap.mean.G1 <- function(fm, iter, alpha, baseline=-999) {
  
  require(nlme)
  boots <- data.frame(onset=numeric(), recovery=numeric(), duration=numeric())
  for (i in 1:iter) {
    displayProgressBar(i,iter,10)
    bs <- get.bootstrap.sample.G(fm)
    dep.var <- all.vars(getCall(fm))[1]
    bootstrap.df <<- getData(fm)
    bootstrap.df[dep.var] <- bs
    try({
      fm.new <- update(fm, data=bootstrap.df)
      #ORD <- duration.G1(fm.new, baseline)
      boots <- rbind(boots, predict.G1(fm.new,level="group")$Pred)
    })
  }
  
  # Clean up
  rm(list=c("bootstrap.df"), pos=".GlobalEnv")
  
    intervals <- t(apply(boots,2,function(x) quantile(x, probs=c(0.25,0.975))))
    colnames(intervals) <- c("lower","upper")
    DPI.xx <- seq(min(getData(fit.full_A)$DPI),max(getData(fit.full_A)$DPI),by=1)
    #polygon(c(DPI.xx,rev(DPI.xx)),c(intervals[,"lower"],rev(intervals[,"upper"])),col=rgb(0,0,1,alpha=.15),border=F)
    polygon(c(DPI.xx,rev(DPI.xx)),c(intervals[,"lower"],rev(intervals[,"upper"])),col=rgb(1,0.65,0,alpha=.15),border=F)
  
  # Calculate bootstrapped SE's and CI's and return in a list along with number of successful fits
#     colnames(boots) <- c("B","I","mu","z")
#     results <- data.frame(lower=numeric(), upper=numeric())
#     results <- rbind(results, quantile(boots$B,probs=c(0.025,0.975)))
#     results <- rbind(results, quantile(boots$I,probs=c(0.025,0.975)))
#     results <- rbind(results, quantile(boots$mu,probs=c(0.025,0.975)))
#     results <- rbind(results, quantile(boots$z,probs=c(0.025,0.975)))
#     rownames(results) <- c("B","I","mu","z")
#     colnames(results) <- c("Lower","Upper")
#     return(list(table=results, fits=nrow(boots)))
#  return(intervals)
}

# Progress bar
  displayProgressBar <- function(i,nI,increments=30)
  {
    if (nI < increments) increments <- nI
    if (i==1)
    {
      cat("Progress:\n")
      cat("|",rep(" ",increments),"|\n ",sep="")
    }
    display.iterates <- ceiling((1:increments) * (nI/increments))
    if (is.element(i,display.iterates)) cat("*")
    if (Sys.info()["sysname"]=="Windows") flush.console()
    if (i==nI) cat("\n")
  }

