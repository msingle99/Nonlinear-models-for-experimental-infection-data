# For testing
test <- function() {
  df<-myData; true.params.A<-group.A.params; true.params.B<-group.B.params; max.iter=5000
}

fit.GIRMM.S2 <- function(df, true.params.A, true.params.B, max.iter, plot.mean=TRUE) {
  
  df.g <- groupedData(Response~DPI|Subject, data=df)
  A <- true.params.A; B <- true.params.B
  Z.start.2 <- ifelse(B$s1$mean==A$s1$mean,1,log(B$s1$mean - A$s1$mean))
  fit <- NULL; fit.type=0;
  good.fit = FALSE
  
  try.fit.full <- try({
    nlme(Response ~ xGauss.log(DPI,mu_,Z,B,I), method="ML",
         data=df.g,
         #groups=~Subject,
         fixed=B+I+mu_+Z~Group,
         random=pdDiag(B+I+mu_+Z~1),
         #start = c(B = c(20,0), I = c(30,10), mu_ = c(45,10), sigma=c(log(3),1)),
         start = c(B = c(A$B$mean,B$B$mean-A$B$mean), I = c(A$I$mean,B$I$mean-A$I$mean), mu_ = c(A$mu$mean,B$mu$mean-A$mu$mean), Z=c(log(A$s1$mean), Z.start.2)),
         #start=c(B=20,I=35,mu=50,Z=1),
         control=nlmeControl(maxIter=max.iter))
  })
 
   if (class(try.fit.full) != "try-error") {
     try.int.full <- try({intervals(try.fit.full)})
     if (class(try.int.full) != "try-error") {
       vc <- VarCorr(try.fit.full)            # Check for unrealistically small random effect variances ("degenerate dist's" for RE's)
       storage.mode(vc) <- "numeric"
       if (!(TRUE %in% (vc[1:4,"StdDev"]<0.01))) {
         print("full model good fit")
         fit <- try.fit.full; fit.type=1
         good.fit = TRUE
       }
     } 
   }
   
  if (good.fit == FALSE) {
     try.fit.BIM <- try({
       nlme(Response ~ xGauss.log(DPI,mu_,Z,B,I), method="ML",
            data=df.g,
            fixed=B+I+mu_+Z~Group,
            #random=list(Group=pdDiag(B+I+mu_~1), Subject=pdDiag(B+I+mu_~1)),
            random=pdDiag(B+I+mu_~1),
            start = c(B = c(A$B$mean,B$B$mean-A$B$mean), I = c(A$I$mean,B$I$mean-A$I$mean), mu_ = c(A$mu$mean,B$mu$mean-A$mu$mean), Z=c(log(A$s1$mean),Z.start.2)),
            control=nlmeControl(maxIter=max.iter))
     })
     if (class(try.fit.BIM) != "try-error") {
       try.int.BIM <- try({intervals(try.fit.BIM)})
       if (class(try.int.BIM) != "try-error") {
         vc <- VarCorr(try.fit.BIM)            # Check for unrealistically small random effect variances
         storage.mode(vc) <- "numeric"
         if (!(TRUE %in% (vc[1:2,"StdDev"]<0.01))) {
           print("B, I and M good fit")
           fit <- try.fit.BIM; fit.type=2
           good.fit = TRUE
         }
       }
     } 
   }
    
   if (good.fit == FALSE) {
     try.fit.BI <- try({
       nlme(Response ~ xGauss.log(DPI,mu_,Z,B,I), method="ML",
            data=df.g,
            fixed=B+I+mu_+Z~Group,
            random=pdDiag(B+I~1),
            start = c(B = c(A$B$mean,B$B$mean-A$B$mean), I = c(A$I$mean,B$I$mean-A$I$mean), mu_ = c(A$mu$mean,B$mu$mean-A$mu$mean), Z=c(log(A$s1$mean),Z.start.2)),
            control=nlmeControl(maxIter=max.iter))
     })
     if (class(try.fit.BI) != "try-error") {
       try.int.BI <- try({intervals(try.fit.BI)})
       if (class(try.int.BI) != "try-error") {
         vc <- VarCorr(try.fit.BI)            # Check for unrealistically small random effect variances
         storage.mode(vc) <- "numeric"
         if (!(TRUE %in% (vc[1:2,"StdDev"]<0.01))) {
           print("B and I good fit")
           fit <- try.fit.BI; fit.type=3
           good.fit = TRUE
         }
       }
     } 
   }
   
   if (good.fit == FALSE) {
     try.fit.I <- try({
       nlme(Response ~ xGauss.log(DPI,mu_,Z,B,I), method="ML",
            data=df.g,
            fixed=B+I+mu_+Z~Group,
            random=pdDiag(I~1),
            start = c(B = c(A$B$mean,B$B$mean-A$B$mean), I = c(A$I$mean,B$I$mean-A$I$mean), mu_ = c(A$mu$mean,B$mu$mean-A$mu$mean), Z=c(log(A$s1$mean),Z.start.2)),
            control=nlmeControl(maxIter=max.iter))
     })
     if (class(try.fit.I) != "try-error") {
       try.int.I <- try({intervals(try.fit.I)})
       if (class(try.int.I) != "try-error") {
         vc <- VarCorr(try.fit.I)            # Check for unrealistically small random effect variances
         storage.mode(vc) <- "numeric"
         if (!(TRUE %in% (vc[1,"StdDev"]<0.01))) {
           print("I only good fit")
           fit <- try.fit.I; fit.type=4
           good.fit = TRUE
         }
       } 
     }
   }
   
#   if (good.fit == FALSE) {
#     try.fit.nls <- try({
#       nls(Response ~ xGauss.log(DPI,mu_,Z,B,I),
#           data=df.g, 
#           start = c(B = A$B$mean, I = A$I$mean, mu_ = A$mu$mean, Z=log(A$s1$mean)),
#           control=nlmeControl(maxIter=max.iter))
#     })
#     if (class(try.fit.nls) != "try-error") {
#       print("No random effects good fit")
#       fit <- try.fit.nls; fit.type=5
#       good.fit = TRUE
#       fit.reduced <- try({nls(Response ~ xGauss.log(DPI,0,0,B,0), data=df.g, start=c(B=A$B$mean))})
#     }
#   } 
#   
   if (good.fit == FALSE) {
     # Both fits failed; save "Failed" result but not fit results
     print("All four models failed")
     fit.type=6
  }
   
   ## EXTRACT RESULTS FROM FITTED MODEL AND RETURN
 
     if (class(fit)[1]=="nlme") {
         
       FE <- fixef(fit)
#       lines(z, xGauss(z, FE["mu_"], exp(FE["Z"]), FE["B"], FE["I"]), type="l", lwd=2, col="red")
     # Standard errors
       SE.B.Intercept <- summary(fit)$tTable["B.(Intercept)","Std.Error"]
       SE.I.Intercept <- summary(fit)$tTable["I.(Intercept)","Std.Error"]
       SE.B.GroupB <- summary(fit)$tTable["B.GroupB","Std.Error"]
       SE.I.GroupB <- summary(fit)$tTable["I.GroupB","Std.Error"]

     # Significant response?
       p.value <- anova(fit, Terms=c("I.Group","mu_.Group","Z.Group"))$p
       reject1 <- (p.value < 0.05)
       p.value <- anova(fit, Terms=c("I.Group"))$p
       reject2 <- (p.value < 0.05)
     # Full width at half maximum
#       #FWHM <- FWHM.nlme(fit)
       FWHM=1
       r <- c(fixef(fit),FWHM,reject1,reject2, SE.B.Intercept, SE.I.Intercept, SE.B.GroupB, SE.I.GroupB, fit.type)
       
#     } else if(class(fit)[1]=="nls") {
#        FE <- coef(fit)
#        lines(z, xGauss(z, FE["mu_"], exp(FE["Z"]), FE["B"], FE["I"]), type="l", lwd=2, col="red")
#       # Standard errors
#         SE.B <- sqrt(diag(vcov(fit))["B"])
#         SE.I <- sqrt(diag(vcov(fit))["I"])
#       # Significant response?
#         if (class(fit.reduced) != "try-error") {
#           reject1 <- as.numeric((anova(fit.reduced, fit)[2,"Pr(>F)"] <0.05))
#           reject2 <- NA 
#         } else { reject1=NA; reject2=NA }
#       # Full width at half maximum
#         FWHM <- FWHM.nls(fit)
#       r <- c(coef(fit),FWHM,reject1,reject2,SE.B.Intercept,SE.I.Intercept, SE.B.GroupB, SE.I.GroupB,fit.type)
#         
     } else {r <- NULL}
 
     return(r)
  
}

fit.2wayANOVA.lme <- function(df) {
  
  # Fit the model
    df$DPI.f <- factor(df$DPI)
    fit <- lme(Response~Group*DPI.f, random=~1|Subject, data=df)
    # aov(Response~Group*DPI.f + Error(Subject), df)
  
  # Plot the estimated mean function for ANOVA
     #pred.ANOVA <- predict(fit, level=0)[1:length(t)]
     #lines(t, pred.ANOVA, col="darkgreen", type="b", lwd=2)
  
  # "Estimates" baseline, intensity and time-to-peak, plus SE's for B and I
    x <- summary(fit)$tTable   # Estimated deviations from Intercept (baseline?) and SE's
    B.Intercept <- x["(Intercept)","Value"]; SE.B.Intercept <- x["(Intercept)","Std.Error"];
    B.GroupB <- x["GroupB","Value"] + B.Intercept; SE.B.GroupB <- x["GroupB","Std.Error"];
    GroupA.estimates <- x[3:(length(t)+1),]
    GroupB.estimates <- x[(length(t)+2):nrow(x),]
    I.Intercept <- max(GroupA.estimates[,"Value"]); SE.I.Intercept <- max(GroupA.estimates[,"Std.Error"])
    I.GroupB <- max(GroupB.estimates[,"Value"] + GroupA.estimates[,"Value"]); 
    SE.I.GroupB <- max(GroupB.estimates[,"Std.Error"])
    mu.Intercept <- t[which.max(GroupA.estimates[,"Value"])+1]    
    mu.GroupB <- t[which.max(GroupB.estimates[,"Value"])+1]
 
  # Full width at half maximum
    #FWHM <- FWHM.anova(x, unique(getData(fit)$DPI), B, I, mu)
    FWHM <- 1
         
  # Test for difference in group responses 
    reject.post <- (anova(fit)[["p-value"]][4] < 0.05)
    
   return(c(B.Intercept, B.GroupB, I.Intercept, I.GroupB, mu.Intercept, mu.GroupB, FWHM, as.numeric(reject.post),SE.B.Intercept, SE.B.GroupB, SE.I.Intercept,SE.I.GroupB))
    
}

####################
## REPORT RESULTS ##
####################

summarize <- function(true, est) {
  true.mean <- true
  error <- est - true.mean
  error2 <- error^2
  bias <- mean(error, na.rm=TRUE)
  var <- var(est, na.rm=TRUE)
  #MSE <- mean(error^2)
  MSE <- bias^2 + var
  SE.bias <- sd(error, na.rm=TRUE)/sqrt(length(est))
  summary <- list(est=mean(est), bias=bias, SE.bias=SE.bias, var=var, MSE=MSE)
  return(summary)
}

summarize.SE <- function(est) {
  return(list(est=mean(est), SE=sd(est)/sqrt(sum(!is.na(est)))))
}

write.report <- function(true.params.A, true.params.B, nlme.est, ANOVA.est, rpt.name) {
	# Write results out to a spreadsheet
  
  # nlme results
    B.Intercept.summary <- summarize(true.params.A$B$mean, nlme.est$B.Intercept)
    B.GroupB.summary <- summarize(true.params.B$B$mean, nlme.est$B.Intercept+nlme.est$B.GroupB)   
    I.Intercept.summary <- summarize(true.params.A$I$mean, nlme.est$I.Intercept)
    I.GroupB.summary <- summarize(true.params.B$I$mean, nlme.est$I.Intercept+nlme.est$I.GroupB)    
    mu.Intercept.summary <- summarize(true.params.A$mu$mean, nlme.est$mu_.Intercept)
    mu.GroupB.summary <- summarize(true.params.B$mu$mean, nlme.est$mu_.Intercept+nlme.est$mu_.GroupB)    
    s.Intercept.summary <- summarize(true.params.A$s1$mean, exp(nlme.est$Z.Intercept))
    s.GroupB.summary <- summarize(true.params.B$s1$mean, exp(nlme.est$Z.Intercept+nlme.est$Z.GroupB))    
    FWHM.summary <- summarize(2.355*true.params.A$s1$mean, nlme.est$FWHM)
    SE.B.Intercept.summary <- summarize.SE(nlme.est$SE.B.Intercept)
    SE.B.GroupB.summary <- summarize.SE(nlme.est$SE.B.GroupB)   
    SE.I.Intercept.summary <- summarize.SE(nlme.est$SE.I.Intercept)
    SE.I.GroupB.summary <- summarize.SE(nlme.est$SE.I.GroupB)    
    nlme.summary <- list(B.Intercept=B.Intercept.summary, B.GroupB=B.GroupB.summary, I.Intercept=I.Intercept.summary, I.GroupB=I.GroupB.summary, mu.Intercept=mu.Intercept.summary, mu.GroupB=mu.GroupB.summary, s.Intercept=s.Intercept.summary, s.GroupB=s.GroupB.summary, FWHM=FWHM.summary, SE.B.Intercept=SE.B.Intercept.summary, SE.B.GroupB=SE.B.GroupB.summary, SE.I.Intercept=SE.I.Intercept.summary, SE.I.GroupB=SE.I.GroupB.summary, reject1=mean(nlme.est$reject1, na.rm=TRUE),reject2=mean(nlme.est$reject2, na.rm=TRUE))
    
    cat("\n*** NLME summary ***\n\n")
    print(unlist(nlme.summary))
    nlme.vector <- unlist(nlme.summary) # Convert list to numeric vector
     write.table(nlme.vector, file="/Users/MikeZhen/Dropbox/DissertationMS/eav-ms/response-curves/Chapter 2/SIMS/test.csv", sep=",")
    
  # ANOVA results
    B.Intercept.summary <- summarize(true.params.A$B$mean, ANOVA.est$B.Intercept)
    B.GroupB.summary <- summarize(true.params.B$B$mean, ANOVA.est$B.GroupB)
    I.Intercept.summary <- summarize(true.params.A$I$mean, ANOVA.est$I.Intercept)
    I.GroupB.summary <- summarize(true.params.B$I$mean, ANOVA.est$I.GroupB)
    mu.Intercept.summary <- summarize(true.params.A$mu$mean, ANOVA.est$mu_.Intercept)
    mu.GroupB.summary <- summarize(true.params.B$mu$mean, ANOVA.est$mu_.GroupB)
    FWHM.summary <- summarize(true.params.A$I$mean, ANOVA.est$FWHM)
    SE.B.Intercept.summary <- summarize.SE(ANOVA.est$SE.B.Intercept)
    SE.B.GroupB.summary <- summarize.SE(ANOVA.est$SE.B.GroupB)
    SE.I.Intercept.summary <- summarize.SE(ANOVA.est$SE.I.Intercept)
    SE.I.GroupB.summary <- summarize.SE(ANOVA.est$SE.I.GroupB)
    ANOVA.summary <- list(B.Intercept=B.Intercept.summary, B.GroupB=B.GroupB.summary, I.Intercept=I.Intercept.summary, I.GroupB=I.GroupB.summary, mu.Intercept=mu.Intercept.summary, mu.GroupB=mu.GroupB.summary, FWHM=FWHM.summary, SE.B.Intercept=SE.B.Intercept.summary, SE.B.GroupB=SE.B.GroupB.summary, SE.I.Intercept=SE.I.Intercept.summary, SE.I.GroupB=SE.I.GroupB.summary,  reject.post=mean(ANOVA.est$reject.post, na.rm=TRUE))
    
    cat("\n*** ANOVA summary ***\n\n")
    print(unlist(ANOVA.summary))
    
}


FWHM.nlme <- function(fm) {
  
  half.max <- fixef(fm)["B"] + fixef(fm)["I"]/2
  DPI <- unique(getData(fm)$DPI)
  
  # Estimated mean curve from fitted model
  DPI.xx <- seq(min(DPI), max(DPI), by=0.1)           
  df.xx <- data.frame(DPI=DPI.xx)    
  pred.fit <- predict(fm, newdata=df.xx, level=0)
  # Split DPI.xx into pre- and post-peak parts
  DPI.pre.mu <- DPI.xx[DPI.xx <= fixef(fm)["mu_"]]
  DPI.post.mu <- DPI.xx[DPI.xx > fixef(fm)["mu_"]]
  
  # Split predicted values into pre- and post-peak parts
  pred.pre.mu <- pred.fit[1:length(DPI.pre.mu)]
  pred.post.mu <- pred.fit[(length(DPI.pre.mu)+1):length(DPI.xx)]
  
  # Find DPI.xx values on either side of peak where prediction is closest to the half max value
  FWHM.lower <- DPI.pre.mu[which.min(abs(pred.pre.mu-half.max))]
  FWHM.upper <- DPI.post.mu[which.min(abs(pred.post.mu-half.max))]
    
  # Return estimated FWHM
  return(FWHM.upper-FWHM.lower)
  
}

FWHM.nls <- function(fm) {
  
  half.max <- coef(fm)["B"] + coef(fm)["I"]/2
  DPI <- unique(getData(fm)$DPI)
  
  # Estimated mean curve from fitted model
  DPI.xx <- seq(min(DPI), max(DPI), by=0.1)           
  df.xx <- data.frame(DPI=DPI.xx)       
  pred.fit <- predict(fm, newdata=df.xx)
  
  # Split DPI.xx into pre- and post-peak parts
  DPI.pre.mu <- DPI.xx[DPI.xx <= coef(fm)["mu_"]]
  DPI.post.mu <- DPI.xx[DPI.xx > coef(fm)["mu_"]]
  
  # Split predicted values into pre- and post-peak parts
  pred.pre.mu <- pred.fit[1:length(DPI.pre.mu)]
  pred.post.mu <- pred.fit[(length(DPI.pre.mu)+1):length(DPI.xx)]
  
  # Find DPI.xx values on either side of peak where prediction is closest to the half max value
  FWHM.lower <- DPI.pre.mu[which.min(abs(pred.pre.mu-half.max))]
  FWHM.upper <- DPI.post.mu[which.min(abs(pred.post.mu-half.max))]
  
  # Return estimated FWHM
  return(FWHM.upper-FWHM.lower)
  
}

FWHM.anova <- function(x, DPI, B, I, mu) {
  
  if (mu > min(DPI) & mu < max(DPI)) {
    
    half.max <- B + I/2
    
    # Estimates at each time point based on fitted model
    pred.fit <- rep(NaN, length(DPI))
    pred.fit[1] <- x[1,"Value"]
    pred.fit[2:length(DPI)] <- pred.fit[1] + x[2:length(DPI),"Value"]
    
    # Split DPI into pre- and post-peak parts
    DPI.pre.mu <- DPI[DPI <= mu]
    DPI.post.mu <- DPI[DPI > mu]
    
    # Split predicted values into pre- and post-peak parts
    pred.pre.mu <- pred.fit[1:length(DPI.pre.mu)]
    pred.post.mu <- pred.fit[(length(DPI.pre.mu)+1):length(DPI)]
    
    # Find DPI values on either side of peak where estimate is closest to the half max value
    FWHM.lower <- DPI.pre.mu[which.min(abs(pred.pre.mu-half.max))]
    FWHM.upper <- DPI.post.mu[which.min(abs(pred.post.mu-half.max))]
    
    # Return estimated FWHM
    return(FWHM.upper-FWHM.lower)
    
  } else return(NA)
  
}

parabolic <- function(DPI,h,a,B,k) {
  
  # For baseline exactly 0 for all subjects - NOT USED
  #  return(ifelse(DPI<(h-sqrt(-k/a)) | DPI>(h+sqrt(-k/a)), 0, a*(DPI-h)^2+k))
  
  # For possibly nonzero baseline for all subjects; a and k either both - or both +
  # Note: to get "onset" at time t for a fixed h and k, set a=k/(t-h)^2
  return(ifelse(DPI<(h-sqrt(k/a)) | DPI>(h+sqrt(k/a)), B, a*(DPI-h)^2+(B-k)))
  
}

# Function that generates a list of N parabolic response curve parameter vectors 
# from an underlying population parabolic response curve parameter vector. The list 
# corresponds to a sample of N stallions selected at random from the population
# of all stallions of a particular breed (who are infected with a particular
# type, strain and dose of virus)

get.two.parabolic.samples <- function(N, t, A, B) {
  
  time.points <- length(t)
  
  # Make vector for Stallion variable
  s <- rep(1,time.points)
  for (i in 2:N) s <- c(s,rep(i,time.points))
  
  # Var
  var.A <- NULL; var.B <- NULL
  var.A.all <- NULL; var.B.all <- NULL
  
  # Draw N fake horses from population A and simulate their data
  group.A <- vector("list", N)
  for (i in 1:N) {
    # Simulate a horse (i.e. its parameters drawn from the population values)
    # Note: don't forget to require that B_i + I_i > 0 if variable is a sperm response
    #group.A[[i]] <- mvrnorm(1,theta.A,diag(sd.theta.A))
    if (A[["a"]][["mean"]]==0 & A[["h"]][["mean"]]==0 & A[["k"]][["mean"]]==0) {
      repeat {
        group.A[[i]][3] <- rnorm(1, A[["B"]][["mean"]], A[["B"]][["sd"]])
        if (group.A[[i]][3] >= A[["B"]][["min"]] & group.A[[i]][3] <= A[["B"]][["max"]]) break
      }
      group.A[[i]][1]=0; group.A[[i]][2]=0; group.A[[i]][4]=0;
    } 
    else {
      repeat { 
        group.A[[i]][1] <- rnorm(1, A[["h"]][["mean"]], A[["h"]][["sd"]])
        if (group.A[[i]][1] >= A[["h"]][["min"]] & group.A[[i]][1] <= A[["h"]][["max"]]) break
      }
      repeat { 
        group.A[[i]][2] <- rnorm(1, A[["a"]][["mean"]], A[["a"]][["sd"]])
        if (group.A[[i]][2] >= A[["a"]][["min"]] & group.A[[i]][2] <= A[["a"]][["max"]]) break
      }
      repeat { 
        group.A[[i]][3] <- rnorm(1, A[["B"]][["mean"]], A[["B"]][["sd"]])
        if (group.A[[i]][3] >= A[["B"]][["min"]] & group.A[[i]][3] <= A[["B"]][["max"]]) break
      }
      repeat { 
        group.A[[i]][4] <- rnorm(1, A[["k"]][["mean"]], A[["k"]][["sd"]])
        if (group.A[[i]][4] >= A[["k"]][["min"]] & group.A[[i]][4] <= A[["k"]][["max"]]) break
      }
    }
    # Simulate a specific response for this horse from its underlying response curve 
    if (A[["a"]][["mean"]]==0 & A[["h"]][["mean"]]==0 & A[["k"]][["mean"]]==0) {
      for (j in 1:time.points) {
        repeat { 
          var.A[j] <- group.A[[i]][3] + rnorm(1, sd=A[["Y"]][["sd"]])
          if (var.A[j] >= A[["Y"]][["min"]] & var.A[j] <= A[["Y"]][["max"]]) break
        }
      }
    }
    
    if (A[["a"]][["mean"]]!=0 | A[["h"]][["mean"]]!=0 | A[["k"]][["mean"]]!=0) {            
      for (j in 1:time.points) {
        repeat { 
          var.A[j] <- parabolic(t[j],group.A[[i]][1],group.A[[i]][2],group.A[[i]][3],group.A[[i]][4]) + rnorm(1, sd=A[["Y"]][["sd"]])
          if (var.A[j] >= A[["Y"]][["min"]] & var.A[j] <= A[["Y"]][["max"]]) break
        }
      }
    }
    var.A.all <- c(var.A.all,var.A)
  }  
  
  # Draw N fake horses from population B and simulate their data
  group.B <- vector("list", N)
  for (i in 1:N) {
    #group.B[[i]] <- mvrnorm(1,theta.B,diag(sd.theta.B))
    #var <- c(var, xGauss(t,group.B[[i]][1],group.B[[i]][2],group.B[[i]][3],group.B[[i]][4])+rnorm(time.points,sd=sd.var.B))
    # Simulate a horse (i.e. its parameters drawn from the population values)
    # Note: don't forget to require that B_i + I_i > 0 if variable is a sperm response
    #group.A[[i]] <- mvrnorm(1,theta.A,diag(sd.theta.A))
    repeat {
      group.B[[i]][1] <- rnorm(1, B[["h"]][["mean"]], B[["h"]][["sd"]])
      if (group.B[[i]][1] >= B[["h"]][["min"]] & group.B[[i]][1] <= B[["h"]][["max"]]) break
    }
    repeat {
      group.B[[i]][2] <- rnorm(1, B[["a"]][["mean"]], B[["a"]][["sd"]])
      if (group.B[[i]][2] >= B[["a"]][["min"]] & group.B[[i]][2] <= B[["a"]][["max"]]) break
    }
    repeat {
      group.B[[i]][3] <- rnorm(1, B[["B"]][["mean"]], B[["B"]][["sd"]])
      if (group.B[[i]][3] >= B[["B"]][["min"]] & group.B[[i]][3] <= B[["B"]][["max"]]) break
    }
    repeat {
      group.B[[i]][4] <- rnorm(1, B[["k"]][["mean"]], B[["k"]][["sd"]])
      if (group.B[[i]][4] >= B[["k"]][["min"]] & group.B[[i]][4] <= B[["k"]][["max"]]) break
    }
    
    # Simulate a response curve for this horse  
    for (j in 1:time.points) {
      repeat {
        var.B[j] <- parabolic(t[j],group.B[[i]][1],group.B[[i]][2],group.B[[i]][3],group.B[[i]][4]) + rnorm(1, sd=B[["Y"]][["sd"]])
        if (var.B[j] >= B[["Y"]][["min"]] & var.B[j] <= B[["Y"]][["max"]]) break
      }
    }
    var.B.all <- c(var.B.all,var.B)
  }
  
  # Create data frame
  df <- data.frame(Group=c(rep('A',N*time.points), rep('B',N*time.points)),
                   Stallion=s,
                   DPI=rep(t,2*N),
                   Response=c(var.A.all, var.B.all))
  return(list(data=df, A=group.A, B=group.B))
}

# Function that generates a list of N extended Gaussian parameter vectors from
# an underlying population extended Guassian parameter vector. The list 
# corresponds to a sample of N stallions selected at random from the population
# of all stallions of a particular breed (who are infected with a particular
# type, strain and dose of virus)

get.two.Gaussian.samples <- function(N, t, A, B, sym) {
  
  time.points <- length(t)
  
  # Make vector for Subject variable
  s <- rep(1,time.points)
  for (i in 2:(2*N)) s <- c(s,rep(i,time.points))
  
  # Var
  var.A <- NULL; var.B <- NULL
  var.A.all <- NULL; var.B.all <- NULL
  
  # Draw N fake horses from population A and simulate their data
  group.A <- vector("list", N)
  for (i in 1:N) {
    # Simulate a horse (i.e. its parameters drawn from the population values)
    # Note: don't forget to require that B_i + I_i > 0 if variable is a sperm response
    #group.A[[i]] <- mvrnorm(1,theta.A,diag(sd.theta.A))
    repeat {
      group.A[[i]][1] <- rnorm(1, A[["mu"]][["mean"]], A[["mu"]][["sd"]])
      if (group.A[[i]][1] >= A[["mu"]][["min"]] & group.A[[i]][1] <= A[["mu"]][["max"]]) break
    }
    repeat {
      group.A[[i]][2] <- rlnorm(1, A[["s1"]][["mean"]], A[["s1"]][["sd"]])
      if (group.A[[i]][2] >= A[["s1"]][["min"]] & group.A[[i]][2] <= A[["s1"]][["max"]]) break
    }
    repeat {
      group.A[[i]][3] <- rlnorm(1, A[["s2"]][["mean"]], A[["s2"]][["sd"]])
      if (group.A[[i]][3] >= A[["s2"]][["min"]] & group.A[[i]][3] <= A[["s2"]][["max"]]) break
    }
    repeat {
      group.A[[i]][4] <- rnorm(1, A[["B"]][["mean"]], A[["B"]][["sd"]])
      if (group.A[[i]][4] >= A[["B"]][["min"]] & group.A[[i]][4] <= A[["B"]][["max"]]) break
    }
    repeat {
      group.A[[i]][5] <- rnorm(1, A[["I"]][["mean"]], A[["I"]][["sd"]])
      if (group.A[[i]][5] >= A[["I"]][["min"]] & group.A[[i]][5] <= A[["I"]][["max"]]) break
    }
    
    # Simulate a response curve for this horse  
    for (j in 1:time.points) {
      if (sym=="FALSE") {
        repeat {
          var.A[j] <- xGauss2(t[j],group.A[[i]][1],group.A[[i]][2],group.A[[i]][3],group.A[[i]][4],group.A[[i]][5]) + rnorm(1, sd=A[["Y"]][["sd"]])
          if (var.A[j] >= A[["Y"]][["min"]] & var.A[j] <= A[["Y"]][["max"]]) break
        }
      }
      if (sym=="TRUE") {
        repeat {
          var.A[j] <- xGauss(t[j],group.A[[i]][1],log(group.A[[i]][2]),group.A[[i]][4],group.A[[i]][5]) + rnorm(1, sd=A[["Y"]][["sd"]])
          if (var.A[j] >= A[["Y"]][["min"]] & var.A[j] <= A[["Y"]][["max"]]) break
        }
      }
    }
    var.A.all <- c(var.A.all,var.A)
  }  
  
  # Draw N fake horses from population B and simulate their data
  group.B <- vector("list", N)
  for (i in 1:N) {
    #group.B[[i]] <- mvrnorm(1,theta.B,diag(sd.theta.B))
    #var <- c(var, xGauss(t,group.B[[i]][1],group.B[[i]][2],group.B[[i]][3],group.B[[i]][4])+rnorm(time.points,sd=sd.var.B))
    # Simulate a horse (i.e. its parameters drawn from the population values)
    # Note: don't forget to require that B_i + I_i > 0 if variable is a sperm response
    #group.A[[i]] <- mvrnorm(1,theta.A,diag(sd.theta.A))
    repeat {
      group.B[[i]][1] <- rnorm(1, B[["mu"]][["mean"]], B[["mu"]][["sd"]])
      if (group.B[[i]][1] >= B[["mu"]][["min"]] & group.B[[i]][1] <= B[["mu"]][["max"]]) break
    }
    repeat {
      group.B[[i]][2] <- rlnorm(1, B[["s1"]][["mean"]], B[["s1"]][["sd"]])
      if (group.B[[i]][2] >= B[["s1"]][["min"]] & group.B[[i]][2] <= B[["s1"]][["max"]]) break
    }
    repeat {
      group.B[[i]][3] <- rlnorm(1, B[["s2"]][["mean"]], B[["s2"]][["sd"]])
      if (group.B[[i]][3] >= B[["s2"]][["min"]] & group.B[[i]][3] <= B[["s2"]][["max"]]) break
    }
    repeat {
      group.B[[i]][4] <- rnorm(1, B[["B"]][["mean"]], B[["B"]][["sd"]])
      if (group.B[[i]][4] >= B[["B"]][["min"]] & group.B[[i]][4] <= B[["B"]][["max"]]) break
    }
    repeat {
      group.B[[i]][5] <- rnorm(1, B[["I"]][["mean"]], B[["I"]][["sd"]])
      if (group.B[[i]][5] >= B[["I"]][["min"]] & group.B[[i]][5] <= B[["I"]][["max"]]) break
    }
    
    # Simulate a response curve for this horse  
    for (j in 1:time.points) {
      if (sym=="FALSE") {
        repeat {
          var.B[j] <- xGauss2(t[j],group.B[[i]][1],group.B[[i]][2],group.B[[i]][3],group.B[[i]][4],group.B[[i]][5]) + rnorm(1, sd=B[["Y"]][["sd"]])
          if (var.B[j] >= B[["Y"]][["min"]] & var.B[j] <= B[["Y"]][["max"]]) break
        }
      }
      if (sym=="TRUE") {
        repeat {
          var.B[j] <- xGauss(t[j],group.B[[i]][1],log(group.B[[i]][2]),group.B[[i]][4],group.B[[i]][5]) + rnorm(1, sd=B[["Y"]][["sd"]])
          if (var.B[j] >= B[["Y"]][["min"]] & var.B[j] <= B[["Y"]][["max"]]) break
        }
      }
    }
    var.B.all <- c(var.B.all,var.B)
  }
  
  # Create data frame
  df <- data.frame(Group=c(rep('A',N*time.points), rep('B',N*time.points)),
                   Subject=s,
                   DPI=rep(t,2*N),
                   Response=c(var.A.all, var.B.all))
  return(list(data=df, A=group.A, B=group.B))
}