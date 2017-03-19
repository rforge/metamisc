logit <- function(x) { if(is.numeric(x))  log(x/(1-x)) else stop("x is not numeric!") }

inv.logit <- function(x) {  if(is.numeric(x)) 1/(1+exp(-x)) else stop("x is not numeric!") }

generateMCMCinits <- function(n.chains, model.pars)
{
  inits <- list()
  for (i in 1:n.chains) {
    inits.i <- list()
    for (j in 1:length(model.pars)) {
      parname <- model.pars[[j]]$param
      fprior <- model.pars[[j]]$param.f
      fargs <- model.pars[[j]]$param.args
      inits.i[[parname]] = do.call(fprior, fargs)
    }
    inits[[i]] <- inits.i
  }
  return(inits)
}


#Update SE(c.index) using Method 4 of Newcombe
restore.c.var<- function(cstat, N.subjects, N.events, restore.method="Newcombe.4", model="normal/logit") {
  n <- N.events #Number of events
  m <- N.subjects-N.events #Number of non-events
  
  if (missing(restore.method)) {
    restore.method <- "Newcombe.4"
  }
  
  if (restore.method=="Hanley" | restore.method=="Newcombe.2") {
    mstar <- m-1
    nstar <- n-1
  } else if (restore.method=="Newcombe.4") {
    mstar <- nstar <- N.subjects/2-1
  } else {
    stop ("Method not implemented yet!")
  }
  
  if (model=="normal/logit") {
    out <- (((1+nstar*(1-cstat)/(2-cstat) + mstar*cstat/(1+cstat)))/(m*n*cstat*(1-cstat)))
  } else if (model=="normal/identity") {
    out <- ((cstat*(1-cstat)*(1+nstar*(1-cstat)/(2-cstat) + mstar*cstat/(1+cstat)))/(m*n))
  } else {
    stop ("Meta-analysis model not implemented!")
  }
  
  return(out)
}

restore.oe.var <- function(citl, citl.se, Po) {
  nom <- ((Po-1)**2)*((Po**2)+1)*((exp(Po+citl))**2)*(citl.se**2)
  denom <- (Po*(-exp(citl))+Po+exp(citl))**2
  out <- nom/denom
  return(out)
}

extrapolateOE <- function(Po, Pe, var.Po, t.val, t.ma, N, model="normal/log") {
  Po.new <- 1-exp(t.ma*log(1-Po)/t.val)
  Pe.new <- 1-exp(t.ma*log(1-Pe)/t.val)
  theta <- theta.var <- NA
  
  if (model=="normal/identity") {
    theta <- Po.new/Pe.new
    
    if (missing(var.Po)) {
      # Approximate SE of Kaplan-Meier using binomial distribution
      theta.var <- ((t.ma**2)*exp(2*t.ma*log(1-Po)/t.val)*Po)/((t.val**2)*N*(1-Po)*(1-exp(t.ma*log(1-Pe)/t.val))**2)
    } else {
      # Use error variance of Po, which is equal to the error variance of the KM estimate
      theta.var <- ((t.ma**2)*exp(2*t.ma*log(1-Po)/t.val)*var.Po)/((t.val**2)*((1-Po)**2)*(1-exp(t.ma*log(1-Pe)/t.val))**2)
    }
  } else if (model=="normal/log" | model=="poisson/log") {
    theta <- log(Po.new/Pe.new)
    
    if (missing(var.Po)) {
      # Approximate SE of Kaplan-Meier using binomial distribution
      theta.var <- ((t.ma**2)*Po*exp(2*t.ma*log(1-Po)/t.val))/((t.val**2)*N*(1-Po)*((1-exp(t.ma*log(1-Po)/t.val))**2))
    } else {
      # Use error variance of Po, which is equal to the error variance of the KM estimate
      theta.var <- ((t.ma**2)*exp(2*t.ma*log(1-Po)/t.val)*var.Po)/((t.val**2)*((1-Po)**2)*((1-exp(t.ma*log(1-Po)/t.val))**2))
    }
  } else {
    stop ("Scale not implemented")
  }
  out <- c(theta, theta.var)
  return(out)
}

generateOEdata <- function(O, E, Po, Pe, OE, OE.se, OE.95CI, citl, citl.se, N, model="normal/log") {
  
  # Derive O or E from OE where possible
  O <- ifelse(is.na(O), OE*E, O)
  O <- ifelse(is.na(O), Po*N, O)
  E <- ifelse(is.na(E), O/OE, E)
  E <- ifelse(is.na(E), Pe*N, E)
  Po <- ifelse(is.na(Po), O/N, Po)
  Po <- ifelse(is.na(Po), OE*Pe, Po)
  Pe <- ifelse(is.na(Pe), E/N, Pe)
  Pe <- ifelse(is.na(Pe), Po/OE, Pe)
  
  # Derive
  
  #TODO: allow confidence intervals of OE ratio
  #TODO: allow E/O ratio
  # Apply necessary data transformations
  if (model == "normal/identity") {
    theta <- OE
    theta <- ifelse(is.na(theta), O/E, theta)
    theta <- ifelse(is.na(theta), Po/Pe, theta)
    theta <- ifelse(is.na(theta), -(exp(citl)*(O/N)-exp(citl)-(O/N)), theta) #derive from CITL
    theta.var <- OE.se**2
    theta.cil <- OE.95CI[,1]
    theta.ciu <- OE.95CI[,2]
    theta.var <- ifelse(is.na(theta.var), ((theta.ciu - theta.cil)/(2*qnorm(0.975)))**2, theta.var) #Derive from 95% CI
    theta.var <- ifelse(is.na(theta.var), (((O/N)**2)+1)*((exp(citl))**2)*(citl.se**2), theta.var)
    theta.var <- ifelse(is.na(theta.var), O*(1-Po)/(E**2), theta.var) #BMJ eq 20 (binomial var)
    theta.var <- ifelse(is.na(theta.var), (O/(E**2)), theta.var) #BMJ eq 30 (Poisson var)
    #Extrapolate theta 
    #if (!missing(t.ma) & !missing(t.val)) {
    #  thetaE <- extrapolateOE(Po=Po, Pe=Pe, var.Po=var.Po, t.val=t.val, t.ma=t.ma, N=N, model=model)
    #}
    
    #Check if continuitiy corrections are needed
    Ecc <- E
    Ncc <- N
    Occ <- O
    cc <- which(E==0 & is.na(theta.var))
    Ecc[cc] <- 0.5
    Ncc[cc] <- N[cc]+0.5
    Occ[cc] <- O[cc]+0.5
    theta <- ifelse(theta==Inf, log(Occ/Ecc), theta)
    theta.var <- ifelse(theta.var==Inf, Occ*(1-(Occ/Ncc))/(Ecc**2), theta.var) #BMJ eq 20 (binomial var)
    theta.var <- ifelse(theta.var==Inf, (Occ/(Ecc**2)), theta.var) #BMJ eq 30 (Poisson var)
  } else if (model == "normal/log" | model == "poisson/log") {
    theta <- log(OE)
    theta <- ifelse(is.na(theta), log(O/E), theta)
    theta <- ifelse(is.na(theta), log(Po/Pe), theta)
    theta <- ifelse(is.na(theta), log(-(exp(citl)*(O/N)-exp(citl)-(O/N))), theta) #derive from CITL
    theta.var <- (OE.se/theta)**2
    theta.cil <- log(OE.95CI[,1])
    theta.ciu <- log(OE.95CI[,2])
    theta.var <- ifelse(is.na(theta.var), ((theta.ciu - theta.cil)/(2*qnorm(0.975)))**2, theta.var)
    theta.var <- ifelse(is.na(theta.var),  restore.oe.var(citl=citl, citl.se=citl.se, Po=Po), theta.var)
    theta.var <- ifelse(is.na(theta.var), (1-Po)/O, theta.var) #BMJ eq 27 (binomial var)
    theta.var <- ifelse(is.na(theta.var), (1/O), theta.var) #BMJ eq 36 (Poisson var)
    
    #Check if continuitiy corrections are needed
    Ecc <- E
    Ncc <- N
    Occ <- O
    cc <- which((O==0 | E==0) & is.na(theta.var))
    Ecc[cc] <- 0.5
    Ncc[cc] <- N[cc]+0.5
    Occ[cc] <- O[cc]+0.5
    theta <- ifelse(theta==Inf, log(Occ/Ecc), theta)
    theta.var <- ifelse(theta.var==Inf, (1-(Occ/Ncc))/Occ, theta.var) #BMJ eq 27 (binomial var)
    theta.var <- ifelse(theta.var==Inf, (1/Occ), theta.var) #BMJ eq 36 (Poisson var)
  } else {
    stop(paste("No appropriate meta-analysis model defined: '", model, "'", sep=""))
  }
  
  #Only calculate 95% CI for which no original values were available
  theta.cil[is.na(theta.cil)] <- (theta+qnorm(0.025)*sqrt(theta.var))[is.na(theta.cil)]
  theta.ciu[is.na(theta.ciu)] <- (theta+qnorm(0.975)*sqrt(theta.var))[is.na(theta.ciu)]
  
  ds <- cbind(theta, sqrt(theta.var), theta.cil, theta.ciu, F)
  colnames(ds) <- c("theta", "theta.se", "theta.95CIl", "theta.95CIu", "cont.corr")
  ds[cc,"cont.corr"] <- T
  ds <- as.data.frame(ds)
  return(ds)
}

plotForest <- function(vmasum, xlab, refline, ...) {
  inv.logit <- function(x) {1/(1+exp(-x)) }
  
  if (!is.null(vmasum$rma)) {
    # Forest plot
    if (vmasum$model=="normal/logit") {
      forest(vmasum$rma, transf=inv.logit, xlab=xlab, addcred=T, refline=refline, ...)
    } else if (vmasum$model == "normal/log") {
      forest(vmasum$rma, transf=exp, xlab=xlab, addcred=T, refline=refline,...)
    } else if (vmasum$model=="normal/identity") {
      forest(vmasum$rma, transf=NULL, xlab=xlab, addcred=T, refline=refline, ...)
    } else {
      stop ("Invalid meta-analysis model!")
    }
  } else {
    col <- c("black", "gray50")
    border <- "black"
    lty <- c("solid", "dotted", "solid")
    cex <- 0.8
    efac <- 1
    efac <- rep(efac, 2)
    
    par.usr <- par("usr")
    height <- par.usr[4] - par.usr[3]
    
    k <- dim(vmasum$data)[1]
    slab <- c(as.character(vmasum$slab), "RE Model")
    yi <- vmasum$data[,"theta"]
    ci.lb <- vmasum$data[,"theta.95CIl"]
    ci.ub <- vmasum$data[,"theta.95CIu"]
    
    if (vmasum$model=="normal/logit") {
      xlim <- c(-0.5, 1.5)
      yi <- sapply(yi, inv.logit)
      ci.lb <- sapply(ci.lb, inv.logit)
      ci.ub <- sapply(ci.ub, inv.logit)
      refline <- NA
    } else if (vmasum$model == "normal/log" | vmasum$model == "poisson/log") {
      yi <- sapply(yi, exp)
      ci.lb <- sapply(ci.lb, exp)
      ci.ub <- sapply(ci.ub, exp)
      refline <- 1
      
      plot.multp.l <- 1.2
      plot.multp.r <- 1.2
      rng <- max(ci.ub, na.rm = TRUE) - min(ci.lb, na.rm = TRUE)
      xlim <- c(min(ci.lb, na.rm = TRUE) - rng * plot.multp.l, 
                max(ci.ub, na.rm = TRUE) + rng * plot.multp.r)
      xlim <- round(xlim, 2)
    }
    
    ylim <- c(-1.5, k + 3)
    
    
    
    #Add the meta-analysis summary to the results
    #Note that no transormations are needed here, as summaries are always presented on original scale
    b <- vmasum$results["estimate"]
    yi <- c(yi, b)
    b.ci.lb <- vmasum$results["95CIl"]
    b.ci.ub <- vmasum$results["95CIu"]
    b.cr.lb <- vmasum$results["95PIl"]
    b.cr.ub <- vmasum$results["95PIu"]
    ci.lb <- c(ci.lb, b.ci.lb)
    ci.ub <- c(ci.ub, b.ci.ub)
    
    
    rows <- c(seq(k,1),-1)
    
    annotext <- round(cbind(yi, ci.lb, ci.ub), 2)
    annotext <- matrix(apply(annotext, 2, format, nsmall = 2), ncol = 3)
    annotext <- paste(annotext[,1], "[", annotext[,2], ",", annotext[,3], "]")
    
    
    par.mar <- par("mar")
    par.mar.adj <- par.mar - c(0, 3, 1, 1)
    par.mar.adj[par.mar.adj < 0] <- 0
    par(mar = par.mar.adj)
    on.exit(par(mar = par.mar))
    
    plot(NA, NA, xlim=xlim, ylim=ylim, ylab="", xlab=xlab, yaxt = "n", xaxt = "n", xaxs = "i", bty = "n", ...)
    
    par.usr <- par("usr")
    height <- par.usr[4] - par.usr[3]
    lheight <- strheight("O")
    cex.adj <- ifelse(k * lheight > height * 0.8, height/(1.25 * k * lheight), 1)
    cex <- par("cex") * cex.adj
    
    for (i in 1:k) {
      points(yi[i], rows[i], pch = 15, ...)
      
      segments(ci.lb[i], rows[i], ci.ub[i], rows[i], ...)
      
      segments(ci.lb[i], rows[i] - (height/150) * cex * 
                 efac[1], ci.lb[i], rows[i] + (height/150) * cex * 
                 efac[1], ...)
      
      segments(ci.ub[i], rows[i] - (height/150) * cex * 
                 efac[1], ci.ub[i], rows[i] + (height/150) * cex * 
                 efac[1], ...)
    }
    text(xlim[1], rows, slab, pos = 4, cex = cex, ...)
    text(x = xlim[2], rows, labels = annotext, pos = 2, cex = cex, ...)
    
    # Add prediction interval
    segments(b.cr.lb, -1, b.cr.ub, -1, lty = lty[2], col = col[2], ...)
    segments(b.cr.lb, -1 - (height/150) * cex * efac[1], 
             b.cr.lb, -1 + (height/150) * cex * efac[1], 
             col = col[2], ...)
    segments(b.cr.ub, -1 - (height/150) * cex * efac[1], 
             b.cr.ub, -1 + (height/150) * cex * efac[1], 
             col = col[2], ...)
    
    # Add diamond for summary estimate
    polygon(x = c(b.ci.lb, b, b.ci.ub, b), y = c(-1, -1 + 
                                                   (height/100) * cex * efac[2], -1, -1 - (height/100) * 
                                                   cex * efac[2]), col = col[1], border = border, ...)
    
    # Add refline
    if (is.numeric(refline)) {
      segments(refline, ylim[1] - 5, refline, ylim[2] - 2, lty = "dotted", ...)
    }
    
    # Add separation line between forest plot and meta-analysis results
    abline(h = 0, lty = 1, ...)
    
    # Add separation line on top of the figure
    abline(h = ylim[2] - 2, lty = lty[3], ...)
    
    if (vmasum$model=="normal/logit") {
      axis(side = 1, at = c(0,0.2,0.4,0.6,0.8,1), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis = 1, ...)
    } else if (vmasum$model == "normal/log" | vmasum$model == "poisson/log") {
      axis(side = 1, at = c(0:ceiling(max(exp(vmasum$data[,"theta.95CIu"]), na.rm=T))), labels = c(0:ceiling(max(exp(vmasum$data[,"theta.95CIu"]), na.rm=T))), cex.axis = 1, ...)
    }
    
  }

}
