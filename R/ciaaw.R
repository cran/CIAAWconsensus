######################################################################

## ciaaw.R
##
## AUTHOR:            Juris Meija & Antonio Possolo
## MODIFICATION DATE: 2018 Sep 8

######################################################################

at.weight = function (ratio, ratio.cov, element, ref.isotope, data=NULL)
{
    ratio.cov <- as.matrix(ratio.cov)
    
    if(is.null(data)){
      e <- new.env(parent = emptyenv())
      data(list='ciaaw.mass.2016', envir = e)
      data <- e$ciaaw.mass.2016
    }

    mass.raw = data[which(data$element==element),]
    isotope.list = as.vector(mass.raw$isotope)
    mass.x = as.numeric(mass.raw$mass)
    mass.u = as.numeric(mass.raw$uncertainty)

    ref = match(ref.isotope, isotope.list)
    p = length(isotope.list) - 1
    max.iter = 1e6
    mass.sequence = c(ref, seq(from = 1, to = (1 + p))[-ref])
    ratios.B = abs(cbind(1, mvtnorm::rmvnorm(max.iter, mean=ratio, sigma=ratio.cov, method = "svd")))
    mass.B = mvtnorm::rmvnorm(max.iter, mean=mass.x, sigma=diag(mass.u^2), method="svd")[,mass.sequence]
    atomicWeight = rowSums(ratios.B * mass.B)/rowSums(ratios.B)
    abundances = ratios.B/rowSums(ratios.B)
    colnames(abundances) = mass.x[mass.sequence]
    aw.U95 = diff(quantile(atomicWeight, probs = c(0.025, 0.975)))/2
    names(aw.U95) = NULL
    abundances.U95 = apply(abundances, 2, function(x) diff(quantile(x, probs = c(0.025, 0.975)))/2)
    aw.digits = round(2 - log10(sd(atomicWeight)), 0)
    x.digits = round(2 - log10(min(apply(abundances, 2, sd))), 0)
    x.mean = colMeans(abundances)
    x.u = sqrt(diag(cov(abundances)))
    cat("\n")
    cat("Atomic weight", formatC(mean(atomicWeight), format = "f", digits = aw.digits), "\n")
    cat("            u", formatC(sd(atomicWeight), format = "f", digits = aw.digits), "\n")
    cat("        U 95%", formatC(aw.U95, format = "f", digits = aw.digits), "\n\n")
    cat("Isotopic abundances", "\n")
    cat("  mass/Da", colnames(abundances), "\n")
    cat("abundance", formatC(x.mean, format = "f", digits = x.digits), "\n")
    cat("        u", formatC(x.u, format = "f", digits = x.digits), "\n")
    cat("    U 95%", formatC(abundances.U95, format = "f", digits = x.digits), "\n")
    result = list(aw = mean(atomicWeight), aw.u = sd(atomicWeight),  
        aw.U95 = aw.U95, abundances = x.mean,
        abundances.cov = cov(abundances),  
        abundances.u = x.u, abundances.U95 = abundances.U95)
}

## Atomic weight and isotopic abundances of iridium which correspond
## to the isotope ratio 191Ir/193Ir = 0.59471(13)
# at.weight(0.59471, 0.00013^2, "iridium", "193Ir")

## Atomic weight and isotopic abundances of silicon which correspond
## to isotope ratios 28Si/29Si = 1.074(69) and 30Si/29Si = 260(11)
## with a correlation of 0.80 between the two isotope ratios
# ratios = c(1.074,260)
# sigma = matrix(c(0.069^2,0.80*0.069*11,0.80*0.069*11,11^2),ncol=2,byrow=TRUE)
# at.weight(ratios, sigma, "silicon", "29Si")

normalize.ratios = function (dat, element, ref.isotope, expand = FALSE) 
{
  e <- new.env(parent = emptyenv())
  data(list='ciaaw.mass.2016', envir = e)
  data <- e$ciaaw.mass.2016
  
  mass.raw = data[which(data$element==element),]
  isotope.list = as.vector(mass.raw$isotope)
  mass.x = as.numeric(mass.raw$mass)
  mass.u = as.numeric(mass.raw$uncertainty)
  
    ## helper functions from miscTools
    insertRow <- function (m, r, v = NA, rName = "") {
      if (r == as.integer(r)) { r <- as.integer(r) }
      else { stop("argument 'r' must be an integer") }
      nr <- nrow(m)
      nc <- ncol(m)
      rNames <- rownames(m)
      if (is.null(rNames) & rName != "") { rNames <- rep("", nr) }
      if (r == 1) { m2 <- rbind(matrix(v, ncol = nc), m)
        if (!is.null(rNames)) { rownames(m2) <- c(rName, rNames) } }
      else if (r == nr + 1) { m2 <- rbind(m, matrix(v, ncol = nc))
        if (!is.null(rNames)) { rownames(m2) <- c(rNames, rName) } }
      else { m2 <- rbind(m[1:(r - 1), , drop = FALSE], matrix(v, ncol = nc), m[r:nr, , drop = FALSE])
        if (!is.null(rNames)) { rownames(m2) <- c(rNames[1:(r - 1)], rName, rNames[r:nr]) } }
      return(m2)
    }
    
    insertCol <- function (m, c, v = NA, cName = "") {
      if (c == as.integer(c)) {c <- as.integer(c)}
      else {stop("argument 'c' must be an integer")}
      nr <- nrow(m)
      nc <- ncol(m)
      cNames <- colnames(m)
      if (is.null(cNames) & cName != "") { cNames <- rep("", nc) }
      if (c == 1) { m2 <- cbind(matrix(v, nrow = nr), m)
        if (!is.null(cNames)) { colnames(m2) <- c(cName, cNames) } }
      else if (c == nc + 1) { m2 <- cbind(m, matrix(v, nrow = nr))
        if (!is.null(cNames)) {colnames(m2) <- c(cNames, cName) } }
      else { m2 <- cbind(m[, 1:(c - 1), drop = FALSE], matrix(v,nrow = nr), m[, c:nc, drop = FALSE])
        if (!is.null(cNames)) { colnames(m2) <- c(cNames[1:(c - 1)], cName, cNames[c:nc]) } }
      return(m2)
    }
    
    ###################

    isotopes = matrix(unlist(strsplit(as.character(dat$Outcome), "\\/")), ncol = 2, byrow = TRUE)
    dat$nominator = isotopes[, 1]
    dat$denominator = isotopes[, 2]
    ref = match(ref.isotope, isotope.list)
    f.J.r.norm = function(x) x/x[ref]
    denominator.list = lapply(split(dat$denominator, dat$Study), 
        function(x) match(unique(x), isotope.list))
    r = mapply(function(x, y) append(x, 1, after = y - 1), split(dat$Value, 
        dat$Study), denominator.list, SIMPLIFY = FALSE)
    cov = mapply(function(x, y)
        insertCol(insertRow(diag(x, nrow = length(x)), y, 0), y, 0)^2, 
        split(dat$Unc, dat$Study), 
        denominator.list, SIMPLIFY = FALSE)
    J.r.norm = mapply(function(x) numDeriv::jacobian(f.J.r.norm, 
        x), r, SIMPLIFY = FALSE)
    cov.norm = mapply(function(x, y) x %*% y %*% t(x), J.r.norm, 
        cov, SIMPLIFY = FALSE)
    r.norm = mapply(function(x, y) x/x[y], r, ref, SIMPLIFY = FALSE)
    n.studies = length(r)
    p = length(r[[1]]) - 1
    if (expand) {
        cov.norm.red = mapply(function(x, y)
            as.matrix(x^2 * y[-ref, -ref]), 
            lapply(split(dat$k_extra, dat$Study), unique), 
            cov.norm, SIMPLIFY = FALSE)
    } else {
        cov.norm.red = mapply(function(y) as.matrix(y[-ref, -ref]), 
                              cov.norm, SIMPLIFY = FALSE) }
    cov.norm.red1 = lapply(lapply(cov.norm.red, diag), function(x) diag(x, 
        nrow = length(x)))
    S = cov.norm.red1
    R = matrix(do.call(rbind, lapply(r.norm, matrix, ncol = p + 
        1, byrow = TRUE))[, -ref], nrow = n.studies)
    Svars = do.call(rbind, lapply(lapply(S, diag), matrix, ncol = p, 
        byrow = TRUE))
    colnames(R) = colnames(Svars) = paste(isotope.list[-ref], 
        "/", ref.isotope, sep = "")
    rownames(R) = rownames(Svars) = unlist(lapply(split(paste(dat$Year, 
        "-", dat$Author, sep = ""), dat$Study), head, 1))
    names(S) = rownames(Svars)
    results = list(R = R, u.R = sqrt(Svars), cov.R = S)
    list(R = R, u.R = sqrt(Svars))
}

## Normalize all iridium isotope data to iridium-193
# normalize.ratios(iridium.data, "iridium", "193Ir")

mmm <- function(y, uy, knha=TRUE, verbose=TRUE) {
  
  y.l <- split(y, col(y))
  var.l <- split(uy^2, col(uy))
  w.l <- split(1/uy^2, col(uy))
  
  n = ncol(y) # no missing outcomes
  n.studies = nrow(y)
  
  y.w = mapply(function(w,y) sum(w*y)/sum(w), w.l, y.l)
  Q = mapply(function(w,y,mu) sum(w*(y-mu)^2), w.l, y.l, y.w)
  tau2 = mapply(function(w,q) max(0, (q - (n.studies-1))/(sum(w)-sum(w^2)/sum(w))), w.l, Q)
  w.s = 1/(uy^2 + matrix(rep(tau2, n.studies), ncol=n, byrow=TRUE))
  
  # variance of the consensus values
  var.beta = 1/colSums(w.s)
  
  # consensus values
  beta = mapply(function(ws, y) sum(ws*y)/sum(ws), split(w.s,col(w.s)), y.l, SIMPLIFY=TRUE)
  
  # covariance estimate
  cov.beta = matrix(0, ncol=n, nrow=n)
  if (n!=1) {
    for (i in 1:n){
      for (j in 1:n){
        if(i!=j) cov.beta[i,j] <- sum(w.s[,i]*w.s[,j]*(y[,i]-beta[i])*(y[,j]-beta[j])/(sum(w.s[,i])*sum(w.s[,j])))
      }}
    diag(cov.beta) <- var.beta
  } else { cov.beta = as.matrix(unlist(var.beta)) }
  
  # spectral decomposition of cov.beta to avoid possible negative variances
  tol=-1e-18
  if (min(eigen(cov.beta)$values)<tol){
    U = eigen(cov.beta)$vectors
    Dpsd = diag(as.vector(eigen(cov.beta)$values), nrow=n)
    Dpsd[which(Dpsd<0)] = 0
    cov.beta = U %*% Dpsd %*% t(U) }
  
  #Knapp-Hartung adjustment for the variances and covariances of each outcome (Birge ratio)
  if(knha==FALSE) H2 <- rep(1, n)
  if(knha==TRUE) {H2 <- mapply(function(w,y,mu) sum(w*(y-mu)^2), split(w.s,col(w.s)), y.l, as.list(beta))/(n.studies-1); H2[which(H2 < 1)] <- 1}
  
  # expand uncertainties by H
  cov.betaH <- outer(sqrt(H2),sqrt(H2)) * cov.beta
  
  # H-adjusted standard uncertainties 
  u.betaH = sqrt(diag(cov.betaH))
  
  # correlation  matrix
  cor.betaH <- if(n==1) as.matrix(1) else {stats::cov2cor(cov.betaH)} 
  
  # 95% confidence intervals for isotope ratios based on t-distribution
  mt = mvtnorm::qmvt(0.975, corr=cor.betaH, df=(n.studies-1)*n, type="shifted")$quantile
  ci.betaH = mt*u.betaH
  
  # Relative total variability due to heterogeneity (in percent)
  I2 = round(100*(Q - (n.studies-1))/Q, 1)
  I2[which(I2 < 0)] <- 0
  
  # determine the number of decimal digits for printout
  d = round(-log10(min(u.betaH))+2,0)
  
  names(beta) <- names(u.betaH) <- names(ci.betaH) <- names(H2) <- names(I2) <- colnames(y)
  colnames(cov.betaH) <- rownames(cov.betaH) <- colnames(y)
  colnames(cor.betaH) <- rownames(cor.betaH) <- colnames(y)
  
  # annotated summary output
  if (verbose ==TRUE) {cat("Multivariate Marginal Method of Moments with working independence assumption", "\n")
    cat("Chen et al (2016) DOI: 10.1002/sim.6789", "\n", "\n")
    cat("Number of studies: ", n.studies, "\n")
    cat("Isotope ratios", "\n", names(beta), "\n", "\n")
    cat("Consensus values", "\n", formatC(unlist(beta), format='f', digits=d), "\n", "\n")
    cat("Standard uncertainty", "\n", formatC(u.betaH, format='f', digits=d), "\n", "\n")
    cat("Expanded uncertainty, 95 %", "\n", formatC(ci.betaH, format='f', digits=d), "\n", "\n")
    cat("Birge ratio (Knapp-Hartung adjustment) for each outcome", "\n", formatC(sqrt(H2), format='f', digits=2), "\n", "\n")
    cat("Relative total variability due to heterogeneity (%) for each outcome", "\n", formatC(I2, format='f', digits=1), "\n", "\n")
  }
  
  # gather the results
  result = list( studies=n.studies,
                 beta=unlist(beta), 
                 beta.u = u.betaH, 
                 beta.U95 = ci.betaH, 
                 beta.cov=cov.betaH, 
                 beta.cor = cor.betaH, 
                 I2 = I2, 
                 H = sqrt(H2), 
                 df=n*(n.studies-1))
}

## Consensus isotope amount ratios for platinum
# df = normalize.ratios(platinum.data, "platinum", "195Pt")
# mmm(df$R, df$u.R)

## A function to convert a given set of isotopic abundances and their uncertainties to the corresponding isotope ratios
abundances2ratios <- function(x,ux,ref=1,iterations=1e4){
  n.isotopes = length(x) # number of isotopes
  cov.R = list()
  u.R = list()
  it.n = 0
  
  # setup
  cor = diag(1,nrow=n.isotopes) #identity
  # assign N random correlations (N = n.isotopes(n.isotopes-3)/2)
  z=upper.tri(cor)&!cor[1,]
  z[,which(apply(z,2,sum)==1)]<-FALSE
  
  cor[z] <- runif(n.isotopes*(n.isotopes-3)/2,min=-1,max=+1)
  cor[lower.tri(cor)] <- t(cor)[lower.tri(cor)]
  cov = cor*outer(ux,ux)
  
  # find columns with a single zero entry (and the corresponding rows)
  cov1<-cov
  # find columns with a single zero entry (and the corresponding rows)
  col1=(1:n.isotopes)[apply(cov1==0,2,sum)==1]
  row1=apply(as.matrix(col1),1,function(x) which(cov1[,x]==0))
  for(i in 1:length(col1)) cov1[row1[i],col1[i]] <- -sum(cov1[,col1[i]])
  cov1[lower.tri(cov1)] <- t(cov1)[lower.tri(cov1)]
  
  cov2<-cov1
  # find missing covariances and assign a variable symbol to them
  for(i in 1:n.isotopes){
    for(j in 1:n.isotopes){ifelse(cov1[i,j]==0,cov2[i,j]<-paste("a",min(i,j),max(i,j),sep=""),cov2[i,j]<-"") }
  }
  variables=unique(c(cov2))[-1]
  d<-matrix(0,ncol=length(variables),nrow=length(variables))
  
  # rows in the design matrix for zero-column and zero-row sum
  for(i in 1:length(variables)){
    d[i,]<-as.numeric(1:length(variables) %in% stats::na.exclude(match(cov2[,i], variables)))
  }
  
  # Monte Carlo
  for(i in 1:iterations){
    cor = diag(1,nrow=n.isotopes) #identity
    cor[z] <- runif(n.isotopes*(n.isotopes-3)/2,min=-1,max=+1)
    cor[lower.tri(cor)] <- t(cor)[lower.tri(cor)]
    cov1 = cor*outer(ux,ux)
    
    for(i in 1:length(col1)) cov1[row1[i],col1[i]] <- -sum(cov1[,col1[i]])
    cov1[lower.tri(cov1)] <- t(cov1)[lower.tri(cov1)]
    
    # the corresponding entries in the response vector
    y=matrix(-apply(cov1,2,sum))[!apply(cov2=="",2,all),]
    
    # find missing covariances by Gauss elimination
    coef<-solve(d)%*%y
    for (i in 1:length(variables)){ cov1[cov2==variables[i]]<-coef[i] }
    
    # zero column-sum and zero row-sum covariance matrix: cov1
    cor <- stats::cov2cor(cov1)
    
    # rejection of isotopic abundance correlations that are outside [-1,+1]
    if(max(abs(cor))<=1){
      # jacobian to get isotope ratio uncertainties
      J <- diag(x[ref]^(-1), nrow=n.isotopes)
      J[,ref] <- -x*x[ref]^(-2)
      J[ref,ref] <- 0
      cov.Ratios <- J %*% cov1 %*% t(J) # isotope ratio covariance matrix
      # further rejection of isotope ratio correlations that are outside [-1,+1]
      if (max(abs(stats::cov2cor(cov.Ratios[-ref,-ref])))<=1){
        it.n = it.n+1
        cov.R[[it.n]] <- cov.Ratios
        u.R[[it.n]] <- sqrt(diag(cov.Ratios))
      }
      else{next}
    }
    else{next}
  }
  
  if(it.n==0) stop("The provided isotopic abundance uncertainties are inconsistent")
  # RESULTS
  # find the average isotope ratio covariance matrix
  cov.R.result <- cov.R[[which.max(unlist(lapply(u.R, mean, na.rm = FALSE)))]]
  R <- x[-ref]/x[ref] # output ratios
  uR <- sqrt(diag(cov.R.result[-ref,-ref])) # uncertainties
  list(R=R, R.u=uR, R.cov = cov.R.result[-ref,-ref], N=it.n)
}
## Example: Convert the isotopic abundances of zinc to the corresponding isotope ratios
#  x = c(0.48630, 0.27900, 0.04100, 0.18750, 0.00620)
# ux = c(0.00091, 0.00076, 0.00031, 0.00135, 0.00010)
# z=abundances2ratios(x, ux, ref=2)
# at.weight(z$R, z$R.cov, "zinc", "66Zn")