# Script to generate case-mother control-mother data with environmental exposure depending on maternal genotype 
# through a linear model and analyze the data with Spmlficmcm, CCMO.dep and CCMO.na

# The results of this simulation are presented in Table 1 of the manuscript
# Methods and software to analyze gene-environment interactions under a case-mother control-mother design with partially missing child genotype
# by Alexandre Bureau, Yuang Tian, Patrick Levallois, Yves Giguere, Jinbo Chen, Hong Zhang

library(doParallel)
RNGkind("L'Ecuyer-CMRG")

# This sets the RNG seed used for the simulations reported in the manuscript
set.seed(149)

# Setting the number of replicates
nrep = 500

# Setting the number of cores for parallel computation
cores <- 20

# Loading CCMO functions
source("https://github.com/yatian20/CCMO.na/blob/main/R/CCMO.na.R?raw=true")
source("https://github.com/yatian20/CCMO.na/blob/main/R/CCMO.dep.R?raw=true")

N <- 2e5
theta <- 0.2

#beta
beta1 <- log(1.8)    #gm
beta2 <- log(1.5)    #gc
beta3 <- log(1)      #no imprinting
beta4 <- log(1.2)    #X
beta5 <- log(1)    #gm X
beta6 <- log(1)    #gc X
beta <- c(0,beta1,beta2,beta3,beta4,beta5,beta6)

n0 <- 1000
n1 <- 1000
n <- n0 + n1
#lambda <- n1 / (n*f) - n0 / (n * (1-f))

    fl <- Y ~ gm + gc + X + X:gm + X:gc

#simulation

fitCMCM = function(...)
{
  library(SPmlficmcm)
  library(Haplin)
  source("genDataRead.R")

#genotype
genotype <- rbind(c(0,0,(1-theta)^2),c(0,1,theta*(1-theta)),c(1,0,theta*(1-theta)),c(1,1,theta^2))
mdip_N <- genotype[sample(1:4, N, replace = TRUE, prob = genotype[,3]),1:2]
gm_N <- apply(mdip_N,1,sum)
# father
pdip_N <- genotype[sample(1:4, N, replace = TRUE, prob = genotype[,3]),1:2]
gp_N <- apply(pdip_N,1,sum)
# child
m_tran <- runif(N) < 0.5
p_tran <- runif(N) < 0.5
gcm_N <- m_tran * mdip_N[,1] + (1 - m_tran) * mdip_N[,2]
gcp_N <- p_tran * pdip_N[,1] + (1 - p_tran) * pdip_N[,2]
gc_N <- gcm_N + gcp_N

#X
e <- rnorm(N)
X_N <- round(log(1.5) *(gm_N - mean(gm_N)) + e)
X_N = pmin(pmax(X_N,-2),2)

# Disease prevalence
f <- 0.01
Z <- cbind(gm_N,gc_N,gcm_N-gcp_N,X_N,gm_N*X_N,gc_N*X_N)
fn <- function(beta0){
  tmp <- exp(beta0 + as.vector(Z %*% beta[-1]))
  p <- tmp/(1+tmp)
  return(mean(p)-f)
}
beta[1] <- uniroot(fn,c(-10,0))$root

#phenotype
tmp <- exp(beta[1] + as.vector(Z %*% beta[-1]))
p <- tmp/(1+tmp)
r <- runif(N)
Y_N <- ifelse(r<p,1,0)

#generate NA
r <- runif(N)
gc_N[r<0.05] <- NA

N1 <- sum(Y_N == 1)
N0 <- N-N1

  case <- sample(1:N1, n1, replace = FALSE)
  control <- sample(1:N0, n0, replace = FALSE)
  gm <- c(gm_N[Y_N == 1][case],gm_N[Y_N == 0][control])
  gc <- c(gc_N[Y_N == 1][case],gc_N[Y_N == 0][control])
  Y <- c(Y_N[Y_N == 1][case],Y_N[Y_N == 0][control])
  X <- c(X_N[Y_N == 1][case],X_N[Y_N == 0][control])
  
  dat=data.frame(Y,X,gm,gc)
  
  fit <- CCMO.na(Y,gm,gc,X,X,X,f,FALSE)
  fit1 <- CCMO.dep(Y,gm,gc,X,X,X,X,f,TRUE)
  fit2 <- try(Spmlficmcm(fl, c(N0,N1), "gm", "gc", dat, 2))
  if (inherits(fit2, "try-error")) fit2.est = fit2.sd = rep(NA,7)
   else 
   {
   	if (is.na(fit2$Value_loglikh)) fit2.est = fit2.sd = rep(NA,7)
   	else 
   	{
   		fit2.est = fit2[["MatR"]][, 1]
   		fit2.sd = fit2[["MatR"]][, 2]
   	}
   }
  # Processing data for Haplin
  gmf = factor(c(NA,gm))
  levels(gmf)=c("AA","AC","CC")
  gcf = factor(c(NA,gc))
  levels(gcf)=c("AA","AC","CC")
  gff = rep(NA,length(gcf))

  dath = data.frame(Y=c(0,Y),X=c(0,X+2),gm=gmf,gf=gff,gc=gcf)
  datf = file()
  write.table(dath,datf,row.names=F,col.names=F,quote=F)

  haplin.data.test <- genDataRead( file.in = datf, file.out = "test",format = "haplin", allele.sep = "", n.vars = 2, cov.header = c("Y","X"),dir.out = tempdir( check = TRUE ),overwrite=T )

  haplin.pre = genDataPreprocess( haplin.data.test, design="cc.triad",overwrite=T)

  haplin.res = haplinStrat( haplin.pre, strata = 2, design = "cc.triad", maternal=T, ccvar = 1, use.missing=T, response="mult", reference=1,max.EM.iter = 100)
  haplin.est = sapply(haplin.res,function(obj)coef(obj$result$result)) 
  haplin.sd = sapply(haplin.res,function(obj) sqrt(diag(obj$var.cov$var.cov)))
  haplin.gxe = gxe(haplin.res)

  return(c(fit$est,fit$sd,fit1$est,fit1$sd,fit2.est,fit2.sd,haplin.est,haplin.sd,fit$est.log,fit$sd.log,haplin.gxe=haplin.gxe$gxe.test[,4]))
}
EST = mclapply(1:nrep,fitCMCM,mc.cores = cores)
CMCM.dep01null.res <- do.call(cbind, EST)
