# Script to compute the statistics in Table 1 of the manuscript
# Methods and software to analyze gene-environment interactions under a case-mother control-mother design with partially missing child genotype
# by Alexandre Bureau, Yuang Tian, Patrick Levallois, Yves Giguere, Jinbo Chen, Hong Zhang

# It assumes the objects CMCM.dep01.res and CMCM.indep01.res are already present in the R environment.
# These objects are created when executing the scripts Data_Generation_dep.R and Data_Generation_indep.R respectively

# Coefficients which are estimated by prospective and retrospective likelihood
beta.vec = beta[-4]
# Coefficients which are estimated by Haplin
hbeta.vec = rep(beta.vec[2:3],5) + rep(-2:2,rep(2,5))*beta.vec[5:6]

# Indices of beta estimates
bind = c(1:6,14:19,28:33,102:107)
# Indices of standard error estimates
seind = c(8:13,22:27,35:40,108:113)
# Indices of Haplin estimates and standard error estimates
hind = c(51:50,56:55,61:60,66:65,71:70)
hseind = hind + 30

# Dependence case

#MCAR 5%

# Proportion of missing beta before removing outliers
missing.beta.dep01 = apply(CMCM.dep01.res[bind,],1,function(vec) mean(is.na(vec)))

    ## Identifying outlying estimates
    par.mad <- apply(CMCM.dep01.res[bind,], 1, mad, na.rm=T)
    par.med <- apply(CMCM.dep01.res[bind,], 1, median, na.rm=T) 
    outliers.mat <- CMCM.dep01.res[bind,] < par.med - 3*sqrt(par.mad) | CMCM.dep01.res[bind,] > par.med + 3*sqrt(par.mad)
    apply(outliers.mat, 1, sum, na.rm=T)
    # Outliers are present only for prospective likelihood
    outliers.vec <- apply(outliers.mat, 2, any, na.rm=T)
    mean(outliers.vec)
    # We keep outlying estimates and their SE

bias.dep01.vec = apply(CMCM.dep01.res[bind,],1,mean,na.rm=T) - beta.vec
emp_se.dep01.vec = apply(CMCM.dep01.res[bind,],1,sd,na.rm=T)
mean_se.dep01.vec = apply(CMCM.dep01.res[seind,],1,mean,na.rm=T)
lower.bound = CMCM.dep01.res[bind,] - qnorm(0.975)*CMCM.dep01.res[seind,]
upper.bound = CMCM.dep01.res[bind,] + qnorm(0.975)*CMCM.dep01.res[seind,]
coverage.dep01.vec = apply(lower.bound<beta.vec & upper.bound>beta.vec,1,mean,na.rm=T)
power.dep01.vec = apply(lower.bound[c(5:6,11:12,17:18,23:24),] > 0,1,mean,na.rm=T)

# Haplin
bias.dep01.haplin.vec = apply(CMCM.dep01.res[hind,],1,mean,na.rm=T) - hbeta.vec
emp_se.dep01.haplin.vec = apply(CMCM.dep01.res[hind,],1,sd,na.rm=T)
mean_se.dep01.haplin.vec = apply(CMCM.dep01.res[hseind,],1,mean,na.rm=T)
lower.bound = CMCM.dep01.res[hind,] - qnorm(0.975)*CMCM.dep01.res[hseind,]
upper.bound = CMCM.dep01.res[hind,] + qnorm(0.975)*CMCM.dep01.res[hseind,]
coverage.dep01.haplin.vec = apply(lower.bound<hbeta.vec & upper.bound>hbeta.vec,1,mean,na.rm=T)
power.dep01.haplin.vec = apply(CMCM.dep01.res[nrow(CMCM.dep01.res)-0:1,]<0.05,1,mean,na.rm=T)

#MCAR 20%

# Proportion of missing beta before removing outliers
missing.beta.dep02 = apply(CMCM.dep02.res[bind,],1,function(vec) mean(is.na(vec)))

    ## Identifying outlying estimates
    par.mad <- apply(CMCM.dep02.res[bind,], 1, mad, na.rm=T)
    par.med <- apply(CMCM.dep02.res[bind,], 1, median, na.rm=T) 
    outliers.mat <- CMCM.dep02.res[bind,] < par.med - 3*sqrt(par.mad) | CMCM.dep02.res[bind,] > par.med + 3*sqrt(par.mad)
    apply(outliers.mat, 1, sum, na.rm=T)
    # Outliers are present only for prospective likelihood
    outliers.vec <- apply(outliers.mat, 2, any, na.rm=T)
    mean(outliers.vec)
    # We keep outlying estimates and their SE

bias.dep02.vec = apply(CMCM.dep02.res[bind,],1,mean,na.rm=T) - beta.vec
emp_se.dep02.vec = apply(CMCM.dep02.res[bind,],1,sd,na.rm=T)
mean_se.dep02.vec = apply(CMCM.dep02.res[seind,],1,mean,na.rm=T)
lower.bound = CMCM.dep02.res[bind,] - qnorm(0.975)*CMCM.dep02.res[seind,]
upper.bound = CMCM.dep02.res[bind,] + qnorm(0.975)*CMCM.dep02.res[seind,]
coverage.dep02.vec = apply(lower.bound<beta.vec & upper.bound>beta.vec,1,mean,na.rm=T)
power.dep02.vec = apply(lower.bound[c(5:6,11:12,17:18,23:24),] > 0,1,mean,na.rm=T)

# Haplin
bias.dep02.haplin.vec = apply(CMCM.dep02.res[hind,],1,mean,na.rm=T) - hbeta.vec
emp_se.dep02.haplin.vec = apply(CMCM.dep02.res[hind,],1,sd,na.rm=T)
mean_se.dep02.haplin.vec = apply(CMCM.dep02.res[hseind,],1,mean,na.rm=T)
lower.bound = CMCM.dep02.res[hind,] - qnorm(0.975)*CMCM.dep02.res[hseind,]
upper.bound = CMCM.dep02.res[hind,] + qnorm(0.975)*CMCM.dep02.res[hseind,]
coverage.dep02.haplin.vec = apply(lower.bound<hbeta.vec & upper.bound>hbeta.vec,1,mean,na.rm=T)
power.dep02.haplin.vec = apply(CMCM.dep02.res[nrow(CMCM.dep02.res)-0:1,]<0.05,1,mean,na.rm=T)

# MNAR 20%
# Proportion of missing beta before removing outliers
missing.beta.dep02MNAR = apply(CMCM.dep02MNAR.res[bind,],1,function(vec) mean(is.na(vec)))

## Identifying outlying estimates
par.mad <- apply(CMCM.dep02MNAR.res[bind,], 1, mad, na.rm=T)
par.med <- apply(CMCM.dep02MNAR.res[bind,], 1, median, na.rm=T) 
outliers.mat <- CMCM.dep02MNAR.res[bind,] < par.med - 3*sqrt(par.mad) | CMCM.dep02MNAR.res[bind,] > par.med + 3*sqrt(par.mad)
apply(outliers.mat, 1, sum, na.rm=T)
# Outliers are present only for prospective likelihood
outliers.vec <- apply(outliers.mat, 2, any, na.rm=T)
mean(outliers.vec)
# We keep outlying estimates and their SE

bias.dep02MNAR.vec = apply(CMCM.dep02MNAR.res[bind,],1,mean,na.rm=T) - beta.vec
emp_se.dep02MNAR.vec = apply(CMCM.dep02MNAR.res[bind,],1,sd,na.rm=T)
mean_se.dep02MNAR.vec = apply(CMCM.dep02MNAR.res[seind,],1,mean,na.rm=T)
lower.bound = CMCM.dep02MNAR.res[bind,] - qnorm(0.975)*CMCM.dep02MNAR.res[seind,]
upper.bound = CMCM.dep02MNAR.res[bind,] + qnorm(0.975)*CMCM.dep02MNAR.res[seind,]
coverage.dep02MNAR.vec = apply(lower.bound<beta.vec & upper.bound>beta.vec,1,mean,na.rm=T)
power.dep02MNAR.vec = apply(lower.bound[c(5:6,11:12,17:18,23:24),] > 0,1,mean,na.rm=T)

# Haplin
bias.dep02MNAR.haplin.vec = apply(CMCM.dep02MNAR.res[hind,],1,mean,na.rm=T) - hbeta.vec
emp_se.dep02MNAR.haplin.vec = apply(CMCM.dep02MNAR.res[hind,],1,sd,na.rm=T)
mean_se.dep02MNAR.haplin.vec = apply(CMCM.dep02MNAR.res[hseind,],1,mean,na.rm=T)
lower.bound = CMCM.dep02MNAR.res[hind,] - qnorm(0.975)*CMCM.dep02MNAR.res[hseind,]
upper.bound = CMCM.dep02MNAR.res[hind,] + qnorm(0.975)*CMCM.dep02MNAR.res[hseind,]
coverage.dep02MNAR.haplin.vec = apply(lower.bound<hbeta.vec & upper.bound>hbeta.vec,1,mean,na.rm=T)
power.dep02MNAR.haplin.vec = apply(CMCM.dep02MNAR.res[nrow(CMCM.dep02MNAR.res)-1:0,]<0.05,1,mean,na.rm=T)

# Independence case

#MCAR 5%

# Proportion of missing beta before removing outliers
missing.beta.indep01 = apply(CMCM.indep01.res[bind,],1,function(vec) mean(is.na(vec)))

    ## Identifying outlying estimates
    par.mad <- apply(CMCM.indep01.res[bind,], 1, mad, na.rm=T)
    par.med <- apply(CMCM.indep01.res[bind,], 1, median, na.rm=T) 
    outliers.mat <- CMCM.indep01.res[bind,] < par.med - 3*sqrt(par.mad) | CMCM.indep01.res[bind,] > par.med + 3*sqrt(par.mad)
    apply(outliers.mat, 1, sum, na.rm=T)
    # Outliers are present only for prospective likelihood
    outliers.vec <- apply(outliers.mat, 2, any, na.rm=T)
    mean(outliers.vec)
    # No outlying estimates 

bias.indep01.vec = apply(CMCM.indep01.res[bind,],1,mean,na.rm=T) - beta.vec
emp_se.indep01.vec = apply(CMCM.indep01.res[bind,],1,sd,na.rm=T)
mean_se.indep01.vec = apply(CMCM.indep01.res[seind,],1,mean,na.rm=T)
lower.bound = CMCM.indep01.res[bind,] - qnorm(0.975)*CMCM.indep01.res[seind,]
upper.bound = CMCM.indep01.res[bind,] + qnorm(0.975)*CMCM.indep01.res[seind,]
coverage.indep01.vec = apply(lower.bound<beta.vec & upper.bound>beta.vec,1,mean,na.rm=T)
power.indep01.vec = apply(lower.bound[c(5:6,11:12,17:18,23:24),] > 0,1,mean,na.rm=T)

# Haplin
bias.indep01.haplin.vec = apply(CMCM.indep01.res[hind,],1,mean,na.rm=T) - hbeta.vec
emp_se.indep01.haplin.vec = apply(CMCM.indep01.res[hind,],1,sd,na.rm=T)
mean_se.indep01.haplin.vec = apply(CMCM.indep01.res[hseind,],1,mean,na.rm=T)
lower.bound = CMCM.indep01.res[hind,] - qnorm(0.975)*CMCM.indep01.res[hseind,]
upper.bound = CMCM.indep01.res[hind,] + qnorm(0.975)*CMCM.indep01.res[hseind,]
coverage.indep01.haplin.vec = apply(lower.bound<hbeta.vec & upper.bound>hbeta.vec,1,mean,na.rm=T)
power.indep01.haplin.vec = apply(CMCM.indep01.res[nrow(CMCM.indep01.res)-0:1,]<0.05,1,mean,na.rm=T)


# Gathering the results in a matrix and writing to a file in a format to facilitate creating Table 1
res = round(cbind(c(bias.indep01.vec,bias.indep01.haplin.vec),c(emp_se.indep01.vec,emp_se.indep01.haplin.vec),c(mean_se.indep01.vec,mean_se.indep01.haplin.vec),c(coverage.indep01.vec,coverage.indep01.haplin.vec),c(bias.dep01.vec,bias.dep01.haplin.vec),c(emp_se.dep01.vec,emp_se.dep01.haplin.vec),c(mean_se.dep01.vec,mean_se.dep01.haplin.vec),c(coverage.dep01.vec,coverage.dep01.haplin.vec)),3)
# We remove results for the intercept
res.sub = cbind(c(rep(c("$\\beta_{G^M}$","$\\beta_{G^C}$","$\\beta_{X}$","$\\beta_{X x G^M}$","$\\beta_{X x G^C}$"),4),"$\\beta_{G^C} -2 \\beta_{X x G^M}$","$\\beta_{G^M} -2 \\beta_{X x G^C}$","$\\beta_{G^C} - \\beta_{X x G^M}$","$\\beta_{G^M} - \\beta_{X x G^C}$","$\\beta_{G^M}$","$\\beta_{G^C}$","$\\beta_{G^C} + \\beta_{X x G^M}$","$\\beta_{G^M} + \\beta_{X x G^C}$","$\\beta_{G^C} +2 \\beta_{X x G^M}$","$\\beta_{G^M} +2 \\beta_{X x G^C}$"),res[-c(1,7,13,19),])
write.table(res.sub,"res_all.tex",row.names=F,quote=F,sep="&")

res.haplin = round(cbind(bias.indep01.haplin.vec,emp_se.indep01.haplin.vec,mean_se.indep01.haplin.vec,coverage.indep01.haplin.vec,bias.dep01.haplin.vec,emp_se.dep01.haplin.vec,mean_se.dep01.haplin.vec,coverage.dep01.haplin.vec),3)
res.haplin = cbind(rep(-2:2,rep(2,5)),rep(c("$G^M$","$G^C$"),5),res.haplin) 
write.table(res.haplin,"res_haplin.tex",row.names=F,quote=F,sep="&")


