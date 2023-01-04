# Script to compute the statistics in Table 1 of the manuscript
# Methods and software to analyze gene-environment interactions under a case-mother control-mother design with partially missing child genotype
# by Alexandre Bureau, Yuang Tian, Patrick Levallois, Yves Giguere, Jinbo Chen, Hong Zhang

# It assumes the objects CMCM.dep01null.res and CMCM.indep01null.res are already present in the R environment.
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

#MCAR

# Proportion of missing beta before removing outliers
missing.beta.dep01null = apply(CMCM.dep01null.res[bind,],1,function(vec) mean(is.na(vec)))

    ## Identifying outlying estimates
    par.mad <- apply(CMCM.dep01null.res[bind,], 1, mad, na.rm=T)
    par.med <- apply(CMCM.dep01null.res[bind,], 1, median, na.rm=T) 
    outliers.mat <- CMCM.dep01null.res[bind,] < par.med - 3*sqrt(par.mad) | CMCM.dep01null.res[bind,] > par.med + 3*sqrt(par.mad)
    apply(outliers.mat, 1, sum, na.rm=T)
    # Outliers are present only for prospective likelihood
    outliers.vec <- apply(outliers.mat, 2, any, na.rm=T)
    mean(outliers.vec)
    # Setting outlying estimates and their SE to NA
    CMCM.dep01null.res[28:41,outliers.vec] = NA

bias.dep01null.vec = apply(CMCM.dep01null.res[bind,],1,mean,na.rm=T) - beta.vec
emp_se.dep01null.vec = apply(CMCM.dep01null.res[bind,],1,sd,na.rm=T)
mean_se.dep01null.vec = apply(CMCM.dep01null.res[seind,],1,mean,na.rm=T)
lower.bound = CMCM.dep01null.res[bind,] - qnorm(0.975)*CMCM.dep01null.res[seind,]
upper.bound = CMCM.dep01null.res[bind,] + qnorm(0.975)*CMCM.dep01null.res[seind,]
coverage.dep01null.vec = apply(lower.bound<beta.vec & upper.bound>beta.vec,1,mean,na.rm=T)
power.dep01null.vec = apply(lower.bound[c(5:6,11:12,17:18,23:24),] > 0| upper.bound[c(5:6,11:12,17:18,23:24),] < 0,1,mean,na.rm=T)

# Haplin
bias.dep01null.haplin.vec = apply(CMCM.dep01null.res[hind,],1,mean,na.rm=T) - hbeta.vec
emp_se.dep01null.haplin.vec = apply(CMCM.dep01null.res[hind,],1,sd,na.rm=T)
mean_se.dep01null.haplin.vec = apply(CMCM.dep01null.res[hseind,],1,mean,na.rm=T)
lower.bound = CMCM.dep01null.res[hind,] - qnorm(0.975)*CMCM.dep01null.res[hseind,]
upper.bound = CMCM.dep01null.res[hind,] + qnorm(0.975)*CMCM.dep01null.res[hseind,]
coverage.dep01null.haplin.vec = apply(lower.bound<hbeta.vec & upper.bound>hbeta.vec,1,mean,na.rm=T)
power.dep01null.haplin.vec = apply(CMCM.dep01null.res[nrow(CMCM.dep01null.res)-0:1,]<0.05,1,mean,na.rm=T)

# MNAR
# Proportion of missing beta before removing outliers
missing.beta.dep01nullMNAR = apply(CMCM.dep01nullMNAR.res[bind,],1,function(vec) mean(is.na(vec)))

## Identifying outlying estimates
par.mad <- apply(CMCM.dep01nullMNAR.res[bind,], 1, mad, na.rm=T)
par.med <- apply(CMCM.dep01nullMNAR.res[bind,], 1, median, na.rm=T) 
outliers.mat <- CMCM.dep01nullMNAR.res[bind,] < par.med - 3*sqrt(par.mad) | CMCM.dep01nullMNAR.res[bind,] > par.med + 3*sqrt(par.mad)
apply(outliers.mat, 1, sum, na.rm=T)
# Outliers are present only for prospective likelihood
outliers.vec <- apply(outliers.mat, 2, any, na.rm=T)
mean(outliers.vec)
# Setting outlying estimates and their SE to NA
CMCM.dep01nullMNAR.res[28:41,outliers.vec] = NA

bias.dep01nullMNAR.vec = apply(CMCM.dep01nullMNAR.res[bind,],1,mean,na.rm=T) - beta.vec
emp_se.dep01nullMNAR.vec = apply(CMCM.dep01nullMNAR.res[bind,],1,sd,na.rm=T)
mean_se.dep01nullMNAR.vec = apply(CMCM.dep01nullMNAR.res[seind,],1,mean,na.rm=T)
lower.bound = CMCM.dep01nullMNAR.res[bind,] - qnorm(0.975)*CMCM.dep01nullMNAR.res[seind,]
upper.bound = CMCM.dep01nullMNAR.res[bind,] + qnorm(0.975)*CMCM.dep01nullMNAR.res[seind,]
coverage.dep01nullMNAR.vec = apply(lower.bound<beta.vec & upper.bound>beta.vec,1,mean,na.rm=T)
power.dep01nullMNAR.vec = apply(lower.bound[c(5:6,11:12,17:18,23:24),] > 0| upper.bound[c(5:6,11:12,17:18,23:24),] < 0,1,mean,na.rm=T)

# Haplin
bias.dep01nullMNAR.haplin.vec = apply(CMCM.dep01nullMNAR.res[hind,],1,mean,na.rm=T) - hbeta.vec
emp_se.dep01nullMNAR.haplin.vec = apply(CMCM.dep01nullMNAR.res[hind,],1,sd,na.rm=T)
mean_se.dep01nullMNAR.haplin.vec = apply(CMCM.dep01nullMNAR.res[hseind,],1,mean,na.rm=T)
lower.bound = CMCM.dep01nullMNAR.res[hind,] - qnorm(0.975)*CMCM.dep01nullMNAR.res[hseind,]
upper.bound = CMCM.dep01nullMNAR.res[hind,] + qnorm(0.975)*CMCM.dep01nullMNAR.res[hseind,]
coverage.dep01nullMNAR.haplin.vec = apply(lower.bound<hbeta.vec & upper.bound>hbeta.vec,1,mean,na.rm=T)
power.dep01nullMNAR.haplin.vec = apply(CMCM.dep01nullMNAR.res[nrow(CMCM.dep01nullMNAR.res)-0:1,]<0.05,1,mean,na.rm=T)

# Independence case

# Proportion of missing beta before removing outliers
missing.beta.indep01null = apply(CMCM.indep01null.res[bind,],1,function(vec) mean(is.na(vec)))

    ## Identifying outlying estimates
    par.mad <- apply(CMCM.indep01null.res[bind,], 1, mad, na.rm=T)
    par.med <- apply(CMCM.indep01null.res[bind,], 1, median, na.rm=T) 
    outliers.mat <- CMCM.indep01null.res[bind,] < par.med - 3*sqrt(par.mad) | CMCM.indep01null.res[bind,] > par.med + 3*sqrt(par.mad)
    apply(outliers.mat, 1, sum, na.rm=T)
    # Outliers are present only for prospective likelihood
    outliers.vec <- apply(outliers.mat, 2, any, na.rm=T)
    mean(outliers.vec)
    # Setting outlying estimates and their SE to NA
    CMCM.indep01null.res[28:41,outliers.vec] = NA

bias.indep01null.vec = apply(CMCM.indep01null.res[bind,],1,mean,na.rm=T) - beta.vec
emp_se.indep01null.vec = apply(CMCM.indep01null.res[bind,],1,sd,na.rm=T)
mean_se.indep01null.vec = apply(CMCM.indep01null.res[seind,],1,mean,na.rm=T)
lower.bound = CMCM.indep01null.res[bind,] - qnorm(0.975)*CMCM.indep01null.res[seind,]
upper.bound = CMCM.indep01null.res[bind,] + qnorm(0.975)*CMCM.indep01null.res[seind,]
coverage.indep01null.vec = apply(lower.bound<beta.vec & upper.bound>beta.vec,1,mean,na.rm=T)
power.indep01null.vec = apply(lower.bound[c(5:6,11:12,17:18,23:24),] > 0| upper.bound[c(5:6,11:12,17:18,23:24),] < 0,1,mean,na.rm=T)

# Haplin
bias.indep01null.haplin.vec = apply(CMCM.indep01null.res[hind,],1,mean,na.rm=T) - hbeta.vec
emp_se.indep01null.haplin.vec = apply(CMCM.indep01null.res[hind,],1,sd,na.rm=T)
mean_se.indep01null.haplin.vec = apply(CMCM.indep01null.res[hseind,],1,mean,na.rm=T)
lower.bound = CMCM.indep01null.res[hind,] - qnorm(0.975)*CMCM.indep01null.res[hseind,]
upper.bound = CMCM.indep01null.res[hind,] + qnorm(0.975)*CMCM.indep01null.res[hseind,]
coverage.indep01null.haplin.vec = apply(lower.bound<hbeta.vec & upper.bound>hbeta.vec,1,mean,na.rm=T)
power.indep01null.haplin.vec = apply(CMCM.indep01null.res[nrow(CMCM.indep01null.res)-0:1,]<0.05,1,mean,na.rm=T)

# Gathering the results in a matrix and writing to a file in a format to facilitate creating Table 1
res = round(cbind(c(bias.indep01null.vec,bias.indep01null.haplin.vec),c(emp_se.indep01null.vec,emp_se.indep01null.haplin.vec),c(mean_se.indep01null.vec,mean_se.indep01null.haplin.vec),c(coverage.indep01null.vec,coverage.indep01null.haplin.vec),c(bias.dep01null.vec,bias.dep01null.haplin.vec),c(emp_se.dep01null.vec,emp_se.dep01null.haplin.vec),c(mean_se.dep01null.vec,mean_se.dep01null.haplin.vec),c(coverage.dep01null.vec,coverage.dep01null.haplin.vec)),3)
# We remove results for the intercept
res.sub = cbind(c(rep(c("$\\beta_{G^M}$","$\\beta_{G^C}$","$\\beta_{X}$","$\\beta_{X x G^M}$","$\\beta_{X x G^C}$"),4),"$\\beta_{G^C} -2 \\beta_{X x G^M}$","$\\beta_{G^M} -2 \\beta_{X x G^C}$","$\\beta_{G^C} - \\beta_{X x G^M}$","$\\beta_{G^M} - \\beta_{X x G^C}$","$\\beta_{G^M}$","$\\beta_{G^C}$","$\\beta_{G^C} + \\beta_{X x G^M}$","$\\beta_{G^M} + \\beta_{X x G^C}$","$\\beta_{G^C} +2 \\beta_{X x G^M}$","$\\beta_{G^M} +2 \\beta_{X x G^C}$"),res[-c(1,7,13,19),])
write.table(res.sub,"res_all.tex",row.names=F,quote=F,sep="&")

res.haplin = round(cbind(bias.indep01null.haplin.vec,emp_se.indep01null.haplin.vec,mean_se.indep01null.haplin.vec,coverage.indep01null.haplin.vec,bias.dep01null.haplin.vec,emp_se.dep01null.haplin.vec,mean_se.dep01null.haplin.vec,coverage.dep01null.haplin.vec),3)
res.haplin = cbind(rep(-2:2,rep(2,5)),rep(c("$G^M$","$G^C$"),5),res.haplin) 
write.table(res.haplin,"res_haplin.tex",row.names=F,quote=F,sep="&")


# Power and Type I error to detect gxe interaction terms

power.res = cbind(rep(c("$X_1 x G^M$","$X_1 x G^C$"),5),round(100*c(power.indep01null.vec,power.indep01null.haplin.vec)),round(100*c(power.indep01.vec,power.indep01.haplin.vec)),round(100*c(power.dep01null.vec,power.dep01null.haplin.vec)),round(100*c(power.dep01.vec,power.dep01.haplin.vec)))
write.table(power.res,"res_power.tex",row.names=F,quote=F,sep="&")

