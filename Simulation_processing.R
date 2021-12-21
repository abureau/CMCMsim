# Script to compute the statistics in Table 1 of the manuscript
# Methods and software to analyze gene-environment interactions under a case-mother control-mother design with partially missing child genotype
# by Alexandre Bureau, Yuang Tian, Patrick Levallois, Yves Giguere, Jinbo Chen, Hong Zhang

# It assumes the objects CMCM.dep01.res and CMCM.indep01.res are already present in the R environment.
# These objects are created when executing the scripts Data_Generation_dep.R and Data_Generation_indep.R respectively

beta.vec = beta[-4]
# Indices of beta estimates
bind = c(1:6,14:19,28:33,42:47)
# Indices of standard error estimates
seind = c(8:13,22:27,35:40,48:53)

# Dependence case

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
    # Setting outlying estimates and their SE to NA
    CMCM.dep01.res[28:41,outliers.vec] = NA

bias.dep01.vec = apply(CMCM.dep01.res[bind,],1,mean,na.rm=T) - beta.vec
emp_se.dep01.vec = apply(CMCM.dep01.res[bind,],1,sd,na.rm=T)
mean_se.dep01.vec = apply(CMCM.dep01.res[seind,],1,mean,na.rm=T)
lower.bound = CMCM.dep01.res[bind,] - qnorm(0.975)*CMCM.dep01.res[seind,]
upper.bound = CMCM.dep01.res[bind,] + qnorm(0.975)*CMCM.dep01.res[seind,]
coverage.dep01.vec = apply(lower.bound<beta.vec & upper.bound>beta.vec,1,mean,na.rm=T)

# Independence case

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
    # Setting outlying estimates and their SE to NA
    CMCM.indep01.res[28:41,outliers.vec] = NA

bias.indep01.vec = apply(CMCM.indep01.res[bind,],1,mean,na.rm=T) - beta.vec
emp_se.indep01.vec = apply(CMCM.indep01.res[bind,],1,sd,na.rm=T)
mean_se.indep01.vec = apply(CMCM.indep01.res[seind,],1,mean,na.rm=T)
lower.bound = CMCM.indep01.res[bind,] - qnorm(0.975)*CMCM.indep01.res[seind,]
upper.bound = CMCM.indep01.res[bind,] + qnorm(0.975)*CMCM.indep01.res[seind,]
coverage.indep01.vec = apply(lower.bound<beta.vec & upper.bound>beta.vec,1,mean,na.rm=T)

# Gathering the results in a matrix and writing to a file in a format to facilitate creating Table 1
res = round(cbind(bias.indep01.vec,emp_se.indep01.vec,mean_se.indep01.vec,coverage.indep01.vec,bias.dep01.vec,emp_se.dep01.vec,mean_se.dep01.vec,coverage.dep01.vec),3)
# We remove results for the intercept
res.sub = cbind(rep(c("$G^M$","$G^C$","$X_1$","$X_1 x G^M$","$X_1 x G^C$"),4),res[-c(1,7,13,19),])
write.table(res.sub,"res.tex",row.names=F,quote=F,sep="&")
