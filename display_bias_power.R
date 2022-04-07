# Power to detect gxe interaction terms

power.res = cbind(rep(c("$X_1 x G^M$","$X_1 x G^C$"),5),round(100*c(power.indep01null.vec,power.indep01null.haplin.vec)),round(100*c(power.indep01.vec,power.indep01.haplin.vec)),round(100*c(power.indep02.vec,power.indep02.haplin.vec)),round(100*c(power.dep01null.vec,power.dep01null.haplin.vec)),round(100*c(power.dep01.vec,power.dep01.haplin.vec)),round(100*c(power.dep02.vec,power.dep02.haplin.vec)))
write.table(power.res,"res_power_raw.tex",row.names=F,quote=F,sep="&")

# Bias plots

bias.indep.GMX = cbind(bias.indep01null.vec,bias.indep01.vec,bias.indep02.vec)[seq(5,23,by=6),]
bias.indep.GCX = cbind(bias.indep01null.vec,bias.indep01.vec,bias.indep02.vec)[seq(6,24,by=6),]
colnames(bias.indep.GMX) = c("Indep, GMxX = 0, 5% m", "Indep, GMxX = log(1.2), 5% m", "Indep, GMxX = log(1.2), 20% m")
colnames(bias.indep.GCX) = c("Indep, GCxX = 0, 5% m", "Indep, GCxX = log(1.2), 5% m", "Indep, GMxX = log(1.2), 20% m")

pdf("bias_indep.pdf")
par(mfrow=c(2,1),cex=0.8)
barplot(bias.indep.GMX,beside=T,main="GM x X")
barplot(bias.indep.GCX,beside=T,main="GC x X",legend.text=c("CCMO.na","CCMO.dep","Spmlficmcm","glm"))
dev.off()

bias.dep.GMX = cbind(bias.dep01null.vec,bias.dep01.vec,bias.dep02.vec)[seq(5,23,by=6),]
bias.dep.GCX = cbind(bias.dep01null.vec,bias.dep01.vec,bias.dep02.vec)[seq(6,24,by=6),]
colnames(bias.dep.GMX) = c("dep, GMxX = 0, 5% m", "dep, GMxX = log(1.2), 5% m", "dep, GMxX = log(1.2), 20% m")
colnames(bias.dep.GCX) = c("dep, GCxX = 0, 5% m", "dep, GCxX = log(1.2), 5% m", "dep, GMxX = log(1.2), 20% m")

pdf("bias_dep.pdf")
par(mfrow=c(2,1),cex=0.8)
barplot(bias.dep.GMX,beside=T,main="GM x X")
barplot(bias.dep.GCX,beside=T,main="GC x X",legend.text=c("CCMO.na","CCMO.dep","Spmlficmcm","glm"))
dev.off()
