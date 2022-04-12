# Power to detect gxe interaction terms

power.res = cbind(rep(c("$X_1 x G^M$","$X_1 x G^C$"),5),round(100*c(power.indep01null.vec,power.indep01null.haplin.vec),1),round(100*c(power.indep01.vec,power.indep01.haplin.vec),1),round(100*c(power.indep02.vec,power.indep02.haplin.vec),1),round(100*c(power.dep01null.vec,power.dep01null.haplin.vec),1),round(100*c(power.dep01.vec,power.dep01.haplin.vec),1),round(100*c(power.dep02.vec,power.dep02.haplin.vec),1))
write.table(power.res,"res_power_raw.tex",row.names=F,quote=F,sep="&")

# Bias plots

# bias.indep.GMX = cbind(bias.indep01null.vec,bias.indep01.vec,bias.indep02.vec)[seq(5,23,by=6),]
# bias.indep.GCX = cbind(bias.indep01null.vec,bias.indep01.vec,bias.indep02.vec)[seq(6,24,by=6),]
# colnames(bias.indep.GMX) = c("Indep, GMxX = 0, 5% m", "Indep, GMxX = log(1.2), 5% m", "Indep, GMxX = log(1.2), 20% m")
# colnames(bias.indep.GCX) = c("Indep, GCxX = 0, 5% m", "Indep, GCxX = log(1.2), 5% m", "Indep, GMxX = log(1.2), 20% m")

# pdf("bias_indep.pdf")
# par(mfrow=c(2,1),cex=0.8)
# barplot(bias.indep.GMX,beside=T,main="GM x X")
# barplot(bias.indep.GCX,beside=T,main="GC x X",legend.text=c("CCMO.na","CCMO.dep","Spmlficmcm","glm"))
# dev.off()

bias.dep.GM = cbind(bias.dep01.vec,bias.dep02.vec,bias.dep02MNAR.vec)[seq(2,20,by=6),]
bias.dep.GC = cbind(bias.dep01.vec,bias.dep02.vec,bias.dep02MNAR.vec)[seq(3,21,by=6),]
bias.dep.X = cbind(bias.dep01.vec,bias.dep02.vec,bias.dep02MNAR.vec)[seq(4,22,by=6),]
bias.dep.GMX = cbind(bias.dep01.vec,bias.dep02.vec,bias.dep02MNAR.vec)[seq(5,23,by=6),]
bias.dep.GCX = cbind(bias.dep01.vec,bias.dep02.vec,bias.dep02MNAR.vec)[seq(6,24,by=6),]
colnames(bias.dep.GM) = colnames(bias.dep.GC) = colnames(bias.dep.X) =colnames(bias.dep.GMX) = colnames(bias.dep.GCX) = c("5% MCAR", "20% MCAR", "20% MNAR")

pdf("bias_dep.pdf")
par(mfrow=c(5,1),cex=0.8)
par(mar=c(3,4,1,2)+0.1)
barplot(bias.dep.GM,beside=T,main="GM",ylim=c(-0.25,0.25))
barplot(bias.dep.GC,beside=T,main="GC",ylim=c(-0.25,0.25))
barplot(bias.dep.X,beside=T,main="X",ylim=c(-0.25,0.25))
barplot(bias.dep.GMX,beside=T,main="GM x X",ylim=c(-0.25,0.25))
barplot(bias.dep.GCX,beside=T,main="GC x X",ylim=c(-0.25,0.25),legend.text=c("CCMO.na","CCMO.dep","Spmlficmcm","glm"),bty="n")
dev.off()

# Empirical SE plots
emp_se.dep.GM = cbind(emp_se.dep01.vec,emp_se.dep02.vec,emp_se.dep02MNAR.vec)[seq(2,20,by=6),]
emp_se.dep.GC = cbind(emp_se.dep01.vec,emp_se.dep02.vec,emp_se.dep02MNAR.vec)[seq(3,21,by=6),]
emp_se.dep.X = cbind(emp_se.dep01.vec,emp_se.dep02.vec,emp_se.dep02MNAR.vec)[seq(4,22,by=6),]
emp_se.dep.GMX = cbind(emp_se.dep01.vec,emp_se.dep02.vec,emp_se.dep02MNAR.vec)[seq(5,23,by=6),]
emp_se.dep.GCX = cbind(emp_se.dep01.vec,emp_se.dep02.vec,emp_se.dep02MNAR.vec)[seq(6,24,by=6),]
colnames(emp_se.dep.GM) = colnames(emp_se.dep.GC) = colnames(emp_se.dep.X) =colnames(emp_se.dep.GMX) = colnames(emp_se.dep.GCX) = c("5% MCAR", "20% MCAR", "20% MNAR")

pdf("emp_se_dep.pdf")
par(mfrow=c(5,1),cex=0.8)
par(mar=c(3,4,1,2)+0.1)
barplot(emp_se.dep.GM,beside=T,main="GM")
barplot(emp_se.dep.GC,beside=T,main="GC")
barplot(emp_se.dep.X,beside=T,main="X")
barplot(emp_se.dep.GMX,beside=T,main="GM x X")
barplot(emp_se.dep.GCX,beside=T,main="GC x X",legend.text=c("CCMO.na","CCMO.dep","Spmlficmcm","glm"),bty="n")
dev.off()

