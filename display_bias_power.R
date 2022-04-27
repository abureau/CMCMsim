#---- Libraries ----
library(ggplot2) #For plots
library(tidyr)   #For data pivot and %>%
library(dplyr)   #For mutate and %>%

# Power to detect gxe interaction terms

power.res = cbind(rep(c("$X_1 x G^M$","$X_1 x G^C$"),5),round(100*c(power.indep01null.vec,power.indep01null.haplin.vec),1),round(100*c(power.indep01.vec,power.indep01.haplin.vec),1),round(100*c(power.dep01null.vec,power.dep01null.haplin.vec),1),round(100*c(power.dep01.vec,power.dep01.haplin.vec),1))
write.table(power.res,"res_power_raw.tex",row.names=F,quote=F,sep="&")

dataForPlot <- function(data = list(), subset, group){
  obj <- cbind(data[[1]], data[[2]], data[[3]])[subset,] %>% as.data.frame()
  colnames(obj) <- c("5% MCAR", "20% MCAR", "20% MNAR")
  out <- pivot_longer(data = obj, cols = c("5% MCAR", "20% MCAR", "20% MNAR"))
  out$group <- rep(group, times = nrow(out))
  out$ingroup <- c(rep("CCMO.na", 3), rep("CCMO.dep", 3), rep("Spmlficmcm", 3), rep("glm", 3))
  out <- out %>% mutate_at(c("group", "ingroup", "name"), as.factor)
  return(out)
}

#Bias
bias.dep.GM  <- dataForPlot(data = list(bias.dep01.vec, bias.dep02.vec, bias.dep02MNAR.vec), subset = seq(2, 20, by=6), group = "GM")
bias.dep.GC  <- dataForPlot(data = list(bias.dep01.vec, bias.dep02.vec, bias.dep02MNAR.vec), subset = seq(3, 21, by=6), group = "GC")
bias.dep.X   <- dataForPlot(data = list(bias.dep01.vec, bias.dep02.vec, bias.dep02MNAR.vec), subset = seq(4, 22, by=6), group = "X")
bias.dep.GMX <- dataForPlot(data = list(bias.dep01.vec, bias.dep02.vec, bias.dep02MNAR.vec), subset = seq(5, 23, by=6), group = "GM x X")
bias.dep.GCX <- dataForPlot(data = list(bias.dep01.vec, bias.dep02.vec, bias.dep02MNAR.vec), subset = seq(6, 24, by=6), group = "GC x X")
bias <- rbind(bias.dep.GM, bias.dep.GC, bias.dep.X, bias.dep.GMX, bias.dep.GCX)
bias$ingroup <- factor(bias$ingroup, levels = c("CCMO.na", "CCMO.dep", "Spmlficmcm", "glm"))
bias$name <- factor(bias$name, levels = c("5% MCAR", "20% MCAR", "20% MNAR"))
bias$group <- factor(bias$group, levels = c("GM", "GC", "X", "GM x X", "GC x X"))

#Empirical SE
emp_se.dep.GM  <- dataForPlot(data = list(emp_se.dep01.vec, emp_se.dep02.vec, emp_se.dep02MNAR.vec), subset = seq(2, 20, by=6), group = "GM")
emp_se.dep.GC  <- dataForPlot(data = list(emp_se.dep01.vec, emp_se.dep02.vec, emp_se.dep02MNAR.vec), subset = seq(3, 21, by=6), group = "GC")
emp_se.dep.X   <- dataForPlot(data = list(emp_se.dep01.vec, emp_se.dep02.vec, emp_se.dep02MNAR.vec), subset = seq(4, 22, by=6), group = "X")
emp_se.dep.GMX <- dataForPlot(data = list(emp_se.dep01.vec, emp_se.dep02.vec, emp_se.dep02MNAR.vec), subset = seq(5, 23, by=6), group = "GM x X")
emp_se.dep.GCX <- dataForPlot(data = list(emp_se.dep01.vec, emp_se.dep02.vec, emp_se.dep02MNAR.vec), subset = seq(6, 24, by=6), group = "GC x X")
emp_se <- rbind(emp_se.dep.GM, emp_se.dep.GC, emp_se.dep.X, emp_se.dep.GMX, emp_se.dep.GCX)
emp_se$ingroup <- factor(emp_se$ingroup, levels = c("CCMO.na", "CCMO.dep", "Spmlficmcm", "glm"))
emp_se$name <- factor(emp_se$name, levels = c("5% MCAR", "20% MCAR", "20% MNAR"))
emp_se$group <- factor(emp_se$group, levels = c("GM", "GC", "X", "GM x X", "GC x X"))

#---- Plots ----
#Biais
pdf(paste0(path, "bias_dep1.pdf"))
ggplot(data = bias) +
  geom_bar(mapping = aes(x = name, y = value, fill = ingroup), stat="identity", position = "dodge") +
  ylim(c(min(bias$value), 0.25)) +
  geom_hline(yintercept = 0) +
  labs(y = "", x = "", title = "", fill = "") +
  facet_grid(vars(group)) +
  scale_fill_grey() +
  theme_classic()
dev.off()


#Empirical SE
pdf(paste0(path, "emp_se_dep1.pdf"))
ggplot(data = emp_se) +
  geom_bar(mapping = aes(x = name, y = value, fill = ingroup), stat="identity", position = "dodge") +
  ylim(c(0, max(emp_se$value))) +
  geom_hline(yintercept = 0) +
  labs(y = "", x = "", title = "", fill = "") +
  facet_grid(vars(group)) +
  scale_fill_grey() +
  theme_classic()
dev.off()
