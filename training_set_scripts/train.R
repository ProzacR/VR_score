#need 'clean' score matrix first only numbers col scores, row componds
#R --no-save --args result_4c.csv < train.R
#
#VR


library(leaps);
cat("-- reading arguments\n", sep = "");
cmd_args = commandArgs();
#for (arg in cmd_args) cat("  ", arg, "\n", sep="");
cmd_args[4]
data<-as.matrix(read.table(cmd_args[4], sep=' '))
stulp=ncol(data)
data[1,]
stulp
data[,stulp]
dataf=data.frame(lKd=log(data[,stulp]))
dataf$score<-as.matrix(data[,c(1:(stulp-1))])
leaps<-regsubsets(dataf$lKd~dataf$score, data=dataf)
plot(leaps, scale="r2")
