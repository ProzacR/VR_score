#TODO? make Y auto and mix-in
> library(cvq2)
> datafc=data.frame(lKd=log(data[,stulp]))
> datafc$Y<-Y
> result<-cvq2(datafc,lKd ~ Y)
> result
