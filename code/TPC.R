rm(list=ls())
dev.off()

library(minpack.lm)

d1 = read.csv("../data/TPC/Danica_TPC.csv")
d2 = read.csv("../data/TPC/Danica_TPC_1.csv")
d3 = read.csv("../data/TPC/Danica_TPC_2.csv")

TPC = rbind(d1,d2,d3)
TPC = TPC[!is.na(TPC$X.13),] # removing rows without data
TPC = TPC[,-1] # removing numbering col
TPC = TPC[,-ncol(TPC)] # removing wavelength col

### for timing each relative abundance with OD
TPC = TPC[4:9]
TPC = TPC[-seq(1, nrow(TPC), 8),]
TPC = TPC[-seq(7, nrow(TPC), 7),]
rownames(TPC) = seq(1:nrow(TPC))
TPC[,seq(1,6,2)] = TPC[,seq(1,6,2)] - TPC[,seq(2,6,2)]
TPC = TPC[,-seq(2,6,2)]
temp = c(12, 15, 18, 22, 25, 28)
time = c(0, 10.5, 14.5, 18.5, 22.5, 35, 38.5, 42.5, 46.5, 58.5, 62.5, 66.5, 70.5, 83, 
         86.5, 90.5, 94.5, 107.5, 119)
sp = c("S18", "W02", "W03")
names(TPC) = sp
TPC$temp = rep(sort(rep(temp, 6)), 19)
TPC$time = sort(rep(time, 36))
TPC$Rep = rep(1:6,(19*6))

mean = data.frame()
sd = data.frame()
for(i in temp){
  subset = subset(TPC, temp == i)
  for(j in time){
    mean = rbind(mean, colMeans(subset[subset$time == j,1:3], na.rm = T))
    sd = rbind(sd, apply(subset[subset$time == j,1:3], 2, sd))
  }
}
names(mean) = sp
names(sd) = sp
mean$temp = sort(rep(temp, 19))
sd$temp = sort(rep(temp, 19))
mean$time = rep(time, 6)
sd$time = rep(time, 6)

pdf("../results/TPC/TPC_single.pdf")
for(i in temp){
  m_subset = subset(mean, temp == i)
  sd_subset = subset(sd, temp == i)
  for(j in 1:length(unique(sp))){
    # png(filename = paste(("../results/TPC/"),sp[j],"_",i,".png", sep=""), width = 480, height = 480)
    plot(1, type="n", xlab="Hour", ylab = "OD",
     main = paste(sp[j],"_", i, "C" ,sep=""),
     xlim = c(time[1],time[length(time)]),
     ylim =c(min((m_subset[,j]-sd_subset[,j]), na.rm = T),
             max(m_subset[,j], na.rm = T)+max(sd_subset[,j],na.rm = T)))
    points(time,m_subset[,j], pch = 1)
    # lines(time,m_subset[,j], type="b", pch = 1)
    arrows(x0=time, y0=m_subset[,j]-sd_subset[,j], x1=time, y1=m_subset[,j]+sd_subset[,j],
           code=3, angle=90, lwd=0.6, length = 0.05)
    # graphics.off()
  }
}
graphics.off()

pdf("../results/TPC/TPC_single_allreps.pdf")
for(i in temp){
  for(rep in 1:6){
    subset = subset(TPC, temp == i & Rep == rep)
    subset[,1:3] = log(subset[,1:3])
    subset[subset == -Inf| subset == Inf] = NaN
    for(s in 1:length(sp)){
      plot(time,subset[,s], xlab = "Hour", ylab = "log(OD)",main = paste(sp[s],"_",i,"_",rep,sep="") , pch = 1)
    }
  }
}
graphics.off()
################################# fitting the growth model###################################
logistic_model <- function(t, r_max, K, N_0){ # The classic logistic equation
  return(N_0 * K * exp(r_max * t)/(K + N_0 * (exp(r_max * t) - 1)))
}

gompertz_model = function(t, r_max, K, N_0, t_lag){ # Modified gompertz growth model (Zwietering 1990)
  return(N_0 + (K - N_0) * exp(-exp(r_max * exp(1) * (t_lag - t)/((K - N_0) * log(10)) + 1)))
}   

t = seq(time[1],time[length(time)] , 0.1)

pdf("../results/TPC/TPCfitting_single_allreps.pdf")
all_r = data.frame()
for(i in temp){
  for(rep in 1:6){
    subset = subset(TPC, temp == i & Rep == rep)
    subset[,1:3] = log(subset[,1:3])
    subset[subset == -Inf| subset == Inf] = NaN
    if(i == 28){
      subset[6:9,1] = NaN
    }
    N_0_get = apply(subset[1:3], 2, min, na.rm =T) # minimum 
    K_get = apply(subset[1:3], 2, max, na.rm =T)
    slopes = apply(subset[1:3], 2, diff)
    slope_max = apply(slopes[-1,], 2, max, na.rm =T)
    r = c()
    for(s in 1:length(sp)){
      while(any(is.na(subset[,s]))){
        subset[which(is.na(subset[,s])),s] = subset[(which(is.na(subset[,s]))+1),s]
      }
      N_0_start = N_0_get[s]
      K_start = K_get[s]
      r_max_start = slope_max[s] 
      t_lag_start = time[which.max(diff(diff(subset[,s])))] - time[1]
      if(t_lag_start>40 || t_lag_start == 0){t_lag_start = 10}
      sub = data.frame(subset[,s], subset$time)
      names(sub) = c("N", "time")
      # if(any(is.na(subset[,s]))){sub = sub[-is.na(sub),]}
      model_fit = nlsLM(N~gompertz_model(t = time, r_max, K, N_0, t_lag), sub,
                           list(t_lag=t_lag_start, r_max=r_max_start,
                                N_0 = N_0_start, K = K_start), control = list(maxiter = 500))
      gr = gompertz_model(t, summary(model_fit)$coefficients[2], 
                          summary(model_fit)$coefficients[4], 
                          summary(model_fit)$coefficients[3], 
                          summary(model_fit)$coefficients[1])
      r_value = round(summary(model_fit)$coefficients[2], 4)
      Model = "Gompertz"
      if(r_value > 2){
        # sub = sub[-1,]
        model_fit = nlsLM(N~logistic_model(t = time, r_max, K, N_0), sub,
                             list(r_max=r_max_start, N_0 = N_0_start, K = K_start),
                          control = list(maxiter = 500))
        gr = logistic_model(t, summary(model_fit)$coefficients[1],
                            summary(model_fit)$coefficients[3],
                            summary(model_fit)$coefficients[2])
        r_value = round(summary(model_fit)$coefficients[1], 4)
        Model = "Logistic"
      }
      plot(time,subset[,s], xlab = "Hour", ylab = "log(OD)",
           main = paste(sp[s],"_",i,"_",rep,sep="") , pch = 1)
      text(60, ((max(subset[,s], na.rm = T)+min(subset[,s], na.rm = T))/2),
           paste("r=",r_value, "\nAIC=",
                 round(AIC(model_fit),3), "\nModel= ",
                 Model, sep=""))
      lines(t, gr)
      r = c(r, r_value)
    }
    all_r = rbind(all_r, r)
  }
}
# S18, 12C
graphics.off()

names(all_r) = sp
all_r$temp = sort(rep(temp,6))
all_r$Rep = rep(1:6,6)

# write.csv(all_r, "../results/TPC/Growth_rates.csv", row.names = F)
