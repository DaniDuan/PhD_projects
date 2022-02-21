library(ggplot2)
library(minpack.lm)

raw = read.csv('../data/Bacteria_growth.csv', header = F)
W02 = data.frame(Time = as.numeric(raw[2,2:9]), S1 = as.numeric(raw[3,2:9]), S2 = as.numeric(raw[4, 2:9]), S3 = as.numeric(raw[5,2:9]))
Mean = rowMeans(W02[,2:4],na.rm = T)
#W02 = cbind(W02, Mean)
#plot(W02$Time, W02$Mean)

w02_M = melt(W02, id.vars = "Time")
colnames(w02_M) = c("Time", "Strain", "Cell_Count")
#ggplot(w02_M, aes(Time,Cell_Count,color = Strain))+ geom_point()

growth = function(t, r_max, K, N_0){
  return(N_0 * K * exp(r_max * t)/(K + N_0 * (exp(r_max * t) - 1)))
}

t = seq(0,60,1)

pdf(file = "../results/W02_growth.pdf") 

p = ggplot(w02_M, aes(Time,Cell_Count,color = Strain))+
  scale_color_brewer(palette = "Dark2") + geom_point()
for(i in unique(w02_M$Strain)){
# i = unique(w02_M$Strain)[1]
  sub = subset(w02_M, Strain == i)
#  plot(sub$Time, sub$Cell_Count)
  r_get = lm(log(Cell_Count)~Time, sub[sub$Time >= 19 & sub$Time <= 39,])
  r_max_start = summary(r_get)$coefficients[2]
  #r_max_start = (log(sub$Cell_Count[sub$Time == 39])-log(sub$Cell_Count[sub$Time == 19]))/20
  N_0_start = min(na.omit(sub$Cell_Count))
  K_start = sub$Cell_Count[sub$Time == 39]
  fit = nlsLM(Cell_Count ~ growth(t = Time, r_max, K, N_0), sub,
                        list(r_max=r_max_start, N_0 = N_0_start, K = K_start), control = list(maxiter = 500))

  r_max = coef(fit)["r_max"]
  N_0 = coef(fit)["N_0"]
  K = coef(fit)["K"]
  forline = growth(t, r_max, K, N_0)
  forp = data.frame(t, forline, Strain = i)
#  plot(W02$Time, W02$Mean)
#  lines(t, growth(t, r_max, K, N_0))
  p = p+geom_line(data=forp, aes(x = t, y = forline, color = Strain)) + labs(x= "Time (hr)", y = "Cell Count/ml")
}  
p

dev.off()
