raw = read.csv('../data/Beechtea+R2.csv', header = F)
data = data.frame()
for(i in 1:10){
  data = rbind(data,raw[((i-1)*7+3):((i-1)*7+5),2:13])
}
data <- data.frame(sapply(data, function(x) as.numeric(as.character(x))))
data[16,1] = NA # high value caused by bubbles
data[1:3,12] = data[1:3,10] # the control for samples is the same with empty control here

S18 = data.frame()
W02 = data.frame()
W03 = data.frame()
blank_sample = data.frame()
blank = data.frame()
for(i in 1:10){
  S18 = rbind(S18, c(rowMeans(data[((i-1)*3+1):((i-1)*3+3),1:3], na.rm = T),apply(data[((i-1)*3+1):((i-1)*3+3),1:3], 1, sd, na.rm=TRUE) ))
  W02 = rbind(W02, c(rowMeans(data[((i-1)*3+1):((i-1)*3+3),4:6]),apply(data[((i-1)*3+1):((i-1)*3+3),4:6], 1, sd) ))
  W03 = rbind(W03,c(rowMeans(data[((i-1)*3+1):((i-1)*3+3),7:9]),apply(data[((i-1)*3+1):((i-1)*3+3),7:9], 1, sd) ))
  blank_sample = rbind(blank_sample, c(mean(data[((i-1)*3+1):((i-1)*3+3),10]),sd(data[((i-1)*3+1):((i-1)*3+3),10])))
  blank = rbind(blank, c(mean(data[((i-1)*3+1):((i-1)*3+3),12]),sd(data[((i-1)*3+1):((i-1)*3+3),12])))
}
names = c("10_mean","20_mean","30_mean","10_sd","20_sd","30_sd")
colnames(S18) = colnames(W02) = colnames(W03) = names
colnames(blank) = colnames(blank_sample) = c("mean", "sd")

Time = c(0, 18.33333, 21.83333, 25.91667, 29.91667, 41.98333, 45.81667, 50.06667,53.95000, 66.06667)
color = c("blue", "blueviolet", "darkorange", "grey", "bisque4")

################ plotting by strains ########################

pdf(file = "../results/by_strains.pdf") 

plot(1, type="n", xlab="Time (hours)", ylab="OD", xlim=c(0, 70), ylim=c(min(S18[,1:3])-max(S18[,4:6]),max(S18[,1:3])+max(S18[,4:6])), main = "S18")
points(Time, S18$`10_mean`, col = color[1], pch = 19)
arrows(Time, S18$`10_mean`-S18$`10_sd`, Time, S18$`10_mean`+S18$`10_sd`, col = color[1], length=0.05, angle=90, code=3)
points(Time, S18$`20_mean`, col = color[2], pch = 19)
arrows(Time, S18$`20_mean`-S18$`20_sd`, Time, S18$`20_mean`+S18$`20_sd`, col = color[2], length=0.05, angle=90, code=3)
points(Time, S18$`30_mean`, col = color[3], pch = 19)
arrows(Time, S18$`30_mean`-S18$`30_sd`, Time, S18$`30_mean`+S18$`30_sd`, col = color[3], length=0.05, angle=90, code=3)
points(Time, blank_sample$mean, col = color[4], pch = 19)
arrows(Time, blank_sample$mean - blank_sample$sd,Time, blank_sample$mean + blank_sample$sd, col = color[4], length=0.05, angle=90, code=3)
points(Time, blank$mean, col = color[5], pch = 19)
arrows(Time, blank$mean - blank$sd,Time, blank$mean + blank$sd, col = color[5], length=0.05, angle=90, code=3)
legend(0,max(S18[,1:3])+max(S18[,4:6]), c("10C","20C","30C","Blank","Medium"), col = color, pch = 19)

plot(1, type="n", xlab="Time (hours)", ylab="OD", xlim=c(0, 70), ylim=c(min(W02[,1:3])-max(W02[,4:6]),max(W02[,1:3])+max(W02[,4:6])), main = "W02")
points(Time, W02$`10_mean`, col = color[1], pch = 19)
arrows(Time, W02$`10_mean`-W02$`10_sd`, Time, W02$`10_mean`+W02$`10_sd`, col = color[1], length=0.05, angle=90, code=3)
points(Time, W02$`20_mean`, col = color[2], pch = 19)
arrows(Time, W02$`20_mean`-W02$`20_sd`, Time, W02$`20_mean`+W02$`20_sd`, col = color[2], length=0.05, angle=90, code=3)
points(Time, W02$`30_mean`, col = color[3], pch = 19)
arrows(Time, W02$`30_mean`-W02$`30_sd`, Time, W02$`30_mean`+W02$`30_sd`, col = color[3], length=0.05, angle=90, code=3)
points(Time, blank_sample$mean, col = color[4], pch = 19)
arrows(Time, blank_sample$mean - blank_sample$sd,Time, blank_sample$mean + blank_sample$sd, col = color[4], length=0.05, angle=90, code=3)
points(Time, blank$mean, col = color[5], pch = 19)
arrows(Time, blank$mean - blank$sd,Time, blank$mean + blank$sd, col = color[5], length=0.05, angle=90, code=3)
legend(0,max(W02[,1:3])+max(W02[,4:6]), c("10C","20C","30C","Blank","Medium"), col = color, pch = 19)

plot(1, type="n", xlab="Time (hours)", ylab="OD", xlim=c(0, 70), ylim=c(min(W03[,1:3])-max(W03[,4:6]),max(W03[,1:3])+max(W03[,4:6])), main = "W03")
points(Time, W03$`10_mean`, col = color[1], pch = 19)
arrows(Time, W03$`10_mean`-W03$`10_sd`, Time, W03$`10_mean`+W03$`10_sd`, col = color[1], length=0.05, angle=90, code=3)
points(Time, W03$`20_mean`, col = color[2], pch = 19)
arrows(Time, W03$`20_mean`-W03$`20_sd`, Time, W03$`20_mean`+W03$`20_sd`, col = color[2], length=0.05, angle=90, code=3)
points(Time, W03$`30_mean`, col = color[3], pch = 19)
arrows(Time, W03$`30_mean`-W03$`30_sd`, Time, W03$`30_mean`+W03$`30_sd`, col = color[3], length=0.05, angle=90, code=3)
points(Time, blank_sample$mean, col = color[4], pch = 19)
arrows(Time, blank_sample$mean - blank_sample$sd,Time, blank_sample$mean + blank_sample$sd, col = color[4], length=0.05, angle=90, code=3)
points(Time, blank$mean, col = color[5], pch = 19)
arrows(Time, blank$mean - blank$sd,Time, blank$mean + blank$sd, col = color[5], length=0.05, angle=90, code=3)
legend(0,max(W03[,1:3])+max(W03[,4:6]), c("10C","20C","30C","Blank","Medium"), col = color, pch = 19)

dev.off()

################ plotting by temperature ########################
temps = c(10, 20, 30)

pdf(file = "../results/by_temp.pdf") 

for(i in 1:3){
  plot(1, type="n", xlab="Time (hours)", ylab="OD", xlim=c(0, 70), ylim=c(min(S18[,i])-max(S18[,i+3]),max(c(W03[,i], W02[,i]))+max(c(W03[,i+3], W02[,i+3]))), main = paste("T =", temps[i], "C"))
  points(Time, S18[,i], col = color[1], pch = 19)
  arrows(Time, S18[,i]-S18[,i+3], Time, S18[,i]+S18[,i+3], col = color[1], length=0.05, angle=90, code=3)
  points(Time, W02[,i], col = color[2], pch = 19)
  arrows(Time, W02[,i]-W02[,i+3], Time, W02[,i]+W02[,i+3], col = color[2], length=0.05, angle=90, code=3)
  points(Time, W03[,i], col = color[3], pch = 19)
  arrows(Time, W03[,i]-W03[,i+3], Time, W03[,i]+W03[,i+3], col = color[3], length=0.05, angle=90, code=3)
  points(Time, blank_sample$mean, col = color[4], pch = 19)
  arrows(Time, blank_sample$mean - blank_sample$sd,Time, blank_sample$mean + blank_sample$sd, col = color[4], length=0.05, angle=90, code=3)
  points(Time, blank$mean, col = color[5], pch = 19)
  arrows(Time, blank$mean - blank$sd,Time, blank$mean + blank$sd, col = color[5], length=0.05, angle=90, code=3)
  legend(0,max(c(W03[,i], W02[,i]))+max(c(W03[,i+3], W02[,i+3])), c("S18","W02","W03","Blank","Medium"), col = color, pch = 19)
}

dev.off()
