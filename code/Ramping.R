rm(list=ls())
dev.off()

r1 = read.csv("../data/Ramping/OD_readings/Danica_Ramping_rr1.csv")
r2 = read.csv("../data/Ramping/OD_readings/Danica_Ramping_rr2.csv")
r3 = read.csv("../data/Ramping/OD_readings/Danica_Ramping_rr3.csv")
r4 = read.csv("../data/Ramping/OD_readings/Danica_Ramping_rr4.csv")

Ramping = rbind(r1,r2,r3,r4)
Ramping = Ramping[!is.na(Ramping$X.11),] # removing rows without data
Ramping = Ramping[,-1] # removing numbering col
Ramping = Ramping[,-ncol(Ramping)] # removing wavelength col
names(Ramping) = seq(1,ncol(Ramping))

############################## Mean and standard deviation ################
data = data.frame()
sd = data.frame()
for(j in 4:11){
  d = c()
  e = c()
  for(i in 1:(35*5)){ # 35 plate readings * 5 treatments
  d = c(d,mean(Ramping[seq((i-1)*8+2, i*8-1),j] - Ramping[seq((i-1)*8+2, i*8-1),5], na.rm = T))
  e = c(e,sd(Ramping[seq((i-1)*8+2, i*8-1),j] - Ramping[seq((i-1)*8+2, i*8-1),5], na.rm = T))
  }
  data = rbind(data,d)
  sd = rbind(sd, e)
}
data = t(data)
row.names(data) = seq(1,(35*5))
data = data[,-2]
sd = t(sd)
row.names(sd) = seq(1,(35*5))
sd = sd[,-2]

############################ Plotting OD ##############################
time = c(0,1,2,3,4,4,5,6,6,7,8,8,9,10,10,11,12,12,13,14,14,15,16,16,17,18,18,19,20,20,21,22,22,23,24)
sp = c("S18", "W02", "W03", "S18+W02", "S18+W03", "W02+W03", "All")
temp = c("10C","10_30C","20C","30C","Ramping")
color = c("blue","yellow", "darkorange", "brown", "black")

# pdf("../results/Ramping_OD.pdf")
# for(i in 1:7){
#   for(j in 1:5){
#     plot(time, data[seq(j,(174+j),5),(8-i)], main = paste(sp[i],"_",temp[j]), xlab = "Day", ylab = "OD", ylim =c(min(data[seq(j,(174+j),5),(8-i)]-sd[seq(j,(174+j),5),(8-i)]),max(data[seq(j,(174+j),5),(8-i)]+sd[seq(j,(174+j),5),(8-i)])))
#     arrows(x0=time, y0=data[seq(j,(174+j),5),(8-i)]-sd[seq(j,(174+j),5),(8-i)], x1=time, y1=data[seq(j,(174+j),5),(8-i)]+sd[seq(j,(174+j),5),(8-i)], code=3, angle=90, lwd=0.6, length = 0.05)
#   }
# }
# graphics.off()

pdf("../results/Ramping_OD_bytemp.pdf")
for(i in 1:7){
  plot(1, type="n", xlab="Day", ylab = "OD",
       main = sp[i], xlim = c(time[1],time[length(time)]),
       ylim =c(min(data[1:175,(8-i)]-sd[1:175,(8-i)]),max(data[1:175,(8-i)]+sd[1:175,(8-i)])))
  for(j in 1:5){
  points(time, data[seq(j,(174+j),5),(8-i)], col = color[j])
  arrows(x0=time, y0=data[seq(j,(174+j),5),(8-i)]-sd[seq(j,(174+j),5),(8-i)],
           x1=time, y1=data[seq(j,(174+j),5),(8-i)]+sd[seq(j,(174+j),5),(8-i)],
           code=3, angle=90, lwd=1, length = 0.05, col = color[j])
  }
  abline(v = 13, lty = 3)
  legend("topleft", c("10C","10_30C", "20C", "30C", "Ramping"), cex = 1, col = color, pch = 1)
}
graphics.off()

######################### Initial growth ############################
# pdf("../results/Ramping_initial.pdf")
# for(i in 1:7){
#   for(j in 1:5){
#     if(j == 1|j == 3| j == 4){
#       plot(time[1:5], data[seq(j,(24+j),5),(8-i)], main = paste(sp[i],"_",temp[j]), xlab = "Day", ylab = "OD", ylim =c(min(data[seq(j,(24+j),5),(8-i)]-sd[seq(j,(24+j),5),(8-i)]),max(data[seq(j,(24+j),5),(8-i)]+sd[seq(j,(24+j),5),(8-i)])))
#       arrows(x0=time[1:5], y0=data[seq(j,(24+j),5),(8-i)]-sd[seq(j,(24+j),5),(8-i)], x1=time[1:5], y1=data[seq(j,(24+j),5),(8-i)]+sd[seq(j,(24+j),5),(8-i)], code=3, angle=90, lwd=0.6, length = 0.05)
#     }
#   }
# }
# graphics.off()

pdf("../results/Ramping_initial_bytemp.pdf")
for(i in 1:7){
  plot(1, type="n", xlab="Day", ylab = "OD",
       main = sp[i], xlim = c(time[1],time[5]),
       ylim =c(min(data[1:30,(8-i)]-sd[1:30,(8-i)]),max(data[1:30,(8-i)]+sd[1:30,(8-i)])))
  for(j in 1:5){
    if(j == 1|j == 3| j == 4){
      points(time[1:5], data[seq(j,(24+j),5),(8-i)], col = color[j])
      arrows(x0=time[1:5], y0=data[seq(j,(24+j),5),(8-i)]-sd[seq(j,(24+j),5),(8-i)],
             x1=time[1:5], y1=data[seq(j,(24+j),5),(8-i)]+sd[seq(j,(24+j),5),(8-i)],
             code=3, angle=90, lwd=1, length = 0.05, col = color[j])
    }
    legend("topleft", c("10C", "20C", "30C"), cex = 1, col = color, pch = 1)
  }
}
graphics.off()

####################### Relative abundance ###########################
ra = read.csv("../data/Ramping/Relative_abundance.csv")

ra_p = data.frame(matrix(NA, nrow = nrow(ra), ncol = 0))
for(i in 1:3){
  ra_p = cbind(ra_p, ra[i*2-1]/(ra[i*2-1]+ra[i*2]), ra[i*2]/(ra[i*2-1]+ra[i*2]))
}
ra_p = cbind(ra_p, ra[7:9]/apply(ra[7:9],1,sum), ra[10:12])

d = unique(ra_p$Day)
t = unique(ra_p$temp)

ra_mean = data.frame()
for(i in 1:(length(t)-1)){
  sub = rbind(ra_p[1:6,],subset(ra_p, temp == t[i+1]))
  for(j in 1:length(d)){
  mean = colMeans(sub[sub$Day==d[j],1:9], na.rm = T)
  ra_mean = rbind(ra_mean, mean)
  }
}
names(ra_mean) = colnames(ra_p)[1:9]
# sd_int = c(apply(ra[1:3,1:6], 2,sd), apply(ra[1:3, 7:9],2,sd))
ra_time = seq(0,24,2)[-12][-2]

pdf("../results/Ramping_relative_abundances.pdf")
for(i in 1:5){
  # S18+W02
  plot(1, type="n", xlab="Day", ylab = "OD",
       main = paste("S18+W02_", temp[i]), xlim = c(ra_time[1],ra_time[11]),
       ylim =c(min(col_mean[(11*(i-1)+1):(11*i),1],col_mean[(11*(i-1)+1):(11*i),2], na.rm = T),
               max(col_mean[(11*(i-1)+1):(11*i),1],col_mean[(11*(i-1)+1):(11*i),2], na.rm = T)))
  points(ra_time,col_mean[(11*(i-1)+1):(11*i),1], col = "darkgreen")
  points(ra_time,col_mean[(11*(i-1)+1):(11*i),2], col = "blue")
  legend("right", c("S18", "W02"), cex = 1, col = c("darkgreen", "blue"), pch = 1, 
         box.lty = 3)
  
  # S18+W03
  plot(1, type="n", xlab="Day", ylab = "OD",
       main = paste("S18+W03_", temp[i]), xlim = c(ra_time[1],ra_time[11]),
       ylim =c(min(col_mean[(11*(i-1)+1):(11*i),3],col_mean[(11*(i-1)+1):(11*i),4], na.rm = T),
               max(col_mean[(11*(i-1)+1):(11*i),3],col_mean[(11*(i-1)+1):(11*i),4], na.rm = T)))
  points(ra_time,col_mean[(11*(i-1)+1):(11*i),3], col = "darkgreen")
  points(ra_time,col_mean[(11*(i-1)+1):(11*i),4], col = "chocolate2")
  legend("right", c("S18", "W03"), cex = 1, col = c("darkgreen", "chocolate2"), pch = 1, box.lty = 3)

  # W02+W03
  plot(1, type="n", xlab="Day", ylab = "OD",
       main = paste("W02+W03_", temp[i]), xlim = c(ra_time[1],ra_time[11]),
       ylim =c(min(col_mean[(11*(i-1)+1):(11*i),5],col_mean[(11*(i-1)+1):(11*i),6], na.rm = T),
               max(col_mean[(11*(i-1)+1):(11*i),5],col_mean[(11*(i-1)+1):(11*i),6], na.rm = T)))
  points(ra_time,col_mean[(11*(i-1)+1):(11*i),5], col = "blue")
  points(ra_time,col_mean[(11*(i-1)+1):(11*i),6], col = "chocolate2")
  legend("right", c("W02", "W03"), cex = 1, col = c("blue", "chocolate2"), pch = 1, box.lty = 3)

  # S18+W02+W03
  plot(1, type="n", xlab="Day", ylab = "OD",
       main = paste("All_", temp[i]), xlim = c(ra_time[1],ra_time[11]),
       ylim =c(min(col_mean[(11*(i-1)+1):(11*i),7],col_mean[(11*(i-1)+1):(11*i),8],
                   col_mean[(11*(i-1)+1):(11*i),9], na.rm = T),
               max(col_mean[(11*(i-1)+1):(11*i),7],col_mean[(11*(i-1)+1):(11*i),8],
                   col_mean[(11*(i-1)+1):(11*i),9], na.rm = T)))
  points(ra_time,col_mean[(11*(i-1)+1):(11*i),7], col = "darkgreen")
  points(ra_time,col_mean[(11*(i-1)+1):(11*i),8], col = "blue")
  points(ra_time,col_mean[(11*(i-1)+1):(11*i),9], col = "chocolate2")
  legend("right", c("S18","W02", "W03"), cex = 1, 
         col = c("darkgreen","blue", "chocolate2"), pch = 1, box.lty = 3)
}
graphics.off()

