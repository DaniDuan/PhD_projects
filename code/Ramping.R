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

### for timing each relative abundance with OD
Ramping = Ramping[4:11]
Ramping = Ramping[-seq(1, nrow(Ramping), 8),]
Ramping = Ramping[-seq(7, nrow(Ramping), 7),]
rownames(Ramping) = seq(1:1050)
Ramping[which(is.na(Ramping[,2])),2] = Ramping[(which(is.na(Ramping[,2]))+1),2]
Ramping[which(is.na(Ramping[,2])),2] = Ramping[(which(is.na(Ramping[,2]))+1),2]
Ramping[,1:8] = Ramping[,1:8] - Ramping[,2]
Ramping = Ramping[,-2]
temp = c("10C","10_30C","20C","30C","Ramping")
time = c(0,1,2,3,4,4,5,6,6,7,8,8,9,10,10,11,12,12,13,14,14,15,16,16,17,18,18,19,20,20,21,22,22,23,24)
Ramping$temp = rep(temp[sort(rep(1:5,6))], 35)
Ramping$Day = sort(rep(time, 6*5))
Ramping$Rep = rep(1:6,(35*5))
######################## Mean and standard deviation ####################
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
data = as.data.frame(t(data))
row.names(data) = seq(1,(35*5))
data = data[,-2]
sd = as.data.frame(t(sd))
row.names(sd) = seq(1,(35*5))
sd = sd[,-2]

############################ Plotting OD ##############################
time = c(0,1,2,3,4,4,5,6,6,7,8,8,9,10,10,11,12,12,13,14,14,15,16,16,17,18,18,19,20,20,21,22,22,23,24)
sp = c("S18", "W02", "W03", "S18+W02", "S18+W03", "W02+W03", "All")
temp = c("10C","10_30C","20C","30C","Ramping")
color = c("blue","yellow", "darkorange", "brown", "black")

names(data) = sp
names(sd) = sp
data$temp = rep(temp, 35)
sd$temp = rep(temp,35)
data$day = sort(rep(time, 5))
sd$day = sort(rep(time,5))

data$sp_group = sort(rep(1:35,5))
# pdf("../results/Ramping_OD.pdf")
# for(i in 1:7){
#   for(j in 1:5){
#     plot(time, data[seq(j,175,5),(8-i)], main = paste(sp[i],"_",temp[j]), xlab = "Day", ylab = "OD", ylim =c(min(data[seq(j,175,5),(8-i)]-sd[seq(j,175,5),(8-i)]),max(data[seq(j,175,5),(8-i)]+sd[seq(j,175,5),(8-i)])))
#     arrows(x0=time, y0=data[seq(j,175,5),(8-i)]-sd[seq(j,175,5),(8-i)], x1=time, y1=data[seq(j,175,5),(8-i)]+sd[seq(j,175,5),(8-i)], code=3, angle=90, lwd=0.6, length = 0.05)
#   }
# }
# graphics.off()

pdf("../results/Ramping_OD_bytemp_20J.pdf")
for(i in 1:7){ # Species
  # png(filename = paste(("../results/Ramping_OD_bytemp/OD_temp_"),sp[i],".png", sep=""), width = 480, height = 480)
  plot(1, type="n", xlab="Day", ylab = "OD",
       main = sp[i], xlim = c(time[1],time[length(time)]),
       ylim =c(min(data[1:175,(8-i)]-sd[1:175,(8-i)]),max(data[1:175,(8-i)]+sd[1:175,(8-i)])))
  for(j in 1:5){ # temperature treatment
    # if(j == 3| j == 2){
      lines(time, data[seq(j,175,5),(8-i)], col = color[j], type="b", pch = 1)
      arrows(x0=time, y0=data[seq(j,175,5),(8-i)]-sd[seq(j,175,5),(8-i)],
             x1=time, y1=data[seq(j,175,5),(8-i)]+sd[seq(j,175,5),(8-i)],
             code=3, angle=90, lwd=1, length = 0.05, col = color[j])
    # }
  }
  abline(v = 13, lty = 3)
  legend("topleft", c(temp[3],temp[2]), cex = 1, col = c(color[3],color[2]), pch = 1, lwd = 1)
  # graphics.off()
}
graphics.off()

######################## Initial growth ############################
pdf("../results/Ramping_initial.pdf")
for(i in 1:7){
  for(j in 1:5){
    if(j == 1|j == 3| j == 4){
      plot(time[1:5], data[seq(j,25,5),(8-i)], main = paste(sp[i],"_",temp[j]),
           xlab = "Day", ylab = "OD",
           ylim =c(min(data[seq(j,25,5),(8-i)]-sd[seq(j,25,5),(8-i)]),max(data[seq(j,25,5),(8-i)]+sd[seq(j,25,5),(8-i)])),
           type = "o")
      arrows(x0=time[1:5], y0=data[seq(j,25,5),(8-i)]-sd[seq(j,25,5),(8-i)], x1=time[1:5], y1=data[seq(j,25,5),(8-i)]+sd[seq(j,25,5),(8-i)], code=3, angle=90, lwd=0.6, length = 0.05)
    }
  }
}
graphics.off()

pdf("../results/Ramping_initial_bytemp.pdf")
for(i in 1:7){
  # png(filename = paste(("../results/Ramping_initial_bytemp/initial_temp_"),sp[i],".png", sep=""), width = 480, height = 480)
  plot(1, type="n", xlab="Day", ylab = "OD",
       main = sp[i], xlim = c(time[1],time[5]),
       ylim =c(min(data[1:30,(8-i)]-sd[1:30,(8-i)]),max(data[1:30,(8-i)]+sd[1:30,(8-i)])))
  for(j in 1:5){
    if(j == 1|j == 3| j == 4){
      lines(time[1:5], data[seq(j,25,5),(8-i)], col = color[j], type="b", pch = 1)
      arrows(x0=time[1:5], y0=data[seq(j,25,5),(8-i)]-sd[seq(j,25,5),(8-i)],
             x1=time[1:5], y1=data[seq(j,25,5),(8-i)]+sd[seq(j,25,5),(8-i)],
             code=3, angle=90, lwd=1, length = 0.05, col = color[j])
    }
    legend("topleft", c("10C", "20C", "30C"), cex = 1, col = color, pch = 1, lwd = 1)
  }
  # graphics.off()
}
graphics.off()

####################### Initial growth rate ###########################
log_od = data.frame()

for(i in 1:25){
  log_od= rbind(log_od, log(Ramping[seq((i-1)*8+2, i*8-1), 4:11] - mean(Ramping[seq((i-1)*8+2, i*8-1), 5], na.rm = T)))
}
# log_od[is.na(log_od)] = 0
# log_od[log_od == -Inf] = 0
log_od = log_od[-2]

# for(s in 1:7){
#   for (reps in 1:30) { # 6 reps * 5 treatments
#       plot(1, type="n", xlab="Day", ylab = "log(OD)",
#            main = paste(sp[s], "_", reps, sep = ""), xlim = c(time[1],time[5]),
#            ylim =c(-9, 0))
#       for(t in 1:5){
#         points(time[t], log_od[(30*(t - 1)+reps), s])
#     }
#   }
# }

r_data = data.frame()
for(s in 1:7){
  for (reps in 1:30) { # 6 reps * 5 treatments
    dp = c()
    for(t in 1:5){
      dp = c(dp, log_od[(30*(t - 1)+reps), s])
    }
    r_data = rbind(r_data, dp)
  }
}
names(r_data) = time[1:5]

pdf("../results/ramping_initial_all_reps.pdf")
data_123 = data.frame()
for(s in 1:7){
  for(i in 1:5){
    if(i == 1 || i == 3|| i == 4){
      r_123 = r_data[((s-1)*30+(i-1)*6+1) :((s-1)*30+6*i),]
      data_123 = rbind(data_123, r_123) # 6 replicates in 3 temperatures in 7 species groups
    for(reps in 1:6){
      plot(time[1:5], r_123[reps,], xlab = "Day", ylab = "log(OD)",
           main = paste(sp[s], "_", temp[i],"_",reps, sep = ""), type = "o")
      }
    }
  }
}
graphics.off()

# getting the highest slope as growth rate
slopes = data.frame(matrix(NA, nrow = nrow(data_123), ncol = 0))
for(i in 1:4){
  slopes = cbind(slopes, data_123[,i+1] - data_123[,i])
}
r = apply(slopes, 1, max, na.rm =T)

all_means = c()
all_sd = c()
for(i in 1:(3*7)){
  all_means = c(all_means, mean(r[(6*(i-1)+1):(6*i)], na.rm = T))
  all_sd = c(all_sd, sd(r[(6*(i-1)+1):(6*i)], na.rm = T))
}
all_sd[is.na(all_sd)] = 0

pdf("../results/Ramping_initial_gr_bytemp.pdf")
color_by_sp = c("darkgreen", "blue", "chocolate2", "darkblue", "burlywood4", "blueviolet", "black")
pch_by_sp = c(1,1,1,2,4,5,3)
plot(1, type="n", xlab="Temps", ylab = "Growth Rate (Day)",
     main = "Initial Growth Rates", xlim = c(10,30),
     ylim =c(min(all_means),max(all_means)))
for(s in 1:7){
  lines(c(10,20,30),all_means[(3*(s-1)+1) : (3*(s-1)+3)], col = color_by_sp[s],
        pch = pch_by_sp[s] ,lwd = 1, type="b")
  # arrows(x0=c(10,20,30), 
  #        y0=all_means[(3*(s-1)+1) : (3*(s-1)+3)] - all_sd[(3*(s-1)+1) : (3*(s-1)+3)],
  #        x1=c(10,20,30), 
  #        y1=all_means[(3*(s-1)+1) : (3*(s-1)+3)] + all_sd[(3*(s-1)+1) : (3*(s-1)+3)],
  #        code=3, angle=90, lwd=1, length = 0.05, col = color_by_sp[s])
}
legend("topleft", sp, cex = 1, col = color_by_sp, pch = pch_by_sp, bty="n",
       ncol = 2, lwd = 1, inset=c(-0.12,-0.17), xpd = T)
graphics.off()

######################### Plotting growth rate for all #######################
### initial growth rate for all treatments
data_in = data.frame()
for(s in 1:7){
  for(i in 1:5){
    r_all = r_data[((s-1)*30+(i-1)*6+1) :((s-1)*30+6*i),]
    data_in = rbind(data_in, r_all) # 6 replicates in 5 treatments in 7 species groups
  }
}

# getting the highest slope as growth rate
slopes_in = data.frame(matrix(NA, nrow = nrow(data_in), ncol = 0))
for(i in 1:4){
  slopes_in = cbind(slopes_in, data_in[,i+1] - data_in[,i])
}
r_in = apply(slopes_in, 1, max, na.rm =T)

in_means = c()
in_sd = c()
for(i in 1:(5*7)){
  in_means = c(in_means, mean(r_in[(6*(i-1)+1):(6*i)], na.rm = T))
  in_sd = c(in_sd, sd(r_in[(6*(i-1)+1):(6*i)], na.rm = T))
}
in_sd[is.na(in_sd)] = 0

in_gr = data.frame()
for(s in 1:7){ # 5 treatments, 7 species groups
  in_gr = rbind(in_gr,in_means[(5*(s-1)+1) : (5*(s-1)+5)])
}

### all growth rate close to carrying capacity
log_all = data.frame()
for(i in 26:175){
  log_all= rbind(log_all, log(Ramping[seq((i-1)*8+2, i*8-1), 4:11] - mean(Ramping[seq((i-1)*8+2, i*8-1), 5], na.rm = T)))
}
log_all = log_all[-2]

rall_data = data.frame()
for(s in 1:7){
  for (reps in 1:30) { # 6 reps * 5 treatments
    dp = c()
    for(t in 1:30){ # 35 - 5 = 30
      dp = c(dp, log_all[(30*(t - 1)+reps), s])
    }
    rall_data = rbind(rall_data, dp)
  }
}
names(rall_data) = time[6:35]

data_5 = data.frame()
for(s in 1:7){
  for(i in 1:5){
    r_5 = rall_data[((s-1)*30+(i-1)*6+1) :((s-1)*30+6*i),]
    data_5 = rbind(data_5, r_5) # 6 replicates in 5 treatments in 7 species groups
    for(reps in 1:6){
      plot(time[6:35], r_5[reps,], xlab = "Day", ylab = "log(OD)",
           main = paste(sp[s], "_", temp[i],"_",reps, sep = ""), type = "o")
    }
  }
}

# getting the highest slope as growth rate
slopes_all = data.frame(matrix(NA, nrow = nrow(data_5), ncol = 0))
for(i in 1:29){
  slopes_all = cbind(slopes_all, data_5[,i+1] - data_5[,i])
}
slopes_all = slopes_all[,-seq(3,29,3)]

all_mx_gr = data.frame(matrix(NA, nrow = nrow(slopes_all), ncol = 0))
for(i in 1:10){ # 10 groups of growth rates, compare every pair for maximum
  subset = cbind(slopes_all[2*(i-1)+1], slopes_all[2*i])
  all_mx_gr = cbind(all_mx_gr, apply(subset, 1, max, na.rm =F))
}

pdf("../results/Ramping_gr_bytemp.pdf")
all_gr_means = data.frame()
n = 0
for(i in 1:7){
  plot(1, type="n", xlab = "Day", ylab = "Growth rate",
       main = sp[i], xlim = c(4,24),
       ylim =c(0,1))
  for(j in 1:5){
    n = 5*(i-1)+j
    mean_gr = colMeans(all_mx_gr[(6*(n-1)+1):(6*n),], na.rm = T)
    all_gr_means = rbind(all_gr_means, mean_gr)
    mean_gr[mean_gr>1] = NaN
    lines(seq(5,23,2),mean_gr, col = color[j], lwd = 1, type="b", pch = 1)
  }
  legend("topright", temp, cex = 1, col = color, pch = 1, box.lty = 2, ncol = 3, lwd = 1)
  abline(v = 13, lty = 3)
}
graphics.off()
names(all_gr_means) = seq(5,23,2)
all_gr_means$temp = rep(temp,7)
all_gr_means$sp = sp[sort(rep(1:7, 5))]

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
  print(d[j])
  ra_mean = rbind(ra_mean, mean)
  }
}
names(ra_mean) = colnames(ra_p)[1:9]
# sd_int = c(apply(ra[1:3,1:6], 2,sd), apply(ra[1:3, 7:9],2,sd))
ra_time = seq(0,24,2)[-12][-2]

ra_mean$day = rep(ra_time, 5)
ra_mean$temp = temp[sort(rep(1:5,11))]

sp_ra = c(sp[sort(rep(4:7,2))],sp[7])
spinside = c("S18", "W02", "S18", "W03", "W02", "W03", "S18", "W02", "W03")
colors_ra = c("darkgreen", "blue", "darkgreen", "chocolate2", "blue", "chocolate2", "darkgreen", "blue", "chocolate2")
pdf("../results/Ramping_relative_abundances.pdf")
for(i in 1:length(temp)){
  for(j in 1:length(unique(sp_ra))){
    sp_ra_v = which(sp_ra == unique(sp_ra)[j])
    plot(1, type="n", xlab="Day", ylab = "%",
       main = paste(unique(sp_ra)[j],"_", temp[i],sep=""), xlim = c(ra_time[1],ra_time[11]),
       ylim =c(min(ra_mean[(11*(i-1)+1):(11*i),sp_ra_v], na.rm = T),
               max(ra_mean[(11*(i-1)+1):(11*i),sp_ra_v], na.rm = T)))
    lines(ra_time,ra_mean[(11*(i-1)+1):(11*i),sp_ra_v[1]], col = colors_ra[sp_ra_v][1], type="b", pch = 1)
    lines(ra_time,ra_mean[(11*(i-1)+1):(11*i),sp_ra_v[2]], col = colors_ra[sp_ra_v][2], type="b", pch = 1)
    if(length(sp_ra_v) == 3){
      lines(ra_time,ra_mean[(11*(i-1)+1):(11*i),sp_ra_v[3]], col = colors_ra[sp_ra_v][3], type="b", pch = 1)
    }
    legend("right", spinside[sp_ra_v], cex = 1, col = colors_ra[sp_ra_v], pch = 1, 
           box.lty = 3, lwd = 1)
  }
}
graphics.off()


##################### relative abundance with growth OD #################
# Initial
ra_data = data[,4:10]
ra_data_sd = sd[,4:9]
ra_data = data.frame(replicate(2, ra_data[,1],simplify = F), 
                     replicate(2, ra_data[,2],simplify = F),
                     replicate(2, ra_data[,3],simplify = F),
                     replicate(3, ra_data[,4],simplify = F),data[,8:10])
ra_data_sd = data.frame(replicate(2, ra_data_sd[,1],simplify = F), 
                     replicate(2, ra_data_sd[,2],simplify = F),
                     replicate(2, ra_data_sd[,3],simplify = F),
                     replicate(3, ra_data_sd[,4],simplify = F),data[,8:10])

ra_OD = data.frame()
ra_OD_sd = data.frame()
for(i in 1:length(unique(ra_mean$day))){
  ra_sub = subset(ra_mean, ra_mean$day == unique(ra_mean$day)[i])
  if(i<10){
    sub = subset(ra_data, ra_data$sp_group == (3*(i -1)+1)|ra_data$sp_group == (3*i-1)|ra_data$sp_group == (3*i))
    sub_sd = subset(ra_data_sd, ra_data_sd$sp_group == (3*(i -1)+1)|ra_data_sd$sp_group == (3*i-1)|ra_data_sd$sp_group == (3*i))
    rep_sub = do.call("rbind",replicate(3, ra_sub,simplify = F))
  }
  else{
    sub = subset(ra_data, ra_data$sp_group == (4*(i -3))|ra_data$sp_group == (4*(i -3)+1)|ra_data$sp_group == (4*(i -3)+2)|ra_data$sp_group == (4*(i -3)+3))
    sub_sd = subset(ra_data_sd, ra_data_sd$sp_group ==  (4*(i -3))|ra_data_sd$sp_group == (4*(i -3)+1)|ra_data_sd$sp_group == (4*(i -3)+2)|ra_data_sd$sp_group == (4*(i -3)+3))
    rep_sub = do.call("rbind",replicate(4, ra_sub,simplify = F))
  }
  ra_OD = rbind(ra_OD, sub[,1:9]*rep_sub[1:9])
  ra_OD_sd = rbind(ra_OD_sd, sub_sd[,1:9]*rep_sub[1:9])
}
names(ra_OD) = names(ra_mean[,1:9])
names(ra_OD_sd) = names(ra_mean[,1:9])
ra_OD = cbind(ra_OD, data[,8:10])
ra_OD_sd = cbind(ra_OD_sd, data[,8:10])

pdf("../results/Ramping_ra_OD.pdf")
for(i in 1:length(temp)){
  for(j in 1:length(unique(sp_ra))){
    sp_ra_v = which(sp_ra == unique(sp_ra)[j])
    sub = ra_OD[seq(i,175,5),sp_ra_v]
    plot(1, type="n", xlab="Day", ylab = "OD",
         main = paste(unique(sp_ra)[j],"_", temp[i],sep=""),
         xlim = c(time[1],time[length(time)]),
         ylim =c(min(sub, na.rm = T), max(sub, na.rm = T)))
    lines(time,sub[,1], col = colors_ra[sp_ra_v][1], type="b", pch = 1)
    lines(time,sub[,2], col = colors_ra[sp_ra_v][2], type="b", pch = 1)
    if(length(sp_ra_v) == 3){
      lines(time,sub[,3], col = colors_ra[sp_ra_v][3], type="b", pch = 1)
    }
    legend("topright", spinside[sp_ra_v], cex = 1, col = colors_ra[sp_ra_v], pch = 1, 
           box.lty = 3, lwd = 1)
  }
}
graphics.off()

################################### OD * RA (reps) ######################################################
Ramping_fat = data.frame(replicate(2, Ramping[,1],simplify = F), 
                     replicate(2, Ramping[,2],simplify = F),
                     replicate(2, Ramping[,3],simplify = F),
                     replicate(3, Ramping[,4],simplify = F),Ramping[,8:10])
names(Ramping_fat) = names(ra_p)
Ramping_fat$group = c(sort(rep(seq(1,6,1), (30*3))), rep(6,30), rep(7,60), sort(rep(8:9,(30*3))), sort(rep(10:11,(30*4))))

# group 1
ra_rep = colMeans(ra_p[1:6,1:9], na.rm = T)
ODra_rep = rbind(do.call("rbind", replicate(5, colMeans((ra_rep* Ramping_fat[Ramping_fat$Day == 0,1:9]), na.rm = T), simplify = FALSE)), 
                 do.call("rbind", replicate(5, colMeans((ra_rep* Ramping_fat[Ramping_fat$Day == 1,1:9]), na.rm = T), simplify = FALSE)),
                 do.call("rbind", replicate(5, colMeans((ra_rep* Ramping_fat[Ramping_fat$Day == 2,1:9]), na.rm = T), simplify = FALSE)))# the first part of the data frame
ODra_rep_sd = rbind(do.call("rbind", replicate(5, apply((ra_rep* Ramping_fat[Ramping_fat$Day == 0,1:9]), 2, sd), simplify = FALSE)),
                    do.call("rbind", replicate(5, apply((ra_rep* Ramping_fat[Ramping_fat$Day == 1,1:9]), 2, sd), simplify = FALSE)),
                    do.call("rbind", replicate(5, apply((ra_rep* Ramping_fat[Ramping_fat$Day == 2,1:9]), 2, sd), simplify = FALSE)))# the first part of the data frame

# groups 2-10
list = do.call("rbind", replicate(3,list(c(1,2),c(3,4),c(5,6))))
reps_data = data.frame()
ra_all = data.frame()
for(i in 1:9){
  reps_data = rbind(reps_data, 
                    subset(Ramping_fat, group == (i+1) & Rep == list[i,]))
  subset = subset(ra_p, Day == (2*(i+1)))
  b_subset = do.call("rbind",replicate((length(unique(Ramping_fat$Day[Ramping_fat$group == (i+1)]))+1), subset, simplify = FALSE))
  ra_all = rbind(ra_all, b_subset)
}

abs_od = reps_data[1:9]*ra_all[1:9]
abs_od = data.frame(abs_od, reps_data[10:13])

mean = data.frame()
sd = data.frame()
for(i in 1:28){
  subset = abs_od[((i-1)*10+1):(i*10),]
  for(j in temp){
    mean = rbind(mean, colMeans(subset[subset$temp == j,1:9],na.rm = T))
    sd = rbind(sd, apply(subset[subset$temp == j,1:9], 2, sd))
  }
}
colnames(mean) = colnames(ODra_rep)
ODra_rep = rbind(ODra_rep, mean)
colnames(sd) = colnames(ODra_rep_sd)
ODra_rep_sd = rbind(ODra_rep_sd, sd)

# group 11
ra_rep = do.call("rbind", replicate(4, subset(ra_p, Day == 24), simplify = FALSE))
data = cbind(Ramping_fat[Ramping_fat$group == 11,1:9]*ra_rep[1:9], Ramping_fat[Ramping_fat$group == 11,10:13])
mean = data.frame()
sd = data.frame()
for(i in 1:4){
  subset = data[((i-1)*30+1):(i*30),]
  for(j in temp){
    mean = rbind(mean, colMeans(subset[subset$temp == j,1:9],na.rm = T))
    sd = rbind(sd, apply(subset[subset$temp == j,1:9], 2, sd))
  }
}
colnames(mean) = colnames(ODra_rep)
ODra_rep = rbind(ODra_rep, mean)
colnames(sd) = colnames(ODra_rep_sd)
ODra_rep_sd = rbind(ODra_rep_sd, sd)

ODra_rep$day = sort(rep(time[0:35],5))
ODra_rep_sd$day = sort(rep(time[0:35],5))
ODra_rep$temp = rep(temp,35)
ODra_rep_sd$temp = rep(temp,35)

sp_ra = c(sp[sort(rep(4:7,2))],sp[7])
spinside = c("S18", "W02", "S18", "W03", "W02", "W03", "S18", "W02", "W03")
colors_ra = c("darkgreen", "blue", "darkgreen", "chocolate2", "blue", "chocolate2", "darkgreen", "blue", "chocolate2")

# pdf("../results/Ramping_ra_OD_errbar.pdf")
for(i in temp){
  for(j in unique(sp_ra)){
    sp_ra_v = which(sp_ra == j)
    sub = ODra_rep[ODra_rep$temp == i,sp_ra_v]
    sub_sd = ODra_rep_sd[ODra_rep_sd$temp == i,sp_ra_v]
    png(filename = paste(("../results/Ramping_ra_OD_errbar/"),j,"_", i,".png",sep=""), width = 1400, height = 1400)
    plot(1, type="n", xlab="Day", ylab = "OD",
         main = paste(j,"_", i,sep=""),
         xlim = c(time[1],time[length(time)]),
         ylim =c(min((sub-sub_sd), na.rm = T), max(sub, na.rm = T)+max(sub_sd,na.rm = T)),
         cex.lab=3, cex.axis=3, cex.main=4)
    lines(time,sub[,1], col = colors_ra[sp_ra_v][1], type="b", pch = 1, lwd=3)
    arrows(x0=time, y0=sub[,1]-sub_sd[,1], x1=time, y1=sub[,1]+sub_sd[,1],
           code=3, angle=90, lwd=3, length = 0.05, colors_ra[sp_ra_v][1])
    lines(time,sub[,2], col = colors_ra[sp_ra_v][2], type="b", pch = 1, lwd=3)
    arrows(x0=time, y0=sub[,2]-sub_sd[,2], x1=time, y1=sub[,2]+sub_sd[,2], 
           code=3, angle=90, lwd=3, length = 0.05, colors_ra[sp_ra_v][2])
    if(length(sp_ra_v) == 3){
      lines(time,sub[,3], col = colors_ra[sp_ra_v][3], type="b", pch = 1, lwd=3)
      arrows(x0=time, y0=sub[,3]-sub_sd[,3], x1=time, y1=sub[,3]+sub_sd[,3],
             code=3, angle=90, lwd=3, length = 0.05, col = colors_ra[sp_ra_v][3])
    }
    op <- par(cex = 4)
    legend("topright", spinside[sp_ra_v], cex = 1, col = colors_ra[sp_ra_v], pch = 1, 
           box.lty = 3, lwd = 3)
    graphics.off()
  }
}
# graphics.off()
