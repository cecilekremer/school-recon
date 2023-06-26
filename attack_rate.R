

### Bias in attack rate calculated as (obs - true)/true

ndatasets = 100

for(c in c(25,50,75,100)){

  for(is in 1:ndatasets){

    datID <- paste0('/Sens4/',c,'/InfNetw_Sim',is,'.csv')

    if(file.exists(datID)){

      dat <- read.csv(datID, header = T)

      dat <- dat[, -1]

      ## True attack rate

      dat$Case <- ifelse(!is.na(dat$InfectionTime), 1, 0)

      dat$InfDay <- floor(dat$InfectionTime) + 1

      dat$PosDay <- floor(dat$DayPosTest) + 1

      ## Weekly true attack rate

      dates.true <- plyr::count(dat$InfDay[!is.na(dat$InfDay)])

      inc.true <- c()
      for(i in 1:100){
        if(i %in% dates.true$x){
          inc.true[i] <- dates.true$freq[dates.true$x == i]
        }else{
          inc.true[i] <- 0
        }
      }

      AR.true <- c()
      for(j in seq(2,100,7)){
        if(j == 100){
          AR.true <- AR.true
        }else if(j == 93){
          # }else if(j == 89){
          # AR.true <- c(AR.true, sum(inc.true[j:(j+8)])/dim(dat)[1])
          AR.true <- c(AR.true, sum(inc.true[j:(j+7)])/dim(dat)[1])
        }else{
          # AR.true <- c(AR.true, sum(inc.true[j:(j+6)])/dim(dat)[1])
          AR.true <- c(AR.true, sum(inc.true[j:(j+6)])/dim(dat)[1])
        }
      }

      inc.t <- c()
      for(j in seq(2,100,7)){
        if(j == 100){
          inc.t <- inc.t
        }else if(j == 93){
          inc.t <- c(inc.t, sum(inc.true[j:(j+7)]))
        }else{
          inc.t <- c(inc.t, sum(inc.true[j:(j+6)]))
        }
      }

      ## Observed attack rate

      data <- dat[dat$Compliance==1, ]

      dates.obs <- plyr::count(data$PosDay[data$PosDay != 'Inf'])

      inc.obs <- c()
      for(i in 1:100){
        if(i %in% dates.obs$x){
          inc.obs[i] <- dates.obs$freq[dates.obs$x == i]
        }else{
          inc.obs[i] <- 0
        }
      }

      AR.obs <- c()
      for(j in seq(2,100,7)){
        if(j == 100){
          AR.obs <- AR.obs
        }else if(j == 93){
          AR.obs <- c(AR.obs, sum(inc.obs[j:(j+7)])/dim(data)[1])
        }else{
          AR.obs <- c(AR.obs, sum(inc.obs[j:(j+6)])/dim(data)[1])
        }
      }

      inc.o <- c()
      for(j in seq(2,100,7)){
        if(j == 100){
          inc.o <- inc.o
        }else if(j == 93){
          inc.o <- c(inc.o, sum(inc.obs[j:(j+7)])*(1/c*100))
        }else{
          inc.o <- c(inc.o, sum(inc.obs[j:(j+6)])*(1/c*100))
        }
      }

      ## Bias in attack rate

      AR.bias <- (AR.obs - AR.true) / AR.true
      mean.bias <- mean(AR.bias)

      ## Write output

      for(k in 1:14){
        write(c(c/100, is, k, AR.bias[k], AR.obs[k], AR.true[k], inc.o[k], inc.t[k]),
              file = paste0("/AR_Sens4_Compl",c,"_week",k,"_net",is,".txt"), ncolumns = 8, append = F)
      }


    }


  }

  for(i in 1:ndatasets){

    for(k in 1:14){

      input.name = paste0("/AR_weekly_230523/AR_Sens4_Compl",c,"_week",k,"_net",i,".txt")
      if(file.exists(input.name)){
        input1 = read.table(file=input.name)
        write.table(input1, file=paste0("AR_weekly_Sens4_Compl",c,"_net_full.csv"), append=T, quote=F, row.names=F, col.names=F, sep=",")

      }
    }

  }

}

#########################################################################################################
#########################################################################################################

### Plot bias

AR_Base_compl0.25 <- read.csv('./AR_weekly_230523/AR_weekly_Base_Compl25_net_full.csv', header = F)
AR_Base_compl0.5 <- read.csv('./AR_weekly_230523/AR_weekly_Base_Compl50_net_full.csv', header = F)
AR_Base_compl0.75 <- read.csv('./AR_weekly_230523/AR_weekly_Base_Compl75_net_full.csv', header = F)
AR_Base_compl1 <- read.csv('./AR_weekly_230523/AR_weekly_Base_Compl100_net_full.csv', header = F)

AR_Base_all <- rbind(AR_Base_compl0.25, AR_Base_compl0.5, AR_Base_compl0.75, AR_Base_compl1)
head(AR_Base_all)
AR_Base_all$prop.sampled <- AR_Base_all$V1
AR_Base_all$netID <- AR_Base_all$V2
AR_Base_all$week <- AR_Base_all$V3
AR_Base_all$AR.bias <- AR_Base_all$V4
AR_Base_all$AR.obs <- AR_Base_all$V5
AR_Base_all$AR.true <- AR_Base_all$V6
AR_Base_all$inc.obs <- AR_Base_all$V7
AR_Base_all$inc.true <- AR_Base_all$V8

AR_Sens2_compl0.25 <- read.csv('./AR_weekly_230523/AR_weekly_Sens2_Compl25_net_full.csv', header = F)
AR_Sens2_compl0.5 <- read.csv('./AR_weekly_230523/AR_weekly_Sens2_Compl50_net_full.csv', header = F)
AR_Sens2_compl0.75 <- read.csv('./AR_weekly_230523/AR_weekly_Sens2_Compl75_net_full.csv', header = F)
AR_Sens2_compl1 <- read.csv('./AR_weekly_230523/AR_weekly_Sens2_Compl100_net_full.csv', header = F)

AR_Sens2_all <- rbind(AR_Sens2_compl0.25, AR_Sens2_compl0.5, AR_Sens2_compl0.75, AR_Sens2_compl1)
head(AR_Sens2_all)
AR_Sens2_all$prop.sampled <- AR_Sens2_all$V1
AR_Sens2_all$netID <- AR_Sens2_all$V2
AR_Sens2_all$week <- AR_Sens2_all$V3
AR_Sens2_all$AR.bias <- AR_Sens2_all$V4
AR_Sens2_all$AR.obs <- AR_Sens2_all$V5
AR_Sens2_all$AR.true <- AR_Sens2_all$V6
AR_Sens2_all$inc.obs <- AR_Sens2_all$V7
AR_Sens2_all$inc.true <- AR_Sens2_all$V8
# 
AR_Sens1_compl0.25 <- read.csv('./AR_weekly_230523/AR_weekly_Sens1_Compl25_net_full.csv', header = F)
AR_Sens1_compl0.5 <- read.csv('./AR_weekly_230523/AR_weekly_Sens1_Compl50_net_full.csv', header = F)
AR_Sens1_compl0.75 <- read.csv('./AR_weekly_230523/AR_weekly_Sens1_Compl75_net_full.csv', header = F)
AR_Sens1_compl1 <- read.csv('./AR_weekly_230523/AR_weekly_Sens1_Compl100_net_full.csv', header = F)

AR_Sens1_all <- rbind(AR_Sens1_compl0.25, AR_Sens1_compl0.5, AR_Sens1_compl0.75, AR_Sens1_compl1)
head(AR_Sens1_all)
AR_Sens1_all$prop.sampled <- AR_Sens1_all$V1
AR_Sens1_all$netID <- AR_Sens1_all$V2
AR_Sens1_all$week <- AR_Sens1_all$V3
AR_Sens1_all$AR.bias <- AR_Sens1_all$V4
AR_Sens1_all$AR.obs <- AR_Sens1_all$V5
AR_Sens1_all$AR.true <- AR_Sens1_all$V6
AR_Sens1_all$inc.obs <- AR_Sens1_all$V7
AR_Sens1_all$inc.true <- AR_Sens1_all$V8

AR_Sens3_compl0.25 <- read.csv('./AR_weekly_230523/AR_weekly_Sens3_Compl25_net_full.csv', header = F)
AR_Sens3_compl0.5 <- read.csv('./AR_weekly_230523/AR_weekly_Sens3_Compl50_net_full.csv', header = F)
AR_Sens3_compl0.75 <- read.csv('./AR_weekly_230523/AR_weekly_Sens3_Compl75_net_full.csv', header = F)
AR_Sens3_compl1 <- read.csv('./AR_weekly_230523/AR_weekly_Sens3_Compl100_net_full.csv', header = F)

AR_Sens3_all <- rbind(AR_Sens3_compl0.25, AR_Sens3_compl0.5, AR_Sens3_compl0.75, AR_Sens3_compl1)
head(AR_Sens3_all)
AR_Sens3_all$prop.sampled <- AR_Sens3_all$V1
AR_Sens3_all$netID <- AR_Sens3_all$V2
AR_Sens3_all$week <- AR_Sens3_all$V3
AR_Sens3_all$AR.bias <- AR_Sens3_all$V4
AR_Sens3_all$AR.obs <- AR_Sens3_all$V5
AR_Sens3_all$AR.true <- AR_Sens3_all$V6
AR_Sens3_all$inc.obs <- AR_Sens3_all$V7
AR_Sens3_all$inc.true <- AR_Sens3_all$V8

AR_Sens4_compl0.25 <- read.csv('./AR_weekly_230523/AR_weekly_Sens4_Compl25_net_full.csv', header = F)
AR_Sens4_compl0.5 <- read.csv('./AR_weekly_230523/AR_weekly_Sens4_Compl50_net_full.csv', header = F)
AR_Sens4_compl0.75 <- read.csv('./AR_weekly_230523/AR_weekly_Sens4_Compl75_net_full.csv', header = F)
AR_Sens4_compl1 <- read.csv('./AR_weekly_230523/AR_weekly_Sens4_Compl100_net_full.csv', header = F)

AR_Sens4_all <- rbind(AR_Sens4_compl0.25, AR_Sens4_compl0.5, AR_Sens4_compl0.75, AR_Sens4_compl1)
head(AR_Sens4_all)
AR_Sens4_all$prop.sampled <- AR_Sens4_all$V1
AR_Sens4_all$netID <- AR_Sens4_all$V2
AR_Sens4_all$week <- AR_Sens4_all$V3
AR_Sens4_all$AR.bias <- AR_Sens4_all$V4
AR_Sens4_all$AR.obs <- AR_Sens4_all$V5
AR_Sens4_all$AR.true <- AR_Sens4_all$V6
AR_Sens4_all$inc.obs <- AR_Sens4_all$V7
AR_Sens4_all$inc.true <- AR_Sens4_all$V8

AR_Base_all$Scenario <- 'Baseline'
AR_Sens2_all$Scenario <- 'Symptomatic isolation'
AR_Sens1_all$Scenario <- 'Testing 2x / week'
AR_Sens3_all$Scenario <- 'Higher symp. prop.'
AR_Sens4_all$Scenario <- '50% immune adults'


AR_all <- rbind(AR_Base_all, AR_Sens1_all, AR_Sens2_all, AR_Sens3_all, AR_Sens4_all)

## Absolute difference?
AR_all2 <- AR_all
AR_all2$week <- as.factor(AR_all2$week)

for(i in 1:14){
  levels(AR_all2$week)[as.character(levels(AR_all2$week))==i] <- paste0('Week ', i)
}
head(AR_all2)

AR_all2$AR.dev <- (AR_all2$inc.obs/477*100)-(AR_all2$inc.true/477*100)
AR_all2$AR.dev <- ifelse((AR_all2$Scenario=='Symptomatic isolation' )#| AR_all2$Scenario == 'Reactive screening')
                         & AR_all2$prop.sampled!=1,
                         NA,
                         AR_all2$AR.dev)
AR_all2$prop.sampled <- ifelse((AR_all2$Scenario=='Symptomatic isolation' ),#| AR_all2$Scenario == 'Reactive screening'),
                         0,
                         AR_all2$prop.sampled)

## Select scenarios to plot
table(AR_all2$Scenario)
AR_all_plot <- AR_all2[AR_all2$Scenario == "Baseline" | 
                         AR_all2$Scenario == "Higher symp. prop." | 
                         AR_all2$Scenario == "50% immune adults",]

library(ggplot2)
jpeg('AR_weekly_sens.jpeg', width = 45, height = 25, units = 'cm', res = 300)
# ggplot(df.summary2, aes(x = as.character(prop.sampled), y = mean, ymin = min, ymax = max)) +
# geom_errorbar(width = 0.1) + geom_point(size = 1) +
ggplot(AR_all_plot, aes(x = as.character(prop.sampled), y = AR.dev, fill = Scenario)) +
  geom_violin(position = position_dodge(0.9), trim = T) +
  stat_summary(fun = median, geom = "point", shape = 4, size = 1, position = position_dodge(1))+
  facet_wrap(~ week, ncol = 7) +
  geom_hline(yintercept = 0, color = 'red', linetype = 'dotted') +
  ylab('Absolute deviation in weekly positivity rate (%)') +
  xlab('Proportion of school sampled') +
  theme(legend.position = "top") +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 20)) +
  theme(legend.text = element_text(size = 15)) +
  theme(legend.title = element_text(size = 20)) +
  theme(strip.text.x = element_text(size = 15)) # facet legend
dev.off()
