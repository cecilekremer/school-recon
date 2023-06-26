

load("WS_fullData.RData")

library(epicontacts)
epi.dat = make_epicontacts(data.full, contacts, id = "CaseID", from = "from", to = "to", directed = F)

#############################

## Create outbreaker data object

library(outbreaker2)
library(epicontacts)
library(adegenet)
library(igraph)
library(EpiEstim)


test.dates = sample.dates$date
names(test.dates) = sample.dates$ID

## generation time distribution 

mu = 4.7; sigma = 2.5
w = dunif(seq(1,10), 1, 10)

## time to sample collection
mu.i = 7; s.i = 2
f = dunif(seq(1,10), 1, 10)

out = outbreaker_data(dates = test.dates, ctd = epi.dat, w_dens = w, f_dens = f
                      # dna = seq.dat
                      )

## Outbreak analysis

# Default settings
# api <- get_cpp_api()
# ls(api)
# args(api$cpp_prior_mu)

cases <- c(1:59)

## define clusters
cluster <- c(rep(1,38),
             rep(2,9),
             rep(3,12))

## configuration outbreaker2
ob_config_sens2 <- list(n_iter = 1000000, # number of iterations
                        sample_every = 100, # thinning factor
                        prior_mu = 100000, # rate for exponential prior of mu (default is 1000)
                        move_kappa = TRUE, # look for missing cases
                        init_pi = 0.7, # prob of observation
                        prior_pi = c(1, 1), # flat prior
                        move_pi = TRUE,
                        init_kappa = 1,
                        init_tree = "star", # type of tree as starting point
                        min_date = -11, # 02/10/2021
                        init_mu = 0.00001,
                        prior_eps = c(1,1), # contact reporting coverage
                        prior_lambda = c(1,1), # eps and lambda supporting a wider range of values (i.e. not as much weight put on reported contacts)
                        init_eps = 0.5,
                        init_lambda = 0.1,
                        pb = T
)

set.seed(123)
res = outbreaker(out,
                 config = ob_config_sens2
)

res
burnin = 10000
plot(res, burn = burnin) # trace plot of posterios
sum(res$post == "-Inf" | is.na(res$post)) / length(res$post)
library(ggplot2)

## convergence
# mutation rate
p1 <- plot(res, "mu", burn = burnin)

# prob of observation
p2 <- plot(res, "pi", burn = burnin)
# plot(res, "pi", type = "density", burn = burnin)

p3 <- plot(res, "eps", burn = burnin)

p4 <- plot(res, "lambda", burn = burnin)

library(ggpubr)
conv <- ggarrange(p1,p2,p3,p4, ncol = 2, nrow = 2)

## consensus based on maximum posterior ancestry (may contain cycles)
res_mpa = summary(res, burnin = burnin, method = "decycle", min_support = 0) # method = "decycle" to remove cycles in the consensus transmission tree
hist(res_mpa$tree$support)

### Shannon entropy to quantify the uncertainty associated with an inferred ancestry
resnew = res
entropy <- c()
for(j in 1:dim(res_mpa$tree)[1]){
  
  resnew[[8+j]] <- ifelse(is.na(resnew[[8+j]]), 0, resnew[[8+j]])
  
  ## if kappa > 1 --> infector = unsampled case = 60
  resnew[[8+j]] <- ifelse(resnew[[8+j+59+59]] > 1, 60, resnew[[8+j]])
  
  freqs <- as.numeric(table(resnew[[8+j]])/sum(table(resnew[[8+j]])))
  entr <- c()
  for(i in 1:length(freqs)){
    entr <- c(entr, freqs[i]*log(freqs[i]))
  }
  
  entropy <- c(entropy, -sum(entr))
}
mean(entropy)
entropy
plot(entropy)
write(entropy, file = 'entropy_OutNoSeq.txt', append = F, ncolumns = 59)
length(entropy)


smry_res = summary(res, burnin = burnin, method = "decycle") # method = "decycle" to remove cycles in the consensus transmission tree
smry_res
quantile(res$mu, c(0.025,0.5,0.975))
quantile(res$pi, c(0.025,0.5,0.975))
quantile(res$eps, c(0.025,0.5,0.975))
quantile(res$lambda, c(0.025,0.5,0.975))

## Check clusters
ok = c()
for(i in 1:dim(smry_res$tree)[1]){
  if(!is.na(smry_res$tree$from[i])){
    if(cluster[smry_res$tree$from[i]] == cluster[smry_res$tree$to[i]]) {
      ok = c(ok, 1)
    }else{
      ok = c(ok, 0)
    }
  }else{
    ok = c(ok, NA)
  }
}
table(ok)
smry_res$tree[which(ok != 1),]


## Define the type of contact (school (1) vs household (2))
setting = matrix(nrow=59, ncol=59)
for(i in 1:59){
  for(j in 1:59){
    setting[i,j] = ifelse(data.full$FOYER[data.full$CaseID==i]==data.full$FOYER[data.full$CaseID==j],2,1) # HH contact if cases belong to the same Foyer
  }
}
# network <- cons_tree$contacts$from
network <- smry_res$tree$from
network[is.na(network)] <- 0
kappa <- smry_res$tree$generations
ncases = 59
for(i in 1:59){
  if(!is.na(kappa[i]) & kappa[i]>1){
    tmp = network[i] # infector of intermediate case
    network[i] = ncases+1
    network[ncases+1] = tmp
    ncases = ncases + 1
  }
}
Infector = network
Infectee = 1:length(network)
dat=data.frame(Infectee,Infector)
type = c()
for(i in 1:ncases){
  if(Infector[i]==0) type[i]=NA
  else if(Infector[i]>59 | i>59) type[i] = 3 # unobserved case
  else if(Infector[i]!=0){
    type[i] = setting[i,Infector[i]]
  }
}
tree=cbind(Infector,Infectee,type); tree=tree[which(Infector!=0),]
g = graph_from_edgelist(tree[,c(1,2)])
edge.mat = tree
# assign adult/child
adult = c()
case = c(1:max(Infector))
for(i in 1:length(case)){
  if(i <= 59) {adult[i] = ifelse(data.full$child[data.full$CaseID==i]==1, 2, 1)} # 2 = child, 1 = adult
  else{adult[i] = 3}
}
study.id = c()
for(i in 1:length(case)){
  if(i<60) {study.id[i] = as.character(data.full$ID[data.full$CaseID==i])}
  else{study.id[i] <- "X"}
}

ncases = length(case)
table(type)
## Transmission between adult-child or child-adult?
age.pair=c()
for(i in 1:length(Infector)){
  if(Infector[i]>0 & Infector[i]<=59 & Infectee[i]<=59){
    if(data.full$child[Infector[i]]==1 & data.full$child[Infectee[i]]==1){
      age.pair[i] = "child-child"
    }
    if(data.full$child[Infector[i]]==1 & data.full$child[Infectee[i]]==0){
      age.pair[i] = "child-adult"
    }
    if(data.full$child[Infector[i]]==0 & data.full$child[Infectee[i]]==1){
      age.pair[i] = "adult-child"
    }
    if(data.full$child[Infector[i]]==0 & data.full$child[Infectee[i]]==0){
      age.pair[i] = "adult-adult"
    }
  }
}
table(age.pair)

vertex.mat = data.frame(case, adult, study.id)

net = graph.data.frame(edge.mat, vertex.mat, directed=T)
library(extrafont)
colr=c("darkblue","red","grey")
colr2 = c("lightblue","pink","grey")
v.shape = c("circle","square","none")
# V(net)$color <- colr2[V(net)$adult]
V(net)$color <- NA
# E(net)$color <- colr[E(net)$type]
edge.col = c("black","black","grey")
E(net)$color <- edge.col[E(net)$type]
V(net)$shape <- v.shape[V(net)$adult]
linetype = c(1,2,3)
E(net)$lty <- linetype[E(net)$type]
par(mfrow=c(1,1))
# jpeg("FigS3.eps", width = 18, height = 15, units = 'cm', res = 300)
# pdf("network.pdf", width=18*0.3937, height=15*0.3937)
jpeg("~/Sequencing/FINAL/Outbreaker/network_noGen.jpeg", width = 20, height = 20, units = 'cm', res=300)
plot(net, vertex.size=10, edge.color=E(net)$color, vertex.shape=V(net)$shape,# vertex.label.col = "black",#vertex.label.col=V(net)$color, #vertex.label=V(net)$study.id,#  vertex.shape="none",  #layout=layout.fruchterman.reingold,
     # vertex.label.family="Times New Roman Uni",  #vertex.label.font=V(net2)$font, vertex.label.color=V(net2)$color,
     vertex.label.cex=0.6, edge.arrow.size=0.2)
legend("topleft",c("School","Household","Unobserved","Adult","Child","Unobserved"),col=c(1,1,1,1,1,"grey"), lty=c(1,2,3,NA,NA,NA), pch=c(rep(NA,3),1,0,NA),title="Transmission setting",bty="n",cex=0.7)
dev.off()


# network posterior probabilities
infector.mat = matrix(0, 59+1, 59)
for(i in 1:59){
  inf.post = (as.data.frame(res))[,i+8]
  inf.post[is.na(inf.post)] = 0
  tab = data.frame(table(inf.post)/dim((as.data.frame(res)))[1])
  for(k in 1:length(tab$inf.post)){
    infector.mat[as.numeric(as.character(tab$inf.post[k]))+1, i] = tab$Freq[k]
  }
}
library(plot.matrix)
library(RColorBrewer)
rownames(infector.mat) = c(0:59)
colnames(infector.mat) = c(1:59)

# cols <- brewer.pal(10, "RdGy")
cols <- brewer.pal(5, "Greys")
newcol <- colorRampPalette(cols)
ncols <- 50
cols2 <- newcol(ncols)

# library(psych)
# jpeg("Fig3.eps", width = 20, height = 12, units = 'cm', res = 1200)
# postscript("Fig3.eps", width=20*0.3937, height=12*0.3937)
# pdf("Fig3.pdf", width=20*0.3937, height=12*0.3937)
jpeg("~/Sequencing/FINAL/Outbreaker/kappa_v2.jpeg", width = 60, height = 30, units = 'cm', res = 300)
par(mfrow=c(1,2))
par(mar=c(5.1, 4.1, 4.1, 4.1))
# plot(infector.mat, xlab="Case", ylab="Infector", main="", fmt.key="%.3f", fmt.cell='%.3f', digits=3, text.cell=list(cex=0.5),
#      col=colorRampPalette(brewer.pal(5, "Oranges")))
plot(infector.mat, xlab="Case", ylab="Infector", main="a", fmt.key="%.2f", cex.axis=1, cex.lab = 1.5,
     col=cols2, breaks = 50, border = NA)
# kappa posterior probabilities
kappa.mat = matrix(0, 59, 5)
for(i in 1:59){
  kappa.post = (as.data.frame(res))[,i+126]
  tab = data.frame(table(kappa.post)/dim((as.data.frame(res)))[1])
  for(k in 1:length(tab$kappa.post)){
    kappa.mat[i, as.numeric(as.character(tab$kappa.post[k]))] = tab$Freq[k]
  }
}
rownames(kappa.mat) = c(1:59)
colnames(kappa.mat) = c(0,1,2,3,4)
# library(psych)
par(mar=c(5.1, 5.1, 4.1, 5.1))
# plot(infector.mat, xlab="Case", ylab="Infector", main="", fmt.key="%.3f", fmt.cell='%.3f', digits=3, text.cell=list(cex=0.5),
#      col=colorRampPalette(brewer.pal(5, "Oranges")))
plot(kappa.mat, xlab="# intermediate cases", ylab="Case", main="b", fmt.key="%.2f", cex.axis=1, cex.lab = 1.5,
     col=cols2, breaks = 50, border = NA)
dev.off()

###

age.pair.mat <- matrix(NA, nrow = 59, ncol = 4)
for(i in 1:59){
  # posterior distribution of infector for each case (direct transmission only)
  inf.post = (as.data.frame(res))[,i+8]
  inf.post[is.na(inf.post)] = 0
  
  kappa.post = (as.data.frame(res))[,i+126]
  kappa.post[is.na(kappa.post)] = 0 # k = NA if infector is NA (0)
  
  inf.post = inf.post[kappa.post == 1]
  
  # tab = data.frame(table(inf.post)/length(inf.post))
  
  # infectors = as.numeric(as.character(tab$inf.post))
  infectors = inf.post
  
  age.pair.post = c()
  k <- 1
  for(j in infectors){
    # posterior distribution of age-pair for each case
    if(data.full$child[j] == 1 & data.full$child[i]==1){
      age.pair.post[k] = "child-child"
    }else if(data.full$child[j] == 1 & data.full$child[i]==0){
      age.pair.post[k] = "child-adult"
    }else if(data.full$child[j] == 0 & data.full$child[i]==1){
      age.pair.post[k] = "adult-child"
    }else if(data.full$child[j] == 0 & data.full$child[i]==0){
      age.pair.post[k] = "adult-adult"
    }
    k <- k + 1
    
  }
  
  age.pair.mat[i, 1] = sum(age.pair.post == 'child-child')
  age.pair.mat[i, 2] = sum(age.pair.post == 'child-adult')
  age.pair.mat[i, 3] = sum(age.pair.post == 'adult-child')
  age.pair.mat[i, 4] = sum(age.pair.post == 'adult-adult')
}

colnames(age.pair.mat) = c('child-child','child-adult','adult-child','adult-adult')
age.pair.freqs = as.data.frame(t(apply(age.pair.mat, 1, function(x) x/sum(x))))
# sum(age.pair.freqs$`child-child`!=0)
# sum(age.pair.freqs$`child-adult`!=0)
rownames(age.pair.freqs) = c(1:59)

cols <- brewer.pal(5, "Greys")
newcol <- colorRampPalette(cols)
ncols <- 50
cols2 <- newcol(ncols)

jpeg("~/Sequencing/FINAL/Outbreaker/pairs_freq.jpeg", width = 20, height = 25, units = 'cm', res = 300)
par(mar=c(5.1, 4.1, 4.1, 5.1))
# plot(infector.mat, xlab="Case", ylab="Infector", main="", fmt.key="%.3f", fmt.cell='%.3f', digits=3, text.cell=list(cex=0.5),
#      col=colorRampPalette(brewer.pal(5, "Oranges")))
plot(as.matrix(age.pair.freqs), xlab="Transmission from-to", ylab="Case", main="", fmt.key="%.2f", cex.axis=1,
     col=cols2, breaks = 50, border = NA)
dev.off()


##################################
head(smry_res$tree) # consensus tree from the posterior samples of trees (ie most frequent ancestor for each case)
hist(smry_res$tree$support, col = "grey", border = "white",
     main = "With genetic information", xlim = c(0,1), xlab = 'Support for consensus ancestry')
hist(smry_res$tree$generations, col = "grey", border = "white",
     main = "Consensus ancestry: number of generations", xlim = c(0,5))
smry_res$pi

## convert infection time estimates from integers to dates
epi.dat$linelist$t_inf_est <- as.Date(smry_res$tree$time,
                                      # origin = min(sample.dates$date.org, na.rm=T))
                                      origin = minDate)


#######################
### PLOT NETWORK

# create color palettes which will be used to display information on the final graph
support_pal <- colorRampPalette(
  c("#918D98", "#645877", "#423359", "#281449", "#1A0340")
)

outcome_pal <- colorRampPalette(
  c("#3288BD", "#ABDDA4")
)

?vis_epicontacts

smry_res$tree$from = ifelse(is.na(smry_res$tree$from), 0, smry_res$tree$from)
smry_res$tree$generations = ifelse(is.na(smry_res$tree$generations), 1, smry_res$tree$generations)

smry_res$tree <- smry_res$tree[smry_res$tree$from != 0, ]
smry_res$tree$school = NA
for(i in 1:dim(smry_res$tree)[1]){
  if(setting[smry_res$tree$from[i], smry_res$tree$to[i]] == 2){
    smry_res$tree$school[i] = 'Household' # HH
  }else{
    smry_res$tree$school[i] = 'School' # school
  }
}
table(smry_res$tree$school)

smry_res$tree$observed <- NA
for(i in 1:dim(smry_res$tree)[1]){
  if(smry_res$tree$generations[i] > 1){
    smry_res$tree$observed[i] <- 'Indirect transmission'
  }else{
    smry_res$tree$observed[i] <- 'Direct transmission'
  }
}
table(smry_res$tree$observed)

# epi.dat$linelist$child = ifelse(epi.dat$linelist$child == 0, 'Adult', 'Child')


cons_tree <- make_epicontacts(epi.dat$linelist,
                              smry_res$tree,
                              directed = TRUE)

###### PLOT NETWORK #####

## network over time
jpeg('NetTime.jpg')
plot(cons_tree,
     # x_axis = 't_inf_est',
     node_color = 'child',
     edge_color = 'school',
     edge_linetype = 'observed',
     width = '100%',
     label = F,
     arrow_size = 0.5,
     node_size = 5,
     network_shape = 'rectangle')
dev.off()

library(ggplot2)

sub_oct <- cons_tree %>%
  subset(node_attribute = list(t_inf_est = c(as.Date(c('2020-10-01','2020-11-12'))))) %>%
  thin('contacts')

sub_dec <- cons_tree %>%
  subset(node_attribute = list(t_inf_est = c(as.Date(c('2020-12-01','2020-12-31'))))) %>%
  thin('contacts')

sub_march <- cons_tree %>%
  subset(node_attribute = list(t_inf_est = c(as.Date(c('2021-03-01','2021-03-31'))))) %>%
  thin('contacts')

### Unobserved: color = black ??
pNetO <- vis_temporal_static(sub_oct,
                             x_axis = 't_inf_est',
                             node_color = 'child',
                             # node_label = 'id',
                             # arrow_size = 0,
                             edge_linetype = 'observed',
                             edge_color = 'school',
                             root_order = 't_inf_est',
                             position_dodge = T,
                             node_size = 50,
                             edge_width = 7,
                             # node_size = 5,
                             # type = 'ttree',
                             rank_contact = 'from',
                             unlinked_pos = 'top',
                             thin = T
) + xlab('Estimated time of infection') + ylab('') +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 200)) +
  theme(legend.position = 'top') +
  theme(legend.key.size = unit(5, 'in')) +
  theme(legend.key.height = unit(2, 'in')) +
  theme(axis.text = element_text(size = 200)) +
  theme(axis.title = element_text(size = 200)) +
  scale_x_continuous(breaks = seq(as.Date('2020-10-01'),as.Date('2020-11-01'), 7))

pNetD <- vis_temporal_static(sub_dec,
                             x_axis = 't_inf_est',
                             node_color = 'child',
                             # node_label = 'id',
                             # arrow_size = 0,
                             edge_linetype = 'observed',
                             edge_color = 'school',
                             root_order = 't_inf_est',
                             position_dodge = T,
                             node_size = 50,
                             edge_width = 7,
                             # node_size = 5,
                             # type = 'ttree',
                             rank_contact = 'from',
                             unlinked_pos = 'top',
                             thin = T
) + xlab('Estimated time of infection') + ylab('') +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 200)) +
  theme(legend.position = 'none') +
  theme(axis.text = element_text(size = 200)) +
  theme(axis.title = element_text(size = 200)) +
  scale_x_continuous(breaks = seq(as.Date('2020-10-01'),as.Date('2020-12-31'), 7))

pNetM <- vis_temporal_static(sub_march,
                             x_axis = 't_inf_est',
                             node_color = 'child',
                             # node_label = 'id',
                             # arrow_size = 0,
                             edge_linetype = 'observed',
                             edge_color = 'school',
                             root_order = 't_inf_est',
                             position_dodge = T,
                             node_size = 50,
                             edge_width = 7,
                             # node_size = 5,
                             # type = 'ttree',
                             rank_contact = 'from',
                             unlinked_pos = 'top',
                             thin = T
) + xlab('Estimated time of infection') + ylab('') +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 200)) +
  theme(legend.position = 'none') +
  theme(axis.text = element_text(size = 200)) +
  theme(axis.title = element_text(size = 200)) +
  scale_x_continuous(breaks = seq(as.Date('2021-03-01'),as.Date('2021-03-19'), 3))

# pNet

# jpeg('network_march_270722.jpg',  width = 500, height = 300, units = 'cm', res = 100)
# pNet
# dev.off()

library(ggpubr)
library(gridExtra)
jpeg('network_all_150922.jpeg', width = 400, height = 500, units = 'cm', res = 100)
grid.arrange(pNetO, pNetD, pNetM, 
             widths = c(1,1),
             layout_matrix = rbind(c(1,1), c(2,3)))
dev.off()

jpeg('network_oct_150922.jpeg', width = 400, height = 350, units = 'cm', res = 100)
pNetO
dev.off()

jpeg('network_dec_mar_150922.jpeg', width = 400, height = 200, units = 'cm', res = 100)
grid.arrange(pNetD, pNetM, ncol = 2)
dev.off()


pNetAll <- vis_temporal_static(cons_tree,
                               type = 'ttree',
                               x_axis = 't_inf_est',
                               # network_shape = 'rectangle',
                               node_color = 'child',
                               edge_linetype = 'observed',
                               # node_label = 'id',
                               edge_color = 'school',
                               # edge_alpha = 'support',
                               root_order = 't_inf_est',
                               # node_order = 't_inf_est',
                               position_dodge = T,
                               rank_contact = 'from',
                               unlinked_pos = 'top',
                               # igraph_type = 'fr',
                               # lineend = 'square',
                               # annot = F,
                               arrow_size = 0,
                               # editor = T,
                               edge_width = 0.5,
                               node_size = 5,
                               # font_size = 15,
                               col_pal = c('Adult' = 'grey', 'Child' = 'white')#,'Household' = 'grey', 'School' = 'black')
                               # x_axis = 't_inf_est',
                               # network_shape = 'rectangle',
                               # node_color = 'child',
                               # col_pal = c('Adult' = 'grey', 'Child' = 'white'),
                               # edge_linetype = 'observed',
                               # arrow_size = 0.1,
                               # node_size = 10
) + xlab('Estimated time of infection') + ylab('') +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 20)) +
  scale_x_continuous(breaks = seq(as.Date('2020-10-01'),as.Date('2021-03-19'), 14))
pNetAll
# ggsave('~/Sequencing/FINAL/Outbreaker/network_all.jpg', pNetAll, dpi = 600, width = 120, height = 70, units = 'cm')

jpeg('network_all_1506.jpg',  width = 150, height = 100, units = 'cm', res = 300)
pNetAll
dev.off()


pNet1 <- vis_temporal_static(sub_fall,
                             x_axis = 't_inf_est',
                             network_shape = 'rectangle',
                             node_color = 'child',
                             col_pal = c('Adult' = 'grey', 'Child' = 'white'),
                             edge_linetype = 'observed',
                             arrow_size = 0.1,
                             node_size = 10
) + xlab('Estimated time of infection') +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 30)) +
  theme(axis.text = element_text(size = 30)) +
  theme(axis.title = element_text(size = 30)) +
  scale_x_continuous(breaks = seq(as.Date('2020-10-01'),as.Date('2020-12-29'), 14))
jpeg('network_fall_0606.jpg', width = 120, height = 70, units = 'cm', res = 300)
pNet1
dev.off()
# ggsave('network_fall_0606.jpg', pNet1, dpi = 300, width = 120, height = 70, units = 'cm')


pNet2 <- vis_temporal_static(sub_march,
                             x_axis = 't_inf_est',
                             network_shape = 'rectangle',
                             node_color = 'child',
                             col_pal = c('Adult' = 'grey', 'Child' = 'white'),
                             edge_linetype = 'observed',
                             arrow_size = 0.1,
                             node_size = 10
) + xlab('Estimated time of infection') +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 30)) +
  theme(axis.text = element_text(size = 30)) +
  theme(axis.title = element_text(size = 30)) +
  scale_x_continuous(breaks = seq(as.Date('2021-03-01'),as.Date('2021-03-17'), 3))
ggsave('network_march_sens2.jpg', pNet2, dpi = 300, width = 60, height = 50, units = 'cm')


