library(GillespieSSA); library(sensobol); library(circlize); library(sensitivity); library(tidyr); library(reshape2); library(ggplot2); library(gridExtra)

SobolModel <- function(tf,iteration, phi, Q, N_H, N_M, N_V){
  epsilon <- 5.88*phi-8.62*phi^2+2.73*phi^3
  
  Area_matrix=phi
  Area_core=1-phi
  eta=0.174
  Area_edge <- epsilon*eta # edge area
  Area_deforest=(1-Area_edge)*Area_matrix/(Area_matrix+Area_core) # deforest area
  Area_forest=(1-Area_edge)*Area_core/(Area_matrix+Area_core) # forest area
  
  N_HC=N_H*.0 # number of human individuals in forest 
  N_HE=N_H*Area_edge/(Area_edge + Area_deforest) # number of human individuals in edge
  N_HM=N_H*Area_deforest/(Area_edge + Area_deforest) # number of human individuals in deforest
  
  N_MC=N_M*Area_forest/(Area_forest + Area_edge) # number of macaque individuals in forest
  N_ME=N_M*Area_edge/(Area_forest + Area_edge) # number of macaque individuals in edge
  N_MM=N_M*.0 # number of macaque individuals in deforest
  
  N_VC=N_V*Area_forest # number of vectors in forest
  N_VE=N_V*Area_edge # number of vectors in edge
  N_VM=N_V*Area_deforest # number of vectors in deforest
  
  f=1/3 ## frequency of biting; 1/gonotrophic cycle (day)
  q_E <- (N_HE*Q)/(N_HE*Q+N_ME*(1-Q)) # proportion of bites taken on humans in edge
  
  alpha_MC <- f # vectors on macaques in forest
  alpha_HE <- f*q_E # vectors on humans in  edge
  alpha_ME <- f*(1-q_E) # vectors on macaques in  edge
  alpha_HM <- f # vectors on humans in  deforest
  
  parms=c(C_VH=0.01, # probability of transmission from vector to human per infectious bite
          C_HV=0.45, # probability of transmission from human to vector per infectious bite
          C_VM=0.53, # probability of transmission from vector to macaque per infectious bite
          C_MV=0.27, # probability of transmission from macaque to vector per infectious bite
          gamma_H=1/14, # human recovery rate
          gamma_M=1/2132, # macaque recovery rate
          mu_H=1/29200, # human death rate
          mu_M=1/3650, # macaque death rate
          mu_V=0.15, #death rate of mosquitoes in forest
          zeta = 0.1, # 1/ duration of sporogony (Extrinsic incubation period)
          alpha_MC = alpha_MC, # include the locally defined alpha_MC
          alpha_HE = alpha_HE, # include the locally defined alpha_HE
          alpha_ME = alpha_ME, # include the locally defined alpha_ME
          alpha_HM = alpha_HM,  # include the locally defined alpha_HM
          N_HE=N_HE,
          N_HM=N_HM,
          N_MC=N_MC,
          N_ME=N_ME,
          N_VC=N_VC,
          N_VE=N_VE,
          N_VM=N_VM)
  
  Prevalence_H=0 # prevalence in humans
  Prevalence_M=0.97 # prevalence in macaques
  Prevalence_VC=0.039 # prevalence of vectors in forest
  Prevalence_VE=0.027 # prevalence of vectors in edge
  Prevalence_VM=0.0026 # prevalence of vectors in deforest
  
  x0 <- round(c(S_H=N_H*(1-Prevalence_H), # susceptible humans
                I_H=N_H*Prevalence_H, # infected humans
                R_H=0, # recovered humans
                S_M=N_M*(1-Prevalence_M), # susceptible macaques
                I_M=N_M*Prevalence_M, # infected macaques
                R_M=0, # recovered macaques
                S_VC=N_VC*(1-Prevalence_VC), # susceptible vectors in forest
                E_VC=N_VC*0, # exposed vectors in forest
                I_VC=N_VC*Prevalence_VC, # infected vectors in forest
                S_VE=N_VE*(1-Prevalence_VE), # susceptible vectors in edge
                E_VE=N_VE*0, # exposed vectors in edge
                I_VE=N_VE*Prevalence_VE, # infected vectors in edge
                S_VM=N_VM*(1-Prevalence_VM), # susceptible vectors in deforest
                E_VM=N_VM*0, # exposed vectors in deforest
                I_VM=N_VM*Prevalence_VM # infected vectors in deforest
  ))
  
  a <- c(	"((alpha_HE*C_VH*I_VE/N_HE) + (alpha_HM*C_VH*I_VM/N_HM))*S_H",	"gamma_H*I_H",	"mu_H*S_H",	"mu_H*I_H",	"mu_H*R_H",	"mu_H*(S_H+I_H+R_H)",	
          "((alpha_MC*C_VM*I_VC/N_MC) + (alpha_ME*C_VM*I_VE/N_ME))*S_M",	"gamma_M*I_M",	"mu_M*S_M",	"mu_M*I_M",	"mu_M*R_M",	"mu_M*(S_M+I_M+R_M)",	
          "alpha_MC*C_MV*I_M*S_VC/N_VC",	"zeta*E_VC",	"mu_V*S_VC",	"mu_V*E_VC",	"mu_V*I_VC",	"mu_V*(S_VC+E_VC+I_VC)",	
          "(alpha_HE*C_HV*I_H + alpha_ME*C_MV*I_M)*S_VE/N_VE",	"zeta*E_VE",	"mu_V*S_VE",	"mu_V*E_VE",	"mu_V*I_VE",	"mu_V*(S_VE+E_VE+I_VE)",	
          "alpha_HM*C_HV*I_H*S_VM/N_VM",	"zeta*E_VM",	"mu_V*S_VM",	"mu_V*E_VM",	"mu_V*I_VM",	"mu_V*(S_VM+E_VM+I_VM)")
  
  nu <- matrix(c(-1,	0 ,	-1,	0 ,	0 ,	1 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,
                 1 ,	-1,	0 ,	-1,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,
                 0 ,	1 ,	0 ,	0 ,	-1,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,
                 0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	-1,	0 ,	-1,	0 ,	0 ,	1 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,
                 0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	1 ,	-1,	0 ,	-1,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,
                 0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	1 ,	0 ,	0 ,	-1,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,
                 0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	-1,	0 ,	-1,	0 ,	0 ,	1 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,
                 0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	1 ,	-1,	0 ,	-1,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,
                 0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	1 ,	0 ,	0 ,	-1,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,
                 0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	-1,	0 ,	-1,	0 ,	0 ,	1 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,
                 0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	1 ,	-1,	0 ,	-1,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,
                 0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	1 ,	0 ,	0 ,	-1,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,
                 0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	-1,	0 ,	-1,	0 ,	0 ,	1 ,
                 0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	1 ,	-1,	0 ,	-1,	0 ,	0 ,
                 0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	1 ,	0 ,	0 ,	-1,	0 ),nrow=15,byrow=TRUE)
  
  FinalI=numeric()
  for(i in 1:iteration){
    out <- ssa(x0,a,nu,parms,tf=tf,method=ssa.d(),simName="Malaria model")
    FinalI[i] <- tail(as.data.frame(out$data),1)$I_H  
  }
  MeanFI=mean(FinalI)
  
  return(list(MeanFI,FinalI))
}
#####

N <- 2^7 # sample size
params <- c("phi","Q","N_H","N_M","N_V")
matrices <- c("A", "B", "AB", "BA")
first <- total <- "azzini" # estimator
order <- "second"
R <- round(N/10)
type <- "percent"
conf <- 0.95

X <- sobol_matrices(matrices=matrices, N=N, params=params, order=order, type="LHS"); dim(X)

X[,1] <- qunif(X[,1], min=0.01, max=0.99)
X[,2] <- qunif(X[,2], min=0.01, max=0.99)
X[,3] <- qunif(X[,3], min=1, max=144)
X[,4] <- qunif(X[,4], min=1, max=144)
X[,5] <- qunif(X[,5], min=1, max=3593)

X=data.frame(X)
y=rep(0, nrow(X))

iteration=100

for(s in 1:nrow(X)){
  ww = with(X[s,], SobolModel(1,iteration,phi,Q,N_H,N_M,N_V))
  y[s] <- ww[[1]]; print(s)}
#save(N, params, matrices, first, total, order, R, type, conf, X, y, file="RDataFiles/SobolSA_tf1.RData")

for(s in 1:nrow(X)){
  ww = with(X[s,], SobolModel(14,iteration,phi,Q,N_H,N_M,N_V))
  y[s] <- ww[[1]]; print(s)}
#save(N, params, matrices, first, total, order, R, type, conf, X, y, file="RDataFiles/SobolSA_tf14.RData")

###################################################################################################################################
# plot_scatter
load(file="RDataFiles/SobolSA_6params_tf1.RData")
XX1=cbind(melt(X),Y=rep(y,dim(X)[2])); XX1$variable <- factor(XX1$variable, levels = c("phi", "Q", "N_H", "N_M", "N_V"))
levels(XX1$variable)[levels(XX1$variable) == "N_H"] <- 'N[H]'; levels(XX1$variable)[levels(XX1$variable) == "N_M"] <- 'N[M]'; levels(XX1$variable)[levels(XX1$variable) == "N_V"] <- 'N[V]'

load(file="RDataFiles/SobolSA_6params_tf14.RData") 
XX14=cbind(melt(X),Y=rep(y,dim(X)[2])); XX14$variable <- factor(XX14$variable, levels = c("phi", "Q", "N_H", "N_M", "N_V"))
levels(XX14$variable)[levels(XX14$variable) == "N_H"] <- 'N[H]'; levels(XX14$variable)[levels(XX14$variable) == "N_M"] <- 'N[M]'; levels(XX14$variable)[levels(XX14$variable) == "N_V"] <- 'N[V]'

XX1=subset(XX1, Y<0.3)
XX=cbind(rbind(XX1,XX14),Day=c(rep(1,nrow(XX1)),rep(14,nrow(XX14))))
XX$Y[XX$Y > 1.5] <- NA
XX$Day <- ifelse(XX$Day==1, '1-day', '14-day')

# Figure 4
ggplot(XX,aes(value,Y)) + geom_hex() + facet_grid(Day~variable,scales="free",labeller = label_parsed)+
  scale_fill_continuous(trans = 'log',breaks=c(1,3,10,30,100),name="Frequency") + labs(x = "", y = "Number of human infections")+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent", color = NA), legend.key = element_rect(fill = "transparent", color = NA), strip.background = element_rect(fill = "white"),strip.text = element_text(size=12),legend.position = "bottom", axis.text.x=element_text(angle=70,hjust=1)) +
  stat_summary_bin(fun = "mean", geom = "point", colour = "red", size = 1.3)

##############################################################################################################
# 1 day
load(file="RDataFiles/SobolSA_6params_tf1.RData"); X1=X; y1=y; dim(X); length(y); log(N,2); iteration
ind1 <- sobol_indices(matrices = matrices, Y = y1, N = N, params = params,first = first, total = total, order = order, boot = TRUE, R = R, parallel = "no", type = type, conf = conf)
cols <- colnames(ind1$results)[1:5]
ind1$results[, (cols) := round(.SD, 3), .SDcols = (cols)]
results<-ind1$results
#ind1.dummy <- sobol_dummy(Y = y, N = N, params = params, boot = TRUE,R = R)
ind1$results$parameters  <- factor(ind1$results$parameters, levels = c("phi", "Q", "N_H", "N_M", "N_V"))
F1 <- plot(ind1) + scale_x_discrete(labels=c(expression(phi,Q,N[H],N[M],N[V])))+labs(x='Parameter')+theme(legend.position='none', axis.text.x=element_text(size=14,face="bold"))+ylim(-0.03,0.86)+ggtitle('(A)')

# Interaction
param_order <- c("phi", "eta", "Q", "N_H", "N_M", "N_V")
param_labels <- c("phi" = expression(phi), "eta" = expression(eta), "Q"="Q","N_H" = expression(N[H]), "N_M" = expression(N[M]), "N_V" = expression(N[V]))
interaction_data <- results[results$sensitivity == "Sij", ]
interaction_data <- dcast(interaction_data, parameters ~ sensitivity, value.var = "original")
interaction_data <- interaction_data %>% separate(parameters, into = c("Parameter1", "Parameter2"), sep = "\\.")
interaction_data$Parameter1 <- factor(interaction_data$Parameter1, levels = names(param_labels))
interaction_data$Parameter2 <- factor(interaction_data$Parameter2, levels = names(param_labels))

fixed_sector_colors <- c("phi" = "#FF6347", "eta" = "#4682B4", "Q" = "#32CD32","N_H"='yellow3', "N_M"="cyan3", "N_V"="magenta3")
fixed_link_colors <- c("#FF000080", "#0000FF80", "#00FF0080","gray","lightyellow2")

# Figure 6
layout(matrix(1:2, nrow = 1))
chordDiagram(interaction_data,grid.col = fixed_sector_colors,col=fixed_link_colors,transparency = 0.5,annotationTrack = "grid",preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {sector.name = get.cell.meta.data("sector.index")
  circos.text(x = get.cell.meta.data("xcenter"),y = get.cell.meta.data("ylim")[1] + 1,labels = param_labels[sector.name],facing = "clockwise",niceFacing = TRUE,adj = c(0, 0.5))}, bg.border = NA)
title(main="(A)")
################################################################################################################################
# 14 days
load(file="RDataFiles/SobolSA_6params_tf14.RData"); X14=X; y14=y; dim(X); length(y); log(N,2); iteration

ind14 <- sobol_indices(matrices = matrices, Y = y14, N = N, params = params,first = first, total = total, order = order, boot = TRUE, R = R, parallel = "no", type = type, conf = conf)
cols <- colnames(ind14$results)[1:5]
ind14$results[, (cols) := round(.SD, 3), .SDcols = (cols)]
results14<-ind14$results
#ind14.dummy <- sobol_dummy(Y = y, N = N, params = params, boot = TRUE,R = R)
ind14$results$parameters  <- factor(ind14$results$parameters, levels = c("phi", "Q", "N_H", "N_M", "N_V"))
F14 <- plot(ind14) + scale_x_discrete(labels=c(expression(phi,Q,N[H],N[M],N[V])))+labs(x='Parameter')+theme(legend.position='bottom', axis.text.x=element_text(size=14,face="bold"))+ylim(-0.03,0.86)+ggtitle("(B)")

interaction_data <- results14[results14$sensitivity == "Sij", ]
interaction_data <- dcast(interaction_data, parameters ~ sensitivity, value.var = "original")
interaction_data <- interaction_data %>% separate(parameters, into = c("Parameter1", "Parameter2"), sep = "\\.")
interaction_data$Parameter1 <- factor(interaction_data$Parameter1, levels = names(param_labels))
interaction_data$Parameter2 <- factor(interaction_data$Parameter2, levels = names(param_labels))

fixed_sector_colors <- c("phi" = "#FF6347", "eta" = "#4682B4", "Q" = "#32CD32","N_H"='yellow3', "N_M"="cyan3", "N_V"="magenta3")
fixed_link_colors <- c("#FF000080", "#0000FF80", "#00FF0080","gray","lightyellow2")
chordDiagram(interaction_data,grid.col = fixed_sector_colors,col=fixed_link_colors,transparency = 0.5,annotationTrack = "grid",preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {sector.name = get.cell.meta.data("sector.index")
  circos.text(x = get.cell.meta.data("xcenter"),y = get.cell.meta.data("ylim")[1] + 1,labels = param_labels[sector.name],facing = "clockwise",niceFacing = TRUE,adj = c(0, 0.5))}, bg.border = NA)
title(main="(B)")

# Figure 5
grid.arrange(F1, F14, nrow=2, heights=c(1,1.2))
