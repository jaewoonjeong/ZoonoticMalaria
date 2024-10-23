library(GillespieSSA)

S_model <- function(tf, phi, scenario){
  Q=0.5
  eta=0.17
  epsilon <- 5.88*phi-8.62*phi^2+2.73*phi^3
  
  Area_matrix=phi
  Area_core=1-phi
  
  Area_edge <- epsilon*eta # edge area
  Area_deforest=(1-Area_edge)*Area_matrix/(Area_matrix+Area_core) # deforest area
  Area_forest=(1-Area_edge)*Area_core/(Area_matrix+Area_core) # forest area
  
  Full_pop=144 # maximum host population size
  
  if(scenario==1){
    N_H=144 
    N_M=144
  } else if(scenario==2){
    N_H=144
    N_M=ceiling(Full_pop*(Area_edge + Area_forest)) # number of macaques (Assumption)    
  } else if(scenario==3){
    N_H <- ceiling(Full_pop*(Area_edge + Area_deforest)) # number of humans (Assumption)
    N_M=144
  } else {
    N_H <- ceiling(Full_pop*(Area_edge + Area_deforest)) # number of humans (Assumption)
    N_M <- ceiling(Full_pop*(Area_edge + Area_forest)) # number of macaques (Assumption)
  }
  
  N_HC=N_H*.0 # number of human individuals in forest 
  N_HE=N_H*Area_edge/(Area_edge + Area_deforest) # number of human individuals in edge
  N_HM=N_H*Area_deforest/(Area_edge + Area_deforest) # number of human individuals in deforest
  
  N_MC=N_M*Area_forest/(Area_forest + Area_edge) # number of macaque individuals in forest
  N_ME=N_M*Area_edge/(Area_forest + Area_edge) # number of macaque individuals in edge
  N_MM=N_M*.0 # number of macaque individuals in deforest
  
  N_V=3593
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
  
  out <- ssa(x0,a,nu,parms,tf=tf,method=ssa.d(),simName="Malaria model")
  FinalI <- tail(as.data.frame(out$data),1)$I_H  
  
  return(FinalI)
}

####################################################
iteration=2000

N=15 # how many values per parameter
phi_variation=seq(0.01, 0.99, length=N)

tf=c(1,7,14) # simulation duration in days

# Scenario 1: constant human pop & constant macaque pop
# Scenario 2: constant human pop & area-proportional macaque pop
# Scenario 3: area-proportional human pop & constant macaque pop
# Scenario 4: area-proportional human pop & area-proportional macaque pop

Output=matrix(NA,iteration,length(tf)*length(phi_variation)*4); dim(Output)=c(iteration,length(phi_variation),length(tf),4)
for(i in 1:iteration){for(p in 1:length(phi_variation)){for(t in 1:length(tf)){for(s in 1:4){
  Output[i,p,t,s]=S_model(tf[t],phi_variation[p],s)
}}};print(i)}

#save(Output,N,tf,iteration,phi_variation, file='RDataFiles/Output_scenarios.RData')
#load(file='RDataFiles/Output_scenarios.RData')

# Define sequence for phi values
phi <- seq(0.01, 0.99, length = N)

# Initialize vectors for storing human (N_H) and macaque (N_M) populations
N_H <- numeric(length(phi))
N_M <- numeric(length(phi))

# Loop through each phi value and calculate corresponding populations
for (p in seq_along(phi)) {
  # Compute epsilon based on phi
  epsilon <- 5.88 * phi[p] - 8.62 * phi[p]^2 + 2.73 * phi[p]^3
  
  # Define areas for matrix and core
  Area_matrix <- phi[p]
  Area_core <- 1 - phi[p]
  
  # Constants
  eta <- 0.17  # Edge effect constant
  
  # Calculate areas
  Area_edge <- epsilon * eta                 # Edge area
  Area_deforest <- (1 - Area_edge) * Area_matrix / (Area_matrix + Area_core)  # Deforested area
  Area_forest <- (1 - Area_edge) * Area_core / (Area_matrix + Area_core)      # Forest area
  
  # Full population size assumption
  Full_pop <- 144
  
  # Calculate populations based on areas
  N_H[p] <- ceiling(Full_pop * (Area_edge + Area_deforest))  # Human population
  N_M[p] <- ceiling(Full_pop * (Area_edge + Area_forest))    # Macaque population
}

# Load required package
library(reshape2)

# Compute the mean of 'Output' array along the specified dimensions and reshape it into a data frame
df <- melt(apply(Output, c(2, 3, 4), mean))

# Rename columns for better clarity
colnames(df) <- c('Phi', 'Time', 'Scenario', 'Infected')

# Update the 'Phi' and 'Time' columns with corresponding values from 'phi_variation' and 'tf'
df$Phi <- phi_variation[df$Phi]
df$Time <- tf[df$Time]
df$Scenario <- paste("Scenario", df$Scenario)

# Format the 'Time' column as "days" and set as a factor with unique levels to preserve order
for(i in 1:nrow(df)){df$Time[i]=ifelse(df$Time[i]==1,'1 day',ifelse(df$Time[i]==7,'7 days','14 days'))}
df$Time <- factor(df$Time, levels = unique(df$Time))

# Create 'Human' and 'Macaque' columns based on specific conditions
total_phi_time <- length(phi_variation) * length(tf)
df$Human <- c(rep(144, total_phi_time), rep(144, total_phi_time), rep(N_H, length(tf)), rep(N_H, length(tf)))
df$Macaque <- c(rep(144, total_phi_time), rep(N_M, length(tf)), rep(144, total_phi_time), rep(N_M, length(tf)))

# Calculate infections per person
df$PerPerson <- df$Infected / df$Human

library(ggplot2); library(gridExtra); library(gtable)

## by simulation periods
ggplot(data = df, aes(x = Phi, y = Infected, color = Scenario)) +
  geom_line(linewidth = 1, aes(linetype = Scenario)) +geom_point(size = 2, aes(shape=Scenario)) +
  facet_wrap(. ~ Time, ncol = 1, dir = "h", scales = 'free_y') +
  scale_color_manual(values = c("Scenario 1" = "skyblue3", "Scenario 2" = "blue3", "Scenario 3" = "red3", "Scenario 4" = "orange3")) +
  theme_bw() +
  theme(legend.position = 'bottom', axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(color = NULL, shape = NULL, linetype = NULL, x = expression("Proportion of deforestation (" ~ phi ~ ")"), y = 'Number of human infections')
