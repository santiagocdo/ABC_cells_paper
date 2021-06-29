#### R-W for ABC paper ####
rm(list = ls()) # errase previous Global Environment, import manually the file

# read csv with simulation parameters
sim <- read.csv(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/conditionsToSimuate_v2.csv"))

# trial types
condNames <- unique(sim$conditionName)
nCondType <- length(condNames)



#### SIMULATION ####
# n simulations per condition 
nSim <- 50

#### free parameters ####
# alpha or learning rate
aRW <- 0.2
# number of blocks
nBlock <- 50

# parameters in one same vector
paramRW <- c(aRW,nBlock)

source("C:/Users/scastiello/OneDrive - Nexus365/1.- Contingency Learning (OBE-TM)/ABC paper/simulations/f_RW_model.R")

for (c in 1:nCondType) {
  message(paste0("Condition: ",condNames[c]))
  
  trialType <- sim[sim$conditionName == condNames[c],]
  
  # number of trials per block
  nTrialBlock <- nrow(trialType)
  
  # Rescorla-Wagner
  INPUT <- cbind(trialType$ctx,trialType$cue)
  OUTPUT <- matrix(trialType[,c("outcome")],ncol=1,nrow=nTrialBlock)
  for (s in 1:nSim) {
    message(paste0("sim: ",s))
    
    aRW <- f_RW_model(paramRW,INPUT,OUTPUT)
    if (s == 1) {
      b <- data.frame(s,rep(1:nBlock,each=nTrialBlock),aRW[,1,])
    } else {
      b <- rbind(b,data.frame(s,rep(1:nBlock,each=nTrialBlock),aRW[,1,]))
    }
  } # end sim
  
  if (c == 1) {
    actO_rw <- data.frame(trialType[1,1],trialType[1,2],trialType[1,3],condNames[c],b)
  } else
    actO_rw <- rbind(actO_rw,data.frame(trialType[1,1],trialType[1,2],trialType[1,3],condNames[c],b))
  
}; rm(b,aRW,INPUT,OUTPUT,nTrialBlock,trialType) # end conditionType

# what do we need?
# Fig A: normal, x frequency or duration
# Fig B: normal and adjusted, x frequency
# Fig C: normal and adjusted, x duration
# Fig D: normal and overadjusted, x frequency


#actO_rw <- actO_rw[,1:(ncol(actO_rw)-1)]
names(actO_rw) <- c(colnames(sim)[1:4],"subj","block","ctx","cue")
if (!require(reshape2)) {install.packages("reshape2")}; library(reshape2) # melt()
actO_rw <- melt(actO_rw, measure.vars = c("ctx","cue")); names(actO_rw) <- c(colnames(actO_rw)[1:6],"cue","V")
actO_rw$cue <- factor(actO_rw$cue, levels = c("cue","ctx"))
# actO_rw$netV <- actO_rw$ctx + actO_rw$cue

for_plot <- actO_rw[actO_rw$block == max(actO_rw$block) & actO_rw$cue == "cue",]
for_plot$condition <- factor(for_plot$condition,levels = c("few","many"))
for_plot$general <- factor(for_plot$general, levels = c("Normal","Adjusted","Overadjusted"))

figAdb <- for_plot[for_plot$general == "Normal",]
figBdb <- rbind(for_plot[for_plot$general == "Normal",],for_plot[for_plot$general == "Adjusted",])
figCdb <- rbind(for_plot[for_plot$general == "Normal",],for_plot[for_plot$general == "Adjusted",])
figDdb <- rbind(for_plot[for_plot$general == "Normal",],for_plot[for_plot$general == "Overadjusted",])


#### visualization #### 3600 x 2400 pixeles
if (!require(ggplot2)) {install.packages("ggplot2")}; library(ggplot2) # ggplot()
if (!require(viridis)) {install.packages("viridis")}; library(viridis) # viridis()

#### figure per paper ####
fig2A <- ggplot(figAdb,aes(x=condition, y=V, shape=condition)) +
  geom_hline(yintercept = 0, color = "gray50", size = 0.5) +
  geom_smooth(method = "lm", aes(group = interaction(cell)), colour = "black", size = 1.3, se = F) +
  labs(title = "E1: Frequency or Duration effect", x = "Frequency or Duration", y = "Associative Strength (V)",
       linetype = "") +
  scale_y_continuous(limits = c(-0.45, 0.45), breaks = c(-0.2, 0, 0.2)) +
  scale_linetype_manual(values = c("solid","dotdash"), labels = c("Non-adjusted","Adjusted or Overadjusted")) +
  facet_wrap(. ~ cell) + 
  theme_classic() + theme(axis.text.x = element_blank(), axis.title.y = element_blank(),
                          axis.ticks.x = element_line(size = 0),
                          axis.ticks.length.x = unit(0, "pt"))

fig2C <- ggplot(figBdb, aes(x=condition, y=V, shape=condition, linetype=general)) +
  geom_hline(yintercept = 0, color = "gray50", size = 0.5) +
  geom_smooth(method = "lm", aes(group = interaction(cell,general)), colour = "black", size = 1.3, se = F) +
  labs(title = "E2 & E4: Adjusted for Duration", x = "Frequency", y = "Associative Strength (V)",
       linetype = "") +
  scale_y_continuous(limits = c(-0.45, 0.45), breaks = c(-0.2, 0, 0.2)) +
  scale_linetype_manual(values = c("solid","dotdash"), labels = c("Non-adjusted","Adjusted or Overadjusted")) +
  facet_wrap(. ~ cell) + 
  theme_classic() + theme(axis.text.x = element_blank(), axis.title.y = element_blank(),
                          axis.ticks.x = element_line(size = 0),
                          axis.ticks.length.x = unit(0, "pt"))

fig2B <- ggplot(figCdb, aes(x=condition, y=V, shape=condition, linetype=general)) +
  geom_hline(yintercept = 0, color = "gray50", size = 0.5) +
  geom_smooth(method = "lm", aes(group = interaction(cell,general)), colour = "black", size = 1.3, se = F) +
  labs(title = "E3: Adjusted for Frequency", x = "Duration", y = "Associative Strength (V)",
       linetype = "") +
  scale_y_continuous(limits = c(-0.45, 0.45), breaks = c(-0.2, 0, 0.2)) +
  scale_linetype_manual(values = c("solid","dotdash"), labels = c("Non-adjusted","Adjusted or Overadjusted")) +
  facet_wrap(. ~ cell) + 
  theme_classic() + theme(axis.text.x = element_blank(), axis.title.y = element_blank(),
                          axis.ticks.x = element_line(size = 0),
                          axis.ticks.length.x = unit(0, "pt"))

fig2D <- ggplot(figDdb, aes(x=condition, y=V, shape=condition, linetype=general)) +
  geom_hline(yintercept = 0, color = "gray50", size = 0.5) +
  geom_smooth(method = "lm", aes(group = interaction(cell,general)), colour = "black", size = 1.3, se = F) + 
  labs(title = "E5: Overadjusted for Duration ", x = "Frequency", y = "Associative Strength (V)",
       linetype = "") +
  scale_y_continuous(limits = c(-0.45, 0.45), breaks = c(-0.2, 0, 0.2)) +
  scale_linetype_manual(values = c("solid","dotdash"), labels = c("Non-adjusted","Adjusted or Overadjusted")) +
  facet_wrap(. ~ cell) + 
  theme_classic() + theme(axis.text.x = element_blank(), axis.title.y = element_blank(),
                          axis.ticks.x = element_line(size = 0),
                          axis.ticks.length.x = unit(0, "pt"))

if (!require(ggpubr)) {install.packages("ggpubr")}; library(ggpubr) # annotate_figure()
fig2 <- annotate_figure(ggarrange(fig2A,fig2C,fig2B,fig2D,nrow = 2,ncol = 2,labels = c("a","b","c","d"),common.legend = T),
                left =  text_grob("Associative Strength (V)", size = 15, face = "bold", rot = 90, x = 0.5, y = 0.5),
                top =  text_grob("Simulations", size = 15, face = "bold", rot = 0, x = 0.5, y = 0.5))



ggsave(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/figure2.png"),
       plot = fig2, width = 18, height = 18, dpi = 900, units = "cm", limitsize = T)

