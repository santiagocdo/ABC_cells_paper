rm(list = ls()) # erase previous Global Environment, import manually the file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set workspace

# print csv files
write_csv <- 0
# print figures?
print_png <- 0

if(!require(lmerTest)) {install.packages("lmerTest")}; library(lmerTest) # lmer() + step()
if(!require(report)) {install.packages("report")}; library(report) # report()
if(!require(reshape2)) {install.packages("reshape2")}; library(reshape2) # melt()
if(!require(dplyr)) {install.packages("dplyr")}; library(dplyr) # %>%
if(!require(plyr)) {install.packages("plyr")}; library(plyr) # revalue()
if(!require(ggplot2)) {install.packages("ggplot2")}; library(ggplot2) # ggplot()
if(!require(ggpubr)) {install.packages("ggpubr")}; library(ggpubr) # ggarrange()
if(!require(xtable)) {install.packages("xtable")}; library(xtable) # print()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # Experiment 1a and 1b # # # # # # # # # # # # # # # ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

exp01_lf <- read.csv("exp01_lf.csv")
exp01_lf$subject <- paste0(exp01_lf$subject, exp01_lf$location)

# general characteristics see ABCcell_raw_data_pooler.R

######## STATISTICAL ANALYSIS ######## 
# reorder factor order, and change levels order, where baseline is factor #1
exp01_lf$freq_A <- factor(exp01_lf$freq_A, levels = c("36", "9", "144"))
exp01_lf$freq_B <- factor(exp01_lf$freq_B, levels = c("36", "9", "144"))
exp01_lf$freq_C <- factor(exp01_lf$freq_C, levels = c("36", "9", "144"))
exp01_lf$freq_D <- factor(exp01_lf$freq_D, levels = c("36", "9", "144"))
exp01_lf$dur_A <- factor(exp01_lf$dur_A, levels = c("800", "200", "3200"))
exp01_lf$dur_B <- factor(exp01_lf$dur_B, levels = c("800", "200", "3200"))
exp01_lf$dur_C <- factor(exp01_lf$dur_C, levels = c("800", "200", "3200"))
exp01_lf$dur_D <- factor(exp01_lf$dur_D, levels = c("800", "200", "3200"))


plotANDstat_all <- 1
# if == mt, then experiment 1b, elseif == nz, then experiment 1a, otherwise, pooled experiment 1.
exp01_lf <- exp01_lf[exp01_lf$location == "mt",]
### exclude data for visualization
if (plotANDstat_all == 1) {
  # plot and analysis with 0.
  lf_noB <- exp01_lf[!is.na(exp01_lf$rating),]
  
  # graphs data
  limits_vec_fre <- c(0,150)
  break_vec_fre <- c(9,36,144)
  limits_vec_dur <- c(0,3300)
  break_vec_dur <- c(200,800,3200)
} else {
  # plot and analysis without 0.
  lf_noB <- exp01_lf[!is.na(exp01_lf$rating) & exp01_lf$condition != "baseline",]
  
  # graphs data
  limits_vec_fre <- c(0,150)
  break_vec_fre <- c(9,36,144)
  limits_vec_dur <- c(0,3300)
  break_vec_dur <- c(200,800,3200)
}
###

# for non convergence models
#relgrad <- with(mFreq@optinfo$derivs,solve(Hessian,gradient))
#max(abs(relgrad))

mFull <- lmer(rating ~ 1 + freq_A + freq_B + freq_C + freq_D + dur_A + dur_B + dur_C + dur_D + 
                (1 | subject), 
              #control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)),
              REML = FALSE, data = lf_noB); anova(mFull)
mFreq <- lmer(rating ~ 1 + freq_A + freq_B + freq_C + freq_D + (1 | subject), 
              control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)),
              REML = FALSE, data = lf_noB); anova(mFreq)
mDur <- lmer(rating ~ 1 + dur_A + dur_B + dur_C + dur_D + (1 | subject), 
             #control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)),
             REML = FALSE, data = lf_noB); anova(mDur)



# Full against reduced
anova(mFull,mFreq, test="LRT"); exp((BIC(mFull) - BIC(mFreq))/2); AIC(mFull) - AIC(mFreq)
anova(mFull,mDur,test="LRT"); exp((BIC(mFull) - BIC(mDur))/2); AIC(mFull) - AIC(mDur)

# Freq against Dur
anova(mDur,mFreq,test="LRT"); exp((BIC(mDur) - BIC(mFreq))/2); AIC(mDur) - AIC(mFreq)

BIC(mFull,mFreq,mDur)
AIC(mFull,mFreq,mDur)

# print table with LMM
csv_lmm <- cbind(c(rep("Full", 17),rep("Freq", 9),rep("Dur", 9)),
                 rbind(table_long(report(mFull))[1:17,],
                       table_long(report(mFreq))[1:9,],
                       table_long(report(mDur))[1:9,]))
colnames(csv_lmm) <- c("Model",colnames(csv_lmm)[2:14])

# repeat previous code for NZ (exp 1a), MT (exp 1b), and for overall 1.
csv_lmm_all <- csv_lmm
csv_lmm_nz <- csv_lmm
csv_lmm_mt <- csv_lmm

if (write_csv == 1) {
  write.csv(csv_lmm_all,"C:/Users/scastiello/OneDrive - Nexus365/1.- Contingency Learning (OBE-TM)/ABC paper/stats/exp1_gitHub_LMM_pooled.csv",
            row.names = F)
  write.csv(csv_lmm_nz,"C:/Users/scastiello/OneDrive - Nexus365/1.- Contingency Learning (OBE-TM)/ABC paper/stats/exp1a_gitHub_LMM.csv",
            row.names = F)
  write.csv(csv_lmm_mt,"C:/Users/scastiello/OneDrive - Nexus365/1.- Contingency Learning (OBE-TM)/ABC paper/stats/exp1b_gitHub_LMM.csv",
            row.names = F)
}



######## VISUALIZATION ######## 
# prepare data base for plots
# add baseline for all cells (B, C, and D)
nSubject <- length(unique(lf_noB$subject))
baseline <- rbind(lf_noB[lf_noB$condition == "baseline",], lf_noB[lf_noB$condition == "baseline",], lf_noB[lf_noB$condition == "baseline",])
baseline$cell <- c(rep("B",nSubject),rep("C",nSubject),rep("D",nSubject))

db_plot <- rbind(lf_noB,baseline)
remove(baseline)

# plot size in pixels 1500 x 1000. PDF 15 x 10 (actual paper is 2200 x 1600 pixels)
sf <- 1 # scaling factor
plot_title_size <- 48#24
strip_text_size <- 44#22
axis_title_size <- 40#20
axis_text_size <- 36#18
legend_title_size <- 40#20
legend_text_size <- 36#18

# stratified by place (nz or mt)
# if == mt, then experiment 1b, elseif == nz, then experiment 1a, otherwise, pooled experiment 1.
db_plot <- db_plot#[db_plot$location == "mt", ] 
# change condition name label
db_plot$condition_simple <- as.factor(db_plot$condition_simple)
levels(db_plot$condition_simple) <- c("Baseline","Few","Long","Many","Short")

pos <- position_dodge(10)

fig2left <- ggplot(db_plot[db_plot$effect != "duration",], aes(x = frequency, y = rating, shape = condition_simple)) +
  theme_bw() +
  geom_hline(yintercept = 0, color = "gray50", size = 1.2) +
  geom_line(aes(group = subject), alpha = 0.1, size = 0.8, col = "gray30", position = pos) +
  geom_jitter(width = 2, height = 0.2, alpha = 0.2, size = 1.2, col = "gray30") +
  geom_smooth(method = "lm", aes(group = cell), colour = "black", size = 2, se = F) + 
  stat_summary(fun.data = "mean_cl_boot", size = 2, position = pos) +
  labs(title = "Frequency effect", #subtitle = "Frequency effect",
       shape = "Condition") + xlab("Frequency") + #ylab("Rating") +
  scale_shape_manual(values = c(15, 17, 19), breaks = c("Few","Baseline","Many")) +  
  scale_x_continuous(limits = limits_vec_fre, breaks = break_vec_fre, minor_breaks = NULL) + 
  scale_y_continuous(limits = c(-10.5, 10.5), breaks = c(-10, 0, 10)) +
  facet_wrap(. ~ cell) + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = sf*strip_text_size),
        plot.title = element_text(size = sf*plot_title_size, hjust = 0,face = "bold"),
        plot.subtitle = element_text(size = sf*plot_title_size/2, hjust = 0, face = "italic", color = "black"),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = sf*axis_text_size, colour="black"),
        axis.title.x = element_text(size = sf*axis_title_size, face = "bold"),
        axis.text.x = element_text(size = sf*axis_text_size, colour="black"),
        axis.line = element_line(colour = 'black', size = 1.6),
        axis.ticks = element_line(size = 1.6),
        axis.ticks.length = unit(10, "pt"),
        legend.title = element_text(size = sf*legend_title_size), 
        legend.text = element_text(size = sf*legend_text_size),
        legend.position = "bottom",
        legend.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),
        panel.border = element_rect(color = "black"),
        panel.spacing = unit(1, "lines")
  )
#fig2left

pos <- position_dodge(120)

fig2right <- ggplot(db_plot[db_plot$effect != "frequency",], aes(x = duration, y = rating, shape = condition_simple)) +
  theme_bw() +
  geom_hline(yintercept = 0, color = "gray50", size = 1.2) +
  geom_line(aes(group = subject), alpha = 0.1, size = 0.8, col = "gray30", position = pos) +
  geom_jitter(width = 2, height = 0.2, alpha = 0.2, size = 1.2, col = "gray30") +
  geom_smooth(method = "lm", aes(group = cell), colour = "black", size = 2, se = F) + 
  stat_summary(fun.data = "mean_cl_boot", size = 2, position = pos) +
  labs(title = "Duration effect", #subtitle = "Duration effect effect",
       shape = "Condition") + xlab("Duration (ms)") + #ylab("Rating") +
  scale_shape_manual(values = c(15, 17, 19), breaks = c("Short","Baseline","Long")) +
  scale_x_continuous(limits = limits_vec_dur, breaks = break_vec_dur, minor_breaks = NULL) + 
  scale_y_continuous(limits = c(-10.5, 10.5), breaks = c(-10, 0, 10)) +
  facet_wrap(. ~ cell) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = sf*strip_text_size),
        plot.title = element_text(size = sf*plot_title_size, hjust = 0,face = "bold"),
        plot.subtitle = element_text(size = sf*plot_title_size/2, hjust = 0, face = "italic", color = "black"),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = sf*axis_text_size, colour="black"),
        axis.title.x = element_text(size = sf*axis_title_size, face = "bold"),
        axis.text.x = element_text(size = sf*axis_text_size*0.8, colour="black"),
        axis.line = element_line(colour = 'black', size = 1.6),
        axis.ticks = element_line(size = 1.6),
        axis.ticks.length = unit(10, "pt"),
        legend.title = element_text(size = sf*legend_title_size), 
        legend.text = element_text(size = sf*legend_text_size),
        legend.position = "bottom",
        legend.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),
        panel.border = element_rect(color = "black"),
        panel.spacing = unit(1, "lines")
  )
#fig2right

results_plot <- annotate_figure(
  ggarrange(fig2left, fig2right, ncol = 2),
  left = text_grob("Rating",
                   color = "black", 
                   face = "bold", 
                   size = plot_title_size, 
                   rot = 90)
  # top = text_grob(paste("Pooled data (NZ+MT); n =",length(levels(as.factor(db_plot$subject)))),
  #                 color = "black",
  #                 face = "bold", 
  #                 size = plot_title_size*1.2)
)
# ggpubr::ggarrange(p_freq,p_dur, plotlist = NULL, ncol = 2, nrow = NULL,
#                   labels = NULL, label.x = 0, label.y = 1, hjust = -0.5,
#                   vjust = 1.5, font.label = list(size = 52, color = "black", face = "bold", family = NULL), 
#                   align = c("none", "h", "v", "hv"), widths = c(1,0.87), heights = 1, legend = "bottom", common.legend = TRUE)
# 2200 x 1600 pixels
results_plot

if (print_png == 1) {
  ggsave("C:/Users/scastiello/OneDrive - Nexus365/1.- Contingency Learning (OBE-TM)/ABC paper/figures and tables/fig2_github.png",
         plot = results_plot, height = 16, width = 22)
  ggsave("C:/Users/scastiello/OneDrive - Nexus365/1.- Contingency Learning (OBE-TM)/ABC paper/figures and tables/fig3_github.png",
         plot = results_plot, height = 16, width = 22)
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # Experiment 2 # # # # # # # # # # # # # # # # # # # ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

exp02_lf <- read.csv("exp02_lf.csv")

# general characteristics
nSubj <- length(unique(exp02_lf$subject)); nSubj
sex <- colSums(table(exp02_lf$subject,exp02_lf$sex)>0); sex
meanAge <- mean(exp02_lf$age[!duplicated(exp02_lf$subject)], na.rm = T); meanAge
sdAge <- sd(exp02_lf$age[!duplicated(exp02_lf$subject)], na.rm = T); sdAge
rangeAge <- range(exp02_lf$age[!duplicated(exp02_lf$subject)], na.rm = T); rangeAge

######## STATISTICAL ANALYSIS ######## 
# reorder factor order, and change levels order, where baseline is factor #1
exp02_lf$frequency <- factor(exp02_lf$frequency, levels = c("36", "9", "144"))
exp02_lf$duration <- factor(exp02_lf$duration, levels = c("800", "200", "3200"))



### filter data and created plot limits vectors
# plot and analysis with 0.
lf_noB <- exp02_lf[!is.na(exp02_lf$rating),]

# graphs data
break_vec_fre <- c(9,36,144)
limits_vec_fre <- c(0,144+break_vec_fre[3]*0.05)
###

# for non convergence models
# if max(abs(relgrad)) is <0.001 then things might be ok
# relgrad <- with(mFreq@optinfo$derivs,solve(Hessian,gradient)); max(abs(relgrad))

# FULL BASELINE MODEL 
mFull <- lmer(rating ~ (frequency * adjusted * cell) + cond_order + (1 | subject),
              #control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = F),
              REML = F, data = lf_noB); anova(mFull)

# FULL NO-BASELINE MODEL
mFull_nb <- lmer(rating ~ (frequency * adjusted * cell) + cond_order + (1 | subject),
                 #control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = F),
                 REML = F, data = lf_noB[lf_noB$condition != "Baseline",]); anova(mFull_nb)

# REDUCED BASELINE MODEL
mFull_red <- lmer(rating ~ (frequency * cell) + cond_order + (1 | subject),
                  #control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = F),
                  REML = F, data = lf_noB); anova(mFull_red)

# REDUCED NO-BASELINE MODEL
mFull_nb_red <- lmer(rating ~ frequency * cell + cond_order + (1 | subject),
                     #control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = F),
                     REML = F, data = lf_noB[lf_noB$condition != "Baseline",]); anova(mFull_nb_red)


# Full against reduced
anova(mFull,mFull_red, test="LRT"); exp((BIC(mFull) - BIC(mFull_red))/2); AIC(mFull) - AIC(mFull_red)
anova(mFull_nb,mFull_nb_red, test="LRT"); exp((BIC(mFull_nb) - BIC(mFull_nb_red))/2); AIC(mFull_nb) - AIC(mFull_nb_red)


## post hoc exploration ##
mA <- lmer(rating ~ (frequency * adjusted) + cond_order + (1 | subject), 
           REML = F, data = lf_noB[lf_noB$cell == "A",]); anova(mA)
mA_red <- lmer(rating ~ (frequency) + cond_order + (1 | subject), 
               REML = F, data = lf_noB[lf_noB$cell == "A",]); anova(mA_red)
anova(mA,mA_red, test="LRT"); exp((BIC(mA) - BIC(mA_red))/2); AIC(mA) - AIC(mA_red)

mB <- lmer(rating ~ (frequency * adjusted) + cond_order + (1 | subject), 
           REML = F, data = lf_noB[lf_noB$cell == "B",]); anova(mB)
mB_red <- lmer(rating ~ (frequency) + cond_order + (1 | subject), 
               REML = F, data = lf_noB[lf_noB$cell == "B",]); anova(mB_red)
anova(mB,mB_red, test="LRT"); exp((BIC(mB) - BIC(mB_red))/2); AIC(mB) - AIC(mB_red)

mC <- lmer(rating ~ (frequency * adjusted) + cond_order + (1 | subject), 
           REML = F, data = lf_noB[lf_noB$cell == "C",]); anova(mC)
mC_red <- lmer(rating ~ (frequency) + cond_order + (1 | subject), 
               REML = F, data = lf_noB[lf_noB$cell == "C",]); anova(mC_red)
anova(mC,mC_red, test="LRT"); exp((BIC(mC) - BIC(mC_red))/2); AIC(mC) - AIC(mC_red)





csv_print <- cbind("duplicated baseline",table_long(report(mFull))[1:13,]); colnames(csv_print) <- c("model",colnames(csv_print)[2:8])
csv_print_mod <- table_long(report(mFull))
csv_print_nb <- cbind("no baseline",table_long(report(mFull_nb))[1:13,]); colnames(csv_print_nb) <- c("model",colnames(csv_print_nb)[2:8])
to_print <- rbind(csv_print,csv_print_nb)
to_print_posthoc <- cbind(c(rep("A",4),rep("B",4),rep("C",4)),
                          rbind(table_long(report(anova(mA)))[1:4,],
                                table_long(report(anova(mB)))[1:4,],
                                table_long(report(anova(mC)))[1:4,]))

if (write_csv == 1) {
  write.csv(csv_print_mod,"C:/Users/scastiello/OneDrive - Nexus365/1.- Contingency Learning (OBE-TM)/ABC paper/stats/exp2_github_LMM.csv",
            row.names = F)
  write.csv(to_print_posthoc,"C:/Users/scastiello/OneDrive - Nexus365/1.- Contingency Learning (OBE-TM)/ABC paper/stats/exp2_github_posthoc_LMM.csv",
            row.names = F)
}


######## VISUALIZATION ######## 
# reorder factor data
lf_noB$frequency <- as.integer(as.character(factor(lf_noB$frequency, levels = c("9","36","144"))))
lf_noB$Adjusted <- lf_noB$adjusted
lf_noB$Adjusted <- as.factor(lf_noB$Adjusted)
levels(lf_noB$Adjusted) <- c("none", "with duration")

# 5000 and 4000 with sf = 4.  With 1.25 times larger (5000/4000)
sf <- 1 # scaling factor. paper = 4
plot_title_size <- 42#24
strip_text_size <- 44#22
axis_title_size <- 40#20
axis_text_size <- 36#18
legend_title_size <- 40#20
legend_text_size <- 36#18
pos <- position_dodge(limits_vec_fre[2]*0.1)


fig4sup <- ggplot(lf_noB, aes(x = frequency, y = rating, col = Adjusted, shape = Adjusted, fll = Adjusted, group = Adjusted)) +
  theme_bw() +
  geom_hline(yintercept = 0, color = "gray50", size = sf*1.2) +
  geom_line(aes(group = interaction(subject, Adjusted)), alpha = 0.05, size = sf*0.8) +
  geom_point(alpha = 0.1, size = sf*1.2, aes(y = jitter(rating, 0, 0.25), shape = Adjusted, fill = Adjusted), position = pos) +
  #geom_jitter(width = 2, height = 0.2, alpha = 0.2, size = 1.2, col = "gray30") +
  geom_smooth(method = "lm", size = sf*2, se = F, aes(linetype = Adjusted)) +
  stat_summary(fun.data = "mean_cl_boot", size = sf*2, position = pos, stroke = sf*2, aes(shape = Adjusted, fill = Adjusted)) + # aes(fill = Adjusted)
  labs(title = "Exp 2: Frequency and adjusted for duration") + #subtitle = "Same baseline for Adjused and Cells") +
  xlab("Frequency") + ylab("Rating") +
  scale_shape_manual(values = c(19, 24), breaks = c("none", "with duration")) +
  scale_fill_manual(values = c("black", "white"), breaks = c("none", "with duration")) +
  #scale_linetype_manual(values=c("solid", "dotted")) +
  scale_colour_manual(values=c("#440154FF","#218F8DFF")) +
  scale_x_continuous(limits = limits_vec_fre, breaks = break_vec_fre, minor_breaks = NULL) + 
  scale_y_continuous(limits = c(-10.5, 10.5), breaks = c(-10, 0, 10)) +
  facet_grid(. ~ cell) + 
  coord_fixed(ratio = limits_vec_fre[2]/12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = sf*strip_text_size),
        plot.title = element_text(size = sf*plot_title_size, hjust = 0,face = "bold"),
        plot.subtitle = element_text(size = sf*plot_title_size/2, hjust = 0, face = "italic", color = "black"),
        plot.margin = unit(c(0.2,1,0.2,0.2), "cm"),
        axis.title.y = element_text(size = sf*axis_title_size, face = "bold"),
        axis.text.y = element_text(size = sf*axis_text_size, colour="black"),
        axis.title.x = element_text(size = sf*axis_title_size, face = "bold"),
        axis.text.x = element_text(size = sf*axis_text_size, colour="black"),
        axis.line = element_line(colour = 'black', size = sf*1.6),
        axis.ticks = element_line(size = sf*1.6),
        axis.ticks.length = unit(sf*10, "pt"),
        legend.title = element_text(size = sf*legend_title_size), 
        legend.text = element_text(size = sf*legend_text_size),
        legend.position = "bottom",
        legend.margin = margin(0.2, 0.2, 1, 0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),
        panel.border = element_rect(color = "black"),
        panel.spacing = unit(sf*1, "lines")
  )
fig4sup

if (print_png == 1) {
  ggsave("C:/Users/scastiello/OneDrive - Nexus365/1.- Contingency Learning (OBE-TM)/ABC paper/figures and tables/Fig4sup_github.png",
         plot = fig4sup, height = 16, width = 22)
}



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # Experiment 3 # # # # # # # # # # # # # # # # # # # ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

exp03_lf <- read.csv("exp03_lf.csv")
exp03_wf <- read.csv("exp03_wf.csv")

# general characteristics
nSubj <- nrow(exp03_wf); nSubj
sex <- table(exp03_wf$sex); sex
meanAge <- mean(2019-as.integer(as.character(exp03_wf$birthYear[substr(exp03_wf$birthYear,1,1)==1]))); meanAge
sdAge <- sd(2019-as.integer(as.character(exp03_wf$birthYear[substr(exp03_wf$birthYear,1,1)==1]))); sdAge
rangeAge <- range(2019-as.integer(as.character(exp03_wf$birthYear[substr(exp03_wf$birthYear,1,1)==1]))); rangeAge

######## STATISTICAL ANALYSIS ######## 
# reorder factor order, and change levels order, where baseline is factor #1
exp03_lf$freq <- factor(exp03_lf$freq, levels = c("12", "4", "36"))
exp03_lf$dur <- factor(exp03_lf$dur, levels = c("1800", "600", "5400"))



### filter data and created plot limits vectors
# plot and analysis with 0.
lf_noB <- exp03_lf[!is.na(exp03_lf$rating),]

# graphs data
break_vec_dur <- c(600,1800,5400)
limits_vec_dur <- c(0,5400+break_vec_dur[3]*0.05)
###

# for non convergence models
# if max(abs(relgrad)) is <0.001 then things might be ok
# relgrad <- with(mFreq@optinfo$derivs,solve(Hessian,gradient)); max(abs(relgrad))

# FULL BASELINE MODEL
mFull <- lmer(rating ~ (dur * adjusted * cell) + cond_order + (1 | subject),
              #control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = F),
              REML = F, data = lf_noB); anova(mFull)

# FULL NO-BASELINE MODEL
mFull_nb <- lmer(rating ~ (dur * adjusted * cell) + cond_order + (1 | subject),
                 #control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = F),
                 REML = F, data = lf_noB[lf_noB$condition != "baseline",]); anova(mFull_nb)

# REDUCED BASELINE MODEL
mFull_red <- lmer(rating ~ (dur * cell) + cond_order + (1 | subject),
                  #control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = F),
                  REML = F, data = lf_noB); anova(mFull_red)

# REDUCED NO-BASELINE MODEL
mFull_nb_red <- lmer(rating ~ (dur * cell) + cond_order + (1 | subject),
                     #control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = F),
                     REML = F, data = lf_noB[lf_noB$condition != "baseline",]); anova(mFull_nb_red)


# Full against reduced
anova(mFull,mFull_red, test="LRT"); exp((BIC(mFull) - BIC(mFull_red))/2); AIC(mFull) - AIC(mFull_red)
anova(mFull_nb,mFull_nb_red, test="LRT"); exp((BIC(mFull_nb) - BIC(mFull_nb_red))/2); AIC(mFull_nb) - AIC(mFull_nb_red)


## post hoc exploration ##
mA <- lmer(rating ~ (dur * adjusted) + cond_order + (1 | subject), 
           REML = F, data = lf_noB[lf_noB$cell == "A",]); anova(mA)
mA_red <- lmer(rating ~ (dur) + cond_order + (1 | subject), 
               REML = F, data = lf_noB[lf_noB$cell == "A",]); anova(mA_red)
anova(mA,mA_red, test="LRT"); exp((BIC(mA) - BIC(mA_red))/2); AIC(mA) - AIC(mA_red)

mB <- lmer(rating ~ (dur * adjusted) + cond_order + (1 | subject), 
           REML = F, data = lf_noB[lf_noB$cell == "B",]); anova(mB)
mB_red <- lmer(rating ~ (dur) + cond_order + (1 | subject), 
               REML = F, data = lf_noB[lf_noB$cell == "B",]); anova(mB_red)
anova(mB,mB_red, test="LRT"); exp((BIC(mB) - BIC(mB_red))/2); AIC(mB) - AIC(mB_red)

mC <- lmer(rating ~ (dur * adjusted) + cond_order + (1 | subject), 
           REML = F, data = lf_noB[lf_noB$cell == "C",]); anova(mC)
mC_red <- lmer(rating ~ (dur) + cond_order + (1 | subject), 
               REML = F, data = lf_noB[lf_noB$cell == "C",]); anova(mC_red)
anova(mC,mC_red, test="LRT"); exp((BIC(mC) - BIC(mC_red))/2); AIC(mC) - AIC(mC_red)

mD <- lmer(rating ~ (dur * adjusted) + cond_order + (1 | subject), 
           REML = F, data = lf_noB[lf_noB$cell == "D",]); anova(mD)
mD_red <- lmer(rating ~ (dur) + cond_order + (1 | subject), 
               REML = F, data = lf_noB[lf_noB$cell == "D",]); anova(mD_red)
anova(mD,mD_red, test="LRT"); exp((BIC(mD) - BIC(mD_red))/2); AIC(mD) - AIC(mD_red)





csv_print <- cbind("duplicated baseline",table_long(report(mFull)))
csv_print_mod <- table_long(report(mFull))
csv_print_nb <- cbind("no baseline",table_long(report(mFull_nb)))
to_print <- rbind(csv_print,csv_print_nb)
to_print_posthoc <- cbind(c(rep("A",4),rep("B",4),rep("C",4),rep("D",4)),
                          rbind(table_long(report(mA))[1:4,],
                                table_long(report(mB))[1:4,],
                                table_long(report(mC))[1:4,],
                                table_long(report(mD))[1:4,]))

if (write_csv == 1) {
  write.csv(csv_print_mod,"C:/Users/scastiello/OneDrive - Nexus365/1.- Contingency Learning (OBE-TM)/ABC paper/stats/exp3_github_LMM.csv",
            row.names = F)
  write.csv(to_print_posthoc,"C:/Users/scastiello/OneDrive - Nexus365/1.- Contingency Learning (OBE-TM)/ABC paper/stats/exp3_github_posthoc_LMM.csv",
            row.names = F)
}



######## VISUALIZATION ######## 
# reorder factor data
lf_noB$freq <- as.integer(as.character(factor(lf_noB$freq, levels = c("4","12","36"))))
lf_noB$dur <- as.integer(as.character(factor(lf_noB$dur, levels = c("600","1800","5400"))))
lf_noB$Adjusted <- lf_noB$adjusted
lf_noB$Adjusted <- as.factor(lf_noB$Adjusted)
levels(lf_noB$Adjusted) <- c("none", "with frequency")

# x and 4000 with sf = 4.  With 1.6 times larger (40/25)
sf <- 1 # scaling factor. paper = 4
plot_title_size <- 42#24
strip_text_size <- 44#22
axis_title_size <- 40#20
axis_text_size <- 36#18
legend_title_size <- 40#20
legend_text_size <- 36#18
pos <- position_dodge(limits_vec_dur[2]*0.1)


fig4inf <- ggplot(lf_noB, aes(x = dur, y = rating, col = Adjusted, shape = Adjusted, fll = Adjusted, group = Adjusted)) +
  theme_bw() +
  geom_hline(yintercept = 0, color = "gray50", size = sf*1.2) +
  geom_line(aes(group = interaction(subject, Adjusted)), alpha = 0.05, size = sf*0.8) +
  geom_point(alpha = 0.1, size = sf*1.2, aes(y = jitter(rating, 0, 0.25), shape = Adjusted, fill = Adjusted), position = pos) +
  #geom_jitter(width = 2, height = 0.2, alpha = 0.2, size = 1.2, col = "gray30") +
  geom_smooth(method = "lm", size = sf*2, se = F, aes(linetype = Adjusted)) + 
  stat_summary(fun.data = "mean_cl_boot", size = sf*2, position = pos, stroke = sf*2, aes(shape = Adjusted, fill = Adjusted)) +
  labs(title = "Exp 3: Duration and adjusted for frequency") + #subtitle = "Same baseline for Adjused and Cells") + 
  xlab("Duration (ms)") + ylab("Rating") +
  scale_fill_manual(values = c("black", "white"), breaks = c("none", "with frequency")) +
  scale_shape_manual(values = c(19, 24), breaks = c("none", "with frequency")) +
  #scale_linetype_manual(values=c("solid", "dotted")) +
  scale_colour_manual(values=c("#440154FF","#218F8DFF")) +
  scale_x_continuous(limits = limits_vec_dur, breaks = break_vec_dur, minor_breaks = NULL) + 
  scale_y_continuous(limits = c(-10.5, 10.5), breaks = c(-10, 0, 10)) +
  facet_grid(. ~ cell) + 
  coord_fixed(ratio = limits_vec_dur[2]/12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = sf*strip_text_size),
        plot.title = element_text(size = sf*plot_title_size, hjust = 0,face = "bold"),
        plot.subtitle = element_text(size = sf*plot_title_size/2, hjust = 0, face = "italic", color = "black"),
        plot.margin = unit(c(0.2,1,0.2,0.2), "cm"),
        axis.title.y = element_text(size = sf*axis_title_size, face = "bold"),
        axis.text.y = element_text(size = sf*axis_text_size, colour="black"),
        axis.title.x = element_text(size = sf*axis_title_size, face = "bold"),
        axis.text.x = element_text(size = sf*axis_text_size*0.6, colour="black"),
        axis.line = element_line(colour = 'black', size = sf*1.6),
        axis.ticks = element_line(size = sf*1.6),
        axis.ticks.length = unit(sf*10, "pt"),
        legend.title = element_text(size = sf*legend_title_size), 
        legend.text = element_text(size = sf*legend_text_size),
        legend.position = "bottom",
        legend.margin = margin(0.2, 0.2, 1, 0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),
        panel.border = element_rect(color = "black"),
        panel.spacing = unit(sf*1, "lines")
  )
fig4inf

if (print_png == 1) {
  ggsave("C:/Users/scastiello/OneDrive - Nexus365/1.- Contingency Learning (OBE-TM)/ABC paper/figures and tables/Fig4inf_github.png",
         plot = fig4inf, height = 16, width = 22)
}



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # Experiment 2 and 3 # # # # # # # # # # # # # # # # ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
fig4sup # 5000 x 4000 pixels
fig4inf # 6400 x 4000 pixels
ggpubr::ggarrange(fig4sup,fig4inf, plotlist = NULL, ncol = NULL, nrow = 2,
                  labels = NULL, label.x = 0, label.y = 1, hjust = -0.5,
                  vjust = 1.5, font.label = list(size = 52, color = "black", face = "bold", family = NULL), 
                  align = c("none", "h", "v", "hv"), widths = c(1.1), heights = 1, legend = "bottom", common.legend = TRUE)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # Experiment 4 # # # # # # # # # # # # # # # # # # # # ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

exp04_lf <- read.csv("exp04_lf.csv")
exp04_wf <- read.csv("exp04_wf.csv")

# general characteristics
nSubj <- nrow(exp04_wf); nSubj
sex <- table(exp04_wf$sex); sex
meanAge <- mean(2019-as.integer(as.character(exp04_wf$birthYear[substr(exp04_wf$birthYear,1,1)==1]))); meanAge
sdAge <- sd(2019-as.integer(as.character(exp04_wf$birthYear[substr(exp04_wf$birthYear,1,1)==1]))); sdAge
rangeAge <- range(2019-as.integer(as.character(exp04_wf$birthYear[substr(exp04_wf$birthYear,1,1)==1]))); rangeAge


######## STATISTICAL ANALYSIS ######## 
# reorder factor order, and change levels order, where baseline is factor #1
exp04_lf$frequency <- factor(exp04_lf$frequency, levels = c("16", "4", "64"))
exp04_lf$duration <- factor(exp04_lf$duration, levels = c("1200", "300", "4800"))



### filter data and created plot limits vectors
# plot and analysis with 0.
lf_noB <- exp04_lf[!is.na(exp04_lf$rating),]

# graphs data
break_vec_fre <- c(4,16,64)
limits_vec_fre <- c(0,64+break_vec_fre[3]*0.05)
###

# for non convergence models
# if max(abs(relgrad)) is <0.001 then things might be ok
# relgrad <- with(mFreq@optinfo$derivs,solve(Hessian,gradient)); max(abs(relgrad))

# FULL BASELINE MODEL 
mFull <- lmer(rating ~ (frequency * adjusted * cell) + cond_order + (1 | subject),
              #control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = F),
              REML = F, data = lf_noB); anova(mFull)

# FULL NO-BASELINE MODEL
mFull_nb <- lmer(rating ~ (frequency * adjusted * cell) + cond_order + (1 | subject),
                 #control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = F),
                 REML = F, data = lf_noB[lf_noB$condition != "baseline",]); anova(mFull_nb)

# REDUCED BASELINE MODEL
mFull_red <- lmer(rating ~ (frequency * cell) + cond_order + (1 | subject),
                  #control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = F),
                  REML = F, data = lf_noB); anova(mFull_red)

# REDUCED NO-BASELINE MODEL
mFull_nb_red <- lmer(rating ~ frequency * cell + cond_order + (1 | subject),
                     #control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = F),
                     REML = F, data = lf_noB[lf_noB$condition != "baseline",]); anova(mFull_nb_red)


# Full against reduced
anova(mFull,mFull_red, test="LRT"); exp((BIC(mFull) - BIC(mFull_red))/2); AIC(mFull) - AIC(mFull_red)
anova(mFull_nb,mFull_nb_red, test="LRT"); exp((BIC(mFull_nb) - BIC(mFull_nb_red))/2); AIC(mFull_nb) - AIC(mFull_nb_red)


## post hoc exploration ##
mA <- lmer(rating ~ (frequency * adjusted) + cond_order + (1 | subject), 
           REML = F, data = lf_noB[lf_noB$cell == "A",]); anova(mA)
mA_red <- lmer(rating ~ (frequency) + cond_order + (1 | subject), 
               REML = F, data = lf_noB[lf_noB$cell == "A",]); anova(mA_red)
anova(mA,mA_red, test="LRT"); exp((BIC(mA) - BIC(mA_red))/2); AIC(mA) - AIC(mA_red)

mB <- lmer(rating ~ (frequency * adjusted) + cond_order + (1 | subject), 
           REML = F, data = lf_noB[lf_noB$cell == "B",]); anova(mB)
mB_red <- lmer(rating ~ (frequency) + cond_order + (1 | subject), 
               REML = F, data = lf_noB[lf_noB$cell == "B",]); anova(mB_red)
anova(mB,mB_red, test="LRT"); exp((BIC(mB) - BIC(mB_red))/2); AIC(mB) - AIC(mB_red)

mC <- lmer(rating ~ (frequency * adjusted) + cond_order + (1 | subject), 
           REML = F, data = lf_noB[lf_noB$cell == "C",]); anova(mC)
mC_red <- lmer(rating ~ (frequency) + cond_order + (1 | subject), 
               REML = F, data = lf_noB[lf_noB$cell == "C",]); anova(mC_red)
anova(mC,mC_red, test="LRT"); exp((BIC(mC) - BIC(mC_red))/2); AIC(mC) - AIC(mC_red)





csv_print <- cbind("duplicated baseline",table_long(report(mFull))); colnames(csv_print) <- c("model",colnames(csv_print)[2:8])
csv_print_mod <- table_long(report(mFull))
csv_print_nb <- cbind("no baseline",table_long(report(mFull_nb))); colnames(csv_print_nb) <- c("model",colnames(csv_print_nb)[2:8])
to_print <- rbind(csv_print,csv_print_nb)
to_print_posthoc <- cbind(c(rep("A",4),rep("B",4),rep("C",4)),
                          rbind(table_long(report(mA))[1:4,],
                                table_long(report(mB))[1:4,],
                                table_long(report(mC))[1:4,]))

if (write_csv == 1) {
  write.csv(csv_print_mod,"C:/Users/scastiello/OneDrive - Nexus365/1.- Contingency Learning (OBE-TM)/ABC paper/stats/exp4_github_LMM.csv",
            row.names = F)
  write.csv(to_print_posthoc,"C:/Users/scastiello/OneDrive - Nexus365/1.- Contingency Learning (OBE-TM)/ABC paper/stats/exp4_github_posthoc.csv",
            row.names = F)
}


######## VISUALIZATION ######## 
# reorder factor data
lf_noB$frequency <- as.integer(as.character(factor(lf_noB$frequency, levels = c("4","16","64"))))
lf_noB$Adjusted <- lf_noB$adjusted
lf_noB$Adjusted <- as.factor(lf_noB$Adjusted)
levels(lf_noB$Adjusted) <- c("none", "with duration")

# 5000 and 4000 with sf = 4.  With 1.25 times larger (5000/4000)
sf <- 1 # scaling factor
plot_title_size <- 42#24
strip_text_size <- 44#22
axis_title_size <- 40#20
axis_text_size <- 36#18
legend_title_size <- 40#20
legend_text_size <- 36#18
pos <- position_dodge(limits_vec_fre[2]*0.1)


fig5sup <- ggplot(lf_noB, aes(x = frequency, y = rating, col = Adjusted, shape = Adjusted, group = Adjusted)) +
  theme_bw() +
  geom_hline(yintercept = 0, color = "gray50", size = sf*1.2) +
  geom_line(aes(group = interaction(subject, Adjusted)), alpha = 0.05, size = sf*0.8) +
  geom_point(alpha = 0.1, size = sf*1.2, aes(y = jitter(rating, 0, 0.25), shape = Adjusted, fill = Adjusted), position = pos) +
  #geom_jitter(width = 2, height = 0.2, alpha = 0.2, size = 1.2, col = "gray30") +
  geom_smooth(method = "lm", size = sf*2, se = F, aes(linetype = Adjusted)) + 
  stat_summary(fun.data = "mean_cl_boot", size = sf*2, position = pos, stroke = sf*2, aes(shape = Adjusted, fill = Adjusted)) + # aes(fill = Adjusted)
  labs(title = "Exp 4: Frequency and adjusted for duration") + #subtitle = "Same baseline for Adjused and Cells") +
  xlab("Frequency") + ylab("Rating") +
  scale_fill_manual(values = c("black", "white"), breaks = c("none", "with duration")) +
  scale_shape_manual(values = c(19, 24), breaks = c("none", "with duration")) +
  #scale_linetype_manual(values=c("solid", "dotted")) +
  scale_colour_manual(values=c("#440154FF","#218F8DFF")) +
  scale_x_continuous(limits = limits_vec_fre, breaks = break_vec_fre, minor_breaks = NULL) + 
  scale_y_continuous(limits = c(-10.5, 10.5), breaks = c(-10, 0, 10)) +
  facet_grid(. ~ cell) + 
  coord_fixed(ratio = limits_vec_fre[2]/12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = sf*strip_text_size),
        plot.title = element_text(size = sf*plot_title_size, hjust = 0,face = "bold"),
        plot.subtitle = element_text(size = sf*plot_title_size/2, hjust = 0, face = "italic", color = "black"),
        plot.margin = unit(c(0.2,1,0.2,0.2), "cm"),
        axis.title.y = element_text(size = sf*axis_title_size, face = "bold"),
        axis.text.y = element_text(size = sf*axis_text_size, colour="black"),
        axis.title.x = element_text(size = sf*axis_title_size, face = "bold"),
        axis.text.x = element_text(size = sf*axis_text_size, colour="black"),
        axis.line = element_line(colour = 'black', size = sf*1.6),
        axis.ticks = element_line(size = sf*1.6),
        axis.ticks.length = unit(sf*10, "pt"),
        legend.title = element_text(size = sf*legend_title_size), 
        legend.text = element_text(size = sf*legend_text_size),
        legend.position = "bottom",
        legend.margin = margin(0.2, 0.2, 1, 0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),
        panel.border = element_rect(color = "black"),
        panel.spacing = unit(sf*1, "lines")
  )
fig5sup


if (print_png == 1) {
  ggsave("C:/Users/scastiello/OneDrive - Nexus365/1.- Contingency Learning (OBE-TM)/ABC paper/figures and tables/Fig5sup_github.png",
         plot = fig5sup, height = 16, width = 22)
}



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # Experiment 5 # # # # # # # # # # # # # # # # # # # # ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

exp05_lf <- read.csv("exp05_lf.csv")
exp05_wf <- read.csv("exp05_wf.csv")

# general characteristics
nSubj <- nrow(exp05_wf); nSubj
sex <- table(exp05_wf$sex); sex
meanAge <- mean(2020-exp05_wf$birthYear,na.rm=T); meanAge
sdAge <- sd(2020-exp05_wf$birthYear,na.rm=T); sdAge
rangeAge <- range(2020-exp05_wf$birthYear,na.rm=T); rangeAge

######## STATISTICAL ANALYSIS ######## 
# reorder factor order, and change levels order, where baseline is factor #1
exp05_lf$freq <- factor(exp05_lf$frequency, levels = c("18", "6", "54")); exp05_lf$frequency <- NULL
exp05_lf$dur <- factor(exp05_lf$duration, levels = c("900", "150", "5400")); exp05_lf$duration <- NULL



### filter data and created plot limits vectors
# plot and analysis with 0.
lf_noB <- exp05_lf[!is.na(exp05_lf$rating),]

# graphs data
break_vec_dur <- c(6,18,54)
limits_vec_dur <- c(0,54+break_vec_dur[3]*0.05)
###

# for non convergence models
# if max(abs(relgrad)) is <0.001 then things might be ok
# relgrad <- with(mFreq@optinfo$derivs,solve(Hessian,gradient)); max(abs(relgrad))

# FULL BASELINE MODEL
mFull <- lmer(rating ~ (freq * adjusted * cell) + cond_order + (1 | subject),
              #control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = F),
              REML = F, data = lf_noB); anova(mFull)

# FULL NO-BASELINE MODEL
mFull_nb <- lmer(rating ~ (freq * adjusted * cell) + cond_order + (1 | subject),
                 #control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = F),
                 REML = F, data = lf_noB[lf_noB$condition != "baseline",]); anova(mFull_nb)

# REDUCED BASELINE MODEL
mFull_red <- lmer(rating ~ (freq * cell) + cond_order + (1 | subject),
                  #control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = F),
                  REML = F, data = lf_noB); anova(mFull_red)

# REDUCED NO-BASELINE MODEL
mFull_nb_red <- lmer(rating ~ (freq * cell) + cond_order + (1 | subject),
                     #control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = F),
                     REML = F, data = lf_noB[lf_noB$condition != "baseline",]); anova(mFull_nb_red)


# Full against reduced
anova(mFull,mFull_red, test="LRT"); exp((BIC(mFull) - BIC(mFull_red))/2); AIC(mFull) - AIC(mFull_red)
anova(mFull_nb,mFull_nb_red, test="LRT"); exp((BIC(mFull_nb) - BIC(mFull_nb_red))/2); AIC(mFull_nb) - AIC(mFull_nb_red)


## post hoc exploration ##
mA <- lmer(rating ~ (freq * adjusted) + cond_order + (1 | subject), 
           REML = F, data = lf_noB[lf_noB$cell == "A",]); anova(mA)
mA_red <- lmer(rating ~ (freq) + cond_order + (1 | subject), 
               REML = F, data = lf_noB[lf_noB$cell == "A",]); anova(mA_red)
anova(mA,mA_red, test="LRT"); exp((BIC(mA) - BIC(mA_red))/2); AIC(mA) - AIC(mA_red)

mB <- lmer(rating ~ (freq * adjusted) + cond_order + (1 | subject), 
           REML = F, data = lf_noB[lf_noB$cell == "B",]); anova(mB)
mB_red <- lmer(rating ~ (freq) + cond_order + (1 | subject), 
               REML = F, data = lf_noB[lf_noB$cell == "B",]); anova(mB_red)
anova(mB,mB_red, test="LRT"); exp((BIC(mB) - BIC(mB_red))/2); AIC(mB) - AIC(mB_red)

mC <- lmer(rating ~ (freq * adjusted) + cond_order + (1 | subject), 
           REML = F, data = lf_noB[lf_noB$cell == "C",]); anova(mC)
mC_red <- lmer(rating ~ (freq) + cond_order + (1 | subject), 
               REML = F, data = lf_noB[lf_noB$cell == "C",]); anova(mC_red)
anova(mC,mC_red, test="LRT"); exp((BIC(mC) - BIC(mC_red))/2); AIC(mC) - AIC(mC_red)

mD <- lmer(rating ~ (freq * adjusted) + cond_order + (1 | subject), 
           REML = F, data = lf_noB[lf_noB$cell == "D",]); anova(mD)
mD_red <- lmer(rating ~ (freq) + cond_order + (1 | subject), 
               REML = F, data = lf_noB[lf_noB$cell == "D",]); anova(mD_red)
anova(mD,mD_red, test="LRT"); exp((BIC(mD) - BIC(mD_red))/2); AIC(mD) - AIC(mD_red)


csv_print <- cbind("duplicated baseline",table_long(report(mFull))); colnames(csv_print) <- c("model",colnames(csv_print)[2:8])
csv_print_mod <- table_long(report(mFull))
csv_print_nb <- cbind("no baseline",table_long(report(mFull_nb))); colnames(csv_print_nb) <- c("model",colnames(csv_print_nb)[2:8])
to_print <- rbind(csv_print,csv_print_nb)
to_print_posthoc <- cbind(c(rep("A",4),rep("B",4),rep("C",4),rep("D",4)),
                          rbind(table_long(report(mA))[1:4,],
                                table_long(report(mB))[1:4,],
                                table_long(report(mC))[1:4,],
                                table_long(report(mD))[1:4,]))

if (write_csv == 1) {
  write.csv(csv_print_mod,"C:/Users/scastiello/OneDrive - Nexus365/1.- Contingency Learning (OBE-TM)/ABC paper/stats/exp5_github_LMM.csv",
            row.names = F)
  write.csv(to_print_posthoc,"C:/Users/scastiello/OneDrive - Nexus365/1.- Contingency Learning (OBE-TM)/ABC paper/stats/exp5_github_posthoc.csv",
            row.names = F)
}



######## VISUALIZATION ######## 
# reorder factor data
lf_noB$freq <- as.integer(as.character(factor(lf_noB$freq, levels = c("6","18","54"))))
lf_noB$dur <- as.integer(as.character(factor(lf_noB$dur, levels = c("150","900","5400"))))
lf_noB$Adjusted <- as.factor(lf_noB$adjusted)
levels(lf_noB$Adjusted) <- c("none", "with duration")

# x and 4000 with sf = 4.  With 1.6 times larger (40/25)
sf <- 1 # scaling factor. paper = 4
plot_title_size <- 42#24
strip_text_size <- 44#22
axis_title_size <- 40#20
axis_text_size <- 36#18
legend_title_size <- 40#20
legend_text_size <- 36#18
pos <- position_dodge(limits_vec_dur[2]*0.1)


fig5inf <- ggplot(lf_noB, aes(x = freq, y = rating, col = Adjusted, shape = Adjusted, fll = Adjusted, group = Adjusted)) +
  theme_bw() +
  geom_hline(yintercept = 0, color = "gray50", size = sf*1.2) +
  geom_line(aes(group = interaction(subject, Adjusted)), alpha = 0.05, size = sf*0.8) +
  geom_point(alpha = 0.1, size = sf*1.2, aes(y = jitter(rating, 0, 0.25), shape = Adjusted, fill = Adjusted), position = pos) +
  #geom_jitter(width = 2, height = 0.2, alpha = 0.2, size = 1.2, col = "gray30") +
  geom_smooth(method = "lm", size = sf*2, se = F, aes(linetype = Adjusted)) + 
  stat_summary(fun.data = "mean_cl_boot", size = sf*2, position = pos, stroke = sf*2, aes(shape = Adjusted, fill = Adjusted)) +
  labs(title = "Exp 5: Frequency and overadjusting for duration") + #subtitle = "Same baseline for Adjused and Cells") + 
  xlab("Frequency") + ylab("Rating") +
  scale_fill_manual(values = c("black", "white"), breaks = c("none", "with duration")) +
  scale_shape_manual(values = c(19, 24), breaks = c("none", "with duration")) +
  #scale_linetype_manual(values=c("solid", "dotted")) +
  scale_colour_manual(values=c("#440154FF","#218F8DFF")) +
  scale_x_continuous(limits = limits_vec_dur, breaks = break_vec_dur, minor_breaks = NULL) + 
  scale_y_continuous(limits = c(-10.5, 10.5), breaks = c(-10, 0, 10)) +
  facet_grid(. ~ cell) + 
  coord_fixed(ratio = limits_vec_dur[2]/12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = sf*strip_text_size),
        plot.title = element_text(size = sf*plot_title_size, hjust = 0,face = "bold"),
        plot.subtitle = element_text(size = sf*plot_title_size/2, hjust = 0, face = "italic", color = "black"),
        plot.margin = unit(c(0.2,1,0.2,0.2), "cm"),
        axis.title.y = element_text(size = sf*axis_title_size, face = "bold"),
        axis.text.y = element_text(size = sf*axis_text_size, colour="black"),
        axis.title.x = element_text(size = sf*axis_title_size, face = "bold"),
        axis.text.x = element_text(size = sf*axis_text_size, colour="black"),
        axis.line = element_line(colour = 'black', size = sf*1.6),
        axis.ticks = element_line(size = sf*1.6),
        axis.ticks.length = unit(sf*10, "pt"),
        legend.title = element_text(size = sf*legend_title_size), 
        legend.text = element_text(size = sf*legend_text_size),
        legend.position = "bottom",
        legend.margin = margin(0.2, 0.2, 1, 0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),
        panel.border = element_rect(color = "black"),
        panel.spacing = unit(sf*1, "lines")
  )
fig5inf # 6400 x 4000 pixels


if (print_png == 1) {
  ggsave("C:/Users/scastiello/OneDrive - Nexus365/1.- Contingency Learning (OBE-TM)/ABC paper/figures and tables/Fig5inf_github.png",
         plot = fig5inf, height = 16, width = 22)
}