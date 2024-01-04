#######################################################################
#######  Making heatmaps for Paper LE8 and ADRDs mediated by proteomics
#######  Mission: Plot the betas and mediation effects  
#######  Programmer: Yi-Han Hu
#######  Date: Dec. 27 2023
#######################################################################

op <- options(nwarnings = 10000)
# --------------------------------------
# Specify working directory where the script and data files are
# --------------------------------------
WorkingDirectory = "C:/Users/huy13/Box/Projects/May/Paper_LE8_ADRDs_proteomics/"

# --------------------------------------
# Set working directory
# --------------------------------------
setwd(WorkingDirectory)

# --------------------------------------
# simple function
# --------------------------------------
floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)

# --------------------------------------
# Turn off scientific notation
# --------------------------------------
options(scipen=999)

# --------------------------------------
# Install/load the packages
# --------------------------------------
library(haven)
library(tidyr) 
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape)
library(purrr)
library(data.table)
library(sjmisc)
library(RColorBrewer)
library(ggnewscale)
library(scico)
# ---------------------------------------------------------------------------------------------#
# ---------------------------------------------------------------------------------------------#
# ---------------------------------- Part 1 Data preprocess -----------------------------------#
# ---------------------------------------------------------------------------------------------#
# ---------------------------------------------------------------------------------------------#
# --------------------------------------
# Load data:
# --------------------------------------

## 1.	No pure mediation: non-significant PIE
NO_MEDIATION_GROUP <- read_dta("Data/NO_MEDIATION_GROUPA_wide.dta")
dim(NO_MEDIATION_GROUP)
head(NO_MEDIATION_GROUP)
colnames(NO_MEDIATION_GROUP)

NO_MEDIATION_GROUP <- NO_MEDIATION_GROUP %>% 
  dplyr::select(protein, betate, pte, betacde, pcde, betaintref, pintref, betaintmed, pintmed, betapie, ppie) %>% 
  dplyr::rename(Protein = protein, TE = betate, p_TE = pte, CDE = betacde, p_CDE = pcde,
                INTREF = betaintref, p_INT = pintref, INTMED = betaintmed, p_INTMED = pintmed, PIE = betapie, p_PIE = ppie) %>% 
  mutate(p_mediated = ifelse(!is.na(p_TE) & p_TE < 0.05 & !is.na(p_PIE) & p_PIE < 0.05, 0.01, 0.1)) %>% 
  mutate_all(~ ifelse(is.na(.) | . == "#DIV/0!" | . == "", NA, .))

NO_MEDIATION_GROUP.long <- to_long(NO_MEDIATION_GROUP, keys = 'term',
                                  values = c('estimate','p'), 
                                  c('TE','CDE','INTREF', 'INTMED', 'PIE'),
                                  c('p_TE','p_CDE','p_INT', 'p_INTMED', 'p_PIE'))

## 2.	Inconsistent mediation: PIE is significant but CDE>TE because PIE goes in the opposite direction to TE.
INCONSISTENT_MEDIATION_GROUP <- read_dta("Data/INCONSISTENT_MEDIATION_GROUPB_wide.dta")
dim(INCONSISTENT_MEDIATION_GROUP)
head(INCONSISTENT_MEDIATION_GROUP)
colnames(INCONSISTENT_MEDIATION_GROUP)


INCONSISTENT_MEDIATION_GROUP <- INCONSISTENT_MEDIATION_GROUP %>% 
  dplyr::select(protein, betate, pte, betacde, pcde, betaintref, pintref, betaintmed, pintmed, betapie, ppie) %>% 
  dplyr::rename(Protein = protein, TE = betate, p_TE = pte, CDE = betacde, p_CDE = pcde,
                INTREF = betaintref, p_INT = pintref, INTMED = betaintmed, p_INTMED = pintmed, PIE = betapie, p_PIE = ppie) %>% 
  mutate(p_mediated = ifelse(!is.na(p_TE) & p_TE < 0.05 & !is.na(p_PIE) & p_PIE < 0.05, 0.01, 0.1)) %>% 
  mutate_all(~ ifelse(is.na(.) | . == "#DIV/0!" | . == "", NA, .))

INCONSISTENT_MEDIATION_GROUP.long <- to_long(INCONSISTENT_MEDIATION_GROUP, keys = 'term',
                                  values = c('estimate','p'), 
                                  c('TE','CDE','INTREF', 'INTMED', 'PIE'),
                                  c('p_TE','p_CDE','p_INT', 'p_INTMED', 'p_PIE'))

## 3.	Consistent mediation: PIE is significant and CDE<TE because PIE goes in the same direction as TE
CONSISTENT_MEDIATION_GROUP <- read_dta("Data/CONSISTENT_MEDIATION_GROUPC_wide.dta")
dim(CONSISTENT_MEDIATION_GROUP)
head(CONSISTENT_MEDIATION_GROUP)
colnames(CONSISTENT_MEDIATION_GROUP)

CONSISTENT_MEDIATION_GROUP <- CONSISTENT_MEDIATION_GROUP %>% 
  dplyr::select(protein, betate, pte, betacde, pcde, betaintref, pintref, betaintmed, pintmed, betapie, ppie, betapct_pie) %>% 
  dplyr::rename(Protein = protein, TE = betate, p_TE = pte, CDE = betacde, p_CDE = pcde,
                INTREF = betaintref, p_INT = pintref, INTMED = betaintmed, p_INTMED = pintmed, PIE = betapie, p_PIE = ppie, percent_mediated = betapct_pie) %>% 
  mutate(p_mediated = ifelse(!is.na(p_TE) & p_TE < 0.05 & !is.na(p_PIE) & p_PIE < 0.05, 0.01, 0.1)) %>% 
  mutate_all(~ ifelse(is.na(.) | . == "#DIV/0!" | . == "", NA, .))

CONSISTENT_MEDIATION_GROUP.long <- to_long(CONSISTENT_MEDIATION_GROUP, keys = 'term',
                                  values = c('estimate','p'), 
                                  c('TE','CDE','INTREF', 'INTMED', 'PIE','percent_mediated'),
                                  c('p_TE','p_CDE','p_INT', 'p_INTMED', 'p_PIE', 'p_mediated'))

# ---------------------------------------------------------------------------------------------#
# ---------------------------------------------------------------------------------------------#
# ------------------------------------- Part 2 Heat maps --------------------------------------#
# ---------------------------------------------------------------------------------------------#
# ---------------------------------------------------------------------------------------------#
# --------------------------------------
# Heatmap for all outcomes
# --------------------------------------
colnames(NO_MEDIATION_GROUP)
# combine all outcome

# List of datasets
datasets <- list(NO_MEDIATION_GROUP.long, INCONSISTENT_MEDIATION_GROUP.long, CONSISTENT_MEDIATION_GROUP.long)

# List of modality labels
modalities <- c("No mediation", "Inconsistent mediation", "Consistent mediation")

# Function to append modality column
append_modality <- function(data, modality) {
  data %>% mutate(modality = modality)
}

# Apply function to each dataset
datasets_new <- map2(datasets, modalities, append_modality)

# Combine datasets
all_mediation.long <- bind_rows(datasets_new)



heatmap <- function(data, outcome = "No mediation"){
  # color for mediation
  pal <- scico(7, palette = 'acton')
  
  # set limits based on the maximum absolute value across all exposure-outcome pairs to ensure consistency.
  # and only based on those significant and not mediation
  temp_data <- data %>% filter(term != "percent_mediated" & p < 0.05 & !is.na(estimate))
  max_abs_value <- max(abs(c(range(as.numeric(temp_data$estimate)))))
  
  data.long <- data %>% 
    mutate(ap=ifelse(p < 0.01, 1,
                     ifelse(p >= 0.01 & p < 0.05, 2, 3)),
           ap=ifelse(is.na(p), 3, ap),
           aq=ifelse(p < 0.05 , "Pass","insig"),
           aq.mediation=ifelse(p < 0.05 , "Pass.mediation", "insig"),
           bg.line=ifelse(term=="percent_mediated", "White", "Dark Grey"),
           bg.color=ifelse(term=="percent_mediated", "Dark Grey", "White"),
           break.mediation = ifelse(term != "percent_mediated", NA,
                                    ifelse(abs(estimate) <= 0.10, 1,
                                           ifelse(abs(estimate) > 0.10 & abs(estimate) <= 0.20, 2,
                                                  ifelse(abs(estimate) > 0.20 & abs(estimate) <= 0.30, 3,
                                                         ifelse(abs(estimate) > 0.30 & abs(estimate) <= 0.40, 4,
                                                                ifelse(abs(estimate) > 0.40 & abs(estimate) <= 0.50, 5, 6))))))) %>% 
    arrange(factor(Protein, levels = unique(data$Protein)), factor(term, levels = c('TE', 'CDE','INTREF', 'INTMED', 'PIE','percent_mediated'))) %>% 
    mutate(term = ifelse(term == "percent_mediated", "%mediated", term),
           break.mediation = factor(break.mediation, levels = c("1", "2", "3", "4", "5", "6", "NaN"))) %>% 
    filter(modality == outcome)
  
  p.plot <- ggplot(data = data.long, aes(x = forcats::fct_rev(factor(Protein, levels = unique(Protein))), y = factor(term, levels = c('TE', 'CDE','INTREF', 'INTMED', 'PIE','%mediated'))))+
    geom_tile(color = data.long$bg.line, fill = data.long$bg.color)+
    geom_point(data= subset(data.long, term %in% c('TE', 'CDE','INTREF', 'INTMED', 'PIE')),
               aes(shape=factor(aq),
                   size=factor(ap), 
                   fill=estimate), na.rm = FALSE)+
    scale_fill_scico(palette = "vik", midpoint = 0, 
                     limits = c(-max_abs_value, max_abs_value+0.05),
                     aesthetics = c("colour","fill")) +
    scale_size_manual(values=c(4, 2.5, 0.5), labels = c("< .01", "< .05", "\u2265 .05"))+
    new_scale_color() +
    geom_point(data = subset(data.long, term %in% c("%mediated")), size=6,
               aes(shape=factor(aq.mediation),
                   color=break.mediation), na.rm = FALSE)+
    scale_shape_manual(values=c('insig'=1, 'Pass.mediation'=16, 'Pass'=21), guide = "none")+
    scale_color_manual(breaks=c(1, 2, 3, 4, 5, 6), drop = FALSE, labels = c("\U2264 10%", "10% - 20%", "20% - 30%", "30% - 40%", "40% - 50%", "> 50%"), values=pal)+
    guides(size = guide_legend(override.aes = list(shape = c(21, 21, 1), fill = c("black", "black", "white")))) +
    labs(title=paste("Heatmap (Four-way decomposition models of selected proteins - ",  outcome, ")", sep = ""),
         # subtitle=paste("", sep = ""),
         x=paste("Proteins (", sum(!is.na(unique(data.long$Protein)))," out of 144 selected proteins)", sep = ""),
         y="",
         size=paste("p-value\nsolid circle: p<.05"), fill=(expression(paste(beta," coefficients"))), color=("% mediated"),
         caption="TE: Total effect; CDE: Controlled direct effect; INTREF: Interaction referent;\nINTMED: Mediated interaction; PIE: Pure indirect effect; \n% mediated is the percent of total effect that is pure indirect effect. No p-values were generated.") +
    theme(plot.title = element_text(color="Dark blue", size=13, face="bold.italic", hjust = 0.5),
          plot.subtitle=element_text(size=10, hjust=0.5, face="italic", color="Dark blue"),
          plot.caption=element_text(size=9, hjust=0.5, color="Dark grey"),
          axis.title.x = element_text(color="deepskyblue", size=11, face="bold"),
          axis.text.x = element_text(angle = 90, size = 9, vjust = 0.5),
          aspect.ratio=2/7)+
    coord_fixed()
  
  return(p.plot)
}


dir.create(paste(WorkingDirectory,"Output/plot",sep=""), recursive = TRUE)
plot.out.folder <- paste(WorkingDirectory,"Output/plot/",sep="")


outcome.list <- c("No mediation", "Inconsistent mediation", "Consistent mediation")
for (outcome in outcome.list){
  heatmap.hide.unsig <- heatmap(data = all_mediation.long, outcome = outcome)
  ggsave(paste(plot.out.folder, outcome, "_LE8_ADRDs_proteins_heatmap_hide_unsig.jpeg",sep=""), heatmap.hide.unsig, width = 12, height = 6, units = "in", dpi = 300)
} 