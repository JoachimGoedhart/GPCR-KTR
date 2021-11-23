library(tidyverse)
library(scales)
library(patchwork)
library(data.table)

heatmap_plot <- function(df, x_var, y_var, group_var, Compound) {
  
  x_var <- enquo(x_var)
  y_var <- enquo(y_var)
  group_var <- enquo(group_var)
  
  plot <- ggplot(df, aes(!!x_var,!!y_var,fill=!!group_var)) + geom_tile() +
  xlab("Time post stimulation (min)") + 
    ylab("") + 
    scale_x_continuous(breaks=seq(0,60, by=15), labels=seq(0,60, by=15), limits=c(-8,60))  +
    scale_fill_gradientn(colours=c("white", "pink", "red", "violet", "blue", "blue4"), values = rescale(c(-0.3, 0.1, 0.4, 0.7, 1.3, 1.8)), na.value="white", guide="colorbar", limits = range(-0.3,1.82)) +  
    labs(fill = "C/N ERK") + labs(title=Compound) +
    theme(text = element_text(size=24),
          plot.title = element_text(hjust = 0.5),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          legend.key.size = unit(1.5, "cm"),
          axis.ticks.length=unit(0.25, "cm"),
          panel.background = element_blank(),
          axis.ticks.y = element_blank())
  
  return(plot)
}



### Define path and process files ----
dir <- "Data/"


## Experiments or unstimulated control
filelist <- list.files(path=dir, pattern = "DMSO.csv", recursive = FALSE)
# files <- c(files,"no_ligand_interpolated.csv")

df_input_list <- map(paste0(dir,filelist), fread)
names(df_input_list) <- gsub(filelist, pattern="\\..*", replacement="")
df_all <- bind_rows(df_input_list, .id = "id")

#Subtraction to set the baseline activity to zero
df_sub <- df_all %>% group_by(Unique_Object,Ligand,Condition) %>% arrange(Time_in_min) %>% mutate(Aktn=CN_AktRB-mean(CN_AktRB[5:6]), ERKn=CN_ERK-mean(CN_ERK[5:6])) %>% ungroup()


df_sub <- df_sub %>% dplyr::select(id, Time_in_min, Inhibitor, Ligand, Condition, Unique_Object, Aktn, ERKn) 

#Sort the data of single cell traces, based on the maximal value
df_sub$Unique_Object <- reorder(df_sub$Unique_Object, df_sub$ERKn, max, na.rm = TRUE)

#Set the time: define t=0 as the timepoint of addition
df_sub$Time_in_min <- df_sub$Time_in_min-23

################## UK ################

Compound <- "UK"

#Filter UK data, maximal concentration
df_new <- df_sub %>% filter(Ligand==Compound)
df_new <- df_new %>% filter(Condition==max(df_new$Condition))

#Save the data as PNG
p2 <- heatmap_plot(df_new, Time_in_min, Unique_Object, ERKn, Compound)


################## His ################

Compound <- "His"

#Filter UK data, maximal concentration
df_new <- df_sub %>% filter(Ligand==Compound) 
df_new <- df_new %>% filter(Condition==max(df_new$Condition))

#Save the data as PNG
p1 <- heatmap_plot(df_new, Time_in_min, Unique_Object, ERKn, Compound)

################## S1P ################

Compound <- "S1P"

#Filter UK data, maximal concentration
df_new <- df_sub %>% filter(Ligand==Compound)
df_new <- df_new %>% filter(Condition==max(df_new$Condition))

#Save the data as PNG
p3 <- heatmap_plot(df_new, Time_in_min, Unique_Object, ERKn, Compound)

################## none ################

Compound <- "none"

#Filter UK data, maximal concentration
df_new <- df_sub %>% filter(Ligand==Compound)
df_new <- df_new %>% filter(Condition==max(df_new$Condition))

#Save the data as PNG
p4 <- heatmap_plot(df_new, Time_in_min, Unique_Object, ERKn, Compound)

###### Save ######

png(file=paste0("Figure_2.png"), width = 1600, height = 600, units = "px")
  (p1|p2|p3|p4)+plot_layout(guides = 'collect')
dev.off()
