#### Packages and functions ----
library(tidyverse)
library(patchwork)
library(stringr)
library(data.table)
library(scales)

#### Define type of plot and other variables ----

Okabe_Ito <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
coulz <- c("black", "orange3", "red2", "blue", "brown", "purple3", "green4", "yellow3")

newColors <- coulz


##### READ THE DATA AND RENUMBER THE CLUSTERS ########

  ## Read csv files ----
  cluster1 = fread("Data/ALL_ERK_Akt_cluster_series.csv", header = TRUE)
  cluster1$Time_in_min=cluster1$Time_in_min-23
  
  # Remove this column
  cluster1$Cluster <- NULL
  
  #Calculate number of clusters for ERK/Akt
  ncla = max(cluster1$clu_A)
  ncle = max(cluster1$clu_E)
  
  
  ## Rename clusters by given order - in this case the average response
  cluaa = cluster1 %>% group_by(clu_A) %>% dplyr::summarise_all(funs(mean)) %>% arrange(Akt)
  cluee = cluster1 %>% group_by(clu_E) %>% dplyr::summarise_all(funs(mean)) %>% arrange(ERK)
  
  mapa = setNames(c(cluaa$clu_A), c(1:ncla))
  cluster1$clu_A <- names(mapa)[match(cluster1$clu_A, mapa)]
  mape = setNames(c(cluee$clu_E), c(1:ncle))
  cluster1$clu_E <- names(mape)[match(cluster1$clu_E, mape)]

  ##### Determine the maximum number of cells
  nnn <- cluster1 %>% group_by(Time_in_min) %>% dplyr::summarise(n=n())
  ntot <- max(nnn$n)
  
  
  ##### Calculate summary statistics and store in dataframe
  df_traces <- cluster1 %>% select(Time_in_min,ERK,Cluster=clu_E) %>% group_by(Time_in_min, Cluster) %>% dplyr::summarise(ERK_mean = mean(ERK))
  p_traces <- ggplot(df_traces, aes(x=Time_in_min)) + 
    geom_line(aes(x=Time_in_min, y=ERK_mean), color="grey", size=2) +
    facet_wrap(~Cluster, ncol=8) +
    geom_text(data=df_traces%>%filter(Time_in_min==21.5),aes(label=Cluster, x=Time_in_min), y=1.1, color="black", size=6)+
    ylim(-0.1, 1.2) +
    # geom_text(data=df_text, aes(x = -Inf, y = Inf, label = text, hjust = -.1, vjust = 2)) +
    geom_rect(xmin=-20,xmax=-5,ymin=1.0,ymax=1.2, aes(fill=Cluster)) +
    theme_void(base_size = 16) +
    # labs(title = "A:  Average ERK C/N ratio per cluster") +
    scale_fill_manual(values=newColors) + 
    labs(x="",y="") +
    # coord_cartesian(clip = 'off') +
    # theme(strip.background = element_rect(color="grey", fill=NA, size=1)) +
    theme(strip.text = element_blank()) +
    # geom_rect(xmin=-20,xmax=-5,ymin=1.30,ymax=1.42, aes(fill=Cluster)) +
    theme(legend.position = "none", plot.title.position = "plot", plot.title = element_text(size=24)) +
    
    NULL
  
  
  ##### PREPARE & PLOT DATA OF CLUSTER CONTRIBUTIONS ########
  
  df_contribution <- cluster1 %>% select (cond_inh, clu_E, clu_A, Unique_Unique_Object) %>% distinct()
  
  ## Select a subset of the data and rename the conditions
  df_contribution_subset <- df_contribution[grepl("0.0 DMSO|100.0 His DMSO|100.0 UK DMSO|1300.0 S1P DMSO", df_contribution$cond_inh),]
  
  df_contribution_subset$cond_inh[df_contribution_subset$cond_inh=="0.0 DMSO"] <- "DMSO"
  # df_contribution_subset$cond_inh[df_contribution_subset$cond_inh=="0.0 YM"] <- "YM"
  # df_contribution_subset$cond_inh[df_contribution_subset$cond_inh=="0.0 PTx"] <- "PTx"
  df_contribution_subset$cond_inh[df_contribution_subset$cond_inh=="100.0 His DMSO"] <- "His"
  df_contribution_subset$cond_inh[df_contribution_subset$cond_inh=="100.0 UK DMSO"] <- "UK"
  df_contribution_subset$cond_inh[df_contribution_subset$cond_inh=="1300.0 S1P DMSO"] <- "S1P"
  
  order <- c("DMSO", "His", "UK", "S1P")
  
  df_contribution <- df_contribution  %>% 
    separate(cond_inh, c("Concentration","Ligand","Inhibitor"), sep=" ") %>%
    mutate(Concentration = round(as.numeric(Concentration), 1))
  
  df_contribution$Inhibitor[df_contribution$Inhibitor=="ymptx"] <- "YM+PTx"
  
  df_contribution_S1P <- df_contribution[grepl("S1P", df_contribution$Ligand),]
  
  
  
  p_His <- df_contribution %>% filter(Ligand=="His") %>% ggplot(aes(x=(as.factor(Concentration)), fill=as.factor(clu_E))) + scale_fill_manual(values=coulz)+
    geom_bar(position = position_fill(reverse = TRUE), width = 0.8) +
    labs(x = "[µM]", y = "Percentage of cells") +
    facet_wrap(~Inhibitor, ncol=2)+
    coord_flip() + 
    theme_light(base_size = 16) + 

    scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1), labels = percent) +
    labs(fill="Cluster") +
    labs(title = "A:  Histamine") +
    guides(fill = guide_legend(nrow = 1)) +
    theme(plot.title.position = "plot", plot.title = element_text(size=24)) +
    theme(strip.text = element_text(color = "black", face = "bold"), strip.background = element_rect(fill = NA, colour = NA)) + 
    scale_x_discrete(limits = rev) +
    NULL
  
  p_UK <- df_contribution %>% filter(Ligand=="UK") %>% ggplot(aes(x=(as.factor(Concentration)), fill=as.factor(clu_E))) + scale_fill_manual(values=coulz)+
    geom_bar(position = position_fill(reverse = TRUE), width = 0.8) +
    labs(x = "[pM]", y = "Percentage of cells") +
    facet_wrap(~Inhibitor, ncol=2)+
    coord_flip() + 
    theme_light(base_size = 16) + 
    scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1), labels = percent) +
    labs(fill="Cluster") +
    labs(title = "B:  UK") +
    guides(fill = guide_legend(nrow = 1)) +
    theme(plot.title.position = "plot", plot.title = element_text(size=24)) +
    theme(strip.text = element_text(color = "black", face = "bold"), strip.background = element_rect(fill = NA, colour = NA)) + 
    scale_x_discrete(limits = rev) +
    NULL
  
  
  
  p_S1P <- df_contribution %>% filter(Ligand=="S1P") %>% ggplot(aes(x=(as.factor(Concentration)), fill=as.factor(clu_E))) + scale_fill_manual(values=coulz)+
    geom_bar(position = position_fill(reverse = TRUE), width = 0.8) +
    labs(x = "[nM]", y = "Percentage of cells") +
    facet_wrap(~Inhibitor, ncol=2)+
    coord_flip() + 
    theme_light(base_size = 16) + 
    scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1), labels = percent) +
    labs(fill="Cluster") +
    labs(title = "C:  Sphingosine 1-phosphate") +
    guides(fill = guide_legend(nrow = 1)) +
    theme(plot.title.position = "plot", plot.title = element_text(size=24)) +
    theme(strip.text = element_text(color = "black", face = "bold"), strip.background = element_rect(fill = NA, colour = NA)) + 
    scale_x_discrete(limits = rev) +
    NULL

  
  

  
  ##### SAVE THE MULTIPANEL PLOT ########
  
  png(file=paste0("Figure_8.png"), width = 3200, height = 5000, units = "px", res = 400)
  (p_traces/p_His/p_UK/p_S1P)  +
    # plot_layout(design = layout) +
    plot_layout(heights = c(0.8,2, 0.8,2))+
    plot_layout(guides = 'collect') & theme(legend.position = 'none')

  
  dev.off()
  
  

  
