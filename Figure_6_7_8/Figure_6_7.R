#### Packages and functions ----
library(tidyverse)
library(patchwork)
library(stringr)
library(grid)
library(data.table)
library(gridExtra)
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
  
  
  ##### PREPARE & PLOT DATA OF DYNAMICS PER CLUSTER ########
  
  ##### Generate a dataframe with number of cells and percentage per cluster
  # df_Erk <- cluster1 %>% group_by(Time_in_min, clu_E) %>%
  #   dplyr::summarise(n=n()) %>% ungroup() %>%
  #   select(clu_E, n) %>% distinct() %>%
  # mutate(Cluster=paste0(clu_E,": ",n, " cells (",round(n/ntot*100),"%)"))
  
  df_Erk <- cluster1 %>% select(cond_inh,Unique_Unique_Object,Time_in_min,ERK,Cluster=clu_E)

  df_Akt <- cluster1 %>% select(cond_inh,Unique_Unique_Object,Time_in_min,Akt,Cluster=clu_A)
  ##### Bind the new cluster names (with the info on number of cells) to the raw daat
  # cluster1 <- cluster1 %>% left_join(df_text, by="clu_E")
  
  ##### Calculate summary statistics and store in dataframe
  summary_ERK <- df_Erk %>% group_by(Time_in_min, Cluster) %>% dplyr::summarise(ERK_mean = mean(ERK), ERK_sd = sd(ERK), n=n())
  summary_Akt <- df_Akt %>% group_by(Time_in_min, Cluster) %>% dplyr::summarise(Akt_mean = mean(Akt), Akt_sd = sd(Akt), n=n())
  
    
  p_ERK <- ggplot(summary_ERK, aes(x=Time_in_min)) + 
    # geom_line(data=cluster1, aes(x=Time_in_min, y=ERK, group=Unique_Unique_Object), color="lightblue1", alpha=0.03) +
    
    geom_line(data = summary_ERK, aes(x=Time_in_min, y=ERK_mean), color="blue", size=1.7) + geom_ribbon(data = summary_ERK, aes(x=Time_in_min, ymin=ERK_mean-ERK_sd, ymax=ERK_mean+ERK_sd), fill="blue", alpha=0.3) +
    facet_wrap(~Cluster, ncol=4, labeller = label_both) +
    ylim(-0.1, 1.2) +
    # geom_text(data=df_text, aes(x = -Inf, y = Inf, label = text, hjust = -.1, vjust = 2)) +
    geom_rect(xmin=-20,xmax=-5,ymin=1.0,ymax=1.2, aes(fill=Cluster)) +
    theme_light(base_size = 16) +
    xlab("time post stimulation (min)") + 
    ylab("C/N ratio change")  +
    labs(title = "B:  Average ERK C/N ratio per cluster") +
    scale_fill_manual(values=newColors) + 
    # coord_cartesian(clip = 'off') +
    # theme(strip.background = element_rect(color="grey", fill=NA, size=1)) +
    theme(strip.text = element_text(color = "black", face = "bold")) + 
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
  
  
  p_cERK <- ggplot(data=df_contribution_subset, aes(x=factor(cond_inh, levels = rev(levels(factor(cond_inh, levels=order)))), fill=as.factor(clu_E))) + scale_fill_manual(values=coulz)+
    geom_bar(position = position_fill(reverse = TRUE), width = 0.8) +
    labs(x = "", y = "Percentage of cells") +
    coord_flip() + 
    theme_light(base_size = 16) + 
    theme(legend.position="right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
          # axis.line = element_line(color="black", size = 1), axis.ticks.length = unit(0.3, "cm"), 
          NULL) +
    scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1), labels = percent) +
    labs(fill="Cluster") +
    labs(title = "A:  Contributions of the ERK response clusters") +
    theme(legend.position = "right", plot.title.position = "plot", plot.title = element_text(size=24)) +
    
    NULL

  

  
  ##### SAVE THE MULTIPANEL PLOT ########
  
  png(file=paste0("Figure_6.png"), width = 4000, height = 3600, units = "px", res = 400)
  (p_cERK/p_ERK) +   
    # plot_layout(design = layout) +
    plot_layout(heights = c(1, 3))+
    NULL

  dev.off()
  
  
  ########## Figure 8 ############
  
  p_Akt <- ggplot(summary_Akt, aes(x=Time_in_min)) + 
    # geom_line(data=cluster1, aes(x=Time_in_min, y=Akt, group=Unique_Unique_Object), color="lightblue1", alpha=0.03) +
    
    geom_line(data = summary_Akt, aes(x=Time_in_min, y=Akt_mean), color="blue", size=1.7) + geom_ribbon(data = summary_Akt, aes(x=Time_in_min, ymin=Akt_mean-Akt_sd, ymax=Akt_mean+Akt_sd), fill="blue", alpha=0.3) +
    facet_wrap(~Cluster, ncol=3, labeller = label_both) +
    ylim(-0.1, 0.4) +
    # geom_text(data=df_text, aes(x = -Inf, y = Inf, label = text, hjust = -.1, vjust = 2)) +
    geom_rect(xmin=-20,xmax=-5,ymin=0.32,ymax=0.4, aes(fill=Cluster)) +
    theme_light(base_size = 16) +
    xlab("time post stimulation (min)") + 
    ylab("C/N ratio change")  +
    labs(title = "B:  Average Akt C/N ratio per cluster") +
    scale_fill_manual(values=newColors) + 
    # coord_cartesian(clip = 'off') +
    # theme(strip.background = element_rect(color="grey", fill=NA, size=1)) +
    theme(strip.text = element_text(color = "black", face = "bold")) + 
    # geom_rect(xmin=-20,xmax=-5,ymin=1.30,ymax=1.42, aes(fill=Cluster)) +
    theme(legend.position = "none", plot.title.position = "plot", plot.title = element_text(size=24)) +
    
    NULL
  
  
  p_cAkt <- ggplot(data=df_contribution_subset, aes(x=factor(cond_inh, levels = rev(levels(factor(cond_inh, levels=order)))), fill=as.factor(clu_A))) + scale_fill_manual(values=coulz)+
    geom_bar(position = position_fill(reverse = TRUE), width = 0.8) +
    labs(x = "", y = "Percentage of cells") +
    coord_flip() + 
    theme_light(base_size = 16) + 
    theme(legend.position="right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
          # axis.line = element_line(color="black", size = 1), axis.ticks.length = unit(0.3, "cm"), 
          NULL) +
    scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1), labels = percent) +
    labs(fill="Cluster") +
    labs(title = "A:  Contributions of the Akt response clusters") +
    theme(legend.position = "right", plot.title.position = "plot", plot.title = element_text(size=24)) +
    
    NULL
  
  
  
  ##################### Plot the count for eco-occurence of ERK and Akt clusters #######
  
  bubble_plot <-  ggplot(data=df_contribution_subset, aes(x=clu_E,y=clu_A))+geom_count()+facet_wrap(~cond_inh) +
    labs(x = "ERK cluster", y = "Akt cluster") +
    theme_light(base_size = 16) + 
    theme(legend.position="bottom", panel.grid.minor = element_blank(), panel.background = element_blank(), 
          # panel.grid.minor = element_blank(),
          # axis.line = element_line(color="black", size = 1), axis.ticks.length = unit(0.3, "cm"), 
          NULL) +
    theme(strip.text = element_text(color = "black", face = "bold")) + 
    labs(size="Count") +
    labs(title = "C: Co-occurence of ERK and Akt clusters") +
    theme(legend.position = "right", plot.title.position = "plot", plot.title = element_text(size=24)) +
    
    NULL
  
  png(file=paste0("Bubble.png"), width = 3200, height = 2000, units = "px", res = 500)
  bubble_plot
  
  dev.off()
  
  
  
  ##### SAVE THE MULTIPANEL PLOT ########
  
  png(file=paste0("Figure_7.png"), width = 4000, height = 4400, units = "px", res = 400)
  (p_cAkt/p_Akt/bubble_plot) +   
    # plot_layout(design = layout) +
    plot_layout(heights = c(1, 2, 1.5))+
    NULL
  
  dev.off()
  
  
  

  
  
