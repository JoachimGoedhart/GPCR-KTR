#### Packages and functions ----
require(tidyverse)
require(gridExtra)
require(grid)
library(data.table)
library(patchwork)
library(stringr)

Confidence_level = 0.95

## Start ratio limits
limitsAkt <- c(0.25, 0.40, 0.55, 0.80)
limitsErk <- c(0.15, 0.35, 0.60, 0.85)


################# Read the data

dir <- "Data/"

filelist <- list.files(path=dir, pattern = ".csv", recursive = FALSE)
# files <- c(files,"no_ligand_interpolated.csv")

df_input_list <- map(paste0(dir,filelist), fread)
names(df_input_list) <- gsub(filelist, pattern=".csv", replacement="")
alln <- bind_rows(df_input_list, .id = "id") %>% dplyr::rename(filename=id)

# ligands = c("S1P", "UK", "His")
alln$Time_in_min=alln$Time_in_min-23

#Generate a unique id for each trace from each cell
alln <- alln %>% unite("Unique_id",c(Ligand, Condition, Unique_Object),remove = FALSE)

#Use time point 6&7 to determine the average start ratio for each trace.
alln_baseline <- alln %>% arrange(Ligand, Condition, Time_in_min) %>% group_by(Unique_id) %>% slice(6:7) %>% ungroup()
alln_baseline_avg <- alln_baseline %>% group_by(Unique_id) %>% dplyr::summarise(start_ERK=mean(CN_ERK), start_Akt=mean(CN_AktRB))

#Plot a histogram with start values
p_start_ERK <- ggplot(data=alln_baseline_avg,aes(x=start_ERK))+
  geom_histogram(binwidth = 0.05, fill="grey80", color="black")+xlim(0.1,1) +

  theme_light(base_size = 16) + xlab("ERK C/N start ratio") + ylab("Frequency") +
  theme(legend.position = "none", plot.title.position = "plot", plot.title = element_text(size=24)) +
  labs(title="A") +
  NULL

#Plot a histogram with start values
p_start_Akt <- ggplot(data=alln_baseline_avg,aes(x=start_Akt))+
  geom_histogram(binwidth = 0.05, fill="grey80", color="black")+xlim(0.1,1) +
  
  theme_light(base_size = 16) + xlab("Akt C/N start ratio") + ylab("Frequency") +
  theme(legend.position = "none", plot.title.position = "plot", plot.title = element_text(size=24)) +
  
  NULL

p_ERK_vs_Akt <- ggplot(data=alln_baseline_avg,aes(x=start_ERK, y=start_Akt))+
  geom_point(fill="black", alpha=0.05, size=2, shape=16)+xlim(0.1,.9) +ylim(0.1,.9) +
  
  theme_light(base_size = 16) + xlab("ERK C/N start ratio") + ylab("Akt C/N start ratio") +
  theme(legend.position = "none", plot.title.position = "plot", plot.title = element_text(size=24)) +
  labs(title="B") +
  NULL

#Group the start values in three bins
alln_baseline_avg <- alln_baseline_avg %>%
  mutate(bin_ERK = case_when(start_ERK < 0.35 ~ "<0.35",
                             start_ERK >=0.35 & start_ERK <0.6  ~ "0.35-0.60",
                             start_ERK >=0.6 ~"0.60-1.00"),
         bin_Akt = case_when(start_Akt < 0.4 ~ "<0.40", 
                             start_Akt >=0.4 & start_Akt <0.55  ~ "0.40-0.55",
                             start_Akt >=0.55 ~"0.55-1.00")
         )

alln <- alln %>% left_join(alln_baseline_avg, by="Unique_id")
alln <- alln %>% mutate(ERK_norm=CN_ERK-start_ERK, Akt_norm = CN_AktRB-start_Akt)

alln_summary_ERK <- alln %>% group_by(Time_in_min, bin_ERK, Ligand, Condition) %>%
  dplyr::summarise(n=n(),mean_ERK = mean(ERK_norm), sd_ERK=sd(ERK_norm)) %>%
  mutate(sem=sd_ERK/sqrt(n),
         CI = qt((1-Confidence_level)/2, n - 1) * sem)

alln_summary_Akt <- alln %>% group_by(Time_in_min, bin_Akt, Ligand, Condition) %>%
  dplyr::summarise(n=n(),mean_Akt = mean(Akt_norm), sd_Akt=sd(Akt_norm)) %>%
  mutate(sem=sd_Akt/sqrt(n),
         CI = qt((1-Confidence_level)/2, n - 1) * sem)

#Generate an index for each ligand concentration series to make sure that the colors are re-used for the different ligands
alln_summary_ERK <- alln_summary_ERK %>% group_by(Ligand,Time_in_min, bin_ERK) %>% dplyr::mutate(ligand_conc = row_number(Condition))%>%ungroup()

#Generate an index for each ligand concentration series to make sure that the colors are re-used for the different ligands
alln_summary_Akt <- alln_summary_Akt %>% group_by(Ligand,Time_in_min, bin_Akt) %>% dplyr::mutate(ligand_conc = row_number(Condition))%>%ungroup()

p_ERK <- ggplot(data=alln_summary_ERK %>% 
                  ### This is used to filter responses for the highest ligand concentration
                  group_by(Ligand) %>% 
                  filter(Condition==max(Condition)),
                aes(x=Time_in_min))+
  facet_grid(Ligand~bin_ERK)+
  geom_line(aes(x=Time_in_min, y=mean_ERK), color="blue", size=1)+
  geom_ribbon(aes(x=Time_in_min, ymin=mean_ERK-sd_ERK, ymax=mean_ERK+sd_ERK, color=as.factor(ligand_conc)), fill="blue", alpha=0.2, size=0.1) +
  theme_light(base_size = 16) + xlab("time post stimulation (min)") + ylab("ERK C/N ratio change") +
  theme(legend.position = "none", plot.title.position = "plot", plot.title = element_text(size=24)) +
  labs(title="C") +

  NULL

p_Akt <- ggplot(data=alln, aes(x=Time_in_min))+
  # geom_line(data=alln, aes(y=ERK_norm, group=Unique_id),alpha=0.02,color="black")+
  facet_grid(Ligand~bin_Akt)+
  geom_line(data=alln_summary_Akt, aes(x=Time_in_min, y=mean_Akt, color=as.factor(ligand_conc)))+
  geom_ribbon(data=alln_summary_Akt, aes(x=Time_in_min, ymin=mean_Akt-sd_Akt, ymax=mean_Akt+sd_Akt, fill=as.factor(ligand_conc), color=as.factor(ligand_conc)), alpha=0.05, size=0.1) +
  theme_light(base_size = 16) + xlab("time post stimulation (min)") + ylab("Akt C/N ratio change") + 
  theme(legend.position = "none", plot.title.position = "plot", plot.title = element_text(size=24)) +
  
  NULL


layout <- "
AB
CC
CC
"

png(file=paste0("Figure_6.png"), width = 4000, height = 4800, units = "px", res = 400)
p_start_ERK + p_ERK_vs_Akt + p_ERK +   plot_layout(design = layout)
  # plot_layout(guides = 'collect')+ 

dev.off()


pdf(file=paste0("Figure_6.pdf"), width = 4000/12, height = 4800/12)
p_start_ERK + p_ERK_vs_Akt + p_ERK +   plot_layout(design = layout)
# plot_layout(guides = 'collect')+ 

dev.off()


