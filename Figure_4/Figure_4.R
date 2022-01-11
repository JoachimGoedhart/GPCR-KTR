## Packages ----
library(tidyverse)
library(plyr)
library(data.table)
library(patchwork)

### Define path and process files ----
dir <- "Data/"

## Experiments or unstimulated control
filelist <- list.files(path=dir, pattern = ".csv", recursive = FALSE)
# files <- c(files,"no_ligand_interpolated.csv")

df_input_list <- map(paste0(dir,filelist), fread)
names(df_input_list) <- gsub(filelist, pattern=" nac.csv", replacement="")
alln <- bind_rows(df_input_list, .id = "id") %>% dplyr::rename(Ligand=id)

## Ligand data ----
## Remove experiments excluded from clustering
alln = alln[!grepl("20200611 P1|20200701 P1|20200526 P1", alln$Experiment),]

alln$Condition=factor(alln$Condition)
# Remove no stimulation data
allnn = alln %>% filter(Concen!=0.00)
levels(allnn$Condition) <- c(levels(allnn$Condition), "YM+PTx")
allnn$Condition[allnn$Condition=="ymptx"] <- "YM+PTx"
allnn$Condition = droplevels(allnn$Condition)

# cors <- ddply(allnn, c("Ligand", "Condition"), summarize, intercept=round(coef(lm(ERKnac ~ Aktnac))[1],2), slope=round(coef(lm(ERKnac ~ Aktnac))[2],2), R=round(sqrt(as.numeric(summary(lm(ERKnac ~ Aktnac))[8])),2))

pltv1 = ggplot(allnn, aes(x=Aktnac, y = ERKnac)) + ylab("AUC ERK") + xlab("AUC Akt") + theme_light(base_size = 16) + labs(colour = "")  + facet_grid(Ligand ~ Condition, scales = "fixed") + 
  stat_density_2d(aes(fill = ..level..), geom = "polygon", n=100, contour_var = "ndensity", alpha=0.5, bins=20) + scale_fill_gradientn(name="Density", colours = rainbow(10), breaks=c(0,0.25,0.5,0.75,1), limits=c(0,1)) + 
  theme(plot.title.position = "plot", strip.text = element_text(size = 22, colour = "black"),  axis.text = element_text(size=18), axis.title = element_text(size=20), legend.key.size = unit(1.8,"line"), legend.title = element_text(size = 18), legend.text = element_text(size = 13)) +
  # geom_label(data=cors, aes(x = 0, y = 10, label=paste("slope =",slope,", R =",R)), size=5) + stat_smooth(method=lm, color="grey50", alpha=0.3) +
  xlim(-1.5,3) + ylim(-1,10) + labs(title="A") +
  NULL

png(file=paste0("AUC Akt vs ERK.png"), width = 1200, height = 1000, units = "px")

plot(pltv1)
dev.off()


## No stimulation data ----
nolig = alln %>% filter(Concen==0.00, Condition=="DMSO") %>% mutate(Ligand="none")

pltnl = ggplot(nolig, aes(x=Aktnac, y = ERKnac)) + ylab("AUC ERK") + xlab("AUC Akt") + theme_light(base_size = 16) + facet_grid(Ligand ~ Condition, scales = "fixed") +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", n=100, contour_var = "ndensity", alpha=0.5, bins=20) + scale_fill_gradientn(name="Density", colours = rainbow(10), breaks=c(0,0.25,0.5,0.75,1), limits=c(0,1))+
  scale_x_continuous(limits=c(layer_scales(pltv1)$x$range$range[1], layer_scales(pltv1)$x$range$range[2])) + scale_y_continuous(limits=c(layer_scales(pltv1)$y$range$range[1], layer_scales(pltv1)$y$range$range[2])) + 
  theme(plot.title.position = "plot", strip.text = element_text(size = 22, colour = "black"), axis.text = element_text(size=18), axis.title = element_text(size=20), legend.key.size = unit(1.8,"line"), legend.title = element_text(size = 18), legend.text = element_text(size = 13)) +
  # geom_label(data=cors, aes(x = 0, y = 10, label=paste("slope =",slope,", R =",R)), size=5) + stat_smooth(method=lm, color="grey50", alpha=0.3) +
  xlim(-1.5,3) + ylim(-1,10) + labs(title="B") +
  NULL

png(file=paste0("AUC no ligand.png"), width = 400, height = 350, units = "px")

plot(pltnl)
dev.off()

layout <- "
AAAA
AAAA
AAAA
B###
"

png(file=paste0("Figure_4.png"), width = 1000, height = 1000, units = "px")
((pltv1+ theme(legend.position="none"))/(pltnl+ theme(legend.position="none")))+
  # plot_layout(guides = 'collect')+ 
  plot_layout(design = layout)+ theme(legend.position= c(1.5, 0.8))
dev.off()







