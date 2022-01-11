#### Packages and functions ----
require(dplyr)
require(ggplot2)
require(gridExtra)
require(grid)
library(data.table)

get_CI_half_width <- function(x) {
  n <- length(x)
  z_t <- qt(1 - (1 - 0.95) / 2, df = n - 1)
  z_t * sd(x) / sqrt(n) }

lower <- function(x) {mean(x) - get_CI_half_width(x)}
upper <- function(x) {mean(x) + get_CI_half_width(x)}

sem <- function(x){
  n <- length(x)
  sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)]))
}
upper_sem <- function(x){mean(x) + sem(x)}
lower_sem <- function(x){mean(x) - sem(x)}
upper_sd <- function(x){mean(x) + sd(x)}
lower_sd <- function(x){mean(x) - sd(x)}

#### Define variables ----
### norm: Normalize data?
### inhi: Read big files and split by inhibitor?
### conc: Read big files and split by concentration?
### none: Read big files and split by experiment?
norms = c(0, 1)
inhis = c(0, 1)
concs = c(0, 1)
nor = 1
inhi = 1
conc = 0
fs=1.5

## Single XY axis for all plots? Position?
sharedaxis = "yes"

x_pos = 0.20
y_pos = 0.25


## With labels (A,B,C...) per plot
withlabs = "yes"
dir <- "Data"

#### Code ----
      ligands = c("S1P", "UK", "His")

        for (lig in ligands){

          # gc()
          ## Find data files ----
          all <- list.files(path = dir , pattern = paste0("ALL_", lig, sep = ""), recursive = FALSE, full.names=TRUE)
          

          ### Define path and process files ----


          ## Loop through files ----
          for (i in all){
            #### Define ligand information ----
            if (grepl("His", i)){ligand = "uM His"}
            if (grepl("S1P", i)){ligand = "nM S1P"}
            if (grepl("UK", i)){ligand = "pM UK"}
            title = substr(i, 10, nchar(i)-4)
            exp = fread((i))
            # if (inhi < 1 && conc < 1){title = substr(i, 10, nchar(i)-4)} else {title = substr(i, 3, nchar(i)-4)}
            # if(any(grepl('Object', colnames(exp)))){object_type = "Object"}
            # if(any(grepl('Unique_Object', colnames(exp)))){object_type = "Unique_Object"}
            if(any(grepl('Unique_Unique_Object', colnames(exp)))){object_type = "Unique_Unique_Object"}
            
            exp[ ,c("Position", "meanArea", "Date", "Slide", "Ligand", "X", "Y")] <- list(NULL)
            exp$Time_in_min=exp$Time_in_min-23
            
              #  Normalization related
              i3 = -0.05
              exp <- exp %>% arrange(Time_in_min) %>% group_by(get(object_type)) %>% mutate(CN_AktRBc_norm=CN_AktRB-mean(CN_AktRB[5:6]), CN_ERKc_norm=CN_ERK-mean(CN_ERK[5:6]))
              ### The following line should be executed to use drift corrected data:
              # exp <- exp %>% arrange(Time_in_min) %>% group_by(get(object_type)) %>% mutate(CN_AktRBc_norm=CN_AktRBc-mean(CN_AktRBc[5:6]), CN_ERKc_norm=CN_ERKc-mean(CN_ERKc[5:6]))
              exp$CN_AktRBc = NULL
              exp$CN_ERKc = NULL
              exp <- exp %>% mutate(CN_AktRBc = CN_AktRBc_norm, CN_ERKc = CN_ERKc_norm)
              i4 = "_norm"

  
            #### Define basic plot layout ----
            plt <- ggplot(exp, aes(x=Time_in_min)) + theme_light(base_size = 16) + xlab("time (min)") + theme(text = element_text(size=20), axis.text.x = element_text(size = rel(fs*1.2)),
                        axis.text.y = element_text(size = rel(fs*1.2)), axis.title.y = element_text(size = rel(fs)), axis.title.x = element_text(size = rel(fs)), 
                       legend.text = element_text(size = rel(fs*0.8)), legend.title = element_text(size=rel(fs)), plot.margin=unit(c(1.5,1,1.5,1),"cm"), legend.spacing.x = unit(0.5, 'cm'),
                       legend.spacing.y = unit(0.5, 'cm'), legend.key.height = unit(1, "cm"), plot.title = element_text(hjust = 0.5, size=rel(fs*1.2)), plot.tag = element_text(size = rel(fs*1.2)))
            # gc()
            
            #### Split data by inhibitor/condition/experiment ----
            if (inhi > 0){exps <- split(exp, f = exp$Inhibitor)} else if (conc > 0) {exps <- split(exp, f = exp$Condition)} else {exps <- split(exp, f = exp$Experiment)}
            expss <- lapply(seq_along(exps), function(x) as.data.frame(exps[[x]])) 
            nexp <- length(exps)
            
            for (con in seq(1:nexp)){
              data = as.data.frame(expss[con])
              x <- (paste0("Set",con))
              eval(call("<-", as.name(x), data))
            }
            allsets=ls(pattern="Set")
            listlabs=c("A", "B", "C", "D", "E", "F", "G", "H")
            
            #### Generate individual plots ----
            for (i2 in allsets){
              ## By inhibitor ----
                exp <- get(i2) %>% group_by(Condition) %>% arrange(Condition) %>% mutate(n=length(unique(get(object_type))))
                exp <- exp %>% mutate(Condition2 = paste(Condition, "   (", n, ")", sep=""))
                exp$Condition2 <- as.factor(exp$Condition2)
                exp$Condition2 <- factor(exp$Condition2, levels = unique(exp$Condition2))
                data_summary <- exp %>% group_by(Time_in_min, Condition2, Inhibitor) %>% select(CN_AktRBc, CN_ERKc, Time_in_min, Condition2, Inhibitor) %>% summarise_all(funs(mean, upper, lower, upper_sem, lower_sem, upper_sd, lower_sd))
                pllt <- plt + ggtitle(as.character(unique(get(i2)$Inhibitor))) + theme(plot.title = element_text(hjust = 0.5))
                if(unique(get(i2)$Inhibitor)=="ymptx"){pllt <- plt + ggtitle("YM + PTx")}
                alfa = 0.02
 
              
              ## Generate ERK and Akt plots ----
              pltx <- pllt   + ylab("CN ERK") + geom_line(data=get(i2), aes(x=Time_in_min, y=CN_ERKc, group=get(object_type)), color="grey", alpha=alfa) +
                geom_line(data = data_summary, aes(x=Time_in_min, y=CN_ERKc_mean, group=Condition2, colour = Condition2), size=3, alpha=0.8) +
                geom_ribbon(data = data_summary, aes(x=Time_in_min, ymin=CN_ERKc_lower_sd, ymax=CN_ERKc_upper_sd, colour = Condition2), alpha=0) +
                coord_cartesian(ylim = c(i3, i3+0.85)) + labs(colour=ligand)
                # geom_segment(aes(x=22.8, xend=86.5, y=i3+0.85, yend=i3+0.85), size=3) + annotate("text", x=30, y=(i3+0.3+0.025), label= "ligand", size=8)
                # annotate("text", x=54.7, y=(i3+0.85+0.025), label= "ligand addition", size=8)

              plty <- pllt   + ylab("CN Akt") + geom_line(data=get(i2), aes(x=Time_in_min, y=CN_AktRBc, group=get(object_type)), color="grey", alpha=alfa) +
                geom_line(data = data_summary, aes(x=Time_in_min, y=CN_AktRBc_mean, group=Condition2, colour = Condition2), size=3, alpha=0.8) +
                geom_ribbon(data = data_summary, aes(x=Time_in_min, ymin=CN_AktRBc_lower_sd, ymax=CN_AktRBc_upper_sd, colour = Condition2), alpha=0) +
                coord_cartesian(ylim = c(i3, i3+0.3)) + labs(colour=ligand)
              # geom_segment(aes(x=22.8, xend=86.5, y=i3+0.3, yend=i3+0.3), size=3) + annotate("text", x=30, y=(i3+0.3+0.025), label= "ligand", size=8)
              # annotate("text", x=54.7, y=(i3+0.3+0.025), label= "ligand addition", size=8)
              
                ## For shared X/Y ----
              if (sharedaxis=="yes"){
                pltx <- pltx + theme(plot.margin=unit(c(1,0.5,0.3,1),"cm"), axis.title = element_text(size=0))
                plty <- plty + theme(plot.margin=unit(c(1,0.5,0.3,1),"cm"), axis.title = element_text(size=0))}
              
              if(withlabs=="yes"){
                pltx <- pltx + theme(plot.margin=unit(c(0.5,0.5,0,0.5),"cm"), axis.title = element_text(size=0)) + labs(tag = listlabs[which(i2==allsets)])
                plty <- plty + theme(plot.margin=unit(c(0.5,0.5,0,0.5),"cm"), axis.title = element_text(size=0)) + labs(tag = listlabs[which(i2==allsets)])}
            
              ## Assign plots to objects ----
              x <- (paste0("plte",which(i2==allsets)))
              eval(call("<-", as.name(x), pltx))
              
              x <- paste0("plta",which(i2==allsets))
              eval(call("<-", as.name(x), plty))
              
              # ## Plot ERK and Akt plots together per inhibitor ----
              # if (inhi > 0){
              #   pl1 <- grid.arrange(pltx, plty, nrow=2)
              #   pl1
              #   png(file=paste(title, "_", unique(exp$Inhibitor), "_sd", i4, ".png", sep=""), width = 1200, height = 1500, units = "px")
              #   plot(pl1)
              #   dev.off()
              # }
            }
            #### Generate grouped plots per kinase ----
            lista <- mget(paste0("plta", seq(1:nexp)))
            liste <- mget(paste0("plte", seq(1:nexp)))
            if (inhi+conc > 0) {colfac = 2; rowfac = 1} else {colfac=1; rowfac=2}
            
              #### With shared X/Y axis ----
            if (sharedaxis=="yes"){
              # top <- textGrob("TITLE", gp = gpar(fontsize = 35))
              bottom <- textGrob("time post stimulation (min)", gp = gpar(fontsize = 35), x=x_pos)
              
              png(bg = "white", file=paste(title, "_per_inh_sd", "_Akt.png", sep=""), width = 1000*nexp/colfac, height = 1100, units = "px")
              grid.arrange(grobs=lista, ncol = nexp/colfac, bottom = bottom, left = textGrob("Akt CN ratio change", rot = 90, gp = gpar(fontsize = 35), y=y_pos))
              dev.off()
            
              png(bg = "white", file=paste(title, "_per_inh_sd", "_ERK.png", sep=""), width = 1000*nexp/colfac, height = 1100, units = "px")
              grid.arrange(grobs=liste, ncol = nexp/colfac, bottom = bottom, left = textGrob("ERK CN ratio change", rot = 90, gp = gpar(fontsize = 35), y=y_pos))
              dev.off()
            }
            
#               #### Without shared X/Y axis ----
# p1 <- grid.arrange(grobs = lista, ncol = nexp/colfac, as.table = TRUE)
# p2 <- grid.arrange(grobs = liste, ncol = nexp/colfac, as.table = TRUE)
# pl1 <- grid.arrange(p2, p1, ncol=colfac, nrow=rowfac, top=textGrob(title, gp=gpar(fontsize=45,font=8)))
# 
# if (inhi > 0) {png(file=paste(title, "_per_inh_sd", i4, "_Akt.png", sep=""), width = 1000*nexp/colfac, height = 1500, units = "px")
# } else if (conc > 0) {png(file=paste(title, "_per_con_sd", i4, "_Akt.png", sep=""), width = 1000*nexp/colfac, height = 1500, units = "px")
# } else {png(file=paste(title, "_per_exp_sd", i4, "_Akt.png", sep=""), width = 1000*nexp/colfac, height = 750, units = "px")}
# plot(p1)
# dev.off()
# 
# if (inhi > 0) {png(file=paste(title, "_per_inh_sd", i4, "_ERK.png", sep=""), width = 1000*nexp/colfac, height = 1500, units = "px")
# } else if (conc > 0) {png(file=paste(title, "_per_con_sd", i4, "_ERK.png", sep=""), width = 1000*nexp/colfac, height = 1500, units = "px")
# } else {png(file=paste(title, "_per_exp_sd", i4, "_ERK.png", sep=""), width = 1000*nexp/colfac, height = 750, units = "px")}
# plot(p2)
# dev.off()
# 
# if (inhi > 0) {png(file=paste(title, "_per_inh_sd", i4, ".png", sep=""), width = 1000*nexp, height = 1500, units = "px")
# } else if (conc > 0) {png(file=paste(title, "_per_con_sd", i4, ".png", sep=""), width = 1000*nexp, height = 1500, units = "px")
# } else {png(file=paste(title, "_per_exp_sd", i4, ".png", sep=""), width = 1000*nexp, height = 1500, units = "px")}
# plot(pl1)
# dev.off()

            #### Remove objects ----
            # rm(list = ls(pattern="pl"))
            # rm(list = ls(pattern="xp"))
            # rm(list = ls(pattern="Set"))
            # gc()
          }
        }
#       }
#     }
#   }
# }


