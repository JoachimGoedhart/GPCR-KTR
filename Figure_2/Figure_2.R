### Load packages and functions ----
library(dplyr)
library(tidyverse)
library(broom)
library(drc)
library(modelr)
library(scales)

myspread <- function(df, key, value) {
  keyq <- rlang::enquo(key)
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% gather(variable, value, !!!s) %>%
    unite(temp, !!keyq, variable) %>%
    spread(temp, value)
}

### Define path and process files ----
dir <- "Data/"

#Get a list with all csv files
filelist = list.files(path=dir, pattern="nac.csv")

#read all csv files and put in df_input_list
df_input_list <- lapply(paste0(dir,filelist), read.csv)

#get the filenames, remove extension and use as "id"
names(df_input_list) <- gsub(filelist, pattern="\\..*", replacement="")

#Merge all the dataframes and use the filenames as id
df <- bind_rows(df_input_list, .id = "id")

# Exclude experiments excluded from clustering
df = df[!grepl("20200611 P1|20200701 P1|20200526 P1", df$Experiment),]

#Set zero values to a 'low' value, to enable a log10 scale on x-axis
df <- df %>%
  mutate(Concen = ifelse((Concen == 0 & id == 'His nac'),
                         yes = 0.02,
                         no = Concen),
         Concen = ifelse((Concen == 0 & id == 'S1P nac'),
                         yes = 4,
                         no = Concen),
         Concen = ifelse((Concen == 0 & id == 'UK nac'),
                         yes = 0.04,
                         no = Concen)
  ) #%>% filter(Condition!='ymptx')

df$id[df$id=="His nac"] <- "Histamine [µM]"
df$id[df$id=="S1P nac"] <- "S1P [nM]"
df$id[df$id=="UK nac"] <- "UK [pM]"
levels(df$Condition)[match("ymptx",levels(df$Condition))] <- "YM+PTx"

### Split data
df0 = df %>% filter(Concen %in% c(0.02, 0.04, 4))
df1 = df %>% filter(!Concen %in% c(0.02, 0.04, 4))

# Get replicate ID within group and summary per replicate
df1 <- df1 %>% group_by(id, Condition) %>% group_split() %>% lapply(., function(x){transform(replicate = as.numeric(factor(Experiment)), x)}) %>% bind_rows()

df_sum0 <- df0 %>% group_by(Condition, Concen, id) %>% summarise(mean_ERK=mean(ERKnac), mean_Akt=mean(Aktnac))
df_sum0 <- do.call("rbind", list(df_sum0 %>% mutate(replicate = "1"), df_sum0 %>% mutate(replicate = "2"), df_sum0 %>% mutate(replicate = "3"), df_sum0 %>% mutate(replicate = "4")))
df_sum1 <- df1 %>% group_by(Condition, Concen, id, replicate=as.character(replicate)) %>% summarise(mean_ERK=mean(ERKnac), mean_Akt=mean(Aktnac))
df_summary <- rbind(df_sum0,df_sum1)

# Trying to color the cells (does not work)
# df1 <- df1 %>% mutate(replicate=as.character(replicate))
# df0 <- df0 %>% mutate(replicate="1")
# dfn <- rbind(df0, df1)
# dfn$replicate <- factor(dfn$replicate)

# Summary of all the cells (ignoring replicates)
df_sum1_1 <- df1 %>% group_by(Condition, Concen, id) %>% summarise(mean_ERK=mean(ERKnac), mean_Akt=mean(Aktnac))
df_summary1 <- rbind(df_sum0[-6],df_sum1_1)

## Plot the ERK data ----
# https://stackoverflow.com/questions/36780357/plotting-dose-response-curves-with-ggplot2-and-drc/37043751

ERKcrc = ggplot(data = df, aes(x = Concen, y = ERKnac)) + 
  geom_jitter(alpha=0.1, width=0.1,size=0.4) + 
  geom_point(data=df_summary, aes(x = Concen, y = mean_ERK, fill=replicate), size=8, shape=21, alpha=0.5) +
  
  # to get a SINGLE CURVE with average of all cells:
  # geom_smooth(data=df_summary1, aes(x = Concen, y = mean_ERK),color='black',method = drm, method.args = list(fct = L.4()), se = FALSE) +
  # to get a SINGLE CURVE with average of average replicate:
  geom_smooth(data=df_summary, aes(x = Concen, y = mean_ERK),color='black',method = drm, method.args = list(fct = L.4()), se = FALSE) +
  # to get curves PER REPLICATE:
  # geom_smooth(data=df_summary, aes(x = Concen, y = mean_ERK, fill=replicate), color="black", method = drm, method.args = list(fct = L.4()), se = FALSE) +
  
  scale_x_log10()+ylim(-1,10)+annotation_logticks(sides="b")+
  facet_grid(Condition~id, scales="free_x") +theme_light()+ ylab("AUC ERK") + xlab("ligand concentration")+
  theme(strip.text = element_text(size = 25, colour = "black"),  axis.text = element_text(size=20), axis.title = element_text(size=22), legend.key.size = unit(1.8,"line"), legend.title = element_text(size = 21), legend.text = element_text(size = 18))

# png(file=paste0("mean CRC ERK avg of all cells.png"), width = 1200, height = 1200, units = "px")
png(file=paste0("Figure_2.png"), width = 1200, height = 1200, units = "px")
# png(file=paste0("mean CRC ERK per replicate.png"), width = 1200, height = 1200, units = "px")
plot(ERKcrc)
dev.off()

### Plot the Akt data ----
# Aktcrc = ggplot(data = df, aes(x = Concen, y = Aktnac)) +
#   geom_jitter(alpha=0.1, width=0.1,size=0.4) + 
#   geom_point(data=df_summary, aes(x = Concen, y = mean_Akt, fill=replicate), size=8, shape=21, alpha=0.5) +
#   
#   # to get a SINGLE CURVE with average of all cells:
#   geom_smooth(data=df_summary1, aes(x = Concen, y = mean_Akt),color='black',method = drm, method.args = list(fct = L.4()), se = FALSE) +
#   # to get a SINGLE CURVE with average of average replicate:
#   # geom_smooth(data=df_summary, aes(x = Concen, y = mean_ERK),color='black',method = drm, method.args = list(fct = L.4()), se = FALSE) +
#   # to get curves PER REPLICATE:
#   # geom_smooth(data=df_summary, aes(x = Concen, y = mean_ERK, fill=replicate), color="black", method = drm, method.args = list(fct = L.4()), se = FALSE) +
#   
#   scale_x_log10()+ylim(-1,3)+annotation_logticks(sides="b")+
#   facet_grid(Condition~id, scales="free_x") +theme_light()+ ylab("AUC Akt") + xlab("ligand concentration")+
#   theme(strip.text = element_text(size = 25, colour = "black"),  axis.text = element_text(size=20), axis.title = element_text(size=22), legend.key.size = unit(1.8,"line"), legend.title = element_text(size = 21), legend.text = element_text(size = 18))
# 
# # png(file=paste0("mean CRC Akt avg of all cells.png"), width = 1200, height = 1200, units = "px")
# png(file=paste0("mean CRC Akt avg of avg.png"), width = 1200, height = 1200, units = "px")
# # png(file=paste0("mean CRC Akt per replicate.png"), width = 1200, height = 1200, units = "px")
# plot(Aktcrc)
# dev.off()

### Scatter plot Erk vs. Akt ----
# ggplot(data = df, aes(x = Aktnac, y = ERKnac))+ 
#   # geom_point(alpha=0.1)+
#   geom_hex() + scale_fill_viridis_c() +
#   # geom_bin2d() +
#   facet_grid(id~Condition)+xlim(-2,5)+ylim(-2,12)

## Define drm function (for Akt or ERK) ----
# https://gist.github.com/angelovangel/c69d78a27c4360df3057b40f4a705837
drm.func <- function(x) {
  drm(mean_ERK ~ Concen, 
      fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), 
      data = x)
}

# drm.func <- function(x) {
#   drm(mean_Akt ~ Concen, 
#       fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), 
#       data = x)
# }

drm.ci <- function(x) {confint(drm.func(x))}

coefs.fun <- function(x) {coef(x) %>% tidy}

## Define data for calculation of coefficients ----
## Coefficients for average of all cells
# coefs <- df_summary1 %>% group_by(id, Condition) %>% nest() %>%
#   mutate(drmod = map(data, drm.func), co = map(drmod, coefs.fun), ci = map(data, drm.ci)
#   )

## Coefficients for average of average with CI
coefs <- df_summary %>% group_by(id, Condition) %>% nest() %>% 
  mutate(drmod = map(data, drm.func), co = map(drmod, coefs.fun), ci = map(data, drm.ci)
  )

## Coefficients per replicate (execute without excluding the 3 experiments)
# df_summarym <- rbind(filter(df_summary, !replicate=="4"), filter(df_summary, replicate=="4", Condition=="ymptx", id=="S1P nac"))
# coefs <- df_summarym %>% group_by(id, Condition, replicate) %>% nest() %>%
#   mutate(drmod = map(data, drm.func), co = map(drmod, coefs.fun), ci = map(data, drm.ci)
#   )

## Transform CI intervals to tibble
coefs$ci = lapply(coefs$ci, as_tibble, rownames="names")

## Summarize coefficients ----
coefficients <- coefs %>% unnest(co) %>% spread(names,x)
coefficients <- data.frame(lapply(coefficients, as.character), stringsAsFactors=FALSE)
coefficients2 <- coefs %>% unnest(ci) %>% myspread(names, c('2.5 %', '97.5 %'))
coefficients2 <- data.frame(lapply(coefficients2, as.character), stringsAsFactors=FALSE)

allcoef <- merge(coefficients, coefficients2)
allcoefred <- subset(allcoef, select = -c(data, drmod, co, ci))
colnames(allcoefred) <- sub("Intercept", "", colnames(allcoefred))
allcoefred[3:14] <- as.numeric(unlist(allcoefred[3:14]))
allcoefred[3:14] <- round(allcoefred[3:14],2)

## Save coefficients ----
# write.csv(allcoef, "EC50 ERK coefficients avg of avg.csv", row.names=FALSE)
write.csv(allcoefred, "EC50 ERK coefficients avg of avg red.csv", row.names=FALSE)
# write.csv(allcoef, "EC50 ERK coefficients avg of all cells.csv", row.names=FALSE)
# write.csv(allcoef, "EC50 ERK coefficients per replicate.csv", row.names=FALSE)