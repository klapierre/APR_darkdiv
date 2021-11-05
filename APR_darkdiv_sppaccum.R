################################################################################
##  APR_darkdiv_sppaccum.R: Species accumulation curves for Dark Diversity Network plots at APR.
##
##  Author: Kimberly Komatsu
##  Date created: February 5, 2020
################################################################################

library(codyn)
library(vegan)
library(tidyverse)

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=24, vjust=-0.35), axis.text.x=element_text(size=20, color = "black"),
             axis.title.y=element_text(size=24, angle=90, vjust=0.7), axis.text.y=element_text(size=20, color= "black"),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=18), panel.border=element_rect(color="black", fill = NA, size = 1)) #axis.line=element_line(color="black")

###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))
barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  

setwd('C:\\Users\\komatsuk\\Dropbox (Smithsonian)\\American Prairie Reserve\\APR_darkdiv') #desktop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\American Prairie Reserve\\APR_darkdiv') #laptop
setwd("~/Dropbox (Smithsonian)/APR_darkdiv") #shelley
setwd("~/Documents/r_stuff/APR_darkdiv") #skye's mac
setwd("C:/Users/hrusk/Dropbox (Smithsonian)/APR_darkdiv") #amy's laptop
setwd("C:/Users/Sarah Alley/Dropbox (Smithsonian)/APR_darkdiv") #Sarah's laptop

 
###read in data
spp <- read.csv('DarkDivNet_Species_Comp_Traits_2019_forAnalysis.csv')%>%
  mutate(management=ifelse(location=='BLM', 'cattle', 'bison'))

spp10 <- spp%>%
  filter(!darkdiv_plot_id %in% c('A1S', 'N1S')) #filter down to just 10x10 m plots
  
#calculate relative cover
totCover <- spp10%>%
  group_by(APR_plot_id)%>%
  summarise(total_cover=sum(cover))%>%
  ungroup()

relCover <- spp10%>%
  left_join(totCover)%>%
  mutate(rel_cover=100*(cover/total_cover))



###percent native
provenance <- relCover%>%
  filter(provenance!='')%>% #filter out unknown species
  group_by(management, APR_plot_id, provenance)%>%
  summarise(cover=sum(rel_cover))%>%
  ungroup()

#figure - percent native cover
ggplot(data=barGraphStats(data=provenance, variable="cover", byFactorNames=c("management", "provenance")), aes(x=management, y=mean, fill=provenance)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  scale_fill_manual(values=c("grey40", "grey"), labels=c("Non-native", "Native")) +
  ylab('Relative Percent Cover\n')+ xlab(element_blank())+
  expand_limits(y=80)+
  theme(legend.position=c(0.21,.93), axis.text.x=element_text(size=24, color = "black"))+
  scale_x_discrete(labels = c("Bison", "Cattle"))+
  annotate("text", x= 2.23, y = 82, label= "a", size = 7)+ 
  annotate("text", x= 1.77, y = 45, label= "b", size = 7)+
  annotate("text", x= 0.77, y = 68, label= "ab", size = 7)+
  annotate("text", x= 1.23, y = 60, label= "ab", size = 7)
#export at 600x600

#anova for percent cover
anova1 <- aov(cover~management*provenance, data = provenance)
summary(anova1)
#Effect of management depends on effect of species (native or non-native), only the interaction was significant (P value = .00706, F 1,36 = 8.163)
anova1posthoc <- TukeyHSD(anova1)
#a above cattle non-native, b above cattle native, ab above bison bars (or the other way)


###functional group
growthForm <- relCover%>%
  filter(provenance!='')%>% #filter out unknown species
  group_by(management, APR_plot_id, growth_form)%>%
  summarise(cover=sum(rel_cover))%>%
  ungroup()

#figure - percent native cover
ggplot(data=barGraphStats(data=growthForm, variable="cover", byFactorNames=c("management", "growth_form")), aes(x=management, y=mean, fill=growth_form)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  #scale_fill_manual(values=c("#FF9900", "#009900")) + ##need 2 more values##
  ylab('Relative Percent Cover') + xlab('Management')
#export at 600x600




###species richness, evenness
management <- relCover%>%
  select(management, APR_plot_id)%>%
  unique()

communityStructure <- community_structure(relCover, time.var=NULL, abundance.var='rel_cover', replicate.var='APR_plot_id')%>%
  left_join(management)

#figure - richness
ggplot(data=barGraphStats(data=communityStructure, variable="richness", byFactorNames=c("management")), aes(x=management, y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  ylab('Plant Species Richness') + xlab('Management')
#export at 400x600

#figure - evenness
ggplot(data=barGraphStats(data=communityStructure, variable="Evar", byFactorNames=c("management")), aes(x=management, y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  ylab('Plant Species Evenness') + xlab('Management')
#export at 400x600



###community differences
#make matrix
sppMatrix <- relCover%>%
  select(location, APR_plot_id, genus_species, rel_cover)%>%
  spread(key=genus_species, value=rel_cover, fill=0)

sppBC <- metaMDS(sppMatrix[,3:85])

plots <- 1:nrow(sppMatrix)
plotData <- sppMatrix[,1:2]
plot(sppBC$points,col=as.factor(plotData$location))
ordiellipse(sppBC, groups = as.factor(plotData$Watershed), kind = "sd", display = "plots", label = T)

#Use the vegan ellipse function to make ellipses           
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

BC_NMDS = data.frame(MDS1 = sppBC$points[,1], MDS2 = sppBC$points[,2],group=plotData$location)
BC_NMDS_Graph <- cbind(plotData,BC_NMDS)
BC_Ord_Ellipses<-ordiellipse(sppBC, plotData$location, display = "sites",
                             kind = "se", conf = 0.95, label = T)               





test

