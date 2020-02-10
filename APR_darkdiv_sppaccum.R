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
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=15))

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
  scale_fill_manual(values=c("#FF9900", "#009900")) +
  ylab('Relative Percent Cover') + xlab('Management')
#export at 600x600



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
  # scale_fill_manual(values=c("#FF9900", "#009900")) +
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

