################################################################################
##  APR_darkdiv_sppaccum.R: Species accumulation curves for Dark Diversity Network plots at APR.
##
##  Author: Kimberly Komatsu
##  Date created: February 5, 2020
################################################################################


##For loading your GIT auth token### library(gitcreds) then gitcreds::gitcreds_set()


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
setwd("~/Dropbox (University of Michigan)/APR_darkdiv") #eb mac
setwd("C://Users/alyou/Dropbox/APR_darkdiv") # Alyssa's laptop


 
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
#Effect of management depends on effect of species (native or non-native), only the interaction was significant (P value = .00706, F (1,36) = 8.163)
anova1posthoc <- TukeyHSD(anova1)
#a above cattle non-native, b above cattle native, ab above bison bars (or the other way)


###functional group
growthForm <- relCover%>%
  filter(provenance!='')%>% #filter out unknown species
  group_by(management, APR_plot_id, growth_form)%>%
  summarise(cover=sum(rel_cover))%>%
  ungroup()


#figure - percent native cover
ggplot(data=barGraphStats(data=growthForm, variable="cover", byFactorNames=c("management", "growth_form")), aes(x=growth_form, y=mean, fill=management)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  ylab('Relative Percent Cover\n') + xlab("")+
  scale_fill_manual(values=c("grey40", "grey"), labels=c("Bison", "Cattle"))+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  scale_x_discrete(limits = c("graminoid", "forb", "woody", "succulent"),breaks= c("graminoid", "forb", "woody", "succulent"),labels = c("Graminoid", "Forb", "Woody", "Succulent"))+
  theme(legend.position=c(.87,.9), axis.text.x=element_text(size=24, color = "black"))+
  expand_limits(y=80) +
  annotate("text", x= .77, y = 73, label= "a", size = 7)+ #gram bison
  annotate("text", x= 1.23, y = 75, label= "a", size = 7)+ #gram cattle
  annotate("text", x= 1.77, y = 40, label= "b", size = 7)+ #forb bison
  annotate("text", x= 2.23, y = 20, label= "bc", size = 7)+ #forb cattle
  annotate("text", x= 2.77, y = 20, label= "bc", size = 7)+ #woody bison
  annotate("text", x= 3.23, y = 32, label= "b", size = 7)+ #woody cattle
  annotate("text", x= 3.77, y = 6, label= "c", size = 7)+ #succ bison
  annotate("text", x= 4.23, y = 8, label= "c", size = 7) #succ cattle
#export at 600x600

anovafunctionalgroup <- aov(cover~management*growth_form, data = growthForm)
summary(anovafunctionalgroup)
#Growth form alone (P value = <.001, F(3, 71)= 79.732 and growth form's interaction with management (P value = 0.015, F(3, 71) = 3.731) are significant

anovafunctionalgroupposthoc <- TukeyHSD(anovafunctionalgroup)
#a for bison graminoid, cattle graminoid 
#b for bison forb, cattle woody
#c for bison succulent, cattle succulent 
#bc for cattle forb, bison woody 

#This is where we left off on 11/5/21

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
  ylab('Plant Species Richness\n') + xlab("") + expand_limits(y=40) +
  theme(axis.text.x=element_text(size=24, color = "black")) +
  scale_x_discrete(labels = c("Bison", "Cattle"))+
  annotate("text", x= 1, y = 33, label= "a", size = 7)+ #bison
  annotate("text", x= 2, y = 37, label= "b", size = 7) #cattle
#export at 400x600

##Richness t-test
richnessttest <- t.test(communityStructure$richness ~ communityStructure$management)
richnessttest
#t = -2.2123, df = 16.294, p-value = 0.04155
#mean bison: 28.2 
#mean cattle: 32.6 

#Richness by Provenance

provmanagment <- relCover%>%
  filter(provenance== "introduced"|provenance == "native")%>%
  select(management, APR_plot_id, provenance)%>%
  unique()

communityStructureprov <- provmanagement %>%
  group_by(APR_plot_id, provenance)%>%
  community_structure(relCover, time.var=NULL, abundance.var='rel_cover', replicate.var= 'APR_plot_id')


##Richness and Provenance pt2
Richnesses<-relCover%>%
  group_by(APR_plot_id,provenance)%>%
  mutate(richness=(count=unique(genus_species)))%>%
  group_by(APR_plot_id,provenance)%>%
  community_structure(time.var=NULL,abundance.var="rel_cover",replicate.var="APR_plot_id")


#anova richness and provenance
anovarichprov <- aov(richness~management*provenance, data = provmanagment)


#figure - evenness
ggplot(data=barGraphStats(data=communityStructure, variable="Evar", byFactorNames=c("management")), aes(x=management, y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  ylab('Plant Species Evenness\n') + xlab('') +  theme(axis.text.x=element_text(size=24, color = "black")) + expand_limits(y=.085) +
  scale_x_discrete(labels = c("Bison", "Cattle"))
#export at 400x600


#T-test for eveness 
evenessttest <- t.test(communityStructure$Evar ~ communityStructure$management)
evenessttest
#t = -1.5805, df = 12.847, p-value = 0.1383
#Bison mean: 0.06226895 
#Cattle mean:0.07236698 

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

