################################################################################
##  APR_darkdiv_sppaccum.R: Species accumulation curves for Dark Diversity Network plots at APR.
##
##  Author: Kimberly Komatsu, Shelley Bennett, Sarah Alley, Elizabeth Blake, Kathryn Bloodworth, Rachael Brenneman, Alyssa Young
##  Date created: February 5, 2020
################################################################################


##For loading your GIT auth token### library(gitcreds) then gitcreds::gitcreds_set()


library(codyn)
library(vegan)
library(grid)
library(tidyverse)
library(nlme)
library(factoextra)
library(ggfortify)
library(fuzzySim) #spCode function
library(ggord) #plotting ordinations

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
setwd("~/Dropbox/APR_darkdiv") #rachael's mac


 
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

#made new data frame to change genus and species to codes without disrupting initial analyses
relCover2 <-relCover %>%
  separate(genus_species, c("genus", "species"), sep = "_")%>%
  mutate(sppName = paste(genus,species,sep=" "))
  
relCover2$sppName <- spCodes(relCover2$sppName, nchar.gen = 2, nchar.sp = 2, nchar.ssp = 0, sep.species = " ", sep.spcode = "", verbosity = 2)

relCover3 <-relCover %>%
  separate(genus_species, c("genus", "species"), sep = "_")%>%
  mutate(sppName = paste(genus,species,sep=" "))

#Histograms to check for outliers
#Skewed to the right but no obvious outliers 
#Alyssum desertorum oddly high? 
hist(relCover$cover)
hist(relCover$rel_cover)


###percent native
provenance <- relCover%>%
  filter(provenance!='')%>% #filter out unknown species
  group_by(management, APR_plot_id, provenance)%>%
  summarise(cover=sum(rel_cover))%>%
  ungroup()

hist(provenance$cover)

#figure - percent native cover
ggplot(data=barGraphStats(data=provenance, variable="cover", byFactorNames=c("management", "provenance")), aes(x=management, y=mean, fill=provenance)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  scale_fill_manual(values=c("grey", "grey40"), labels=c("Non-native", "Native")) +
  ylab('Relative Percent Cover\n')+ xlab(element_blank())+
  expand_limits(y=80)+
  theme(legend.position=c(0.21,.93), axis.text.x=element_text(size=24, color = "black"))+
  scale_x_discrete(labels = c("Sun Prairie", "BLM"))+
  annotate("text", x= 2.23, y = 82, label= "a", size = 7)+ 
  annotate("text", x= 1.77, y = 45, label= "b", size = 7)+
  annotate("text", x= 0.77, y = 68, label= "ab", size = 7)+
  annotate("text", x= 1.23, y = 60, label= "ab", size = 7)
#export at 600x600


# ###percent native by growth form
# provenanceGrowthForm <- relCover%>%
#   filter(provenance!='')%>% #filter out unknown species
#   group_by(management, APR_plot_id, provenance, growth_form)%>%
#   summarise(cover=sum(rel_cover))%>%
#   ungroup()
# 
# #figure - percent native cover by growth form
# ggplot(data=barGraphStats(data=provenanceGrowthForm, variable="cover", byFactorNames=c("management", "provenance", "growth_form")), aes(x=growth_form, y=mean, fill=provenance)) +
#   geom_bar(stat='identity', position=position_dodge()) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
#   # scale_fill_manual(values=c("grey40", "grey"), labels=c("Non-native", "Native")) +
#   ylab('Relative Percent Cover\n')+ xlab(element_blank())+
#   # expand_limits(y=80)+
#   # theme(legend.position=c(0.21,.93), axis.text.x=element_text(size=24, color = "black"))+
#   # scale_x_discrete(labels = c("Bison", "Cattle"))+
#   # annotate("text", x= 2.23, y = 82, label= "a", size = 7)+ 
#   # annotate("text", x= 1.77, y = 45, label= "b", size = 7)+
#   # annotate("text", x= 0.77, y = 68, label= "ab", size = 7)+
#   # annotate("text", x= 1.23, y = 60, label= "ab", size = 7) +
#   facet_wrap(~management)
# #export at 600x600
# 
# #anova for percent cover by provenance and growth form
# anova1 <- aov(cover~management*provenance*growth_form, data = provenanceGrowthForm)
# summary(anova1) #no three-way interaction so no need to pursue this


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

hist(sqrt(growthForm$cover))
hist(growthForm$cover)

#figure - percent native cover
ggplot(data=barGraphStats(data=growthForm, variable="cover", byFactorNames=c("management", "growth_form")), aes(x=growth_form, y=mean, fill=management)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  ylab('Relative Percent Cover\n') + xlab("")+
  scale_fill_manual(values=c("grey40", "grey"), labels=c("Sun Prairie", "BLM"))+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  scale_x_discrete(limits = c("graminoid", "forb", "woody", "succulent"),breaks= c("graminoid", "forb", "woody", "succulent"),labels = c("Graminoid", "Forb", "Woody", "Succulent"))+
  theme(legend.position=c(.82,.9), axis.text.x=element_text(size=24, color = "black"))+
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

anovafunctionalgroup <- aov(sqrt(cover)~management*growth_form, data = growthForm)
summary(anovafunctionalgroup)
#Growth form alone (P value = <.001, F(3, 71)= 79.732 and growth form's interaction with management (P value = 0.015, F(3, 71) = 3.731) are significant
#sqrt more normal
#sqrt growth form (P value = < 2e-16, F(3, 71) = 80.218 and growth form's interaction with management (P value = 0.00129, F(3, 71 = 5.825) are significant

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

hist(communityStructure$richness)
hist(communityStructure$Evar)

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

### AY & KB figured out how to separate richness by provenance and ran the stats and made figure 
#Richness by Provenance and Management
relCoverProv <- spp10%>%
  left_join(totCover) %>%
  filter(provenance== "introduced"|provenance == "native") %>%
  mutate(ID_provenance = paste(APR_plot_id,provenance,sep="-")) %>%
  mutate(rel_cover=100*(cover/total_cover))

provmanagment <- relCoverProv%>%
  select(ID_provenance,management)%>%
  unique()

communityStructureprov<-community_structure(relCoverProv, time.var=NULL, abundance.var='rel_cover', 
  replicate.var= 'ID_provenance') %>%
  left_join(provmanagment) %>%
  separate(ID_provenance,c("APR_plot_id","provenance"),sep="-")



#figure - richness, management, and provenance
ggplot(data=barGraphStats(data=communityStructureprov, variable="richness", byFactorNames=c("management","provenance")), aes(x=management, y=mean,fill=provenance)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  theme(legend.position=c(.2,.9), axis.text.x=element_text(size=24, color = "black"))+
  ylab('Plant Species Richness\n') + xlab("") + expand_limits(y=30)+
  scale_fill_manual(values=c("grey", "grey40"), labels=c("Non-native", "Native"))+
  theme(axis.text.x=element_text(size=24, color = "black")) +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  annotate("text", x= 0.77, y = 11, label= "a", size = 7)+ #bison introduced
  annotate("text", x= 1.23, y = 23, label= "b", size = 7)+ #bison native
  annotate("text", x= 1.77, y = 10, label= "a", size = 7)+ #cattle introduced
  annotate("text", x= 2.23, y = 28, label= "c", size = 7) #cattle native
#export at 600x600


#anova richness and provenance
anovarichprov <- aov(richness~management*provenance, data = communityStructureprov)
summary(anovarichprov)
#                        Df Sum Sq Mean Sq F value   Pr(>F)    
# management             1   48.4    48.4   4.324  0.04477 *  
# provenance             1 1988.1  1988.1 177.597 1.73e-15 ***
# management:provenance  1  129.6   129.6  11.577  0.00165 ** 
# Residuals             36  403.0    11.2           

anovaRichProvManageposthoc <- TukeyHSD(anovarichprov)
#a for bison introduced, cattle introduced
#b for bison native
#c for cattle native




#figure - evenness
ggplot(data=barGraphStats(data=communityStructure, variable="Evar", byFactorNames=c("management")), aes(x=management, y=mean)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  ylab('Plant Species Evenness\n') + xlab('') +  theme(axis.text.x=element_text(size=24, color = "black")) + expand_limits(y=.085) +
  scale_x_discrete(labels = c("Sun\nPrairie", "BLM"))
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

sppBC <- metaMDS(sppMatrix[,3:87])

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

ord3 <- data.frame(plotData,scores(sppBC,display="sites"))%>%
  group_by(location)

BC_Ord_Ellipses<-ordiellipse(sppBC, plotData$location, display = "sites",
                            kind = "se", conf = 0.95, label = T)
#Make a new empty data frame called BC_Ellipses                
BC_Ellipses <- data.frame()
#Generate ellipses points
for(g in unique(BC_NMDS$group)){
  BC_Ellipses <- rbind(BC_Ellipses, cbind(as.data.frame(with(BC_NMDS[BC_NMDS$group==g,], 
  veganCovEllipse(BC_Ord_Ellipses[[g]]$cov,BC_Ord_Ellipses[[g]]$center,BC_Ord_Ellipses[[g]]$scale)))
  ,group=g))
}



#Make a data frame called BC_NMDS and at a column using the first set of "points" in BC_Data and a column using the second set of points.  Group them by watershed
#BC_NMDS = data.frame(MDS1 = BC_Data$points[,1], MDS2 = BC_Data$points[,2],group=BC_Meta_Data$Watershed)

ggplot(BC_NMDS_Graph, aes(x=MDS1, y=MDS2, color=group,linetype = group, shape = group))+
  geom_point(size=6)+ 
  geom_path(data = BC_Ellipses, aes(x = NMDS1, y = NMDS2), size = 3)+
  labs(color="", linetype = "", shape = "")+
  scale_colour_manual(values=c("grey40", "grey"), labels = c("Bison", "Cattle"), name = "")+
  scale_linetype_manual(values = c("twodash", "solid"), labels = c("Bison", "Cattle"), name = "")+
  scale_shape_manual(values = c(15, 16),labels = c("Bison", "Cattle"), name = "")+
  xlab("NMDS1")+ 
  ylab("NMDS2")+ 
  theme(axis.text.x=element_text(size=24, color = "black"), axis.text.y = element_text(size = 24, color = "black"), legend.text = element_text(size = 24))

##permanova data
PermanovaData <- sppMatrix %>%
  select(-APR_plot_id, -location)
Permanova <- adonis2(formula = PermanovaData~ location, data = plotData, permutations = 999, method= "bray")
print(Permanova)
#F= 5.7566   DF = 1,19  P = .001

#betadisper
Veg <- vegdist(PermanovaData, method = "bray")
Dispersion <- betadisper(Veg, plotData$location)
permutest(Dispersion, pairwise = TRUE, permutations = 999) 

#species accumulation

coversection <- spp%>%
  mutate(plotIDsize = paste(APR_plot_id, plot_section, sep = "::"))

plotsection <- coversection%>%
  select(plotIDsize)%>%
  unique()

communitystructuresection<- community_structure(coversection, time.var=NULL, abundance.var='cover', replicate.var='plotIDsize')%>%
  left_join(plotsection)%>%
  separate(plotIDsize, c("APR_plot_id", "plot_section"), sep = "::")

specaccum <- read.csv("DarkDivSpecAccum.csv")

ggplot(data=barGraphStats(data=specaccum, variable="richness", byFactorNames=c("plot_size", "management")), aes(x=plot_size, y=mean, fill = management, colour = management))+
  geom_line(aes(linetype = management), size = 1.5) +
  geom_point(size = 2, color = "black")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9), color = "black")+
  ylab('Average Species Richness\n')+ xlab(bquote('Plot Size'~(m^2)))+
  theme(axis.title=element_text(size=24, color = "black"),legend.position=c(.2,.9))+
  scale_colour_manual(values=c("grey", "grey40"), labels=c("BLM", "Sun Prairie"))+
  expand_limits(y= 35)



                                         

###rank abundance curves
rankAbundance <- relCover3%>%
  group_by(location, provenance, growth_form, sppName)%>%
  summarize(avg_cover=mean(rel_cover))%>%
  ungroup()%>%
  filter(sppName!='unknown:small thin')%>%
  mutate(provenance=ifelse(sppName=='Penstemon albicula', 'native', as.character(provenance)))%>%
  arrange(location, -avg_cover)%>%
  group_by(location)%>%
  mutate(rank=seq_along(location))%>%
  ungroup()

BLMrank <- ggplot(data=subset(rankAbundance, location=='BLM', avg_cover>0), aes(x=rank, y=avg_cover)) +
  geom_line() +
  geom_point(aes(colour=provenance, shape=growth_form), size=3) +
  scale_color_manual(values=c("grey50", "grey20"), labels=c("Non-native", "Native"))+
  scale_shape_discrete(labels=c("Forb","Graminoid", "Succulent", "Woody"),name="Growth Form")+
  xlab('') +
  ylab('BLM\nRelative Percent Cover\n') +
  # scale_x_continuous(expand=c(0,0), limits=c(0.5,17), breaks=seq(0,17,5)) +
  # scale_y_continuous(expand=c(0,0), limits=c(0,60), breaks=seq(0,60,10)) +
  geom_text(aes(y=avg_cover+1.2, x=rank+0.1, label=sppName), hjust='left', vjust='center', angle=90, size=4)+
  expand_limits(y=40)
BLMrank

BDrank <- ggplot(data=subset(rankAbundance, location=='BD', avg_cover>0), aes(x=rank, y=avg_cover)) +
  geom_line() +
  geom_point(aes(colour=provenance, shape=growth_form), size=3) +
  scale_color_manual(values=c("grey50", "grey20"), labels=c("Non-native", "Native"))+
  scale_shape_discrete(labels=c("Forb","Graminoid", "Succulent", "Woody"),name="Growth Form")+
  xlab('Species Rank') +
  ylab('Sun Prairie\nRelative Percent Cover\n') +
  # scale_x_continuous(expand=c(0,0), limits=c(0.5,17), breaks=seq(0,17,5)) +
  # scale_y_continuous(expand=c(0,0), limits=c(0,60), breaks=seq(0,60,10)) +
  geom_text(aes(y=avg_cover+1.5, x=rank+0.1, label=sppName), hjust='left', vjust='center', angle=90, size=4)+
  expand_limits(y=50)
BDrank
pushViewport(viewport(layout=grid.layout(2,1)))
print(BLMrank, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(BDrank, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))





###soil and precipitation data
covar <- read.csv('Kim_covar.csv') %>%
  mutate(location = ifelse(Treatment == "Bison_d", "BD", "BLM")) %>%
  mutate(APR_plot_id = paste(location,id,sep="_"))

covarRelCov <- relCover %>%
  left_join(covar) %>%
  filter(!is.na(elevation_)) %>%
  select(APR_plot_id, growth_form, provenance, management, rel_cover, soil_wat_2, soil_san_2, soil_Cla_2, soil_bul_2, slope_WGS8, elevation_, aspect_WGS)


#PCA
RelCovPCA<- prcomp(covarRelCov[,6:12]) #prcomp is the function that does a PCA on plantWide dataframe, removing column 1 (which was the plot ids)
                                               # PC1      PC2     PC3     PC4     PC5     PC6     PC7
                        #Standard deviation     76.8056 15.94035 2.89597 2.09064 1.38459 0.84450 0.55555
                        #Proportion of Variance  0.9563  0.04119 0.00136 0.00071 0.00031 0.00012 0.00005
                        #Cumulative Proportion   0.9563  0.99746 0.99882 0.99952 0.99983 0.99995 1.00000
RelCovAxes <- predict(RelCovPCA, newdata=covarRelCov)%>% #makes a new dataframe called plantAxes and puts together the plantPCA output with the plantWide dataframe
  cbind(covarRelCov[,1:5])%>% #binds on the plot ids (column 1 of the plantWide dataframe)
  select(APR_plot_id, growth_form, provenance, management, rel_cover, PC1, PC2) #keeps only the plot ids, and eigenvalues from the first two PC axes (could modify this to keep more axes too)


##all soil, elevation, slope, aspect, precip, temp into PCA, animal stuff separate
##or everything into one big PCA
##out of pca, see what % variation explained by first several axes. if 90% use first axes, or 50, 30 use 2 pca numbers
## anovas for richness and evenness, do mixed model instead and include pca as random factor using lme4

summary(APRModel <- lme(rel_cover~as.factor(management)*as.factor(provenance),data=RelCovAxes, random=~1|PC1))
                           
var <- get_pca_var(RelCovPCA)
var
head(var$contrib,7)

#PCA Graph

g<-autoplot(RelCovPCA,data=covarRelCov,scale=0,
            colour="Site",loadings=TRUE,loadings.colour="black",size=3,
            loadings.label=TRUE,loadings.label.colour="black",loadings.label.size=6)

arrow_ends <- layer_data(g, 2)[,c(2,4)]

autoplot(RelCovPCA,data=covarRelCov,scale=0,
         colour="management",loadings=TRUE,loadings.colour="black",size=3,
         loadings.label=TRUE,loadings.label.colour="black",loadings.label.size=5,
         loadings.label.vjust = 1.5,frame=T,frame.colour = 'management') +
  geom_point(size = 2) + 
  theme(plot.background=element_blank(),
        panel.background=element_rect(fill='transparent',color='black',size=1),
        legend.key=element_blank())


###soil and precip 2 with all data
covar2 <- read.csv('Kim_covar03312022.csv') %>%
  mutate(location = ifelse(Treatment == "Bison_d", "BD", "BLM")) %>%
  mutate(APR_plot_id = paste(location,id,sep="_")) 

covarRelCov2 <- relCover %>%
  left_join(covar2) %>%
  filter(!is.na(elevation_)) %>%
  select(APR_plot_id, growth_form, provenance, management, rel_cover, soil_wat_2, soil_san_2, soil_Cla_2, soil_bul_2, slope_WGS8, elevation_, aspect_WGS, prcp_reduc, distwaterS, distwaterP)

#standardize covariates into z-scores
covarRelCov2standard <- covarRelCov2 %>%
  mutate(soilwat2_stand = (soil_wat_2 - mean(soil_wat_2, na.rm = TRUE)) / sd(soil_wat_2, na.rm = TRUE)) %>%
  mutate(soilsan2_stand = (soil_san_2 - mean(soil_san_2, na.rm = TRUE)) / sd(soil_san_2, na.rm = TRUE)) %>%
  mutate(soilCla2_stand = (soil_Cla_2 - mean(soil_Cla_2, na.rm = TRUE)) / sd(soil_Cla_2, na.rm = TRUE)) %>%
  mutate(soilbul2_stand = (soil_bul_2 - mean(soil_bul_2, na.rm = TRUE)) / sd(soil_bul_2, na.rm = TRUE)) %>%
  mutate(slope_stand = (slope_WGS8 - mean(slope_WGS8, na.rm = TRUE)) / sd(slope_WGS8, na.rm = TRUE))  %>%
  mutate(elevation_stand = (elevation_ - mean(elevation_, na.rm = TRUE)) / sd(elevation_, na.rm = TRUE)) %>%
  mutate(aspect_stand = (aspect_WGS - mean(aspect_WGS, na.rm = TRUE)) / sd(aspect_WGS, na.rm = TRUE)) %>%
  mutate(prcp_stand = (prcp_reduc - mean(prcp_reduc, na.rm = TRUE)) / sd(prcp_reduc, na.rm = TRUE)) %>%
  mutate(distwaterS_stand = (distwaterS - mean(distwaterS, na.rm = TRUE)) / sd(distwaterS, na.rm = TRUE)) %>%
  mutate(distwaterP_stand = (distwaterP - mean(distwaterP, na.rm = TRUE)) / sd(distwaterP, na.rm = TRUE))
 
RelCov2StandardPCA <- prcomp(covarRelCov2standard[,16:25])
RelCov2StandardAxes <- predict(RelCov2StandardPCA, newdata=covarRelCov2standard)%>% #makes a new dataframe called plantAxes and puts together the plantPCA output with the plantWide dataframe
  cbind(covarRelCov2standard[,1:5])%>% #binds on the plot ids (column 1 of the plantWide dataframe)
  select(APR_plot_id, growth_form, provenance, management, rel_cover, PC1, PC2, PC3)

g3<-autoplot(RelCov2StandardPCA,data=covarRelCov2standard,scale=0,
             colour="management",loadings=TRUE,loadings.colour="black",size=3,
             loadings.label=TRUE,loadings.label.colour="black",loadings.label.size=6)

arrow_ends <- layer_data(g3, 3)[,c(2,4)]

autoplot(RelCov2StandardPCA,data=covarRelCov2standard,scale=0,
         colour="management",loadings=TRUE,loadings.colour="black",size=3,
         loadings.label=TRUE,loadings.label.colour="black",loadings.label.size=5,
         loadings.label.vjust = 1.5,frame=T,frame.colour = 'management') +
  geom_point(size = 2) + 
  theme(plot.background=element_blank(),
        panel.background=element_rect(fill='transparent',color='black',size=1),
        legend.key=element_blank())

##ask Alyssa how to change plot to include 3 PC columns


#PCA
RelCovPCA2<- prcomp(covarRelCov2[,6:15])

RelCovAxes2 <- predict(RelCovPCA2, newdata=covarRelCov2)%>% #makes a new dataframe called plantAxes and puts together the plantPCA output with the plantWide dataframe
  cbind(covarRelCov[,1:5])%>% #binds on the plot ids (column 1 of the plantWide dataframe)
  select(APR_plot_id, growth_form, provenance, management, rel_cover, PC1, PC2)

g2<-autoplot(RelCovPCA2,data=covarRelCov2,scale=0,
             colour="management",loadings=TRUE,loadings.colour="black",size=3,
             loadings.label=TRUE,loadings.label.colour="black",loadings.label.size=6)

arrow_ends <- layer_data(g2, 2)[,c(2,4)]

autoplot(RelCovPCA2,data=covarRelCov2,scale=0,
         colour="management",loadings=TRUE,loadings.colour="black",size=3,
         loadings.label=TRUE,loadings.label.colour="black",loadings.label.size=5,
         loadings.label.vjust = 1.5,frame=T,frame.colour = 'management') +
  geom_point(size = 2) + 
  theme(plot.background=element_blank(),
        panel.background=element_rect(fill='transparent',color='black',size=1),
        legend.key=element_blank())




##RDA
# The matrix of species composition (sample x species)--> in our case plot x species
plots <- read.csv("plot.csv")
plotsenv <- read.csv("plotenv.csv")

spe <- relCover2%>%
  merge(plots, by = ("APR_plot_id"))%>%
  select(plot, sppName, rel_cover)%>%
  spread(key=sppName, value=rel_cover, fill=0)

spe[-1] <- log1p(spe[-1]) # species data are in percentage scale which is strongly rightskewed, better to transform them

spe.hell <- decostand (spe[-1], 'hell') #but they're already relative so we just need to do a square root: you've only taken a subset of the species, then yes, you can just apply a square root transformation to the data you are using and it would have been the same if you'd done the entire Hellinger transformation on the entire data set and then thrown out some of the species.


#and the matrix of environmental variables (sample x env.variables, for simplicity containing only one env. variable in the illustration below) needs to be available --> in our case plot x env. variables (covar2) --> covarRDA

envRDA <- covar2 %>%
  filter(!is.na(elevation_)) %>%
  merge(plotsenv, by = ("APR_plot_id"))%>%
  mutate(location = ifelse(Treatment == "Bison_d", "Sun Prairie", "BLM"))%>%
  select(plot, location, soil_wat_2, soil_san_2, soil_Cla_2, soil_bul_2, slope_WGS8, elevation_, aspect_WGS, prcp_reduc, distwaterS, distwaterP)%>%
  mutate(Soil.Water = (soil_wat_2 - mean(soil_wat_2, na.rm = TRUE)) / sd(soil_wat_2, na.rm = TRUE)) %>%
  mutate(Sand = (soil_san_2 - mean(soil_san_2, na.rm = TRUE)) / sd(soil_san_2, na.rm = TRUE)) %>%
  mutate(Clay = (soil_Cla_2 - mean(soil_Cla_2, na.rm = TRUE)) / sd(soil_Cla_2, na.rm = TRUE)) %>%
  mutate(Bulk.Density = (soil_bul_2 - mean(soil_bul_2, na.rm = TRUE)) / sd(soil_bul_2, na.rm = TRUE)) %>%
  mutate(Slope = (slope_WGS8 - mean(slope_WGS8, na.rm = TRUE)) / sd(slope_WGS8, na.rm = TRUE))  %>%
  mutate(Elev. = (elevation_ - mean(elevation_, na.rm = TRUE)) / sd(elevation_, na.rm = TRUE)) %>%
  mutate(Aspect = (aspect_WGS - mean(aspect_WGS, na.rm = TRUE)) / sd(aspect_WGS, na.rm = TRUE)) %>%
  mutate(Precip. = (prcp_reduc - mean(prcp_reduc, na.rm = TRUE)) / sd(prcp_reduc, na.rm = TRUE)) %>%
  mutate(Dist.Water = (distwaterP - mean(distwaterP, na.rm = TRUE)) / sd(distwaterP, na.rm = TRUE))

env <- envRDA %>%
  select(plot, location, Soil.Water, Sand, Clay , Bulk.Density, Slope, Elev., Aspect, Precip., Dist.Water)

#run the tbRDA

tbRDA <- rda(spe.hell ~ Soil.Water + Sand + Clay + Bulk.Density + Slope + Elev. + Aspect + Precip. + Dist.Water, data = env)
tbRDA
summary(tbRDA)


# Call:
#   rda(formula = spe.hell ~ Soil.Water + Sand + Clay + BD + Slope +      Elev + Aspect + Precip + DW, data = env) 
# 
# Partitioning of variance:
#   Inertia Proportion
# Total          0.4549     1.0000
# Constrained    0.2836     0.6235
# Unconstrained  0.1713     0.3765
# 
# Eigenvalues, and their contribution to the variance 
# 
# Importance of components:
#   RDA1    RDA2    RDA3    RDA4    RDA5    RDA6    RDA7     RDA8
# Eigenvalue            0.1229 0.04156 0.02955 0.02761 0.02010 0.01513 0.01148 0.009667
# Proportion Explained  0.2701 0.09136 0.06496 0.06068 0.04419 0.03326 0.02523 0.021250
# Cumulative Proportion 0.2701 0.36150 0.42646 0.48714 0.53133 0.56459 0.58982 0.611070
# RDA9     PC1     PC2     PC3     PC4     PC5     PC6      PC7
# Eigenvalue            0.005636 0.04935 0.03707 0.02162 0.01857 0.01606 0.01431 0.008698
# Proportion Explained  0.012388 0.10847 0.08148 0.04752 0.04082 0.03531 0.03146 0.019119
# Cumulative Proportion 0.623459 0.73193 0.81341 0.86093 0.90175 0.93706 0.96852 0.987642
# PC8
# Eigenvalue            0.005622
# Proportion Explained  0.012358
# Cumulative Proportion 1.000000
# 
# Accumulated constrained eigenvalues
# Importance of components:
#   RDA1    RDA2    RDA3    RDA4    RDA5    RDA6    RDA7     RDA8
# Eigenvalue            0.1229 0.04156 0.02955 0.02761 0.02010 0.01513 0.01148 0.009667
# Proportion Explained  0.4333 0.14655 0.10419 0.09734 0.07088 0.05335 0.04047 0.034084
# Cumulative Proportion 0.4333 0.57983 0.68402 0.78135 0.85223 0.90558 0.94605 0.980130
# RDA9
# Eigenvalue            0.005636
# Proportion Explained  0.019870
# Cumulative Proportion 1.000000
# 
# Scaling 2 for species and site scores
# * Species are scaled proportional to eigenvalues
# * Sites are unscaled: weighted dispersion equal on all dimensions
# * General scaling constant of scores:  1.667611 
# 
# 
# Species scores
# 
# RDA1       RDA2       RDA3       RDA4       RDA5       RDA6
# Acmi  0.0110733 -2.843e-02  0.0115349  0.0151611  0.0156498  0.0127787
# Agcr -0.4320868 -5.441e-02  0.0412808  0.0130168 -0.0411117  0.0515534
# Alde -0.3073684  9.432e-02 -0.0287679 -0.0363444  0.0324119 -0.0526226
# Alte  0.0052866 -3.012e-03 -0.0054221  0.0012174  0.0044188  0.0012623
# Amps -0.0015772  9.698e-04 -0.0016251 -0.0014441  0.0011822 -0.0027385
# Anoc -0.0009110 -6.181e-04  0.0004925 -0.0006903 -0.0019193  0.0011783
# Anpa  0.0244064 -4.715e-04  0.0479356  0.0338196  0.0008608 -0.0257131
# Arca -0.0151944 -9.411e-02  0.0304930 -0.0349998 -0.0498455  0.0537468
# Arfr  0.1934118 -6.549e-02 -0.1194745 -0.0167128  0.0538539  0.0300606
# Artr  0.2844432 -1.182e-01  0.0098994 -0.0223958  0.0070874 -0.0363175
# Aser  0.0174739 -5.331e-03  0.0372337  0.0239434 -0.0125237 -0.0215693
# Aslo -0.0016387 -2.226e-03  0.0019214 -0.0009317 -0.0001139 -0.0050379
# Asmi  0.0302153 -1.577e-03  0.0261905 -0.0181154 -0.0156068 -0.0011792
# Atca  0.0204771  7.427e-03  0.0133913 -0.0067413  0.0319732  0.0346867
# Basc -0.0028443 -2.845e-04 -0.0019641 -0.0028757  0.0012806 -0.0012527
# Bogr  0.1282145 -1.915e-01 -0.0068221 -0.0726381 -0.0206911  0.0872940
# Brar  0.0832499 -1.722e-02 -0.1720797 -0.0500185  0.0824663 -0.0261518
# Brte -0.0649315 -8.353e-02  0.0403815  0.1357085  0.0410827  0.0271055
# Buda  0.0011470 -2.715e-04 -0.0037777  0.0019176 -0.0034334 -0.0015002
# Cadu  0.0021484  2.632e-03  0.0023912  0.0020323  0.0125362 -0.0066429
# Cafi  0.0003532  3.397e-04 -0.0001984  0.0013383  0.0001992  0.0020244
# Cami  0.0015265  7.066e-04  0.0002438 -0.0020955  0.0013699  0.0007210
# Cast  0.0088188 -4.919e-02  0.0273245  0.0283111  0.0067779  0.0148905
# Chal -0.0464607  4.539e-02  0.0141791 -0.0253006 -0.0046894  0.0212380
# Chvi  0.1530222  1.522e-02  0.0616704  0.0137624  0.0529698 -0.0091018
# Coum  0.1359084  9.154e-02 -0.0522847  0.0453866 -0.0866489 -0.0081283
# Covi -0.0058359 -3.376e-03  0.0008743 -0.0083795  0.0023129  0.0206099
# Crmo -0.0043470  3.972e-03 -0.0007122 -0.0033249  0.0008090 -0.0014124
# Daca  0.0044537  2.756e-03 -0.0006230  0.0067324  0.0058740  0.0019056
# Deso -0.0007804 -1.200e-03  0.0013066 -0.0004435  0.0004592 -0.0031055
# Dica  0.0315015 -1.735e-02  0.0259086  0.0334407  0.0275512  0.0036896
# Drre -0.0042570  5.223e-03 -0.0017084 -0.0037585  0.0004927 -0.0023955
# Elel -0.0443862 -3.976e-02  0.0659183 -0.0197695 -0.1015448  0.0018038
# Eltr  0.1980812  4.235e-02  0.1453907 -0.1230869  0.1106047  0.0003799
# Erci  0.0047155 -9.313e-03  0.0014543 -0.0083823  0.0120805 -0.0153879
# Ersp -0.0709847  7.614e-02  0.0167061 -0.0481423 -0.0093145  0.0263326
# Eunu -0.0037732  3.836e-03 -0.0026022 -0.0023655 -0.0009088 -0.0021630
# Grsq  0.0064540  8.795e-03  0.0033494 -0.0028180  0.0011917  0.0045244
# Gusa  0.0017796  1.552e-03 -0.0009878  0.0050326  0.0018740  0.0062360
# Gypa -0.0062565  5.663e-03  0.0012702  0.0009603  0.0038333  0.0025827
# Heco  0.0274211 -1.543e-01  0.1631251  0.0269522  0.0258884 -0.1122077
# Hehi  0.0131073 -3.809e-03  0.0031431 -0.0177019  0.0122477  0.0062048
# Hisc -0.0055214  5.984e-03  0.0018199 -0.0037493 -0.0007438  0.0026435
# Koma  0.1140091 -1.308e-02  0.1333860 -0.0501406 -0.0644581 -0.0273213
# Laoc -0.0140497  1.289e-02  0.0042528 -0.0027206 -0.0005921  0.0048411
# Lase  0.0143285  1.647e-02 -0.0032989 -0.0026116  0.0080281  0.0127969
# Lede  0.0219477  1.847e-03 -0.0191166 -0.0031671 -0.0350652 -0.0175737
# Lepe -0.0031514  1.938e-03 -0.0032472 -0.0028854  0.0023621 -0.0054719
# Liri  0.0530082  3.674e-03 -0.0493963  0.0386980 -0.0115947 -0.0044204
# Loar -0.0367709  2.121e-02  0.0264622  0.0493476 -0.0343843  0.0027052
# Mahi -0.0433128  4.431e-02  0.0135195 -0.0299829 -0.0055659  0.0216165
# Meof  0.1186287  5.526e-02 -0.0478626  0.0969903  0.0953983  0.0232838
# Mesa -0.0078872 -2.345e-02 -0.0570854 -0.0116781 -0.0150840 -0.0081716
# Mili  0.0021544  9.972e-04  0.0003440 -0.0029574  0.0019334  0.0010176
# Monu -0.0549865  5.392e-02  0.0172490 -0.0370078 -0.0096233  0.0284700
# Mucu -0.0035606  3.583e-03 -0.0039694 -0.0063665  0.0045956 -0.0098584
# Navi  0.2428542  1.792e-01 -0.0225329  0.1157971 -0.0944245  0.0115717
# Oppo  0.1161374  1.759e-02  0.0597526 -0.0581420 -0.0782130  0.0021564
# Pasm  0.0209003 -2.096e-01 -0.1009685  0.0014265 -0.0562755 -0.0191628
# Peal -0.0022382  3.838e-04  0.0002988 -0.0001432 -0.0012742  0.0018098
# Pear  0.0333027  9.449e-03  0.0385098 -0.0197211 -0.0115229 -0.0131869
# Peni  0.0025238  1.912e-03 -0.0013792  0.0040498  0.0042292  0.0016891
# Phho  0.0030957  1.161e-03 -0.0052941  0.0083295 -0.0010074  0.0042709
# Phse  0.0023456  2.685e-03 -0.0010228  0.0005180 -0.0025440  0.0004738
# Plel -0.0237543 -2.156e-02 -0.0049733 -0.0184663  0.0056340  0.0254695
# Plpa -0.0440016  1.552e-02  0.0618102 -0.0479858  0.0079685  0.1802133
# Poar -0.0018027 -1.784e-03 -0.0002646 -0.0017364 -0.0006646  0.0008380
# Poav -0.0331728  3.595e-02  0.0109339 -0.0225261 -0.0044690  0.0158825
# Poco  0.0166166  7.691e-03  0.0026535 -0.0228105  0.0149122  0.0078484
# Pora -0.0130047  1.735e-02  0.0019500 -0.0112878  0.0040490  0.0100229
# Pose  0.0157740  3.169e-02  0.0720463  0.1906255  0.0742368  0.0613076
# Rhar -0.0070132 -6.938e-03  0.0009733 -0.0039911 -0.0099920 -0.0084838
# Satr -0.0045306 -3.898e-05 -0.0009095  0.0085193  0.0024628 -0.0007222
# Save -0.1323259 -3.460e-02  0.0424834  0.0796487  0.0455062 -0.0550823
# Scpa  0.1569186  9.769e-02  0.0834196 -0.0142448 -0.0575349 -0.0116176
# Silo -0.1889523  7.009e-02 -0.0441175 -0.0951262 -0.0283315 -0.0600301
# Somi -0.0015833  2.715e-04  0.0002114 -0.0001013 -0.0009014  0.0012803
# Sosp  0.0007644  3.538e-04  0.0001221 -0.0010494  0.0006860  0.0003611
# Spco  0.0161917 -9.999e-02 -0.0310751  0.0642244 -0.0828786  0.0609806
# Ta   -0.1301761 -5.055e-02  0.0086026  0.1015969 -0.0354144 -0.0322077
# Thar -0.0074151 -3.264e-03 -0.0023830 -0.0055484  0.0026137 -0.0072032
# Trdu  0.1165785  1.625e-02 -0.0319115  0.0531406 -0.1189781 -0.0246137
# unth  0.0011740  1.344e-03 -0.0005119  0.0002592 -0.0012733  0.0002372
# Viam  0.0632107  7.240e-02 -0.0013313  0.0010743 -0.0002311  0.0429232
# Vuoc -0.0278597 -2.018e-02  0.0330798  0.0633414  0.0089125  0.0047970
# 
# 
# Site scores (weighted sums of species scores)
# 
# RDA1      RDA2     RDA3     RDA4     RDA5      RDA6
# row1  -0.51695 -0.191236 -0.11158 -0.15706 -0.26013  0.916797
# row2  -0.32544 -0.232054 -0.10960 -0.26978  0.22115  0.487030
# row3  -0.70776  0.503643 -0.54007 -0.34494  0.08268 -0.856428
# row4  -0.69933  0.352425  0.39356  1.07806  0.35718 -0.016569
# row5  -0.79744  0.848996  0.24035 -0.38350 -0.23304  0.282211
# row6  -0.09155 -0.462302  0.12188 -0.03437 -0.20293 -0.366262
# row7   0.13906 -0.615317  0.21155 -0.22847  0.37766 -1.013284
# row8  -0.17123 -0.700303  0.19584  0.25868 -0.01182  0.332969
# row9  -0.16514 -0.399237  0.16482 -0.14501 -0.79775 -0.006202
# row10  0.33592 -0.008838 -0.90700  0.39797 -0.82448 -0.325053
# row11  0.53222  0.448152  0.83096 -0.16277 -0.25519 -0.146759
# row12  0.36176  0.149018 -0.43307  0.80838  0.44394 -0.013135
# row13  0.26404 -0.016235  0.63223  0.39104 -0.17957 -0.200551
# row14  0.25358 -0.508041 -0.50682 -0.12508  0.26722  0.173080
# row15  0.33776  0.052066  0.16198 -0.56698  0.24076  0.183650
# row16  0.44193  0.460978 -0.02863 -0.11541 -0.34078  0.132488
# row17  0.41819  0.541892 -0.39076  0.03633  0.27732  0.250040
# row18  0.39037 -0.223608  0.07435 -0.43707  0.83778  0.185978
# 
# 
# Site constraints (linear combinations of constraining variables)
# 
# RDA1     RDA2     RDA3     RDA4     RDA5      RDA6
# row1  -0.414635  0.07109  0.05536 -0.02653 -0.23606  0.335275
# row2  -0.369227 -0.36548 -0.09877 -0.41717  0.02868  0.432985
# row3  -0.465170  0.28602 -0.47931 -0.42591  0.34867 -0.807703
# row4  -0.679558  0.23220  0.14238  1.01080  0.18923 -0.081712
# row5  -0.490201  0.53126  0.16157 -0.33287 -0.06604  0.234699
# row6  -0.179116 -0.17718  0.02486 -0.10193 -0.25519 -0.216676
# row7  -0.232660 -0.35777  0.38957 -0.13222  0.13692 -0.925898
# row8  -0.001852 -0.98947  0.25118  0.36921  0.21595  0.326766
# row9  -0.347318 -0.23565  0.18775 -0.26318 -0.73173  0.449210
# row10  0.261745 -0.06196 -0.86209  0.43760 -0.78352 -0.342347
# row11  0.680266  0.19301  0.78663 -0.40284 -0.23538 -0.269366
# row12  0.318281  0.24114 -0.17394  0.51074  0.53336  0.213019
# row13  0.424224  0.02783  0.62819  0.48395 -0.04386 -0.226698
# row14  0.490209 -0.61147 -0.62267 -0.23109  0.19955  0.009301
# row15  0.217304  0.10058  0.03470 -0.29830  0.19501  0.102638
# row16  0.426797  0.48849 -0.18611  0.09425 -0.46290  0.086219
# row17  0.179820  0.54919 -0.08638 -0.03466  0.18280  0.435316
# row18  0.181091  0.07816 -0.15291 -0.23983  0.78452  0.244972
# 
# 
# Biplot scores for constraining variables
# 
# RDA1     RDA2     RDA3     RDA4      RDA5      RDA6
# Soil.Water -0.15895  0.12212 -0.10192 -0.78020 -0.089523  0.222853
# Sand       -0.11797  0.05516 -0.71734  0.11047 -0.171528  0.009569
# Clay        0.03056  0.13557  0.65635 -0.25070  0.007679 -0.075348
# BD          0.11871  0.56879  0.13816  0.74388 -0.092174 -0.222743
# Slope       0.46114 -0.05615 -0.62734 -0.04008  0.175404 -0.119593
# Elev        0.86874  0.30514  0.03160  0.10736 -0.011033 -0.296267
# Aspect      0.19481  0.08033 -0.17268  0.24851 -0.551621  0.183638
# Precip      0.89504  0.33204 -0.04332  0.15898  0.092733  0.049884
# DW         -0.77848 -0.31971  0.27499 -0.22577 -0.124285  0.241172


#plot rda 1 vs rda 2 with species labels
library(ggord)

ggord(tbRDA, env$location, poly = FALSE, ptslab = TRUE, size = 3,veclsz = 0.4, arrow = 0.3, addcol = "grey10", grp_title = "Site", repel = TRUE, alpha = 2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_linetype_manual(values = c('solid', 'solid'))+
  labs(color="", linetype = "", shape = "")+
  scale_colour_manual(values=c("grey60", "grey80"), labels = c("BLM", "Sun Prairie"), name = "")+
  scale_linetype_manual(values = c("twodash", "solid"), labels = c("BLM", "Sun Prairie"), name = "")+
  scale_shape_manual(values = c(17, 16),labels = c("BLM", "Sun Prairie"), name = "")
  #coord_cartesian(xlim=c(-1.1,1.2), ylim=c(-1.1,1.2))


##final graph
ggord(tbRDA, env$location, poly = FALSE, ptslab = TRUE, size = 4,veclsz = 0.4, arrow = 0.3, addcol = "grey10",addsize = 2.5, grp_title = "Site", repel = TRUE, alpha =2, ext = 1.3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_linetype_manual(values = c('solid', 'solid'))+
  labs(color="", linetype = "", shape = "")+
  scale_colour_manual(values=c("grey60", "grey80"), labels = c("BLM", "Sun Prairie"), name = "")+
  scale_shape_manual(values = c(17, 16),labels = c("BLM", "Sun Prairie"), name = "")+
  theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16, color = "black"),
               axis.title.y=element_text(size=20, angle=90, vjust=0.7), axis.text.y=element_text(size=16, color= "black"),
               panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
               legend.title=element_blank(), legend.text=element_text(size=18), panel.border=element_rect(color="black", fill = NA, size = 1))

#poly=FALSE no polygon fill
#ptslab=TRUE add species names to points
#size env variables point size
#repel=TRUE puts labels inside plot
#ext distance from env variable label to arrow
#addsize = the size of species labels
#addcol = color of species labels
# grp_in = NULL,
# cols = NULL,
# facet = FALSE,
# nfac = NULL,
# addpts = NULL,
# obslab = FALSE,
# ptslab = FALSE,
# ellipse = TRUE,
# ellipse_pro = 0.95,
# poly = TRUE,
# polylntyp = "solid",
# hull = FALSE,
# arrow = 0.4,
# labcol = "black",
# veccol = "black",
# vectyp = "solid",
# veclsz = 0.5,
# ext = 1.2,
# repel = FALSE,
# vec_ext = 1,
# vec_lab = NULL,
# size = 4,
# sizelab = NULL,
# addsize = size/2,
# addcol = "blue",
# addpch = 19,
# txt = 4,
# alpha = 1,
# alpha_el = 0.4,
# xlims = NULL,
# ylims = NULL,
# var_sub = NULL,
# coord_fix = TRUE,
# parse = TRUE,
# grp_title = "Groups",
# force = 1,
# max.overlaps = 10,
# exp = c(0, 0)


#check how much each rda/pca explains
constrained_eig <- tbRDA$CCA$eig/tbRDA$tot.chi*100
unconstrained_eig <- tbRDA$CA$eig/tbRDA$tot.chi*100
expl_var <- c(constrained_eig, unconstrained_eig)
barplot (expl_var[1:20], col = c(rep ('red', length (constrained_eig)), rep ('black', length (unconstrained_eig))),
         las = 2, ylab = '% variation')


ordiplot (tbRDA)

ordiplot (tbRDA, display = 'species', choices = c(1,2), type = 'n')
orditorp (tbRDA, display = 'species', choices = c(1,2), pcol = 'grey', pch = '+')








#standardize values for envt conditions, make into z scores
#RDA or CCA, layers on plant community, which environmental variables are driving different plant communities
#change previous graphs to BLM, Sun Prairie
#describing plant communities in two sites - apply to other northern mixed grass prairies, conservation/restoration using certain species
