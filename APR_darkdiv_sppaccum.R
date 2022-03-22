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
  scale_fill_manual(values=c("grey40", "grey"), labels=c("Non-Native", "Native"))+
  theme(axis.text.x=element_text(size=24, color = "black")) +
  scale_x_discrete(labels = c("Bison", "Cattle"))+
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

specaccummean <- ggplot(data=barGraphStats(data=specaccum, variable="richness", byFactorNames=c("plot_size", "management")), aes(x=plot_size, y=mean, fill = management, colour = management))+
  geom_line(aes(linetype = management), size = 1.5) +
  geom_point(size = 2, color = "black")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9), color = "black")+
  ylab('Average Species Richness\n')+ xlab(bquote('Plot Size'~(m^2)))+
  theme(axis.title=element_text(size=24, color = "black"),legend.position=c(.2,.9))+
  scale_colour_manual(values=c("grey40", "grey"))+
  expand_limits(y= 35)

theme(legend.position=c(.87,.9), axis.text.x=element_text(size=24, color = "black"))

                                         

###rank abundance curves
rankAbundance <- relCover%>%
  group_by(location, provenance, growth_form, genus_species)%>%
  summarize(avg_cover=mean(rel_cover))%>%
  ungroup()%>%
  filter(genus_species!='unknown:small_thin_leaf')%>%
  mutate(provenance=ifelse(genus_species=='Penstemon_albicula', 'native', as.character(provenance)))%>%
  arrange(location, -avg_cover)%>%
  group_by(location)%>%
  mutate(rank=seq_along(location))%>%
  ungroup()

BLMrank <- ggplot(data=subset(rankAbundance, location=='BLM', avg_cover>0), aes(x=rank, y=avg_cover)) +
  geom_line() +
  geom_point(aes(colour=provenance, shape=growth_form), size=3) +
  xlab('Species Rank') +
  ylab('Relative Percent Cover') +
  # scale_x_continuous(expand=c(0,0), limits=c(0.5,17), breaks=seq(0,17,5)) +
  # scale_y_continuous(expand=c(0,0), limits=c(0,60), breaks=seq(0,60,10)) +
  geom_text(aes(y=avg_cover+0.5, x=rank+0.1, label=genus_species), hjust='left', vjust='bottom', angle=25, size=5)

BDrank <- ggplot(data=subset(rankAbundance, location=='BD', avg_cover>0), aes(x=rank, y=avg_cover)) +
  geom_line() +
  geom_point(aes(colour=provenance, shape=growth_form), size=3) +
  xlab('Species Rank') +
  ylab('Relative Percent Cover') +
  # scale_x_continuous(expand=c(0,0), limits=c(0.5,17), breaks=seq(0,17,5)) +
  # scale_y_continuous(expand=c(0,0), limits=c(0,60), breaks=seq(0,60,10)) +
  geom_text(aes(y=avg_cover+0.5, x=rank+0.1, label=genus_species), hjust='left', vjust='bottom', angle=25, size=5)

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


