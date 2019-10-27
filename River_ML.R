#### import libraries ####
# data wrangling

library(tidyverse)
library(reshape)
library(data.table)
library(xlsx)
library(gridExtra)
library(grid)
library(chron)
library(devtools)
library(SOfun)
library(usethis)

# data visualization 

library(GGally)
library(RColorBrewer)
library(proj4)
library(leaflet)
library(leaflet.minicharts)
library(RColorBrewer)
library(mapview)
library(htmlwidgets)
library(corrplot)
library(mice)
library(VIM)
library(ggmosaic)
library(esquisse)

# data analysis - machine learning - statistics

library(Rcpp)
library(vegan)
library(cluster)
library(MuMIn) # R2c and R2m
library(nlme)
library(gamm4)
library(lme4)
library(dunn.test)  # Kruskal Wallis test
library(car)
library(psych)
library(psycho)
library(remotes)
library(mlr)
library(mlrMBO)
library(DiceKriging)
library(rgenoud)
library(randomForest)
library(Metrics)
library(Hmisc)
library(xgboost)
library(checkmate)
library(ranger)
library(parallel)
library(parallelMap)


#### import data ####

river <- read_csv("River.csv")
# change some character columns to factor

river$Date <- as.factor(river$Date)
river$Location <- as.factor(river$Location)
river$River <- as.factor(river$River)
river$LB <- as.factor(river$LB)
river$RB <- as.factor(river$RB)
river$Shading <- as.factor(river$Shading)
river$Erosion <- as.factor(river$Erosion)
river$Flow_variation <- as.factor(river$Flow_variation)

# change the colnames

colnames(river)[6:38] <- c("Water temperature", "pH", "DO", "EC", "TDS", "Salinity", "Turbility", "Chlorophyll", "Left Bank", 
                           "Right Bank", "Shading", "Erosion", "Flow variability", "Average depth", "Average velocity",
                           "Pool Class", "BOD", "COD", "TN" ,"NH4", "NO2", "NO3", "TP", "PO4", "Air temperature", "Wind velocity",
                           "Rain", "Solar radiation", "Latitude", "Longitude", "Dissolved N2O", "Dissolved CH4", "Dissolved CO2"
                           )
#### Correlation coefficients #### 
variable_river <- cbind(river[,6:13], river[,22:33]) %>% select(-Rain)

corr_river <- cor(variable_river, use = 'pairwise')
p.mat <- cor.mtest(variable_river)$p
colnames(p.mat) <- colnames(corr_river)
row.names(p.mat) <- colnames(corr_river)

# GGally Not really nice
ggsave("Corr_coeff.tiff", ggpairs(variable_river,
                                  lower = list(continuous = wrap("smooth", color = "deepskyblue")),
                                  upper = list(continuous = wrap("cor", size = 3, color = "tomato"))
) + theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()),
units = 'cm', height = 50, width = 50, dpi = 300)


# corrplot nicer but cannot handle the categorical variables

tiff("corr_coeff_2.tiff",units = 'cm',height = 50,width = 50,res = 300, pointsize = 12)
corrplot(corr_river, p.mat = p.mat, method = "circle", type = "upper",
         sig.level = 0.05, insig = "blank", order = "alphabet")
dev.off()

# Using mosaic plot to represent the relationship among two or more categorical variables 
# in this case, only for river and hydromorphological data

ggsave("Mosaic_river_LB.tiff", ggplot(river)+
           geom_mosaic(aes(x= product(River), fill = "Lelf Bank"))+
           labs(x ="", y = ""),
       units = 'cm', height = 30, width = 40, dpi = 300)
ggsave("Mosaic_river_RB.tiff", ggplot(river)+
           geom_mosaic(aes(x= product(River), fill = "Right Bank"))+
           labs(x ="", y = ""),
       units = 'cm', height = 30, width = 40, dpi = 300)
ggsave("Mosaic_river_FV.tiff",ggplot(river)+
           geom_mosaic(aes(x= product(River), fill = "Flow variability"))+
           labs(x ="", y = ""),
       units = 'cm', height = 30, width = 40, dpi = 300)
ggsave("Mosaic_river_shading.tiff",ggplot(river)+
           geom_mosaic(aes(x= product(River), fill = Shading))+
           labs(x ="", y = ""),
       units = 'cm', height = 30, width = 40, dpi = 300)
ggsave("Mosaic_river_erosion.tiff",ggplot(river)+
           geom_mosaic(aes(x= product(River), fill = Erosion))+
           labs(x ="", y = ""),
       units = 'cm', height = 30, width = 40, dpi = 300)
ggsave("Mosaic_river_pool.tiff",ggplot(river)+
           geom_mosaic(aes(x= product(River), fill = "Pool Class"))+
           labs(x ="", y = ""),
       units = 'cm', height = 30, width = 40, dpi = 300)


#### Spatio-temporal variability ####
#*** Mixed model ####

# log transform and data standardization 

river_sta <- river
river_sta$log_N2O <- log(river_sta$Dis_N2O_cor)
river_sta$sta_N2O <- standardize(river_sta$log_N2O) 
river_sta$log_CH4 <- log(river_sta$Dis_CH4_cor)
river_sta$sta_CH4 <- standardize(river_sta$log_CH4) 
river_sta$log_CO2 <- log(river_sta$Dis_CO2_cor)
river_sta$sta_CO2 <- standardize(river_sta$log_CO2) 

#***** N2O #####
# diagnostic outliers 

# using box plot 

boxplot.stats(river_sta$sta_N2O)$out
# both result in two outliers whose values are 3.939392 and 2.722888 --> remove them or not --> yes!!!

river_sta_N2O <- river_sta %>% filter (sta_N2O < 2.7)


# using cleveland dotplot

ggsave("Cleveland_N2O.tiff", ggplot(river_sta_N2O) +
    aes(x = sta_N2O, y = No) +
    geom_point(size = 3L, colour = "#0c4c8a") +
    xlab(bquote("Standardized Dissolved "*N[2]*"O")) +
    ylab("Order of the data")+
    theme_bw()+
    theme(axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          text = element_text(size = 14)),
    units = 'cm', height = 20, width = 20, dpi = 300)

# mixed model

set.seed(1)
river_lmm_N2O <- lmer(sta_N2O~1 + (1|River/Date), data = river_sta_N2O)
r.squaredGLMM(river_lmm_N2O)
summary(river_lmm_N2O)
vcov(river_lmm_N2O)

# using built-in function
river_lmm_N2O_fort <- fortify(river_lmm_N2O)
diagPlot<-function(model){
    p1<-ggplot(model, aes(.fitted, .resid))+geom_point()
    p1<-p1+stat_smooth(method="loess")+geom_hline(yintercept=0, col="red", linetype="dashed")
    p1<-p1+xlab("Fitted values")+ylab("Residuals")
    p1<-p1+ggtitle("Residual vs Fitted Plot")+theme_bw()
    
    p2<-ggplot(model, aes(qqnorm(.scresid)[[1]], .scresid))+geom_point(na.rm = TRUE)
    p2<-p2+geom_line(aes(qqline(.scresid)))+xlab("Theoretical Quantiles")+ylab("Standardized Residuals")
    p2<-p2+ggtitle("Normal Q-Q")+theme_bw()
    
    p3<-ggplot(model, aes(.fitted, sqrt(abs(.scresid))))+geom_point(na.rm=TRUE)
    p3<-p3+stat_smooth(method="loess", na.rm = TRUE)+xlab("Fitted Value")
    p3<-p3+ylab(expression(sqrt("|Standardized residuals|")))
    p3<-p3+ggtitle("Scale-Location")+theme_bw()
    
    return(list(rvfPlot=p1, 
                qqPlot=p2, 
                sclLocPlot=p3))
}

diagPlts <- diagPlot(river_lmm_N2O)
ggsave("Diagnostic_plot1.tiff",diagPlts[[1]],
       units = 'cm', height = 20, width = 20, dpi = 300)
ggsave("Diagnostic_plot2.tiff",diagPlts[[2]],
       units = 'cm', height = 20, width = 20, dpi = 300)
ggsave("Diagnostic_plot3.tiff",diagPlts[[3]],
       units = 'cm', height = 20, width = 20, dpi = 300)

# ICC values !!! river should be date and river
ICC_river_N2O <- as.numeric(VarCorr(river_lmm2_N2O)[[2]])/
    (as.numeric(VarCorr(river_lmm2_N2O)[[2]]) + as.numeric(VarCorr(river_lmm2_N2O)[[4]])+
         as.numeric(VarCorr(river_lmm2_N2O)[[5]]))
ICC_date_N2O <- (as.numeric(VarCorr(river_lmm2_N2O)[[2]])+ as.numeric(VarCorr(river_lmm2_N2O)[[4]]))/
    (as.numeric(VarCorr(river_lmm2_N2O)[[2]]) + as.numeric(VarCorr(river_lmm2_N2O)[[4]])+
         as.numeric(VarCorr(river_lmm2_N2O)[[5]]))

# ICC river and date are low meaning low spatiotemporal variability

# As such, can applied Kruskal-Wallis

#***** CH4 #####
# diagnostic outliers 

# using box plot 

boxplot.stats(river_sta$sta_CH4)$out

# using cleveland dotplot

ggsave("Cleveland_CH4.tiff",ggplot(river_sta) +
    aes(x = sta_CH4, y = No) +
    geom_point(size = 3L, colour = "#0c4c8a") +
    xlab(bquote("Standardized Dissolved "*CH[4]*"")) +
    ylab("Order of the data")+
    theme_bw()+
    theme(axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          text = element_text(size = 14)),
    units = 'cm', height = 20, width = 20, dpi = 300)


# mixed model

set.seed(1)

river_lmm_CH4 <- lmer(sta_CH4~1 + (1|River/Date), data = river_sta)
r.squaredGLMM(river_lmm_CH4)
summary(river_lmm_CH4)
vcov(river_lmm_CH4)

# using built-in function
river_lmm_CH4_fort <- fortify(river_lmm_CH4)


diagPlts <- diagPlot(river_lmm_CH4)
ggsave("Diagnostic_plot1.tiff",diagPlts[[1]],
       units = 'cm', height = 20, width = 20, dpi = 300)
ggsave("Diagnostic_plot2.tiff",diagPlts[[2]],
       units = 'cm', height = 20, width = 20, dpi = 300)
ggsave("Diagnostic_plot3.tiff",diagPlts[[3]],
       units = 'cm', height = 20, width = 20, dpi = 300)

# ICC values 
ICC_river_CH4 <- 0.2713/ (0.2713 + 0.2504 + 0.5999)
ICC_date_CH4 <- (0.2713 + 0.2504)/ (0.2713 + 0.2504 + 0.5999)

# ICC river and date are low meaning low spatiotemporal variability

# As such, can applied Kruskal-Wallis


#***** CO2 #####
# diagnostic outliers 

# using box plot 

boxplot.stats(river_sta$sta_CO2)$out

# using cleveland dotplot

ggsave("Cleveland_CO2.tiff",ggplot(river_sta) +
           aes(x = sta_CO2, y = No) +
           geom_point(size = 3L, colour = "#0c4c8a") +
           xlab(bquote("Standardized Dissolved "*CH[4]*"")) +
           ylab("Order of the data")+
           theme_bw()+
           theme(axis.title.y = element_text(size = 14),
                 axis.title.x = element_text(size = 14),
                 text = element_text(size = 14)),
       units = 'cm', height = 20, width = 20, dpi = 300)


# mixed model

set.seed(1)

river_lmm_CO2 <- lmer(sta_CO2~1 + (1|River/Date), data = river_sta)
r.squaredGLMM(river_lmm_CO2)
summary(river_lmm_CO2)
vcov(river_lmm_CO2)

# using built-in function
river_lmm_CO2_fort <- fortify(river_lmm_CO2)


diagPlts <- diagPlot(river_lmm_CO2)
ggsave("Diagnostic_plot1.tiff",diagPlts[[1]],
       units = 'cm', height = 20, width = 20, dpi = 300)
ggsave("Diagnostic_plot2.tiff",diagPlts[[2]],
       units = 'cm', height = 20, width = 20, dpi = 300)
ggsave("Diagnostic_plot3.tiff",diagPlts[[3]],
       units = 'cm', height = 20, width = 20, dpi = 300)

# ICC values 
ICC_river_CO2 <- 0.2155/ (0.2155 + 0.2635  + 0.6616 )
ICC_date_CO2 <- (0.2155 + 0.2635 )/ (0.2155 + 0.2635  + 0.6616 )

# ICC river and date are low meaning low spatiotemporal variability

# As such, can applied Kruskal-Wallis


#*** Kruskal-wallis with Bonferroni correction ####

river_dun_CO2 <- as.data.frame(dunn.test(river$Dis_CO2_cor, river$River, method=c('bonferroni'))[4:5])
river_dun_CH4 <- as.data.frame(dunn.test(river$Dis_CH4_cor, river$River, method=c('bonferroni'))[4:5])
river_dun_N2O <- as.data.frame(dunn.test(river$Dis_N2O_cor, river$River, method=c('bonferroni'))[4:5])
river_dun <- bind_rows(river_dun_CO2, river_dun_CH4, river_dun_N2O, .id = "Dissolved gases")
write_csv(river_dun, "river_dun.csv")

river_dun_CO2_D <- as.data.frame(dunn.test(river$Dis_CO2_cor, river$Date, method=c('bonferroni'))[4:5])
river_dun_CH4_D <- as.data.frame(dunn.test(river$Dis_CH4_cor, river$Date, method=c('bonferroni'))[4:5])
river_dun_N2O_D <- as.data.frame(dunn.test(river$Dis_N2O_cor, river$Date, method=c('bonferroni'))[4:5])
river_dun_D <- bind_rows(river_dun_CO2_D, river_dun_CH4_D, river_dun_N2O_D, .id = "Dissolved gases")
write_csv(river_dun_D, "river_dun_D.csv")

#*** Friedmann test ####

# Different river, block date 




# Different date, block river



#### RF - TRY ####

# choose the correct dataset. Remove correlated parameters.

river_RF <- variable_river %>% select(-c(Salinity, EC, TDS, BOD, TN, TP))
river_RF$N2O <- river_sta$sta_N2O
river_RF$CH4 <- river_sta$sta_CH4
river_RF$CO2 <- river_sta$sta_CO2

# make a random model for N2O

set.seed(1234)
regressor_N2O <- randomForest(river_RF[,-c(14:16)], river_RF$N2O, ntree = 500, importance = TRUE)

imp_N2O <- randomForest::importance(regressor_N2O, type=1)
featureImportance_N2O <- data.frame(Feature=row.names(imp_N2O), Importance=imp_N2O[,1])
ggsave("RF_N2O_try.tiff",ggplot(featureImportance_N2O[1:10,], aes(x=reorder(Feature, Importance), y=Importance)) +
    geom_bar(stat="identity", fill="tomato") +
    coord_flip() + 
    theme_bw(base_size=20) +
    xlab("") +
    ylab("Importance") + 
    ggtitle("Random Forest Feature Importance") +
    theme(plot.title=element_text(size=18)),
    units = 'cm', height = 20, width = 20, dpi = 300)


# make a random model for CH4

set.seed(1234)
regressor_CH4 <- randomForest(river_RF[,-c(14:16)], river_RF$CH4, ntree = 500, importance = TRUE)

imp_CH4 <- randomForest::importance(regressor_CH4, type=1)
featureImportance_CH4 <- data.frame(Feature=row.names(imp_CH4), Importance=imp_CH4[,1])
ggsave("RF_CH4_try.tiff",ggplot(featureImportance_CH4[1:10,], aes(x=reorder(Feature, Importance), y=Importance)) +
    geom_bar(stat="identity", fill="tomato") +
    coord_flip() + 
    theme_bw(base_size=20) +
    xlab("") +
    ylab("Importance") + 
    ggtitle("Random Forest Feature Importance") +
    theme(plot.title=element_text(size=18)),
    units = 'cm', height = 20, width = 20, dpi = 300)


# make a random model for CO2

set.seed(1234)
regressor_CO2 <- randomForest(river_RF[,-c(14:16)], river_RF$CO2, ntree = 1000, importance = TRUE)

imp_CO2 <- randomForest::importance(regressor_CO2, type=1)
featureImportance_CO2 <- data.frame(Feature=row.names(imp_CO2), Importance=imp_CO2[,1])
ggsave("RF_CO2_try.tiff",ggplot(featureImportance_CO2[1:10,], aes(x=reorder(Feature, Importance), y=Importance)) +
    geom_bar(stat="identity", fill="tomato") +
    coord_flip() + 
    theme_bw(base_size=20) +
    xlab("") +
    ylab("Importance") + 
    ggtitle("Random Forest Feature Importance") +
    theme(plot.title=element_text(size=18)),
units = 'cm', height = 20, width = 20, dpi = 300)
#### RF - TRY 2 #### Try to do CV on this # START FROM HERE 22:43 20/08/19

# choose the correct dataset. Remove correlated parameters. But you take into account the categorical data 

river_RF_2 <- river[,-c(1:5)] %>% select(-c(Salinity, EC, TDS, BOD, TN, TP, Rain, Latitude, Longitude, 
                                            `Dissolved N2O`, `Dissolved CH4`, `Dissolved CO2`))
river_RF_2$N2O <- river_sta$sta_N2O
river_RF_2$CH4 <- river_sta$sta_CH4
river_RF_2$CO2 <- river_sta$sta_CO2

# make a random model for N2O

set.seed(1234)
regressor_N2O_2 <- randomForest(river_RF_2[,-c(22:24)], river_RF_2$N2O, ntree = 500, importance = TRUE)

imp_N2O_2 <- randomForest::importance(regressor_N2O_2, type=1)
featureImportance_N2O_2 <- data.frame(Feature=row.names(imp_N2O_2), Importance=imp_N2O_2[,1])
ggsave("RF_2_N2O_try_2.tiff",ggplot(featureImportance_N2O_2[,], aes(x=reorder(Feature, Importance), y=Importance)) +
           geom_bar(stat="identity", fill="tomato") +
           coord_flip() + 
           theme_bw(base_size=20) +
           xlab("") +
           ylab("Importance") + 
           ggtitle("Random Forest Feature Importance") +
           theme(plot.title=element_text(size=18)),
       units = 'cm', height = 20, width = 20, dpi = 300)


# make a random model for CH4

set.seed(1234)
regressor_CH4_2 <- randomForest(river_RF_2[,-c(22:24)], river_RF_2$CH4, ntree = 500, importance = TRUE)

imp_CH4_2 <- randomForest::importance(regressor_CH4_2, type=1)
featureImportance_CH4_2 <- data.frame(Feature=row.names(imp_CH4_2), Importance=imp_CH4_2[,1])
ggsave("RF_2_CH4_try_2.tiff",ggplot(featureImportance_CH4_2[,], aes(x=reorder(Feature, Importance), y=Importance)) +
           geom_bar(stat="identity", fill="tomato") +
           coord_flip() + 
           theme_bw(base_size=20) +
           xlab("") +
           ylab("Importance") + 
           ggtitle("Random Forest Feature Importance") +
           theme(plot.title=element_text(size=18)),
       units = 'cm', height = 20, width = 20, dpi = 300)


# make a random model for CO2

set.seed(1234)
regressor_CO2_2 <- randomForest(river_RF_2[,-c(22:24)], river_RF_2$CO2, ntree = 1000, importance = TRUE)

imp_CO2_2 <- randomForest::importance(regressor_CO2_2, type=1)
featureImportance_CO2_2 <- data.frame(Feature=row.names(imp_CO2_2), Importance=imp_CO2_2[,1])
ggsave("RF_2_CO2_try_2.tiff",ggplot(featureImportance_CO2_2[,], aes(x=reorder(Feature, Importance), y=Importance)) +
           geom_bar(stat="identity", fill="tomato") +
           coord_flip() + 
           theme_bw(base_size=20) +
           xlab("") +
           ylab("Importance") + 
           ggtitle("Random Forest Feature Importance") +
           theme(plot.title=element_text(size=18)),
       units = 'cm', height = 20, width = 20, dpi = 300)
#### RF - MLR ####

colnames(river_RF)[1] <- "Water_temperature"
colnames(river_RF)[11] <- "Air_temperature"
colnames(river_RF)[12] <- "Wind_velocity"
colnames(river_RF)[13] <- "Solar_radiation"

#*** N2O ####
# creat mlr task and convert factors to dummy features

trainTask_N2O = createDummyFeatures(river_RF[,-c(15:16)],target = "N2O")
trainTask_N2O <- makeRegrTask(data = river_RF[,-c(15:16)],target = "N2O")

# create mlr learner
set.seed(1234)
lrn_N2O <- makeLearner("regr.ranger")
# lrn_N2O$par.vals <- list(ntree = 100L, importance = TRUE)
cv_N2O <- makeResampleDesc(method = "LOO")
# set parallel backend
parallelStartSocket(cpus = detectCores()-1)
res_N2O <- resample(learner = lrn_N2O, task = trainTask_N2O, resampling = cv_N2O)
# Tuning hyperparameters

# Parameter Tuning: Mainly, there are three parameters in the random forest algorithm which you should look at (for tuning):
# ntree - The number of trees to grow. Larger the tree, it will be more computationally expensive to build models.
# mtry - It refers to how many variables we should select at a node split. 
# Also as mentioned above, the default value is p/3 for regression and sqrt(p) for classification. 
# We should always try to avoid using smaller values of mtry to avoid overfitting.
# nodesize - It refers to how many observations we want in the terminal nodes.
# This parameter is directly related to tree depth. Higher the number, lower the tree depth. 
# With lower tree depth, the tree might even fail to recognize useful signals from the data.

# To know which hyperparameter can be tuned using getParamSet
getParamSet(lrn_N2O)

# Only tune three abovementioned hyperparameters

params_N2O <- makeParamSet(makeIntegerParam("mtry", lower = 2, upper = 10),
                           makeIntegerParam("min.node.size", lower = 2, upper = 25))
tc_N2O <- makeTuneControlMBO(budget = 100)
tr_N2O = tuneParams(learner = lrn_N2O, task = trainTask_N2O, resampling = cv_N2O,
                par.set = params_N2O, control = tc_N2O, show.info = T)

parallelStop()
tr_N2O
# Tune result:
# Op. pars: mtry=10; min.node.size=2
# mse.test.mean=0.3917445

# Apply the optimal RF

set.seed(1234)
regressor_N2O_opt <- randomForest(river_RF[,-c(14:16)], river_RF$N2O, ntree = 500, importance = TRUE, mtry = 10, nodesize = 2)

imp_N2O_opt <- randomForest::importance(regressor_N2O_opt, type=1)
featureImportance_N2O_opt <- data.frame(Feature=row.names(imp_N2O_opt), Importance=imp_N2O_opt[,1])

labels_N2O <- c("NH[4]^+{}", "Water ~~ temperature", "PO[4]^3^-{}", "Turbility", "DO", "NO[2]^-{}", "NO[3]^-{}", 
                "COD", "Wind ~ velocity", "Solar ~ radiation", "pH", "Chlorophyll *~alpha", "Air ~ temperature")
labels_N2O <- rev(labels_N2O)
labels_N2O_parse <- parse(text = labels_N2O)

ggsave("RF_N2O_opt.tiff", ggplot(featureImportance_N2O_opt[1:13,], aes(x=reorder(Feature, Importance), y=Importance)) +
           geom_bar(stat="identity", fill="tomato") +
           coord_flip() + 
           theme_bw(base_size=20) +
           scale_x_discrete(labels = labels_N2O_parse) +
           labs(x = NULL, y = "Importance") + 
           ggtitle("Random Forest Feature Importance") +
           theme(plot.title=element_text(size=18)),
       units = 'cm', height = 20, width = 20, dpi = 300)


#*** CH4 ####
# creat mlr task and convert factors to dummy features

trainTask_CH4 = createDummyFeatures(river_RF %>% select(-N2O,-CO2),target = "CH4")
trainTask_CH4 <- makeRegrTask(data = river_RF %>% select(-N2O,-CO2),target = "CH4")

# create mlr learner
set.seed(1234)
lrn_CH4 <- makeLearner("regr.ranger")
# lrn_CH4$par.vals <- list(ntree = 100L, importance = TRUE)
cv_CH4 <- makeResampleDesc(method = "LOO")
# set parallel backend
parallelStartSocket(cpus = detectCores()-1)
res_CH4 <- resample(learner = lrn_CH4, task = trainTask_CH4, resampling = cv_CH4)
# Tuning hyperparameters

# Only tune three abovementioned hyperparameters

params_CH4 <- makeParamSet(makeIntegerParam("mtry", lower = 2, upper = 10),
                           makeIntegerParam("min.node.size", lower = 2, upper = 25))
tc_CH4 <- makeTuneControlMBO(budget = 100)
tr_CH4 = tuneParams(learner = lrn_CH4, task = trainTask_CH4, resampling = cv_CH4,
                par.set = params_CH4, control = tc_CH4, show.info = F)

parallelStop()
tr_CH4
# Tune result:
# mtry=5; min.node.size=2
# mse.test.mean=0.3260113

# Apply the optimal RF

set.seed(1234)
regressor_CH4_opt <- randomForest(river_RF[,-c(14:16)], river_RF$CH4, ntree = 500, importance = TRUE, mtry = 5, nodesize = 2)

imp_CH4_opt <- randomForest::importance(regressor_CH4_opt, type=1)
featureImportance_CH4_opt <- data.frame(Feature=row.names(imp_CH4_opt), Importance=imp_CH4_opt[,1])

labels_CH4 <- c("DO", "Turbility",  "NO[2]^-{}", "NO[3]^-{}",  "PO[4]^3^-{}", "Water ~~ temperature", "COD" , "Air ~ temperature",
                "NH[4]^+{}","pH",  "Wind ~ velocity", "Solar ~ radiation", "Chlorophyll*~alpha")
labels_CH4 <- rev(labels_CH4)
labels_CH4_parse <- parse(text = labels_CH4)

ggsave("RF_CH4_opt.tiff",ggplot(featureImportance_CH4_opt[1:13,], aes(x=reorder(Feature, Importance), y=Importance)) +
           geom_bar(stat="identity", fill="tomato") +
           coord_flip() + 
           theme_bw(base_size=20) +
           scale_x_discrete(labels = labels_CH4_parse) +
           labs(x = NULL, y = "Importance") + 
           ggtitle("Random Forest Feature Importance") +
           theme(plot.title=element_text(size=18)),
       units = 'cm', height = 20, width = 20, dpi = 300)


#*** CO2 ####
# creat mlr task and convert factors to dummy features

trainTask_CO2 = createDummyFeatures(river_RF %>% select(-N2O,-CH4),target = "CO2")
trainTask_CO2 <- makeRegrTask(data = river_RF %>% select(-N2O,-CH4),target = "CO2")

# create mlr learner
set.seed(1234)
lrn_CO2 <- makeLearner("regr.ranger")
# lrn_CO2$par.vals <- list(ntree = 100L, importance = TRUE)
cv_CO2 <- makeResampleDesc(method = "LOO")
# set parallel backend
parallelStartSocket(cpus = detectCores()-1)
res_CO2 <- resample(learner = lrn_CO2, task = trainTask_CO2, resampling = cv_CO2)
# Tuning hyperparameters

# Only tune three abovementioned hyperparameters

params_CO2 <- makeParamSet(makeIntegerParam("mtry", lower = 2, upper = 10),
                           makeIntegerParam("min.node.size", lower = 2, upper = 25))
tc_CO2 <- makeTuneControlMBO(budget = 100)
tr_CO2 = tuneParams(learner = lrn_CO2, task = trainTask_CO2, resampling = cv_CO2,
                    par.set = params_CO2, control = tc_CO2, show.info = F)

parallelStop()
tr_CO2
# Tune result:
# mtry=9; min.node.size=2
# mse.test.mean=0.2113096

# Apply the optimal RF

set.seed(1234)
regressor_CO2_opt <- randomForest(river_RF[,-c(14:16)], river_RF$CO2, ntree = 500, importance = TRUE, mtry = 9, nodesize = 2)

imp_CO2_opt <- randomForest::importance(regressor_CO2_opt, type=1)
featureImportance_CO2_opt <- data.frame(Feature=row.names(imp_CO2_opt), Importance=imp_CO2_opt[,1])

labels_CO2 <- c("DO", "NH[4]^+{}", "Water ~~ temperature", "NO[2]^-{}","pH", "PO[4]^3^-{}", "NO[3]^-{}", "COD" , 
                "Turbility", "Air ~ temperature", "Chlorophyll*~alpha", "Wind ~ velocity", "Solar ~ radiation")
labels_CO2 <- rev(labels_CO2)
labels_CO2_parse <- parse(text = labels_CO2)

ggsave("RF_CO2_opt.tiff", ggplot(featureImportance_CO2_opt[1:13,], aes(x=reorder(Feature, Importance), y=Importance)) +
           geom_bar(stat="identity", fill="tomato") +
           coord_flip() + 
           theme_bw(base_size=20) +
           scale_x_discrete(labels = labels_CO2_parse) +
           labs(x = NULL, y = "Importance") + 
           ggtitle("Feature Importance") +
           theme(plot.title=element_text(size=18)),
       units = 'cm', height = 20, width = 20, dpi = 300)
#### RF - MLR-2 ####

colnames(river_RF_2)[1] <- "Water_temperature"
colnames(river_RF_2)[6:7] <- c("Left_bank", "Right_bank")
colnames(river_RF_2)[10:13] <- c("Flow_variability", "Average_depth", "Average_velocity", "Pool_class")
colnames(river_RF_2)[19] <- "Air_temperature"
colnames(river_RF_2)[20] <- "Wind_velocity"
colnames(river_RF_2)[21] <- "Solar_radiation"

#*** N2O-2 ####
# creat mlr task and convert factors to dummy features

trainTask_N2O_2 = createDummyFeatures(river_RF_2[,-c(23:24)],target = "N2O")
trainTask_N2O_2 <- makeRegrTask(data = river_RF_2[,-c(23:24)],target = "N2O")

# create mlr learner
set.seed(1234)
lrn_N2O_2 <- makeLearner("regr.ranger")
# lrn_N2O_2$par.vals <- list(ntree = 100L, importance = TRUE)
cv_N2O_2 <- makeResampleDesc(method = "LOO")
# set parallel backend
parallelStartSocket(cpus = detectCores()-1)
res_N2O_2 <- resample(learner = lrn_N2O_2, task = trainTask_N2O_2, resampling = cv_N2O_2)
# Tuning hyperparameters

# Parameter Tuning: Mainly, there are three parameters in the random forest algorithm which you should look at (for tuning):
# ntree - The number of trees to grow. Larger the tree, it will be more computationally expensive to build models.
# mtry - It refers to how many variables we should select at a node split. 
# Also as mentioned above, the default value is p/3 for regression and sqrt(p) for classification. 
# We should always try to avoid using smaller values of mtry to avoid overfitting.
# nodesize - It refers to how many observations we want in the terminal nodes.
# This parameter is directly related to tree depth. Higher the number, lower the tree depth. 
# With lower tree depth, the tree might even fail to recognize useful signals from the data.

# To know which hyperparameter can be tuned using getParamSet
getParamSet(lrn_N2O_2)

# Only tune three abovementioned hyperparameters

params_N2O_2 <- makeParamSet(makeIntegerParam("mtry", lower = 2, upper = 10),
                           makeIntegerParam("min.node.size", lower = 2, upper = 25))
tc_N2O_2 <- makeTuneControlMBO(budget = 100)
tr_N2O_2 <- tuneParams(learner = lrn_N2O_2, task = trainTask_N2O_2, resampling = cv_N2O_2,
                    par.set = params_N2O_2, control = tc_N2O_2, show.info = T)

parallelStop()
tr_N2O_2
# Tune result:
# Op. pars: mtry=10; min.node.size=2
# mse.test.mean=0.4196192
# Apply the optimal RF


# using randomForest package
set.seed(1234)
regressor_N2O_2_opt <- randomForest(river_RF_2[,-c(22:24)], river_RF$N2O, ntree = 1000, importance = TRUE, mtry = 10, nodesize = 2)
imp_N2O_2_opt <- randomForest::importance(regressor_N2O_2_opt, type=1)
featureImportance_N2O_2_opt <- data.frame(Feature=row.names(imp_N2O_2_opt), Importance=imp_N2O_2_opt[,1])

ggsave("RF_N2O_2_opt.tiff", ggplot(featureImportance_N2O_2_opt[,], aes(x=reorder(Feature, Importance), y=Importance)) +
           geom_bar(stat="identity", fill="tomato") +
           coord_flip() + 
           theme_bw(base_size=20) +
           # scale_x_discrete(labels = labels_N2O_2_parse) +
           labs(x = NULL, y = "Scale Importance") + 
           ggtitle("Random Forest Feature Importance") +
           theme(plot.title=element_text(size=18)),
       units = 'cm', height = 20, width = 20, dpi = 300)
# using ranger package with permutation feature importance
set.seed(1234)
regressor_N2O_2_opt_ranger <- ranger(formula = N2O~., data = river_RF_2[,-c(23:24)], num.trees =  1000, mtry = 10, 
                                     importance = 'permutation', min.node.size = 2) 
imp_N2O_2_opt_per <- importance_pvalues(x = regressor_N2O_2_opt_ranger, method = 'janitza', num.permutations = 100)

featureImportance_N2O_2_opt_per <- data.frame(Feature=row.names(imp_N2O_2_opt_per), Importance=imp_N2O_2_opt_per[,1], 
                                              pvalue=imp_N2O_2_opt_per[,2])
# remove the variables with pvalue >0.05
featureImportance_N2O_2_opt_per <- featureImportance_N2O_2_opt_per %>% filter(pvalue <=0.05)

labels_N2O_2_per <- c("NH[4]^+{}", "DO", "Turbidity", "Water ~ temperature", "NO[2]^-{}", "PO[4]^3^-{}", "COD", "NO[3]^-{}", 
                 "Average ~ velocity", "Solar ~ radiation", "Average ~ depth", "Wind ~ velocity", "Air ~ temperature")
labels_N2O_2_per <- rev(labels_N2O_2_per)
labels_N2O_2_per_parse <- parse(text = labels_N2O_2_per)

ggsave("RF_N2O_2_opt_per.tiff", ggplot(featureImportance_N2O_2_opt_per[,], aes(x=reorder(Feature, Importance), y=Importance)) +
           geom_bar(stat="identity", fill="tomato") +
           coord_flip() + 
           theme_bw(base_size=20) +
           scale_x_discrete(labels = labels_N2O_2_per_parse) +
           labs(x = NULL, y = "Scale Importance") + 
           ggtitle(bquote(""~N[2]~"O")) +
           theme(plot.title=element_text(size=18)),
       units = 'cm', height = 20, width = 20, dpi = 300)

#*** CH4-2 ####
# creat mlr task and convert factors to dummy features

trainTask_CH4_2 = createDummyFeatures(river_RF_2 %>% select(-N2O,-CO2),target = "CH4")
trainTask_CH4_2 <- makeRegrTask(data = river_RF_2 %>% select(-N2O,-CO2),target = "CH4")

# create mlr learner
set.seed(1234)
lrn_CH4_2 <- makeLearner("regr.ranger")
# lrn_CH4_2$par.vals <- list(ntree = 100L, importance = TRUE)
cv_CH4_2 <- makeResampleDesc(method = "LOO")
# set parallel backend
parallelStartSocket(cpus = detectCores()-1)
res_CH4_2 <- resample(learner = lrn_CH4_2, task = trainTask_CH4_2, resampling = cv_CH4_2)
# Tuning hyperparameters

# Only tune three abovementioned hyperparameters

params_CH4_2 <- makeParamSet(makeIntegerParam("mtry", lower = 2, upper = 10),
                           makeIntegerParam("min.node.size", lower = 2, upper = 25))
tc_CH4_2 <- makeTuneControlMBO(budget = 100)
tr_CH4_2 = tuneParams(learner = lrn_CH4_2, task = trainTask_CH4_2, resampling = cv_CH4_2,
                    par.set = params_CH4_2, control = tc_CH4_2, show.info = F)

parallelStop()
tr_CH4_2
# Tune result:
# mtry=6; min.node.size=2
# mse.test.mean=0.2276998

lrn_CH4_2_tuned <- setHyperPars(learner = lrn_CH4_2,
                      par.vals = tr_CH4_2$x)
resample(learner = lrn_CH4_2_tuned, task = trainTask_CH4_2, resampling = cv_CH4_2)

# Apply the optimal RF

set.seed(1234)
regressor_CH4_2_opt <- randomForest(river_RF_2[,-c(22:24)], river_RF$CH4, ntree = 1000, importance = TRUE, mtry = 6, nodesize = 2)

imp_CH4_2_opt <- randomForest::importance(regressor_CH4_2_opt, type=1)
featureImportance_CH4_2_opt <- data.frame(Feature=row.names(imp_CH4_2_opt), Importance=imp_CH4_2_opt[,1])
# remove the negative importance values which shows the unimportant variable

labels_CH4_2 <- c("DO", "Turbidity",  "NO[2]^-{}", "NO[3]^-{}",  "PO[4]^3^-{}", "Water ~ temperature", "COD" , "Air ~ temperature",
                "NH[4]^+{}","pH",  "Wind ~ velocity", "Solar ~ radiation", "Chlorophyll*~alpha")
labels_CH4_2 <- rev(labels_CH4_2)
labels_CH4_2_parse <- parse(text = labels_CH4_2)

ggsave("RF_CH4_2_opt.tiff",ggplot(featureImportance_CH4_2_opt[,], aes(x=reorder(Feature, Importance), y=Importance)) +
           geom_bar(stat="identity", fill="tomato") +
           coord_flip() + 
           theme_bw(base_size=20) +
           # scale_x_discrete(labels = labels_CH4_2_parse) +
           labs(x = NULL, y = "Importance") + 
           ggtitle("Random Forest Feature Importance") +
           theme(plot.title=element_text(size=18)),
       units = 'cm', height = 20, width = 20, dpi = 300)

# using ranger package with permutation feature importance
set.seed(1234)
regressor_CH4_2_opt_ranger <- ranger(formula = CH4~., data = river_RF_2[,-c(22,24)], num.trees =  1000, mtry = 6, 
                                     importance = 'permutation', min.node.size = 2) 
imp_CH4_2_opt_per <- importance_pvalues(x = regressor_CH4_2_opt_ranger, method = 'janitza', num.permutations = 100)

featureImportance_CH4_2_opt_per <- data.frame(Feature=row.names(imp_CH4_2_opt_per), Importance=imp_CH4_2_opt_per[,1], 
                                              pvalue=imp_CH4_2_opt_per[,2])
# remove the variables with pvalue >0.05
featureImportance_CH4_2_opt_per <- featureImportance_CH4_2_opt_per %>% filter(pvalue <=0.05)

labels_CH4_2_per <- c("DO", "NO[2]^-{}", "Turbidity", "Average ~ velocity", "Average ~ depth", "NO[3]^-{}", "PO[4]^3^-{}", "Water ~ temperature", 
                      "COD", "Air ~ temperature", "NH[4]^+{}",  "Solar ~ radiation", "Left ~ bank", "Wind ~ velocity")
labels_CH4_2_per <- rev(labels_CH4_2_per)
labels_CH4_2_per_parse <- parse(text = labels_CH4_2_per)

ggsave("RF_CH4_2_opt_per.tiff", ggplot(featureImportance_CH4_2_opt_per[,], aes(x=reorder(Feature, Importance), y=Importance)) +
           geom_bar(stat="identity", fill="tomato") +
           coord_flip() + 
           theme_bw(base_size=20) +
           scale_x_discrete(labels = labels_CH4_2_per_parse) +
           labs(x = NULL, y = "Scale Importance") + 
           ggtitle(bquote("CH"[4])) +
           theme(plot.title=element_text(size=18)),
       units = 'cm', height = 20, width = 20, dpi = 300)

#*** CO2-2 ####
# creat mlr task and convert factors to dummy features

trainTask_CO2_2 <- createDummyFeatures(river_RF_2 %>% select(-N2O,-CH4),target = "CO2")
trainTask_CO2_2 <- makeRegrTask(data = river_RF_2 %>% select(-N2O,-CH4),target = "CO2")

# create mlr learner
set.seed(1234)
lrn_CO2_2 <- makeLearner("regr.ranger")
# lrn_CO2_2$par.vals <- list(ntree = 100L, importance = TRUE)
cv_CO2_2 <- makeResampleDesc(method = "LOO")
# set parallel backend
parallelStartSocket(cpus = detectCores()-1)
res_CO2_2 <- resample(learner = lrn_CO2_2, task = trainTask_CO2_2, resampling = cv_CO2_2)
# Tuning hyperparameters

# Only tune three abovementioned hyperparameters

params_CO2_2 <- makeParamSet(makeIntegerParam("mtry", lower = 2, upper = 10),
                           makeIntegerParam("min.node.size", lower = 2, upper = 25))
tc_CO2_2 <- makeTuneControlMBO(budget = 100)
tr_CO2_2 <- tuneParams(learner = lrn_CO2_2, task = trainTask_CO2_2, resampling = cv_CO2_2,
                    par.set = params_CO2_2, control = tc_CO2_2, show.info = F)

parallelStop()
tr_CO2_2
# Tune result:
# mtry=10; min.node.size=2
# mse.test.mean=0.221592

# Apply the optimal RF

set.seed(1234)
regressor_CO2_2_opt <- randomForest(river_RF_2[,-c(22:24)], river_RF$CO2, ntree = 500, importance = TRUE, mtry = 9, nodesize = 2)

imp_CO2_2_opt <- randomForest::importance(regressor_CO2_2_opt, type=1)
featureImportance_CO2_2_opt <- data.frame(Feature=row.names(imp_CO2_2_opt), Importance=imp_CO2_2_opt[,1])

labels_CO2_2 <- c("DO", "NH[4]^+{}", "Water ~~ temperature", "NO[2]^-{}","pH", "PO[4]^3^-{}", "NO[3]^-{}", "COD" , 
                "Turbidity", "Air ~ temperature", "Chlorophyll*~alpha", "Wind ~ velocity", "Solar ~ radiation")
labels_CO2_2 <- rev(labels_CO2_2)
labels_CO2_2_parse <- parse(text = labels_CO2_2)

ggsave("RF_CO2_2_opt.tiff", ggplot(featureImportance_CO2_2_opt[,], aes(x=reorder(Feature, Importance), y=Importance)) +
           geom_bar(stat="identity", fill="tomato") +
           coord_flip() + 
           theme_bw(base_size=20) +
           # scale_x_discrete(labels = labels_CO2_2_parse) +
           labs(x = NULL, y = "Importance") + 
           ggtitle("Feature Importance") +
           theme(plot.title=element_text(size=18)),
       units = 'cm', height = 20, width = 20, dpi = 300)


# using ranger package with permutation feature importance
set.seed(1234)
regressor_CO2_2_opt_ranger <- ranger(formula = CO2~., data = river_RF_2[,-c(22,23)], num.trees =  1000, mtry = 6, 
                                     importance = 'permutation', min.node.size = 2) 
imp_CO2_2_opt_per <- importance_pvalues(x = regressor_CO2_2_opt_ranger, method = 'janitza', num.permutations = 100)

featureImportance_CO2_2_opt_per <- data.frame(Feature=row.names(imp_CO2_2_opt_per), Importance=imp_CO2_2_opt_per[,1], 
                                              pvalue=imp_CO2_2_opt_per[,2])
# remove the variables with pvalue >0.05
featureImportance_CO2_2_opt_per <- featureImportance_CO2_2_opt_per %>% filter(pvalue <=0.05)

labels_CO2_2_per <- c("DO", "NH[4]^+{}",  "Water ~ temperature", "NO[2]^-{}", "PO[4]^3^-{}", "Average ~ velocity", "Average ~ depth", 
                      "pH ", "NO[3]^-{}", "Turbidity",  "COD", "Right ~ bank", "Chlorophyll*~alpha")
labels_CO2_2_per <- rev(labels_CO2_2_per)
labels_CO2_2_per_parse <- parse(text = labels_CO2_2_per)

ggsave("RF_CO2_2_opt_per.tiff", ggplot(featureImportance_CO2_2_opt_per[,], aes(x=reorder(Feature, Importance), y=Importance)) +
           geom_bar(stat="identity", fill="tomato") +
           coord_flip() + 
           theme_bw(base_size=20) +
           scale_x_discrete(labels = labels_CO2_2_per_parse) +
           labs(x = NULL, y = "Scale Importance") + 
           ggtitle(bquote("CO"[2])) +
           theme(plot.title=element_text(size=18)),
       units = 'cm', height = 20, width = 20, dpi = 300)
