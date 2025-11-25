# Project: Puerto Rican Solar Power Predictability Analysis
# Code by: Renee Obringer
# Last Run: 17 November 2025

# ORGANIZATION: 
# This code is organized into sections, the start of each is denoted by multiple #
# The sections can be run independently by loading the rdata files at the beginning of each section
# Each section is described below
#
# LOAD DATA: load data
# DATA PRE-PROCESSING: pre-processing to integrate solar power and weather data
# SENSITIVITY ANALYSIS: senstivity analysis for some of the assumptions in the PVGIS simulation
# MODELING: run random forest for each municipality 
#           **NOTE: this runs 5-fold cross validation for 15 years of daily data for 228 model runs 
#             (76 municipalities x 3 panels), so it takes some time (~3 hours)
# INTERPRETATION: interpretation of results (significance tests, etc.)
# FIGURES AND TABLES: plotting figures and creating tables included in manuscript 

rm(list=ls())
options(scipen = 999)

# libraries
library(stringr)           # for working with strings
library(readxl)            # for reading in excel files
library(tidyverse)         # for data organization
library(reshape2)          # for data organization
library(randomForest)      # random forest algorithm
library(ggplot2)           # for general plotting
library(cowplot)           # for plotting multiple panels
library(sf)                # for plotting maps with ggplot
library(ggspatial)         # for plotting maps with ggplot
library(lubridate)         # for working with dates
library(humidity)          # for calculating humidity
library(corrplot)          # for creating correlation plots
library(transport)         # for calculating Wasserstein Distance

# directories
datadir1 <- '/Users/rqo5125/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/data/solarpower'
datadir2 <- '/Users/rqo5125/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/data/daynet_1kmweather/PuertoRico'
rdatadir <- '/Users/rqo5125/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/2025_26/papers/PRsolarpower/rdatafiles'
outputdir <- '/Users/rqo5125/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/2025_26/papers/PRsolarpower/results'
figuredir <- '/Users/rqo5125/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/2025_26/papers/PRsolarpower/figures'

################ LOAD DATA ################

# Weather data
setwd(datadir2)

# initialize variables
allnarrdata <- list()

# loop through each municipality
for (i in 1:76) {
  
  # set municipality number
  if (i < 10) {
    municipnum <- paste(0, i, sep = '')
  } else {
    municipnum <- as.character(i)
  }
  
  # get filenames for a given municipality
  filenames <- list.files(pattern = municipnum, full.names=F)
  
  # initialize data
  narrdata <- list(); narrdatanames <- c()
  
  # loop through each file
  for (j in 1:length(filenames)) {
    # read in data
    narrdata[[j]] <- read.csv(paste(datadir2, '/', filenames[j], sep = ""), header=T, fill = TRUE, stringsAsFactors = F)
    narrdatanames[j] <- names(narrdata[[j]])[3]
  }
  
  # reorganize data
  narrdataframe <- data.frame(narrdata[[1]][,2:3], narrdata[[2]][,3], narrdata[[3]][,3], narrdata[[4]][,3], narrdata[[5]][,3], narrdata[[6]][,3])
  names(narrdataframe) <- c('date', narrdatanames)
  
  # store as list
  allnarrdata[[i]] <- narrdataframe
}

# solar power data
setwd(datadir1)

# get filenames
filenames <- list.files(pattern = '*.xlsx', full.names=F)

# initialize variables
solardata <- list()

# loop through all filenames
for (i in 1:length(filenames)) {
  
  # get solar panel name
  solarpanel <- str_remove(filenames[i],'_35deg_10percloss_power_2005-20.xlsx')
  
  # read in data
  solardata[[solarpanel]] <- read_excel(paste(datadir1, '/', filenames[i], sep = ""))
}

# solar power data (sensitivity analysis)
setwd(paste(datadir1, '/sensitivityanalysis', sep = ''))

# get filenames
filenames <- list.files(pattern = '*.xlsx', full.names=F)

# initialize variables
sensitivityanalysis <- list()

# loop through all filenames
for (i in 1:length(filenames)) {
  
  # get experiment name
  experiment <- str_remove(filenames[i],'_power_2015-20.xlsx')
  
  # read in data
  sensitivityanalysis[[i]] <- read_excel(paste(datadir1, '/sensitivityanalysis/', filenames[i], sep = ""))
  sensitivityanalysis[[i]]$experiment <- experiment
}

setwd(rdatadir)
save(solardata, allnarrdata, sensitivityanalysis, file = 'rawdata.rdata')

################ DATA PRE-PROCESSING ################

setwd(rdatadir)
load('rawdata.rdata')

# aggregate to daily narr values

# initialize data
dailynarrdata <- list()

# loop through all municipalities
for (i in 1:length(allnarrdata)) {
  municdata <- allnarrdata[[i]]
  
  # convert units
  municdata$temp_degC <- municdata$temp_K - 273.15

  # calculate relative humidity from clausius-clapeyron - https://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
  #vaporpressure <- municdata$pressure_Pa*municdata$sphum_kgPkg/(0.622 + municdata$sphum_kgPkg)
  #satvaporpressure <- 611 * exp(17.67 * (municdata$temp_K - 273.16)/(municdata$temp_K - 29.65))
  #municdata$rhum <- (vaporpressure/satvaporpressure)*100
  
  municdata$rhum_perc <- SH2RH(municdata$sphum_kgPkg, municdata$temp_K, municdata$pressure_Pa, isK = T)
  
  # remove unnecessary columns
  municdata <- municdata[,-c(3,5,6)]
  
  # remove time from date column
  municdata$date <- sub(' .*', '', municdata$date)
  
  # aggregate some variables by min, mean, and max
  minmeanmax <- municdata %>% group_by(date) %>% summarize_at(vars(windspeed_mPs:rhum_perc), list(min = min, mean = mean, max = max))
  minmeanmax$date <- as.Date(minmeanmax$date)
  
  # aggregate some variables by total
  total <- municdata %>% group_by(date) %>% summarize_at(vars(precip_mmPs, rad_wPm2), list(total = sum))
  total$date <- as.Date(total$date)
  
  # add empty rows for missing dates (February 29th)
  ts <- seq(ymd("2005-01-01"), ymd("2019-12-31"), by="day")
  total <- tidyr::complete(total, date = ts, fill = list(windspeed_mPs_min = NA, temp_degC_min = NA, rhum_perc_min = NA,
                                                         windspeed_mPs_mean = NA, temp_degC_mean = NA, rhum_perc_mean = NA,
                                                         windspeed_mPs_max = NA, temp_degC_max = NA, rhum_perc_max = NA, 
                                                         precip_mmPs_total = NA, rad_wPm2_total = NA))
  minmeanmax <- tidyr::complete(minmeanmax, date = ts, fill = list(windspeed_mPs_min = NA, temp_degC_min = NA, rhum_perc_min = NA,
                                                                   windspeed_mPs_mean = NA, temp_degC_mean = NA, rhum_perc_mean = NA,
                                                                   windspeed_mPs_max = NA, temp_degC_max = NA, rhum_perc_max = NA, 
                                                                   precip_mmPs_total = NA, rad_wPm2_total = NA))
  
  # combine to final narr dataset
  dailynarrdata[[i]] <- reduce(list(minmeanmax, total), left_join, by = 'date')
}

# merge datasets

# initialize variables
allpaneldata <- list()

# loop through each panel 
for (i in 1:length(solardata)) {
  panel <- solardata[[i]]
  
  # initialize variables
  mergeddata <- list()
  
  # loop through each municipality
  for (j in 1:length(dailynarrdata)) {
    # merge solar power and narr data for each panel-municipality combination
    mergeddata[[j]] <- cbind(dailynarrdata[[j]], panel[1:5478,j+1]) # remove 2020 data
    names(mergeddata[[j]])[length(mergeddata[[j]])] <- 'power_watts'
  }
  
  # store results
  allpaneldata[[names(solardata)[i]]] <- mergeddata
}

setwd(rdatadir)
save(dailynarrdata, allpaneldata, file = 'processeddata.rdata')

################ SENSITIVITY ANALYSIS ################

setwd(rdatadir)
load('rawdata.rdata')

# Figure showing the difference between time series

# organize data
plotdata <- c()

for (i in 1:length(sensitivityanalysis)) {
  expdata <- sensitivityanalysis[[i]]
  expdata <- separate(expdata, 'experiment', into = c('panel', 'tilt', 'loss'), sep = '_')
  names(expdata)[2:77] <- c(1:76)
  expdata <- melt(expdata, id = c('Day', 'panel', 'tilt', 'loss'), variable.name = 'municipality', value.name = 'energy')
  
  plotdata <- rbind(plotdata, expdata)
}

setwd(figuredir)
pdf('CdTe_sensitivityanalysis.pdf', width = 15, height = 15)
ggplot(plotdata[which(plotdata$panel == 'CdTe'),]) + 
  geom_line(aes(x = Day, y = energy, group = municipality), color = '#D3D3D3') +
  xlab('Day') + ylab('Energy Generated (Wh/day)') + ggtitle('CdTe') +
  facet_wrap(~ tilt + loss, nrow = 4, scales = 'free') + theme_light(base_size = 12) 
dev.off()

setwd(figuredir)
pdf('CIS_sensitivityanalysis.pdf', width = 15, height = 15)
ggplot(plotdata[which(plotdata$panel == 'CIS'),]) + 
  geom_line(aes(x = Day, y = energy, group = municipality), color = '#D3D3D3') +
  xlab('Day') + ylab('Energy Generated (Wh/day)') + ggtitle('CIS') +
  facet_wrap(~ tilt + loss, nrow = 4, scales = 'free') + theme_light(base_size = 12) 
dev.off()

setwd(figuredir)
pdf('CrystSI_sensitivityanalysis.pdf', width = 15, height = 15)
ggplot(plotdata[which(plotdata$panel == 'crystSi'),]) + 
  geom_line(aes(x = Day, y = energy, group = municipality), color = '#D3D3D3') +
  xlab('Day') + ylab('Energy Generated (Wh/day)') + ggtitle('CrystSI') +
  facet_wrap(~ tilt + loss, nrow = 4, scales = 'free') + theme_light(base_size = 12) 
dev.off()

# stability tests for reduced form model (Wasserstein-Fourier Distance; Cazelles et al., 2021)

# get averages
averages <- plotdata %>% group_by(Day, panel, tilt, loss) %>% summarize(mean(energy))
names(averages)[5] <- 'energy'

panellist <- c('CdTe', 'CIS', 'crystSi')
tiltlist <- c('15deg', '25deg', '35deg', '45deg')
losslist <- c('1percloss', '5percloss', '10percloss', '15percloss')
allwd <- data.frame(panel = c(), tilt = c(), loss1 = c(), loss2 = c(), wd = c())

# loop through panels
for (p in 1:length(panellist)) {
  # loop through tilt values
  for (t in 1:length(tiltlist)) {
    allfreq <- data.frame(percloss1 = numeric(937), percloss5 = numeric(937), percloss10 = numeric(937), percloss15 = numeric(937))
    allnpsd <- data.frame(percloss1 = numeric(937), percloss5 = numeric(937), percloss10 = numeric(937), percloss15 = numeric(937))
    # loop through loss values
    for (l in 1:length(losslist)) {
      # calculate normalized power spectral density
      psd <- spec.pgram(averages[which(averages$panel == panellist[p] & averages$tilt == tiltlist[t] & averages$loss == losslist[l]),5], plot = FALSE)$spec
      allfreq[,l] <- spec.pgram(averages[which(averages$panel == panellist[p] & averages$tilt == tiltlist[t] & averages$loss == losslist[l]),5], plot = FALSE)$freq
      allnpsd[,l] <- psd/sum(psd)
    }
    
    wdvalues1vA <- data.frame(panel = panellist[p], tilt = tiltlist[t], loss1 = names(allfreq)[1], loss2 = names(allfreq)[2:4], 
                              wd = c(wasserstein1d(a = allfreq[,1], b = allfreq[,2], wa = allnpsd[,1], wb = allnpsd[,2], p = 2),
                                     wasserstein1d(a = allfreq[,1], b = allfreq[,3], wa = allnpsd[,1], wb = allnpsd[,3], p = 2),
                                     wasserstein1d(a = allfreq[,1], b = allfreq[,4], wa = allnpsd[,1], wb = allnpsd[,4], p = 2)))
    
    wdvalues5vA <- data.frame(panel = panellist[p], tilt = tiltlist[t], loss1 = names(allfreq)[2], loss2 = names(allfreq)[3:4], 
                              wd = c(wasserstein1d(a = allfreq[,2], b = allfreq[,3], wa = allnpsd[,2], wb = allnpsd[,3], p = 2),
                                     wasserstein1d(a = allfreq[,2], b = allfreq[,4], wa = allnpsd[,2], wb = allnpsd[,4], p = 2)))
    
    wdvalues10vA <- data.frame(panel = panellist[p], tilt = tiltlist[t], loss1 = names(allfreq)[3], loss2 = names(allfreq)[4], 
                               wd = wasserstein1d(a = allfreq[,3], b = allfreq[,4], wa = allnpsd[,3], wb = allnpsd[,4], p = 2))
    
    # store each panel-tilt-loss combination
    allwd <- rbind(allwd, wdvalues1vA, wdvalues5vA, wdvalues10vA)
  }
}

setwd(outputdir)
write.csv(allwd, file = 'allwassersteindistances.csv')

# plot wasserstein distances

allwd_long <- allwd
allwd_long$rows <- sub("percloss", "", allwd_long$loss1)
allwd_long$cols <- sub("percloss", "", allwd_long$loss2)
allwd_long$tilt <- as.numeric(sub('deg', '', allwd_long$tilt))

allwd_long_reverse <- allwd_long
allwd_long_reverse$rows <- allwd_long$cols
allwd_long_reverse$cols <- allwd_long$rows

combinedwd <- rbind(allwd_long, allwd_long_reverse)

combinedwd$rows <- as.factor(combinedwd$rows)
combinedwd$cols <- as.factor(combinedwd$cols)

p1 <- ggplot(data = combinedwd[which(combinedwd$panel == 'CdTe'),], aes(x = rows, y = cols)) +
  geom_tile(aes(fill = wd), colour = "white") + geom_text(aes(label = sprintf("%1.2f",wd/10e-7)), vjust = 1) +
  scale_fill_gradient(low = "#d0d1e6", high = "#023858") + theme_light() + theme(legend.position = "none") +
  facet_wrap(~tilt, nrow = 1) + ylab('Loss Value (%)') + xlab('Loss Value (%)') + ggtitle('CdTe') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

p2 <- ggplot(data = combinedwd[which(combinedwd$panel == 'CIS'),], aes(x = rows, y = cols)) +
  geom_tile(aes(fill = wd), colour = "white") + geom_text(aes(label = sprintf("%1.2f",wd/10e-7)), vjust = 1) +
  scale_fill_gradient(low = "#d0d1e6", high = "#023858") + theme_light() + theme(legend.position = "none") +
  facet_wrap(~tilt, nrow = 1) + ylab('Loss Value (%)') + xlab('Loss Value (%)') + ggtitle('CIS') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

p3 <- ggplot(data = combinedwd[which(combinedwd$panel == 'crystSi'),], aes(x = rows, y = cols)) +
  geom_tile(aes(fill = wd), colour = "white") + geom_text(aes(label = sprintf("%1.2f",wd/10e-7)), vjust = 1) +
  scale_fill_gradient(low = "#d0d1e6", high = "#023858") + theme_light() + theme(legend.position = "none") +
  facet_wrap(~tilt, nrow = 1) + ylab('Loss Value (%)') + xlab('Loss Value (%)') + ggtitle('CrystSi') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

setwd(figuredir)
pdf('wassersteindistances.pdf', width = 10, height = 8)
plot_grid(p1, p2, p3, nrow = 3)
dev.off()


################ MODELING ################

setwd(rdatadir)
load('processeddata.rdata')

# variable names
seasonnames <- c('Rainy', 'Dry')

# initialize variables
panelrsq <- list(); panelrmse <- list(); panelyhat <- list(); paneltest <- list(); panelmodel <- list()

# loop through each panel
for (i in 1:length(allpaneldata)) {
  panel <- allpaneldata[[i]]
  
  # initialize variables
  municrsq <- list(); municnrmse <- list(); municyhat <- list(); munictest <- list(); municmodel <- list()
  
  # for tracking progress of loop
  print(paste('panel: ', i, sep = ''))
  
  # loop through each municipality
  for (j in 1:length(panel)) {
    municipality <- panel[[j]]
    
    # extract rainy season months
    municipality$wetseason <- (month(municipality$date) >= 5 & month(municipality$date) <= 11)
    
    # for tracking progress of loop
    print(paste('municipality: ', j, sep = ''))
    print(paste('time: ', Sys.time(), sep = ''))
    
    # set up seasonal analysis & remove date/season columns
    rainyseason <- municipality[which(municipality$wetseason == T), -c(1,14)]
    dryseason <- municipality[which(municipality$wetseason == F), -c(1,14)]
    
    seasons <- list(rainyseason, dryseason)
    
    # initialize variables
    seasonrsq <- list(); seasonnrmse <- list(); seasonyhat <- list(); seasontest <- list(); seasonmodel <- list()
    
    for (s in 1:length(seasons)) {
      seasondata <- na.omit(seasons[[s]])
      
      # set up 5-fold cross validation with selected variables
      k <- 5
      n <- nrow(seasondata)
      set.seed(11)
      seasondata <- seasondata[sample(n),]
      folds <- cut(seq(1,n), breaks = k, labels = FALSE)
      
      # initialize variables
      cvrsq <- c(); cvnrmse <- c(); rfyhat_all <- c(); testdata_all <- c()
      
      # run model
      for (l in 1:k) {
        # split into training and test dataset
        testIndex <- which(folds == l, arr.ind = TRUE)
        testData <- seasondata[testIndex,]
        trainData <- seasondata[-testIndex,]
        
        # train model
        rfmodel <- randomForest(power_watts ~ ., data = trainData)
        
        # test model
        rfyhat <- predict(rfmodel, newdata = testData)
        
        # calculate error
        cvrsq[l] <- cor(rfyhat, testData$power_watts)^2
        cvnrmse[l] <- sqrt(sum((rfyhat-testData$power_watts)^2)/nrow(testData))/(range(testData$power_watts)[2]-range(testData$power_watts)[1])
        
        # store important variables
        rfyhat_all <- append(rfyhat_all, rfyhat)
        testdata_all <- append(testdata_all, testData$power_watts)
      }
      
      # plot partial dependence
      imp <- importance(rfmodel)
      impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
      
      setwd(paste(figuredir, '/PDPplots', sep = ''))
      for (v in 1:4) {
        pdf(paste('pdp_', names(allpaneldata)[i], '_m_', j, '_', impvar[v], '_', seasonnames[s], '.pdf', sep = ''))
        partialPlot(rfmodel, trainData, impvar[v], xlab=impvar[v],
                    main=paste("Partial Dependence on", impvar[v]))
        dev.off()
      }
      
      # store data for each season
      seasonrsq[[s]] <- cvrsq
      seasonnrmse[[s]] <- cvnrmse
      seasonyhat[[s]] <- rfyhat_all
      seasontest[[s]] <- testdata_all
      seasonmodel[[s]] <- rfmodel
    }
    
    # store data for each municipality
    municrsq[[j]] <- seasonrsq
    municnrmse[[j]] <- seasonnrmse
    municyhat[[j]] <- seasonyhat
    munictest[[j]] <- seasontest
    municmodel[[j]] <- seasonmodel
  }
  
  # store data for each panel
  panelrsq[[names(allpaneldata)[i]]] <- municrsq
  panelrmse[[names(allpaneldata)[i]]] <- municnrmse
  panelyhat[[names(allpaneldata)[i]]] <- municyhat
  paneltest[[names(allpaneldata)[i]]] <- munictest
  panelmodel[[names(allpaneldata)[i]]] <- municmodel
}

setwd(rdatadir)
save(panelrsq, panelrmse, panelyhat, paneltest, panelmodel, file = 'modelresults.rdata')

################ INTERPRETATION ################

setwd(rdatadir)
load('modelresults.rdata')

allnrmse_r <- data.frame('CdTe_min' = sapply(1:76, function(x) min(panelrmse[[1]][[x]][[1]])), 'CdTe_mean' = sapply(1:76, function(x) mean(panelrmse[[1]][[x]][[1]])), 'CdTe_max' = sapply(1:76, function(x) max(panelrmse[[1]][[x]][[1]])),
                         'CIS_min' = sapply(1:76, function(x) min(panelrmse[[2]][[x]][[1]])), 'CIS_mean' = sapply(1:76, function(x) mean(panelrmse[[2]][[x]][[1]])), 'CIS_max' = sapply(1:76, function(x) max(panelrmse[[2]][[x]][[1]])),
                         'crystSI_min' = sapply(1:76, function(x) min(panelrmse[[3]][[x]][[1]])), 'crystSI_mean' = sapply(1:76, function(x) mean(panelrmse[[3]][[x]][[1]])), 'crystSI_max' = sapply(1:76, function(x) max(panelrmse[[3]][[x]][[1]])),
                         'Region' = latlon$region, 'Municipality' = latlon$name)

allnrmse_r$bestpanel_mean <- case_when(abs(allnrmse_r$CdTe_mean) < abs(allnrmse_r$CIS_mean) & abs(allnrmse_r$CdTe_mean) < abs(allnrmse_r$crystSI_mean) ~ 'CdTe',
                                       abs(allnrmse_r$CIS_mean) < abs(allnrmse_r$CdTe_mean) & abs(allnrmse_r$CIS_mean) < abs(allnrmse_r$crystSI_mean) ~ 'CIS',
                                       abs(allnrmse_r$crystSI_mean) < abs(allnrmse_r$CdTe_mean) & abs(allnrmse_r$crystSI_mean) < abs(allnrmse_r$CIS_mean) ~ 'CrystSi')

allnrmse_r$bestpanel_min <- case_when(abs(allnrmse_r$CdTe_min) < abs(allnrmse_r$CIS_min) & abs(allnrmse_r$CdTe_min) < abs(allnrmse_r$crystSI_min) ~ 'CdTe',
                                      abs(allnrmse_r$CIS_min) < abs(allnrmse_r$CdTe_min) & abs(allnrmse_r$CIS_min) < abs(allnrmse_r$crystSI_min) ~ 'CIS',
                                      abs(allnrmse_r$crystSI_min) < abs(allnrmse_r$CdTe_min) & abs(allnrmse_r$crystSI_min) < abs(allnrmse_r$CIS_min) ~ 'CrystSi')

allnrmse_r$bestpanel_max <- case_when(abs(allnrmse_r$CdTe_max) < abs(allnrmse_r$CIS_max) & abs(allnrmse_r$CdTe_max) < abs(allnrmse_r$crystSI_max) ~ 'CdTe',
                                      abs(allnrmse_r$CIS_max) < abs(allnrmse_r$CdTe_max) & abs(allnrmse_r$CIS_max) < abs(allnrmse_r$crystSI_max) ~ 'CIS',
                                      abs(allnrmse_r$crystSI_max) < abs(allnrmse_r$CdTe_max) & abs(allnrmse_r$crystSI_max) < abs(allnrmse_r$CIS_max) ~ 'CrystSi')



allnrmse_d <- data.frame('CdTe_min' = sapply(1:76, function(x) min(panelrmse[[1]][[x]][[2]])), 'CdTe_mean' = sapply(1:76, function(x) mean(panelrmse[[1]][[x]][[2]])), 'CdTe_max' = sapply(1:76, function(x) max(panelrmse[[1]][[x]][[2]])),
                         'CIS_min' = sapply(1:76, function(x) min(panelrmse[[2]][[x]][[2]])), 'CIS_mean' = sapply(1:76, function(x) mean(panelrmse[[2]][[x]][[2]])), 'CIS_max' = sapply(1:76, function(x) max(panelrmse[[2]][[x]][[2]])),
                         'crystSI_min' = sapply(1:76, function(x) min(panelrmse[[3]][[x]][[2]])), 'crystSI_mean' = sapply(1:76, function(x) mean(panelrmse[[3]][[x]][[2]])), 'crystSI_max' = sapply(1:76, function(x) max(panelrmse[[3]][[x]][[2]])),
                         'Region' = latlon$region, 'Municipality' = latlon$name)

allnrmse_d$bestpanel_mean <- case_when(abs(allnrmse_d$CdTe_mean) < abs(allnrmse_d$CIS_mean) & abs(allnrmse_d$CdTe_mean) < abs(allnrmse_d$crystSI_mean) ~ 'CdTe',
                                       abs(allnrmse_d$CIS_mean) < abs(allnrmse_d$CdTe_mean) & abs(allnrmse_d$CIS_mean) < abs(allnrmse_d$crystSI_mean) ~ 'CIS',
                                       abs(allnrmse_d$crystSI_mean) < abs(allnrmse_d$CdTe_mean) & abs(allnrmse_d$crystSI_mean) < abs(allnrmse_d$CIS_mean) ~ 'CrystSi')

allnrmse_d$bestpanel_min <- case_when(abs(allnrmse_d$CdTe_min) < abs(allnrmse_d$CIS_min) & abs(allnrmse_d$CdTe_min) < abs(allnrmse_d$crystSI_min) ~ 'CdTe',
                                      abs(allnrmse_d$CIS_min) < abs(allnrmse_d$CdTe_min) & abs(allnrmse_d$CIS_min) < abs(allnrmse_d$crystSI_min) ~ 'CIS',
                                      abs(allnrmse_d$crystSI_min) < abs(allnrmse_d$CdTe_min) & abs(allnrmse_d$crystSI_min) < abs(allnrmse_d$CIS_min) ~ 'CrystSi')

allnrmse_d$bestpanel_max <- case_when(abs(allnrmse_d$CdTe_max) < abs(allnrmse_d$CIS_max) & abs(allnrmse_d$CdTe_max) < abs(allnrmse_d$crystSI_max) ~ 'CdTe',
                                      abs(allnrmse_d$CIS_max) < abs(allnrmse_d$CdTe_max) & abs(allnrmse_d$CIS_max) < abs(allnrmse_d$crystSI_max) ~ 'CIS',
                                      abs(allnrmse_d$crystSI_max) < abs(allnrmse_d$CdTe_max) & abs(allnrmse_d$crystSI_max) < abs(allnrmse_d$CIS_max) ~ 'CrystSi')

# t-tests
t.test(allnrmse_r$crystSI_mean, allnrmse_r$CIS_mean)
t.test(allnrmse_d$crystSI_mean, allnrmse_d$CIS_mean)

# best panel counts
nrow(allnrmse_r[allnrmse_r$bestpanel_min == 'CrystSi',])
nrow(allnrmse_r[allnrmse_r$bestpanel_max == 'CrystSi',])
nrow(allnrmse_r[allnrmse_r$bestpanel_mean == 'CrystSi',])

nrow(allnrmse_r[allnrmse_r$bestpanel_min == 'CIS',])
nrow(allnrmse_r[allnrmse_r$bestpanel_max == 'CIS',])
nrow(allnrmse_r[allnrmse_r$bestpanel_mean == 'CIS',])

nrow(allnrmse_r[allnrmse_r$bestpanel_min == 'CdTe',])
nrow(allnrmse_r[allnrmse_r$bestpanel_max == 'CdTe',])
nrow(allnrmse_r[allnrmse_r$bestpanel_mean == 'CdTe',])

nrow(allnrmse_d[allnrmse_d$bestpanel_min == 'CrystSi',])
nrow(allnrmse_d[allnrmse_d$bestpanel_max == 'CrystSi',])
nrow(allnrmse_d[allnrmse_d$bestpanel_mean == 'CrystSi',])

nrow(allnrmse_d[allnrmse_d$bestpanel_min == 'CIS',])
nrow(allnrmse_d[allnrmse_d$bestpanel_max == 'CIS',])
nrow(allnrmse_d[allnrmse_d$bestpanel_mean == 'CIS',])

nrow(allnrmse_d[allnrmse_d$bestpanel_min == 'CdTe',])
nrow(allnrmse_d[allnrmse_d$bestpanel_max == 'CdTe',])
nrow(allnrmse_d[allnrmse_d$bestpanel_mean == 'CdTe',])


# best panel by deviation counts
alldeviations <- list()

for (i in 1:3) {
  # initialize variable
  municdeviations <- list()
  
  # calculate deviation for each municipality and season
  for (j in 1:76) {
    # extract actual and predicted data
    predicted <- panelyhat[[i]][[j]]
    actual <- paneltest[[i]][[j]]
    
    # initialize variable
    deviation <- list()
    
    for (s in 1:2) {
      deviation[[s]] <- predicted[[s]] - actual[[s]]
    }
    municdeviations[[j]] <- deviation
  }
  alldeviations[[i]] <- municdeviations
}

latlon <- read_excel('/Users/rqo5125/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/2025_26/papers/PRsolarpower/location_coordinates_allmunicipalties.xlsx')

alldeviations2_r <- data.frame('CdTe_min' = sapply(1:76, function(x) min(alldeviations[[1]][[x]][[1]])), 'CdTe_mean' = sapply(1:76, function(x) mean(alldeviations[[1]][[x]][[1]])), 'CdTe_max' = sapply(1:76, function(x) max(alldeviations[[1]][[x]][[1]])),
                               'CIS_min' = sapply(1:76, function(x) min(alldeviations[[2]][[x]][[1]])), 'CIS_mean' = sapply(1:76, function(x) mean(alldeviations[[2]][[x]][[1]])), 'CIS_max' = sapply(1:76, function(x) max(alldeviations[[2]][[x]][[1]])),
                               'crystSI_min' = sapply(1:76, function(x) min(alldeviations[[3]][[x]][[1]])), 'crystSI_mean' = sapply(1:76, function(x) mean(alldeviations[[3]][[x]][[1]])), 'crystSI_max' = sapply(1:76, function(x) max(alldeviations[[3]][[x]][[1]])),
                               'Region' = latlon$region, 'Municipality' = latlon$name)

alldeviations2_r$bestpanel_mean <- case_when(abs(alldeviations2_r$CdTe_mean) < abs(alldeviations2_r$CIS_mean) & abs(alldeviations2_r$CdTe_mean) < abs(alldeviations2_r$crystSI_mean) ~ 'CdTe',
                                             abs(alldeviations2_r$CIS_mean) < abs(alldeviations2_r$CdTe_mean) & abs(alldeviations2_r$CIS_mean) < abs(alldeviations2_r$crystSI_mean) ~ 'CIS',
                                             abs(alldeviations2_r$crystSI_mean) < abs(alldeviations2_r$CdTe_mean) & abs(alldeviations2_r$crystSI_mean) < abs(alldeviations2_r$CIS_mean) ~ 'CrystSi')

alldeviations2_r$bestpanel_min <- case_when(abs(alldeviations2_r$CdTe_min) < abs(alldeviations2_r$CIS_min) & abs(alldeviations2_r$CdTe_min) < abs(alldeviations2_r$crystSI_min) ~ 'CdTe',
                                            abs(alldeviations2_r$CIS_min) < abs(alldeviations2_r$CdTe_min) & abs(alldeviations2_r$CIS_min) < abs(alldeviations2_r$crystSI_min) ~ 'CIS',
                                            abs(alldeviations2_r$crystSI_min) < abs(alldeviations2_r$CdTe_min) & abs(alldeviations2_r$crystSI_min) < abs(alldeviations2_r$CIS_min) ~ 'CrystSi')

alldeviations2_r$bestpanel_max <- case_when(abs(alldeviations2_r$CdTe_max) < abs(alldeviations2_r$CIS_max) & abs(alldeviations2_r$CdTe_max) < abs(alldeviations2_r$crystSI_max) ~ 'CdTe',
                                            abs(alldeviations2_r$CIS_max) < abs(alldeviations2_r$CdTe_max) & abs(alldeviations2_r$CIS_max) < abs(alldeviations2_r$crystSI_max) ~ 'CIS',
                                            abs(alldeviations2_r$crystSI_max) < abs(alldeviations2_r$CdTe_max) & abs(alldeviations2_r$crystSI_max) < abs(alldeviations2_r$CIS_max) ~ 'CrystSi')

alldeviations2_d <- data.frame('CdTe_min' = sapply(1:76, function(x) min(alldeviations[[1]][[x]][[2]])), 'CdTe_mean' = sapply(1:76, function(x) mean(alldeviations[[1]][[x]][[2]])), 'CdTe_max' = sapply(1:76, function(x) max(alldeviations[[1]][[x]][[2]])),
                               'CIS_min' = sapply(1:76, function(x) min(alldeviations[[2]][[x]][[2]])), 'CIS_mean' = sapply(1:76, function(x) mean(alldeviations[[2]][[x]][[2]])), 'CIS_max' = sapply(1:76, function(x) max(alldeviations[[2]][[x]][[2]])),
                               'crystSI_min' = sapply(1:76, function(x) min(alldeviations[[3]][[x]][[2]])), 'crystSI_mean' = sapply(1:76, function(x) mean(alldeviations[[3]][[x]][[2]])), 'crystSI_max' = sapply(1:76, function(x) max(alldeviations[[3]][[x]][[2]])),
                               'Region' = latlon$region, 'Municipality' = latlon$name)

alldeviations2_d$bestpanel_mean <- case_when(abs(alldeviations2_d$CdTe_mean) < abs(alldeviations2_d$CIS_mean) & abs(alldeviations2_d$CdTe_mean) < abs(alldeviations2_d$crystSI_mean) ~ 'CdTe',
                                             abs(alldeviations2_d$CIS_mean) < abs(alldeviations2_d$CdTe_mean) & abs(alldeviations2_d$CIS_mean) < abs(alldeviations2_d$crystSI_mean) ~ 'CIS',
                                             abs(alldeviations2_d$crystSI_mean) < abs(alldeviations2_d$CdTe_mean) & abs(alldeviations2_d$crystSI_mean) < abs(alldeviations2_d$CIS_mean) ~ 'CrystSi')

alldeviations2_d$bestpanel_min <- case_when(abs(alldeviations2_d$CdTe_min) < abs(alldeviations2_d$CIS_min) & abs(alldeviations2_d$CdTe_min) < abs(alldeviations2_d$crystSI_min) ~ 'CdTe',
                                            abs(alldeviations2_d$CIS_min) < abs(alldeviations2_d$CdTe_min) & abs(alldeviations2_d$CIS_min) < abs(alldeviations2_d$crystSI_min) ~ 'CIS',
                                            abs(alldeviations2_d$crystSI_min) < abs(alldeviations2_d$CdTe_min) & abs(alldeviations2_d$crystSI_min) < abs(alldeviations2_d$CIS_min) ~ 'CrystSi')

alldeviations2_d$bestpanel_max <- case_when(abs(alldeviations2_d$CdTe_max) < abs(alldeviations2_d$CIS_max) & abs(alldeviations2_d$CdTe_max) < abs(alldeviations2_d$crystSI_max) ~ 'CdTe',
                                            abs(alldeviations2_d$CIS_max) < abs(alldeviations2_d$CdTe_max) & abs(alldeviations2_d$CIS_max) < abs(alldeviations2_d$crystSI_max) ~ 'CIS',
                                            abs(alldeviations2_d$crystSI_max) < abs(alldeviations2_d$CdTe_max) & abs(alldeviations2_d$crystSI_max) < abs(alldeviations2_d$CIS_max) ~ 'CrystSi')


nrow(alldeviations2_r[alldeviations2_r$bestpanel_min == 'CrystSi',])
nrow(alldeviations2_r[alldeviations2_r$bestpanel_max == 'CrystSi',])
nrow(alldeviations2_r[alldeviations2_r$bestpanel_mean == 'CrystSi',])

nrow(alldeviations2_r[alldeviations2_r$bestpanel_min == 'CIS',])
nrow(alldeviations2_r[alldeviations2_r$bestpanel_max == 'CIS',])
nrow(alldeviations2_r[alldeviations2_r$bestpanel_mean == 'CIS',])

nrow(alldeviations2_r[alldeviations2_r$bestpanel_min == 'CdTe',])
nrow(alldeviations2_r[alldeviations2_r$bestpanel_max == 'CdTe',])
nrow(alldeviations2_r[alldeviations2_r$bestpanel_mean == 'CdTe',])

nrow(alldeviations2_d[alldeviations2_d$bestpanel_min == 'CrystSi',])
nrow(alldeviations2_d[alldeviations2_d$bestpanel_max == 'CrystSi',])
nrow(alldeviations2_d[alldeviations2_d$bestpanel_mean == 'CrystSi',])

nrow(alldeviations2_d[alldeviations2_d$bestpanel_min == 'CIS',])
nrow(alldeviations2_d[alldeviations2_d$bestpanel_max == 'CIS',])
nrow(alldeviations2_d[alldeviations2_d$bestpanel_mean == 'CIS',])

nrow(alldeviations2_d[alldeviations2_d$bestpanel_min == 'CdTe',])
nrow(alldeviations2_d[alldeviations2_d$bestpanel_max == 'CdTe',])
nrow(alldeviations2_d[alldeviations2_d$bestpanel_mean == 'CdTe',])

sum(alldeviations2_r$crystSI_mean)
sum(alldeviations2_d$crystSI_mean)

sum(alldeviations2_r$crystSI_max) * (7*30) + sum(alldeviations2_d$crystSI_max) * (5*30)

(sum(alldeviations2_r$crystSI_max) * (7*30) + sum(alldeviations2_d$crystSI_max) * (5*30)) / 1000

################ FIGURES AND TABLES ################

setwd(rdatadir)
load('modelresults.rdata')
load('processeddata.rdata')

# TABLES - NRMSE & R^2 (By Municipality)

# organize data
latlon <- read_excel('/Users/rqo5125/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/2025_26/papers/PRsolarpower/location_coordinates_allmunicipalties.xlsx')

allnrmse_r <- data.frame('CdTe_min' = sapply(1:76, function(x) min(panelrmse[[1]][[x]][[1]])), 'CdTe_mean' = sapply(1:76, function(x) mean(panelrmse[[1]][[x]][[1]])), 'CdTe_max' = sapply(1:76, function(x) max(panelrmse[[1]][[x]][[1]])),
                         'CIS_min' = sapply(1:76, function(x) min(panelrmse[[2]][[x]][[1]])), 'CIS_mean' = sapply(1:76, function(x) mean(panelrmse[[2]][[x]][[1]])), 'CIS_max' = sapply(1:76, function(x) max(panelrmse[[2]][[x]][[1]])),
                         'crystSI_min' = sapply(1:76, function(x) min(panelrmse[[3]][[x]][[1]])), 'crystSI_mean' = sapply(1:76, function(x) mean(panelrmse[[3]][[x]][[1]])), 'crystSI_max' = sapply(1:76, function(x) max(panelrmse[[3]][[x]][[1]])),
                         'Region' = latlon$region, 'Season' = rep('Rainy',76))
allnrmse_d <- data.frame('CdTe_min' = sapply(1:76, function(x) min(panelrmse[[1]][[x]][[2]])), 'CdTe_mean' = sapply(1:76, function(x) mean(panelrmse[[1]][[x]][[2]])), 'CdTe_max' = sapply(1:76, function(x) max(panelrmse[[1]][[x]][[2]])),
                         'CIS_min' = sapply(1:76, function(x) min(panelrmse[[2]][[x]][[2]])), 'CIS_mean' = sapply(1:76, function(x) mean(panelrmse[[2]][[x]][[2]])), 'CIS_max' = sapply(1:76, function(x) max(panelrmse[[2]][[x]][[2]])),
                         'crystSI_min' = sapply(1:76, function(x) min(panelrmse[[3]][[x]][[2]])), 'crystSI_mean' = sapply(1:76, function(x) mean(panelrmse[[3]][[x]][[2]])), 'crystSI_max' = sapply(1:76, function(x) max(panelrmse[[3]][[x]][[2]])),
                         'Region' = latlon$region, 'Season' = rep('Dry',76))
allnrmse <- rbind(allnrmse_r, allnrmse_d)

allrsq_r <- data.frame('CdTe_min' = sapply(1:76, function(x) min(panelrsq[[1]][[x]][[1]])), 'CdTe_mean' = sapply(1:76, function(x) mean(panelrsq[[1]][[x]][[1]])), 'CdTe_max' = sapply(1:76, function(x) max(panelrsq[[1]][[x]][[1]])),
                       'CIS_min' = sapply(1:76, function(x) min(panelrsq[[2]][[x]][[1]])), 'CIS_mean' = sapply(1:76, function(x) mean(panelrsq[[2]][[x]][[1]])), 'CIS_max' = sapply(1:76, function(x) max(panelrsq[[2]][[x]][[1]])),
                       'crystSI_min' = sapply(1:76, function(x) min(panelrsq[[3]][[x]][[1]])), 'crystSI_mean' = sapply(1:76, function(x) mean(panelrsq[[3]][[x]][[1]])), 'crystSI_max' = sapply(1:76, function(x) max(panelrsq[[3]][[x]][[1]])),
                       'Region' = latlon$region, 'Season' = rep('Rainy',76))
allrsq_d <- data.frame('CdTe_min' = sapply(1:76, function(x) min(panelrsq[[1]][[x]][[2]])), 'CdTe_mean' = sapply(1:76, function(x) mean(panelrsq[[1]][[x]][[2]])), 'CdTe_max' = sapply(1:76, function(x) max(panelrsq[[1]][[x]][[2]])),
                       'CIS_min' = sapply(1:76, function(x) min(panelrsq[[2]][[x]][[2]])), 'CIS_mean' = sapply(1:76, function(x) mean(panelrsq[[2]][[x]][[2]])), 'CIS_max' = sapply(1:76, function(x) max(panelrsq[[2]][[x]][[2]])),
                       'crystSI_min' = sapply(1:76, function(x) min(panelrsq[[3]][[x]][[2]])), 'crystSI_mean' = sapply(1:76, function(x) mean(panelrsq[[3]][[x]][[2]])), 'crystSI_max' = sapply(1:76, function(x) max(panelrsq[[3]][[x]][[2]])),
                       'Region' = latlon$region, 'Season' = rep('Dry',76))
allrsq <- rbind(allrsq_r, allrsq_d)

# write to csv
setwd(outputdir)
write.csv(allnrmse, 'nrmse_municipalities.csv')
write.csv(allrsq, 'rsq_municipalities.csv')

# TABLES - NRMSE & R^2 (By Region)

# organize data
regionnrmse <- aggregate(allnrmse[,1:9], list(allnrmse$Region, allnrmse$Season), FUN = mean)
regionrsq <- aggregate(allrsq[,1:9], list(allrsq$Region, allrsq$Season), FUN = mean)

# write to csv
setwd(outputdir)
write.csv(regionnrmse, 'nrmse_regions.csv')
write.csv(regionrsq, 'rsq_regions.csv')

# TABLES - NRMSE & R^2 (By Panel)

# organize data
panelnrmse <- aggregate(allnrmse[,1:9], list(allnrmse$Season), FUN = mean)
panelrsq <- aggregate(allrsq[,1:9], list(allrsq$Season), FUN = mean)

# write to csv
setwd(outputdir)
write.csv(panelnrmse, 'nrmse_panels.csv')
write.csv(panelrsq, 'rsq_panels.csv')

# FIGURES - Box Plot of NRMSE (By Panel)

# organize data
plotdata <- data.frame('mins' = c(panelnrmse[,2], panelnrmse[,5], panelnrmse[,8]), 
                       'means' = c(panelnrmse[,3], panelnrmse[,6], panelnrmse[,9]), 
                       'maxs' = c(panelnrmse[,4], panelnrmse[,7], panelnrmse[,10]), 
                       'panel' = c('CdTe', 'CdTe', 'CIS', 'CIS', 'CrystalSI', 'CrystalSI'),
                       'season' = c('Dry', 'Rainy', 'Dry', 'Rainy', 'Dry', 'Rainy'))

# create & store figure
setwd(figuredir)
pdf('nrmse_panels.pdf', width = 6, height = 6)
ggplot(plotdata, aes(x = panel, ymin = mins, ymax = maxs, lower = mins, upper = maxs, middle = means)) + 
  geom_boxplot(stat = 'identity', aes(fill = panel)) + theme_light() + guides(fill="none") +
  xlab('Panel') + ylab('NRMSE') + theme(legend.position = 'bottom', text = element_text(size = 20)) +
  scale_fill_manual(values = c('#66c2a5', '#fc8d62', '#8da0cb')) +
  facet_wrap(~season, nrow = 2)
dev.off()

# FIGURES - Box Plots of NRMSE (By Region)

# organize data
cdtedata <- data.frame('mins' = regionnrmse$CdTe_min, 'means' = regionnrmse$CdTe_mean, 
                       'maxs' = regionnrmse$CdTe_max, 'region' = regionnrmse$Group.1, 'season' = regionnrmse$Group.2)
cisdata <- data.frame('mins' = regionnrmse$CIS_min, 'means' = regionnrmse$CIS_mean, 
                       'maxs' = regionnrmse$CIS_max, 'region' = regionnrmse$Group.1, 'season' = regionnrmse$Group.2)
crystSIdata <- data.frame('mins' = regionnrmse$crystSI_min, 'means' = regionnrmse$crystSI_mean, 
                       'maxs' = regionnrmse$crystSI_max, 'region' = regionnrmse$Group.1, 'season' = regionnrmse$Group.2)

paneldata <- list(cdtedata, cisdata, crystSIdata)

# create & store figures
setwd(figuredir)
panelfigurenames <- c('nrmse_regions_CdTe.pdf', 'nrmse_regions_CIS.pdf', 'nrmse_regions_crystSI.pdf')

# loop through panels
for (i in 1:3) {
  # extract plot data
  plotdata <- paneldata[[i]]
  
  # create & store figure
  pdf(panelfigurenames[i], width = 6, height = 6)
  p <- ggplot(plotdata, aes(x = region, ymin = mins, ymax = maxs, lower = mins, upper = maxs, middle = means)) + 
        geom_boxplot(stat = 'identity', aes(fill = region)) + theme_light() + guides(fill = "none") +
        xlab('Region') + ylab('NRMSE') + theme(legend.position = 'bottom', text = element_text(size = 20)) +
    scale_fill_manual(name = 'Region', values = c('#33a02c', '#a6cee3', '#1f78b4', '#ffff99', '#b2df8a', '#fff200'),
                      labels = c('Central', 'East', 'Metro', 'North', 'South', 'West')) +
    scale_x_discrete(labels = c('Central', 'East', 'Metro', 'North', 'South', 'West')) +
    facet_wrap(~season, nrow = 2)
  print(p)
  dev.off()
}

# FIGURES - Box Plots of NRMSE (By Municipality)

# organize data
allnrmse$Municipality <- latlon$name

# organize data
cdtedata <- data.frame('mins' = allnrmse$CdTe_min, 'means' = allnrmse$CdTe_mean, 'maxs' = allnrmse$CdTe_max, 'region' = allnrmse$Region, 'municipality' = allnrmse$Municipality, 'season' = allnrmse$Season)
cisdata <- data.frame('mins' = allnrmse$CIS_min, 'means' = allnrmse$CIS_mean, 'maxs' = allnrmse$CIS_max, 'region' = allnrmse$Region, 'municipality' = allnrmse$Municipality, 'season' = allnrmse$Season)
crystSIdata <- data.frame('mins' = allnrmse$crystSI_min, 'means' = allnrmse$crystSI_mean, 'maxs' = allnrmse$crystSI_max, 'region' = allnrmse$Region, 'municipality' = allnrmse$Municipality, 'season' = allnrmse$Season)

paneldata <- list(cdtedata, cisdata, crystSIdata)
panelnames <- c('CdTe', 'CIS', 'CrystSI')
regionnames <- c('Central', 'East', 'Metro', 'North', 'South', 'West')

# create & store figures
setwd(figuredir)

# loop through panels
for (i in 1:3) {
  # extract panel of interest
  municipalitydata <- paneldata[[i]]
  
  # separate into regions
  centraldata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'central')], 'means' = municipalitydata$means[which(municipalitydata$region == 'central')], 
                            'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'central')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'central')],
                            'season' = municipalitydata$season[which(municipalitydata$region == 'central')])
  eastdata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'east')], 'means' = municipalitydata$means[which(municipalitydata$region == 'east')], 
                         'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'east')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'east')],
                         'season' = municipalitydata$season[which(municipalitydata$region == 'east')])
  metrodata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'metro')], 'means' = municipalitydata$means[which(municipalitydata$region == 'metro')], 
                          'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'metro')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'metro')],
                          'season' = municipalitydata$season[which(municipalitydata$region == 'metro')])
  northdata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'north')], 'means' = municipalitydata$means[which(municipalitydata$region == 'north')], 
                          'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'north')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'north')],
                          'season' = municipalitydata$season[which(municipalitydata$region == 'north')])
  southdata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'south')], 'means' = municipalitydata$means[which(municipalitydata$region == 'south')], 
                          'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'south')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'south')],
                          'season' = municipalitydata$season[which(municipalitydata$region == 'south')])
  westdata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'west')], 'means' = municipalitydata$means[which(municipalitydata$region == 'west')], 
                         'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'west')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'west')],
                         'season' = municipalitydata$season[which(municipalitydata$region == 'west')])
  
  regiondata <- list(centraldata, eastdata, metrodata, northdata, southdata, westdata)
  
  for (j in 1:6) {
    # extract plot data
    plotdata <- regiondata[[j]]
    
    # create & store figure
    pdf(paste('nrmse_municipalities_', regionnames[j], '_', panelnames[i], '.pdf', sep = ''), width = 6, height = 6)
    p <- ggplot(plotdata, aes(x = municipality, ymin = mins, ymax = maxs, lower = mins, upper = maxs, middle = means)) + 
      geom_boxplot(stat = 'identity', fill = 'darkgrey') + theme_light() + 
      xlab('Municipality') + ylab('NRMSE') + theme(legend.position = 'bottom', text = element_text(size = 20),
                                                   axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_wrap(~season, nrow = 2)
    print(p)
    dev.off()
  }
}

# FIGURES - Box Plot of yhat (By Panel)

allyhat_r <- data.frame('CdTe_min' = sapply(1:76, function(x) min(panelyhat[[1]][[x]][[1]])), 'CdTe_mean' = sapply(1:76, function(x) mean(panelyhat[[1]][[x]][[1]])), 'CdTe_max' = sapply(1:76, function(x) max(panelyhat[[1]][[x]][[1]])),
                        'CIS_min' = sapply(1:76, function(x) min(panelyhat[[2]][[x]][[1]])), 'CIS_mean' = sapply(1:76, function(x) mean(panelyhat[[2]][[x]][[1]])), 'CIS_max' = sapply(1:76, function(x) max(panelyhat[[2]][[x]][[1]])),
                        'crystSI_min' = sapply(1:76, function(x) min(panelyhat[[3]][[x]][[1]])), 'crystSI_mean' = sapply(1:76, function(x) mean(panelyhat[[3]][[x]][[1]])), 'crystSI_max' = sapply(1:76, function(x) max(panelyhat[[3]][[x]][[1]])),
                        'Region' = latlon$region, 'Municipality' = latlon$name, 'Season' = rep('Rainy', 76))
allyhat_d <- data.frame('CdTe_min' = sapply(1:76, function(x) min(panelyhat[[1]][[x]][[2]])), 'CdTe_mean' = sapply(1:76, function(x) mean(panelyhat[[1]][[x]][[2]])), 'CdTe_max' = sapply(1:76, function(x) max(panelyhat[[1]][[x]][[2]])),
                        'CIS_min' = sapply(1:76, function(x) min(panelyhat[[2]][[x]][[2]])), 'CIS_mean' = sapply(1:76, function(x) mean(panelyhat[[2]][[x]][[2]])), 'CIS_max' = sapply(1:76, function(x) max(panelyhat[[2]][[x]][[2]])),
                        'crystSI_min' = sapply(1:76, function(x) min(panelyhat[[3]][[x]][[2]])), 'crystSI_mean' = sapply(1:76, function(x) mean(panelyhat[[3]][[x]][[2]])), 'crystSI_max' = sapply(1:76, function(x) max(panelyhat[[3]][[x]][[2]])),
                        'Region' = latlon$region, 'Municipality' = latlon$name, 'Season' = rep('Dry', 76))
allyhat <- rbind(allyhat_r, allyhat_d)

regionyhat <- aggregate(allyhat[,1:9], list(allyhat$Region, allyhat$Season), FUN = mean)
panelyhat2 <- aggregate(allyhat[,1:9], list(allyhat$Season), FUN = mean)

# organize data
plotdata <- data.frame('mins' = c(panelyhat2[,2], panelyhat2[,5], panelyhat2[,8]), 
                       'means' = c(panelyhat2[,3], panelyhat2[,6], panelyhat2[,9]), 
                       'maxs' = c(panelyhat2[,4], panelyhat2[,7], panelyhat2[,10]), 
                       'panel' = c('CdTe', 'CdTe', 'CIS', 'CIS', 'CrystalSI', 'CrystalSI'),
                       'season' = c('Dry', 'Rainy', 'Dry', 'Rainy', 'Dry', 'Rainy'))

# create & store figure
setwd(figuredir)
pdf('yhat_panels.pdf', width = 6, height = 6)
ggplot(plotdata, aes(x = panel, ymin = mins, ymax = maxs, lower = mins, upper = maxs, middle = means)) + 
  geom_boxplot(stat = 'identity', aes(fill = panel)) + theme_light() + guides(fill = 'none') +
  xlab('Panel') + ylab('Predicted Daily Energy Generated (Wh)') + theme(legend.position = 'bottom', text = element_text(size = 20)) +
  scale_fill_manual(values = c('#66c2a5', '#fc8d62', '#8da0cb')) +
  facet_wrap(~season, nrow = 2)
dev.off()

# FIGURES - Box Plots of yhat (By Region)

# organize data
cdtedata <- data.frame('mins' = regionyhat$CdTe_min, 'means' = regionyhat$CdTe_mean, 
                       'maxs' = regionyhat$CdTe_max, 'region' = regionyhat$Group.1, 'season' = regionyhat$Group.2)
cisdata <- data.frame('mins' = regionyhat$CIS_min, 'means' = regionyhat$CIS_mean, 
                      'maxs' = regionyhat$CIS_max, 'region' = regionyhat$Group.1, 'season' = regionyhat$Group.2)
crystSIdata <- data.frame('mins' = regionyhat$crystSI_min, 'means' = regionyhat$crystSI_mean, 
                          'maxs' = regionyhat$crystSI_max, 'region' = regionyhat$Group.1, 'season' = regionyhat$Group.2)

paneldata <- list(cdtedata, cisdata, crystSIdata)

# create & store figures
setwd(figuredir)
panelfigurenames <- c('yhat_regions_CdTe.pdf', 'yhat_regions_CIS.pdf', 'yhat_regions_crystSI.pdf')

# loop through panels
for (i in 1:3) {
  # extract plot data
  plotdata <- paneldata[[i]]
  
  # create & store figure
  pdf(panelfigurenames[i], width = 6, height = 6)
  p <- ggplot(plotdata, aes(x = region, ymin = mins, ymax = maxs, lower = mins, upper = maxs, middle = means)) + 
    geom_boxplot(stat = 'identity', aes(fill = region)) + theme_light() + guides(fill = 'none') +
    xlab('Region') + ylab('Predicted Daily Energy Generated (Wh)') + theme(legend.position = 'bottom', text = element_text(size = 20)) +
    scale_fill_manual(name = 'Region', values = c('#33a02c', '#a6cee3', '#1f78b4', '#ffff99', '#b2df8a', '#fff200'),
                      labels = c('Central', 'East', 'Metro', 'North', 'South', 'West')) +
    scale_x_discrete(labels = c('Central', 'East', 'Metro', 'North', 'South', 'West')) +
    facet_wrap(~season, nrow = 2, scales = 'free')
  print(p)
  dev.off()
}

# FIGURES - Box Plots of yhat (By Municipality)

# organize data
cdtedata <- data.frame('mins' = allyhat$CdTe_min, 'means' = allyhat$CdTe_mean, 'maxs' = allyhat$CdTe_max, 'region' = allyhat$Region, 'municipality' = allyhat$Municipality, 'season' = allyhat$Season)
cisdata <- data.frame('mins' = allyhat$CIS_min, 'means' = allyhat$CIS_mean, 'maxs' = allyhat$CIS_max, 'region' = allyhat$Region, 'municipality' = allyhat$Municipality, 'season' = allyhat$Season)
crystSIdata <- data.frame('mins' = allyhat$crystSI_min, 'means' = allyhat$crystSI_mean, 'maxs' = allyhat$crystSI_max, 'region' = allyhat$Region, 'municipality' = allyhat$Municipality, 'season' = allyhat$Season)

paneldata <- list(cdtedata, cisdata, crystSIdata)
panelnames <- c('CdTe', 'CIS', 'CrystSI')
regionnames <- c('Central', 'East', 'Metro', 'North', 'South', 'West')

# create & store figures
setwd(figuredir)

# loop through panels
for (i in 1:3) {
  # extract panel of interest
  municipalitydata <- paneldata[[i]]
  
  # separate into regions
  centraldata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'central')], 'means' = municipalitydata$means[which(municipalitydata$region == 'central')], 
                            'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'central')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'central')],
                            'season' = municipalitydata$season[which(municipalitydata$region == 'central')])
  eastdata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'east')], 'means' = municipalitydata$means[which(municipalitydata$region == 'east')], 
                         'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'east')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'east')],
                         'season' = municipalitydata$season[which(municipalitydata$region == 'east')])
  metrodata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'metro')], 'means' = municipalitydata$means[which(municipalitydata$region == 'metro')], 
                          'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'metro')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'metro')],
                          'season' = municipalitydata$season[which(municipalitydata$region == 'metro')])
  northdata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'north')], 'means' = municipalitydata$means[which(municipalitydata$region == 'north')], 
                          'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'north')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'north')],
                          'season' = municipalitydata$season[which(municipalitydata$region == 'north')])
  southdata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'south')], 'means' = municipalitydata$means[which(municipalitydata$region == 'south')], 
                          'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'south')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'south')],
                          'season' = municipalitydata$season[which(municipalitydata$region == 'south')])
  westdata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'west')], 'means' = municipalitydata$means[which(municipalitydata$region == 'west')], 
                         'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'west')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'west')],
                         'season' = municipalitydata$season[which(municipalitydata$region == 'west')])
  
  regiondata <- list(centraldata, eastdata, metrodata, northdata, southdata, westdata)
  
  for (j in 1:6) {
    # extract plot data
    plotdata <- regiondata[[j]]
    
    # create & store figure
    pdf(paste('yhat_municipalities_', regionnames[j], '_', panelnames[i], '.pdf', sep = ''), width = 6, height = 6)
    p <- ggplot(plotdata, aes(x = municipality, ymin = mins, ymax = maxs, lower = mins, upper = maxs, middle = means)) + 
      geom_boxplot(stat = 'identity', fill = 'darkgrey') + theme_light() + 
      xlab('Municipality') + ylab('Predicted Daily Energy Generated (Wh)') + 
      theme(legend.position = 'bottom', text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_wrap(~season, nrow = 2, scales = 'free')
    print(p)
    dev.off()
  }
}

# FIGURES - Box Plot of ytest (By Panel)

allytest_r <- data.frame('CdTe_min' = sapply(1:76, function(x) min(paneltest[[1]][[x]][[1]])), 'CdTe_mean' = sapply(1:76, function(x) mean(paneltest[[1]][[x]][[1]])), 'CdTe_max' = sapply(1:76, function(x) max(paneltest[[1]][[x]][[1]])),
                        'CIS_min' = sapply(1:76, function(x) min(paneltest[[2]][[x]][[1]])), 'CIS_mean' = sapply(1:76, function(x) mean(paneltest[[2]][[x]][[1]])), 'CIS_max' = sapply(1:76, function(x) max(paneltest[[2]][[x]][[1]])),
                        'crystSI_min' = sapply(1:76, function(x) min(paneltest[[3]][[x]][[1]])), 'crystSI_mean' = sapply(1:76, function(x) mean(paneltest[[3]][[x]][[1]])), 'crystSI_max' = sapply(1:76, function(x) max(paneltest[[3]][[x]][[1]])),
                        'Region' = latlon$region, 'Municipality' = latlon$name, 'Season' = rep('Rainy', 76))
allytest_d <- data.frame('CdTe_min' = sapply(1:76, function(x) min(paneltest[[1]][[x]][[2]])), 'CdTe_mean' = sapply(1:76, function(x) mean(paneltest[[1]][[x]][[2]])), 'CdTe_max' = sapply(1:76, function(x) max(paneltest[[1]][[x]][[2]])),
                        'CIS_min' = sapply(1:76, function(x) min(paneltest[[2]][[x]][[2]])), 'CIS_mean' = sapply(1:76, function(x) mean(paneltest[[2]][[x]][[2]])), 'CIS_max' = sapply(1:76, function(x) max(paneltest[[2]][[x]][[2]])),
                        'crystSI_min' = sapply(1:76, function(x) min(paneltest[[3]][[x]][[2]])), 'crystSI_mean' = sapply(1:76, function(x) mean(paneltest[[3]][[x]][[2]])), 'crystSI_max' = sapply(1:76, function(x) max(paneltest[[3]][[x]][[2]])),
                        'Region' = latlon$region, 'Municipality' = latlon$name, 'Season' = rep('Dry', 76))
allytest <- rbind(allytest_r, allytest_d)

regionytest <- aggregate(allytest[,1:9], list(allytest$Region, allytest$Season), FUN = mean)
paneltest2 <- aggregate(allytest[,1:9], list(allytest$Season), FUN = mean)

# organize data
plotdata <- data.frame('mins' = c(paneltest2[,2], paneltest2[,5], paneltest2[,8]), 
                       'means' = c(paneltest2[,3], paneltest2[,6], paneltest2[,9]), 
                       'maxs' = c(paneltest2[,4], paneltest2[,7], paneltest2[,10]), 
                       'panel' = c('CdTe', 'CdTe', 'CIS', 'CIS', 'CrystalSI', 'CrystalSI'),
                       'season' = c('Dry', 'Rainy', 'Dry', 'Rainy', 'Dry', 'Rainy'))

# create & store figure
setwd(figuredir)
pdf('ytest_panels.pdf', width = 6, height = 6)
ggplot(plotdata, aes(x = panel, ymin = mins, ymax = maxs, lower = mins, upper = maxs, middle = means)) + 
  geom_boxplot(stat = 'identity', aes(fill= panel)) + theme_light() + guides(fill = 'none') +
  xlab('Panel') + ylab('Actual Daily Energy Generated (Wh)') + theme(legend.position = 'bottom', text = element_text(size = 20))  +
  scale_fill_manual(values = c('#66c2a5', '#fc8d62', '#8da0cb')) +
  facet_wrap(~season, nrow = 2)
dev.off()

# FIGURES - Box Plots of ytest (By Region)

# organize data
cdtedata <- data.frame('mins' = regionytest$CdTe_min, 'means' = regionytest$CdTe_mean, 
                       'maxs' = regionytest$CdTe_max, 'region' = regionytest$Group.1, 'season' = regionytest$Group.2)
cisdata <- data.frame('mins' = regionytest$CIS_min, 'means' = regionytest$CIS_mean, 
                      'maxs' = regionytest$CIS_max, 'region' = regionytest$Group.1, 'season' = regionytest$Group.2)
crystSIdata <- data.frame('mins' = regionytest$crystSI_min, 'means' = regionytest$crystSI_mean, 
                          'maxs' = regionytest$crystSI_max, 'region' = regionytest$Group.1, 'season' = regionytest$Group.2)

paneldata <- list(cdtedata, cisdata, crystSIdata)

# create & store figures
setwd(figuredir)
panelfigurenames <- c('ytest_regions_CdTe.pdf', 'ytest_regions_CIS.pdf', 'ytest_regions_crystSI.pdf')

# loop through panels
for (i in 1:3) {
  # extract plot data
  plotdata <- paneldata[[i]]
  
  # create & store figure
  pdf(panelfigurenames[i], width = 6, height = 6)
  p <- ggplot(plotdata, aes(x = region, ymin = mins, ymax = maxs, lower = mins, upper = maxs, middle = means)) + 
    geom_boxplot(stat = 'identity', aes(fill = region)) + theme_light() + guides(fill = 'none') +
    xlab('Region') + ylab('Actual Daily Energy Generated (Wh)') + theme(legend.position = 'bottom', text = element_text(size = 20)) +
    scale_fill_manual(name = 'Region', values = c('#33a02c', '#a6cee3', '#1f78b4', '#ffff99', '#b2df8a', '#fff200'),
                      labels = c('Central', 'East', 'Metro', 'North', 'South', 'West')) +
    scale_x_discrete(labels = c('Central', 'East', 'Metro', 'North', 'South', 'West')) +
    facet_wrap(~season, nrow = 2, scales = 'free')
  print(p)
  dev.off()
}

# FIGURES - Box Plots of ytest (By Municipality)

# organize data
cdtedata <- data.frame('mins' = allytest$CdTe_min, 'means' = allytest$CdTe_mean, 'maxs' = allytest$CdTe_max, 'region' = allytest$Region, 'municipality' = allytest$Municipality, 'season' = allytest$Season)
cisdata <- data.frame('mins' = allytest$CIS_min, 'means' = allytest$CIS_mean, 'maxs' = allytest$CIS_max, 'region' = allytest$Region, 'municipality' = allytest$Municipality, 'season' = allytest$Season)
crystSIdata <- data.frame('mins' = allytest$crystSI_min, 'means' = allytest$crystSI_mean, 'maxs' = allytest$crystSI_max, 'region' = allytest$Region, 'municipality' = allytest$Municipality, 'season' = allytest$Season)

paneldata <- list(cdtedata, cisdata, crystSIdata)
panelnames <- c('CdTe', 'CIS', 'CrystSI')
regionnames <- c('Central', 'East', 'Metro', 'North', 'South', 'West')

# create & store figures
setwd(figuredir)

# loop through panels
for (i in 1:3) {
  # extract panel of interest
  municipalitydata <- paneldata[[i]]
  
  # separate into regions
  centraldata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'central')], 'means' = municipalitydata$means[which(municipalitydata$region == 'central')], 
                            'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'central')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'central')],
                            'season' = municipalitydata$season[which(municipalitydata$region == 'central')])
  eastdata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'east')], 'means' = municipalitydata$means[which(municipalitydata$region == 'east')], 
                         'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'east')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'east')],
                         'season' = municipalitydata$season[which(municipalitydata$region == 'east')])
  metrodata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'metro')], 'means' = municipalitydata$means[which(municipalitydata$region == 'metro')], 
                          'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'metro')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'metro')],
                          'season' = municipalitydata$season[which(municipalitydata$region == 'metro')])
  northdata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'north')], 'means' = municipalitydata$means[which(municipalitydata$region == 'north')], 
                          'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'north')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'north')],
                          'season' = municipalitydata$season[which(municipalitydata$region == 'north')])
  southdata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'south')], 'means' = municipalitydata$means[which(municipalitydata$region == 'south')], 
                          'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'south')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'south')],
                          'season' = municipalitydata$season[which(municipalitydata$region == 'south')])
  westdata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'west')], 'means' = municipalitydata$means[which(municipalitydata$region == 'west')], 
                         'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'west')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'west')],
                         'season' = municipalitydata$season[which(municipalitydata$region == 'west')])
  
  regiondata <- list(centraldata, eastdata, metrodata, northdata, southdata, westdata)
  
  for (j in 1:6) {
    # extract plot data
    plotdata <- regiondata[[j]]
    
    # create & store figure
    pdf(paste('ytest_municipalities_', regionnames[j], '_', panelnames[i], '.pdf', sep = ''), width = 6, height = 6)
    p <- ggplot(plotdata, aes(x = municipality, ymin = mins, ymax = maxs, lower = mins, upper = maxs, middle = means)) + 
      geom_boxplot(stat = 'identity', fill = 'darkgrey') + theme_light() + 
      xlab('Municipality') + ylab('Actual Daily Energy Generated (Wh)') + 
      theme(legend.position = 'bottom', text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_wrap(~season, nrow = 2, scales = 'free')
    print(p)
    dev.off()
  }
}

# FIGURE - Study Area Map

# organize data
puertorico <- st_read('/Users/rqo5125/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/data/shapefiles/puerto-rico-municipality-boundaries/puerto-rico-municipality-boundaries.shp')
names(puertorico)[2] <- 'name'
puertorico <- puertorico[which(puertorico$name != 'Vieques'),]

plotdata <- merge(puertorico, latlon, by = 'name', all = TRUE)

# create and store figure
setwd(figuredir)
pdf('studyarea.pdf', width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata, color = "#808588", aes(fill = region)) + annotation_scale(pad_y = unit(-0.003, "cm"), pad_x = unit(12.5, "cm")) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), panel.grid = element_blank(), panel.border = element_blank(), legend.position="bottom") +
  scale_fill_manual(name = 'Region', values = c('#33a02c', '#a6cee3', '#1f78b4', '#ffff99', '#b2df8a', '#fff200'),
                    labels = c('Central', 'East', 'Metro', 'North', 'South', 'West')) 
dev.off() 

# FIGURE - STUDY AREA WITH IDS

setwd(figuredir)
pdf('studyarea_ids.pdf', width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata, color = "#808588", aes(fill = region)) + annotation_scale(pad_y = unit(-0.003, "cm"), pad_x = unit(12.5, "cm")) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), panel.grid = element_blank(), panel.border = element_blank(), legend.position="bottom") +
  scale_fill_manual(name = 'Region', values = c('#33a02c', '#a6cee3', '#1f78b4', '#ffff99', '#b2df8a', '#fff200'),
                    labels = c('Central', 'East', 'Metro', 'North', 'South', 'West')) +
  geom_sf_text(data = plotdata, aes(label = OBJECTID), color = 'black')
dev.off() 

# write csv file for info on municipalities
municinfo <- as.data.frame(plotdata)
write.csv(municinfo[,-10], file = 'municipalityinfo.csv')

# FIGURE - Study Area Map + NRMSE by Panel

# organize data
names(allnrmse)[12] <- 'name'
plotdata <- merge(puertorico, allnrmse, by = 'name', all = TRUE) 

# create and store figures
setwd(figuredir)

pdf(paste('nrmsemap_CdTe_rainy.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$season == 'Rainy'),], color = "#808588", aes(fill = CdTe_mean)) +
    theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
          axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
          panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'NRMSE    ', low = '#f7fcfd', high = '#41ae76') + ggtitle('CdTe - Rainy Season') 
dev.off()

pdf(paste('nrmsemap_CdTe_dry.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$season == 'Dry'),], color = "#808588", aes(fill = CdTe_mean)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'NRMSE    ', low = '#f7fcfd', high = '#41ae76') + ggtitle('CdTe - Dry Season') 
dev.off()

pdf(paste('nrmsemap_CIS_rainy.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$season == 'Rainy'),], color = "#808588", aes(fill = CIS_mean)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'NRMSE    ', low = '#fff7ec', high = '#ef6548') + ggtitle('CIS - Rainy Season')
dev.off()

pdf(paste('nrmsemap_CIS_dry.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$season == 'Dry'),], color = "#808588", aes(fill = CIS_mean)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'NRMSE    ', low = '#fff7ec', high = '#ef6548') + ggtitle('CIS - Dry Season')
dev.off()

pdf(paste('nrmsemap_CrystSI_rainy.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$season == 'Rainy'),], color = "#808588", aes(fill = crystSI_mean)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'NRMSE    ', low = '#fcfbfd', high = '#807dba') + ggtitle('CrystSI - Rainy Season')
dev.off()

pdf(paste('nrmsemap_CrystSI_dry.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$season == 'Dry'),], color = "#808588", aes(fill = crystSI_mean)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'NRMSE    ', low = '#fcfbfd', high = '#807dba') + ggtitle('CrystSI - Dry Season')
dev.off()

# FIGURE - Study Area Map + Best Panel by NRMSE

# organize data
allnrmse$bestpanel_mean <- case_when(allnrmse$CdTe_mean < allnrmse$CIS_mean & allnrmse$CdTe_mean < allnrmse$crystSI_mean ~ 'CdTe',
                                     allnrmse$CIS_mean < allnrmse$CdTe_mean & allnrmse$CIS_mean < allnrmse$crystSI_mean ~ 'CIS',
                                     allnrmse$crystSI_mean < allnrmse$CdTe_mean & allnrmse$crystSI_mean < allnrmse$CIS_mean ~ 'CrystSi')

allnrmse$bestpanel_min <- case_when(allnrmse$CdTe_min < allnrmse$CIS_min & allnrmse$CdTe_min < allnrmse$crystSI_min ~ 'CdTe',
                                    allnrmse$CIS_min < allnrmse$CdTe_min & allnrmse$CIS_min < allnrmse$crystSI_min ~ 'CIS',
                                    allnrmse$crystSI_min < allnrmse$CdTe_min & allnrmse$crystSI_min < allnrmse$CIS_min ~ 'CrystSi')

allnrmse$bestpanel_max <- case_when(allnrmse$CdTe_max < allnrmse$CIS_max & allnrmse$CdTe_max < allnrmse$crystSI_max ~ 'CdTe',
                                    allnrmse$CIS_max < allnrmse$CdTe_max & allnrmse$CIS_max < allnrmse$crystSI_max ~ 'CIS',
                                    allnrmse$crystSI_max < allnrmse$CdTe_max & allnrmse$crystSI_max < allnrmse$CIS_max ~ 'CrystSi')

plotdata <- merge(puertorico, allnrmse, by = 'name', all = TRUE) 

# create and store figures
setwd(figuredir)
pdf(paste('nrmsemap_bestpanel_mean_rainy.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$season == 'Rainy'),], color = "#808588", aes(fill = bestpanel_mean)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
        scale_fill_manual(name = '', values = c('#807dba')) + ggtitle('Average - Rainy Season')
dev.off()

pdf(paste('nrmsemap_bestpanel_mean_dry.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$season == 'Dry'),], color = "#808588", aes(fill = bestpanel_mean)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_manual(name = '', values = c('#66c2a5', '#fc8d62', '#8da0cb')) + ggtitle('Average - Dry Season')
dev.off()

pdf(paste('nrmsemap_bestpanel_min_rainy.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$season == 'Rainy'),], color = "#808588", aes(fill = bestpanel_min)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_manual(name = '', values = c('#807dba')) + ggtitle('Minimum - Rainy Season')
dev.off()

pdf(paste('nrmsemap_bestpanel_min_dry.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$season == 'Dry'),], color = "#808588", aes(fill = bestpanel_min)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_manual(name = '', values = c('#66c2a5', '#fc8d62', '#8da0cb')) + ggtitle('Minimum - Dry Season')
dev.off()

pdf(paste('nrmsemap_bestpanel_max_rainy.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$season == 'Rainy'),], color = "#808588", aes(fill = bestpanel_max)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_manual(name = '', values = c('#fc8d62', '#8da0cb')) + ggtitle('Maximum - Rainy Season')
dev.off()

pdf(paste('nrmsemap_bestpanel_max_dry.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$season == 'Dry'),], color = "#808588", aes(fill = bestpanel_max)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_manual(name = '', values = c('#66c2a5', '#fc8d62', '#8da0cb')) + ggtitle('Maximum - Dry Season')
dev.off()

# FIGURE - Box Plot of Deviations (By Panel)

# calculate deviations
alldeviations <- list()

for (i in 1:3) {
  # initialize variable
  municdeviations <- list()
  
  # calculate deviation for each municipality and season
  for (j in 1:76) {
    # extract actual and predicted data
    predicted <- panelyhat[[i]][[j]]
    actual <- paneltest[[i]][[j]]
    
    # initialize variable
    deviation <- list()
    
    for (s in 1:2) {
      deviation[[s]] <- predicted[[s]] - actual[[s]]
    }
    municdeviations[[j]] <- deviation
  }
  alldeviations[[i]] <- municdeviations
}

alldeviations2_r <- data.frame('CdTe_min' = sapply(1:76, function(x) min(alldeviations[[1]][[x]][[1]])), 'CdTe_mean' = sapply(1:76, function(x) mean(alldeviations[[1]][[x]][[1]])), 'CdTe_max' = sapply(1:76, function(x) max(alldeviations[[1]][[x]][[1]])),
                               'CIS_min' = sapply(1:76, function(x) min(alldeviations[[2]][[x]][[1]])), 'CIS_mean' = sapply(1:76, function(x) mean(alldeviations[[2]][[x]][[1]])), 'CIS_max' = sapply(1:76, function(x) max(alldeviations[[2]][[x]][[1]])),
                               'crystSI_min' = sapply(1:76, function(x) min(alldeviations[[3]][[x]][[1]])), 'crystSI_mean' = sapply(1:76, function(x) mean(alldeviations[[3]][[x]][[1]])), 'crystSI_max' = sapply(1:76, function(x) max(alldeviations[[3]][[x]][[1]])),
                               'Region' = latlon$region, 'Municipality' = latlon$name, 'Season' = rep('Rainy', 76))
alldeviations2_d <- data.frame('CdTe_min' = sapply(1:76, function(x) min(alldeviations[[1]][[x]][[2]])), 'CdTe_mean' = sapply(1:76, function(x) mean(alldeviations[[1]][[x]][[2]])), 'CdTe_max' = sapply(1:76, function(x) max(alldeviations[[1]][[x]][[2]])),
                               'CIS_min' = sapply(1:76, function(x) min(alldeviations[[2]][[x]][[2]])), 'CIS_mean' = sapply(1:76, function(x) mean(alldeviations[[2]][[x]][[2]])), 'CIS_max' = sapply(1:76, function(x) max(alldeviations[[2]][[x]][[2]])),
                               'crystSI_min' = sapply(1:76, function(x) min(alldeviations[[3]][[x]][[2]])), 'crystSI_mean' = sapply(1:76, function(x) mean(alldeviations[[3]][[x]][[2]])), 'crystSI_max' = sapply(1:76, function(x) max(alldeviations[[3]][[x]][[2]])),
                               'Region' = latlon$region, 'Municipality' = latlon$name, 'Season' = rep('Dry', 76))

alldeviations2 <- rbind(alldeviations2_r, alldeviations2_d)
regiondeviations <- aggregate(alldeviations2[,1:9], list(alldeviations2$Region, alldeviations2$Season), FUN = mean)
panelydeviations <- aggregate(alldeviations2[,1:9], list(alldeviations2$Season), FUN = mean)

# organize data
plotdata <- data.frame('mins' = c(panelydeviations[,2], panelydeviations[,5], panelydeviations[,8]), 
                       'means' = c(panelydeviations[,3], panelydeviations[,6], panelydeviations[,9]), 
                       'maxs' = c(panelydeviations[,4], panelydeviations[,7], panelydeviations[,10]), 
                       'panel' = c('CdTe', 'CdTe', 'CIS', 'CIS', 'CrystalSI', 'CrystalSI'),
                       'season' = c('Dry', 'Rainy', 'Dry', 'Rainy', 'Dry', 'Rainy'))

# create & store figure
setwd(figuredir)
pdf('deviations_panels.pdf', width = 6, height = 6)
ggplot(plotdata, aes(x = panel, ymin = mins, ymax = maxs, lower = mins, upper = maxs, middle = means)) + 
  geom_boxplot(stat = 'identity', aes(fill= panel)) + theme_light() + guides(fill = 'none') +
  xlab('Panel') + ylab('Deviation from Actual Daily Generation (Wh)') + theme(legend.position = 'bottom', text = element_text(size = 20))  +
  scale_fill_manual(values = c('#66c2a5', '#fc8d62', '#8da0cb')) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  facet_wrap(~season, nrow = 2, scales = 'free')
dev.off()

# FIGURES - Box Plots of Deviations (By Region)

# organize data
cdtedata <- data.frame('mins' = regiondeviations$CdTe_min, 'means' = regiondeviations$CdTe_mean, 
                       'maxs' = regiondeviations$CdTe_max, 'region' = regiondeviations$Group.1, 'season' = regiondeviations$Group.2)
cisdata <- data.frame('mins' = regiondeviations$CIS_min, 'means' = regiondeviations$CIS_mean, 
                      'maxs' = regiondeviations$CIS_max, 'region' = regiondeviations$Group.1, 'season' = regiondeviations$Group.2)
crystSIdata <- data.frame('mins' = regiondeviations$crystSI_min, 'means' = regiondeviations$crystSI_mean, 
                          'maxs' = regiondeviations$crystSI_max, 'region' = regiondeviations$Group.1, 'season' = regiondeviations$Group.2)

paneldata <- list(cdtedata, cisdata, crystSIdata)

# create & store figures
setwd(figuredir)
panelfigurenames <- c('deviations_regions_CdTe.pdf', 'deviations_regions_CIS.pdf', 'deviations_regions_crystSI.pdf')

# loop through panels
for (i in 1:3) {
  # extract plot data
  plotdata <- paneldata[[i]]
  
  # create & store figure
  pdf(panelfigurenames[i], width = 6, height = 6)
  p <- ggplot(plotdata, aes(x = region, ymin = mins, ymax = maxs, lower = mins, upper = maxs, middle = means)) + 
    geom_boxplot(stat = 'identity', aes(fill = region)) + theme_light() + guides(fill = 'none') +
    xlab('Region') + ylab('Deviation from Actual Daily Generation (Wh)') + theme(legend.position = 'bottom', text = element_text(size = 20)) +
    scale_fill_manual(name = 'Region', values = c('#33a02c', '#a6cee3', '#1f78b4', '#ffff99', '#b2df8a', '#fff200'),
                      labels = c('Central', 'East', 'Metro', 'North', 'South', 'West')) +
    scale_x_discrete(labels = c('Central', 'East', 'Metro', 'North', 'South', 'West'))+
    geom_hline(yintercept = 0, linetype = 'dashed') +
    facet_wrap(~season, nrow = 2, scales = 'free')
  print(p)
  dev.off()
}

# FIGURES - Box Plots of Deviation (By Municipality)

# organize data
cdtedata <- data.frame('mins' = alldeviations2$CdTe_min, 'means' = alldeviations2$CdTe_mean, 'maxs' = alldeviations2$CdTe_max, 'region' = alldeviations2$Region, 'municipality' = alldeviations2$Municipality, 'season' = alldeviations2$Season)
cisdata <- data.frame('mins' = alldeviations2$CIS_min, 'means' = alldeviations2$CIS_mean, 'maxs' = alldeviations2$CIS_max, 'region' = alldeviations2$Region, 'municipality' = alldeviations2$Municipality, 'season' = alldeviations2$Season)
crystSIdata <- data.frame('mins' = alldeviations2$crystSI_min, 'means' = alldeviations2$crystSI_mean, 'maxs' = alldeviations2$crystSI_max, 'region' = alldeviations2$Region, 'municipality' = alldeviations2$Municipality, 'season' = alldeviations2$Season)

paneldata <- list(cdtedata, cisdata, crystSIdata)
panelnames <- c('CdTe', 'CIS', 'CrystSI')
regionnames <- c('Central', 'East', 'Metro', 'North', 'South', 'West')

# create & store figures
setwd(figuredir)

# loop through panels
for (i in 1:3) {
  # extract panel of interest
  municipalitydata <- paneldata[[i]]
  
  # separate into regions
  centraldata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'central')], 'means' = municipalitydata$means[which(municipalitydata$region == 'central')], 
                            'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'central')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'central')],
                            'season' = municipalitydata$season[which(municipalitydata$region == 'central')])
  eastdata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'east')], 'means' = municipalitydata$means[which(municipalitydata$region == 'east')], 
                         'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'east')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'east')],
                         'season' = municipalitydata$season[which(municipalitydata$region == 'east')])
  metrodata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'metro')], 'means' = municipalitydata$means[which(municipalitydata$region == 'metro')], 
                          'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'metro')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'metro')],
                          'season' = municipalitydata$season[which(municipalitydata$region == 'metro')])
  northdata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'north')], 'means' = municipalitydata$means[which(municipalitydata$region == 'north')], 
                          'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'north')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'north')],
                          'season' = municipalitydata$season[which(municipalitydata$region == 'north')])
  southdata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'south')], 'means' = municipalitydata$means[which(municipalitydata$region == 'south')], 
                          'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'south')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'south')],
                          'season' = municipalitydata$season[which(municipalitydata$region == 'south')])
  westdata <- data.frame('mins' = municipalitydata$mins[which(municipalitydata$region == 'west')], 'means' = municipalitydata$means[which(municipalitydata$region == 'west')], 
                         'maxs' = municipalitydata$maxs[which(municipalitydata$region == 'west')], 'municipality' = municipalitydata$municipality[which(municipalitydata$region == 'west')],
                         'season' = municipalitydata$season[which(municipalitydata$region == 'west')])
  
  regiondata <- list(centraldata, eastdata, metrodata, northdata, southdata, westdata)
  
  for (j in 1:6) {
    # extract plot data
    plotdata <- regiondata[[j]]
    
    # create & store figure
    pdf(paste('deviation_municipalities_', regionnames[j], '_', panelnames[i], '.pdf', sep = ''), width = 6, height = 7)
    p <- ggplot(plotdata, aes(x = municipality, ymin = mins, ymax = maxs, lower = mins, upper = maxs, middle = means)) + 
      geom_boxplot(stat = 'identity', fill = 'darkgrey') + theme_light() + 
      xlab('Municipality') + ylab('Deviation from Actual Daily Generation (Wh)') + geom_hline(yintercept = 0, linetype = 'dashed') +
      theme(legend.position = 'bottom', text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_wrap(~season, nrow = 2, scales = 'free')
    print(p)
    dev.off()
  }
}

# FIGURE - Study Area Map + Average Deviations by Panel

# organize data
names(alldeviations2)[11] <- 'name'
plotdata <- merge(puertorico, alldeviations2, by = 'name', all = TRUE) 

# create and store figures
setwd(figuredir)

pdf(paste('avg_deviationsmap_CdTe_rainy.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Rainy'),], color = "#808588", aes(fill = CdTe_mean)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'Deviation (Wh)   ', high = '#99d8c9', low = '#00441b') + ggtitle('CdTe - Rainy Season')
dev.off()

pdf(paste('avg_deviationsmap_CdTe_dry.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Dry'),], color = "#808588", aes(fill = CdTe_mean)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'Deviation (Wh)   ', high = '#99d8c9', low = '#00441b') + ggtitle('CdTe - Dry Season')
dev.off()

pdf(paste('avg_deviationsmap_CIS_rainy.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Rainy'),], color = "#808588", aes(fill = CIS_mean)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'Deviation (Wh)   ', high = '#fdbb84', low = '#7f0000') + ggtitle('CIS - Rainy Season')
dev.off()

pdf(paste('avg_deviationsmap_CIS_dry.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Dry'),], color = "#808588", aes(fill = CIS_mean)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'Deviation (Wh)   ', high = '#fdbb84', low = '#7f0000') + ggtitle('CIS - Dry Season')
dev.off()

pdf(paste('avg_deviationsmap_CrystSI_rainy.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Rainy'),], color = "#808588", aes(fill = crystSI_mean)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'Deviation (Wh)   ', high = '#bcbddc', low = '#3f007d') + ggtitle('CrystSI - Rainy Season')
dev.off()

pdf(paste('avg_deviationsmap_CrystSI_dry.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Dry'),], color = "#808588", aes(fill = crystSI_mean)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'Deviation (Wh)   ', high = '#bcbddc', low = '#3f007d') + ggtitle('CrystSI - Dry Season')
dev.off()

# FIGURE - Study Area Map + Maximum Deviations by Panel

# create and store figures
setwd(figuredir)

pdf(paste('max_deviationsmap_CdTe_rainy.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Rainy'),], color = "#808588", aes(fill = CdTe_max)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'Deviation (Wh)   ', low = '#99d8c9', high = '#00441b') + ggtitle('CdTe - Rainy Season')
dev.off()

pdf(paste('max_deviationsmap_CdTe_dry.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Dry'),], color = "#808588", aes(fill = CdTe_max)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'Deviation (Wh)   ', low = '#99d8c9', high = '#00441b') + ggtitle('CdTe - Dry Season')
dev.off()

pdf(paste('max_deviationsmap_CIS_rainy.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Rainy'),], color = "#808588", aes(fill = CIS_max)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'Deviation (Wh)   ', low = '#fdbb84', high = '#7f0000') + ggtitle('CIS - Rainy Season')
dev.off()

pdf(paste('max_deviationsmap_CIS_dry.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Dry'),], color = "#808588", aes(fill = CIS_max)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'Deviation (Wh)   ', low = '#fdbb84', high = '#7f0000') + ggtitle('CIS - Dry Season')
dev.off()

pdf(paste('max_deviationsmap_CrystSI_rainy.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Rainy'),], color = "#808588", aes(fill = crystSI_max)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'Deviation (Wh)   ', low = '#bcbddc', high = '#3f007d') + ggtitle('CrystSI - Rainy Season')
dev.off()

pdf(paste('max_deviationsmap_CrystSI_dry.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Dry'),], color = "#808588", aes(fill = crystSI_max)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'Deviation (Wh)   ', low = '#bcbddc', high = '#3f007d') + ggtitle('CrystSI - Dry Season')
dev.off()

# FIGURE - Study Area Map + Minimum Deviations by Panel

# create and store figures
setwd(figuredir)

pdf(paste('min_deviationsmap_CdTe_rainy.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Rainy'),], color = "#808588", aes(fill = CdTe_min)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'Deviation (Wh)   ', high = '#99d8c9', low = '#00441b') + ggtitle('CdTe - Rainy Season')
dev.off()

pdf(paste('min_deviationsmap_CdTe_dry.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Dry'),], color = "#808588", aes(fill = CdTe_min)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'Deviation (Wh)   ', high = '#99d8c9', low = '#00441b') + ggtitle('CdTe - Dry Season')
dev.off()

pdf(paste('min_deviationsmap_CIS_rainy.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Rainy'),], color = "#808588", aes(fill = CIS_min)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'Deviation (Wh)   ', high = '#fdbb84', low = '#7f0000') + ggtitle('CIS - Rainy Season')
dev.off()

pdf(paste('min_deviationsmap_CIS_dry.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Dry'),], color = "#808588", aes(fill = CIS_min)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'Deviation (Wh)   ', high = '#fdbb84', low = '#7f0000') + ggtitle('CIS - Dry Season')
dev.off()

pdf(paste('min_deviationsmap_CrystSI_rainy.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Rainy'),], color = "#808588", aes(fill = crystSI_min)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'Deviation (Wh)   ', high = '#bcbddc', low = '#3f007d') + ggtitle('CrystSI - Rainy Season')
dev.off()

pdf(paste('min_deviationsmap_CrystSI_dry.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Dry'),], color = "#808588", aes(fill = crystSI_min)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_gradient(name = 'Deviation (Wh)   ', high = '#bcbddc', low = '#3f007d') + ggtitle('CrystSI - Dry Season')
dev.off()

# FIGURE - Study Area Map + Best Panel by Deviation

# organize data
alldeviations2$bestpanel_mean <- case_when(abs(alldeviations2$CdTe_mean) < abs(alldeviations2$CIS_mean) & abs(alldeviations2$CdTe_mean) < abs(alldeviations2$crystSI_mean) ~ 'CdTe',
                                           abs(alldeviations2$CIS_mean) < abs(alldeviations2$CdTe_mean) & abs(alldeviations2$CIS_mean) < abs(alldeviations2$crystSI_mean) ~ 'CIS',
                                           abs(alldeviations2$crystSI_mean) < abs(alldeviations2$CdTe_mean) & abs(alldeviations2$crystSI_mean) < abs(alldeviations2$CIS_mean) ~ 'CrystSi')

alldeviations2$bestpanel_min <- case_when(abs(alldeviations2$CdTe_min) < abs(alldeviations2$CIS_min) & abs(alldeviations2$CdTe_min) < abs(alldeviations2$crystSI_min) ~ 'CdTe',
                                          abs(alldeviations2$CIS_min) < abs(alldeviations2$CdTe_min) & abs(alldeviations2$CIS_min) < abs(alldeviations2$crystSI_min) ~ 'CIS',
                                          abs(alldeviations2$crystSI_min) < abs(alldeviations2$CdTe_min) & abs(alldeviations2$crystSI_min) < abs(alldeviations2$CIS_min) ~ 'CrystSi')

alldeviations2$bestpanel_max <- case_when(abs(alldeviations2$CdTe_max) < abs(alldeviations2$CIS_max) & abs(alldeviations2$CdTe_max) < abs(alldeviations2$crystSI_max) ~ 'CdTe',
                                          abs(alldeviations2$CIS_max) < abs(alldeviations2$CdTe_max) & abs(alldeviations2$CIS_max) < abs(alldeviations2$crystSI_max) ~ 'CIS',
                                          abs(alldeviations2$crystSI_max) < abs(alldeviations2$CdTe_max) & abs(alldeviations2$crystSI_max) < abs(alldeviations2$CIS_max) ~ 'CrystSi')

plotdata <- merge(puertorico, alldeviations2, by = 'name', all = TRUE) 

# create and store figures
setwd(figuredir)
pdf(paste('deviationmap_bestpanel_mean_rainy.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Rainy'),], color = "#808588", aes(fill = bestpanel_mean)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_manual(name = '', values = c('#66c2a5', '#fc8d62', '#8da0cb')) + ggtitle('Average - Rainy Season')
dev.off()

pdf(paste('deviationmap_bestpanel_mean_dry.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Dry'),], color = "#808588", aes(fill = bestpanel_mean)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_manual(name = '', values = c('#66c2a5', '#fc8d62', '#8da0cb')) + ggtitle('Average - Dry Season')
dev.off()

pdf(paste('deviationmap_bestpanel_min_rainy.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Rainy'),], color = "#808588", aes(fill = bestpanel_min)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_manual(name = '', values = c('#fc8d62', '#8da0cb')) + ggtitle('Minimum - Rainy Season')
dev.off()

pdf(paste('deviationmap_bestpanel_min_dry.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Dry'),], color = "#808588", aes(fill = bestpanel_min)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_manual(name = '', values = c('#66c2a5', '#fc8d62', '#8da0cb')) + ggtitle('Minimum - Dry Season')
dev.off()

pdf(paste('deviationmap_bestpanel_max_rainy.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Rainy'),], color = "#808588", aes(fill = bestpanel_max)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_manual(name = '', values = c('#fc8d62', '#8da0cb')) + ggtitle('Maximum - Rainy Season')
dev.off()

pdf(paste('deviationmap_bestpanel_max_dry.pdf', sep = ''), width = 8, height = 5)
ggplot() + theme_light() + geom_sf(data = plotdata[which(plotdata$Season == 'Dry'),], color = "#808588", aes(fill = bestpanel_max)) +
  theme(text = element_text(size = 20), panel.background = element_blank(), axis.text = element_blank(), 
        axis.title = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
        panel.border = element_blank(), legend.position="bottom", legend.key.width = unit(1.5, "cm")) +
  scale_fill_manual(name = '', values = c('#fc8d62', '#8da0cb')) + ggtitle('Maximum - Dry Season')
dev.off()

# FIGURE - Correlation matrix across all municipalities 

cdtedata <- do.call(rbind, allpaneldata[[1]])
cisdata <- do.call(rbind, allpaneldata[[2]])
crystsidata <- do.call(rbind, allpaneldata[[3]])

alldata <- cbind(cdtedata, cisdata$power_watts, crystsidata$power_watts)
names(alldata)[13:15] <- c('CdTe_generation', 'CIS_generation', 'CrystSI_generation')

M <- cor(alldata[,-1], use = "complete.obs")

setwd(figuredir)
pdf('correlationplot.pdf', width = 11.5, height = 10)
corrplot.mixed(M, lower = 'color', upper = 'number', tl.pos = "lt", order = 'alphabet')
dev.off()

# FIGURE - VARIABLE IMPORTANCE

# extract variable importance

# create variable names
panels <- c('CdTe', 'CIS', 'CrystSi')
municipalities <- municipalitydata$municipality[1:76]
regions <- municipalitydata$region[1:76]
seasons <- c('Rainy', 'Dry')

# initialize variables
varimp <- data.frame()

# loop through panels
for (p in 1:3) {
  
  # loop through the municipality
  for (m in 1:76) {
    
    # loop through the seasons
    for (s in 1:2) {
      vitemp <- data.frame(panel = panels[p], municipality = municipalities[m], region = regions[m], season = seasons[s], t(panelmodel[[p]][[m]][[s]]$importance))
      
      varimp <- rbind(varimp, vitemp)
    }
  }
}

# get average node impurity by region
regionVI <- aggregate(varimp[,5:15], by = list(varimp$panel, varimp$region, varimp$season), FUN = mean)
regionVI_long <- melt(regionVI, id.vars = c('Group.1', 'Group.2', 'Group.3'))

setwd(figuredir)
pdf('variableimportance_regions.pdf', width = 7, height = 7)
ggplot(regionVI_long) +
  geom_point(aes(x = value/100000000, y = variable, color = Group.2)) +
  scale_color_manual(name = 'Region', values = c('#33a02c', '#a6cee3', '#1f78b4', '#ffff99', '#b2df8a', '#fff200'),
                     labels = c('Central', 'East', 'Metro', 'North', 'South', 'West')) +
  theme_light() + theme(legend.position = 'bottom') + xlab('Increased Node Impurity (10^8)') + ylab('Variable') +
  facet_wrap(~Group.3 + Group.1, nrow = 2, scales = "free_x") 
dev.off()

# get average node impurity by panel
panelVI <- aggregate(varimp[,5:15], by = list(varimp$panel, varimp$season), FUN = mean)
panelVI_long <- melt(panelVI, id.vars = c('Group.1', 'Group.2'))

pdf('variableimportance_panels.pdf', width = 5, height = 7)
ggplot(panelVI_long) +
  geom_point(aes(x = value/100000000, y = variable, color = Group.1)) +
  scale_color_manual(name = 'Region', values = c('#66c2a5', '#fc8d62', '#8da0cb')) +
  theme_light() + theme(legend.position = 'bottom') + xlab('Increased Node Impurity (10^8)') + ylab('Variable') +
  facet_wrap(~Group.2, nrow = 2, scales = "free_x") 
dev.off()

# FIGURE - PREDICTED VS ACTUAL
daycount <- c(1:length(panelyhat[[1]][[1]][[1]]), 1:length(panelyhat[[2]][[1]][[1]]), 1:length(panelyhat[[3]][[1]][[1]]))
season <- rep('Rainy', length(daycount))
panel <- c(rep('CdTe', length(panelyhat[[1]][[1]][[1]])), rep('CIS', length(panelyhat[[2]][[1]][[1]])), rep('crystSi', length(panelyhat[[3]][[1]][[1]])))
allyhat_r <- rbind(sapply(1:76, function(x) unlist(panelyhat[[1]][[x]][[1]])), sapply(1:76, function(x) unlist(panelyhat[[2]][[x]][[1]])),sapply(1:76, function(x) unlist(panelyhat[[3]][[x]][[1]])))
allyhat_r <- as.data.frame(cbind(daycount, season, panel, allyhat_r))
names(allyhat_r)[4:79] <- latlon$name
allyhat_r <- melt(allyhat_r, id.vars = c('daycount', 'season', 'panel'), variable.name = 'name', value.name = 'predicted')
allyhat_r <- merge(allyhat_r, latlon, on = c('name'))

daycount <- c(1:length(panelyhat[[1]][[1]][[2]]), 1:length(panelyhat[[2]][[1]][[2]]), 1:length(panelyhat[[3]][[1]][[2]]))
season <- rep('Dry', length(daycount))
panel <- c(rep('CdTe', length(panelyhat[[1]][[1]][[2]])), rep('CIS', length(panelyhat[[2]][[1]][[2]])), rep('crystSi', length(panelyhat[[3]][[1]][[2]])))
allyhat_d <- rbind(sapply(1:76, function(x) unlist(panelyhat[[1]][[x]][[2]])), sapply(1:76, function(x) unlist(panelyhat[[2]][[x]][[2]])),sapply(1:76, function(x) unlist(panelyhat[[3]][[x]][[2]])))
allyhat_d <- as.data.frame(cbind(daycount, season, panel, allyhat_d))
names(allyhat_d)[4:79] <- latlon$name
allyhat_d <- melt(allyhat_d, id.vars = c('daycount', 'season', 'panel'), variable.name = 'name', value.name = 'predicted')
allyhat_d <- merge(allyhat_d, latlon, on = c('name'))

daycount <- 1:length(paneltest[[1]][[1]][[1]])
season <- rep('Rainy', length(daycount))
allytest_r <- as.data.frame(cbind(daycount, season, sapply(1:76, function(x) unlist(paneltest[[1]][[x]][[1]]))))
names(allytest_r)[3:78] <- latlon$name
allytest_r <- melt(allytest_r, id.vars = c('daycount', 'season'), variable.name = 'name', value.name = 'actual')
allytest_r <- merge(allytest_r, latlon, on = c('name'))

daycount <- c(1:length(paneltest[[1]][[1]][[1]]), 1:length(paneltest[[2]][[1]][[1]]), 1:length(paneltest[[3]][[1]][[1]]))
season <- rep('Rainy', length(daycount))
panel <- c(rep('CdTe', length(paneltest[[1]][[1]][[1]])), rep('CIS', length(paneltest[[2]][[1]][[1]])), rep('crystSi', length(paneltest[[3]][[1]][[1]])))
alltest_r <- rbind(sapply(1:76, function(x) unlist(paneltest[[1]][[x]][[1]])), sapply(1:76, function(x) unlist(paneltest[[2]][[x]][[1]])),sapply(1:76, function(x) unlist(paneltest[[3]][[x]][[1]])))
alltest_r <- as.data.frame(cbind(daycount, season, panel, alltest_r))
names(alltest_r)[4:79] <- latlon$name
alltest_r <- melt(alltest_r, id.vars = c('daycount', 'season', 'panel'), variable.name = 'name', value.name = 'actual')
alltest_r <- merge(alltest_r, latlon, on = c('name'))

daycount <- c(1:length(paneltest[[1]][[1]][[2]]), 1:length(paneltest[[2]][[1]][[2]]), 1:length(paneltest[[3]][[1]][[2]]))
season <- rep('Dry', length(daycount))
panel <- c(rep('CdTe', length(paneltest[[1]][[1]][[2]])), rep('CIS', length(paneltest[[2]][[1]][[2]])), rep('crystSi', length(paneltest[[3]][[1]][[2]])))
alltest_d <- rbind(sapply(1:76, function(x) unlist(paneltest[[1]][[x]][[2]])), sapply(1:76, function(x) unlist(paneltest[[2]][[x]][[2]])),sapply(1:76, function(x) unlist(paneltest[[3]][[x]][[2]])))
alltest_d <- as.data.frame(cbind(daycount, season, panel, alltest_d))
names(alltest_d)[4:79] <- latlon$name
alltest_d <- melt(alltest_d, id.vars = c('daycount', 'season', 'panel'), variable.name = 'name', value.name = 'actual')
alltest_d <- merge(alltest_d, latlon, on = c('name'))

allyhat <- rbind(allyhat_r, allyhat_d)
allytest <- rbind(alltest_r, alltest_d)

plotdata <- merge(allyhat, allytest, on = c('name', 'daycount', 'season', 'region', 'panel'))
plotdata$predicted <- as.numeric(plotdata$predicted)
plotdata$actual <- as.numeric(plotdata$actual)

setwd(figuredir)
pdf('reliabilityplot.pdf', height = 7, width = 10)
ggplot(plotdata) +
  geom_point(aes(x = actual, y = predicted, color = region), alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  scale_color_manual(name = 'Region', values = c('#33a02c', '#a6cee3', '#1f78b4', '#ffff99', '#b2df8a', '#fff200')) +
  facet_wrap(~season + panel, nrow = 2, scales = 'free') +
  theme_light() + xlab('Actual Values') + ylab('Predicted Values') +
  theme(legend.position = 'bottom')
dev.off()

# OOB Plot

# extract OOB MSE

# create variable names
panels <- c('CdTe', 'CIS', 'CrystSi')
municipalities <- municipalitydata$municipality[1:76]
regions <- municipalitydata$region[1:76]
seasons <- c('Rainy', 'Dry')

# initialize variables
oobmse <- data.frame()

# loop through panels
for (p in 1:3) {
  
  # loop through the municipality
  for (m in 1:76) {
    
    # loop through the seasons
    for (s in 1:2) {
      oobtemp <- data.frame(panel = panels[p], municipality = municipalities[m], region = regions[m], season = seasons[s], tree = c(1:500), mse = panelmodel[[p]][[m]][[s]]$mse)
      
      oobmse <- rbind(oobmse, oobtemp)
    }
  }
}

# get average OOB by region
regionOOB <- aggregate(oobmse$mse, by = list(oobmse$panel, oobmse$region, oobmse$season, oobmse$tree), FUN = mean)

# plot data
setwd(figuredir)
pdf('oobrmse_regions.pdf', width = 7, height = 10)
ggplot(regionOOB) +
  geom_point(aes(x = Group.4, y = sqrt(x), color = Group.2)) +
  facet_wrap(~Group.3 + Group.1, scales = 'free') + xlab('Tree') + ylab('Out-of-Bag RMSE') +
  scale_color_manual(name = 'Region', values = c('#33a02c', '#a6cee3', '#1f78b4', '#ffff99', '#b2df8a', '#fff200')) +
  theme_light() + theme(legend.position = 'bottom')
dev.off()
