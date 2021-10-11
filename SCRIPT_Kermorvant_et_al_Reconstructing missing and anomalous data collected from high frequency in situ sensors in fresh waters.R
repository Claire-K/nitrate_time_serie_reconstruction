## SCRIPT PUBLICATION Reconstructing missing and anomalous data collected from high-frequency in-situ sensors in fresh waters (KERMORVANT et al.
######################################################################################"

# Load required packages
library(tidyverse)
library(mgcv)
library(mgcViz)
library(tsibble)
library(feasts)
library(lubridate)
library(forecast)
library(modelr)
library(Metrics)
library(zoo)

#####################################
# 1 - Load data from NEON database #
###################################

if(file.exists("waq_arik1.rda")) {
  waq_arik1 <- readRDS("waq_arik1.rda")
} else {waq_arik1 <- neonUtilities::loadByProduct(
  dpID = "DP1.20288.001",
  site = "ARIK",
  startdate = "2018-10",
  enddate = "2020-10",
  package = "expanded",
  token = Sys.getenv("NEON_TOKEN"),
  check.size = FALSE)

saveRDS(waq_arik1, "waq_arik1.rda")
}
# Take only water quality at the down station
water_quality_arik1 <- waq_arik1$waq_instantaneous %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  filter(horizontal_position == "102") %>%
  select(c(5,19,31,33,45,47,59,61,73,75,87,89,101,103,118)) %>%
  mutate(
    start_date_time = ymd_hms(start_date_time) - second(start_date_time)
  ) %>%
  rename(
    spec_cond = specific_conductance,
    label_spec_cond = specific_cond_final_qf,
    oxygen = dissolved_oxygen,
    label_oxygen = dissolved_oxygen_final_qf,
    oxygen_sat = dissolved_oxygen_saturation,
    label_oxygen_sat = dissolved_oxygen_sat_final_qf,
    ph = p_h,
    label_ph = p_h_final_qf,
    chloro = chlorophyll,
    label_chloro = chlorophyll_final_qf,
    label_turbidity = turbidity_final_qf,
    label_f_dom = f_dom_final_qf
  )


# download data of interest - Nitrate in Suface Water
if(file.exists("nsw_arik1.rda")) {
  nsw_arik1 <- readRDS("nsw_arik1.rda")
} else {
  nsw_arik1 <-  neonUtilities::loadByProduct(
    dpID="DP1.20033.001",
    site="ARIK",
    startdate = "2018-10",
    enddate = "2020-10",
    package="expanded",
    token = Sys.getenv("NEON_TOKEN"),
    check.size = FALSE
  )
  saveRDS(nsw_arik1, "nsw_arik1.rda")
}

# Extract nitrate at 15 minute intervals
nitrate_arik1 <- nsw_arik1$NSW_15_minute %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(
    start_date_time = ymd_hms(start_date_time) - second(start_date_time)
  ) %>%
  select(c(start_date_time, surf_water_nitrate_mean, final_qf)) %>%
  rename(
    nitrate_mean = surf_water_nitrate_mean,
    label_nitrate_mean = final_qf
  )

# download data of interest - Temperature Suface Water
if(file.exists("temp_arik1.rda")) {
  temp_arik1 <- readRDS("temp_arik1.rda")
} else {
  temp_arik1 <-  neonUtilities::loadByProduct(
    dpID="DP1.20053.001",
    site="ARIK",
    startdate = "2018-10",
    enddate = "2020-10",
    package="expanded",
    token = Sys.getenv("NEON_TOKEN"),
    check.size = FALSE
  )
  saveRDS(temp_arik1, "temp_arik1.rda")
}

# Extract nitrate at 15 minute intervals
temp_arik1 <- temp_arik1$TSW_5min %>%
  as_tibble() %>%
  janitor::clean_names()  %>%
  filter(horizontal_position == "101") %>%
  mutate(start_date_time = ymd_hms(start_date_time) - second(start_date_time)) %>%
  select(c(start_date_time, surf_water_temp_mean, final_qf)) %>%
  rename(
    temp_mean = surf_water_temp_mean,
    label_temp_mean = final_qf
  )


# download data of interest - ELEVATION
if(file.exists("level_arik1.rda")) {
  level_arik1 <- readRDS("level_arik1.rda")
} else {
  level_arik1 <- neonUtilities::loadByProduct(
    dpID = "DP1.20016.001",
    site = "ARIK",
    startdate = "2018-10",
    enddate = "2020-10",
    package = "expanded",
    token = Sys.getenv("NEON_TOKEN"),
    check.size = FALSE
  )
  saveRDS(level_arik1, "level_arik1.rda")
}

# Take only water quality at the down station
level_arik1 <- level_arik1$EOS_5_min %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  filter(horizontal_position == "101") %>%
  select(c(5,7,37)) %>%
  mutate(
    start_date_time = ymd_hms(start_date_time) - second(start_date_time)
  ) %>%
  rename(
    elev = surfacewater_elev_mean,
    label_elev = s_wat_elev_final_qf,
  )


# Join water_quality and nitrate, and clean up
data_arik1 <- left_join(
  nitrate_arik1,
  temp_arik1,
  by = "start_date_time")

data_arik1 <- left_join(
  data_arik1,
  water_quality_arik1,
  by = "start_date_time") 

data_arik1 <- left_join(
  data_arik1,
  level_arik1,
  by = "start_date_time") 


data_arik1 <- data_arik1 %>%
  distinct(start_date_time, .keep_all=TRUE)

label_nitrate_mean = if_else(data_arik1$nitrate_mean  > 25, 1, data_arik1$label_nitrate_mean)
data_arik1$daily_cycle = as.factor(rep(seq(from = 1, to = 24*15, by = 1), length.out = 73056))
data_arik1$log_turbidity <- log(data_arik1$turbidity +1)
data_arik1$n_1 <- dplyr::lag(data_arik1$nitrate_mean)
data_arik1$n_2 <- dplyr::lag(data_arik1$nitrate_mean, n = 2)

##############################################
# 2 - QUIT abnormal data : with "label" = 1 #
############################################
data_without_abnormal<-data_arik1 %>%
  mutate(
    spec_cond = if_else(label_spec_cond == 1, NA_real_, spec_cond),
    oxygen = if_else(label_oxygen == 1, NA_real_, oxygen),
    oxygen_sat = if_else(label_oxygen_sat == 1, NA_real_, oxygen_sat),
    turbidity = if_else(label_turbidity == 1, NA_real_, turbidity),
    log_turbidity = if_else(label_turbidity == 1, NA_real_, log_turbidity),
    ph = if_else(label_ph == 1, NA_real_, ph),
    chloro = if_else(label_chloro == 1, NA_real_, chloro),
    f_dom = if_else(label_f_dom == 1, NA_real_, f_dom),
    temp_mean = if_else(label_temp_mean == 1, NA_real_, temp_mean),
    elev = if_else(label_elev == 1, NA_real_ , elev),
    nitrate_mean = if_else(label_nitrate_mean == 1, NA_real_ , nitrate_mean),
    nitrate_mean = if_else(nitrate_mean <= 0, NA_real_ , nitrate_mean)
  )  %>%
  distinct(start_date_time, .keep_all=TRUE)


## Plot data
data_without_abnormal %>%
  dplyr::select(start_date_time, nitrate_mean, log_turbidity, temp_mean, spec_cond, oxygen, elev ) %>%
  pivot_longer(-start_date_time, names_to = "variable") %>%
  ggplot(aes(x = start_date_time , y = value, col = variable)) +
  geom_point() +
  facet_grid(rows = vars(variable),  scales = "free", space="free_x") +
  labs(y = "", x = "Time") +
  scale_x_datetime(date_labels ="%Y-%m") +
  theme(legend.position = "none")


##################
# 3 - Framework #
################

#################################
## 3- 1 Ponctual data anomalies #
################################
# example for 30%

data_res <- list()

for(j in 1:100) {
#Create a start (no available data before) a end (strange data after) 
data_model <- data_without_abnormal %>%
  dplyr::select(nitrate_mean, log_turbidity, temp_mean, spec_cond, oxygen, elev, n_1,n_2 )
start<-which(data_model$nitrate_mean != "NA")[1] # Data before are NA
end<- 56976 # data after are strange 
data_model <- data_model[start:end,]


# Create a sorted vector, but keeping the first and the three last data non NA (for ARIMA not to fail) / Choose the percent of keeping data
dt = c(1,sort(sample((2:(nrow(data_model)-3)), (nrow(data_model)-6)*.7, replace = FALSE)),(nrow(data_model)-2),(nrow(data_model)-1),nrow(data_model))
data_model[-dt,1] <- NA

#crete empty objects that will be filled in the loop

se_gam<-vector()
lower_arima<-vector()
upper_arima<-vector()
data_modelled_gam <- data_model
data_modelled_arima <- data_model 

# 1/ predict all posible NA with a GAM

model1 <- mgcv::gam(nitrate_mean ~ s(temp_mean, k = 12) + s(spec_cond, k = 12) + s(oxygen, k = 12) 
                            + s(log_turbidity, k = 12)  + elev + n_1 + n_2, data = data_model)

predicted <- mgcv::predict.gam(model1, newdata = data_modelled_gam, se.fit = TRUE)
    data_modelled_gam$nitrate_mean<- predicted$fit
    data_modelled_gam$n_1 <- dplyr::lag(predicted$fit)
    data_modelled_gam$n_2 <- dplyr::lag(predicted$fit, n = 2)
    se_gam<-predicted$se.fit
# 2/ predict all posible NA with an ARIMA
    
for(i in 500:(nrow(data_modelled_arima)-3)) {

if (is.na(data_modelled_arima$nitrate_mean[i])) {
  
    model <- auto.arima(data_modelled_arima$nitrate_mean[(i - 500):(i-1)])
    na_seq <- which(!is.na(data_modelled_arima$nitrate_mean[i:nrow(data_model)])) # find the lenght of the NA sequence to predict it / na_seq[1] is the first non na data
    predicted <- forecast(model, h = (na_seq[1]-1), level = 95) # na_seq[1]-1 is the last na value
    data_modelled_arima$nitrate_mean[i:(i+na_seq[1]-2)] <- predicted$mean
    data_modelled_arima$n_1[(i+1):(i+na_seq[1]-1)] <- predicted$mean
    data_modelled_arima$n_2[(i+2):(i+na_seq[1])] <- predicted$mean
    lower_arima[i:(i+na_seq[1]-2)] <- predicted$lower  
    upper_arima[i:(i+na_seq[1]-2)] <- predicted$upper

} 

}

      data_res[[j]] <- data_frame("time" = data_without_abnormal$start_date_time[start:end],
                     "real" = data_without_abnormal$nitrate_mean[start:end],
                     "predicted" = ifelse(is.na(data_model$nitrate_mean), ifelse(is.na(data_modelled_gam$nitrate_mean), data_modelled_arima$nitrate_mean,data_modelled_gam$nitrate_mean),  NA),
                     "se_up" =  ifelse(is.na(data_modelled_gam$nitrate_mean),upper_arima,data_modelled_gam$nitrate_mean + se_gam),
                     "se_low" = ifelse(is.na(data_modelled_gam$nitrate_mean),lower_arima,data_modelled_gam$nitrate_mean - se_gam))    
    

print(j) 

}

save(data_res, file = "results_30percent.RData") # save in RData

####################################
## 3- 2 Sequences of abnormal data #
###################################
# example for 10 times one weeks

data_res <- list()

for(j in 1:100) {
  #Create a start (no available data before) a end (strange data after) 
  data_model <- data_without_abnormal %>%
    dplyr::select(nitrate_mean, log_turbidity, temp_mean, spec_cond, oxygen, elev, n_1,n_2 )
  start<-which(data_model$nitrate_mean != "NA")[1] # Data before are NA
  end<- 56976 # data after are strange 
  data_model <- data_model[start:end,]
  
  
  # Create a sorted vector, but keeping the first and the three last data non NA (for ARIMA not to fail) /
  
  starts <-sample(2:(nrow(data_model)-672), 10, replace = FALSE)
  
  na_sequences<-vector()
  for(i in 1:10){
    na_sequences <- c(na_sequences,seq(from = starts[i], to = starts[i] + 672, by = 1 ))
  }
  
  na_sequences<-unique(na_sequences)
  
  data_model[na_sequences,1] <- NA
  
  #crete empty objects that will be filled in the loop
  
  se_gam<-vector()
  lower_arima<-vector()
  upper_arima<-vector()
  data_modelled_gam <- data_model
  data_modelled_arima <- data_model 
  
  # 1/ predict all posible NA with a GAM
  
  model1 <- mgcv::gam(nitrate_mean ~ s(temp_mean, k = 12) + s(spec_cond, k = 12) + s(oxygen, k = 12) 
                      + s(log_turbidity, k = 12)  + elev + n_1 + n_2, data = data_model)
  
  predicted <- mgcv::predict.gam(model1, newdata = data_modelled_gam, se.fit = TRUE)
  data_modelled_gam$nitrate_mean<- predicted$fit
  data_modelled_gam$n_1 <- dplyr::lag(predicted$fit)
  data_modelled_gam$n_2 <- dplyr::lag(predicted$fit, n = 2)
  se_gam<-predicted$se.fit
  # 2/ predict all posible NA with an ARIMA
  
  for(i in 500:(nrow(data_modelled_arima)-3)) {
    
    if (is.na(data_modelled_arima$nitrate_mean[i])) {
      
      model <- auto.arima(data_modelled_arima$nitrate_mean[(i - 500):(i-1)])
      na_seq <- which(!is.na(data_modelled_arima$nitrate_mean[i:nrow(data_model)])) # find the lenght of the NA sequence to predict it / na_seq[1] is the first non na data
      predicted <- forecast(model, h = (na_seq[1]-1), level = 95) # na_seq[1]-1 is the last na value
      data_modelled_arima$nitrate_mean[i:(i+na_seq[1]-2)] <- predicted$mean
      data_modelled_arima$n_1[(i+1):(i+na_seq[1]-1)] <- predicted$mean
      data_modelled_arima$n_2[(i+2):(i+na_seq[1])] <- predicted$mean
      lower_arima[i:(i+na_seq[1]-2)] <- predicted$lower  
      upper_arima[i:(i+na_seq[1]-2)] <- predicted$upper
      
    } 
    
  }
  
  data_res[[j]] <- data_frame("time" = data_without_abnormal$start_date_time[start:end],
                              "real" = data_without_abnormal$nitrate_mean[start:end],
                              "predicted" = ifelse(is.na(data_model$nitrate_mean), ifelse(is.na(data_modelled_gam$nitrate_mean), data_modelled_arima$nitrate_mean,data_modelled_gam$nitrate_mean),  NA),
                              "se_up" =  ifelse(is.na(data_modelled_gam$nitrate_mean),upper_arima,data_modelled_gam$nitrate_mean + se_gam),
                              "se_low" = ifelse(is.na(data_modelled_gam$nitrate_mean),lower_arima,data_modelled_gam$nitrate_mean - se_gam))    
  
  
  print(j) 
  
}

save(data_res, file = "10weeks.RData") # save in RData


###############
# 4 - Results #
##############
# example of code for 30% of ponctual data removed

load(file = "results_30percent.RData") # not possible to give a name directly, it keep the data-res name
results_30percent <- data_res

# So ye have a list of 100 simulation of predictions with 30% of data is NA 



# Delete IC when there is no predicted data - delete prediction when there is a real data
dat30<-list()
for (i in 1:100){
  dat30[[i]]<-data.frame("time" = results_30percent[[i]]$time, 
        "real" = results_30percent[[i]]$real,
        "predicted" =  results_30percent[[i]]$predicted,
        "se_up" = ifelse(is.na(results_30percent[[i]]$predicted),NA, results_30percent[[i]]$se_up),
        "se_low" = ifelse(is.na(results_30percent[[i]]$predicted),NA,                     results_30percent[[i]]$se_low))  

dat30[[i]]$real[dat30[[i]]$real<0]<- NA
}

# calculate RMSE
rmse30percent<-vector()
for(i in 1:100) {
rmse1 <- cbind(dat30[[i]][,2],dat30[[i]][,3])
rmse1 <- rmse1[complete.cases(rmse1),]
rmse30percent[i] <- Metrics::rmse(rmse1[,1], rmse1[,2])  
}



# PWPI
# create sensor precision interval
dat30modif<-list()
for (i in 1:100){
  dat30modif[[i]]<-data.frame("time" = dat30[[i]]$time, 
        "real" = dat30[[i]]$real,
        "predicted" =  dat30[[i]]$predicted,
        "real_up_limit" = ifelse(dat30[[i]]$real >= 20, dat30[[i]]$real + (0.1*dat30[[i]]$real), dat30[[i]]$real + 2),
        "real_down_limit" = ifelse(dat30[[i]]$real >= 20, dat30[[i]]$real - (0.1*dat30[[i]]$real), dat30[[i]]$real - 2),
        "se_up" = dat30[[i]]$se_up,
        "se_low" = dat30[[i]]$se_low)  

}

# proportion of predicted values inside the error confidence of the sonde (10%)
mat <- matrix( ,nrow = 100, ncol = nrow(dat30modif[[i]]))
for(i in 1:100) {
  full <- dat30modif[[i]][complete.cases(dat30modif[[i]]),]
      for (j in 1:nrow(full)) {
        vec_real <- seq(from = full[j,]$real_down_limit, to = full[j,]$real_up_limit, by = 0.01)
          vec_pred <- full[j,]$predicted
          mat[i,j] <- length(intersect(vec_real, round(vec_pred,2))) # number of intersection by line 
      }

}


prop_inside_CI_30_mod <-vector()
for(i in 1 : 100){
prop_inside_CI_30_mod[i] <- length(which(mat[i,] != 0)) / nrow(dat30modif[[i]][complete.cases(dat30modif[[i]]),])
}

boxplot(prop_inside_CI_30_mod, main = "With ajusted CI")

# Code for figure 4
rmse_final <- tibble("rmse" = c(rmse20percent,rmse30percent,rmse40percent),
                             "Percent_NA" = c(rep("20",100),
rep("30",100),rep("40",100)))

gg1<-ggplot(rmse_final, aes(x = Percent_NA, y = rmse, color = Percent_NA)) +
  ylab(expression("RMSE "*mu*"mol/L")) +
  xlab("Proportion of point data removed")+
  geom_boxplot() +
  scale_color_manual(values = c("brown1","maroon3","blue")) +
  theme(legend.position='none',,
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14))



prop_inside_CI_mod <- tibble("Prop" = c(prop_inside_CI_20_mod,prop_inside_CI_30_mod,prop_inside_CI_40_mod),
                             "Percent_NA" = c(rep("20",100),
rep("30",100),rep("40",100)))

gg2<-ggplot(prop_inside_CI_mod, aes(x = Percent_NA, y = Prop, color = Percent_NA)) +
  geom_boxplot() +
  ylab("PWPI (%)") +
  xlab("Proportion of point data removed")+ 
  scale_color_manual(values = c("brown1","maroon3","blue")) +
  theme(legend.position='none',
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14))

require(gridExtra)
grid.arrange(gg1, gg2, ncol=2)



## Figure Nitrate time series
Sys.setlocale(locale = "en_US.UTF-8")
i=5
gg5<-ggplot(data = dat30modif[[i]])+
  geom_point(aes(x = time, y = real)) + 
  geom_ribbon(aes(x = time, y = real,ymin= real-(0.1*real), ymax=real + (0.1*real)),  alpha = 0.5) +
  xlab("") + ylab(expression("Nitrate concentration ("*mu*"mol/L)"))+ 
  theme(legend.position='none',
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14))





######################################
# 5 - Use the method on the ARIK TS #
#####################################


data_model <- data_without_abnormal %>%
  dplyr::select(nitrate_mean, log_turbidity, temp_mean, spec_cond, oxygen, elev, n_1,n_2 )
start<-which(data_model$nitrate_mean != "NA")[1] # Data before are NA
end<- 56976 # data after are strange 
data_model <- data_model[start:end,]


#crete empty objects that will be filled in the loop

se_gam<-vector()
lower_arima<-vector()
upper_arima<-vector()
data_modelled_gam <- data_model
data_modelled_arima <- data_model 

# 1/ predict all posible NA with a GAM

model1 <- mgcv::gam(nitrate_mean ~ s(temp_mean, k = 12) + s(spec_cond, k = 12) + s(oxygen, k = 12) 
                            + s(log_turbidity, k = 12)  + elev + n_1 + n_2, data = data_model)

predicted <- mgcv::predict.gam(model1, newdata = data_modelled_gam, se.fit = TRUE)
    data_modelled_gam$nitrate_mean<- predicted$fit
    data_modelled_gam$n_1 <- dplyr::lag(predicted$fit)
    data_modelled_gam$n_2 <- dplyr::lag(predicted$fit, n = 2)
    se_gam<-predicted$se.fit
# 2/ predict all posible NA with an ARIMA
    
for(i in 500:(nrow(data_modelled_arima)-3)) {

if (is.na(data_modelled_arima$nitrate_mean[i])) {
  
    model <- auto.arima(data_modelled_arima$nitrate_mean[(i - 500):(i-1)])
    na_seq <- which(!is.na(data_modelled_arima$nitrate_mean[i:nrow(data_model)])) # find the lenght of the NA sequence to predict it / na_seq[1] is the first non na data
    predicted <- forecast(model, h = (na_seq[1]-1), level = 95) # na_seq[1]-1 is the last na value
    data_modelled_arima$nitrate_mean[i:(i+na_seq[1]-2)] <- predicted$mean
    data_modelled_arima$n_1[(i+1):(i+na_seq[1]-1)] <- predicted$mean
    data_modelled_arima$n_2[(i+2):(i+na_seq[1])] <- predicted$mean
    lower_arima[i:(i+na_seq[1]-2)] <- predicted$lower  
    upper_arima[i:(i+na_seq[1]-2)] <- predicted$upper

} 

}

      data_res <- data_frame("time" = data_without_abnormal$start_date_time[start:end],
                     "real" = data_without_abnormal$nitrate_mean[start:end],
                     "predicted" = ifelse(is.na(data_model$nitrate_mean), ifelse(is.na(data_modelled_gam$nitrate_mean), data_modelled_arima$nitrate_mean,data_modelled_gam$nitrate_mean),  NA),
                     "se_up" =  ifelse(is.na(data_modelled_gam$nitrate_mean),upper_arima,data_modelled_gam$nitrate_mean + se_gam),
                     "se_low" = ifelse(is.na(data_modelled_gam$nitrate_mean),lower_arima,data_modelled_gam$nitrate_mean - se_gam))    
    



data_res<-data.frame("time" = data_res$time, 
        "real" = data_res$real,
        "predicted" =  ifelse(is.na(data_res$real),data_res$predicted,NA),
        "se_up" = ifelse(is.na(data_res$predicted),NA,data_res$se_up),
        "se_low" = ifelse(is.na(data_res$predicted),NA,data_res$se_low),
        "real_up_limit" = ifelse(data_res$real >= 20, data_res$real + (0.1*data_res$real), data_res$real + 2),
        "real_down_limit" = ifelse(data_res$real >= 20, data_res$real - (0.1*data_res$real), (data_res$real - 2)) ) 


# Code for figure 3

Sys.setlocale(locale = "en_US.UTF-8") #months in english


gg1<-ggplot(data = data_res[2000:3000,]) +
  geom_point(aes(x = time, y = real), size = .2) + 
    geom_ribbon(aes(x = time, y = real, ymin= real_down_limit, ymax=real_up_limit),  alpha = 0.5) +
  geom_point(aes(x = time, y = predicted), col= "green") +
    geom_errorbar(aes(x = time, y = predicted, ymin=se_up, ymax=se_low),  col= "green")+
  xlab("") + ylab(expression("Nitrate ("*mu*"mol/L)"))+ 
  annotate(geom="text", x=data_res[1950,1], y=4.1, label="(a)", size = 5)+ 
  theme(legend.position='none',
        axis.text.y = element_text(size=14))

gg2<-ggplot(data = data_res[26050:26150,]) +
  geom_point(aes(x = time, y = real), size = .2) + 
    geom_ribbon(aes(x = time, y = real,ymin= real_down_limit, ymax=real_up_limit),  alpha = 0.5) +
  geom_point(aes(x = time, y = predicted), col= "green") +
    geom_errorbar(aes(x = time, y = predicted, ymin=se_up, ymax=se_low),  col= "green", alpha = 0.5)+
  xlab("") + ylab(expression("Nitrate ("*mu*"mol/L)"))+
  annotate(geom="text", x=data_res[26050,1], y=18.5, label="(b)", size = 5)+ 
  theme(legend.position='none',
        axis.text.y = element_text(size=14))

gg3<-ggplot(data = data_res[31000:31400,]) +
  geom_point(aes(x = time, y = real), size = .2) + 
    geom_ribbon(aes(x = time, y = real,ymin= real_down_limit, ymax=real_up_limit),  alpha = 0.5) +
  geom_point(aes(x = time, y = predicted), col= "green") +
    geom_errorbar(aes(x = time, y = predicted, ymin=se_up, ymax=se_low),  col= "green", alpha = 0.5)+
  xlab("") + ylab(expression("Nitrate ("*mu*"mol/L)"))+
  annotate(geom="text", x=data_res[30980,1], y=6, label="(c)", size = 5)+ 
  theme(legend.position='none',
        axis.text.y = element_text(size=14))

gg4<-ggplot(data = data_res[6000:14000,]) +
  geom_point(aes(x = time, y = real)) + 
    geom_ribbon(aes(x = time, y = real,ymin= real_down_limit, ymax=real_up_limit),  alpha = 0.5) +
  geom_point(aes(x = time, y = predicted), col= "green") +
    geom_errorbar(aes(x = time, y = predicted, ymin=se_up, ymax=se_low),  col= "green", alpha = 0.5)+
  xlab("") + ylab(expression("Nitrate ("*mu*"mol/L)"))+
  annotate(geom="text", x=data_res[6000,1], y=80, label="(d)", size = 5)+ 
  theme(legend.position='none',
        axis.text.y = element_text(size=14))



ggpubr::ggarrange(gg1,gg2,gg3, gg4,  nrow = 2, ncol = 2)


na_in_ts<-length(data_res$real[which(is.na(data_res$real) == TRUE)])
number_of_still_na <- nrow(data_res[which(data_res$real[which(is.na(data_res$real) == TRUE)] &&
data_res$predicted[which(is.na(data_res$predicted) == FALSE)]),])
sd_range_min<- min(data_res$se_up - data_res$se_low, na.rm = TRUE)
sd_range_max<- max(data_res$se_up - data_res$se_low, na.rm = TRUE)
sd_range_med<- median(data_res$se_up - data_res$se_low, na.rm = TRUE)

dat_table<-data.frame("Missing Value in TS" = na_in_ts,
                      "Non predicted Missing Value" = number_of_still_na,
                      "Range min of prediction interv." = round(sd_range_min,2),
                      "Range max of prediction interv." = round(sd_range_max,2),
                      "Median range of prediction interv." = round(sd_range_med,2))
dat_table


