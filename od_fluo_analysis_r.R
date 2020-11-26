#plate reader analysis script.


#load required packages
library(tidyverse)
library(readxl)
library(plyr)


#1. import dataframe
X20_10_m9_238 <- read_excel("C:/Users/ah2218/Desktop/PHD/plate reader data/20 10 m9 238.xlsx")

#select rows and columns corresponding to OD and fluorescence tables
OD700 <- data.matrix(X20_10_m9_238[17:65, 3:98])
fluo_gain_60 <- data.matrix(X20_10_m9_238[66:114,3:98])
fluo_gain_80 <- data.matrix(X20_10_m9_238[115:163,3:98])



#specify time duration and step size for measurments
time <- seq(0, 12, by = 0.25)


#blank the samples & add time
# subtract blanks: media for OD, wt cells for fluorescence
#add time column

OD_blanked <- cbind(time, OD700 - rowMeans(subset(OD700, select = c(1,2,3,4,5))))
GFP60_blanked <- cbind(time, fluo_gain_60 - rowMeans(subset(fluo_gain_60, select = c(13,14,15,16,17))))
GFP80_blanked <- cbind(time, fluo_gain_80 - rowMeans(subset(fluo_gain_80, select = c(13,14,15,16,17))))


#remove NEGATIVES - replace with zero
OD_blanked[OD_blanked<0] <- 0
GFP60_blanked[GFP60_blanked<0] <- 0
GFP80_blanked[GFP80_blanked<0] <- 0



#GROWTH RATE


#empty matrix to store values
growth_rate <- matrix(nrow = 47, ncol = 96)
#loop through each well and each time point to calculate growth rate
#formula: growth_rate(t2) = (log(od(t3))-log(od(t1)))/(t3/t1)
for(i in 2:97)
{
  for(j in 2:48)
  {
    growth_rate[j-1,i-1] <- (log(OD_blanked[j+1,i])-log(OD_blanked[j-1,i]))/0.5
  }
}

#REMOVE INF, NAN 
growth_rate[!is.finite(growth_rate)] <- 0

#time here excludes the first and last points, as no rates for those
t2 <- seq(0.25,11.75, by=0.25)





# prepare for plotting by makingdfs with all the means for the things we want to plot
# in a ggplot friendly form
prepare_for_plot <- function(c1, timex) {
  BP16 <-  rowMeans(subset(c1, select = c(14,15,16,17,18)))
  x1mM <-rowMeans(subset(c1, select = c(26,27,28,29,30)))
  x2mM <-  rowMeans(subset(c1, select = c(38,39,40,41,42)))
  x0_5mM <- rowMeans(subset(c1, select = c(50,51,52,53,54)))
  x0mM <-  rowMeans(subset(c1, select = c(62,63,64,65,66)))
  
  means <- as.data.frame(cbind(timex, BP16, x1mM, x2mM, x0_5mM, x0mM))
  
  df1x <- means %>%
    select(timex, BP16, x1mM,x2mM,x0_5mM, x0mM) %>%
    gather(key = "variable", value = "value", -timex)
  
  return(df1x)
  
}

OD_means <- prepare_for_plot(OD_blanked,OD_blanked[,1])


# GROWTH RATE
GR_means <- prepare_for_plot(growth_rate,t2)



#GFP60

GFP60_means <- prepare_for_plot(GFP60_blanked,time)


#GFP80
GFP80_means <- prepare_for_plot(GFP80_blanked,time)



# GFP PER OD

# USE GFP80 as an example

GFP_per_od <- cbind(OD_blanked[1],GFP80_blanked[,-1]/OD_blanked[,-1])
GFP_per_od <- data.matrix(GFP_per_od)
GFP_per_od[!is.finite(GFP_per_od)] <- 0


GFP_OD_means <- prepare_for_plot(GFP_per_od,time)

# PLOT
# beware: ORDER OF PLOTTING IS ALPHABETICAL!!!!!

line_plot <- function(df1,timex) {
  ggplot(df1, aes(x = timex, y = value)) + 
    geom_line(aes(color = variable)) + 
    scale_color_manual(values = c("darkred", "steelblue", 'yellow', 'green', 'orange'),limits = c('BP16', 'x1mM', 'x2mM', 'x0_5mM', 'x0mM'), labels=c("BP16", "1mM", "2mM", '0.5mM', '0mM'), name = 'Induction') +
    labs(x="Time (hrs)", y = "Fluorescence (AU)")
}

plot1 <- line_plot(OD_means,time)
plot2 <- line_plot(GR_means,t2)
plot3 <- line_plot(GFP60_means,time)
plot4 <- line_plot(GFP80_means,time)


