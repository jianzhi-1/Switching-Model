# Import libraries
library(ggplot2) # plotting
library(naniar) # fill na
library(visdat)
library(dplyr)
library(ggmosaic) # Mosaic plot
library(gridExtra) # multipanel plots
library(MASS) # mvrnorm
library(mvtnorm)
library(Matrix)
library(glmnet)
library(nortest)
library(nnet) # multinomial logistic regression
library(MASS) # ordinal logistic regression
library(MSwM)

### 0. Set up
set.seed(42)
data = read.csv("230Adata.csv")
df = as.data.frame(data)
n = nrow(df) # 915
colnames(df) # "DATE"     "UNRATE"   "FEDFUNDS" "GDP"      "CPIAUCSL" "BOPTEXP"  "SPY"

subdf = df[,c(1,2,3,4,5)]
colnames(subdf) # "DATE"     "UNRATE"   "FEDFUNDS" "GDP"      "CPIAUCSL"
subdf = na.omit(subdf)
nsub = nrow(subdf) # 834 
# subdf is the data with all populated fields, contiguous

n = nrow(subdf)
n.train = 700
n.test = n - n.train

subdf.train = subdf[1:n.train,]
subdf.test = subdf[(n.train+1):n,]

### 1. Data visualisation
summary(subdf)

# line plot of UNRATE
ggplot() + geom_line(aes(x=1:n.train,y=subdf.train$UNRATE))
ggplot() + geom_line(aes(x=1:(n.train-1),y=diff(subdf.train$UNRATE, 1)))

acf(diff(subdf.train$UNRATE, 1)) # some significant ACF
pacf(diff(subdf.train$UNRATE, 1)) # some significant PACF
# PACF suggests yt should depend on y(t-1), y(t-2), y(t-3), y(t-4), y(t-12), y(t-24)

subsubdf = subdf.train[1:(nrow(subdf.train)-1),]
subsubdf$diff = diff(subdf.train$UNRATE, 1)
cor(subsubdf[,c("diff","UNRATE","FEDFUNDS","GDP","CPIAUCSL")]) # correlation matrix
# diff positively correlated with FEDFUNDS, GDP,  CPIAUCSL
#      negatively correlated with UNRATE (previous step)
#           diff      UNRATE    FEDFUNDS         GDP    CPIAUCSL
# diff      1.00000000 -0.05214588  0.11122681  0.0249548  0.01972545

ggplot() + geom_histogram(colour="black", fill="white", aes(subsubdf$diff))

ggplot(data = subsubdf, aes(x = UNRATE, y = diff)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE) +
  labs(x = "UNRATE", y = "diff") 
# UNRATE and diff has negative correlation (not significant)
cor.test(subsubdf$UNRATE, subsubdf$diff) # t = -1.3786, df = 697, p-value = 0.1685

ggplot(data = subsubdf, aes(x = FEDFUNDS, y = diff)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE) +
  labs(x = "FEDFUNDS", y = "diff")
# FEDFUNDS and diff has positive correlation (significant)
cor.test(subsubdf$FEDFUNDS, subsubdf$diff) # t = 2.9548, df = 697, p-value = 0.003234

ggplot(data = subsubdf, aes(x = GDP, y = diff)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE) +
  labs(x = "GDP", y = "diff")
# GDP and diff has positive correlation (not significant)
cor.test(subsubdf$GDP, subsubdf$diff) # t = 0.65903, df = 697, p-value = 0.5101

ggplot(data = subsubdf, aes(x = CPIAUCSL, y = diff)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE) +
  labs(x = "CPIAUCSL", y = "diff")
# CPI and diff has positive correlation (not significant)
cor.test(subsubdf$CPIAUCSL, subsubdf$diff) # t = 0.52087, df = 697, p-value = 0.6026

calculate_mse <- function(observed_values, predicted_values) {
  squared_diff <- (observed_values - predicted_values)^2
  mse <- mean(squared_diff)
  return(mse)
}

### 2. Baseline Model
# Simply predict using the previous month's unemployment rate
