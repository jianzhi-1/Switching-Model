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

# 2.1 Simply predict using the previous month's unemployment rate
y.pred = c(subdf.train$UNRATE[n.train], subdf.test$UNRATE[1:(n.test - 1)])
calculate_mse(subdf.test$UNRATE, y.pred) # 0.9312687
calculate_mse_cutoff(subdf.test$UNRATE, y.pred) # 0.027
ggplot() + geom_line(aes(x=1:n.test,y=y.pred), color="blue") + 
  geom_line(aes(x=1:n.test,y=subdf.test$UNRATE), color="darkgreen")
# Prediction lags behind by 1

# 2.2 Simply use the previous two data points to predict the next month's unemployment rate
y.pred = 2*c(subdf.train$UNRATE[n.train], subdf.test$UNRATE[1:(n.test - 1)]) - c(subdf.train$UNRATE[n.train-1], subdf.train$UNRATE[n.train], subdf.test$UNRATE[1:(n.test - 2)])
calculate_mse(subdf.test$UNRATE, y.pred) # 1.842836
calculate_mse_cutoff(subdf.test$UNRATE, y.pred) # 0.0568
ggplot() + geom_line(aes(x=1:n.test,y=y.pred), color="blue") + 
  geom_line(aes(x=1:n.test,y=subdf.test$UNRATE), color="darkgreen")

### 3. Ridge Regression
# Idea is to use other covariates to improve prediction error

# 3.1 Use the covariates at the same time: UNRATE    FEDFUNDS         GDP    CPIAUCSL
subdf.train.ridge = subdf.train[1:(n.train - 1),]
subdf.train.ridge$diff = diff(subdf.train$UNRATE, 1)
subdf.train.ridge$DATE = NULL

y = subdf.train.ridge$diff
subdf.train.ridge$diff = NULL
X = as.matrix(subdf.train.ridge)

cv_model = cv.glmnet(X, y, alpha = 0)
best_lambda = cv_model$lambda.min # 0.00212585 (this is pathetic)
ridge_model_best = glmnet(X, y, alpha = 0, lambda = best_lambda)

coef(ridge_model_best)
# s0
# (Intercept)  6.838000e-03
# UNRATE      -7.232122e-03
# FEDFUNDS     9.390256e-03
# GDP          1.498277e-05
# CPIAUCSL    -8.679753e-04

fitted_values = predict(ridge_model_best, newx = X)
calculate_mse(y, fitted_values) # 0.03546902
ggplot() + geom_line(aes(x=1:length(fitted_values),y=fitted_values), color="blue", linewidth=0.2) + 
  geom_line(aes(x=1:length(y),y=y), color="darkgreen", linewidth=0.2)

# Evaluate on test set
subdf.test.ridge = subdf.test[1:(n.test-1),]
subdf.test.ridge$diff = diff(subdf.test$UNRATE, 1)
subdf.test.ridge$DATE = NULL

y = subdf.test.ridge$diff
subdf.test.ridge$diff = NULL
X = as.matrix(subdf.test.ridge)

y.pred = predict(ridge_model_best, newx = X)
calculate_mse(y, y.pred) # 0.9388009
calculate_mse_cutoff(y, y.pred) # 0.03039402
ggplot() + geom_line(aes(x=1:length(y.pred),y=y.pred), color="blue", linewidth=0.2) + 
  geom_line(aes(x=1:length(y),y=y), color="darkgreen", linewidth=0.2)

# Cumulative Plot
y.pred.cumulative = y.pred + subdf.test$UNRATE[1:(n.test-1)]
y = subdf.test[2:n.test,]$UNRATE
ggplot() + geom_line(aes(x=1:length(y.pred.cumulative),y=y.pred.cumulative), color="blue", linewidth=0.2) + 
  geom_line(aes(x=1:length(y),y=y), color="darkgreen", linewidth=0.2)

# 3.2 Use in addition the lags y(t-1), y(t-2), y(t-3), y(t-4), y(t-12), y(t-24) suggested by PACF
subdf.train.ridge.ii = subdf.train
subdf.train.ridge.ii$DATE = NULL
subdf.train.ridge.ii <- subdf.train.ridge.ii %>%
  mutate(
    lag_1 = lag(UNRATE, 1),
    lag_2 = lag(UNRATE, 2),
    lag_3 = lag(UNRATE, 3),
    lag_4 = lag(UNRATE, 4),
    lag_12 = lag(UNRATE, 12),
    lag_24 = lag(UNRATE, 24)
  )
subdf.train.ridge.ii = na.omit(subdf.train.ridge.ii) # omit rows with no lags
nsub = nrow(subdf.train.ridge.ii)
subdf.train.ridge.ii.final = subdf.train.ridge.ii[1:(nsub - 1),]
subdf.train.ridge.ii.final$diff = diff(subdf.train.ridge.ii$UNRATE, 1)

y = subdf.train.ridge.ii.final$diff
subdf.train.ridge.ii.final$diff = NULL
X = as.matrix(subdf.train.ridge.ii.final)

cv_model = cv.glmnet(X, y, alpha = 0)
best_lambda = cv_model$lambda.min # 0.004143626 (this is pathetic)
ridge_model_best = glmnet(X, y, alpha = 0, lambda = best_lambda)

coef(ridge_model_best)
# s0
# (Intercept)  7.108730e-02
# UNRATE       5.511277e-02
# FEDFUNDS     8.260150e-03
# GDP          6.931299e-06
# CPIAUCSL    -2.727976e-04
# lag_1        5.537325e-02
# lag_2       -1.089737e-02
# lag_3       -3.174353e-02
# lag_4       -6.846218e-02
# lag_12      -1.422965e-02
# lag_24      -5.662170e-03

# The small lag coefficients are concerning, because their tick size is 0.1, so a change in 0.1 units correspond to a 0.001 change, which is 2 magnitudes less than a tick size of diff (0.1)
# I don't expect good results here

fitted_values = predict(ridge_model_best, newx = X)
calculate_mse(y, fitted_values) # 0.03068849 lower apparent error!
ggplot() + geom_line(aes(x=1:length(fitted_values),y=fitted_values), color="blue", linewidth=0.2) + 
  geom_line(aes(x=1:length(y),y=y), color="darkgreen", linewidth=0.2)

# Evaluate on test set
subdf.train.ridge.ii.test = subdf[(n.train-23):n,]
subdf.train.ridge.ii.test$DATE = NULL
subdf.train.ridge.ii.test <- subdf.train.ridge.ii.test %>%
  mutate(
    lag_1 = lag(UNRATE, 1),
    lag_2 = lag(UNRATE, 2),
    lag_3 = lag(UNRATE, 3),
    lag_4 = lag(UNRATE, 4),
    lag_12 = lag(UNRATE, 12),
    lag_24 = lag(UNRATE, 24)
  )
subdf.train.ridge.ii.test = na.omit(subdf.train.ridge.ii.test)
  
subdf.train.ridge.ii.test = subdf.train.ridge.ii.test[1:(n.test-1),]
subdf.train.ridge.ii.test$diff = diff(subdf.test$UNRATE, 1)

y = subdf.train.ridge.ii.test$diff
subdf.train.ridge.ii.test$diff = NULL
X = as.matrix(subdf.train.ridge.ii.test)

y.pred = predict(ridge_model_best, newx = X)
calculate_mse(y, y.pred) # 1.012032
calculate_mse_cutoff(y, y.pred) # 0.02651024 # better than naive method
ggplot() + geom_line(aes(x=1:length(y.pred),y=y.pred), color="blue", linewidth=0.2) + 
  geom_line(aes(x=1:length(y),y=y), color="darkgreen", linewidth=0.2)

# Cumulative Plot
y.pred.cumulative = y.pred + subdf.test$UNRATE[1:(n.test-1)]
y = subdf.test[2:n.test,]$UNRATE
ggplot() + geom_line(aes(x=1:length(y.pred.cumulative),y=y.pred.cumulative), color="blue", linewidth=0.2) + 
  geom_line(aes(x=1:length(y),y=y), color="darkgreen", linewidth=0.2)

### 4. Indicators on past drops

### 4. Indicators on past drops

### 5. Markov Switching Model

subdf.train.ridge.ii = subdf.train
subdf.train.ridge.ii$DATE = NULL
subdf.train.ridge.ii <- subdf.train.ridge.ii %>%
  mutate(
    lag_1 = lag(UNRATE, 1),
    lag_2 = lag(UNRATE, 2),
    lag_3 = lag(UNRATE, 3),
    lag_4 = lag(UNRATE, 4),
    lag_12 = lag(UNRATE, 12),
    lag_24 = lag(UNRATE, 24)
  )
subdf.train.ridge.ii = na.omit(subdf.train.ridge.ii) # omit rows with no lags
nsub = nrow(subdf.train.ridge.ii)
subdf.train.ridge.ii.final = subdf.train.ridge.ii[1:(nsub - 1),]
subdf.train.ridge.ii.final$diff = diff(subdf.train.ridge.ii$UNRATE, 1)

lm.model = lm(diff~., data=subdf.train.ridge.ii.final)
summary(lm.model)
mod.mswm = msmFit(lm.model,
                  k=2, # 2 regimes
                  p=0, # number of AR coefficients
                  sw=c(rep(T, 12)))

summary(mod.mswm)
par(mar=c(5.1,4.1,4.1,2.1))
plotProb(mod.mswm,which=1) # plotProb(mod.mswm,which=1)
plotProb(mod.mswm,which=2)

combined_df <- data.frame(XX = as.numeric(rownames(subdf.train.ridge.ii.final)), X = subdf$UNRATE[as.numeric(rownames(subdf.train.ridge.ii.final))], Y = mod.mswm@Fit@filtProb[,1])
ggplot(combined_df, aes(xx = XX, x = X, y = Y)) +
  geom_line(aes(x=XX, y=X), color = "blue") +
  geom_rect(aes(xmin = XX - 0.5, xmax = XX + 0.5, ymin = -Inf, ymax = Inf, fill = Y > 0.1), alpha = 0.3) +
  scale_fill_manual(values = c("TRUE" = "green", "FALSE" = "red")) +
  theme_minimal() +
  labs(title = "Line Graph with Colored Rectangles",
       x = "X Axis",
       y = "Y Axis")
