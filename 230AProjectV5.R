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

### 1.5. Effect of Interest Rate on Unemployment
#        Is interest rate an important covariate?
onedf = subdf
onedf$DATE = NULL

# add in diff 
twodf = onedf[1:(nrow(onedf)-1),]
twodf$diff = diff(onedf$UNRATE, 1)

ggplot() + geom_histogram(colour="black", fill="white", aes(onedf$FEDFUNDS)) # histogram of FEDFUNDS
cor(twodf$diff, twodf$FEDFUNDS) # correlation between UNRATE and FEDFUNDS

# * Regression 1: diff ~ 1 + FEDFUNDS
ols.fit = lm(diff~FEDFUNDS, data=twodf)
summary(ols.fit)
ols.fit.hc0 = sqrt(diag(hccm(ols.fit, type="hc0")))
ols.fit.hc0 # robust beta error

# * Regression 2: diff ~ 1 + FEDFUNDS + GDP + CPI + UNRATE
ols.fit = lm(diff~., data=twodf)
summary(ols.fit)
ols.fit.hc0 = sqrt(diag(hccm(ols.fit, type="hc0")))
ols.fit.hc0 # robust beta error

# adding the lags
onepdf = onedf
onepdf$lag1 = lag(onedf$UNRATE, 1)
onepdf$lag2 = lag(onedf$UNRATE, 2)
onepdf$lag3 = lag(onedf$UNRATE, 3)
onepdf$lag4 = lag(onedf$UNRATE, 4)
onepdf$lag12 = lag(onedf$UNRATE, 12)
onepdf$lag24 = lag(onedf$UNRATE, 24)
one.prime.df = onepdf[1:(nrow(onepdf)-1),]
one.prime.df$diff = diff(onepdf$UNRATE, 1)
one.prime.df = na.omit(one.prime.df)

# * Regression 3: diff ~ 1 + FEDFUNDS + GDP + CPI + UNRATE + lag1 + lag2 + lag3 + lag4 + lag12 + lag24
ols.fit = lm(diff~., data=one.prime.df)
summary(ols.fit)
ols.fit.hc0 = sqrt(diag(hccm(ols.fit, type="hc0")))
ols.fit.hc0 # robust beta error

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

# 2.3 Use AR(4)
arima_model = arima(subdf.train$UNRATE, order = c(4, 0, 0))
summary(arima_model)
arima.coef = coef(arima_model)
arima.pred = numeric(n.test) # AR(4) predictions
for (i in (n.train+1):n){
  arima.pred[i - n.train] = arima.coef%*%c(subdf[(i - 1),]$UNRATE, subdf[(i - 2),]$UNRATE, subdf[(i - 3),]$UNRATE, subdf[(i - 4),]$UNRATE, 0)
}

calculate_mse(subdf.test$UNRATE, arima.pred) # 1.089113
calculate_mse_cutoff(subdf.test$UNRATE, arima.pred) # 0.03116826
ggplot() + geom_line(aes(x=1:n.test,y=arima.pred), color="blue") + 
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

ggplot() + geom_line(aes(x=1:(n.test-1),y=y - y.pred.cumulative), color="blue", linewidth=0.2) # residual plot

### 4. Indicators on past drops

gap = 5
subsubdf = subdf.train[1:(nrow(subdf.train)-1),]
subsubdf$diff = diff(subdf.train$UNRATE, 1)
subsubdf$DATE = NULL
subsubdf.iv = subsubdf[(gap+1):nrow(subsubdf),]
subsubdf.iv$ind = diff(subsubdf$UNRATE, lag=gap) # ind column = UNRATE_t - UNRATE_{t - gap}

summary(subsubdf.iv$ind)

# 4.1 Finding the best gap and threshold lambda
gap_arr = seq(3, 10, 1)
# gap_arr = seq(5, 5, 1)
for (gap in gap_arr){
  
  subsubdf.iv = subsubdf[(gap+1):nrow(subsubdf),]
  subsubdf.iv$ind = diff(subsubdf$UNRATE, lag=gap)
  
  first_quartile = as.numeric(summary(subsubdf.iv$ind)[2])
  third_quartile = as.numeric(summary(subsubdf.iv$ind)[5])
  
  lambda_arr = seq(first_quartile, third_quartile, 0.1) # brute force through -0.7, ..., 0.7
  #lambda_arr = seq(-0.7, 0.7, 0.1)
  
  l = length(lambda_arr)
  res = numeric(l)
  cnt = 1
  for (lambda in lambda_arr){
    # 4.1.1 Separate into two regimes
    regime.one.lambda = subsubdf.iv[subsubdf.iv$ind > lambda,]
    regime.two.lambda = subsubdf.iv[subsubdf.iv$ind <= lambda,]
    
    # 4.1.2 Do separate regression
    y.one = regime.one.lambda$diff
    regime.one.lambda$diff = NULL
    X.one = as.matrix(regime.one.lambda)
    
    cv_model.one = cv.glmnet(X.one, y.one, alpha = 0)
    best_lambda.one = cv_model.one$lambda.min # 0.004143626 (this is pathetic)
    ridge_model_best.one = glmnet(X.one, y.one, alpha = 0, lambda = best_lambda.one)
    
    fitted_values.one = predict(ridge_model_best.one, newx = X.one)
    error.one = calculate_mse(y.one, fitted_values.one)
    
    y.two = regime.two.lambda$diff
    regime.two.lambda$diff = NULL
    X.two = as.matrix(regime.two.lambda)
    
    cv_model.two = cv.glmnet(X.two, y.two, alpha = 0)
    best_lambda.two = cv_model.two$lambda.min # 0.004143626 (this is pathetic)
    ridge_model_best.two = glmnet(X.two, y.two, alpha = 0, lambda = best_lambda.two)
    
    fitted_values.two = predict(ridge_model_best.two, newx = X.two)
    calculate_mse(y.two, fitted_values.two)
    
    error.two = calculate_mse(y.two, fitted_values.two)
    
    error.tot = (error.one*length(y.one) + error.two*length(y.two))/(length(y.one) + length(y.two))
    
    print(paste(lambda, error.tot, sep=" -> "))
    res[cnt] = error.tot
    cnt = cnt + 1
    # 4.1.3 
  }
  idx = which.min(res) # select index of minimum loss
  print(gap)
  print(lambda_arr[idx]) # print corresponding minimum threshold
}

# after analysis (lambda, gap) = (0.1, 5) works well
lambda = 0.1
gap = 5
subsubdf.iv = subsubdf[(gap+1):nrow(subsubdf),]
subsubdf.iv$ind = diff(subsubdf$UNRATE, lag=gap)
regime.one.lambda.star = subsubdf.iv[subsubdf.iv$ind > lambda,]
regime.two.lambda.star = subsubdf.iv[subsubdf.iv$ind <= lambda,]

ggplot() + geom_point(aes(x=as.numeric(rownames(regime.one.lambda.star)),y=regime.one.lambda.star$UNRATE)) +
  geom_point(aes(x=as.numeric(rownames(regime.two.lambda.star)),y=regime.two.lambda.star$UNRATE), color="blue")

y.one.star = regime.one.lambda.star$diff
regime.one.lambda.star$diff = NULL
X.one.star = as.matrix(regime.one.lambda.star)

cv_model.one.star = cv.glmnet(X.one.star, y.one.star, alpha = 0)
best_lambda.one.star = cv_model.one.star$lambda.min # 0.007515883
ridge_model_best.one.star = glmnet(X.one.star, y.one.star, alpha = 0, lambda = best_lambda.one.star)

fitted_values.one.star = predict(ridge_model_best.one.star, newx = X.one.star)
error.one.star = calculate_mse(y.one.star, fitted_values.one.star)
error.one.star # 0.03859857 (apparent error for this section)

y.two.star = regime.two.lambda.star$diff
regime.two.lambda.star$diff = NULL
X.two.star = as.matrix(regime.two.lambda.star)

cv_model.two.star = cv.glmnet(X.two.star, y.two.star, alpha = 0)
best_lambda.two.star = cv_model.two.star$lambda.min # 0.004143626 (this is pathetic)
ridge_model_best.two.star = glmnet(X.two.star, y.two.star, alpha = 0, lambda = best_lambda.two.star)

fitted_values.two.star = predict(ridge_model_best.two.star, newx = X.two.star)
error.two.star = calculate_mse(y.two.star, fitted_values.two.star)
error.two.star # 0.02491733 (apparent error for this section)

ggplot() + geom_point(aes(x=as.numeric(rownames(regime.one.lambda.star)),y=y.one.star)) +
  geom_point(aes(x=as.numeric(rownames(regime.two.lambda.star)),y=y.two.star), color="blue") +
  geom_point(aes(x=as.numeric(rownames(regime.one.lambda.star)),y=fitted_values.one.star), color="green") +
  geom_point(aes(x=as.numeric(rownames(regime.two.lambda.star)),y=fitted_values.two.star), color="red")

mse.up = calculate_mse(y.one.star, fitted_values.one.star)
mse.down = calculate_mse(y.two.star, fitted_values.two.star)

# Here is the train error
mse.tot = (length(y.one.star)*mse.up + length(y.two.star)*mse.down)/(length(y.two.star) + length(y.one.star)) # 0.0292149

indicator.pred = numeric(n.test)
subdf.test.cur = subdf
subdf.test.cur$DATE = NULL
subdf.test.cur$ind = c(rep(NA, gap), diff(subdf.test.cur$UNRATE, lag=gap))
for (i in (n.train+1):n){
  if (subdf.test.cur[i,]$ind > lambda){
    # regime 1
    indicator.pred[i - n.train] = predict(ridge_model_best.one.star, newx = subdf.test.cur[i,])
  } else {
    # regime 2
    indicator.pred[i - n.train] = predict(ridge_model_best.two.star, newx = subdf.test.cur[i,])
  }
}
y.ind.true = subdf.test.cur$UNRATE[(n.train+1):n] - subdf.test.cur$UNRATE[n.train:(n-1)]
calculate_mse(y.ind.true, indicator.pred) # 0.8221356
calculate_mse_cutoff(y.ind.true, indicator.pred) # 0.02453915

ggplot() + geom_line(aes(x=1:n.test,y=subdf.test.cur$UNRATE[n.train:(n-1)]+indicator.pred), color="blue") +
  geom_line(aes(x=1:n.test,y=subdf.test.cur$UNRATE[(n.train+1):n]), color="green") # plot of prediction vs actual

ggplot() + geom_line(aes(x=1:n.test,y=subdf.test.cur$UNRATE[(n.train+1):n] - (subdf.test.cur$UNRATE[n.train:(n-1)]+indicator.pred)), color="black") # residual pplot

ggplot() + geom_line(aes(x=1:75,y=(subdf.test.cur$UNRATE[(n.train+1):n] - (subdf.test.cur$UNRATE[n.train:(n-1)]+indicator.pred))[1:75]), linewidth=0.2,color="black") +
 geom_line(aes(x=1:75,y=c(NA, (y - y.pred.cumulative)[1:74])), linewidth=0.2, color="blue") # residual plot

# Indicator regression beats ridge regression on turning points (?)

### 5. Markov Switching Model

### 5.1
# tried This is from fitting the switching model on the difference. Need more degrees of freedom
subdf.train.v = subdf.train
subdf.train.v$DATE = NULL
subdf.train.v2 = subdf.train.v[1:(nrow(subdf.train.v)-1),]
subdf.train.v2$diff = diff(subdf.train.v$UNRATE, 1)
model = lm(diff~1,data=subdf.train.v2)
mod = msmFit(model, k=2, p=1, sw=rep(TRUE,3))

plotProb(mod, which=3)
plotProb(mod, which=2)
plotProb(mod, which=1)

print(mod@Coef) # coefficients
print(mod@transMat) # translation coefficients
print(mod@iniProb) # initial probability
print(mod@seCoef) # standard error of coefficients
print(mod@Fit@filtProb) # filtered probability
# print(mod@Fit@smoTransMat) # transition matrix array

state.here = numeric(nrow(mod@Fit@filtProb))
for(i in 1:nrow(mod@Fit@filtProb)){
  if (mod@Fit@filtProb[i,1]>mod@Fit@filtProb[i,2]){
    state.here[i] = 1
  } else {
    state.here[i] = 2
  }
}

ggplot() + geom_line(aes(x=1:n.train,y=subdf.train$UNRATE)) +
  geom_line(aes(x=1:n.train,y=c(rep(NA, 2), state.here)),color="orange")+
  labs(x = "Months since 1954 Jul", y = "Unemployment Rate")+
  theme_minimal()

subdf.v = subdf
subdf.v$DATE = NULL
subdf.v2 = subdf.train.v[1:(nrow(subdf.v)-1),]
subdf.v2$diff = diff(subdf.v$UNRATE, 1)
pred.v = numeric(n.test)
print(nrow(subdf.v2) - n.train)
count = 0
for (i in n.train:nrow(subdf.v2)){
  count = count + 1
  if (count > 50){
    break
  }
  print(count)
  model = lm(diff~1,data=subdf.v2[1:(i-1),])
  mod = msmFit(model, k=2, p=1, sw=rep(TRUE,3))
  nnn = nrow(mod@Fit@filtProb)
  if (mod@Fit@filtProb[nnn,1]>mod@Fit@filtProb[nnn,2]){
    pred.v[i - n.train + 1] = as.numeric(as.numeric(mod@Coef[1,])%*%c(1, subdf.v2[(i-1),]$diff))
  } else {
    pred.v[i - n.train + 1] = as.numeric(as.numeric(mod@Coef[2,])%*%c(1, subdf.v2[(i-1),]$diff))
  }
}

ggplot() + geom_line(aes(x=1:n.train,y=subdf.train$UNRATE)) +
  geom_line(aes(x=1:n.train,y=c(rep(NA, 2), state.here)),color="orange")+
  geom_line(aes(x=(n.train+1):(n.train+50),y=pred.v[1:50] + subdf.v[n.train:(n.train+49),]$UNRATE),color="blue")+
  geom_line(aes(x=(n.train+1):(n.train+50),y=subdf.v[(n.train+1):(n.train+50),]$UNRATE),color="green")+
  labs(x = "Months since 1954 Jul", y = "Unemployment Rate")+
  theme_minimal()

calculate_mse_cutoff(subdf.v[(n.train+1):(n.train+50),]$UNRATE, pred.v[1:50] + subdf.v[n.train:(n.train+49),]$UNRATE) # 0.02454834 wow!

### 5.2
subdf.train.v = subdf.train
subdf.train.v$DATE = NULL
subdf.train.v2 = subdf.train.v[1:(nrow(subdf.train.v)-1),]
subdf.train.v2$diff = diff(subdf.train.v$UNRATE, 1)
model = lm(diff~1,data=subdf.train.v2)
mod = msmFit(model, k=2, p=4, sw=rep(TRUE,6))

plotProb(mod, which=3)
plotProb(mod, which=2)
plotProb(mod, which=1)

print(mod@Coef) # coefficients
print(mod@transMat) # translation coefficients
print(mod@iniProb) # initial probability
print(mod@seCoef) # standard error of coefficients
# print(mod@Fit@smoTransMat) # transition matrix array

state.here = numeric(nrow(mod@Fit@filtProb))
for(i in 1:nrow(mod@Fit@filtProb)){
  if (mod@Fit@filtProb[i,1]>mod@Fit@filtProb[i,2]){
    state.here[i] = 1
  } else {
    state.here[i] = 2
  }
}

ggplot() + geom_line(aes(x=1:n.train,y=subdf.train$UNRATE)) +
  geom_line(aes(x=1:n.train,y=c(rep(NA, 5), state.here)),color="orange")+
  labs(x = "Months since 1954 Jul", y = "Unemployment Rate")+
  theme_minimal()

subdf.v = subdf
subdf.v$DATE = NULL
subdf.v2 = subdf.train.v[1:(nrow(subdf.v)-1),]
subdf.v2$diff = diff(subdf.v$UNRATE, 1)
pred.v = numeric(n.test)
print(nrow(subdf.v2) - n.train)
count = 0
for (i in n.train:nrow(subdf.v2)){
  count = count + 1
  if (count > 50){
    break
  }
  print(count)
  model = lm(diff~1,data=subdf.v2[1:(i-1),])
  mod = msmFit(model, k=2, p=4, sw=rep(TRUE,6))
  nnn = nrow(mod@Fit@filtProb)
  if (mod@Fit@filtProb[nnn,1]>mod@Fit@filtProb[nnn,2]){
    pred.v[i - n.train + 1] = as.numeric(as.numeric(mod@Coef[1,])%*%c(1, subdf.v2[(i-1),]$diff, subdf.v2[(i-2),]$diff, subdf.v2[(i-3),]$diff, subdf.v2[(i-4),]$diff))
  } else {
    pred.v[i - n.train + 1] = as.numeric(as.numeric(mod@Coef[2,])%*%c(1, subdf.v2[(i-1),]$diff, subdf.v2[(i-2),]$diff, subdf.v2[(i-3),]$diff, subdf.v2[(i-4),]$diff))
  }
}

ggplot() + geom_line(aes(x=1:n.train,y=subdf.train$UNRATE)) +
  geom_line(aes(x=1:n.train,y=c(rep(NA, 5), state.here)),color="orange")+
  geom_line(aes(x=(n.train+1):(n.train+50),y=pred.v[1:50] + subdf.v[n.train:(n.train+49),]$UNRATE),color="blue")+
  geom_line(aes(x=(n.train+1):(n.train+50),y=subdf.v[(n.train+1):(n.train+50),]$UNRATE),color="green")+
  labs(x = "Months since 1954 Jul", y = "Unemployment Rate")+
  theme_minimal()

calculate_mse_cutoff(subdf.v[(n.train+1):(n.train+50),]$UNRATE, pred.v[1:50] + subdf.v[n.train:(n.train+49),]$UNRATE) # 0.02313789 wow wow!


### 5.3
subdf.v = subdf.train
subdf.v$DATE = NULL
subdf.v2 = subdf.v[1:(nrow(subdf.v)-1),]
subdf.v2$diff = diff(subdf.v$UNRATE, 1)
model = lm(diff~.,data=subdf.v2)
mod = msmFit(model, k=2, p=4, sw=rep(TRUE,10))

plotProb(mod, which=3)
plotProb(mod, which=2)
plotProb(mod, which=1)

print(mod@Coef) # coefficients
print(mod@transMat) # translation coefficients
print(mod@iniProb) # initial probability
print(mod@seCoef) # standard error of coefficients
# print(mod@Fit@smoTransMat) # transition matrix array

state.here = numeric(nrow(mod@Fit@filtProb))
for(i in 1:nrow(mod@Fit@filtProb)){
  if (mod@Fit@filtProb[i,1]>mod@Fit@filtProb[i,2]){
    state.here[i] = 1
  } else {
    state.here[i] = 2
  }
}

ggplot() + geom_line(aes(x=1:n.train,y=subdf.train$UNRATE)) +
  geom_line(aes(x=1:n.train,y=c(rep(NA, 5), state.here)),color="orange")+
  labs(x = "Months since 1954 Jul", y = "Unemployment Rate")+
  theme_minimal()

subdf.v = subdf
subdf.v$DATE = NULL
subdf.v2 = subdf.v[1:(nrow(subdf.v)-1),]
subdf.v2$diff = diff(subdf.v$UNRATE, 1)
pred.v = numeric(n.test)
print(nrow(subdf.v2) - n.train)
count = 0
for (i in n.train:nrow(subdf.v2)){
  count = count + 1
  if (count > 50){
    break
  }
  print(count)
  model = lm(diff~.,data=subdf.v2[1:(i-1),])
  mod = msmFit(model, k=2, p=4, sw=rep(TRUE,10))
  nnn = nrow(mod@Fit@filtProb)
  if (mod@Fit@filtProb[nnn,1]>mod@Fit@filtProb[nnn,2]){
    pred.v[i - n.train + 1] = as.numeric(as.numeric(mod@Coef[1,])%*%c(1, subdf.v2[i,]$UNRATE, subdf.v2[i,]$FEDFUNDS, subdf.v2[i,]$GDP, subdf.v2[i,]$CPIAUCSL, subdf.v2[(i-1),]$diff, subdf.v2[(i-2),]$diff, subdf.v2[(i-3),]$diff, subdf.v2[(i-4),]$diff))
  } else {
    pred.v[i - n.train + 1] = as.numeric(as.numeric(mod@Coef[2,])%*%c(1, subdf.v2[i,]$UNRATE, subdf.v2[i,]$FEDFUNDS, subdf.v2[i,]$GDP, subdf.v2[i,]$CPIAUCSL, subdf.v2[(i-1),]$diff, subdf.v2[(i-2),]$diff, subdf.v2[(i-3),]$diff, subdf.v2[(i-4),]$diff))
  }
}

ggplot() + geom_line(aes(x=1:n.train,y=subdf.train$UNRATE)) +
  geom_line(aes(x=1:n.train,y=c(rep(NA, 5), state.here)),color="orange")+
  geom_line(aes(x=(n.train+1):(n.train+50),y=pred.v[1:50] + subdf.v[n.train:(n.train+49),]$UNRATE),color="blue")+
  geom_line(aes(x=(n.train+1):(n.train+50),y=subdf.v[(n.train+1):(n.train+50),]$UNRATE),color="green")+
  labs(x = "Months since 1954 Jul", y = "Unemployment Rate")+
  theme_minimal()

calculate_mse_cutoff(subdf.v[(n.train+1):(n.train+50),]$UNRATE, pred.v[1:50] + subdf.v[n.train:(n.train+49),]$UNRATE) # 0.02397587 wow wow!



### 4. Indicators on past changes

# 4a. Grid Search
gap_arr = seq(3, 10, 1)
for (gap in gap_arr){
  
  subsubdf.iv = subsubdf[(gap+1):nrow(subsubdf),]
  subsubdf.iv$ind = diff(subsubdf$UNRATE, lag=gap)
  
  first_quartile = as.numeric(summary(subsubdf.iv$ind)[2])
  third_quartile = as.numeric(summary(subsubdf.iv$ind)[5])
  
  lambda_arr = seq(first_quartile, third_quartile, 0.1)

  l = length(lambda_arr)
  res = numeric(l)
  cnt = 1
  for (lambda in lambda_arr){
    # 4.1 Separate into two regimes
    regime.one.lambda = subsubdf.iv[subsubdf.iv$ind > lambda,]
    regime.two.lambda = subsubdf.iv[subsubdf.iv$ind <= lambda,]
    
    # 4.2 Do separate regression
    y.one = regime.one.lambda$diff
    regime.one.lambda$diff = NULL
    X.one = as.matrix(regime.one.lambda)
    
    cv_model.one = cv.glmnet(X.one, y.one, alpha = 0)
    best_lambda.one = cv_model$lambda.min # 0.004143626 (this is pathetic)
    ridge_model_best.one = glmnet(X.one, y.one, alpha = 0, lambda = best_lambda.one)
    
    fitted_values.one = predict(ridge_model_best.one, newx = X.one)
    error.one = calculate_mse(y.one, fitted_values.one)
    
    y.two = regime.two.lambda$diff
    regime.two.lambda$diff = NULL
    X.two = as.matrix(regime.two.lambda)
    
    cv_model.two = cv.glmnet(X.two, y.two, alpha = 0)
    best_lambda.two = cv_model$lambda.min # 0.004143626 (this is pathetic)
    ridge_model_best.two = glmnet(X.two, y.two, alpha = 0, lambda = best_lambda.two)
    
    fitted_values.two = predict(ridge_model_best.two, newx = X.two)
    calculate_mse(y.two, fitted_values.two)
    
    error.two = calculate_mse(y.two, fitted_values.two)
    
    error.tot = (error.one*length(y.one) + error.two*length(y.two))/(length(y.one) + length(y.two))
    
    print(paste(lambda, error.tot, sep=" -> "))
    res[cnt] = error.tot
    cnt = cnt + 1
    # 4.3 
  }
  idx = which.min(res)
  print(gap)
  print(lambda_arr[idx])
}

# 4b. Optimal
lambda = 0.1
gap = 5
regime.one.lambda.star = subsubdf.iv[subsubdf.iv$ind > lambda,]
regime.two.lambda.star = subsubdf.iv[subsubdf.iv$ind <= lambda,]

ggplot() + geom_point(aes(x=as.numeric(rownames(regime.one.lambda.star)),y=regime.one.lambda.star$UNRATE)) +
  geom_point(aes(x=as.numeric(rownames(regime.two.lambda.star)),y=regime.two.lambda.star$UNRATE), color="blue")

### 4.2 Model Selection

### 4.3 Prediction

### 5. Autoregressive Switching Model

# 5.1 Switching model on the unemployment data itself
subdf.train.v = subdf.train
subdf.train.v$DATE = NULL
model = lm(UNRATE~1,data=subdf.train.v)
mod = msmFit(model, k=2, p=2, sw=rep(TRUE,4))

plotProb(mod, which=3)
plotProb(mod, which=2)
plotProb(mod, which=1)

print(mod@Coef) # coefficients
print(mod@transMat) # translation coefficients
print(mod@iniProb) # initial probability
print(mod@seCoef) # standard error of coefficients
print(mod@Fit@smoTransMat) # transition matrix array

# 5.2 More autoregressive coefficients
subdf.train.v = subdf.train
subdf.train.v$DATE = NULL
model = lm(UNRATE~1,data=subdf.train.v)
mod = msmFit(model, k=2, p=4, sw=rep(TRUE,6))

plotProb(mod, which=3)
plotProb(mod, which=2)
plotProb(mod, which=1)

print(mod@Coef) # coefficients

# 5.3 Operating on diff
subdf.train.v = subdf.train
subdf.train.v$DATE = NULL
subdf.train.v2 = subdf.train.v[1:(nrow(subdf.train.v)-1),]
subdf.train.v2$diff = diff(subdf.train.v$UNRATE, 1)
model = lm(diff~1,data=subdf.train.v2)
mod = msmFit(model, k=2, p=1, sw=rep(TRUE,3))

plotProb(mod, which=3)
plotProb(mod, which=1)

# 5.4 Operating on diff
subdf.train.v = subdf.train
subdf.train.v$DATE = NULL
subdf.train.v2 = subdf.train.v[1:(nrow(subdf.train.v)-1),]
subdf.train.v2$diff = diff(subdf.train.v$UNRATE, 1)
model = lm(diff~1,data=subdf.train.v2)
mod = msmFit(model, k=2, p=3, sw=rep(TRUE,5))

plotProb(mod, which=3)
plotProb(mod, which=1)
