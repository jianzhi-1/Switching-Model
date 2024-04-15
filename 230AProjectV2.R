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

### 0. Set up
set.seed(42)
data <- read.csv("UNRATE.csv")
df = as.data.frame(data)
n = length(df$UNRATE) # 915
n.train = 450
n.test = n - n.train
# ggplot() + geom_line(aes(x=1:n),y=df$UNRATE)
df.train = df[1:n.train,]
df.test = df[(n.train+1):n,]
ggplot() + geom_line(aes(x=1:n.train,y=df.train$UNRATE)) # plot training data

df.train.diff = diff(df.train$UNRATE)
ggplot() + geom_line(aes(x=1:length(df.train.diff),y=df.train.diff)) # plot difference array

### 1. Best Linear Approximation, Rolling Window Version

calculate_mse <- function(observed_values, predicted_values) {
  squared_diff <- (observed_values - predicted_values)^2
  mse <- mean(squared_diff)
  return(mse)
}

max.window.size = 10
window.size.seq = seq(2, max.window.size, 1)
mse.arr = numeric(length(window.size.seq))

mse.min = NULL
pred.best = numeric(n.test)

for (window.size in window.size.seq){
  pred = numeric(n.test)
  for (i in 1:n.test){
    x = (n.train + i - window.size):(n.train - 1 + i)
    y = df[x,]$UNRATE
    lm.model = lm(y~x)
    
    pred[i] = predict(lm.model, newdata=data.frame(x=(n.train+i):(n.train+i)))
  }
  #ggplot() + geom_line(aes(x=1:length(pred),y=pred), color="blue", linewidth=0.3) +
  #  geom_line(aes(x=1:n.test,y=df.test$UNRATE), color="green", linewidth=0.3) # plot difference array
  mse.arr[window.size] = calculate_mse(pred, df.test$UNRATE)
  if (is.null(mse.min) || mse.min > mse.arr[window.size]){
    mse.min = mse.arr[window.size]
    pred.best = pred
  }
}

ggplot() + geom_line(aes(x=window.size.seq,y=mse.arr[2:length(mse.arr)]), color="blue", linewidth=0.3) 
ggplot() + geom_line(aes(x=1:length(pred.best),y=pred.best), color="blue", linewidth=0.3) +
  geom_line(aes(x=1:n.test,y=df.test$UNRATE), color="green", linewidth=0.3) # plot difference array
ggplot() + geom_line(aes(x=1:length(pred.best),y=pred.best - df.test$UNRATE), color="red", linewidth=0.3) # residual plot

# 4 seems to be the best rolling window size
# MSE = 0.4671989
# [1] 0.5617849 0.4869654 0.4671989 0.4897314 0.5131216 0.5407904
# [7] 0.5649271 0.5839810 0.5994566

# 2. Improving the model by adding ARIMA

acf(pred.best - df.test$UNRATE) # seems like can model as a MA(4) model
pacf(pred.best - df.test$UNRATE)

# 3. Is there an adaptive way of selecting rolling window? Is 4 always the best?
smoother = function(val, lam = 0.25){
  ret = numeric(length(val))
  for (i in 1:length(val)){
    if (i == 1){
      ret[i] = (1 - lam)*val[i] + lam*val[i + 1]
    } else if (i == length(val)){
      ret[i] = (1 - lam)*val[i] + lam*val[i - 1]
    } else {
      ret[i] = (1 - 2*lam)*val[i] + lam*val[i - 1] + lam*val[i + 1]
    }
  }
  return(ret)
}
# plot the best rolling window size: which provided the best estimate
max.window.size = 20
window.size.seq = seq(2, max.window.size, 1)
nn = length((max.window.size + 1):n.train)
res = numeric(nn)
best.l.arr = numeric(n.train)

X.feature = matrix(nrow = nn, ncol = max.window.size)
y.feature = numeric(nn)

X.feature.smooth = matrix(nrow = nn, ncol = max.window.size)

for (i in (max.window.size + 1):n.train){
  # Task 4 Preparation
  X.feature[(i - max.window.size),] = df$UNRATE[(i - max.window.size):(i-1)]
  X.feature.smooth[(i - max.window.size),] = smoother(df$UNRATE[(i - max.window.size):(i - 1)])
  
  l.min.mse = NULL
  for (l in window.size.seq){
    x = (i - l):(i - 1) # size l window size
    y = df[x,]$UNRATE
    lm.model = lm(y~x)
    cur.pred = predict(lm.model, newdata=data.frame(x=i:i))
    if (is.null(l.min.mse) || l.min.mse > (cur.pred - df[i,]$UNRATE)^2){
      l.min.mse = (cur.pred - df[i,]$UNRATE)^2
      best.l.arr[i] = l
    }
  }
  
  # Task 4 Preparation
  y.feature[(i - max.window.size)] = best.l.arr[i]
}

ggplot() + geom_histogram(colour="black", fill="white", aes(best.l.arr[(max.window.size + 1):n.train])) # concentrate efforts on window size 15 and stuff

# Take the previous max.window.size features and 

#
#logistic.model <- glm(y.feature ~ X.feature, family = binomial(link = "logit"))
categorical.df = data.frame(y=factor(y.feature), X.feature)

### 4.1 Multinomial Regression
multinomial.regression.model <- multinom(y ~ ., data=categorical.df)

pred.iv = numeric(n)
for (i in (n.train+1):n){
  X.feature.cur = df$UNRATE[(i - max.window.size):(i-1)]
  X.feature.cur = as.data.frame(t(X.feature.cur))
  colnames(X.feature.cur) = colnames(categorical.df)[2:length(colnames(categorical.df))]
  wl = predict(multinomial.regression.model, newdata=X.feature.cur, type = "class") # rolling window size
  wl = 1 + as.numeric(wl)
  print(wl)
  X.feature.new = (i - wl):(i - 1)
  y.new = df[X.feature.new,]$UNRATE
  lm.model = lm(y.new~X.feature.new)
  pred.df = data.frame(x=(i:i))
  colnames(pred.df) = "X.feature.new"
  pred.iv[i] = predict(lm.model, newdata=pred.df)
}

ggplot() + geom_line(aes(x=(n.train+1):n,y=pred.iv[(n.train+1):n]), color="blue", linewidth=0.3) +
  geom_line(aes(x=(n.train+1):n,y=df.test$UNRATE), color="green", linewidth=0.3) # plot difference array
calculate_mse(pred.iv[(n.train+1):n], df.test$UNRATE) #0.5951225


### 4.2 Ordinal Regression
ordinal.regression.model <- polr(y ~ ., data = categorical.df )

ordinal.df.smooth = data.frame(y=factor(y.feature), X.feature.smooth)
ordinal.regression.model.smooth <- polr(y ~ ., data = ordinal.df.smooth)

#ordinal.regression.model <- polr(factor(y.feature) ~ as.data.frame(X.feature))

pred.iv = numeric(n)
for (i in (n.train+1):n){
  X.feature.cur = df$UNRATE[(i - max.window.size):(i-1)]
  X.feature.cur = as.data.frame(t(X.feature.cur))
  colnames(X.feature.cur) = colnames(ordinal.df)[2:length(colnames(ordinal.df))]
  wl = predict(ordinal.regression.model, newdata=X.feature.cur, type = "class") # rolling window size
  wl = 1 + as.numeric(wl)
  print(wl)
  X.feature.new = (i - wl):(i - 1)
  y.new = df[X.feature.new,]$UNRATE
  lm.model = lm(y.new~X.feature.new)
  pred.df = data.frame(x=(i:i))
  colnames(pred.df) = "X.feature.new"
  pred.iv[i] = predict(lm.model, newdata=pred.df)
}

# Most of the ordinal predicted 2 i.e. just use two values

ggplot() + geom_line(aes(x=(n.train+1):n,y=pred.iv[(n.train+1):n]), color="blue", linewidth=0.3) +
  geom_line(aes(x=(n.train+1):n,y=df.test$UNRATE), color="green", linewidth=0.3) # plot difference array
calculate_mse(pred.iv[(n.train+1):n], df.test$UNRATE) #0.6123597

#ggplot() + geom_line(aes(x=1:length(pred.best),y=pred.best - df.test$UNRATE), color="red", linewidth=0.3) # residual plot

pred.iv = numeric(n)
for (i in (n.train+1):n){
  X.feature.cur = smoother(df$UNRATE[(i - max.window.size):(i-1)])
  X.feature.cur = as.data.frame(t(X.feature.cur))
  colnames(X.feature.cur) = colnames(ordinal.df)[2:length(colnames(ordinal.df))]
  wl = predict(ordinal.regression.model.smooth, newdata=X.feature.cur, type = "class") # rolling window size
  wl = 1 + as.numeric(wl)
  X.feature.new = (i - wl):(i - 1)
  y.new = df[X.feature.new,]$UNRATE
  lm.model = lm(y.new~X.feature.new)
  pred.df = data.frame(x=(i:i))
  colnames(pred.df) = "X.feature.new"
  pred.iv[i] = predict(lm.model, newdata=pred.df)
}

# Most of the ordinal predicted 2 i.e. just use two values

ggplot() + geom_line(aes(x=(n.train+1):n,y=pred.iv[(n.train+1):n]), color="blue", linewidth=0.3) +
  geom_line(aes(x=(n.train+1):n,y=df.test$UNRATE), color="green", linewidth=0.3) # plot difference array
calculate_mse(pred.iv[(n.train+1):n], df.test$UNRATE) #0.6123597 no difference
