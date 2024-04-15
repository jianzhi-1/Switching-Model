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

### 0. Set up
set.seed(42)
data <- read.csv("UNRATE.csv")
df = as.data.frame(data)
n = length(df$UNRATE) # 915
n.train = 450
# ggplot() + geom_line(aes(x=1:n),y=df$UNRATE)
df.train = df[1:n.train,]
ggplot() + geom_line(aes(x=1:n.train,y=df.train$UNRATE)) # plot training data

df.train.diff = diff(df.train$UNRATE)
ggplot() + geom_line(aes(x=1:length(df.train.diff),y=df.train.diff)) # plot difference array

ggplot() + geom_line(aes(x=1:850,y=diff[1:850]))

nn = length(diff)

### 1. Gaussian Update Formula
