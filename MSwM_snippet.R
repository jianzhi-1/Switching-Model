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
