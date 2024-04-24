library(quantmod)

ticker = "SPY"

start_date = as.Date("1947-12-31")
end_date = Sys.Date() 

getSymbols(ticker, src = "yahoo", from = start_date, to = end_date, periodicity = "monthly")

write.csv(SPY, "SPY.csv", row.names = T)
