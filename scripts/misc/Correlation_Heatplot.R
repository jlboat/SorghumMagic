library(metan)
library(readxl)

df <- read_xlsx("/Users/nkumar/Desktop/file.xlsx")
head(df)

df[,1] <- lapply(df[,1], function(x) as.factor(x))
df[,2:13] <- lapply(df [,2:13], function(x) as.numeric(x))
str(df)

corrl <-corr_coef(df)
plot(corrl)
print(corrl)
sink("corrl.txt")



