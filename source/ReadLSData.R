install.packages('readr')
library(readr)
filename <- '/Users/oeztuerk/Documents/CaTeNa/Susceptibility/RobertLSData.csv'
T <- read_csv(filename)
dim(T)
sort(unique(T$Geology))
