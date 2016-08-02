# Script to generate data in package pcev
logit <- function(t) log(t) - log(1-t)
BLKData <- read.csv("data-raw/BLKData.csv", header = TRUE)

BLKmethyl <- data.frame("position" = BLKData$position,
                        BLKData[,grep("read", names(BLKData))]/
                          BLKData[,grep("coverage", names(BLKData))])

# Remove duplicated rows
BLKmethyl <- BLKmethyl[!duplicated(BLKmethyl[,-1]),]
# Remove 1s and 0s
BLKmethyl[BLKmethyl == 0] <- 0.01
BLKmethyl[BLKmethyl == 1] <- 0.99
# Remove NAs, which correspond to 0 coverage
naRows <- apply(BLKmethyl[,-1], 1, function(row) sum(is.na(row)) > 0)
BLKmethyl <- BLKmethyl[!naRows,]

variances <- apply(logit(as.matrix(BLKmethyl[-1,])), 1, var, na.rm=TRUE)

# Select top B most variable
B = 1000
BLKmethyl <- BLKmethyl[order(variances, decreasing = TRUE)[1:B],]
BLKmethyl <- BLKmethyl[order(BLKmethyl$position),]

# Create data for use in package
methylation2 <- as.matrix(BLKmethyl[,-1])
position2 <- BLKmethyl$position
cell_type <- 1*(grepl("_BC_", names(BLKmethyl[,-1]))) +
  2*(grepl("_TC_", names(BLKmethyl[,-1]))) + 
  3*(grepl("_Mono_", names(BLKmethyl[,-1])))
pheno2 <- c("BC", "TC", "Mono")[cell_type]

devtools::use_data(methylation2, position2, pheno2)