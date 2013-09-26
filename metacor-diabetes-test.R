setwd("~/Dropbox/Datasets/Implicit bias meta analysis")
source("metacor.r")
source("remlcor.r")
library(metaSEM)

# Read in the data file
d <- read.csv("diabetes4.csv", header = T, row.names = 1)

# Replace missing values with numbers indicating a small amount of information
d[, paste0("r", 1:6)] <- lapply(d[, paste0("r", 1:6)], 
                                function(i) ifelse(is.na(i), .0001, i))
d[, paste0("n", 1:6)] <- lapply(d[, paste0("n", 1:6)], 
                                function(i) ifelse(is.na(i), .001, i))

# Create trt, which will contain the calculated study effect sizes
trt <- matrix(nrow = dim(d)[1], ncol = 5)
colnames(trt) <- paste0("trt", 2:6)

# Create helper matrices that will assist in making the column names of S
mat1 <- matrix(rep(2:6, 5), 5, 5, byrow = T)
mat2 <- matrix(rep(2:6, 5), 5, 5)
# Create S, the matrix containing the within-studies variance-covariance matrices
S <- matrix(nrow = dim(d)[1], ncol = (5 * 6) / 2)
colnames(S) <- paste0(vech(mat1), vech(mat2))

# For every treatment besides 1, the reference
for (i in 2:6)
{
  # Get the reference r and n
  r1 <- d[, "r1"]
  n1 <- d[, "n1"]
  
  # Get the information for the current treatment
  trt_name <- paste0("trt", i)
  r_cur_1 <- d[, paste0("r", i)]
  n_cur_1 <- d[, paste0("n", i)]
  
  # Calculate the effect size and variance for the current comparison between
  # the reference and the current treatment
  trt[, trt_name] <- log(r_cur_1 / (n_cur_1 - r_cur_1)) - log(r1 / (n1 - r1))
  S[, paste0(i, i)] <- 1/r_cur_1 + 1/(n_cur_1 - r_cur_1) + 1/r1 + 1/(n1 - r1)
  
  # If we're not on the last treatment
  if (i != 6)
  {
    # For all the remaining treatments, calculate the covariances
    for (j in (i + 1):6)
    {
      r_cur_2 <- d[, paste0("r", j)]
      n_cur_2 <- d[, paste0("n", j)]
      
      S[, paste0(i, j)] <- 1/r1 + 1/(n1 - r1)
    } 
  }
}

con <- diag(rep("0*a", 5))
cor <- matrix(.5, nrow = 5, ncol = 5)
diag(cor) <- 1

mod <- metacor(trt, S, RE.constraints = con, cor.constraints = cor)
summary(mod)

mod.reml <- reml(trt, S)
summary(mod.reml)
varComp <- diag(c(coef(mod.reml), coef(mod.reml)))
mod <- metacor(y, v, RE.constraints = varComp)
summary(mod)