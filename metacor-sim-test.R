setwd("~/Dropbox/Datasets/Implicit bias meta analysis")
source("metacor.r")
library(metaSEM)


# Create data-frame of observed
AB <- c(.5, .3, .2, .2, NA, NA, .1)
AC <- c(NA, .3, .4, NA, .2, .1, .5)
y <- data.frame(AB, AC)

n1 <- 50
n2 <- 50
nT <- 150 # Total N for 3-arm studies in this example

# Create matrix of observed variances of the study effect sizes and covariances between effect sizes.
# These are calculated using the formulas for the variance and covariance of a within-studies effect size
# from Gleser and Olkin (1994) assuming that each study arm has an N of 50 (and therefore all 2-arm
# studies have a total N of 100 and all 3-arm studies have an N of 150)
v <- matrix(c(1/n1 + 1/n2 + .5^2/(2*(n1 + n2)), NA, NA,
              1/n1 + 1/n2 + .3^2/(2*nT), 1/n1 + .3 * .3 / (2*nT), 1/n1 + 1/n2 + .3^2/(2*nT),
              1/n1 + 1/n2 + .2^2/(2*nT), 1/n1 + .2 * .4 / (2*nT), 1/n1 + 1/n2 + .4^2/(2*nT),
              1/n1 + 1/n2 + .2^2/(2*(n1 + n2)), NA, NA,
              NA, NA, 1/n1 + 1/n2 + .2^2/(2*(n1 + n2)),
              NA, NA, 1/n1 + 1/n2 + .1^2/(2*(n1 + n2)),
              1/n1 + 1/n2 + .1^2/(2*nT), 1/n1 + .1 * .5 / (2*nT), 1/n1 + 1/n2 + .5^2/(2*nT)), 
            nrow = 7, ncol = 3, byrow = T)

# Testing: no random effect constraints
mod <- metacor(y, v)
summary(mod)
# The covariance value is negative

# Testing: equal variances, separate covariance
con <- matrix(c("0*a", "0*b",
                "0*b", "0*a"), nrow = 2, ncol = 2, byrow = T)

mod <- metacor(y, v, RE.constraints = con)
summary(mod)
# Covariance value still negative

con <- diag(c("0*a", "0*a"))
cor <- matrix(c(1, .5,
                .5, 1), nrow = 2, ncol = 2, byrow = T)

mod <- metacor(y, v, RE.constraints = con, cor.constraints = cor)
summary(mod)


# Create a second set of data
AB <- rnorm(30, mean = .3, sd = .3)
AC <- rnorm(30, mean = .1, sd = .3)
AD <- rnorm(30, mean = -.2, sd = .3)
AE <- rnorm(30, mean = 0, sd = .3)

# Distribute NAs
AB[sample(1:30, size = 5)] <- NA
AC[sample(1:30, size = 5)] <- NA
AD[sample(1:30, size = 5)] <- NA
AE[sample(1:30, size = 5)] <- NA

# Create data frame
y <- data.frame(AB, AC, AD, AE)
ns <- data.frame(ABn = ifelse(is.na(AB), NA, 50), 
                 ACn = ifelse(is.na(AC), NA, 50), 
                 ADn = ifelse(is.na(AD), NA, 50), 
                 AEn = ifelse(is.na(AE), NA, 50))

# Create variances
ABvar <- ifelse(!is.na(y$AB), 2 * (1/ns$ABn) + y$AB^2 / (2 * rowSums(ns, na.rm = T)), NA)
ACvar <- ifelse(!is.na(y$AC), 2 * (1/ns$ACn) + y$AC^2 / (2 * rowSums(ns, na.rm = T)), NA)
ADvar <- ifelse(!is.na(y$AD), 2 * (1/ns$ADn) + y$AD^2 / (2 * rowSums(ns, na.rm = T)), NA)
AEvar <- ifelse(!is.na(y$AE), 2 * (1/ns$AEn) + y$AE^2 / (2 * rowSums(ns, na.rm = T)), NA)
AB_ACcov <- ifelse(!is.na(y$AB) & !is.na(y$AC), 1/ns$ABn + y$AB * y$AC / (2 * rowSums(ns, na.rm = T)), NA)
AB_ADcov <- ifelse(!is.na(y$AB) & !is.na(y$AD), 1/ns$ABn + y$AB * y$AD / (2 * rowSums(ns, na.rm = T)), NA)
AB_AEcov <- ifelse(!is.na(y$AB) & !is.na(y$AE), 1/ns$ABn + y$AB * y$AE / (2 * rowSums(ns, na.rm = T)), NA)
AC_ADcov <- ifelse(!is.na(y$AC) & !is.na(y$AD), 1/ns$ACn + y$AC * y$AD / (2 * rowSums(ns, na.rm = T)), NA)
AC_AEcov <- ifelse(!is.na(y$AC) & !is.na(y$AE), 1/ns$ACn + y$AC * y$AE / (2 * rowSums(ns, na.rm = T)), NA)
AD_AEcov <- ifelse(!is.na(y$AD) & !is.na(y$AE), 1/ns$ADn + y$AD * y$AE / (2 * rowSums(ns, na.rm = T)), NA)

v <- data.frame(ABvar, AB_ACcov, AB_ADcov, AB_AEcov, ACvar, AC_ADcov, AC_AEcov, ADvar, AD_AEcov, AEvar)

# First test: No constraints on the between-studies variance-covariance matrix
mod <- metacor(y, v)
summary(mod)

# Second test: Variances constrained to be equal, covariances constrained to be equal
con <- matrix(c("0*a", "0*b", "0*b", "0*b",
                "0*b", "0*a", "0*b", "0*b",
                "0*b", "0*b", "0*a", "0*b",
                "0*b", "0*b", "0*b", "0*a"), nrow = 4, ncol = 4, byrow = T)

mod <- metacor(y, v, RE.constraints = con)
summary(mod)

# Create correlation constraints
con <- matrix(0, nrow = 4, ncol = 4)
diag(con) <- paste0("0*", letters[1:4])

# All correlations constrained to equal .5
cor <- matrix(.5, nrow = 4, ncol = 4)
diag(cor) <- 1

# Test the function
mod <- metacor(y, v, RE.constraints = con, cor.constraints = cor)
summary(mod)