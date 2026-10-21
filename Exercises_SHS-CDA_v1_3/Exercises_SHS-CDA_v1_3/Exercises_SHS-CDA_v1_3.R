Sys.setlocale("LC_TIME", "C")
# date_update <- format(Sys.time(), "%A, %B %e, %Y at %I:%M %p")  # AM/PM
date_update <- format(Sys.time(), "%A, %B %e, %Y at %H:%M:%S")  # 24H





# install.packages("mice") # if you need to install the 'mice' package
library(mice)
data(nhanes2)
nhanes2
summary(nhanes2)

# ?nhanes2

getBMIgroup <- function(bmi) {  # numerical BMI
  ### categorization of BMI according to WHO:
  # underweight if BMI < 18.5
  # normal      if BMI in [18.5, 25)
  # overweight  if BMI in [25, 30)
  # obese       if BMI >= 30
  bmibreaks <- c(0, 18.5, 25, 30, Inf)
  bmilevels <- c("underweight", "normal", "overweight", "obese")
  res <- cut(bmi,
             breaks = bmibreaks,
             labels = bmilevels,
             include.lowest = TRUE,
             right = FALSE)
  res <- droplevels(res)
  return(res)
}

nhanes2$bmig <-  getBMIgroup(nhanes2$bmi)

# check:
library(doBy)
summaryBy(bmi ~ bmig, FUN = c(min, max), data = nhanes2)
summary(nhanes2)

### variables missingnes (%):
(miss100 <- sort(100 * colMeans(is.na(nhanes2)), decreasing = TRUE))
### complete cases (%):
(cc100 <- 100 * mean(complete.cases(nhanes2)))

mice::md.pattern(nhanes2, rotate.names = TRUE)

par(las = 1, cex.lab = 1.5, cex.axis = 1.2)
plot(nhanes2$bmi, nhanes2$chl,
     xlab = "BMI (kg/m²)",
     ylab = "Total serum cholesterol (mg/dl)",
     pch = 19, cex = 2, col = "blue")

getBetaBMI <- function(mydata = nhanes2) {
  mod <- lm(chl ~ bmi + hyp + age, data = mydata)
  coefs <- summary(mod)$coefficients
  ci <- confint(mod)
  coefs <- cbind(coefs, ci)
  beta <- coefs["bmi", c("Estimate", "Std. Error", "Pr(>|t|)", "2.5 %", "97.5 %")]
  return(beta)
}

effectCC <- getBetaBMI(nhanes2)
effectCC

imputeMean <- function(x) {
  xnew <- x
  xmean <- mean(x, na.rm = TRUE)
  xnew[is.na(x)] <- xmean
  return(xnew)
}

imputeMode <- function(x, seed = NULL) {
  xnew <- x
  freqs <- table(xnew)
  modepos <- which(freqs == max(freqs))
  # sampling if multimodal:
  l <- length(modepos)
  if (l > 1) {
    if (!is.null(seed))
      set.seed(seed) 
    modepos <- sample(1:l, size = 1)
  }
  xmode <- levels(x)[modepos]
  xnew[is.na(x)] <- xmode
  return(xnew)
 }

imputeMeanMode <- function(x, seed = 666) {
  if (!anyNA(x))
    return(x)
  condNum <- inherits(x, c("numeric", "integer"))
  condFac <- inherits(x, "factor")
  condOther <- !condNum & !condFac
  if (condNum) xnew <- imputeMean(x = x)
  if (condFac) xnew <- imputeMode(x = x, seed = seed)
  if (condOther) xnew <- x
  return(xnew)
 }

nhanes2I <- nhanes2
for (i in 1:dim(nhanes2I)[2])
 nhanes2I[, i] <- imputeMeanMode(nhanes2I[, i])
# for consistence, re-build bmig:
nhanes2I$bmig <- getBMIgroup(nhanes2I$bmi)
nhanes2I

effectI <- getBetaBMI(nhanes2I)
effectI

library(mice)
ini <- mice(nhanes2, maxit = 0)
ini ### note that age doesn't need to be imputed

### Get the imputation method and modify:
meth <- ini$method
### passive imputation for bmi:
meth["bmig"] <- "~getBMIgroup(bmi)"
### get the current predictorMatrix and modify in order "bmig" is not a predictor
### of "bmi":
pred <- ini$predictorMatrix
pred["bmi", "bmig"] <- 0



# Number of imputations:
M <- 6

# Multiple imputation (set printFlag = TRUE to see the progress):
imp <- mice(nhanes2,
            m = M,
            method = meth,
            predictorMatrix = pred,
            printFlag = FALSE,
            seed = 666)

imp$imp

complete(data = imp)

complete(data = imp, action = 2)

### 1) Fit the model to each of the imputed datasets:

modimp <- with(data = imp, exp = lm(chl ~ bmi + hyp + age))
modimp


### 2) Pool of fitted models:

modpool <- pool(modimp)
modpool


### 3) Get pooled estimates and confidence intervals:

coeffsMI <- summary(modpool, conf.int = TRUE, conf.level = 0.95)
coeffsMI

### 4) Get the association of interest:
effectMI <- coeffsMI[coeffsMI$term == "bmi",
                     c("estimate", "std.error", "p.value", "2.5 %", "97.5 %")]
names(effectMI) <- names(effectCC)
effectMI

results <- rbind(effectCC, effectI, effectMI)
rownames(results) <- c("complete cases",
                       "simple mean imputation",
                       "multiple imputation")
results

library(miclust)

data(minhanes)
class(minhanes)
length(minhanes)

### raw data in the first element of "minhanes":
nhanes <- minhanes[[1]]
dim(nhanes)
head(nhanes)
anyNA(nhanes)
summary(nhanes)

### first imputed data:
nhanes_imp1 <- minhanes[[2]]
### first imputed data set:
dim(nhanes_imp1)
head(nhanes_imp1)
anyNA(nhanes_imp1)
summary(nhanes_imp1)









### prepare data for miclust:
library(miclust)
minhanes1 <- miclust::getdata(data = minhanes)
class(minhanes1)
length(minhanes1)
names(minhanes1)


### exploration of the raw data:
head(minhanes1$rawdata)
summary(minhanes1$rawdata)

### exploration of the first imputed data set after transformations
### required to apply the clustering algorithm have been applied by
### the function ""getdata":
head(minhanes1$impdata[[1]])
summary(minhanes1$impdata[[1]])

nhanesmiclust <- miclust(data = minhanes1,
                         search = "backward",
                         ks = 2:3,
                         distance = "manhattan",
                         centpos = "means",
                         initcl = "hc",
                         verbose = TRUE,
                         seed = 666)

### print results:
nhanesmiclust

### plot results:
plot(nhanesmiclust,
     col.nclfreq = "blue",
     col.critcf = "orange",
     col.nvarfreq = "forestgreen")

nhanesmiclust_summary <- summary(nhanesmiclust, quantilevars = 0.5)
nhanesmiclust_summary

### selected variables
y <- getvariablesfrequency(nhanesmiclust)
y

# plot skeleton:
par(las = 1)
plot(y$percfreq, type = "n", main = "",
     xlab = "Variable", ylab = "Percentage of times selected",
     xlim = 0.5 + c(0, length(y$varnames)),
     lwd = 5, col = "blue",
     xaxt = "n", yaxt = "n")  # for customized axes
# customized axes:
axis(1, at = 1:length(y$varnames), labels = y$varnames, las = 2)
aty <- 10 * (0:10)
axis(2, at = aty)
# background grid:
abline(h = aty, lty = 2, col = "gray")
# data bars:
lines(y$percfreq, type = "h", lwd = 5, col = "blue")

myBinomApprox <- function(x, n, conf.level = 0.95) {
  z <- qnorm((1 + conf.level) / 2)
  p <- x / n
  error <- z * sqrt(p * (1 - p) / n)
  lower <- p - error
  upper <- p + error
  res <- data.frame(x, n, proportion = p, lower, upper, conf.level)
  return(res)
}

































# estimateRisks <- function(n, piE, pi0, pi1, propE, propD, seed = NULL) {
#   ...
# }

estimateRisks(n = 20000, piE = 1/5, pi0 = 1/25, pi1 = 1/10, propE = 1/2,
              propD = 1/2, seed = 4321)

# library(Epi)       # ci.lin

########################################################################
########################################################################
########################################################################
###
### simulateData: Function to simulate a database with a binary response
###               (D), two independent normal covariates (X1, X2), and a
###               binaryexposure (E).
###
########################################################################
###
### Input:
###
### n:           sample size
### muX:         vector with the means of X1, X2
### sigmaX:      vector with the sd of X1, X2
### pE:          prevalence of E
### pDmuXE0:     P(D|X=E(X), E=0)
### pDmuXE1:     P(D|X=E(X), E=1)
### ORDeltaX1E0: OR(DeltaX1=1|E=0)
### ORDeltaX1E1: OR(DeltaX1=1|E=1)
### ORDeltaX2E0: OR(DeltaX2=1|E=0)
### ORDeltaX2E1: OR(DeltaX2=1|E=1)
### seed:        seed for reproducibility of results.
###
########################################################################
###
### Output:
###
### a simulated database
###
########################################################################
########################################################################
########################################################################

simulateData <- function(n,
                         muX,
                         sigmaX,
                         pE,
                         pDmuXE0,
                         pDmuXE1,
                         ORDeltaX1E0,
                         ORDeltaX1E1,
                         ORDeltaX2E0,
                         ORDeltaX2E1,
                         seed = NULL) {
  if (!is.null(seed))
   set.seed(seed)
  X1 <- rnorm(n, mean = muX[1], sd = sigmaX[1]) 
  X2 <- rnorm(n, mean = muX[2], sd = sigmaX[2]) 
  E <- rbinom(n = n, size = 1, prob = pE)
  b1 <- log(ORDeltaX1E0)
  b2 <- log(ORDeltaX2E0)
  b0 <- log(pDmuXE0 / (1 - pDmuXE0)) - b1 * muX[1] - b2 * muX[2]
  b4 <- log(ORDeltaX1E1) - b1
  b5 <- log(ORDeltaX2E1) - b2
  numerb3 <- (1 - pDmuXE0) * pDmuXE1
  denomb3 <- pDmuXE0 * (1 - pDmuXE1)
  b3 <- log(numerb3 / denomb3) - b4 * muX[1] - b5 * muX[2]
  L <- b0 + b1 * X1 + b2 * X2 + b3 * E + (b4 * X1 + b5 * X2) * E
  pD = 1 / (1 + exp(-L))
  D <- apply(X = as.matrix(pD),
             MARGIN = 1,
             FUN = function(x) rbinom(n = 1, size = 1, prob = x))
  sim <- data.frame(X1, X2, E, D)
  return(sim)
}

########################################################################
########################################################################
###
### Simulated datasets:
###
########################################################################
########################################################################


########################################################################
### Simulating a database in which D depends only on E. This is equivalent
### to (check it):
# b1 = 0 --> ORDeltaX1E0 = 1
# b2 = 0 --> ORDeltaX2E0 = 1
# b4 = 0 --> ORDeltaX1E1 = 1
# b5 = 0 --> ORDeltaX2E1 = 1
########################################################################

dat0 <- simulateData(n = 5000,
                     muX = c(10, 20),
                     sigmaX = c(1, 3),
                     pE = 0.10,
                     pDmuXE0 = 0.05,
                     pDmuXE1 = 0.25,
                     ORDeltaX1E0 = 1,
                     ORDeltaX1E1 = 1,
                     ORDeltaX2E0 = 1,
                     ORDeltaX2E1 = 1,
                     seed = 666)

########################################################################
### Simulating a database in which D depends only on both E and X1, with
### no interactions. This is equivalent to (check it):
# b2 = 0 --> ORDeltaX2E0 = 1
# b4 = 0 --> ORDeltaX1E1 = 1
# b5 = 0 --> ORDeltaX2E1 = 1
########################################################################

dat1 <- simulateData(n = 5000,
                     muX = c(10, 20),
                     sigmaX = c(1, 3),
                     pE = 0.10,
                     pDmuXE0 = 0.05,
                     pDmuXE1 = 0.25,
                     ORDeltaX1E0 = 1.15,
                     ORDeltaX1E1 = 1,
                     ORDeltaX2E0 = 1,
                     ORDeltaX2E1 = 1,
                     seed = 666)

########################################################################
### Simulating a database in which D depends only on both E and X1, and
### the interaction (E * X1). This is equivalent to (check it):
# b2 = 0 --> ORDeltaX2E0 = 1
# b5 = 0 --> ORDeltaX2E1 = 1
########################################################################

dat2 <- simulateData(n = 5000,
                     muX = c(10, 20),
                     sigmaX = c(1, 3),
                     pE = 0.10,
                     pDmuXE0 = 0.05,
                     pDmuXE1 = 0.25,
                     ORDeltaX1E0 = 1.15,
                     ORDeltaX1E1 = 1.50,
                     ORDeltaX2E0 = 1,
                     ORDeltaX2E1 = 1)

########################################################################
### Simulating a database in which D depends only on E, X1 and X2, with
### no interactions. This is equivalent to (check it):
# b4 = 0 --> ORDeltaX1E1 = 1
# b5 = 0 --> ORDeltaX2E1 = 1
########################################################################

dat3 <- simulateData(n = 5000,
                     muX = c(10, 20),
                     sigmaX = c(1, 3),
                     pE = 0.10,
                     pDmuXE0 = 0.05,
                     pDmuXE1 = 0.25,
                     ORDeltaX1E0 = 1.15,
                     ORDeltaX1E1 = 1,
                     ORDeltaX2E0 = 1.20,
                     ORDeltaX2E1 = 1)

mod0 <- glm(D ~ E, data = dat0, family = binomial)
mod1 <- glm(D ~ E + X1, data = dat1, family = binomial)
mod2 <- glm(D ~ E * X1, data = dat2, family = binomial)
mod3 <- glm(D ~ E + X1 + X2, data = dat3, family = binomial)

summary(mod0)
summary(mod1)
summary(mod2)
summary(mod3)








### simulate a X vector of length 1000:
set.seed(666)
# vector length:
n <- 10000
X <- runif(n = n, min = 0.1, max = 50000)

### use microbenchmark() to compare computation times:
library(microbenchmark)

# number of runs:
m <- 5000
res <- microbenchmark(f1 = f1(x = X),
                      f2 = f2(x = X),
                      f3 = f3(x = X),
                      times = m)
res <- as.data.frame(res)

### convert time from nanoseconds to miliseconds:
res$time <- res$time / 1e6

### summary: mean, median, and 2.5th and 97.5th percentiles:
res_tab <- aggregate(time ~ expr,
                     data = res,
                     FUN = function(x) {
                       c(mean   = mean(x),
                         median = median(x),
                         p2.5   = quantile(x, 0.025),
                         p97.5  = quantile(x, 0.975))
                       })
# convert the aggregate() structure to data.frame:
res_tab <- do.call(what = "data.frame", args = res_tab)
names(res_tab) <- c("function", "mean", "median", "2.5%", "97.5%")

res_tab

compareEstimates <- function(n, p) {
  Xn <- rbinom(n = 1, size = n, prob = p)
  aproxIC <- as.numeric(prop.test(x = Xn, n = n, conf.level = 0.95)$conf.int)
  exactIC <- as.numeric(binom.test(x = Xn, n = n, conf.level = 0.95)$conf.int)
  res <- c(n = n, p = p, aproxIC, exactIC)
  names(res)[c(3, 5)] <- paste0("lower", c("Approx", "Exact"))
  names(res)[c(4, 6)] <- paste0("upper", c("Approx", "Exact"))
  return(res)
 }

nValues <- seq(10, 100, by = 10)
pValues <- c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99)
res <- expand.grid(n = nValues, p = pValues)
aux <- data.frame(matrix(nrow = nrow(res), ncol = 4))
auxnames <-  c("lowerApprox", "upperApprox", "lowerExact", "upperExact")
names(aux) <- auxnames
for (i in 1:nrow(aux))
 aux[i, ] <- compareEstimates(n = res$n[i], p = res$p[i])[auxnames]
myResults <- cbind(res, aux)

# par(mfrow = c(3, 3), las = 1, cex.axis = 0.8)
# xmin <- min(nValues)
# xmax <- max(nValues)
# ymin <- min(myResults[, c("lowerApprox", "lowerExact")])
# ymax <- max(myResults[, c("upperApprox", "upperExact")])
# delta <- 0.5
# legendPos <- c(rep("topright", 5), rep("bottomright", 4))
# for (i in 1:length(pValues)) {
#   aux <- myResults[myResults$p == pValues[i], -2]
#   ymin <- min(aux[, c("lowerApprox", "lowerExact")])
#   ymax <- max(aux[, c("upperApprox", "upperExact")])
#   plot(c(xmin, xmax), c(pValues[i], pValues[i]), type = "l",
#        main = paste("p =", pValues[i]), ylim = c(ymin, ymax),
#        xlab = "n", ylab = expression(hat(p)), lty = 2, col = "blue")
#   segments(aux$n - delta, aux$lowerApprox, aux$n - delta, aux$upperApprox,
#            col = "red", lwd = 2)
#   segments(aux$n + delta, aux$lowerExact, aux$n + delta, aux$upperExact,
#            col = "forestgreen", lwd = 2)
#   abline(h = 0, lty = 3)
#   abline(h = 1, lty = 3)
#   legend(legendPos[i], c("Approx.", "Exact"), col = c("red", "forestgreen"),
#          lwd = 2, bty = "n", cex = 0.8)
#  }

# ########################################################################
# ########################################################################
# ###
# ### 'multinomialSample': auxiliar function to simulate a random sample
# ###                      from a multivariate distribution.
# ###
# ########################################################################
# ###
# ### Input:
# ###
# ### n:     sample size
# ### probs: vector of length 'l' number proportional to the probabilities
# ### seed:  seed for results reproducibility
# ###
# ########################################################################
# ###
# ### Output:
# ###
# ### vector of length 'l' with the frequencies
# ###
# ###
# multinomialSample <- function(n, probs, seed = NULL) {
#   if (!is.null(seed))
#     set.seed(seed)
#   l <- length(probs)
#   res <- rep(0, l)
#   names(res) <- 1:l
#   samp <- rep(1:length(probs), rmultinom(1, size = n, probs))
#   aux <- table(samp)
#   res[names(aux)] <- aux
#   return(res)
# }
# ########################################################################
# ###
# ### END OF  'multinomialSample' FUNCTION
# ###
# ########################################################################

# ########################################################################
# ########################################################################
# ###
# ### 'estimateRisks': function to simulate a 2x2 contingency table under three
# ### different sampling designs and estimate bias of its probabilities
# ###
# #########################################################################
# ###
# ### Input:
# ###
# ### n:     sample size
# ### piE:   prevalence of the exposure E in the population
# ### pi0:   prevalence of the disease D among the non exposed in the population
# ### pi1:   prevalence of the disease D among the exposed in the population
# ### propE: proportion of exposed individuals in the sampling
# ### propD: proportion of cases in the sampling
# ### seed: seed for results reproducibility
# ###
# #########################################################################
# ###
# ### Output:
# ###
# ### vector of length 'l' with the frequencies
# ###
# ###
# estimateRisks <- function(n, piE, pi0, pi1, propE, propD, seed = NULL) {
#   # Notation:
#   #    D:  Case
#   #    ND: Control
#   #    E: Exposed
#   #    NE: Non exposed
#   #    a: contingency table [1, 1]
#   #    b: contingency table [1, 2]
#   #    c: contingency table [2, 1]
#   #    d: contingency table [2, 2]
# 
#   # Fix the seed if any:
#   if (!is.null(seed))
#     set.seed(seed)
# 
#   # List for the contingency tables:
#   tabs <- vector("list", 3)
# 
#   # Sampling design 1:
#   # rows: 1 (NE); 2 (E)
#   # columns: 1 (ND); 2 (D)
#   # Compute the join probabilities:
#   piTab1 <- (1 - piE) * (1 - pi0)
#   piTab2 <- ### TO COMPLETE
#   piTab3 <- ### TO COMPLETE
#   piTab4 <- ### TO COMPLETE
#   piTab <- c(piTab1, piTab2, piTab3, piTab4)
#   # Simulate the sample:
#   tab <- multinomialSample(n = n, probs = piTab)
#   tab <- matrix(tab, nrow = 2, ncol = 2, byrow = TRUE)
#   # Add margins:
#   tab <- addmargins(tab)
#   rownames(tab) <- c("non E", "E", "Total")
#   colnames(tab) <- c("non D", "D", "Total")
#   # Save table:
#   tabs[[1]] <- tab
# 
#   # Sampling design 2:
#   # rows: 1 (NE); 2 (E)
#   # columns: 1 (ND); 2 (D)
#   # Simulate row 1:
#   nE <- round(n * propE)
#   row1col2 <- ### TO COMPLETE
#   row1col1 <- ### TO COMPLETE
#   # Simulate row 2:
#   row2col2 <- ### TO COMPLETE
#   row2col1 <- ### TO COMPLETE
#   tab <- matrix(c(row1col1, row1col2, row2col1, row2col2),
#                 nrow = 2, ncol = 2, byrow = TRUE)
#   # Add margins:
#   tab <- addmargins(tab)
#   rownames(tab) <- c("non E", "E", "Total")
#   colnames(tab) <- c("non D", "D", "Total")
#   # Save table:
#   tabs[[2]] <- tab
# 
#   # Sampling design 3:
#   # rows: 1 (ND); 2 (D)
#   # columns: 1 (NE); 2 (E)
#   # Simulate row 1:
#   nD <- round(n * propD)
#   piD <- pi0 + (pi1 - pi0) * piE
#   row1col2 <- ### TO COMPLETE
#   row1col1 <- ### TO COMPLETE
#   # Simulate row 2:
#   row2col2 <- ### TO COMPLETE
#   row2col1 <- ### TO COMPLETE
#   tab <- matrix(c(row1col1, row1col2, row2col1, row2col2),
#                 nrow = 2, ncol = 2, byrow = TRUE)
#   # Add marginals:
#   tab <- addmargins(tab)
#   rownames(tab) <- c("non D", "D", "Total")
#   colnames(tab) <- c("non E", "E", "Total")
#   # Save table:
#   tabs[[3]] <- tab
# 
#   # Label tables:
#   names(tabs) <- paste0("Design", 1:3)
# 
#   # Summary table for the estimates:
#   est <- matrix(nrow = 4, ncol = 6)
#   rownames(est) <- c(names(tabs), "Population")
#   colnames(est) <- c("P(E)", "P(D)", "P(D|nonE)", "P(D|E)","P(E|nonD)","P(E|D)")
#   for (i in 1:2) {
#     est[i, 1] <- tabs[[i]][2, 3] / n 	
#     est[i, 2] <- tabs[[i]][3, 2] / n 	
#     est[i, 3] <- tabs[[i]][1, 2] / tabs[[i]][1, 3] 	
#     est[i, 4] <- tabs[[i]][2, 2] / tabs[[i]][2, 3] 	
#     est[i, 5] <- tabs[[i]][2, 1] / tabs[[i]][3, 1] 	
#     est[i, 6] <- tabs[[i]][2, 2] / tabs[[i]][3, 2]
#    }
#   i <- 3
#   est[i, 1] <- tabs[[i]][3, 2] / n 	
#   est[i, 2] <- tabs[[i]][2, 3] / n 	
#   est[i, 3] <- tabs[[i]][2, 1] / tabs[[i]][3, 1] 	
#   est[i, 4] <- tabs[[i]][2, 2] / tabs[[i]][3, 2] 	
#   est[i, 5] <- tabs[[i]][1, 2] / tabs[[i]][1, 3] 	
#   est[i, 6] <- tabs[[i]][2, 2] / tabs[[i]][2, 3]
#   est[4, ] <- c(piE, piD, pi0, pi1, (1 - pi1) * piE / (1 - piD), pi1 * piE / piD)
# 
#   # Bias (in percentage):
#   pop <- matrix(rep(est[4, ], 3), nrow = length(tabs), ncol = ncol(est),
#                 byrow = TRUE)
#   bias100 <- round(100 * (est[1:length(tabs), ] / pop - 1), 1)
# 
#   # Result:
#   res <- list(tables2x2 = tabs, estimates = est, percentBias = bias100)
#   return(res)
#  }
# #########################################################################
# ###
# ### END OF  'estimateRisks' FUNCTION
# ###
# #########################################################################
# #########################################################################
