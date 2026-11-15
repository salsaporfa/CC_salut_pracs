## ----libraries, include=FALSE------------------------------------------------------
library(knitr)
library(highlight)  # to highlight R output
library(xtable)     # to create nice LaTeX tables
library(epitools)
library(mice)
library(miclust)




## ----Description_of_missingness, include=TRUE, echo = TRUE-------------------------
load("copd_redux.RData")
n <- length(copd$id) #participant number
v <- ncol(copd) #variables number
M <- n*v #cell number

missing <- colSums(is.na(copd)) 
miss_percentatge <- (missing/n)*100 #missing 

minimum_missing_pp <- min(miss_percentatge) #minimum missing percentage
maximum_missing_pp <- max(miss_percentatge) #maximum missing percentage

miss_under_5 <- sum(miss_percentatge < 5)
pctg_under_5 <- (miss_under_5/v) * 100
#percentage of variables with less than 5% missing

complete <- sum(miss_percentatge == 0)
complete_pctg <- (complete/v)*100
#percentage of patients with complete data on all variables

overall <- sum(is.na(copd))
overall_pctg <- (overall/M)*100 
#Overall percentage of the information is missing


## ----Multiple_imputation_1, include = TRUE, echo = TRUE----------------------------
str(copd)


## ----Multiple_imputation_2, include=TRUE, echo=TRUE--------------------------------
categoricas <- names(copd)[sapply(copd, is.factor)]
categoricas_latex <- gsub("_", "\\\\_", categoricas)
niveles <- sapply(copd[categoricas], nlevels)
binarias <- names(niveles[niveles == 2])
binarias_latex <- gsub("_", "\\\\_", binarias)
no_binarias <- names(niveles[niveles > 2])
no_binarias_latex <- gsub("_", "\\\\_", no_binarias)
level_no_binaria <- niveles[niveles > 2]


## ----Multiple_imputation_3, include=TRUE, echo=TRUE, results='asis'----------------
missing_info <- data.frame(
  variable = names(copd),
  missing_n = colSums(is.na(copd)),
  missing_pct = round(colMeans(is.na(copd)) * 100, 2)
)

completed_variables <- names(copd)[colSums(is.na(copd)) == 0]
completed_variables_latex <- gsub("_", "\\\\_", completed_variables)

xt <- xtable(
  missing_info,
  caption = "Missing description",
  label = "tab:missing_description",
  align = c("l", "l", "r", "r"),
  digits = c(0, 0, 0, 2)
)

print(
  xt,
  include.rownames = FALSE,
  include.colnames = TRUE,
  caption.placement = "top",
  table.placement = "H"
)


## ----Multiple_imputation_4, include=TRUE, echo=TRUE--------------------------------
ini <- mice(copd, maxit = 0)
meth <- ini$method
p_Matrix <- ini$predictorMatrix
meth["id"] <- ""
meth[no_binarias] <- "polr"
for (var in completed_variables) {
  p_Matrix[var, ] <- 0 #the completed variables are not imputed
}
p_Matrix[,"id"] <- 0 #the variable id is not considered as a predictor


## ----Multiple_imputation_5, include=TRUE, echo=TRUE, message=FALSE-----------------
imp <- mice(copd, m = 100, 
            method = meth, predictorMatrix = p_Matrix, 
            seed = 666,  
            printFlag = FALSE)


## ----mclust_1, include=TRUE, echo=TRUE---------------------------------------------
#help(getdata)
bd <- complete(imp, action="all", include=TRUE)
levels(bd[[1]]$graudisn_v1) # revising the variable graudisn_v1 having 6 levels
# bd preparation
bd_preparada <- lapply(bd, function(df) {
  # We do not consider the id
  df$id <- NULL
  # Identifacte the factor variables
  factors <- names(df)[sapply(df, is.factor)]
  # Consider these numerical levels as factors
  df[factors] <- lapply(df[factors], function(x) {
    as.numeric(as.character(x))
  })
  df
})
dades <- getdata(bd_preparada)


## ----mclust_2, include=TRUE, echo=TRUE---------------------------------------------
n <- nrow(dades$rawdata)
p <- ncol(dades$rawdata)
n / p
cluster_results <- miclust(dades, search = "backward",seed = 666)


## ----result1, echo=TRUE, results='hide', fig.width=8, fig.height=6, out.width='\\linewidth'----
plot(cluster_results,
     k = NULL,
     metric = "critcf",
     col.nclfreq = "gray",
     col.critcf = "gray",
     col.nvarfreq = "gray",
     col.varsel = "black")


## ----result2, echo=TRUE, results='hide', fig.width=8, fig.height=6, out.width='\\linewidth'----
frequency <- getvariablesfrequency(cluster_results)

barplot(
  frequency$percfreq,
  names.arg = frequency$varnames
)

