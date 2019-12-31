## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----preliminary steps, results="hide", message=FALSE, warning=FALSE-----

# PRELIMINARY FUNCTIONS -------------------------------------------------------

# Install the development version of the pawnr package
devtools::install_github("arnaldpuy/pawnr", build_vignettes = TRUE)

# Function to read in all required packages in one go:
loadPackages <- function(x) {
  for(i in x) {
    if(!require(i, character.only = TRUE)) {
      install.packages(i, dependencies = TRUE)
      library(i, character.only = TRUE)
    }
  }
}

# Load the packages
loadPackages(c("tidyverse", "data.table", "randtoolbox", "sensitivity", 
               "boot", "parallel", "doParallel", "scales", "cowplot", 
               "overlapping", "pawnr", "sensobol", "sensitivity", "wesanderson", 
               "ggridge"))

# Set checkpoint

dir.create(".checkpoint")
library("checkpoint")

checkpoint("2019-09-22", 
           R.version ="3.6.1", 
           checkpointLocation = getwd())


## ----settings, cache = TRUE----------------------------------------------

# DEFINE SETTINGS -------------------------------------------------------------

N <- seq(500, 10000, 250) # Sample sizes
n <- 15 # Number of conditioning intervals
k <- c(2, 3, 8, 20) # Vector with number of model inputs
R <- 100 # Bootstrap replicas
n_cores <- floor(detectCores() * 0.75) # Use 75% of the cores available
type <- "norm" # Define the confidence interval method
conf <- 0.95 # Define the ci
models <- c("Liu", "Ishigami", "Sobol' G", "Morris")
params <- lapply(k, function(x) paste("X", 1:x, sep = ""))
names(params) <- models

# Function to compute the Liu et al. function
liu <- function(X1, X2) {
  X1 / X2
}

liu_Mapply <- function(X) {
  X[, 1] <- qchisq(X[, 1], df = 10)
  X[, 2] <- qchisq(X[, 2], df = 13.978)
  return(mapply(liu, X[, 1], X[, 2]))
}


## ----sample_matrices, cache=TRUE, dependson="settings"-------------------

# CONSTRUCT SAMPLE MATRICES ---------------------------------------------------

A <- B <- list()
for(i in k) {
  # For Sobol' STi
  A[[i]] <-mclapply(N, function(N) sobol_matrices(n = floor(N / (i + 1)), k = i), mc.cores = n_cores)
  # For PAWN
  B[[i]] <- mclapply(N, function(N) randtoolbox::sobol(n = N, dim = i))
}

A <- A[!sapply(A, is.null)]
B <- B[!sapply(B, is.null)]

names(A) <- models
names(B) <- models

for(i in names(A)) {
  names(A[[i]]) <- N
}

for(i in names(B)) {
  names(B[[i]]) <- N
}


## ----output, cache=TRUE, dependson="sample_matrices"---------------------

# COMPUTE MODEL OUTPUT --------------------------------------------------------

Y <- Y.pawn <- list()
for(i in names(A)) {
  if(i == "Liu") {
    Y[[i]] <- lapply(A[[i]], function(x) liu_Mapply(x))
    Y.pawn[[i]] <- lapply(B[[i]], function(x) liu_Mapply(x))
  } else if(i == "Ishigami") {
    Y[[i]] <- lapply(A[[i]], function(x) sensobol::ishigami_Mapply(x))
    Y.pawn[[i]] <- lapply(B[[i]], function(x) sensobol::ishigami_Mapply(x))
  } else if(i == "Sobol' G") {
    Y[[i]] <- lapply(A[[i]], function(x) sensobol::sobol_Fun(x))
    Y.pawn[[i]] <- lapply(B[[i]], function(x) sensobol::sobol_Fun(x))
  } else {
    Y[[i]] <- lapply(A[[i]], function(x) sensitivity::morris.fun(x))
    Y.pawn[[i]] <- lapply(B[[i]], function(x) sensitivity::morris.fun(x))
  }
}

names(Y) <- models
for(i in names(Y)) {
  names(Y[[i]]) <- N
}

names(Y.pawn) <- models
for(i in names(Y.pawn)) {
  names(Y.pawn[[i]]) <- N
}


## ----model_uncertainty, cache=TRUE, dependson="output", dev="tikz"-------

# PLOT MODEL UNCERTAINTY ------------------------------------------------------

lapply(models, function(models) Y.pawn[[models]]$`10000`) %>%
  do.call(cbind, .) %>%
  data.table() %>%
  setnames(., 1:4, models) %>%
  melt(., measure.vars = 1:4) %>%
  .[, variable:= factor(variable, levels = models)] %>%
  ggplot(., aes(value)) +
  geom_histogram(color = "black",
                 fill = "white") +
  labs(x = "Y",
       y = "Count") +
  facet_wrap(~ variable,
             scales = "free_x",
             ncol = 4) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent",
                                         color = NA),
        legend.key = element_rect(fill = "transparent",
                                  color = NA))


## ----sobol_indices, cache=TRUE, dependson="output"-----------------------

# COMPUTE SOBOL' INDICES AND THEIR CONFIDENCE INTERVALS -----------------------

out <- out.ci <- list()
for(i in names(A)) {
  for(j in names(A[[i]])) {
    out[[i]][[j]] <- sobol_indices(Y[[i]][[j]],
                                   params = params[[i]],
                                   n = floor(as.numeric(j) / (length(params[[i]]) + 1)),
                                   type = "saltelli",
                                   R = R,
                                   parallel = "multicore",
                                   ncpus = n_cores)
    out.ci[[i]][[j]] <- sobol_ci(out[[i]][[j]],
                                 params = params[[i]],
                                 type = type,
                                 conf = conf)
  }
}


## ----sobol_indices_dummy, cache=TRUE, dependson="sobol_indices"----------

# SOBOL INDICES AND CONFIDENCE INTERVALS OF DUMMY PARAMETER -------------------

sobol.dummy <- sobol.dummy.ci <- list()
for(i in names(A)) {
  for(j in names(A[[i]])) {
    sobol.dummy[[i]][[j]] <- sobol_dummy(Y[[i]][[j]],
                                         params = params[[i]],
                                         R = R,
                                         n = floor(as.numeric(j) / (length(params[[i]]) + 1)),
                                         parallel = "multicore",
                                         ncpus = n_cores)
    sobol.dummy.ci[[i]][[j]] <- sobol_ci_dummy(sobol.dummy[[i]][[j]],
                                               type = type,
                                               conf = conf)
  }
}

sobol.dummy.final <- lapply(sobol.dummy.ci, function(x) rbindlist(x, idcol = "N")) %>%
  rbindlist(., idcol = "model") %>%
  .[, model:= factor(model, levels = c("Liu", "Ishigami",
                                       "Sobol' G", "Morris"))]


## ----sobol_convergence_dt, cache=TRUE, dependson="sobol_indices"---------

# SOBOL' CONVERGENCE ---------------------------------------------------------

sobol.convergence <- lapply(out.ci, function(x) rbindlist(x, idcol = "N")) %>%
  rbindlist(., idcol = "model") %>%
  .[, N:= as.numeric(N)] %>%
  .[, diff:= high.ci - low.ci] %>%
  .[, model:= factor(model, levels = c("Liu", "Ishigami",
                                       "Sobol' G", "Morris"))] %>%
  .[, parameters:= factor(parameters,
                          levels = paste("X", 1:20, sep = ""))] %>%
  .[, method:= "$S_{Ti}^*$"] %>%
  .[, .(model, N, parameters, original, low.ci, high.ci, diff, method, sensitivity)]


## ----pawn_indices, cache=TRUE, dependson=c("sample_matrices", "output")----

# COMPUTE PAWN INDICES AND THEIR CONFIDENCE INTERVALS -------------------------

pawn.indices <- pawn.ci <- list()
for(i in names(B)) {
  for(j in names(B[[i]]) ) {
    pawn.indices[[i]][[j]] <- pawn_generic(data = B[[i]][[j]],
                                           Y = Y.pawn[[i]][[j]],
                                           n = n,
                                           test = median,
                                           R = R)
    pawn.ci[[i]][[j]] <- pawn_ci(pawn.indices[[i]][[j]])
  }
}


## ----pawn_dummy, cache=TRUE, dependson=c("sample_matrices", "output", "pawn_indices")----

# PAWN AND CONFIDENCE INTERVALS OF DUMMY PARAMETER ----------------------------

pawn.index.dummy <- list()
for(i in names(Y.pawn)) {
  for(j in names(Y.pawn[[i]]) ) {
    pawn.index.dummy[[i]][[j]] <- pawn_dummy(Y = Y.pawn[[i]][[j]],
                                             n = n,
                                             R = R)
  }
}
pawn.index.dummy <- lapply(pawn.index.dummy, function(x) rbindlist(x, idcol = "N")) %>%
  rbindlist(., idcol = "model") %>%
  .[, model:= factor(model, levels = c(c("Liu", "Ishigami",
                                         "Sobol' G", "Morris")))]


## ----pawn_convergence_dt, cache=TRUE, dependson="pawn_indices"-----------

# PAWN CONVERGENCE ------------------------------------------------------------

pawn.convergence <- lapply(pawn.ci, function(x) rbindlist(x, idcol = "N")) %>%
  rbindlist(., idcol = "model") %>%
  .[, N:= as.numeric(N)] %>%
  .[, diff:= high.ci - low.ci] %>%
  .[, model:= factor(model, levels = c("Liu", "Ishigami",
                                       "Sobol' G", "Morris"))] %>%
  .[, parameters:= gsub("V", "X", parameters)] %>%
  .[, parameters:= factor(parameters,
                          levels = paste("X", 1:20, sep = ""))] %>%
  .[, method:= "PAWN"]


## ----export_sobol_pawn_convergence, cache=TRUE, dependson=c("pawn_convergence_dt", "sobol_convergence_dt")----

# EXPORT SOBOL' AND PAWN CONVERGENCE RATES ------------------------------------

fwrite(sobol.convergence, "sobol.convergence.csv")
fwrite(pawn.convergence, "pawn.convergence.csv")


## ----plot_convergence, cache=TRUE, dependson=c("sobol_convergence_dt", "pawn_convergence_dt"), fig.height=4, fig.width=6.3, dev="tikz"----

# PLOT CONVERGENCE ------------------------------------------------------------

sobol.convergence[sensitivity == "STi"] %>%
  .[, sensitivity:= NULL] %>%
  .[, method:= factor(method, levels = c("PAWN", "$S_{Ti}^*$"))] %>%
  rbind(., pawn.convergence) %>%
  ggplot(., aes(N, diff,
                group = parameters,
                color = parameters)) +
  geom_line() +
  geom_hline(yintercept = 0.05,
             lty = 2) +
  scale_color_discrete(name = "Model inputs") +
  labs(y = expression(Stat[indices]),
       x = "N") +
  facet_grid(method~model) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent",
                                         color = NA),
        legend.key = element_rect(fill = "transparent",
                                  color = NA))


## ----plot_convergence2_samples, cache=TRUE, dependson=c("sobol_convergence_dt", "pawn_convergence_dt"), fig.height=4, fig.width=6.3, dev="tikz"----

# PLOT CONVERGENCE (SHOWING THE RANGE OF SAMPLES USED) ------------------------

sobol.convergence[sensitivity == "STi"] %>%
  .[, sensitivity:= NULL] %>%
  rbind(., pawn.convergence) %>%
  .[, method:= factor(method, levels = c("PAWN", "$S_{Ti}^*$"))] %>%
  ggplot(., aes(N, diff, 
                group = parameters, 
                color = parameters)) +
  geom_line() +
  annotate("rect", 
           xmin = 200, 
           xmax = 2000, 
           ymin = 0, 
           ymax = Inf, 
           alpha = 0.1, 
           fill="red") +
  annotate("rect", xmin = 2500, 
           xmax = 4000, 
           ymin = 0, 
           ymax = Inf, 
           alpha = 0.1, 
           fill="green") +
  geom_hline(yintercept = 0.05, 
             lty = 2) +
  scale_color_discrete(name = "Model inputs") +
  labs(y = expression(Stat[indices]), 
       x = "N") +
  facet_grid(method~model) +
  theme_bw() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent",
                                         color = NA),
        legend.key = element_rect(fill = "transparent",
                                  color = NA))


## ----plot_full_convergence, cache=TRUE, dependson=c("sobol_convergence.dt", "pawn_convergence_dt"), dev="tikz"----

# PLOT SOBOL' AND PAWN INDICES ------------------------------------------------

# Sobol' indices
a <- plot_sobol(sobol.convergence[N==10000],
                dummy = sobol.dummy.final[N==10000]) +
  facet_grid(~model,
             scales = "free_x",
             space = "free_x") +
  labs(x = "",
       y = "Sobol' index") +
  theme(axis.text.x = element_text(size = 6),
        legend.position = "none")
# Get legend
legend <- get_legend(a + theme(legend.position = "top"))

# PAWN indices
b <- pawn.convergence[N==10000] %>%
  plot_pawn(.) +
  geom_rect(data = pawn.index.dummy[N==10000],
            aes(ymin = 0,
                ymax = high.ci,
                xmin = -Inf,
                xmax = Inf),
            fill = "black",
            alpha = 0.1,
            inherit.aes = FALSE) +
  labs(x = "",
       y = "PAWN") +
  facet_grid(~ model,
             scales = "free_x",
             space = "free_x") +
  theme(axis.text.x = element_text(size = 6))
# Merge
bottom <- plot_grid(a, b,
                    ncol = 1,
                    labels = "auto",
                    align = "h")
plot_grid(legend, bottom,
          labels = c("", ""),
          ncol = 1,
          align = "",
          rel_heights = c(0.1, 1))


## ----plot_convergence2, cache=TRUE, dependson="plot_full_convergence", dev="tikz", fig.height=2, fig.width=6.5----

# PLOT SOBOL' AND PAWN INDICES (INDIVIDUAL PLOTS) -----------------------------

a
b


## ----pawn_model, cache=TRUE----------------------------------------------

# THE MODEL -------------------------------------------------------------------

# Function to divide a vector into chunks
chunks <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

# The model
model_pawn <- function(Model, N, n, epsilon, theta) {
  # Check which model to apply to set the number of 
  # parameters
  if(Model == 1) {
    k <- 2
  } else if(Model == 2) {
    k <- 3
  } else if(Model == 3) {
    k <- 8
  } else {
    k <- 20
  }
  # Create the Sobol' matrix
  data <- randtoolbox::sobol(n = N, dim = k) 
  # Transform distribution:
  if(Model == 1) {
    ModelRun <- liu_Mapply
  } else if(Model == 2) {
    ModelRun <- sensobol::ishigami_Mapply
  } else if(Model == 3) {
    ModelRun <- sensobol::sobol_Fun
  } else {
    ModelRun <- sensitivity::morris.fun
  }
  # Run the model
  Y <- ModelRun(data)
  # Set seed to fix the random number generator
  # for the sample function below
  set.seed(epsilon)
  # Sample the unconditional model output
  index <- sample(1:nrow(data),
                  size = floor(nrow(data) / n),
                  # Without replacement
                  replace = FALSE)
  # Bind model inputs and model output
  dt <- data.table::data.table(cbind(data, Y))
  # Subset and obtain the unconditional model output
  Y_unc <- dt[index, Y]
  # Create the intervals
  melted <- data.table::melt(dt,
                             measure.vars = 1:(ncol(dt) - 1),
                             variable.name = "parameters")
  # Compute PAWN indices
  out <- melted[order(parameters, value)][
    , .(chunks(Y, n)), parameters][
    , Y_unc:= .(rep(.(Y_unc), times = n * ncol(data)))][
    , ID:= .I][
    , results:= .(.(mapply(stats::ks.test, Y_unc, V1))), ID][
    , statistic:= sapply(results, function(x) x[, 1]$statistic)]
  if(theta == 1) {
    final <-  out[, mean(statistic), parameters][, V1]
  } else if(theta == 2) {
    final <-  out[, median(statistic), parameters][, V1]
  } else {
    final <-  out[, max(statistic), parameters][, V1]
  }
  return(final)
}


## ----pawn_settings, cache=TRUE-------------------------------------------

# DEFINE SETTINGS --------------------------------------------------------------

# Set sample size
n <- 2 ^ 13

# Define N.min and N.max
N.min <- 200
N.max <- 2000

# Set parameters
parameters <- c("N", "n", "epsilon", "theta")

# Vector with name of functions
models <- c("Liu", "Ishigami", "Sobol' G", "Morris")


## ----pawn_matrix, cache=TRUE, dependson="pawn_settings"------------------

# CREATION OF THE MATRICES ----------------------------------------------------

# Create the A, B and AB matrices, also for the 
# computation of second and third-order indices
tmp <- lapply(1:4, function(x)
  sobol_matrices(n = n, 
                 k = length(parameters), 
                 second = TRUE, 
                 third = TRUE) %>%
    data.table()) 

# Name the slots
names(tmp) <- 1:4

# Rename columns
tmp <- lapply(tmp, setnames, parameters) %>%
  rbindlist(., idcol = "Model")

# Create two copies of the sample matrix and list the 
# original and the copies. One would be to run the 
# calculations in the max in theta setting; the
# other one for the max not in theta setting, 
# and the other in the optimum setting
max <- copy(tmp)
A <- list(tmp, max, copy(tmp))

# Name the slots
names(A) <- c("max", "no.max", "optimum")

# Transform all distributions 
for(i in names(A)) {
  if(i == "max") {
    # where 1=mean, 2=median, 3=max in the model
    A[[i]][, N:= floor(qunif(N, N.min, N.max))]
    A[[i]][, n:= floor(qunif(n, 5, 20))]
    A[[i]][, theta:= floor(theta * (3 - 1 + 1)) + 1]
  } else if(i == "no.max") {
    A[[i]][, N:= floor(qunif(N, N.min, N.max))]
    A[[i]][, n:= floor(qunif(n, 5, 20))]
    A[[i]][, theta:= floor(theta * (2 - 1 + 1)) + 1]
  } else {
    A[[i]][, N:= floor(qunif(N, N.max, 4000))]
    A[[i]][, n:= floor(qunif(n, 15, 20))]
    A[[i]][, theta:= floor(theta * (2 - 1 + 1)) + 1]
  }
}

# Transform all the other distributions
A.pawn <- rbindlist(A, idcol = "setting")[
  , epsilon:= floor(qunif(epsilon, 1, 1000))][
  , Model:= as.numeric(Model)]

print(A.pawn)


## ----run_pawn_model, cache=TRUE, dependson=c("pawn_matrices", "pawn_model", "pawn_settings")----

# RUN MODEL -------------------------------------------------------------------

# Define parallel computing
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Compute
Y.pawn <- foreach(i=1:nrow(A.pawn), 
             .packages = "data.table") %dopar%
  {
    model_pawn(epsilon = A.pawn[[i, "epsilon"]], 
               N = A.pawn[[i, "N"]], 
               n = A.pawn[[i, "n"]], 
               theta = A.pawn[[i, "theta"]],
               Model = A.pawn[[i, "Model"]])
  }

# Stop parallel cluster
stopCluster(cl)


## ----extract_pawn_data, cache=TRUE, dependson="run_pawn_model"-----------

# EXTRACT DATA ----------------------------------------------------------------

rowNumber <- lapply(1:4, function(x) A.pawn[, .I[Model == x]])
names(rowNumber) <- models

out <- list()
for(i in models) {
  out[[i]] <- Y.pawn[rowNumber[[i]]]
}

dt.models <- list()
for(i in seq_along(1:4)) {
  dt.models[[i]] <- cbind(A.pawn[Model == i], data.table(do.call(rbind, out[[i]])))
}


## ----uncertainty_dt, cache=TRUE, dependson="extract_pawn_data"-----------

# DATASET FOR UNCERTAINTY ANALYSIS --------------------------------------------

AB.pawn <- lapply(dt.models, function(x)  {
  x[, .SD[1: (2 * (2 ^ 13))], setting] %>%
    melt(., measure.vars = patterns("V"), 
         variable.name = "parameter")
  }) %>%
  rbindlist() %>%
  .[, Model:= ifelse(Model == 1, models[1], 
                     ifelse(Model == 2, models[2], 
                            ifelse(Model == 3, models[3], models[4])))] %>%
  .[, parameter:= gsub("V", "X", parameter)] %>%
  .[, parameter:= factor(parameter, 
                         levels = paste("X", 1:20, sep = ""))] %>%
  .[, Model:= factor(Model, 
                     levels = c("Liu", "Ishigami", "Sobol' G", "Morris"))] %>%
  .[, setting:= ifelse(setting == "max", "$max \\in \\theta$", 
                      ifelse(setting == "no.max", "$max \\notin \\theta$", "Optimum"))]


## ----pawn_overlap, cache=TRUE, dependson="uncertainty_dt"----------------

# CHECK OVERLAP ---------------------------------------------------------------

overlap.dt <- split(AB.pawn, AB.pawn$setting)

overlap.results <- mclapply(overlap.dt, function(x) {
  split(x, x$Model, drop = TRUE) %>%
    lapply(., function(x) split(x, x$parameter, drop = TRUE)) %>%
    lapply(., function(x) lapply(x, function(y) y[, value])) %>%
    lapply(., function(x) overlap(x))}, 
  mc.cores = n_cores)

tmp <- lapply(overlap.results, function(x) lapply(x, function(y) { 
  cbind(y$OV) %>%
    data.frame() %>%
    setDT(., keep.rownames = TRUE)
  })) 

pawn.overlap.results <- lapply(tmp, function(x) 
  rbindlist(x, idcol = "Model")) %>%
  rbindlist(., idcol = "setting") %>%
  setnames(., ".", "overlap")

par.overlap <- paste("X", 1:6, sep = "")

final.overlap <- lapply(models, function(x) pawn.overlap.results[Model==x, .SD, setting]) %>%
  lapply(., function(x) x[, "overlap":= round(.SD, 3), .SDcols = "overlap"])

final.overlap

# Export results
rbindlist(final.overlap) %>%
  fwrite(., "pawn.overlap.csv")


## ----plot_pawn_uncertainty, cache=TRUE, dependson="uncertainty_dt", dev="tikz"----

# PLOT UNCERTAINTY ------------------------------------------------------------

plot.uncertainty.pawn <- ggplot(AB.pawn, aes(value,
                                             fill = parameter,
                                             color = parameter)) +
  geom_density(alpha = 0.5,
               position = "identity") +
  facet_grid(setting~Model, 
             scales = "free_y") +
  scale_fill_discrete(name = "Model input") +
  scale_color_discrete(guide = FALSE) +
  labs(x = "PAWN",
       y = "Density") +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  theme_bw() +
  theme(legend.position = "top", 
        legend.box = "horizontal",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent", 
                                         color = NA), 
        legend.key = element_rect(fill = "transparent",  
                                  color = NA))

plot.uncertainty.pawn


## ----export_pawn, cache=TRUE, dependson="uncertainty_dt"-----------------

# EXPORT AB MATRIX FOR PAWN ---------------------------------------------------

fwrite(AB.pawn, "AB.pawn.csv")


## ----pawn_sensitivity_dt, cache=TRUE, dependson="extract_pawn_data"------

# DATASET FOR SENSITIVITY ANALYSIS --------------------------------------------

dt.pawn.sens <- lapply(dt.models, function(x) 
  melt(x, measure.vars = patterns("V"), variable.name = "model.input")) %>%
  rbindlist() %>%
  .[, Model:= ifelse(Model == 1, models[1], 
                     ifelse(Model == 2, models[2], 
                            ifelse(Model == 3, models[3], models[4])))] %>%
  .[, model.input:= gsub("V", "X", model.input)] %>%
  .[, model.input:= factor(model.input, 
                         levels = paste("X", 1:20, sep = ""))] %>%
  .[, Model:= factor(Model, 
                     levels = c("Liu", "Ishigami", "Sobol' G", "Morris"))] %>%
  setnames(., "value", "Y") 

# EXPORT AB MATRIX FOR SENSITIVITY --------------------------------------------

fwrite(dt.pawn.sens, "dt.pawn.sens.csv")


## ----pawn_sensitivity, cache=TRUE, dependson="pawn_sensitivity_dt"-------

# SENSITIVITY ANALYSIS --------------------------------------------------------

pawn.sensitivity <- dt.pawn.sens[, sobol_indices(Y, 
                                                 params = parameters,
                                                 R = R, 
                                                 n = 2 ^ 13, 
                                                 parallel = "multicore", 
                                                 second = TRUE, 
                                                 third = TRUE,
                                                 ncpus = n_cores), 
                                 .(setting, Model, model.input)]


## ----pawn_ci, cache=TRUE, dependson="pawn_sensitivity"-------------------

# CONFIDENCE INTERVALS --------------------------------------------------------

# Arrange data
tmp3 <- split(pawn.sensitivity, pawn.sensitivity$setting) %>%
  lapply(., function(x) split(x, x$Model)) %>%
  lapply(., function(x) lapply(x, function(y) split(y, y$model.input, drop = TRUE)))

# Compute confidence intervals
pawn.ci <- list()
for(i in names(tmp3)) {
  for(j in names(tmp3[[i]])) {
    for(k in names(tmp3[[i]][[j]])) {
      pawn.ci[[i]][[j]][[k]] <- sobol_ci(tmp3[[i]][[j]][[k]], 
                                          params = parameters, 
                                          type = type, 
                                          conf = conf, 
                                          second = TRUE, 
                                          third = TRUE)
    }
  }
}

# Rearrange data
final.pawn.ci <- lapply(pawn.ci, function(x) 
  lapply(x, function(y) rbindlist(y, idcol = "model.input"))) %>%
  lapply(., function(x) rbindlist(x, idcol = "model")) %>%
  rbindlist(., idcol = "setting") %>%
  .[, model:= factor(model, levels = c("Liu", "Ishigami", 
                                       "Sobol' G", "Morris"))] %>%
  .[, model.input:= factor(model.input, levels = paste("X", 1:20, sep = ""))] %>%
  .[, parameters:= gsub("epsilon", "$\\\\varepsilon$", parameters)] %>%
  .[, parameters:= gsub("theta", "$\\\\theta$", parameters)] %>%
  .[, setting:= ifelse(setting == "max", "$max \\in \\theta$", 
                       ifelse(setting == "no.max", "$max \\notin \\theta$", "Optimum"))]

# EXPORT DATA -----------------------------------------------------------------

fwrite(final.pawn.ci, "final.pawn.ci.csv")


## ----plot_pawn_sensitivity, cache=TRUE, dependson="pawn_ci", dev="tikz", fig.height=2.7----

# PLOT AGGREGATED SOBOL' INDICES ----------------------------------------------

a <- final.pawn.ci[sensitivity == "Si" | sensitivity == "STi"] %>%
  ggplot(., aes(parameters, original,
                fill = sensitivity)) +
  geom_boxplot(outlier.size = 0.2) +
  labs(x = "", 
       y = "Sobol' index") +
  scale_fill_discrete(name = "Sobol' indices",
                      labels = c(expression(S[italic(i)]),
                                 expression(S[italic(T[i])]))) +
  theme_bw() +
  facet_wrap(~ setting) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent", 
                                         color = NA), 
        legend.key = element_rect(fill = "transparent",  
                                  color = NA),
        legend.position = "none")

legend <- get_legend(a + theme(legend.position = "top"))

# PLOT SUM OF SOBOL' SI ------------------------------------------------------

b <- final.pawn.ci[sensitivity == "Si"][
  , sum(original), .(setting, model, model.input)] %>%
  ggplot(., aes(setting, V1)) +
  geom_boxplot(outlier.size = 0.2) +
  labs(x = "", 
       y = expression(paste("Sum of"~S[i]))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# MERGE AGGREGATED SOBOL' AND SUM OF SOBOL' -----------------------------------

up <- plot_grid(legend, NULL, 
                ncol = 2)

bottom <- plot_grid(a, b, 
                    ncol = 2, 
                    align = "hv", 
                    labels = "auto",
                    rel_widths = c(2.2, 1))

plot_grid(up, bottom, 
          ncol = 1, 
          align = "hv", 
          rel_heights = c(0.21, 1))


## ----sobol_model, cache=TRUE---------------------------------------------

# THE MODEL -------------------------------------------------------------------

# Functions to create A and AB matrices to compute Ti
scrambled_sobol <- function(A, B) {
  X <- rbind(A, B)
  for(i in 1:ncol(A)) {
    AB <- A
    AB[, i] <- B[, i]
    X <- rbind(X, AB)
  }
  AB <- rbind(A, X[((2*nrow(A)) + 1):nrow(X), ])
  return(AB)
}

sobol_matrix <- function(n, k) {
  df <- randtoolbox::sobol(n = n, dim = k * 2)
  A <- df[, 1:k]
  B <- df[, (k + 1) : (k * 2)]
  out <- scrambled_sobol(A = A, B = B)
  return(out)
}

# Functions to estimate Ti
sobol_all <- function(Y_A, Y_AB, type) {
  n <- length(Y_A[!is.na(Y_A)])
  f0 <- (1 / n) * sum(Y_A)
  VY <- 1 / n * sum((Y_A - f0) ^ 2)
  if(type == "jansen") {
    STi <- ((1 / (2 * n)) * sum((Y_A - Y_AB) ^ 2)) / VY
  } else if(type == "homma") {
    STi <- (VY - (1 / n) * sum(Y_A * Y_AB) + f0 ^ 2) / VY
  } else if(type == "sobol") {
    STi <- ((1 / n) * sum(Y_A * (Y_A - Y_AB))) / VY
  } else {
    stop("type should  be either jansen, sobol or homma")
  }
  return(STi)
}

sobol_Ti_Mapply <- function(d, type) {
  return(mapply(sobol_all,
                MoreArgs = list(type = type),
                d[, "Y_A"],
                d[, "Y_AB"]))
}

sobol_Ti <- function(Y, params, type) {
  k <- length(params)
  p <- length(1:(length(Y) / (k + 1)))
  Y_A <- Y[1:p]
  Y_AB <- Y[(p + 1):length(Y)]
  parameters <- rep(params, each = length(Y_A))
  vec <- cbind(Y_A, Y_AB)
  out <- data.table(vec, parameters)
  output <- out[, sobol_Ti_Mapply(.SD, type = type), parameters][, V1]
  return(output)
}

# The model
model_sobol <- function(Model, N, k, Theta) {
  data <- sobol_matrix(n = floor(N / (k + 1)), k = k) 
  if(Model == 1) {
    Y <- liu_Mapply(data)
  } else if(Model == 2) {
    Y <- sensobol::ishigami_Mapply(data)
  } else if(Model == 3) {
    Y <- sensobol::sobol_Fun(data)
  } else {
    Y <- sensitivity::morris.fun(data)
  }
  out <- sobol_Ti(Y, params = paste("X", 1:k, sep = ""), type = Theta)
  return(out)
}


## ----sobol_settings, cache=TRUE------------------------------------------

# DEFINE SETTINGS -------------------------------------------------------------

# Set parameters
parameters.sobol <- c("N", "Theta")


## ----sobol_matrix, cache=TRUE, dependson=c("pawn_settings", "sobol_settings")----

# CREATION OF THE MATRICES ----------------------------------------------------

# Create the A and AB matrices, also for the 
# computation of second and third-order indices
tmp <- lapply(models, function(x)
  sobol_matrices(n = n, k = length(parameters.sobol)) %>%
    data.table()) 

# Rename columns and transform distributions
A <- lapply(tmp, setnames, parameters.sobol) %>%
  rbindlist(., idcol = "Model")

# Create two copies of the sample matrix and list the 
# original and the copies. One would be to run the 
# calculations with uncertainty in N and Theta, 
# the other with uncertainty in N only.
N.only <- copy(A)
A.DT <- list(A, N.only)
names(A.DT) <- c("N.Theta", "N")

A <- rbindlist(A.DT, idcol = "setting")

A.sobol <- A[, k:= ifelse(Model == 1, 2, ifelse(Model == 2, 3, ifelse(Model == 3, 8, 20)))][
    , N:= floor(qunif(N, N.min, N.max))][
    , Model:= as.numeric(Model)][
    , Theta:= floor(Theta * (3 - 1 + 1)) + 1][
    , Theta:= ifelse(Theta == 1, "jansen", ifelse(Theta == 2, "homma", "sobol"))][
    , Theta:= ifelse(setting == "N", "jansen", Theta)]

print(A.sobol)
print(n)




## ----run_sobol_model, cache=TRUE, dependson=c("pawn_matrices", "sobol_matrix", "sobol_settings", "sobol_model")----

# RUN MODEL -------------------------------------------------------------------

# Define parallel computing
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Compute
Y.sobol <- foreach(i=1:nrow(A.sobol), 
                   .packages = "data.table") %dopar%
  {
    model_sobol(N = A.sobol[[i, "N"]], 
                Theta = A.sobol[[i, "Theta"]],
                Model = A.sobol[[i, "Model"]], 
                k = A.sobol[[i, "k"]])
  }

# Stop parallel cluster
stopCluster(cl)


## ----extract_sobol_data, cache=TRUE, dependson="run_sobol_model"---------

# EXTRACT MODEL OUTPUT --------------------------------------------------------

rowNumber <- lapply(1:4, function(x) A.sobol[, .I[Model == x]])
names(rowNumber) <- models

out <- list()
for(i in models) {
  out[[i]] <- Y.sobol[rowNumber[[i]]]
}

dt.models <- list()
for(i in seq_along(1:4)) {
  dt.models[[i]] <- cbind(A[Model == i], data.table(do.call(rbind, out[[i]])))
}


## ----sobol_uncertainty_dt, cache=TRUE, dependson="extract_sobol_data"----

# DATASET FOR UNCERTAINTY ANALYSIS --------------------------------------------

AB.sobol <- lapply(dt.models, function(x) {
  x[, .SD[1: (2 * (2 ^ 13))], setting] %>% 
  melt(., measure.vars = patterns("V"), 
       variable.name = "parameter")}) %>%
  rbindlist(.) %>%
  .[, Model:= ifelse(Model == 1, models[1], 
                     ifelse(Model == 2, models[2], 
                            ifelse(Model == 3, models[3], models[4])))] %>%
  .[, k:= NULL] %>%
  .[, parameter:= gsub("V", "X", parameter)] %>%
  .[, parameter:= factor(parameter, 
                         levels = paste("X", 1:20, sep = ""))] %>%
  .[, Model:= factor(Model, 
                     levels = c("Liu", "Ishigami", "Sobol' G", "Morris"))] %>%
  .[, setting:= ifelse(setting == "N.Theta", "$N,\\theta$", setting)]

# EXPORT AB MATRIX FOR SOBOL' -------------------------------------------------

fwrite(AB.sobol, "AB.sobol.csv")


## ----plot_sobol_uncertainties, cache=TRUE, dependson="sobol_uncertainty_dt", dev="tikz"----

# PLOT UNCERTAINTY ------------------------------------------------------------

AB.sobol %>%
  ggplot(., aes(value,
                fill = parameter,
                color = parameter)) +
  geom_density(alpha = 0.5,
               position = "identity") +
  facet_grid(setting~Model) +
  scale_fill_discrete(name = "Model input") +
  scale_color_discrete(guide = FALSE) +
  labs(x = "Y",
       y = "Count") +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  scale_y_continuous(limits = c(0, 35)) +
  theme_bw() +
  theme(legend.position = "top", 
        legend.box = "horizontal",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent", 
                                         color = NA), 
        legend.key = element_rect(fill = "transparent",  
                                  color = NA))


## ----sobol_overlap, cache=TRUE, dependson="sobol_uncertainty_dt"---------

# CHECK OVERLAP ---------------------------------------------------------------

overlap.dt <- split(AB.sobol, AB.sobol$setting)

overlap.results <- mclapply(overlap.dt, function(x) {
  split(x, x$Model, drop = TRUE) %>%
    lapply(., function(x) split(x, x$parameter, drop = TRUE)) %>%
    lapply(., function(x) lapply(x, function(y) y[, value])) %>%
    lapply(., function(x) overlap(x))}, 
  mc.cores = n_cores)

tmp <- lapply(overlap.results, function(x) lapply(x, function(y) { 
  cbind(y$OV) %>%
    data.frame() %>%
    setDT(., keep.rownames = TRUE)
})) 

sobol.overlap.results <- lapply(tmp, function(x) 
  rbindlist(x, idcol = "Model")) %>%
  rbindlist(., idcol = "setting") %>%
  setnames(., ".", "overlap")

par.overlap <- paste("X", 1:6, sep = "")

final.overlap <- lapply(models, function(x) sobol.overlap.results[Model==x, .SD, setting]) %>%
  lapply(., function(x) x[, "overlap":= round(.SD, 3), .SDcols = "overlap"])

final.overlap

# Export results
rbindlist(final.overlap) %>%
  fwrite(., "sobol.overlap.csv")


## ----sobol_sensitivity_dt, cache=TRUE, dependson="extract_sobol_data"----

# DATASET FOR SENSITIVITY ANALYSIS --------------------------------------------

full.dataset.sobol <- lapply(dt.models, function(x) 
  melt(x, measure.vars = patterns("V"), 
       variable.name = "parameter")) %>%
  rbindlist(.) %>%
  .[, Model:= ifelse(Model == 1, models[1], 
                     ifelse(Model == 2, models[2], 
                            ifelse(Model == 3, models[3], models[4])))] %>%
  .[, k:= NULL] %>%
  .[, parameter:= gsub("V", "X", parameter)] %>%
  .[, parameter:= factor(parameter, 
                         levels = paste("X", 1:20, sep = ""))] %>%
  .[, Model:= factor(Model, 
                     levels = c("Liu", "Ishigami", "Sobol' G", "Morris"))] %>%
  .[, setting:= ifelse(setting == "N.Theta", "$N,\\theta$", setting)]

# EXPORT SENSITIVITY MATRIX ---------------------------------------------------

fwrite(full.dataset.sobol, "full.dataset.sobol.csv")


## ----sobol_sensitivity, cache=TRUE, dependson="sobol_sensitivity_dt"-----

# SENSITIVITY ANALYSIS --------------------------------------------------------

sobol.sensitivity <- full.dataset.sobol[, sobol_indices(value, 
                                                        type = "jansen", 
                                                        params = parameters.sobol, 
                                                        n = 2 ^ 13, 
                                                        R = R, 
                                                        parallel = "multicore", 
                                                        ncpus = n_cores), 
                                        .(Model, parameter, setting)]


## ----sobol_ci, cache=TRUE, dependson="sobol_sensitivity"-----------------

# CONFIDENCE INTERVALS --------------------------------------------------------

# Arrange data
tmp3 <- split(sobol.sensitivity, sobol.sensitivity$setting) %>%
  lapply(., function(x) split(x, x$Model)) %>%
  lapply(., function(x) lapply(x, function(y) split(y, y$parameter, drop = TRUE)))

# Compute confidence intervals
out <- list()
for(i in names(tmp3)) {
  for(j in names(tmp3[[i]])) {
    for(k in names(tmp3[[i]][[j]])) {
      out[[i]][[j]][[k]] <- sobol_ci(tmp3[[i]][[j]][[k]], 
                                         params = parameters.sobol, 
                                         type = type, 
                                         conf = conf)
    }
  }
}

# ARRANGE DATA ----------------------------------------------------------------

final.sobol <- lapply(out, function(x) 
  lapply(x, function(y) rbindlist(y, idcol = "model.input"))) %>%
  lapply(., function(x) rbindlist(x, idcol = "Model")) %>%
  rbindlist(., idcol = "setting") %>%
  .[, Model:= factor(Model, levels = c("Liu", "Ishigami", "Sobol' G", "Morris"))] %>%
  .[, model.input:= factor(model.input, levels = paste("X", 1:20, sep = ""))] %>%
  .[, parameters:= gsub("Theta", "$\\\\theta$", parameters)] 

# EXPORT DATA -----------------------------------------------------------------

fwrite(final.sobol, "final.sobol.csv")


## ----plot_sobol, cache=TRUE, dependson="sobol_ci", dev="tikz", fig.height=2, fig.width=2.5----

# PLOT SOBOL INDICES ----------------------------------------------------------

ggplot(final.sobol, aes(parameters, original,
                        fill = sensitivity)) +
  geom_boxplot(outlier.size = 0.2) +
  labs(x = "", 
       y = "Sobol' index") +
  scale_fill_discrete(name = "Sobol' indices",
                      labels = c(expression(S[italic(i)]),
                                 expression(S[italic(T[i])]))) +
  theme_bw() +
  facet_wrap(~setting) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent", 
                                         color = NA), 
        legend.key = element_rect(fill = "transparent",  
                                  color = NA),
        legend.position = "none")


## ----global_uncertainty, cache=TRUE, fig.height=7, dependson=c("uncertainty_dt", "sobol_uncertainty_dt"), dev="tikz"----

# MERGE UNCERTAINTY IN PAWN AND SOBOL'-----------------------------------------

a <- plot.uncertainty.pawn +
  theme(legend.position = "none")

b <- AB.sobol[!setting == "N"] %>%
  ggplot(., aes(value,
                fill = parameter,
                color = parameter)) +
  geom_density(alpha = 0.5,
               position = "identity") +
  facet_grid(setting~Model) +
  scale_fill_discrete(name = "Model input") +
  scale_color_discrete(guide = FALSE) +
  labs(x = "$S_{Ti}^*$",
       y = "Density") +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  scale_y_continuous(limits = c(0, 35)) +
  theme_bw() +
  theme(legend.position = "none", 
        legend.box = "horizontal",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent", 
                                         color = NA), 
        legend.key = element_rect(fill = "transparent",  
                                  color = NA))

# Get legend
legend <- get_legend(a + theme(legend.position = "top"))

# Merge
bottom <- plot_grid(a, b, 
                    ncol = 1, 
                    labels = "auto", 
                    align = "h", 
                    rel_heights = c(1, 0.46))

plot_grid(legend, bottom, 
          labels = c("", ""), 
          ncol = 1,
          align = "", 
          rel_heights = c(0.2, 1))


## ----global_sensitivities, dev="tikz", cache=TRUE, fig.height=2.5, dependson=c("pawn_ci", "sobol_ci")----

# PLOT AGGREGATED SOBOL' INDICES ----------------------------------------------

a <- final.pawn.ci[sensitivity == "Si" | sensitivity == "STi"] %>%
  ggplot(., aes(parameters, original,
                fill = sensitivity)) +
  geom_boxplot(outlier.size = 0.2) +
  labs(x = "", 
       y = "Sobol' index") +
  scale_fill_discrete(name = "Sobol' indices",
                      labels = c(expression(S[italic(i)]),
                                 expression(S[italic(T[i])]))) +
  theme_bw() +
  facet_wrap(~ setting) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent", 
                                         color = NA), 
        legend.key = element_rect(fill = "transparent",  
                                  color = NA),
        legend.position = "none")

legend <- get_legend(a + theme(legend.position = "top"))


b <- final.sobol[!setting == "N"] %>%
  ggplot(., aes(parameters, original, fill = sensitivity)) +
  geom_boxplot(outlier.size = 0.2) +
  labs(x = "", 
       y = "") +
  scale_fill_discrete(name = expression(paste("Sobol'"~T[i])),
                      labels = c(expression(S[italic(i)]),
                                 expression(S[italic(T[i])]))) +
  facet_wrap(~ setting) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent", 
                                         color = NA), 
        legend.key = element_rect(fill = "transparent",  
                                  color = NA),
        legend.position = "none")

bottom <- plot_grid(a, b, 
                    ncol = 2, 
                    align = "hv", 
                    labels = "auto",
                    rel_widths = c(2.58, 1))

plot_grid(legend, bottom, 
          ncol = 1, 
          align = "hv", 
          rel_heights = c(0.21, 1))


## ----plot_aggregated_Si, cache=TRUE, dev="tikz", fig.height=3, fig.width=3.1----

# PLOT AGGREGATED SUM OF SI ---------------------------------------------------

final.sobol2 <- setnames(final.sobol, "Model", "model")

rbind(final.pawn.ci[sensitivity == "Si"][, type:= "PAWN"], 
      final.sobol2[!setting == "N" & sensitivity == "Si"][, type:= "$S_{Ti}^*$"]) %>%
  .[, sum(original), .(setting, model, model.input, type)] %>%
  .[, setting:= factor(setting, levels = c("$max \\in \\theta$", 
                                           "$max \\notin \\theta$", 
                                           "Optimum", 
                                           "$N,\\theta$"))] %>%
  ggplot(., aes(setting, V1, fill = type)) +
  scale_fill_grey(start = 0.5, end = 0.9, 
                  name = "") +
  geom_boxplot(outlier.size = 0.2) +
  labs(x = "", 
       y = expression(paste("Sum of"~S[i]))) +
  theme_bw() +
  theme(legend.position = "top", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


## ----plot_single_functions, cache=TRUE, dependson="pawn_ci", dev="tikz", fig.height=8, fig.width=7----

# ARRANGE TO PLOT SOBOL' INDICES FOR EACH FUNCTION ---------------------------

tmp <- split(final.pawn.ci, final.pawn.ci$model)
gg <- list()
for(i in names(tmp)) {
  for(j in 1:3) {
    gg[[i]][[j]] <- plot_sobol(tmp[[i]], type = j) +
      scale_y_continuous(breaks = pretty_breaks(n = 3)) +
      facet_grid(model.input ~ setting) +
      labs(x = "", 
           y = "Sobol' index") +
      theme(legend.position = "none")
  }
}
# Extract legend
legend <- get_legend(gg[[1]][[1]] + theme(legend.position = "top"))

# PLOT SOBOL' INDICES FOR LIU, ISHIGAMI AND SOBOL' G --------------------------

all <- lapply(1:3, function(x) {
  left <- plot_grid(gg[[1]][[x]], gg[[2]][[x]],
                    labels = c("a", "b"),
                    align = "h",
                    ncol = 1)
  plot_grid(left, gg[[3]][[x]],
            labels = c("", "c"),
            ncol = 2)
})

plot_grid(legend, all[[1]], 
          ncol = 1, 
          rel_heights = c(0.1, 1))

lapply(2:3, function(x) all[x])


## ----plot_morris, cache=TRUE, dependson="pawn_ci", dev="tikz", fig.height=10, fig.width=5----

# PLOT SOBOL' INDICES FOR THE MORRIS FUNCTION ---------------------------------

plot_grid(legend, gg[[4]][[1]], 
          ncol = 1, 
          rel_heights = c(0.07, 1))

lapply(2:3, function(x) gg[[4]][[x]])


## ----merge_second_third, cahce=TRUE, dependson="pawn_ci", dev="tikz", fig.height=2.7----

# MERGE SECOND AND THIRD-ORDER EFFECTS -----------------------------------------

gg <- list()
second.third <- c("Sij", "Sijk")
for(i in second.third) {
  gg[[i]] <- final.pawn.ci[sensitivity == i] %>%
    ggplot(., aes(parameters, original)) +
    geom_boxplot(outlier.size = 0.2) +
    labs(x = NULL,
         y = "Sobol' index") +
    scale_fill_discrete(name = "Sobol' indices",
                        labels = c(expression(S[italic(i)]),
                                   expression(S[italic(T[i])]))) +
    theme_bw() +
    geom_hline(yintercept = 0,
               lty = 2,
               color = "red") +
    scale_y_continuous(breaks = pretty_breaks(n = 3)) +
    facet_wrap(~ setting) +
    theme(legend.position = "none", 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45,
                                     hjust = 1))
}

# PLOT SECOND AND THIRD-ORDER EFFECTS -----------------------------------------

plot_grid(gg[[1]],
          gg[[2]] + labs(x = "", y = ""),
          ncol = 2,
          labels = "auto",
          align = "hv")


## ----sum_Si_weighted, cache=TRUE, dependson="pawn_ci", dev="tikz", fig.height=2.5----

# PLOT AGGREGATED SOBOL' INDICES AFTER WEIGHTING ------------------------------

a <- final.pawn.ci[sensitivity == "Si" |sensitivity == "STi"] %>%
  # For each function, setting and design parameter, compute 
  # the median value of Si and STi
  .[, .(Median = median(original)), 
    by = .(setting, model, sensitivity, parameters)] %>%
  # Compute the aggregated median and the percentiles
  .[, .(aggregated.median = median(Median), 
        low.ci = quantile(Median, probs = 0.025),
        high.ci = quantile(Median, probs = 0.975)), 
    by = .(setting, sensitivity, parameters)] %>%
  ggplot(., aes(parameters, aggregated.median,
                color = sensitivity)) +
  geom_point(position = position_dodge(0.6)) +
  geom_errorbar(aes(ymin = low.ci, 
                    ymax = high.ci), 
                position = position_dodge(0.6)) +
  labs(x = "", 
       y = "Sobol' index") +
  scale_color_discrete(name = "Sobol' indices",
                       labels = c(expression(S[italic(i)]),
                                  expression(S[italic(T[i])]))) +
  theme_bw() +
  facet_wrap(~ setting) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent", 
                                         color = NA), 
        legend.key = element_rect(fill = "transparent",  
                                  color = NA))

b <- final.sobol[!setting == "N"] %>%
  # For each function, setting and design parameter, compute 
  # the median value of Si and STi
  .[, .(Median = median(original)), 
    by = .(setting, Model, sensitivity, parameters)] %>%
  # Compute the aggregated median and the percentiles
  .[, .(aggregated.median = median(Median), 
        low.ci = quantile(Median, probs = 0.025),
        high.ci = quantile(Median, probs = 0.975)), 
    by = .(setting, sensitivity, parameters)] %>%
  ggplot(., aes(parameters, aggregated.median,
                color = sensitivity)) +
  geom_point(position = position_dodge(0.6)) +
  geom_errorbar(aes(ymin = low.ci, 
                    ymax = high.ci), 
                position = position_dodge(0.6)) +
  labs(x = "", 
       y = NULL) +
  scale_color_discrete(name = "Sobol' indices",
                       labels = c(expression(S[italic(i)]),
                                  expression(S[italic(T[i])]))) +
  theme_bw() +
  facet_wrap(~ setting) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent", 
                                         color = NA), 
        legend.key = element_rect(fill = "transparent",  
                                  color = NA))

all <- plot_grid(a, b, 
                 ncol = 2, 
                 align = "hv",
                 rel_widths = c(2.58, 1),
                 labels = "auto")

plot_grid(legend, all, 
          ncol = 1, 
          align = "hv", 
          rel_heights = c(0.21, 1))


## ----plot_sobol_functions, cache=TRUE, dependson="sobol_ci", dev="tikz", fig.height=8, fig.width=4.5----

# PLOT DESIGN PARAMETERS FOR SOBOL: ALL FUNCTIONS -----------------------------

tmp <- final.sobol[!setting == "N"] %>%
  split(., .$Model)
gg <- list()
for(i in names(tmp)) {
    gg[[i]] <- plot_sobol(tmp[[i]], type = 1) +
      scale_y_continuous(breaks = pretty_breaks(n = 3)) +
      facet_grid(model.input ~.) +
      labs(x = "", 
           y = "Sobol' index") +
      theme(legend.position = "none")
}

# Extract legend
legend <- get_legend(gg[[1]] + theme(legend.position = "top"))

# PLOT SOBOL' INDICES FOR LIU, ISHIGAMI AND SOBOL' G --------------------------

left <- plot_grid(gg[[1]], gg[[2]],
                  labels = c("a", "b"),
                  align = "h",
                  ncol = 1)
all <- plot_grid(left, gg[[3]], 
                 labels = c("", "c"), 
                 ncol = 2)

plot_grid(legend, all, 
          ncol = 1, 
          rel_heights = c(0.1, 1))


## ----plot_morris_sobol, cache=TRUE, dependson="sobol_ci", dev="tikz", fig.height=10, fig.width=3----

# PLOT SOBOL' INDICES FOR THE MORRIS FUNCTION ---------------------------------

plot_grid(legend, gg[[4]], 
          ncol = 1, 
          rel_heights = c(0.07, 1))


## ----session_information-------------------------------------------------

# SESSION INFORMATION ---------------------------------------------------------

sessionInfo()


