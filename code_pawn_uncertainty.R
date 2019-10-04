

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
               "overlapping", "pawnr", "sensobol", "sensitivity", "wesanderson"))

# Set checkpoint

dir.create(".checkpoint")
library("checkpoint")

checkpoint("2019-09-22", 
           R.version ="3.6.1", 
           checkpointLocation = getwd())

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

# SAMPLE MATRICES #############################################################
###############################################################################

# CONSTRUCT SAMPLE MATRICES ---------------------------------------------------

A <- list()
for(i in k) {
  A[[i]] <-mclapply(N, function(N) sobol_matrices(n = N, k = i), mc.cores = n_cores)
}

A <- A[!sapply(A, is.null)]
names(A) <- models

for(i in names(A)) {
  names(A[[i]]) <- N
}

# MODEL OUTPUT ################################################################
###############################################################################

# COMPUTE MODEL OUTPUT --------------------------------------------------------

Y <- list()
for(i in names(A)) {
  if(i == "Liu") {
    Y[[i]] <- lapply(A[[i]], function(x) liu_Mapply(x))
  } else if(i == "Ishigami") {
    Y[[i]] <- lapply(A[[i]], function(x) sensobol::ishigami_Mapply(x))
  } else if(i == "Sobol' G") {
    Y[[i]] <- lapply(A[[i]], function(x) sensobol::sobol_Fun(x))
  } else {
    Y[[i]] <- lapply(A[[i]], function(x) sensitivity::morris.fun(x))
  }
}

names(Y) <- models
for(i in names(Y)) {
  names(Y[[i]]) <- N
}

# PLOT MODEL UNCERTAINTY ------------------------------------------------------

lapply(models, function(models) Y[[models]]$`600`) %>%
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

# SOBOL' INDICES ##############################################################
###############################################################################

# COMPUTE SOBOL' INDICES AND THEIR CONFIDENCE INTERVALS -----------------------

out <- out.ci <- list()
for(i in names(A)) {
  for(j in names(A[[i]])) {
    out[[i]][[j]] <- sobol_indices(Y[[i]][[j]], 
                                   params = params[[i]], 
                                   n = as.numeric(j), 
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

# SOBOL INDICES AND CONFIDENCE INTERVALS OF DUMMY PARAMETER -------------------

sobol.dummy <- sobol.dummy.ci <- list()
for(i in names(A)) {
  for(j in names(A[[i]])) {
    sobol.dummy[[i]][[j]] <- sobol_dummy(Y[[i]][[j]], 
                                         params = params[[i]], 
                                         R = R, 
                                         n = as.numeric(j), 
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

# SOBOL' CONVERGENCE ----------------------------------------------------------

sobol.convergence <- lapply(out.ci, function(x) rbindlist(x, idcol = "N")) %>%
  rbindlist(., idcol = "model") %>%
  .[, N:= as.numeric(N)] %>%
  .[, diff:= high.ci - low.ci] %>%
  .[, model:= factor(model, levels = c("Liu", "Ishigami", 
                                       "Sobol' G", "Morris"))] %>%
  .[, parameters:= factor(parameters, 
                          levels = paste("X", 1:20, sep = ""))] %>%
  .[, method:= "Sobol'"] %>%
  .[, .(model, N, parameters, original, low.ci, high.ci, diff, method, sensitivity)]

# PAWN INDICES ################################################################
###############################################################################

# COMPUTE PAWN INDICES AND THEIR CONFIDENCE INTERVALS -------------------------

# Subset to take only the A matrix and the model output of the A matrix
Y.pawn <- A.pawn <- list()
for(i in names(Y)) {
  for(j in names(Y[[i]])) {
    Y.pawn[[i]][[j]] <- Y[[i]][[j]][1:j]
    A.pawn[[i]][[j]] <- A[[i]][[j]][1:j, ]
  }
}

# Compute PAWN indices and their confidence intervals
pawn.indices <- pawn.ci <- list()
for(i in names(A.pawn)) {
  for(j in names(A.pawn[[i]]) ) {
    pawn.indices[[i]][[j]] <- pawn_generic(data = A.pawn[[i]][[j]], 
                                           Y = Y.pawn[[i]][[j]], 
                                           n = n, 
                                           test = median, 
                                           R = R)
    pawn.ci[[i]][[j]] <- pawn_ci(pawn.indices[[i]][[j]])
  }
}

# PAWN AND CONFIDENCE INTERVALS OF DUMMY PARAMETER ----------------------------

pawn.index.dummy <- list()
for(i in names(Y)) {
  for(j in names(Y[[i]]) ) {
    pawn.index.dummy[[i]][[j]] <- pawn_dummy(Y = Y[[i]][[j]], 
                                             n = n, 
                                             R = R)
  }
}

pawn.index.dummy <- lapply(pawn.index.dummy, function(x) rbindlist(x, idcol = "N")) %>%
  rbindlist(., idcol = "model") %>%
  .[, model:= factor(model, levels = c(c("Liu", "Ishigami", 
                                         "Sobol' G", "Morris")))]

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

# PLOT CONVERGENCE ############################################################
###############################################################################

# PLOT CONVERGENCE ------------------------------------------------------------

sobol.convergence[sensitivity == "STi"] %>%
  .[, sensitivity:= NULL] %>%
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

# PLOT SOBOL' AND PAWN INDICES ------------------------------------------------

# Sobol' indices
a <- plot_sobol(sobol.convergence[N==max(N)], 
                dummy = sobol.dummy.final[N==max(N)]) +
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
b <- pawn.convergence[N==max(N)] %>%
  plot_pawn(.) +
  geom_rect(data = pawn.index.dummy[N==max(N)],
            aes(ymin = 0,
                ymax = high.ci,
                xmin = -Inf,
                xmax = Inf),
            fill = "black",
            alpha = 0.1,
            inherit.aes = FALSE) +
  facet_grid(~model, 
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

# PLOT SOBOL' AND PAWN INDICES (INDIVIDUAL PLOTS) -----------------------------

a
b
