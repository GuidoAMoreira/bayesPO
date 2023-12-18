## ----set_options, include = FALSE---------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(bayesPO)
library(ggplot2)
library(bayesplot)
library(MASS)
theme_set(theme_bw())
color_scheme_set("green")
set.seed(123456789)
temp <- tempfile(fileext = ".rda")
d <- download.file("https://drive.google.com/uc?id=1WoBmyVFj_PI3zGcxIaIvGbRZFUy8_NNU&export=download",
  temp, mode = "wb", quiet = TRUE)
# Try to use downloaded version. If not available, run the model
if (!d) load(temp, verbose = TRUE) else {
  warning("Data failed to download from drive. Please check internet connection and try again.")
}

## ----set_true_params----------------------------------------------------------
if (!d) {
  beta <- c(-1, 2) # Intercept = -1. Only one covariate
delta <- c(3, 4) # Intercept = 3. Only one covariate
lambdaStar <- 1000
} else warning("Data failed to download from drive. Please check internet connection and try again.")


## ----simulation---------------------------------------------------------------
if (!d) {
  # Spread a Poisson amount of points randomly in the random square.
total_points <- rpois(1, lambdaStar)
random_points <- cbind(runif(total_points), runif(total_points))
grid_size <- 50

# Find covariate values to explain the species occurrence.
# We give them a Gaussian spatial structure.
reg_grid <- as.matrix(expand.grid(seq(0, 1, len = grid_size), seq(0, 1, len = grid_size)))
#### Next is commented to save time. Uncomment to replicate results
# Z <- mvrnorm(1, rep(0, total_points + grid_size * grid_size), 3 * exp(-as.matrix(dist(rbind(random_points, reg_grid))) / 0.2))
Z1 <- Z[1:total_points]; Z2 <- Z[-(1:total_points)]

# Thin the points by comparing the retaining probabilities with uniforms
# in the log scale to find the occurrences
#### Next is commented to save time. Uncomment to replicate results
# occurrences <- log(runif(total_points)) <= -log1p(exp(-beta[1] - beta[2] * Z1))
n_occurrences <- sum(occurrences)
occurrences_points <- random_points[occurrences,]
occurrences_Z <- Z1[occurrences]

# Find covariate values to explain the observation bias.
# Additionally create a regular grid to plot the covariate later.
#### Next is commented to save time. Uncomment to replicate results
# W <- mvrnorm(1, rep(0, n_occurrences + grid_size * grid_size), 2 * exp(-as.matrix(dist(rbind(occurrences_points, reg_grid))) / 0.3))
W1 <- W[1:n_occurrences]; W2 <- W[-(1:n_occurrences)]

# Find the presence-only observations.
#### Next is commented to save time. Uncomment to replicate results
# po_sightings <- log(runif(n_occurrences)) <= -log1p(exp(-delta[1] - delta[2] * W1))
n_po <- sum(po_sightings)
po_points <- occurrences_points[po_sightings, ]
po_Z <- occurrences_Z[po_sightings]
po_W <- W1[po_sightings]
} else warning("Data failed to download from drive. Please check internet connection and try again.")


## ----plot_observer_bias, fig.width = 5, fig.height = 3.5, fig.align = 'center'----
if (!d) {
  ggplot(
  data.frame(x = reg_grid[, 1], y = reg_grid[, 2], Covariate = W2),
  aes(x, y)
) +
  geom_tile(aes(fill = Covariate)) +
  geom_point(aes(shape = Occurrences),
             data = data.frame(x = occurrences_points[, 1],
                               y = occurrences_points[, 2],
                               Occurrences = po_sightings), size = 1) +
  geom_point(data = data.frame(x = occurrences_points[!po_sightings, 1],
                               y = occurrences_points[!po_sightings, 2]),
             size = 1, shape = 1, color = "white") +
  scale_shape_manual(values = c(1, 19), labels = c("Not observed", "Observed"))
} else warning("Data failed to download from drive. Please check internet connection and try again.")

## ----prior--------------------------------------------------------------------
if (!d) {
  jointPrior <- prior(
  NormalPrior(rep(0, 2), 10 * diag(2)), # Beta
  NormalPrior(rep(0, 2), 10 * diag(2)), # Delta
  GammaPrior(0.00001, 0.00001) # LambdaStar
)
} else warning("Data failed to download from drive. Please check internet connection and try again.")

## ----create_model-------------------------------------------------------------
if (!d) {
  model <- bayesPO_model(po = cbind(po_Z, po_W),
                       intensitySelection = 1, observabilitySelection = 2,
                       intensityLink = "logit", observabilityLink = "logit",
                       initial_values = 2, joint_prior = jointPrior)
} else warning("Data failed to download from drive. Please check internet connection and try again.")

## ----fit----------------------------------------------------------------------
if (!d) {
  bkg <- cbind(Z2, W2) # Create background
} else warning("Data failed to download from drive. Please check internet connection and try again.")

#### Next is commented to save time. Uncomment to replicate results
# fit <- fit_bayesPO(model, bkg, area = 1, mcmc_setup = list(burnin = 1000, iter = 2000))

## ----print_fit----------------------------------------------------------------
if (!d) {
  summary(fit)
} else warning("Data failed to download from drive. Please check internet connection and try again.")

## ----fit2---------------------------------------------------------------------
#### Next is commented to save time. Uncomment to replicate results
# fit2 <- fit_bayesPO(fit, bkg, mcmc_setup = list(iter = 10000))

## ----print_fit2---------------------------------------------------------------
if (!d) {
  summary(fit2)
} else warning("Data failed to download from drive. Please check internet connection and try again.")

## ----bayesplot, fig.width = 5, fig.width = 5, fig.align = 'center'------------
if (!d) {
  mcmc_trace(fit2)
} else warning("Data failed to download from drive. Please check internet connection and try again.")

## ----marginals, fig.width = 5, fig.width = 5, fig.align = 'center'------------
if (!d) {
  mcmc_dens(fit2)
} else warning("Data failed to download from drive. Please check internet connection and try again.")

