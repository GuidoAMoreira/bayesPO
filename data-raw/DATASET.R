library(MASS)
theme_set(theme_bw())
color_scheme_set("green")
set.seed(123456789)

beta <- c(-1, 2) # Intercept = -1. Only one covariate
delta <- c(3, 4) # Intercept = 3. Only one covariate
lambdaStar <- 1000


# Spread a Poisson amount of points randomly in the random square.
total_points <- rpois(1, lambdaStar)
random_points <- cbind(runif(total_points), runif(total_points))
grid_size <- 50

# Find covariate values to explain the species occurrence.
# We give them a Gaussian spatial structure.
reg_grid <- as.matrix(expand.grid(seq(0, 1, len = grid_size), seq(0, 1, len = grid_size)))
#### Next is commented to save time. Uncomment to replicate results
Z <- mvrnorm(1, rep(0, total_points + grid_size * grid_size), 3 * exp(-as.matrix(dist(rbind(random_points, reg_grid))) / 0.2))
Z1 <- Z[1:total_points]; Z2 <- Z[-(1:total_points)]

# Thin the points by comparing the retaining probabilities with uniforms
# in the log scale to find the occurrences
#### Next is commented to save time. Uncomment to replicate results
occurrences <- log(runif(total_points)) <= pnorm(beta[1] + beta[2] * Z1, log.p = TRUE)
n_occurrences <- sum(occurrences)
occurrences_points <- random_points[occurrences,]
occurrences_Z <- Z1[occurrences]

# Find covariate values to explain the observation bias.
# Additionally create a regular grid to plot the covariate later.
#### Next is commented to save time. Uncomment to replicate results
W <- mvrnorm(1, rep(0, n_occurrences + grid_size * grid_size), 2 * exp(-as.matrix(dist(rbind(occurrences_points, reg_grid))) / 0.3))
W1 <- W[1:n_occurrences]; W2 <- W[-(1:n_occurrences)]

# Find the presence-only observations.
#### Next is commented to save time. Uncomment to replicate results
po_sightings <- log(runif(n_occurrences)) <= pnorm(delta[1] + delta[2] * W1, log.p = TRUE)
n_po <- sum(po_sightings)
po_points <- occurrences_points[po_sightings, ]
po_Z <- occurrences_Z[po_sightings]
po_W <- W1[po_sightings]

bayesPO_sim <- list(po = cbind(po_Z, po_W), bkg = cbind(Z2, W2), grid = reg_grid,
                    po_sightings = po_sightings,
                    occurrences_points = occurrences_points)
usethis::use_data(DATASET, overwrite = TRUE)
