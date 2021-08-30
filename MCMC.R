library(ggplot2)

## MCMC algorithm for Cauchy distribution ##

# following the Roberts, Gelman and Gilks paper on weak convergence (1997),
# we pilot run the Metropolis algorithm up to and including calculating the
# acceptance probability, to find the optimal value for the standard deviation,
# in order to give us 0.234 acceptance rate

# define the metropolis algorithm up to the acceptance rate to return a mean value of the rates
metropolis_optimal <- function(N, sd) {
  # generate the first value
  X0 <- runif(1, min=0, max=1)
  # create a vector to store the states of the chain
  X <- X0
  accept_rates <- NULL
  # start the algorithm here
  for (i in 1:N) {
    # sample from a normal distribution
    Y <- rnorm(1, mean=X, sd=sd)
    # calculate the acceptance probability
    r <- min((1+X^2)/(1+Y^2), 1)
    # store the accept rate
    accept_rates <- c(accept_rates, r)
    # perform the accept-reject step
    if (runif(1,min=0,max=1) < r) {
      X <- Y
    }
  }
  # return the mean accept rate
  return(mean(accept_rates))
}

# define a function that will run metropolis optimal for different values of 
# the standard variation to find approximately the one giving an acceptance rate
# closest to the optimal one of 0.234
optimal_sd <- function(K, N, start_sd, end_sd) {
  # define an interval with changing values for the standard deviation
  interval_sd <- seq(from = start_sd, to = end_sd, length.out = K)
  # set the optimal rate
  optimal_rate <- 0.234
  # run the first chain
  sd <- metropolis_optimal(N, interval_sd[1])
  # calculate the difference
  diff <- abs(sd - optimal_rate)
  
  for (i in 2:length(interval_sd)) {
    # run a chain
    sd_candidate <- metropolis_optimal(N, interval_sd[i])
    # check if the new parameter is better
    if (abs(sd_candidate - optimal_rate) < diff) {
      diff <- abs(sd_candidate - optimal_rate)
      i_best <- i
    }
  }
  return(interval_sd[i_best])
}

# mean average of several pilot runs of optimal_sd function
final_sd <- function(M, K, N, start_sd, end_sd) {
  results <- NULL
  for (s in 1:M) {
    results <- c(results, optimal_sd(K, N, start_sd, end_sd))
  }
  return(mean(results))
}

# calculate the optimal sd
final_optimal_sd <- final_sd(10, 500, 10000, 1, 20)

metropolis_normal <- function(N, sd) {
  # generate the first value
  X0 <- runif(1, min=0, max=1)
  # create a vector to store the states of the chain
  X <- X0
  
  # start the algorithm here
  for (i in 1:N) {
    # sample from a normal distribution
    Y <- rnorm(1, mean=X[i], sd=sd)
    # calculate the acceptance probability
    r <- min((1+X^2)/(1+Y^2), 1)
    # perform the accept-reject step
    if (runif(1,min=0,max=1) < r) {
      X <- c(X, Y)
    } else {
      X <- c(X, X[i])
    }
  }
  # store the chain in a list
  return(list("time" = 1:length(X), "states" = X))
}

# run the MH algorithm
mh <- metropolis_normal(50000, sd=final_optimal_sd)

# store the results in a data frame
output <- data.frame(S = mh$states, T = mh$time)

# plot the results for all three compartments
ggplot(data = output, aes(x=T)) +
  geom_line(aes(y=S)) +
  labs(x = "Time", y = "States")

# produce a histogram of states and overplot the Cauchy distribution on top
hist(mh$states, freq = F, ylim = c(0,0.35), xlim = c(-6,6), breaks = 150, main = "Histogram of the state space of the Markov chain", xlab = "Values")
curve(1/(pi*(1+x^2)), from=-6, to=6, add = T, col="red", lwd=2)
