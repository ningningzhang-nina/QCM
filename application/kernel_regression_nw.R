mNW <- function(x, X, Y, h, K = dnorm) {
  
  # Arguments
  # x: evaluation points
  # X: vector (size n) with the predictors
  # Y: vector (size n) with the response variable
  # h: bandwidth
  # K: kernel
  
  # Matrix of size n x length(x)
  Kx <- sapply(X, function(Xi) K((x - Xi) / h) / h)
  
  # Weights
  W <- Kx / rowSums(Kx) # Column recycling!
  
  # Means at x ("drop" to drop the matrix attributes)
  drop(W %*% Y)
  
}

# Objective function
cvNW <- function(X, Y, h, K = dnorm) {
  
  sum(((Y - mNW(x = X, X = X, Y = Y, h = h, K = K)) /
         (1 - K(0) / colSums(K(outer(X, X, "-") / h))))^2)
  # Beware: outer() is not very memory-friendly!
  
}

# Find optimum CV bandwidth, with sensible grid
bw.cv.grid <- function(X, Y,
                       h.grid = diff(range(X)) * (seq(0.1, 0.5, l = 200))^2,
                       K = dnorm, plot.cv = FALSE) {
  
  obj <- sapply(h.grid, function(h) cvNW(X = X, Y = Y, h = h, K = K))
  h <- h.grid[which.min(obj)]
  if (plot.cv) {
    plot(h.grid, obj, type = "o")
    rug(h.grid)
    abline(v = h, col = 2, lwd = 2)
  }
  h
  
}

# Bandwidth
# example
# hCV <- bw.cv.grid(X = Z1, Y = Y11, plot.cv = TRUE)

