library('lbfgs')
df <- read.csv("SeoulBikeData.csv")
predict <- function(x) {
    y <- df$RentedBikeCount
    yhat <- (df$Temperature * x[1] + df$Humidity * x[2]) * df$Hour *x[4] + x[3] 
    exp(yhat)
}
objective <- function(x) {
    y <- df$RentedBikeCount
    yhat <- (df$Temperature * x[1] + df$Humidity * x[2]) * df$Hour *x[4] + x[3] 
    mu <- exp(yhat)
    sum(mu - y * yhat - y)
}
gradient <- function(x) { ## Gradient of 'fr'
    y <- df$RentedBikeCount
    yhat <- (df$Temperature * x[1] + df$Humidity * x[2]) * df$Hour * x[4] + x[3] 
    mu <- exp(yhat)
    res <- (mu - y)
    g <- c(0, 0, 0, 0)
    g[1] <- sum(df$Temperature * df$Hour * x[4] * res)
    g[2] <- sum(df$Humidity * df$Hour * x[4] * res)
    g[3] <- sum(res)
    g[4] <- sum((df$Temperature * x[1] + df$Humidity * x[2]) * df$Hour * x[4] * res)
    g
}

output <- lbfgs(objective, gradient, c(-0.1, 0.26, 1.58, -0.1), linesearch_algorithm='LBFGS_LINESEARCH_BACKTRACKING_WOLFE')
output
