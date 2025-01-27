# Homework week 2

# Age influences weight through height indirectly and directly through physiological 
# changes. 

# Generative model
weight_sim <- function(N, age, height){
  a <- 3
  ba <- 3
  bh <- 0.3
  #height <- rnorm(N, 20 + log(age)*40, 10)
  mu <- a + ba*log(age) + bh*height
  sigma <- 5
  weight <- rnorm(N, mu, sigma)
  return(weight)
}


# Data creation
N <- 1e4
ages <- rnorm(N, 30, 5)
heights <- rnorm(N, 140, 20) 
weights <- weight_sim(N, ages, heights)

d <- data.frame(ages, heights, weights)

plot(weights ~ ages, data = d)

mAgeHeight <- quap(alist( weights ~ dnorm(mu, sigma),
                    mu <- a + ba*log(ages) + bh*heights,
                    a ~ dnorm(0, 10),
                    ba ~ dnorm(0, 1),
                    bh ~ dnorm(0, 1),
                    sigma ~ dunif(0, 10)
                   ), data = d)
precis(mAgeHeight)


# Using Howell1 consider only people <13 years and estimate total causal effect of each
# year of growth on weight

data("Howell1")
d <- Howell1
str(d)

d13 <- d[d$age<13, ]   # filter out old people
str(d13)

d13$height.s <- (d13$height-mean(d13$height)) / sd(d13$height)
d13$age.s <- (d13$age - mean(d13$age)) / sd(d13$age)
d13$weight.s <- (d13$weight-mean(d13$weight)) / sd(d13$weight)

m13 <- quap(alist( weight.s ~ dnorm(mu, sigma),
                          mu <- a + ba*age.s,
                          a ~ dnorm(0, 1),
                          ba ~ dnorm(0, 1),
                          sigma ~ dunif(0, 10)
), data = d13)

precis(m13)
plot(precis(m13))

# Dont stratify by height because then we only see the effect of age has which does not work through height.