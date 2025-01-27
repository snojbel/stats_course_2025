# Multivariate linear models
# -> Using more than one predictor variable to model outcome.

# Set up   ------------------------
library(rethinking)


# Standerdize predictor ----------------
data("WaffleDivorce")
d <- WaffleDivorce

d$MedianAgeMarriage.s <- with(d,(MedianAgeMarriage-mean(MedianAgeMarriage))/sd(MedianAgeMarriage))

# Fit model

m5.1 <- map(
            alist(
                  Divorce ~ dnorm(mu, sigma),
                  mu <- a + bA * MedianAgeMarriage.s,
                  a ~ dnorm(10, 10), 
                  bA ~ dnorm(0, 1),
                  sigma ~ dunif(0, 10)
            ), data = d)
precis(m5.1)

# Intervals

MAM.seq <- seq(from = -3, to = 3, length.out = 30)
mu <- link(m5.1, data = data.frame(MedianAgeMarriage.s = MAM.seq))
mu.PI <- apply(mu, 2, PI)

# plot

plot(Divorce ~ MedianAgeMarriage.s, data = d, col = rgb(0.1, 0.5, 0.3, 0.2), 
     pch = 16, xlab = " Standerdized Median age at Marriage")
abline(m5.1)
shade(mu.PI, MAM.seq)

# Marriage rate and Divorce

d$Marriage.s <- with(d, (Marriage-mean(Marriage))/sd(Marriage))

m5.2 <- quap(alist(
                  Divorce ~ dnorm(mu, sigma),
                  mu <- a + bR * Marriage.s, 
                  a ~ dnorm(10, 10),
                  bR ~ dnorm(0, 1), 
                  sigma ~ dunif(0, 10))
            , data = d)
 precis(m5.2)

# Hard to compare parameter means between bivariate regressions so we should build
# a multivariate regression and measure partial value of each predictor
 
# What is the predictive value of a variable, once I already know all of the other
# predictor variables?
 
# x. What is the increased value in knowing age at marriage if I already know marriage rate

# Model
 
m5.3 <- quap(alist(
                    Divorce ~ dnorm(mu, sigma),
                    mu <- a + bA * MedianAgeMarriage.s + bR * Marriage.s,
                    a ~ dnorm(10, 10),
                    bA ~ dnorm(0, 1),
                    bR ~ dnorm(0, 1),
                    sigma ~ dunif(0, 10)
                    ), data = d)
precis(m5.3) 
plot(m5.3)

# Shows that once we know median age at marriage marriage rate offers little or 
# no additional predictive power. 

# Predictor residual plots ----------------------------------------------------
# The average prediction error when we use all of the other predictor variables to
# model the predictor of intereset

m5.4 <- quap(
             alist(
                   Marriage.s ~ dnorm(mu, sigma),
                   mu <- a + bA*MedianAgeMarriage.s,
                   a ~ dnorm(0, 10), 
                   bA ~ dnorm(0, 1),
                   sigma ~ dunif(0, 10)
             ), data = d )


# Compute expected value at each state
mu <- coef(m5.4)["a"] + coef(m5.4)["bA"] * d$MedianAgeMarriage.s 
# Compute residuals for each state
m.resid <- d$Marriage.s - mu

plot(Marriage.s ~ MedianAgeMarriage.s, d, col = rangi2, pch = 16)
abline(m5.4)
#loop over states
for(i in 1:length(m.resid)){
  x <- d$MedianAgeMarriage.s[i]  # x coordinate
  y <- d$Marriage.s[i]           # y coordinate 
  # draw line
  lines(c(x, x), c(mu[i], y), lwd = 0.5, col = col.alpha("black", 0.7))
  # each line goes from point (x, mu) to (x, y)
}

# Shows variation in marriage rate that is not explained by the variable age at marriage

# Multivarate models regressions show: remaining association of each predictor after
# already knowing the other


# Counterfactual plots --------------------------------
# Used to see how prediction change as you change only one predictor at a time
# Help understand implications of model

# prepare counterfactual data
A.avg <- mean(d$MedianAgeMarriage.s)
R.seq <- seq(from = -3, to = 3, length.out = 30)
pred.data <- data.frame(
    Marriage.s <- R.seq,
    MedianAgeMarriage.s <- A.avg
)

# compute counterfactual mean divorce
mu <- link(m5.3, data = pred.data)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

# simulate counterfactual outcomes
R.sim <- sim(m5.3, data = pred.data, n = 1e4)
R.PI <- apply(R.sim, 2, PI)

# Display predictions

plot(Divorce ~ Marriage.s, data = d, type = "n")
mtext("MedianAgeMarriage.s = 0")
lines(R.seq, mu.mean)
shade(mu.PI, R.seq)
shade(R.PI, R.seq)

# This uses the model to predict divorce when medianage at marriage is constant at zero
# but marriage rate is allowed to vary.

# We can also predict divorce rate when
# we allow age at marriage to vary and keep marriage rate constant:
R.avg <- mean(d$Marriage.s)
A.seq <- seq(from = -3, to = 3, length.out = 30)
pred.data2 <- data.frame(
              Marriage.s <- R.avg,
              MedianAgeMarriage.s <- A.seq
              )

mu <- link(m5.3, data = pred.data2)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

A.sim <- sim(m5.3, data = pred.data2, n = 1e4)
A.PI <- apply(A.sim, 2, PI)

plot(Divorce ~ MedianAgeMarriage.s, data = d, type = "n")
mtext("Marriage.s = 0")
lines(A.seq, mu.mean)
shade(mu.PI, A.seq)
shade(A.PI, A.seq)

# Posterior prediction plots ---------------------------------------------------

# Call link using orginal data as data:

mu <- link(m5.3)

# Summarize samples across cases
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

# simulate observations

divorce.sim <- sim(m5.3, n = 1e4)
divorce.PI <- apply(divorce.sim, 2, PI)

plot(mu.mean ~ d$Divorce, col = rangi2, ylim = range(mu.PI), xlab = "Observed divorce"
      , ylab = "Predicted Divorce")
abline(a = 0, b = 1, lty = 2)
for(i in 1:nrow(d)){
  lines(rep(d$Divorce[i], 2), c(mu.PI[1, i], mu.PI[2, i]), col = rangi2)
}

# Under predicts with high divorce rates and over predicts for states with low divorce rate

identify(x = d$Divorce, y = mu.mean, labels = d$Loc, cex = 0.8) # Press in plot window to make point get its label

# Compute residuals

divorce.resid <- d$Divorce - mu.mean

# Order
o <- order(divorce.resid)

# Make plot

dotchart(divorce.resid[o], labels = d$Loc[o], xlim = c(-6, 5), cex = 0.6)
abline(v = 0, col = col.alpha("black", 0.2))
for(i in 1:nrow(d)){
  j <- o[i] # Which state in order
  lines(d$Divorce[j]-c(mu.PI[1, j], mu.PI[2, j]), rep(i, 2))
  points(d$Divorce[j]-c(divorce.PI[1,j], divorce.PI[2, j]), rep(i, 2), pch = 3
         , cex = 0.6, col = "gray")
}

# Masked relationship ----------------------------------------------------------

# Sometimes two predictors can be correlated with the outcome and mask the existance of any correlation
# because one is positively correlated and the other negativly correlated

data(milk)
d <- milk
str(d)

m5.5 <- quap(alist(
                   kcal.per.g ~ dnorm(mu, sigma),
                   mu <- a + bn*neocortex.perc,
                   a ~ dnorm(0, 100),
                   bn ~ dnorm(0, 1),
                   sigma ~ dunif(0, 1)
                   ), data = d)
# Doesnt work because of NA's

dcc <- d[complete.cases(d), ]

m5.5 <- quap(alist(
                    kcal.per.g ~ dnorm(mu, sigma),
                    mu <- a + bn*neocortex.perc,
                    a ~ dnorm(0, 100),
                    bn ~ dnorm(0, 1),
                    sigma ~ dunif(0, 1)
                  ), data = dcc)
precis(m5.5, digits = 3)
plot(kcal.per.g ~ neocortex.perc, data = dcc, pch = 16, col = rangi2)

# Not a very strong if any correlation between these two

np.seq <- 0:100
pred.data <- data.frame(neocortex.perc = np.seq)

mu <- link(m5.5, data = pred.data, n = 1e4)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

lines(np.seq, mu.mean)
lines(np.seq, mu.PI[1,], lty = 2)
lines(np.seq, mu.PI[2, ], lty = 2)

# What about mean mass of species?
# Scaling measurements like body mass are often related by magnitudes to other variables
# Taking the log of a measure translates into this magnitude.
# By taking the log of body mass we say that we suggest that the the magnitude
# of a mothers body mass is related to milk energy, linearly. 

dcc$log.mass <- log(dcc$mass)

m5.6 <- quap(alist(
                    kcal.per.g ~ dnorm(mu, sigma),
                    mu <- a + bm * log.mass,
                    a ~ dnorm(0, 100),
                    bm ~ dnorm(0, 1),
                    sigma ~ dunif(0, 1)
                    ), data = dcc)
precis(m5.6, digits = 3)

logmass.seq <- seq(-5, 5 , length.out = 30)
pred.data <- data.frame(log.mass = logmass.seq)
mu <- link(m5.6, data = pred.data)

mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

plot(kcal.per.g ~ log.mass, data = dcc, pch = 16, col = rangi2)
lines(logmass.seq, mu.mean)
lines(logmass.seq, mu.PI[1, ], lty = 2)
lines(logmass.seq, mu.PI[2, ], lty = 2)

# Adding both predictor variables

m5.7 <- quap(alist(
                    kcal.per.g ~ dnorm(mu, sigma),
                    mu <- a + bm*log.mass + bn*neocortex.perc,
                    a ~ dnorm(0, 100),
                    bm ~ dnorm(0, 1),
                    bn ~ dnorm(0, 1),
                    sigma ~ dunif(0, 1)
                    ), data = dcc)
precis(m5.7, digits = 3)

# By incorporating both the correlation of both to kcal increases. 

# Lets look at some counterfactual plots

mean.log.mass <- mean(log(dcc$mass))
np.seq <- 0:100
pred.data <- data.frame(
              neocortex.perc = np.seq,
              log.mass = mean.log.mass)
mu <- link(m5.7, data = pred.data, n = 1e4)

mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

plot(kcal.per.g ~ neocortex.perc, data = dcc, type = "n")
lines(np.seq, mu.mean)
lines(np.seq, mu.PI[1, ], lty = 2)
lines(np.seq, mu.PI[2, ], lty = 2)

# Equivalent plot but with changing body mass

mean.neocortex <- mean(dcc$neocortex.perc)
logmass.seq <- seq(-5, 5 , length.out = 30)
pred.data <- data.frame(
  neocortex.perc = mean.neocortex,
  log.mass = logmass.seq)
mu <- link(m5.7, data = pred.data, n = 1e4)

mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

plot(kcal.per.g ~ log.mass, data = dcc, type = "n")
lines(logmass.seq, mu.mean)
lines(logmass.seq, mu.PI[1, ], lty = 2)
lines(logmass.seq, mu.PI[2, ], lty = 2)

# The two predictors are positively correlated with another, so they tend to cancel each others effect out
# The model asks if species that have high neocortex percent for their body mass also has
# higher milk energy. (and vice versa)

# Multicollinearity ------------------------------------------------------------

# Legs and height

N <- 100
height <- rnorm(N, 10, 2)   # These people are quite short
leg_prop <- runif(N, min = 0.4, max = 0.5)
leg_left <- leg_prop*height + rnorm(N, 0, 0.02) # leg length with some developmental differences
leg_right <- leg_prop*height + rnorm(N, 0, 0.02)

d <- data.frame(height, leg_left, leg_right)

m5.8 <- quap(alist(
                   height ~ dnorm(mu, sigma),
                   mu <- a + bl*leg_left + br*leg_right,
                   a ~ dnorm(10, 100),
                   bl ~ dnorm(2, 10),
                   br ~ dnorm(2, 10),
                   sigma ~ dunif(0, 10)
                   ), data = d)
precis(m5.8)
plot(precis(m5.8))

post <- extract.samples(m5.8)
plot(bl ~ br, post)

# Because leg lengths contain almost the same information is like including the same predictor twice
# which equates to mu <- a + bl*x + br*x. which would be a + (bl + br) *x
# so the some of bl and br is actually the good estimate

sum_blbr <- post$bl + post$br
dens(sum_blbr)

# Only using one leg will produce a much better estimate

m5.9 <- quap(alist(
  height ~ dnorm(mu, sigma),
  mu <- a + bl*leg_left ,
  a ~ dnorm(10, 100),
  bl ~ dnorm(2, 10),
  sigma ~ dunif(0, 10)
), data = d)
precis(m5.9)

# Multicollinerity in MILK

d <- milk

# calories and fat correlation
m5.10 <- quap(alist(
                    kcal.per.g ~ dnorm(mu, sigma),
                    mu <- a + bf * perc.fat,
                    a ~ dnorm(0, 10),
                    bf ~ dnorm(0, 1),
                    sigma ~ dunif(0, 10)
                    ), data = d)


# calories and lactose correlation
m5.11 <- quap(alist(
                    kcal.per.g ~ dnorm(mu, sigma),
                    mu <- a + bl * perc.lactose,
                    a ~ dnorm(0, 10),
                    bl ~ dnorm(0, 1),
                    sigma ~ dunif(0, 10)
                  ), data = d)
precis(m5.10, digits = 3)
precis(m5.11, digits = 3)

# Both

m5.12 <- quap(alist(
                    kcal.per.g ~ dnorm(mu, sigma),
                    mu <- a + bl * perc.lactose + bf*perc.fat,
                    a ~ dnorm(0, 10),
                    bl ~ dnorm(0, 1),
                    bf ~ dnorm(0, 1),
                    sigma ~ dunif(0, 10)
                  ), data = d)
precis(m5.12, digits = 3)

# These two have essentially the same information, so including both, the posterior distribution will
# start describing the long distribution where they are both possible. They are substitues for each other

pairs(~ kcal.per.g + perc.fat + perc.lactose, data = d, col = rangi2, pch = 16)

# after knowing one predictor knowing the other does not add any new additional information.

cor(d$perc.fat, d$perc.lactose)

# Check predictor variables against one another in a pairs plot to see if multicollinearity could
# be a problem from your model

# Post treatment Bias ----------------------------------------------------------
# Including variables that are consequences of other variables

# Fungus and plant example

N <- 100
h0 <- rnorm(N, 10, 2)
treatment <- rep(0:1, each = N/2)
fungus <- rbinom(N, size = 1, prob = 0.5 - treatment*0.4)
h1 <- h0 + rnorm(N, 5-3*fungus)

d <- data.frame(h0, treatment, fungus, h1)
str(d)

m5.13 <- quap(alist(
                    h1 ~ dnorm(mu, sigma),
                    mu <- a + bh*h0 + bt*treatment + bf*fungus,
                    a ~ dnorm(0, 100),
                    c(bh, bt, bf) ~ dnorm(0, 10),
                    sigma ~ dunif(0, 10)
                    ), data = d)
precis(m5.13)

# Treatment has very little effect, but thats because we controlled for fungus in the model. We
# are essentially asking, once we already know if it developed fungus, how much more information
# does knowing treatment give.

m5.14 <- quap(alist(
                    h1 ~ dnorm(mu, sigma),
                    mu <- a + bh*h0 + bt*treatment,
                    a ~ dnorm(0, 100),
                    c(bh, bt) ~ dnorm(0, 10),
                    sigma ~ dunif(0, 10)
                  ), data = d)
precis(m5.14)


# Categorical variables --------------------------------------------------------

data(Howell1)
d <- Howell1
str(d)

m5.15 <- quap(alist(
                    height ~ dnorm(mu, sigma),
                    mu <- a + bi*male,
                    a ~ dnorm(178, 100),
                    bi ~ dnorm(0, 10),
                    sigma ~ dunif(0, 50)
                    ), data = d)
precis(m5.15)

# Get the interval for male height (a will here be interval for female height)

post <- extract.samples(m5.15, n = 1e4)
mu.male <- post$a + post$bi
PI(mu.male)

# MANY categories

data(milk)
d <- milk
unique(d$clade)

# creating dummy data

d$clade.NWM <- ifelse(d$clade=="New World Monkey", 1, 0)
d$clade.OWM <- ifelse(d$clade=="Old World Monkey", 1, 0)
d$clade.S <- ifelse(d$clade=="Strepsirrhine", 1, 0)

m5.16 <- quap(alist(kcal.per.g ~ dnorm(mu, sigma),
                    mu <- a + bn*clade.NWM + bo*clade.OWM + bs*clade.S,
                    a ~ dnorm(0, 10),
                    c(bn, bo, bs) ~ dnorm(0, 1),
                    sigma ~ dunif(0, 10)
                    ), data = d)
precis(m5.16)
# this actually defines four different linear models mu = a, mu = a + bn, mu = a + bo, mu = a +bs
# as the other variables will be turned off when they arent relevant.

# Get averages
post <- extract.samples(m5.16)
str(post)

mu.ape <- post$a
mu.NWM <- post$a + post$bn
mu.OWM <- post$a + post$bo
mu.S   <- post$a + post$bs

#summarize using precis
precis(data.frame(mu.ape, mu.NWM, mu.OWM, mu.S))

diff.NWM.OWM <- mu.NWM - mu.OWM
quantile(diff.NWM.OWM, probs = c(0.025, 0.5, 0.975))

# Using index instead of dummy variables

d$clade_id <- coerce_index(d$clade)   # creates index numbers and indexs based on clade

m5.16_alt <- quap(
                  alist(
                    kcal.per.g ~ dnorm(mu, sigma),
                    mu <- a[clade_id], 
                    a[clade_id] ~ dnorm(0, 10),
                    sigma ~ dunif(0, 10)
                  ), data = d)
precis(m5.16_alt, depth = 2)

# Ordinary Least Squares (lm) models --------------------------------------------
# provided you are okay with flat priors they will achieve the same results as the previous examples

m5.17 <- lm(kcal.per.g ~ 1 + perc.fat, data = d)
m5.18 <- lm(kcal.per.g ~ 1 + clade.S + clade.OWM + clade.NWM, data = d)
m5.18 <- lm(kcal.per.g ~ 1 + clade, data = d)  # No need to split categorical variables into dummy variables, unless you want that controll

# Intercepts are automatically added even without the 1, if you dont want one add - 1 at the end.
# functions like lm are also sensitive to transformations of data within itself, so do these transformations
# yourself before you give it to the model. 

# Practice -------------

# 5E
# 1) Model 2 and model 4,

# 2) Animal diversity is linearly related to latitude, but only after accounting for plant diversity model:
# suggests that there is no direct relationship between latitude and animal diversity but rather plant diversity
# is acting as a collider for the two
mAnimalDiversity <- quap(alist(ad ~ dnorm(mu, sigma),
                               mu <- a + bl * latitude + bp * plant.diversity,
                               a ~ dnorm(0, 10),
                               bl ~ dnorm(0, 10),
                               bp ~ dnorm(0, 10),
                               sigma ~ dunif(0, 10)
                               ), data = dataframe)

# 3) Neither amount of funding nor size of laboratory is by itself a good predictor of time to PhD degree,
# but together these variables are both positively associated with time to degree.
# If they are negatively correlated with one another, then considering each alone may miss the
# positive relationships with the outcome, both parameters positive
mPhD.time <- quap(alist(time ~ dnorm(mu, sigma),
                               mu <- a + bf * funding + bls * lab.space,
                               a ~ dnorm(0, 10),
                               bls ~ dnorm(0, 10),
                               bf ~ dnorm(0, 10),
                               sigma ~ dunif(0, 10)
), data = dataframe)

# 4) 1, 3, 4

# 5M
# 1-2 paper
# 3) More divorce means more singles, which means more opportunities to get married. Could be tested by 
# seeing if there is a correlation between number of singles and divorce rate?

# 4)
data("WaffleDivorce")
d <- WaffleDivorce
mormon <- read.csv("mormon-population-by-state-2024.csv")
mormon <- mormon[, 1:2]
d <- merge(d, mormon, by.x = "Location", by.y = "state")

# Standardize
d$MedianAgeMarriage.s <- with(d,(MedianAgeMarriage-mean(MedianAgeMarriage))/sd(MedianAgeMarriage))
d$MormonPopulation2022.s <- with(d,(MormonPopulation2022-mean(MormonPopulation2022))/sd(MormonPopulation2022))
d$Marriage.s <- with(d, (Marriage-mean(Marriage))/sd(Marriage))

mMormon <- quap(alist(
                Divorce ~ dnorm(mu, sigma),
                mu <- a + bA * MedianAgeMarriage.s + bM * MormonPopulation2022.s + bMa * Marriage.s,
                a ~ dnorm(10, 10),
                bA ~ dnorm(0, 1),
                bM ~ dnorm(0, 1),
                bMa ~ dnorm(0, 1),
                sigma ~ dunif(0, 10)
              ), data = d)

precis(mMormon)

# 5)
# I would start by seeing if the relationship between gas and obesity is changed when filtering for time per
# week of exercise. Then you would isolate the effect through restuarant visits. And then test the opposite.
# I would also try testing with both to see if there is any apperant direct effect of gas on obesity because this
# would suggest possibility of some other confounder as there is little reason to believe that gas price directly
# lowers obesity. 

#5H
data("foxes")
d <- foxes

# 1.1) 
d$area.s <- with(d, (area - mean(area))/sd(area))
mFox.1 <- quap(alist(weight ~ dnorm(mu, sigma),
                     mu <- a + bt * area.s,
                     a ~ dnorm(0, 10),
                     bt ~ dnorm(0, 1), 
                     sigma ~ dunif(0, 10)
                     ), data = d)
plot(precis(mFox.1))

par(mfrow = c(2, 2))

range(d$area.s)
area.seq <- seq(from = -3, to = 3, length.out = 30)

simFox <- sim(mFox.1, data = list(area.s = area.seq))
weight.PI <- apply(simFox, 2, PI)

mu <- link(mFox.1, data = data.frame(area.s = area.seq))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

# plot

plot(weight ~ area.s, data = d, col = rgb(0.1, 0.5, 0.3, 0.2), 
     pch = 16, xlab = " Standardized territory Area")
lines(area.seq, mu.mean)
shade(mu.PI, area.seq)
shade(weight.PI, area.seq)

# 1.2)

d$group.s <- with(d, (group - mean(group))/sd(group))
mFox.2 <- quap(alist(weight ~ dnorm(mu, sigma),
                     mu <- a + bg * group.s,
                     a ~ dnorm(0, 10),
                     bg ~ dnorm(0, 1), 
                     sigma ~ dunif(0, 10)
), data = d)

plot(precis(mFox.2))

range(d$group.s)
group.seq <- seq(from = -3, to = 3, length.out = 30)

simFox <- sim(mFox.2, data = list(group.s = group.seq))
weight.PI <- apply(simFox, 2, PI)

mu <- link(mFox.2, data = data.frame(group.s = group.seq))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

plot(weight ~ group.s, data = d, col = rgb(0.1, 0.5, 0.3, 0.2), 
     pch = 16, xlab = " Standardized group size")
lines(group.seq, mu.mean)
shade(mu.PI, group.seq)
shade(weight.PI, group.seq)

# Doesn't seem like either of them are vary important.

# 2)

mFox.3 <- quap(alist(weight ~ dnorm(mu, sigma),
                     mu <- a + bg * group.s + ba * area.s,
                     a ~ dnorm(0, 10),
                     bg ~ dnorm(0, 1),
                     ba ~ dnorm(0, 1),
                     sigma ~ dunif(0, 10)
), data = d)

plot(precis(mFox.3))

# static area
group.seq <- seq(from = -3, to = 3, length.out = 30)
mean.area <- mean(d$area.s)

simFox <- sim(mFox.3, data = list(group.s = group.seq, area.s = mean.area))
weight.PI <- apply(simFox, 2, PI)

mu <- link(mFox.3, data = data.frame(group.s = group.seq, area.s = mean.area))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

plot(weight ~ group.s, data = d, col = rgb(0.1, 0.5, 0.3, 0.2), 
     pch = 16, xlab = " Standardized group size")
lines(group.seq, mu.mean)
shade(mu.PI, group.seq)
shade(weight.PI, group.seq)

# static group
area.seq <- seq(from = -3, to = 3, length.out = 30)
mean.group <- mean(d$group.s)

simFox <- sim(mFox.3, data = list(area.s = area.seq, group.s = mean.group))
weight.PI <- apply(simFox, 2, PI)

mu <- link(mFox.3, data = data.frame(area.s = area.seq, group.s = mean.group))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

plot(weight ~ area.s, data = d, col = rgb(0.1, 0.5, 0.3, 0.2), 
     pch = 16, xlab = " Standardized area size")
lines(area.seq, mu.mean)
shade(mu.PI, area.seq)
shade(weight.PI, area.seq)

# They both seem to influence but not a huge amount, they are probably correlated.

# 3)

d$avgfood.s <- with(d, (avgfood - mean(avgfood))/sd(avgfood))

mFox.4 <- quap(alist( weight ~ dnorm(mu, sigma),
                      mu <- a + bf*avgfood.s + bg*group.s, 
                      a ~ dnorm(0, 10),
                      bf ~ dnorm(0, 1),
                      bg ~ dnorm(0, 1),
                      sigma ~ dunif(0, 10)
                      ), data = d)
precis(mFox.4)

# Plot
avgfood.seq <- seq(from = -3, to = 3, length.out = 30)
mean.group <- mean(d$group.s)

simFox <- sim(mFox.4, data = list(avgfood.s = avgfood.seq, group.s = mean.group))
weight.PI <- apply(simFox, 2, PI)

mu <- link(mFox.4, data = data.frame(avgfood.s = avgfood.seq, group.s = mean.group))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

plot(weight ~ avgfood.s, data = d, col = rgb(0.1, 0.5, 0.3, 0.2), 
     pch = 16, xlab = " Standardized Average food")
lines(avgfood.seq, mu.mean)
shade(mu.PI, avgfood.seq)
shade(weight.PI, avgfood.seq)

group.seq <- seq(from = -3, to = 3, length.out = 30)
mean.avgfood <- mean(d$avgfood.s)

simFox <- sim(mFox.4, data = list(group.s = group.seq, avgfood.s = mean.avgfood))
weight.PI <- apply(simFox, 2, PI)

mu <- link(mFox.4, data = data.frame(group.s = group.seq, avgfood.s = mean.avgfood))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

plot(weight ~ group.s, data = d, col = rgb(0.1, 0.5, 0.3, 0.2), 
     pch = 16, xlab = " Standardized Group size")
lines(group.seq, mu.mean)
shade(mu.PI, group.seq)
shade(weight.PI, group.seq)

mFox.5 <- quap(alist( weight ~ dnorm(mu, sigma),
                      mu <- a + bf*avgfood.s + bg*group.s + ba*area.s, 
                      a ~ dnorm(0, 10),
                      bf ~ dnorm(0, 1),
                      bg ~ dnorm(0, 1),
                      ba ~ dnorm(0, 1),
                      sigma ~ dunif(0, 10)
), data = d)

precis(mFox.5)

plot(mFox.5)

# Avgfood and area provide the same information, are probably highly correlated


