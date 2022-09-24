
# I.)

# We  wish to fit a bivariate normal distribution with mean m and variance Sigma to the scallop locations (longitude and latitude)
scallops <- read.csv('C:/Users/MAMSL/Downloads/Scallops.csv')
m <- colMeans(scallops[,c("long","lat")]) # Sample mean
Sigma<- var(scallops[,c("long","lat")])
Sigma

x<- seq(-74,-71, length=31)
y<- seq(38, 41, length=31)

dens <- matrix(0,31,31)

for (i in 1:31){
  for (j in 1:31){
    dens[i,j]<- 1/(2*pi*sqrt(det(Sigma)))*exp(-1/2* c(x[i]-m[1], y[j]-m[2])%*%solve(Sigma)%*%c(x[i]-m[1], y[j]-m[2])) 
  }
}

# 3D view of the density 

persp(x,y,dens, phi = 40, theta = 60)

# 2D dimensional viex of the density represnted using the contour function
par(mfrow=c(1,1))
contour(x, y, dens)
points(scallops$long, scallops$lat)
# The BNV seems not to be perfectly adequate as it cannot account for (and gets biased by) the small branch in the NW.
# Relative to the fitted distribution, some points in the NW and the SE appear to be outliers, but
# due to the reason above one has to be careful with this interpretation.



# II.) MAHALANOBIS DISTANCE AND OUTLIER DETECTION
# Here we work on a new dataset called  engine 

engine <-  read.table('C:/Users/MAMSL/Downloads/engine.dat', header  =  TRUE) # download the data to R



# a) Here we compute the sample mean and the sample variance matrix rounded two digits
M    <- round(colMeans(engine),2)
S    <- round(var(engine),2)

M
#  CO   HC  NOX 
# 7.96 0.55 1.33 

S
#       CO    HC   NOX
# CO  27.67  0.80 -1.75
# HC   0.80  0.03 -0.05
# NOX -1.75 -0.05  0.23

# Two have an idea of which observation can be considered as a possible outlier, we produce a pair plot
# of the data

pairs(engine)

plot(engine$CO, engine$HC)
identify(engine$CO, engine$HC)

plot(engine$CO, engine$NOX)
identify(engine$CO, engine$NOX)

plot(engine$HC, engine$NOX)
identify(engine$HC, engine$NOX)

# Observations 34,45, and 39  appear to be possible outliers


# Here, we compute the squared mahalanobis distance for each obersvation of the dataset
d <- mahalanobis(engine, M, S)

# For example, the squared mahalnobis distance for the 5th observation is 2.38544
d[5]
# 2.385444 

# The squared mahalnobis distance can also be computed using its mathematical formula as below: 
eng5 <- as.numeric(engine[5,])
(eng5-M)%*%solve(S)%*%(eng5-M)
#        [,1]
#[1,] 2.385444

# d) For all 46 observations of our data set, we carry out a statistical test of the hypothesis H0: Case i is not outlying. Here we give
#case numbers for which this null hypothesis is rejected at the 5% and 2.5% level of significance,

# This mahalanobis test for outlier detection consists first in assuming that the observations of the variables forming the engine data sets have been 
# generated according to a mutltivariate normal distribution, with unknown population mean and unknown population variance. 
# H0 is rejected for observation i if its squared mahalanobis distance is greater that the critical value of the chi squared
# distribution with 3 degrees of freedom
# 5% of significance
which(d>   qchisq(0.95,3))
# 30 34 35 39 
#  Observation 30 might be surprising.

# 2.5% of significance
which(d>   qchisq(0.975,3))
# 34 35 39

# Observation 34,35,39 seem to be outliers assuming that the data is mutlivariate normal. To check this empirically,
# it is important to notice that the squared mahalnobis distance for some random vector of q variables follow a chi squared distribution with q degrees of freedom,
# if this random vecotr is Multivariate normal. Thus, one can check the assumption of mutlivariate normality of the data, by cheking
# if the mahalnobis distances are generated from a chi square distribution with q degrees of freedom.
# To do so one can make a quantile vs quantile plot.

plot(qchisq( (1:46-0.5)/46,3),sort(d)) 
abline(a=0,b=1)
# clearly some deviation from the straight line,
# especially in the central part (observations are too small)
# and the boundary (observations are too large ==> outliers).
# the d_i  have also some discrete (stepwise) character.


# The above mahalanobis outlier test seen before has some drawbacks, first we must assume that the data is generated from
# a multivariate normal distribution which is not always the case. Estimates for m and Σ are often calculated using same values of x we are now potentially thinking of removing 
#- (potential) outliers are already affecting our choices through using those estimates for computing mahalnobis distances!
# One can show mathematically that the effective type I error rate (i.e., the probability of at least one observation is incorrecly classified
# as an outlier ) is nα where n the number of observations of our dataset and α, the level of significance. Thus, the more observations you have and the more tests you make,
# and the probability of making a mistake increases. Thus, a possible solution is applying a bonferoni correction which replaces α with α/n.

# After Bonferroni correction:
which(d>   qchisq(1-0.05/46,3))
# only one of the observations is identified as outliers
