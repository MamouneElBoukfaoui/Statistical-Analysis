###############################################################################
# 1.1 Manipulating matrices
################################################################################

################################################################################

A <-matrix(c(3,1,4,2,0,0,1,3,4), nrow=3, byrow=TRUE) # 3 by 3 matrix
A


#    [,1] [,2] [,3]
#[1,]    3    1    4
#[2,]    2    0    0
#[3,]    1    3    4

# Determinant.
det(A)

# Transpose.
t(A)

# Inverse
solve(A)

################################################################################
# Write a function tr which calculates the trace of a square matrix, and apply it to matrix A.
# Extract the diagonal of the argument M as a vector, and then sum it.

tr<-function(M){
  Trace <- sum(diag(M))
  return(Trace)
}

tr(A)

################################################################################
# 1.2  Simulations
################################################################################
################################################################################
# Given some density function f(x,y) = 1/2 for -1 <= x <= 1 and  0 <= y <= 1 and
# 0 elsewhere.

# Generate 100 samples each of x and y from this density.
x<-runif(100,-1,1)
y<-runif(100,0,1)

################################################################################
# Produce a scatterplot of the simulated data set..
plot(x,y)

################################################################################

# Combine the x and y values and compute the sample variance matrix.
Z<-cbind(x,y)
var(Z)

# Note that you have to stack x and y next to one another to get the full
# variance matrix. This consists not only of the variances of x and y but the
# covariance between them:

#     x                         y
# x   var(x) = cov(x, x)        cov(x, y) = cov(y, x)
#
# y   cov(y, x) = cov(x, y)     var(y) = cov(y, y)

################################################################################
# Now what if we repeat the experiment for a number of different values of n, both smaller and larger than 100,
# and estimate the variance in each case.

# Generate 10,000 samples each of x and y.
x<-runif(10000,-1,1)
y<-runif(10000,0,1)

# Plot them.
plot(x,y)

# Combine the x and y values and compute the variance matrix.
Z<-cbind(x,y)
var(Z)

# Evidently, the result here is closer to the true value. This is to be
# expected. As the number of samples tends to infinity, the variance of the
# estimate tends to zero, and its value to the true value.



################################################################################
# 1.3  Scallop data: exploration

# Scallops data set information on 127 locations at which scallops
# were collected in a 1990 survey cruise in the Atlantic continental shelf
# off Long Island, New York, USA.
################################################################################

################################################################################

# Read in the scallops data.
scallops <- read.table("C:/Users/HP/Desktop/scallops.dat", header=TRUE)
scallops

#  scallops
#     lat      long     tcatch       y
#1   40.38333 -71.85000      1 0.00000
#2   40.13333 -72.08333      2 0.69315
#3   40.10000 -72.31667      7 1.94591
#4   40.01667 -72.40000     13 2.56495
#....
#....

# The continuous inputs: the  Longitude (long) and latitude (lat) of each location was recorded. 
# The output (tcatch is discrete) and represents the total number of scallops caught at the ith sample location


# Compute the sample mean and display it.
m <- colMeans(scallops[,c("long","lat")])

m
# long       lat
# -72.73215  39.91798

# Compute the sample variance and display it.
Sigma<- var(scallops[,c("long","lat")])
Sigma
#           long       lat
# long 0.2636210 0.2531269
# lat  0.2531269 0.3699811

# Compute the eigenvalues and eigenvectors.
eVs = eigen(Sigma)

eVs$values
# [1] 0.57545402 0.05814811
# Note that all eigenvalues are > 0, meaning Sigma is positive definite.

################################################################################
# Make a scatterplot of the longitude and latitude inputs. 


# Plot them.
plot(scallops$long, scallops$lat)

# The mean can then be added to the scatterplot via
points(m[1],m[2], col=2, pch="+")



################################################################################

# Plot histograms .
hist(scallops$tcatch)
hist(log(scallops$tcatch))

# We see that the distribution of log-abundances is far less skewed, and appears
# to be more 'normal'. Normality of the response is a standard assumption of
# regression models (which we will consider later).
