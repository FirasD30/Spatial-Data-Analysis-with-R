#------------------------------------------------------------------
# Load spatial data and fit a GLM and an autocovariate logistic
# regression
#
# Partly derived and modified from:
# Fletcher & Fortin
# Spatial Ecology and Conservation Monitoring: Applications with R
# Springer (2018)
#------------------------------------------------------------------

library(raster)
library(spdep)
library(vegan)
library(sf)

setwd("c:/users/derek/work/teaching/biol 4062/2023")
source("icorrelogram.R")

# Load the raster
elev<-raster("c:/users/derek/work/teaching/biol 4062/example data/elev.gri")

# Set the projection: Albers conic equal area
elev.crs <- CRS("+proj=aea +lat_1=46 +lat_2=48 +lat_0=44 +lon_0=-109.5
              +x_0=600000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
proj4string(elev) <- elev.crs

# Load the presence/absence point data
# Data are point counts of 100m radius from
# transects of about 3km long of the Varied Thrush
# (Ixoreus naevius) in N. Idaho and W. Montana.
# Data originally from Northern Region Landbird
# Monitoring Program: see Hutto & Young 2002
VariedThrush <- read.csv("c:/users/derek/work/teaching/biol 4062/example data/VariedThrush.csv", header=TRUE)
head(VariedThrush)

# Plot the elevation and the presences/absences
plot(elev)
points(VariedThrush[VariedThrush$VATH==1, c("EASTING","NORTHING")], col="red")      
points(VariedThrush[VariedThrush$VATH==0, c("EASTING","NORTHING")], col="black")

# Extract coordinate data at sampling points
coords <- cbind(VariedThrush$EASTING, VariedThrush$NORTHING)

# Combine bird observations and elevations
VariedThrush <- cbind(VariedThrush, elev = extract(x=elev, y=coords))

# Look at this combined data frame
head(VariedThrush)

# Plot a correlogram of the raw data
VATH.cor <- icorrelogram(locations=coords, z=VariedThrush$VATH, binsize=1000, maxdist=15000)
plot(VATH.cor$dist, VATH.cor$Morans.i, ylim = c(-0.5, 0.5), pch = 19, xlab = "Distance, m", ylab = "Autocorrelation", main = "Raw data")
abline(h=0, lty = "dashed")
lines(VATH.cor$dist, VATH.cor$null.lower, col = "blue")
lines(VATH.cor$dist, VATH.cor$null.upper, col = "blue")



#-----------------------------------------------
# Logistic regression of presence/absence vs. 
# elevation, ignoring spatial dependence
#-----------------------------------------------

# Consider centreing and scaling predictors before modelling
#  a proper analysis to prevent v. high correlation between
# elevation and elevation ^2
# VariedThrush$elev <- scale(VariedThrush$elev, center = T, scale = T)

# Logistic regression with elevation as a linear + quadratic term
VariedThrush.glm <- glm(VATH ~ elev + I(elev^2),family="binomial", data=VariedThrush)
summary(VariedThrush.glm)

# Predict on to new data to plot against elevation
NewElevs<- seq(min(VariedThrush$elev), max(VariedThrush$elev), length=50)
newdata <- data.frame(elev=NewElevs)
glm.pred <- predict(VariedThrush.glm, newdata=newdata, type= "link", se=T) #type=response for predicted probabilities
ucl <- glm.pred$fit + 1.96*glm.pred$se.fit
lcl <- glm.pred$fit - 1.96*glm.pred$se.fit

# Back-transform from link to probability scale
glm.newdata <- data.frame(newdata, pred=plogis(glm.pred$fit), lcl=plogis(lcl), ucl=plogis(ucl))


# Plot predictions vs elevation
dev.new()
plot(glm.newdata$elev, glm.newdata$pred,ylim=c(0,0.5), pch = 19, xlab="Elevation", ylab="Prob. Occurrence", main = "GLM")
lines(glm.newdata$elev, glm.newdata$lcl)
lines(glm.newdata$elev, glm.newdata$ucl)

# Map the model projections
layers <- stack(elev)
names(layers) <- c("elev")
class(layers)
glm.raster <- predict(model=VariedThrush.glm, object=layers,  type="response")
plot(glm.raster, xlab = "Longitude", ylab = "Latitude", main = "Probability of presence")
points(VariedThrush[VariedThrush$VATH==1, c("EASTING","NORTHING")], col="red")      

#------------------------------------------------#
#Inspect spatial dependence of GLM residuals

# Residuals from quadratic elevation model
VATH.elev.res <- residuals(VariedThrush.glm, type="deviance")

# Correlogram on rmodel esiduals
corr.res <- icorrelogram(locations=coords, z=VATH.elev.res, binsize=1000, maxdist=15000)

# Plot correlogram
dev.new()
plot(corr.res$dist, corr.res$Morans.i, ylim = c(-0.5, 0.5), pch = 19, xlab = "Distance", ylab = "Morans I", main = "GLM residuals")
abline(h=0, lty = "dashed")
lines(corr.res$dist, corr.res$null.lower, col = "blue")
lines(corr.res$dist, corr.res$null.upper, col = "blue")

#------------------------------------------
# Fit an autocovariate logistic regression
# Integrates surrounding values

# Create autocovariate with 1km radius
auto1km <- autocov_dist(VariedThrush$VATH, coords, nbs = 1000, type= "one",style="B", zero.policy=T)  #binary weights

# Fit autocovariate model
acv.VariedThrush <- glm(VATH~elev+I(elev^2)+auto1km, family="binomial", data=VariedThrush)
# extra variable added

#inspect
summary(acv.VariedThrush)

#residuals
VATH.auto.res <- residuals(acv.VariedThrush, type="deviance")

#correlogram on residuals
cor.auto.res <- icorrelogram(locations=coords, z=VATH.auto.res, binsize=1000, maxdist=15000)

# Plot correlogram
plot(cor.auto.res$dist, cor.auto.res$Morans.i, ylim = c(-0.5, 0.5), pch = 19, xlab = "Distance", ylab = "Morans I", main = "Autocovariate model residuals")
abline(h=0, lty = "dashed")
lines(cor.auto.res$dist, cor.auto.res$null.lower, col = "blue")
lines(cor.auto.res$dist, cor.auto.res$null.upper, col = "blue")
    

