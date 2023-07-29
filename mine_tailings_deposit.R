#### Create modelgrid ####
library( rgdal)
library( raster )
library(ncdf4)

source( "C:/Users/josej/Escritorio/Universidad 2021/Trabajo de titulo/Scripts/Load_datasets_DGA_AGROMET_DMC.R" )
source("C:/Users/josej/Escritorio/Universidad 2021/Trabajo de titulo/Scripts/Load_datasets_DGA_AGROMET_DMC.R")
source("C:/Users/josej/Escritorio/Universidad 2021/Trabajo de titulo/Scripts/Convert_Time.R")
source("C:/Users/josej/Escritorio/Universidad 2021/Trabajo de titulo/Scripts/Hydrological_Functions.R")
source("C:/Users/josej/Escritorio/Universidad 2021/Trabajo de titulo/Scripts2/Geospatial_Analysis.R")
source( "C:\\Users\\josej\\Escritorio\\Universidad 2021\\Trabajo de titulo\\Scripts2\\Parallel_Computing.R" )
source( "C:\\Users\\josej\\Escritorio\\Universidad 2021\\Trabajo de titulo\\Scripts2\\Charts.R" )

##### Water table map ######

file = "C:/Users/josej/JupyterLab/MODFLOW Flopy/Example1/MF6/Ovalle/Quebrada El Ingenio/Shapefiles/QuebradaElIngenio_Extendida_v2.shp"
shp_bas = readOGRv2( file )

#### Create Hydraulic head map ####
# OBS 
file = "C:/Users/josej/JupyterLab/MODFLOW Flopy/Example1/MF6/Ovalle/Quebrada El Ingenio/Shapefiles clipped/ObservationWells.shp"
shp = readOGRv2( file )

file = "C:/Users/josej/JupyterLab/MODFLOW Flopy/Example1/MF6/Ovalle/Quebrada El Ingenio/Hydrology/Well Levels/Observation_wells_2016-2018.txt"
data = read.csv(file, skip=3, check.names = FALSE)
StationNames = colnames(data)[2:dim(data)[2]]
idx = which( StationNames %in% shp@data$Name )
data = data[,c(1,idx + 1)]

# Data important variables
StationNames = colnames(data)[2:dim(data)[2]]
nstations = length( StationNames )
idx = which( shp@data$Name  %in%  StationNames  )
shp = shp[idx,]

# QA/QC
# sd( station ) != 0

# Filter by incorrect measurements
aux = data[, !names( data ) %in% "ASENTAMIENTO LAS VEGAS" ]
data = aux



# Interpolate missing data as monthly average
time = as.POSIXct( data[,1] , format = "%Y-%m-%d", tz="GMT")

if ( FALSE ){
dataInterpolated = ConvertTime(data = data[,2:dim(data)[2]],
                               time = time,
                               format_time = "%Y-%m-%d",
                               interval_time_output = "month",
                               method="mean",
                               na_accept = 31)
}

dataInterpolated = ConvertTime(data = data[,2:dim(data)[2]],
                               time = time,
                               format_time = "%Y-%m-%d",
                               interval_time_output = "2 months",
                               method="mean",
                               na_accept = 61)

dataWOInterpolation = data

aux = dataInterpolated
data = aux

# Data important variables
nstations
time = data[,1]
ntime = dim( data )[1]
data

#data = dataInterpolated
# Add water table to df
order = match(  StationNames, shp@data$Name   )
ntime = length( data$Time )
for ( itime in 1:ntime){
  newCol = data.frame( t( data[itime, 2:dim(data)[2]] ) )
  #newCol = setNames( newCol, paste0('Time',data$Time[itime]) )
  newCol = setNames( newCol, paste0('Time',sprintf( itime, fmt = '%02d') ) )
  
  shp@data = cbind( shp@data, newCol )
}

# Plot time series
library(RColorBrewer)
# Hexadecimal color specification 
mycol = brewer.pal(n = 8, name = "Set1")

# Time series with values of 0 could mean that water table is over terrain height 
# (i.e. is a drain )
nstations = dim( data )[2] -1
for ( istation in 1:nstations){
  plot(dataInterpolated[,istation + 1], ylim=c(0,500), col = mycol[istation])
  par( new =TRUE)
}


# Grid for interpolation
grid = rgeos:: gCentroid(BBIdomain,byid=TRUE)
plot( grid, pch = 16 )


library(gstat)
library(ggplot2)

#### Kriging ####
# Fitting a variogram
lzn.vgm <- variogram( Time01 ~ 1, shp) # calculates sample variogram values
# How you can generate variograms using colnames with numbers (i.e. dates)
#lzn.vgm <- variogram( x + y ~ 2016-01-01, shp) # calculates sample variogram values 

#### Kriging with trend ####
lm <- function (formula, data, subset, weights, na.action,
                method = "qr", model = TRUE, contrasts = NULL, offset, ...)
{
  cl <- match.call()
  
  ## keep only the arguments which should go into the model frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- quote(stats::model.frame) # was as.name("model.frame"), but
  ##    need "stats:: ..." for non-standard evaluation
  mf <- eval.parent(mf)
  if (method == "model.frame") return(mf)
  
  ## 1) allow model.frame to update the terms object before saving it.
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  
  ## 2) retrieve the weights and offset from the model frame so
  ## they can be functions of columns in arg data.
  w <- model.weights(mf)
  offset <- model.offset(mf)
  x <- model.matrix(mt, mf, contrasts)
  ## if any subsetting is done, retrieve the "contrasts" attribute here.
  
  z <- lm.fit(x, y, offset = offset, ...)
  class(z) <- c(if(is.matrix(y)) "mlm", "lm")
  
  ## 3) return the na.action info
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  
  ## 4) return the contrasts used in fitting: possibly as saved earlier.
  z$contrasts <- attr(x, "contrasts")
  
  ## 5) return the levelsets for factors in the formula
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  if (model)  z$model <- mf
  z
}

# Prepare data
library( reshape2 )
library( dplyr )

time = 1:length( dataWOInterpolation$Time ) 
time = data.frame( "time" = dataWOInterpolation$Time, "t" =time )

data = melt( dataWOInterpolation, id="Time" )
data$x = NA
data$x =  shp$x[ match( data$variable, shp$Name)]
data$y =  shp$y[ match( data$variable, shp$Name)]
data$t =  time$t[ match( data$Time, time$time)]


#data = as_tibble( data )
#data2 = data %>% select( "x","y","t","value")
idx = match( c( "x","y","t","value"), names(data) )
idx = idx[ !is.na(idx) ]
data2 = data[,  idx ]
weights = rep(1, dim( data2 )[1] )
na.action = na.omit
offset = NULL
dataMean = mean( data2$value, na.rm=T )

model_lm = lm( value ~ x + y + t , data2, weights = weights, na.action = na.action,
              method = "qr", model = TRUE, contrasts = NULL, offset = offset)
#trend = 2
aux = predict( model_lm, data2, interval = "confidence", se.fit = T, level = 0.682) # level = 1 sd

data2$lm_fit = aux$fit[,1]
data2$lm_fit_lb = aux$fit[,2]
data2$lm_fit_up = aux$fit[,3]
data2$lm_res = data2$value - data2$lm_fit


model_krig <- krige( 2016-01-01 ~ 1, shp, meuse.grid, model=lzn.fit)

# Creation of grid for new sampling

# Zone array
ext = extent( shp_bas )

dist = 1000
ncols = ceiling( (ext[2] - ext[1])/ dist )
res = ( ext[2] - ext[1] )/ncols
nrows = ceiling( ( ext[4] - ext[3] ) /res )
s
extNew = extent( rbind( c(ext[1], ext[1] + res* ncols), c(ext[3], ext[3] + res* nrows ) ) )
ras  <- raster(ext = extNew,
                  res= res ,
                  crs = crs(shp_bas),
                  vals = 1)
ras =  mask( ras, shp_bas )

pointGrid = rasterToPoints( ras )
zonePoints = SpatialPointsDataFrame( coords = pointGrid, 
                                     proj4string = crs( shp_bas ),
                                     data = as.data.frame( pointGrid ) )
plot( zonePoints )
plot( ras )

shp@data = data2
model_krig <- krige( formula = lm_res ~ x + y + t, data = as(data2, "Spatial"),
                     newdata = zonePoints,
                     na.action = na.pass )

# Try kriging as a function of distance( sqrt(x^2+y^2) ) and time (t)


 # Kriging as function of position and time (x,t) --------------------------

library(gstat)
library(sp)
library(spacetime)
library(raster)
library(rgdal)
library(rgeos) 

# input 
file = "C:/Users/josej/Downloads/ozon_tram1_14102011_14012012(1).csv"

data <- read.table( file, sep=",", header=T)
data$TIME <- as.POSIXlt(as.numeric(substr(paste(data$generation_time), 1, 10)), origin="1970-01-01")
data$LAT <- as.numeric(substr(paste(data$latitude),1,2))+(as.numeric(substr(paste(data$latitude),3,10))/60)
data$LON <- as.numeric(substr(paste(data$longitude),1,1))+(as.numeric(substr(paste(data$longitude),2,10))/60) 

sub <- na.omit(data) 
sub$latitude <-  as.numeric( sub$latitude )
sub$longitude <-  as.numeric( sub$longitude )

sub <- data[data$TIME>=as.POSIXct('2011-12-12 00:00 CET')&data$TIME<=as.POSIXct('2011-12-14 23:00 CET'),]
#Create a SpatialPointsDataFrame
coordinates(sub)=~LON+LAT
projection(sub)=CRS("+init=epsg:4326")

#Transform into Mercator Projection
ozone.UTM <- spTransform(sub,CRS("+init=epsg:3395")) 

#
ozoneSP <- SpatialPoints(ozone.UTM@coords,CRS("+init=epsg:3395")) 

dupl <- zerodist(ozoneSP) 

#
ozoneDF <- data.frame(PPB=ozone.UTM$ozone_ppb[-dupl[,2]]) 
ozoneTM <- as.POSIXct(ozone.UTM$TIME[-dupl[,2]],tz="CET") 
timeDF <- STIDF(ozoneSP,ozoneTM,data=ozoneDF) 
stplot(timeDF) 

var<- variogramST(PPB~1,data=timeDF,tunit="hours",assumeRegular=F,na.omit=T, cores = 4,
                  twindow = 20)
                 
plot( var )
##
## Fit variogram 
{
# input 
file = "C:/Users/josej/Downloads/ozon_tram1_14102011_14012012(1).csv"
  
data <- read.table( file, sep=",", header=T)
paste(data$generation_time[1])
data$TIME <- as.POSIXlt(as.numeric(substr(paste(data$generation_time), 1, 10)), origin="1970-01-01")
as.POSIXlt(as.numeric(substr(paste(data$generation_time[1]), start=1, stop=10)), origin="1970-01-01")
data$LAT <- as.numeric(substr(paste(data$latitude),1,2))+(as.numeric(substr(paste(data$latitude),3,10))/60)
data$LON <- as.numeric(substr(paste(data$longitude),1,1))+(as.numeric(substr(paste(data$longitude),2,10))/60)

View( data)
# data with NAs
dim( data)
# Important step to filter NA values
data <- na.omit(data) 
# data without NAs
dim( data)
sub <- data[data$TIME>=as.POSIXct('2011-12-12 00:00 CET')&data$TIME<=as.POSIXct('2011-12-14 23:00 CET'),]

#Create a SpatialPointsDataFrame
coordinates(sub)=~LON+LAT
projection(sub)=CRS("+init=epsg:4326")

#Transform into Mercator Projection
ozone.UTM <- spTransform(sub,CRS("+init=epsg:3395")) 

ozoneSP <- SpatialPoints(ozone.UTM@coords,CRS("+init=epsg:3395")) 

dupl <- zerodist(ozoneSP) 
ozoneDF <- data.frame(PPB=ozone.UTM$ozone_ppb[-dupl[,2]]) 

ozoneTM <- as.POSIXct(ozone.UTM$TIME[-dupl[,2]],tz="CET") 
timeDF <- STIDF(ozoneSP,ozoneTM,data=ozoneDF) 

var <- variogramST(PPB~1,data=timeDF,tunit="hours",assumeRegular=F,na.omit=T, cores = 4) 
varCorrect = var
plot(var,map=T) 
plot(var,map=F) 
}


# ST Kriging with gstat ---------------------------------------------------
{
  
  # Set data 
  data = dataInterpolated
  
  # Standarization
  {
    meanVar = array( NA, nstations)
    for ( istation in 1:nstations ){
      meanVar[istation] = mean( data[,1+istation], na.rm = T ) 
      
    }
    
    
    # Spatial standarization
    for ( istation in 1:nstations ){
      data[,1+istation]  =   ( data[,1+istation] - mean( meanVar[istation] ) ) /       sd( meanVar, na.rm =T )^2

    }
    
    # Temporal standarization 
    for ( istation in 1:nstations ){
      data[,1+istation] = ( data[,1+istation] - mean( data[,1+istation], na.rm = T ) ) /
        sd( data[,1+istation], na.rm =T )
      
    }
  }
  
  
shp_forVgm =   SpatialPoints(data.frame(x = 0, y = 0))[-1,]
time_forVgm = c()
for (istation in 1:( nstations )){
for (iobs in 1:dim(data)[1]){
  shp_forVgm = bind( shp_forVgm, shp[istation,] )
}
  time_forVgm = cbind( time_forVgm, data[,1])
}


# Convert to correct classes
data_forVgm = data.frame( "value" = unlist( data[,2:nstations] ) )
idx  =  !is.na( data_forVgm )
data_forVgm = subset( data_forVgm, idx)


shp_forVgm = as(shp_forVgm, "SpatialPoints")
shp_forVgm = shp_forVgm[idx]

time_forVgm = as.vector( time_forVgm [idx] )
time_forVgm = as.POSIXct(time_forVgm, origin = '1970-01-01', tz ="GMT")



##
dfIST = STIDF( sp = shp_forVgm,
              time = time_forVgm,#ST( shp, data$Time,  data$Time + days(1) ),
              data = data_forVgm     )

stplot( dfIST )


# Plot
{
cutoff = 25000
width =  1000*1/10
boundaries = seq(0, cutoff, width)

# we are gonna work in days
tunit="days"
twindow = (tail( time, 1 ) -  head( time, 1 ) )
timeWhole = seq( head( time, 1 ), by= tunit, length.out = twindow )
tlags = which(  timeWhole  %in% time  ) - 1
var <- variogramST(value~1,data=dfIST,tunit= tunit, cutoff = cutoff,
                   width = width, boundaries = boundaries,
                   assumeRegular=F,na.omit=T,  twindow = twindow,
                   tlags = tlags) 


plot(var,map=F) 
plot(var,map=T) 

plot(var,wireframe=T) 

}



# Model variogrm fit functions --------------------------------------------
{
  # in gstat we have 5 options: separable, product sum, metric, sum metric, and simple sum metric. 
  # You can find more information to fit these model, including all the equations presented below,
  # in (Gr?ler et al., 2015)
  # lower and upper bounds
  delta = 0.0001
  pars.l <- c(sill.s = 0, range.s = 10, nugget.s = 0 + delta ,
              sill.t = 0, range.t = 1, nugget.t = 0 + delta,
              sill.st = 0, range.st = 10, nugget.st = 0 + delta,
              anis = 0)
  pars.u <- c(sill.s = 1.0, range.s = 1000, nugget.s = 0.4,
              sill.t = 1.0, range.t = 60, nugget.t = 0.4,
              sill.st = 200, range.st = 1000, nugget.st = 100,
              anis = 700) 
  
  ##
  library(rtop)
  
  source( "C:/Users/josej/Escritorio/Universidad 2021/Trabajo de titulo/Scripts/Load_datasets_DGA_AGROMET_DMC.R" )
  source("C:/Users/josej/Escritorio/Universidad 2021/Trabajo de titulo/Scripts/Load_datasets_DGA_AGROMET_DMC.R")
  source("C:/Users/josej/Escritorio/Universidad 2021/Trabajo de titulo/Scripts/Convert_Time.R")
  source("C:/Users/josej/Escritorio/Universidad 2021/Trabajo de titulo/Scripts/Hydrological_Functions.R")
  source("C:/Users/josej/Escritorio/Universidad 2021/Trabajo de titulo/Scripts2/Geospatial_Analysis.R")
  source( "C:\\Users\\josej\\Escritorio\\Universidad 2021\\Trabajo de titulo\\Scripts2\\Parallel_Computing.R" )
  
  
  hydrofunction = function( pars, vgmtype1, vgmtype2, vgmtype3, covtype){
    # Compute daily abstraction with model 1
    # X: Array of dimension as the parameter spatial grid. Daily abstraction in mm.
    # A: Array of dimension as parameter spatial grid. Fraction of area of every cell inside
    # basin.
    
    if ( FALSE ){
      
      vgm_spatial = vgm( pars[1], eval(parse(text=  vgmtype1), envir= .GlobalEnv) ,
                         pars[2], pars[3] )
      vgm_temp = vgm( pars[4], eval(parse(text=  vgmtype2), envir= .GlobalEnv) ,
                      pars[5], pars[6] )
      vgm_joint = vgm( pars[7], eval(parse(text=  vgmtype3), envir= .GlobalEnv) ,
                       pars[8], pars[9] )
    }
    vgm_spatial = vgm( pars[1],  vgmtype1,
                       pars[2], pars[3] )
    vgm_temp = vgm( pars[4],  vgmtype2,
                    pars[5], pars[6] )
    vgm_joint = vgm( pars[7], vgmtype3, 
                     pars[8], pars[9] )
    
    if (covtype == "separable"){
      # 1: Separable
      separable <- vgmST("separable", space = vgm_spatial,time = vgm_temp, sill = vgm_joint$psill[2], stAni=pars[10]) 
      #separable_Vgm <- fit.StVariogram(var, separable,  fit.method = 6 ,method="L-BFGS-B",lower=pars.l,upper=pars.u,tunit="hours")
      attr(separable_Vgm, "MSE")
      Vgm = separable
    } else if (covtype == "sumMetric"){
      
      # 3: Sum metric
      sumMetric <- vgmST("sumMetric", space = vgm_spatial,time = vgm_temp, joint = vgm_joint, stAni=pars[10]) 
      #sumMetric_Vgm <- fit.StVariogram(var, sumMetric,  fit.method = 10 , method="L-BFGS-B",lower=pars.l,upper=pars.u,tunit="hours")
      attr(sumMetric_Vgm, "MSE")
      Vgm = sumMetric
      
    } else if (covtype == "simpleSumMetric"){
      
      # 4: Simple sum metric
      SimplesumMetric <- vgmST("simpleSumMetric",space = vgm_spatial, time = vgm_temp, joint = vgm_joint, nugget=1, stAni=pars[10]) 
      #SimplesumMetric_Vgm <- fit.StVariogram(var, SimplesumMetric,  fit.method = 10 ,method = "L-BFGS-B",lower=pars.l)
      attr(SimplesumMetric_Vgm, "MSE")
      Vgm = SimplesumMetric
    }
    dist_grid = var
    var_model = variogramSurface( Vgm, dist_grid = dist_grid )
    colnames( var_model )[ dim( var_model )[2] ] = "gamma_model"
    sim = var_model$gamma_model
    
    return(sim)
  }
  
  
  MSE = function( sim, obs, mult = 1 ){
    mse = mean( (obs - sim)^2, na.rm = T)*mult
    
    return( mse )
  }
  
  
  ObjFun_MSE = function( pars, obs, mult, vgmtype, covtype ){
    vgmtype1 = vgmtype[1]
    vgmtype2 = vgmtype[2]
    vgmtype3 = vgmtype[3]
    
    sim = hydrofunction( pars, vgmtype1 = vgmtype1, vgmtype2 = vgmtype2, vgmtype3 = vgmtype3,
                         covtype = covtype)
    
    mse = MSE( sim, obs, mult )
    return( mse )
  }
  
  
  
  ##
  models = c("Exp", "Sph", "Gau", "Mat" )
  models = expand.grid( models, models, models, stringsAsFactors = FALSE )
  
  # Covariance models
  ncols = 3
  
  # variograms combinations
  nrows = dim( models )[1]
}



# STFDF
mycrs = crs( shp )
shp_forVgm =   SpatialPoints(data.frame(x = 0, y = 0), proj4string = mycrs)[-1,]
time_forVgm = c()
for (iobs in 1:dim(data)[1]){
  for (istation in 1:( nstations )){
    shp_forVgm = rbind( shp_forVgm, shp[istation,] )
    
  }
  time_forVgm = c( time_forVgm, data[iobs,1])
}

data_forVgm = data.frame( "value" = unlist( data[,2:8]  ) )
#idx  =  !is.na( data_forVgm ) 
#data_forVgm = subset( data_forVgm, idx)

shp_forVgm = as(shp_forVgm, "SpatialPoints")
#shp_forVgm = shp_forVgm[idx]

time_forVgm = as.POSIXct(time_forVgm, origin = '1970-01-01', tz ="GMT")
#time_forVgm = time_forVgm[idx]

#dfST = STFDF(sp  = shp_forVgm, time = time_forVgm, data  = data_forVgm  , endTime = delta(time))
}




##
{
mycrs = crs( shp )
shp_forVgm =   SpatialPoints(data.frame(x = 0, y = 0), proj4string = mycrs)[-1,]
shp_forVgm = rbind( shp_forVgm, shp )
time_forVgm = data[,1]
time_forVgm = as.POSIXct(time_forVgm, format= "%Y-%m-%d", tz ="GMT")


data_forVgm = data.frame( "value" = unlist( data[,2:8]  ) )
shp_forVgm = as(shp, "SpatialPoints")
#shp_forVgm = shp_forVgm[idx]

time_forVgm = as.POSIXct(time_forVgm, origin = '1970-01-01', tz ="GMT")
#time_forVgm = time_forVgm[idx]


data_forVgm = data.frame( "value" = unlist( data[,2:8]  ) )
##
dfST = STIDF( sp = shp_forVgm,
              time = time_forVgm,#ST( shp, data$Time,  data$Time + days(1) ),
        endTime = time_forVgm + days( 1 ),
          data = data_forVgm )


sp = shp_forVgm
time = time_forVgm
endTime = time_forVgm + days( 1 )
data = data_forVgm



gstat:::new

is.nan( time_forVgm )
?ST
stplot(dfST) 
View( dfST )
View( var )

View( data_forVgm )

# histogram of distances
distMatrix = distSpatialPoints( shp )
distMatrix[ !upper.tri( distMatrix ) ]= NA

arr  = na.omit( as.vector(distMatrix) )
arr  = distMatrix[ !is.na( distMatrix )]

h <- hist( arr, breaks = pretty(arr,5), freq = FALSE, col = "gray", ylim = c(0, 1))
h <- hist_def( arr, breaks = pretty(arr,5), freq = FALSE, col = "gray", ylim = c(0, 1))


# Experimental variogram
{
cutoff = 25000
twindow = 30*24
width =  100
boundaries = seq(0, cutoff, width)
var <- variogramST(value~1,data=dfST,tunit="days", cutoff = cutoff,
                   width = width, boundaries = boundaries,
                   assumeRegular=F,na.omit=T,  twindow = twindow) 
stplot( var )


plot(var,map=F) 
plot(var,map=T) 

plot(var,wireframe=T) 
}

}

# Fit model variogram 





# covariograms can be used for spatial (continuous  variable) and temporal ( normally used as
# discrete variable). Auto correlation function also can be used for time and also! spatial
# continuous variable and for 
# 
{
  library( reshape )
  library( rgeos )
  
  plot_directory = "C:/Users/josej/JupyterLab/MT3D USGS/Examples/Mine tailings deposit/images"
  file = "C:/Users/josej/JupyterLab/MODFLOW Flopy/Example1/MF6/Ovalle/Quebrada El Ingenio/model2/model2_modelgrid.shp"
  dsn = "C:/Users/josej/JupyterLab/MODFLOW Flopy/Example1/MF6/Ovalle/Quebrada El Ingenio/model2"
  layer = "model2_modelgridIdomain"
  BBIdomain = readOGR(dsn = dsn, layer = layer)
  
  
  aux = apply( BBIdomain@data, 2, as.numeric)
  for (icol in 1:dim(BBIdomain@data)[2]){
    BBIdomain@data[,icol] = aux[,icol]
  }
  
  plot(BBIdomain)
  
  centroid = gCentroid(BBIdomain, byid=TRUE)
  BBIdomain$x = centroid@coords[,1]
  BBIdomain$y = centroid@coords[,2]
  BBIdomain = BBIdomain[ BBIdomain$idomain_1 == 1, ]
  
  

  # Water table level: variable z
  dfVar = data
  
  # Standarization
  {
    meanVar = array( NA, nstations)
    for ( istation in 1:nstations ){
      meanVar[istation] = mean( dfVar[,1+istation], na.rm = T ) 
      
    }
    
    
    # Spatial standarization
    for ( istation in 1:nstations ){
      dfVar[,1+istation]  =   ( dfVar[,1+istation] - mean( meanVar[istation] ) ) /       sd( meanVar, na.rm =T )^2
      
    }
    
    # Temporal standarization 
    for ( istation in 1:nstations ){
      dfVar[,1+istation] = ( dfVar[,1+istation] - mean( dfVar[,1+istation], na.rm = T ) ) /
        sd( dfVar[,1+istation], na.rm =T )
      
    }
  }
  
  
  
  dfVar$Time = 1:length( data$Time )
  dfVar = melt(  dfVar , id.vars = "Time" )
  x = array( NA, dim( dfVar)[1] )
  y = array( NA, dim( dfVar)[1] )
  
  for ( irow in 1:dim( dfVar )[1] ){
    idx = match( dfVar$variable[irow] ,shp$Name )
    x[ irow ] = shp$x[ idx ]
    y[ irow ] = shp$y[ idx ]
    
  }
  dfVar$x = x
  dfVar$y = y
  
  # Load rch bivariate variable: RCH ( u, t)
  file = "C:/Users/josej/JupyterLab/MODFLOW Flopy/Example1/MF6/Ovalle/Quebrada El Ingenio/model/MODFLOW Variables/RDaily.txt"
  RCH = read.table( file, skip = 3, sep = ",", header = T, check.names = F, stringsAsFactors = F)
  #dfRCH = data.frame( "cell" = 1:(dim(RCH)[2]-1), )
  dfRCH = as.data.frame( RCH )
  colnames( dfRCH )[1] = "Time"
  #dfRCH = melt( dfRCH, id.vars = "Time",   factorsAsStrings = T )
  n =  dim( dfVar)[1] / dim( dfRCH )[1]
  d = dfRCH
  aux = replicate(n, d, simplify = FALSE)
  dfRCH = do.call("rbind", aux ) 
  names( dfRCH) = paste0( "RCH t", 1:dim( dfRCH)[2])
  dfVar = cbind( dfVar, dfRCH)

  # Exploratory Analysis
  arr = dfVar$value[ !is.na( dfVar$value ) ]
  arr_t = arr - min( arr ) + 1
  lambda = seq(-2, 2, 1/10)
  BCtrans = boxcox(lm( arr_t ~ 1) , lambda = lambda)
  lambda = BCtrans$x[ which.max( BCtrans$y ) ]
  arr <- (arr_t ^ lambda - 1) / lambda
  
  main = paste("Histograma de" , "z_t")
  ylab = "Frecuencia relativa"
  xlab = "z_t"
  h = hist( arr, breaks = pretty(arr,5), freq = F )
  h$rel_counts <- h$counts/sum(h$counts)
  h$counts = h$rel_counts
  h$density = h$rel_counts
  plot( h, col = "gray", ylim = c(0, 1),
        main = main, xlab = xlab, ylab = ylab )
  plot( ecdf( arr ), verticals = T, add = T, pch = NA)
  #h <- hist_def( arr, breaks = pretty(arr,5), freq = FALSE, col = "gray", ylim = c(0, 1),
  #               main = main, xlab = xlab, ylab = ylab )
  ks.test(arr, "pnorm", mean( arr, na.rm=T ), sd( arr, na.rm=T ))
  shapiro.test( ( arr - mean( arr, na.rm=T ) )/ ( sd( arr, na.rm=T )^2 ) )
  
  p1 = recordPlot()
  filename = paste( plot_directory,paste0("ExploratoryAnalysis_",xlab,".png"), sep="/")
  png_base(p1,filename,480,480,overwrite=T)
  
  arr_new =  dfVar$value
  arr_new[ !is.na( arr_new ) ] = arr
  #arr_new[ is.na( arr_new ) ] = 0
  
  dfVar$z_t = arr_new

  #dfVar$z_t_decision = ifelse(  dfVar$z_t != 0 , 1, 0 )
  ###
  arr = dfVar$x
  main = main = paste("Histograma de" , "x")
  ylab = "Frecuencia relativa"
  xlab = "x"
  h = hist( arr, breaks = pretty(arr,5), freq = F )
  h$rel_counts <- h$counts/sum(h$counts)
  h$counts = h$rel_counts
  h$density = h$rel_counts
  plot( h, col = "gray", ylim = c(0, 1),
        main = main, xlab = xlab, ylab = ylab )
  plot( ecdf( arr ), verticals = T, add = T, pch = NA)
  
  p1 = recordPlot()
  filename = paste( plot_directory,paste0("ExploratoryAnalysis_",xlab,".png"), sep="/")
  png_base(p1,filename,480,480,overwrite=T)
  
  arr = dfVar$y
  main = main = paste("Histograma de" , "y")
  ylab = "Frecuencia relativa"
  xlab = "y"
  h = hist( arr, breaks = pretty(arr,5), freq = F )
  h$rel_counts <- h$counts/sum(h$counts)
  h$counts = h$rel_counts
  h$density = h$rel_counts
  plot( h, col = "gray", ylim = c(0, 1),
        main = main, xlab = xlab, ylab = ylab )
  plot( ecdf( arr ), verticals = T, add = T, pch = NA)
  
  p1 = recordPlot()
  filename = paste( plot_directory,paste0("ExploratoryAnalysis_",xlab,".png"), sep="/")
  png_base(p1,filename,480,480,overwrite=T)
  
  library(MASS)
  arr =  dfRCH$value[ dfRCH$value > 0 ] 
  arr =  arr[ !is.na ( arr )  ] 
  arr_t = arr - min( arr, na.rm = T ) + 1
  lambda = seq(-2, 2, 1/10)
  BCtrans = boxcox(lm( arr_t ~ 1) , lambda = lambda)
  lambda = BCtrans$x[ which.max( BCtrans$y ) ]
  arr <- (arr_t ^ lambda - 1) / lambda
  
  main = main = paste("Histograma de" , "RCH_t")
  ylab = "Frecuencia relativa"
  xlab = "RCH_t"
  h = hist( arr, breaks = pretty(arr,5), freq = F )
  h$rel_counts <- h$counts/sum(h$counts)
  h$counts = h$rel_counts
  h$density = h$rel_counts
  plot( h, col = "gray", ylim = c(0, 1),
        main = main, xlab = xlab, ylab = ylab )
  plot( ecdf( arr ), verticals = T, add = T, pch = NA)
  ks.test(arr, "pnorm", mean( arr, na.rm=T ), sd( arr, na.rm=T ))
  
  p1 = recordPlot()
  filename = paste( plot_directory,paste0("ExploratoryAnalysis_",xlab,".png"), sep="/")
  png_base(p1,filename,480,480,overwrite=T)
  
  arr_new =  dfRCH$value
  arr_new[ arr_new > 0 &  !is.na( arr_new ) ] = arr
  dfVar$RCH_t = arr_new

  dfVar$RCH_t_decision = ifelse(  dfVar$RCH_t != 0 , 1, 0 )  
  ##
  
  # Add variables
  dfVar$cos_dn = cos( dfVar$Time *2*pi/365 )
  
  # Correlation between variables
  # A) Homogeneous spatial variables
  library( corrplot )
  library(dplyr)
  
  #vars = names( dfVar )
  vars = c("z_t", "x", "y", "RCH_t", "cos_dn")
  
  corvar =   dfVar %>%
    select( vars )

  corvar = corvar[ rowSums(is.na( corvar )) == 0,]
  View( corvar )  
  #corvar = as.matrix( corvar[!is.na( corvar ),] )
  #View ( corvar )  
  corvar2 = cor( corvar , use = "pairwise.complete.obs", method =  "pearson")
  
corrplot( as.matrix( corvar2 ), type = "upper", 
          is.corr = FALSE,
           method = "square", 
           addCoef.col = "black", 
           tl.pos = "lt",
           tl.col = "black")
  
  # B) Satial variables
  
  # variogram
  vars = c("z_t", "x", "y", "RCH_t", "RCH_t_decision", "cos_dn")
  
  corvar =   dfVar %>%
    select( vars )

  coordinates( corvar  ) = ~ x + y 
  #varX = variogram( z_t ~ x + y , locations = x ,data =  dfVar )
  
  alpha = 0
  cutoff = sqrt( ( extent( BBIdomain )[2] - extent( BBIdomain )[2] ) ^2 + 
                   ( extent( BBIdomain )[3] - extent( BBIdomain )[4] ) ^2 )
  width = cutoff/15
  varX = variogram( z_t ~  RCH_t * RCH_t_decision, locations = ~ x + y,  data = corvar )
  
  # Variogram doesn't work since there are little observations
}


  # MLR
 {
 
  library( reshape )
  library( rgeos )
  
  plot_directory = "C:/Users/josej/JupyterLab/MT3D USGS/Examples/Mine tailings deposit/images"
  file = "C:/Users/josej/JupyterLab/MODFLOW Flopy/Example1/MF6/Ovalle/Quebrada El Ingenio/model2/model2_modelgrid.shp"
  dsn = "C:/Users/josej/JupyterLab/MODFLOW Flopy/Example1/MF6/Ovalle/Quebrada El Ingenio/model2"
  layer = "model2_modelgrid"
  BBIdomain = readOGR(dsn = dsn, layer = layer)
  
  
  aux = apply( BBIdomain@data, 2, as.numeric)
  for (icol in 1:dim(BBIdomain@data)[2]){
    BBIdomain@data[,icol] = aux[,icol]
  }
  
  plot(BBIdomain)
  
  centroid = gCentroid(BBIdomain, byid=TRUE)
  BBIdomain$x = centroid@coords[,1]
  BBIdomain$y = centroid@coords[,2]
  BBIdomain = BBIdomain[ BBIdomain$idomain_1 == 1, ]
  
  nrow = max( BBIdomain$row )
  ncol = max( BBIdomain$column )
  
  # Water table level: variable z
  dfVar = data
  flag_value = -999
  
  # Standarization
  {
    meanVar = array( NA, nstations)
    for ( istation in 1:nstations ){
      meanVar[istation] = mean( dfVar[,1+istation], na.rm = T ) 
      
    }
    
    
    # Spatial standarization
    for ( istation in 1:nstations ){
      dfVar[,1+istation]  =   ( dfVar[,1+istation] - mean( meanVar[istation] ) ) /       sd( meanVar, na.rm =T )^2
      
    }
    
    # Temporal standarization 
    for ( istation in 1:nstations ){
      dfVar[,1+istation] = ( dfVar[,1+istation] - mean( dfVar[,1+istation], na.rm = T ) ) /
        sd( dfVar[,1+istation], na.rm =T )
      
    }
  }
  
  
  View( dfVar )
  # standarization
  dfVar$Time = 1:length( data$Time )
  dfVar = melt(  dfVar , id.vars = "Time" )
  
  
  x = array( NA, dim( dfVar)[1] )
  y = array( NA, dim( dfVar)[1] )
  
  for ( irow in 1:dim( dfVar )[1] ){
    idx = match( dfVar$variable[irow] ,shp$Name )
    x[ irow ] = shp$x[ idx ]
    y[ irow ] = shp$y[ idx ]
    
  }
  dfVar$x = x
  dfVar$y = y
  
  # Load rch bivariate variable: RCH ( u, t)
  file = "C:/Users/josej/JupyterLab/MODFLOW Flopy/Example1/MF6/Ovalle/Quebrada El Ingenio/model/MODFLOW Variables/RDaily.txt"
  RCH = read.table( file, skip = 3, sep = ",", header = T, check.names = F, stringsAsFactors = F)
  #dfRCH = data.frame( "cell" = 1:(dim(RCH)[2]-1), )
  dfRCH = as.data.frame( array(NA, c( dim( RCH )[1], prod( c( nrow, ncol) ) ) ) )
  for ( icell in 2:dim( RCH )[2] ){
    cellid = names(RCH)[ icell ]
    cellid = strsplit( cellid, " " )[[1]][4]
    cellid = as.numeric( cellid )
    
    idx_cell = which( BBIdomain$node == cellid )
    jrow = BBIdomain$row[ idx_cell ]
    jcol = BBIdomain$column[ idx_cell ]
    
    dfRCH[ , cellid ] = dfRCH[, icell]
  }
  colnames( dfRCH )[1] = "Time"
  ntime = dim( dfRCH )[1]
  
  #debugonce( fun = ConvertTime )
  
  ## Support for bimonth support
  time = as.POSIXct( "2016-01-01", format="%Y-%m-%d", tz="GMT" ) + days( dfRCH[,1] ) - days(1)
  dfRCH = ConvertTime(data = dfRCH[,2:dim(dfRCH)[2]],
                                   time = time,
                                   format_time = "%Y-%m-%d",
                                   interval_time_output = "2 months",
                                   method="sum",
                                   na_accept = 61)
  ##
  
  n =  dim( dfVar )[1] / dim( dfRCH )[1]
  d = dfRCH
  dfRCH = do.call("rbind", replicate(n, d, simplify = FALSE))
  #names( dfRCH) = paste0( "RCH_t_", 1:dim( dfRCH)[2])
  dfVar = cbind( dfVar, dfRCH[,2:dim(dfRCH)[2]])
  
  # Exploratory Analysis
  # version v2: Support for bimonthly
  version = "v2"
  
  library(MASS)
  arr = dfVar$value[ !is.na( dfVar$value ) ]
  arr_t = arr - min( arr, na.rm=T ) + 1
  lambda = seq(-40, 40, 1/10)
  BCtrans = boxcox(lm( arr_t ~ 1) , lambda = lambda)
  lambda = BCtrans$x[ which.max( BCtrans$y ) ]
  arr <- (arr_t ^ lambda - 1) / lambda
  
  main = paste("Histograma de" , "z_t")
  ylab = "Frecuencia relativa"
  xlab = "z_t"
  h = hist( arr, breaks = pretty(arr,5), freq = F )
  h$rel_counts <- h$counts/sum(h$counts)
  h$counts = h$rel_counts
  h$density = h$rel_counts
  plot( h, col = "gray", ylim = c(0, 1),
        main = main, xlab = xlab, ylab = ylab )
  plot( ecdf( arr ), verticals = T, add = T, pch = NA)
  #h <- hist_def( arr, breaks = pretty(arr,5), freq = FALSE, col = "gray", ylim = c(0, 1),
  #               main = main, xlab = xlab, ylab = ylab )
  ks.test(arr, "pnorm", mean( arr, na.rm=T ), sd( arr, na.rm=T ))
  shapiro.test( ( arr - mean( arr, na.rm=T ) )/ ( sd( arr, na.rm=T )^2 ) )
  
  p1 = recordPlot()
  filename = paste( plot_directory,paste0("ExploratoryAnalysis_",version,"_",xlab,".png"), sep="/")
  png_base(p1,filename,480,480,overwrite=T)
  
  arr_new =  dfVar$value
  arr_new[ !is.na( arr_new ) ] = arr
  #arr_new[ is.na( arr_new ) ] = 0
  
  dfVar$z_t = arr_new
  idx = which( names(dfVar)  == "z_t" )
  dfVar = dfVar[, c(1:(5), idx, 6:( dim( dfVar )[2] ) ) ]
  
  View( dfVar )
  View( arr )
  

  #dfVar$z_t_decision = ifelse(  dfVar$z_t != 0 , 1, 0 )
  ###
  arr = dfVar$x
  main = main = paste("Histograma de" , "x")
  ylab = "Frecuencia relativa"
  xlab = "x"
  h = hist( arr, breaks = pretty(arr,5), freq = F )
  h$rel_counts <- h$counts/sum(h$counts)
  h$counts = h$rel_counts
  h$density = h$rel_counts
  plot( h, col = "gray", ylim = c(0, 1),
        main = main, xlab = xlab, ylab = ylab )
  plot( ecdf( arr ), verticals = T, add = T, pch = NA)
  
  p1 = recordPlot()
  filename = paste( plot_directory,paste0("ExploratoryAnalysis_",version,"_",xlab,".png"), sep="/")
  png_base(p1,filename,480,480,overwrite=T)
  
  arr = dfVar$y
  main = main = paste("Histograma de" , "y")
  ylab = "Frecuencia relativa"
  xlab = "y"
  h = hist( arr, breaks = pretty(arr,5), freq = F )
  h$rel_counts <- h$counts/sum(h$counts)
  h$counts = h$rel_counts
  h$density = h$rel_counts
  plot( h, col = "gray", ylim = c(0, 1),
        main = main, xlab = xlab, ylab = ylab )
  plot( ecdf( arr ), verticals = T, add = T, pch = NA)
  
  p1 = recordPlot()
  filename = paste( plot_directory,paste0("ExploratoryAnalysis_",version,"_",xlab,".png"), sep="/")
  png_base(p1,filename,480,480,overwrite=T)
  
  #for (icol in 6:( 6 + ntime  ) ) {
  idxCol =  7:( 7 + dim(dfRCH)[2] - 1  )
  arr =  dfVar[ dfVar[ , idxCol] > 0, idxCol ] 
  arr =  arr[ arr > 0 &  !is.na( arr )  ] 
  arr_t = arr - min( arr, na.rm = T ) + 1
  lambda = seq(-20, 20, 1/10)
  BCtrans = boxcox(lm( arr_t ~ 1) , lambda = lambda)
  lambda = BCtrans$x[ which.max( BCtrans$y ) ]
  arr <- (arr_t ^ lambda - 1) / lambda
  arr <- array( arr, dim( dfVar[, idxCol] ) )
  #}
  
  main = paste("Histograma de" , "RCH_t")
  ylab = "Frecuencia relativa"
  xlab = "RCH_t"
  h = hist( arr, breaks = pretty(arr,5), freq = F )
  h$rel_counts <- h$counts/sum(h$counts)
  h$counts = h$rel_counts
  h$density = h$rel_counts
  plot( h, col = "gray", ylim = c(0, 1),
        main = main, xlab = xlab, ylab = ylab )
  plot( ecdf( arr ), verticals = T, add = T, pch = NA)
  ks.test(arr, "pnorm", mean( arr, na.rm=T ), sd( arr, na.rm=T ))
  
  p1 = recordPlot()
  filename = paste( plot_directory,paste0("ExploratoryAnalysis_",version,"_",xlab,".png"), sep="/")
  png_base(p1,filename,480,480,overwrite=T)
  
  for ( icol in 2:( dim(dfRCH)[2]  ) ){
  #idxCol = 7:( 7 + ntime  )
  arr_new =  dfRCH[ , icol ]
  #for icol i
  idx = which(  arr_new> 0 &  !is.na( arr_new ) )#, arr.ind =  TRUE  )
  arr_new[ idx ] = arr[ , icol ]
  idx = which(  arr_new == 0 )#, arr.ind =  TRUE  )
  arr_new[ idx ] = NA
  
  idx_dfVar = icol + 7 
  dfVar[ , idx_dfVar ] = arr_new
  
  newName = paste0( "RCH_dec_t_", icol )
  idx_newName = 739 + icol
  dfVar$newCol = ifelse(  dfRCH[ ,icol ] != 0 , 1, 0 )
  
  names(dfVar)[names(dfVar) == "newCol"] <- newName
  #itime = icol - 7 + 1
  #names( dfVar )[idx_newName] = eval( as.name( newName ) )
  }
  ##
  
  # Add variables
  dfVar$cos_dn = cos( dfVar$Time *2*pi/365 )
  
  # Correlation between variables
  # B) spatial variables
  library( corrplot )
  library(dplyr)
  

  corvar =   dfVar %>%
    select( -c("Time", "variable" ) )
  
  corvar =   dfVar %>%
    select( c(3,4,5,6, 7:(7+ntime)) )
  
  #corvar = corvar[ rowSums(is.na( corvar )) == 0,]
  View( corvar )  

  
  corMatrix = cor( corvar , use = "pairwise.complete.obs", method =  "pearson")
  corUnique = corMatrix
  corUnique[ !upper.tri( corMatrix, diag = TRUE ) ] = NA
  idx = unique( which( corUnique > 0.4 , arr.ind = T )[,1] )
  aux = c( 4, idx)
  idx = aux
  
  View( corvar[,idx] )
  #model_mlr =lm( z_t ~ .,  data = corvar[,idx] )

  
  # PCA
  
  # Load data
  # Water table level: variable z
  dfVar = data
  
  aux = data.frame( "time" = dfVar$Time, array(NA, dim= c( length(dfVar$Time),
                                                           dim( BBIdomain )[1] )) )
  
  
  wt_measurements = intersectFeatures( x= BBIdomain, y = shp, arr.ind = T)
  for ( iwell in 1:length( wt_measurements )){
    idx = wt_measurements[ iwell ]
    idx_dfVar = which( names( dfVar ) == shp$Name[ iwell ] )
    aux[, idx + 1 ] = dfVar[ , 1 + iwell]
  }

  dfVar = aux
  
  # Standarization
  
  # Spatial standarization
  ncells = dim( dfVar)[2] - 1
  meanVar = array( NA, ncells  )
  for ( istation in 1:ncells ){
    
    meanVar[istation] = mean( dfVar[,1+istation], na.rm = T )
    dfVar[,1+istation] = ( dfVar[,1+istation] - meanVar[istation] ) /
      sd( dfVar[,1+istation], na.rm =T )
  }
  
  # Temporal standarization 
  for ( istation in 1:ncells ){
    
    
    dfVar[,1+istation]  =   ( dfVar[,1+istation] - mean( meanVar, na.rm = T  ) ) / 
      sd( meanVar, na.rm =T )
    
    #dfVar[,1+istation] = ( dfVar[,1+istation] - mean( dfVar[,1+istation], na.rm = T ) ) /
    #  sd( dfVar[,1+istation], na.rm =T )^2
    
  }
  
  
  # Load rch bivariate variable: RCH ( u, t)
  file = "C:/Users/josej/JupyterLab/MODFLOW Flopy/Example1/MF6/Ovalle/Quebrada El Ingenio/model/MODFLOW Variables/RDaily.txt"
  RCH = read.table( file, skip = 3, sep = ",", header = T, check.names = F, stringsAsFactors = F)
  #dfRCH = data.frame( "cell" = 1:(dim(RCH)[2]-1), )
  dfRCH = as.data.frame( RCH )
  colnames( dfRCH )[1] = "Time"
  ntime = dim( dfRCH )[1]
  
  ## Support for bimonth support
  time = as.POSIXct( "2016-01-01", format="%Y-%m-%d", tz="GMT" ) + days( dfRCH[,1] ) - days(1)
  dfRCH = ConvertTime(data = dfRCH[,2:dim(dfRCH)[2]],
                      time = time,
                      format_time = "%Y-%m-%d",
                      interval_time_output = "2 months",
                      method="sum",
                      na_accept = 61)
  ##
  
  #sstascale  <- scale(SSTaavg)    # Scale each column (i.e., x' = (x - mean(x))/sd(x))
  #zs         <- var(sstascale)    # Compute the variance covariance matrix.
  # Note that it is symetric!
  
  X = as.matrix( dfRCH[, 2:dim( dfRCH )[2]] )
  Y = as.matrix( dfVar[, colSums(is.na(dfVar))<nrow(dfVar)]   )
  
  #View( X )
  #do an Eigen decomposition..
  xs         <- var( X, use = "pairwise.complete.obs")#"all.obs")#"na.or.complete")
  xsvd       <- svd( xs )
  
  #Principal Components...
  SST.PCs        <- X %*%  xsvd$u # Note that zsvd$u is equivalent to matrix E
  #Eigen Values.. - fraction variance 
  SST.lambdas    <- (xsvd$d/sum(xsvd$d))
  
  ############## Plot the eigen spectrum for SST anomalies
  # NpcsSpecSST cant not be higher than the number of observations sites.
  NpcsSpecSST = 4
  first_modes <- NpcsSpecSST
  
  plot(1:first_modes, SST.lambdas[1:first_modes], type="l", xlab="Modes",
       ylab="Frac. Var. explained",las = 1)
  mtext('April-September SST anomalies', side=3, line = 1)
  points(1:first_modes, SST.lambdas[1:first_modes], col="red")

  
  ############## Plot eigen loadings and the time series for the first PCs
  plot( BBIdomain )
  xlong    <-  sort(unique(xygrid[,1]))
  ylat     <-  sort(unique(xygrid[,2]))
  
  par(mfrow = c(4, 2))
  par(mar = c(2, 3, 2, 1))
  for(imode in 1){ # Start loop over first modes
    
    zfull       <- rep(NaN,Nglobe)   #also equal 72*36
    zfull[indexsst.filt]=zsvd$u[,imode]
    zmat = matrix(zfull,nrow=Nx_sst,ncol=Ny_sst)
    image.plot(xlong,ylat,zmat,ylim=range(-60,70),ann=FALSE,las=1)
    mtext(paste("EOF no.",imode,sep=""), side = 3, line = 0.2, cex = 0.8)
    contour(xlong,ylat,(zmat),ylim=range(-60,70),add=TRUE,nlev=6,lwd=1.5)
    map("world2",wrap=c(-360,0),interior=F,add=TRUE,lwd=1)
    box()
    
    plot(Matrix.Rows, scale(SST.PCs[,imode]),type="l",xlab="Year",las=1,ylab="")
    mtext(paste("PC no.",imode,sep=""), side = 3, line = 0.2, cex = 0.8)
    
  }# End loop over first modes
  
  }

{
  dfVar = cbind( )
  corvar
  ##
  library(sp)
  data(meuse)
  coordinates(meuse) = ~x+y
  # let's do some manual fitting of two direct variograms and a cross variogram
  g <- gstat(id = "ln.zinc", formula = log(zinc)~1, data = meuse)
  g <- gstat(g, id = "ln.lead", formula = log(lead)~1, data = meuse)
  # examine variograms and cross variogram:
  plot(variogram(g))
  # enter direct variograms:
  g <- gstat(g, id = "ln.zinc", model = vgm(.55, "Sph", 900, .05))
  g <- gstat(g, id = "ln.lead", model = vgm(.55, "Sph", 900, .05))
  # enter cross variogram:
  g <- gstat(g, id = c("ln.zinc", "ln.lead"), model = vgm(.47, "Sph", 900, .03))
  # examine fit:
  plot(variogram(g), model = g$model, main = "models fitted by eye")
  # see also demo(cokriging) for a more efficient approach
  g["ln.zinc"]
  g["ln.lead"]
  g[c("ln.zinc", "ln.lead")]
  ##
  
  
  
  h <- hist(dfVar$value, freq = FALSE, col = "gray", ylim = c(0, 1))
  lines(h$mids, cumsum(h$counts)/sum(h$counts) )
  # Semivariogram as a lag of functions: x, y and t.
  g = gstat(NULL, "log-zinc", value ~ x,data = dfRCH,  )
  variogram(g)
  varX = variogram( value ~ x , locations = x ,data = na.omit(dfRCH) )
  plot( varX )
  

  
  # Cross correlogram with spatial variable precipitation RCH( u, t) and RCH(u,t-1),..., . u = (x,y)
  
  
  distVgm = variogram( )
  x = data[,2]
  y = data[,3]
  !is.na(x)
  !is.na(y)
  acf(x, y,type = "covariance", lag.max = 60, na.action = na.pass, plot= TRUE)
  ?acf
}
 ##
# Correction for well pumping

# Correction for 
##
if (FALSE){
  # Provisory variograms
  df_s = dataframe( "x"= ozone.UTM@coords[,1],"y" = ozone.UTM@coords[,2], "z"= ozoneDF)
  df_t = dataframe( "x"= ozone.UTM@coords[,1],"y" = ozone.UTM@coords[,2], "t"= ozoneDF)
  ozone_var_s <- gstat::variogram(log(zinc) ~ 1, ozoneDF)
  ozone_var_t <- gstat::variogram(log(zinc) ~ 1, ozoneTM)
  var_spatial <- fit.variogram(PPB~1,data=ozoneDF,tunit="hours",assumeRegular=F,na.omit=T) 
  var_temporal <- fit.variogram(PPB~1,data=ozoneTM,tunit="hours",assumeRegular=F,na.omit=T) 
}


# Plotting the Variogram
var <- variogramST(PPB~1,data=timeDF,tunit="hours",assumeRegular=F,na.omit=T, cores =4) 

system.time({
var2 <- variogramST(PPB~1,data=timeDF,tunit="hours",assumeRegular=F,na.omit=T, cores = 4,
                    twindow =30) 
})


plot(var,map=F) 
plot(var,map=T) 

plot(var,wireframe=T) 


# Spatial correlation and correlation test
{
  library(pheatmap)
  library( lattice )
  
  # Define input data
  data = dataInterpolated[,2:8]
  
  # Spearman
  # Compute correlation matrix
  corrMatrix = round(cor(data, use = "pairwise.complete.obs",method = "spearman" ),2)
  
  # Perform a correlation test
  cor.test_result = array( "", dim= dim(corrMatrix) )
  for ( istation in 1:nstations){
    for ( jstation in  1:nstations){
      
      if ( istation <= jstation){
        x = data[,istation]
        y = data[,jstation]
        
        
        aux =  cor.test( x, y, alternative = "two.sided", method = "spearman", conf.level = 0.45, na.action = na.omit )
        if ( !is.na(aux$p.value ) & ( aux$p.value < 0.05 | aux$p.value > 0.95 ) ){
          value = sprintf( fmt = "%.2f", corrMatrix[ istation, jstation ] ) 
          cor.test_result[istation,jstation]  = value
        }
      }
    }
  }
  
  
  # Plot heatmap
  legend_breaks = c(-1,-0.95,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,0.95,1)
  legend_labels  = legend_breaks
  breaks= legend_breaks
  color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                           "RdBu")))( length(breaks) - 1 )
  pheatmap( corrMatrix, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = cor.test_result,
            number_format = "%.2f", number_color = "black", legend_breaks = legend_breaks, legend_labels  = legend_labels,
            breaks = breaks, color = color)
  
  library(grid)
  grid.text(1:(dim( corrMatrix )), x=rep(seq(0.05, 0.91, length.out=10), 7), 
            y=rep(seq(0, 1, 0.1)+0.05, each=7))
  
  #add text elements
  if ( FALSE ){
  for ( istation in 1:nstations){
    for ( jstation in  1:nstations){
      if ( istation <= jstation){
        x = 1
        y = 1
        text(x=x, y=y, sprintf( "%.2f", corrMatrix[ istation, jstation ] ), cex = 2 )
      }
    }
  }
  }
  

  
  distMatrix = distSpatialPoints( shp )
  corrMatrix2 = corrMatrix
  corrMatrix2[ corrMatrix2 = 1 ] = NA
  plot(  distMatrix, corrMatrix2 )
  abline( h=0 )
}



       
       
tlags =c( 31,28,31,30,31,30,31,31,30,31,30,31)
diag = sqrt( ( extent( shp )[2] - extent( shp )[1]  )^2 + ( extent( shp )[4] - extent( shp )[3]  )^2  )

cutoff = diag
width = cutoff/15
boundaries = seq(0, cutoff, width)
var_exp <- variogramST(value~1,cutoff = cutoff, width = width, boundaries = boundaries,
                       data=dfST,tunit="days",assumeRegular=F,na.omit=T,
                       twindow = 30/3) 
plot(var_exp,map=F) 

var_exp2 = var_exp
idx = !is.na( var_exp2$gamma )
length( idx )
var_exp2 = var_exp2[  idx, ] 
plot(var_exp2,wireframe=T) 


# Exploratory Analysis ----------------------------------------------------
















df_results = array(NA,dim=c( nrows, ncols))
for (irow in 1:nrows){
  vgm_spatial = vgm( (pars.u - pars.l)[1] /2, as.character(models[irow,1]) , (pars.u - pars.l)[2] /2,
                     (pars.u - pars.l)[3] /2)
  vgm_temp = vgm( (pars.u - pars.l)[4] /2, as.character(models[irow,2]) , (pars.u - pars.l)[5] /2,
                  (pars.u - pars.l)[6] /2)
  vgm_joint = vgm( (pars.u - pars.l)[7] /2, as.character(models[irow,3]) , (pars.u - pars.l)[8] /2,
                   (pars.u - pars.l)[9] /2)
  
  if ( FALSE ){
  # 1: Separable
  separable <- vgmST("separable", space = vgm_spatial,time = vgm_temp, sill = vgm_joint$psill[2], stAni=500) 
  separable_Vgm <- fit.StVariogram(var, separable,  fit.method = 6 ,method="L-BFGS-B",lower=pars.l,upper=pars.u,tunit="hours")
  attr(separable_Vgm, "MSE")
  
  
  # 2: Sum metric
  sumMetric <- vgmST("sumMetric", space = vgm_spatial,time = vgm_temp, joint = vgm_joint, stAni=500) 
  sumMetric_Vgm <- fit.StVariogram(var, sumMetric,  fit.method = 10 , method="L-BFGS-B",lower=pars.l,upper=pars.u,tunit="hours")
  attr(sumMetric_Vgm, "MSE")
  
  
  # 3: Simple sum metric
  SimplesumMetric <- vgmST("simpleSumMetric",space = vgm_spatial, time = vgm_temp, joint = vgm_joint, nugget=1, stAni=500) 
  SimplesumMetric_Vgm <- fit.StVariogram(var, SimplesumMetric,  fit.method = 10 ,method = "L-BFGS-B",lower=pars.l)
  attr(SimplesumMetric_Vgm, "MSE")
  }
  
  
  ## Calibration
  if ( TRUE ){
  pars = (pars.u - pars.l) /2
  upper = pars.u
  lower = pars.l
  vgmtype = c( toString( models[irow,1] ),  toString( models[irow,2] ),  toString( models[irow,3] )    )
  covtype = "separable"
  obs = var$gamma
  result = sceuaWTime( ObjFun_MSE, pars = pars, lower, upper, maxn = 20000, kstop = 5, pcento = 5,
                       ngs = 20, npg = 5, nps = 5, nspl = 5, mings = 5, iniflg = 1, iprint = 0, iround = 3,
                       timeout = 60, peps = 0.0001,  obs = obs, vgmtype = vgmtype, covtype = covtype,
                       mult=10^6)
  
  pars = result$par
  vgm_spatial = vgm( pars[1],  vgmtype1,
                     pars[2], pars[3] )
  vgm_temp = vgm( pars[4],  vgmtype2,
                  pars[5], pars[6] )
  vgm_joint = vgm( pars[7], vgmtype3, 
                   pars[8], pars[9] )
  # 1: Separable
  separable_Vgm <- vgmST("separable", space = vgm_spatial,time = vgm_temp, sill = vgm_joint$psill[2], stAni=pars[10]) 
  
  
  pars = (pars.u - pars.l) /2
  upper = pars.u
  lower = pars.l
  vgmtype =  c( toString( models[irow,1] ),  toString( models[irow,2] ),  toString( models[irow,3] )    )
  covtype = "sumMetric"
  obs = var$gamma
  result = sceuaWTime( ObjFun_MSE, pars = pars, lower, upper, maxn = 20000, kstop = 5, pcento = 5,
                       ngs = 20, npg = 5, nps = 5, nspl = 5, mings = 5, iniflg = 1, iprint = 0, iround = 3,
                       timeout = 60*10, peps = 0.0001,  obs = obs, vgmtype = vgmtype, covtype = covtype,
                       mult=10^6)


  pars = result$par
  vgm_spatial = vgm( pars[1],  vgmtype1,
                     pars[2], pars[3] )
  vgm_temp = vgm( pars[4],  vgmtype2,
                  pars[5], pars[6] )
  vgm_joint = vgm( pars[7], vgmtype3, 
                   pars[8], pars[9] )
  # 2: Sum metric
  sumMetric_Vgm <- vgmST("sumMetric", space = vgm_spatial,time = vgm_temp, joint = vgm_joint, stAni=pars[10]) 

  
  
  pars = (pars.u - pars.l) /2
  upper = pars.u
  lower = pars.l
  vgmtype =  c( toString( models[irow,1] ),  toString( models[irow,2] ),  toString( models[irow,3] )    )
  covtype = "simpleSumMetric"
  obs = var$gamma
  result = sceuaWTime( ObjFun_MSE, pars = pars, lower, upper, maxn = 20000, kstop = 5, pcento = 5,
                       ngs = 20, npg = 5, nps = 5, nspl = 5, mings = 5, iniflg = 1, iprint = 0, iround = 3,
                       timeout = 60*10, peps = 0.0001,  obs = obs, vgmtype = vgmtype, covtype = covtype,
                       mult=10^6)
  
  pars = result$par
  vgm_spatial = vgm( pars[1],  vgmtype1,
                     pars[2], pars[3] )
  vgm_temp = vgm( pars[4],  vgmtype2,
                  pars[5], pars[6] )
  vgm_joint = vgm( pars[7], vgmtype3, 
                   pars[8], pars[9] )
  # 3: Simple sum metric
  SimplesumMetric_Vgm <- vgmST("simpleSumMetric",space = vgm_spatial, time = vgm_temp, joint = vgm_joint, nugget=1, stAni=pars[10]) 
  }
  
  # Save results
  if (FALSE){
  df_results[irow,1] = attr(separable_Vgm, "MSE")
  df_results[irow,2] = attr(sumMetric_Vgm, "MSE")
  df_results[irow,3] = attr(SimplesumMetric_Vgm, "MSE")
  }
  
  dist_grid = as.data.frame( var )
  
  if ( TRUE ){
  var_model <- variogramSurface( separable_Vgm, dist_grid = dist_grid )
  colnames( var_model )[ dim( var_model )[2] ] = "gamma_model"
  
  df_results[irow,1] = MSE( var_model$gamma, var_model$gamma_model )
  
  var_model <- variogramSurface( sumMetric_Vgm, dist_grid = dist_grid )
  colnames( var_model )[ dim( var_model )[2] ] = "gamma_model"
  
  df_results[irow,2] = MSE( var_model$gamma, var_model$gamma_model )
  }
  var_model <- variogramSurface( SimplesumMetric_Vgm, dist_grid = dist_grid )
  colnames( var_model )[ dim( var_model )[2] ] = "gamma_model"
 
  df_results[irow,3] = MSE( var_model$gamma, var_model$gamma_model )
  
  print( paste0("row ",irow,"/",nrows) )
}

which( df_results == min( df_results ), arr.ind = TRUE)


# results 
vgm_spatial = vgm( (pars.u - pars.l)[1] /2, as.character(models[irow,1]) , (pars.u - pars.l)[2] /2,
                   (pars.u - pars.l)[3] /2)
vgm_temp = vgm( (pars.u - pars.l)[4] /2, as.character(models[irow,2]) , (pars.u - pars.l)[5] /2,
                (pars.u - pars.l)[6] /2)
vgm_joint = vgm( (pars.u - pars.l)[7] /2, as.character(models[irow,3]) , (pars.u - pars.l)[8] /2,
                 (pars.u - pars.l)[9] /2)

## Calibration
if ( TRUE ){
  pars = (pars.u - pars.l) /2
  upper = pars.u
  lower = pars.l
  vgmtype = c( toString( models[irow,1] ),  toString( models[irow,2] ),  toString( models[irow,3] )    )
  covtype = "separable"
  obs = var$gamma
  result = sceuaWTime( ObjFun_MSE, pars = pars, lower, upper, maxn = 20000, kstop = 5, pcento = 5,
                       ngs = 20, npg = 5, nps = 5, nspl = 5, mings = 5, iniflg = 1, iprint = 0, iround = 3,
                       timeout = 60, peps = 0.0001,  obs = obs, vgmtype = vgmtype, covtype = covtype,
                       mult=10^6)
  
  pars = result$par
  vgm_spatial = vgm( pars[1],  vgmtype1,
                     pars[2], pars[3] )
  vgm_temp = vgm( pars[4],  vgmtype2,
                  pars[5], pars[6] )
  vgm_joint = vgm( pars[7], vgmtype3, 
                   pars[8], pars[9] )
  # 1: Separable
  separable_Vgm <- vgmST("separable", space = vgm_spatial,time = vgm_temp, sill = vgm_joint$psill[2], stAni=pars[10]) 
  
}

#demo(stkrige) 
par( mfrow = c(2,1) )
sp::plot(var,separable_Vgm,all = T, map=F) 
plot(var,map=F) 


# How much distance?; How much time in time window?
# Kriging without function is block kriging



#### IDW ####
# Example of interpolation
surfaceIDW <- idw(formula =  Time01 ~ 1, locations = shp, 
                  newdata = grid, nmax=10, maxdist=50000, idp=2)




#### Validation of model: LOOCV ####
nstations = dim(shp@data)[1]
ntime = dim(shp@data)[]
maxd = 5000
idp = 2

IDWoptim = function( maxd, idp ){
  error = array( NA, dim= c(nstations,1) )
  for ( istation in 1:nstations){
    shpLOOCV = shp[-istation,]
    gridLOOCV = shp[ istation, ]
    
    for (itime in 1){
      surfaceIDW <- idw(formula =  Time01 ~ 1 ,
                        locations = shpLOOCV, 
                        newdata = gridLOOCV, nmax=10, maxdist= maxd, idp= idp)
      
      error[ istation, itime ] = gridLOOCV$Time01 -  surfaceIDW$var1.pred
    }
  }
  error = mean( error, na.rm = T)
  return (error)
}


library(rtop)
ObjFun_MSE = function( pars ){
  result = IDWoptim( maxd = pars[1], idp = pars[2] )
  return( result )
}

#Gives you a labeled distance matrix. You just need the upper or lower triangle:
dist_df <- as.data.frame(st_distance( sf::st_as_sf(shp)) )
range(dist_df[lower.tri(dist_df)])
maxd_lower = min( apply(dist_df, 2, FUN = function(x){ sort(x, decreasing = F)[2]}) )
maxd_upper = max( apply(dist_df, 2, FUN = function(x){ sort(x, decreasing = T)[2]}) )


X = c(5000,2)
lower = c( maxd_lower, 0.1)
upper = c( maxd_upper, 4)

system.time( {
  optim = sceuaWTime( ObjFun_MSE, pars = X, lower, upper, maxn = 50000, kstop = 5, pcento = 0.01,
                      ngs = 5, npg = 5, nps = 5, nspl = 5, mings = 5, iniflg = 1, iprint = 0, iround = 3,
                      peps = 0.0001
  )
})



surfaceIDW <- idw(formula =  Time01 ~ 1, locations = shp, 
                  newdata = grid, nmax=10, maxdist=optim$par[1], idp=optim$par[2])

surfaceIDW <- idw(formula =  Time01 ~ 1, locations = shp, 
                  newdata = grid, nmax=10, maxdist=optim$par[1], idp=optim$par[2])

# BC
file = "C:/Users/josej/JupyterLab/MODFLOW Flopy/Example1/MF6/Ovalle/Quebrada El Ingenio/Shapefiles/QuebradaElIngenio_Extendida_OutFlowBC.shp"
boundary = readOGRv2 ( file )
icells = which( rgeos::gIntersects( BBIdomain, boundary, byid = TRUE ) )


# This is the estimated head at BC
surfaceIDW$var1.pred[icells]

aux = sf::st_as_sf(surfaceIDW)
plot( aux["var1.pred"])

icells =  rgeos::gNearestPoints( surfaceIDW[!is.na(surfaceIDW$var1.pred),], boundary ) 
idx = which.min( st_distance( st_as_sf(surfaceIDW),
                              st_as_sf(icells[1,]), by_element = TRUE) )
surfaceIDW$var1.pred[idx]
#  124.49 [masl]

# Plot
BBIdomain@data$surfaceIDW = surfaceIDW@data$var1.pred
aux = sf::st_as_sf(BBIdomain)


nbreaks = 10
mycol = brewer.pal(n = nbreaks, name = "Blues")
plot( aux["surfaceIDW"], col=mycol)
plot( aux["surfaceIDW"], axes = TRUE, key.pos = 4, nbreaks = nbreaks  +1,
      pal = mycol, key.width = lcm(4.5))

library( )
surfaceIDW@data$var1.pred
aux = raster::po(surfaceIDW)
sf::as_Spatial()

# Kriging
plot(lzn.vgm$dist, lzn.vgm$gamma )
lzn.fit <- fit.variogram(lzn.vgm, model=vgm(1, "exp", 900, 1)) # fit model
#show.vgms()
plot(lzn.vgm, lzn.fit)

colnames( shp@data )
#coordinates(meuse.grid) <- ~ x + y # step 3 above
lzn.kriged <- krige( 2016-01-01 ~ 1, shp, meuse.grid, model=lzn.fit)

lzn.kriged %>% as.data.frame %>%
  ggplot(aes(x=x, y=y)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="red") +
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw()

##### Topographic map ######

file = "C:\\Users\\josej\\JupyterLab\\MODFLOW Flopy\\Example1\\MF6\\Ovalle\\Quebrada El Ingenio\\DEM\\Patch.tif"
ras = raster( file )

file = "C:\\Users\\josej\\JupyterLab\\MODFLOW Flopy\\Example1\\MF6\\Ovalle\\Quebrada El Ingenio\\Shapefiles\\QuebradaElIngenio_Extendida_v2.shp"
shp = readOGRv2 ( file )

file = "C:\\Users\\josej\\JupyterLab\\MT3D USGS\\Examples\\Mine tailings deposit\\Coords\\Proyected_tailings_dam_v2.shp"
shp = readOGRv2 ( file )
#ras2 = clip_raster2( ras, shp )
ras2 = raster::crop( ras, extent(shp) )
plot( ras2 )

hist(ras)
hist(ras2)


plot( ras2 )
DrawTopographicMap(raster = ras2 , dz = 5 )

