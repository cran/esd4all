constructESD <- function(location,plot=TRUE,
                         get.data="data(esdsummary,envir=environment())",
                         mfrow=c(2,2)) {
  
  if (is.character(get.data)) eval(parse(text=get.data)) else
                              esdsummary <- get.data
  X <- esdsummary$X
  yr <- esdsummary$X + esdsummary$X.0

#  print(location)
  imatch <- is.element(esdsummary$locations,location)
#  print(summary(esdsummary)); print(sum(imatch))
  lon <- esdsummary$lons[imatch]
  lat <- esdsummary$lats[imatch]
  alt <- esdsummary$alts[imatch]
  ele <- esdsummary$eles[imatch]
  coefs <- esdsummary$coefs[imatch,,,]
#  print(ele); print(dim(coefs))
  if (plot) par(mfrow=mfrow)
  # print("constructESD"); print(dim(coefs)); print(c(coefs))

  coefs[is.na(coefs)] <- 0
  sesong <- c("December-February","March-May","June-August","September-November")
  coluns <- switch(as.character(ele),"101"="pink","601"="lightblue")
  coltrn <- switch(as.character(ele),"101"="darkred","601"="darkblue")
  ylab <- switch(as.character(ele),"101"="Temperature (deg C)",
                                   "601"="Precipitation (mm/season)")
  ave <- rep(NA,length(X)*4); dim(ave) <- c(length(X),4); q05 <- ave; q95 <- ave
  for (is in 1:4) {
    fit.mean <- coefs[1,1,is] +     coefs[2,1,is]*X +   coefs[3,1,is]*X^2 +
                coefs[4,1,is]*X^3 + coefs[5,1,is]*X^4 + coefs[6,1,is]*X^5
    fit.q05  <- coefs[1,2,is] +     coefs[2,2,is]*X +   coefs[3,2,is]*X^2 +
                coefs[4,2,is]*X^3 + coefs[5,2,is]*X^4 + coefs[6,2,is]*X^5
    fit.q95  <- coefs[1,3,is] +     coefs[2,3,is]*X +   coefs[3,3,is]*X^2 +
                coefs[4,3,is]*X^3 + coefs[5,3,is]*X^4 + coefs[6,3,is]*X^5

    ave[,is] <- fit.mean; q05[,is] <- fit.q05; q95[,is] <- fit.q95
     
    if (plot) {
      plot(range(yr),range(c(fit.q05,fit.q95)),ylim=range(c(fit.q05,fit.q95)),
           type="n",main=paste(location,": ",sesong[is],sep=""),ylab=ylab,xlab="Year")
      grid()
      polygon(c(yr,reverse(yr)),
              c(fit.q05,reverse(fit.q95)),lty=2,lwd=1,col=coluns,border="grey")
      lines(yr,fit.mean,lty=1,lwd=7,col=coltrn)
    }
  }
  results <- list(year=yr,q05=q05,ave=ave,q95=q95)
  invisible(results)
  
}



pdfESD <- function(location,plot=TRUE,get.data="data(esdsummary,envir=environment())",
                   year=2050,ref=NULL,mfrow=c(2,2),what="pdf") {
  if (is.character(get.data)) eval(parse(text=get.data)) else
                              esdsummary <- get.data
#  print(location)
  imatch <- is.element(esdsummary$locations,location)
#  print(summary(esdsummary)); print(sum(imatch))
  lon <- esdsummary$lons[imatch]
  lat <- esdsummary$lats[imatch]
  alt <- esdsummary$alts[imatch]
  ele <- esdsummary$eles[imatch]
  coefs <- esdsummary$coefs[imatch,,,]
    if (is.null(ref)) {
    now <- date()
    ref <- as.integer(substr(now,21,24))
  }
#  print(ele); print(dim(coefs))
  if (plot) par(mfrow=mfrow)
  # print("constructESD"); print(dim(coefs)); print(c(coefs))

  X <- year - esdsummary$X.0; X0 <- ref - esdsummary$X.0
  coefs[is.na(coefs)] <- 0
  sesong <- c("December-February","March-May","June-August","September-November")
  coltrn <- switch(as.character(ele),"101"="darkred","601"="darkblue")
  ylab <- switch(as.character(ele),"101"="Temperature (deg C)",
                                   "601"="Precipitation (mm/season)")
  ave <- rep(NA,length(X)*4); dim(ave) <- c(length(X),4); q05 <- ave; q95 <- ave
  for (is in 1:4) {
    fit.mean.0 <- coefs[1,1,is] +     coefs[2,1,is]*X0 +   coefs[3,1,is]*X0^2 +
                coefs[4,1,is]*X0^3 + coefs[5,1,is]*X0^4 + coefs[6,1,is]*X0^5
    fit.q95.0  <- coefs[1,3,is] +     coefs[2,3,is]*X0 +   coefs[3,3,is]*X0^2 +
                coefs[4,3,is]*X0^3 + coefs[5,3,is]*X0^4 + coefs[6,3,is]*X0^5

    fit.mean.x <- coefs[1,1,is] +     coefs[2,1,is]*X +   coefs[3,1,is]*X^2 +
                coefs[4,1,is]*X^3 + coefs[5,1,is]*X^4 + coefs[6,1,is]*X^5
    fit.q95.x  <- coefs[1,3,is] +     coefs[2,3,is]*X +   coefs[3,3,is]*X^2 +
                coefs[4,3,is]*X^3 + coefs[5,3,is]*X^4 + coefs[6,3,is]*X^5

    SD.0 <- (fit.q95.0 - fit.mean.0)/1.64  # From Table B.1 in Wilks (1995) p. 434
    SD.x <- (fit.q95.x - fit.mean.x)/1.64  # From Table B.1 in Wilks (1995) p. 434
    x.min <- min(c(fit.mean.0 - 3*SD.0,fit.mean.x - 3*SD.x))
    x.max <- max(c(fit.mean.0 + 3*SD.0,fit.mean.x + 3*SD.x))
    x <- seq(x.min,x.max,by=0.1)
    if (plot) {
      if (what=="pdf") {
        plot(x,dnorm(x=x,mean=fit.mean.0,sd=SD.0),type="l",
             lwd=5,col="grey80",main=paste(location,": ",sesong[is],sep=""),
             xlab="Temperature",ylab=ylab,
             sub=paste("Year",X+esdsummary$X.0," &",ref))
        grid()
        lines(x,dnorm(x=x,mean=fit.mean.x,sd=SD.x),
              lty=1,lwd=7,col=coltrn)
      } else if (what=="cdf") {
        plot(x,pnorm(q=x,mean=fit.mean.0,sd=SD.0),type="l",
             lwd=5,col="grey80",main=paste(location,": ",sesong[is],sep=""),
             xlab="Temperature",ylab=ylab,
             sub=paste("Year",X+esdsummary$X.0," &",ref))
        grid()
        lines(x,pnorm(q=x,mean=fit.mean.x,sd=SD.x),
              lty=1,lwd=7,col=coltrn)
      }
    }
  }
  results <- list(ave=fit.mean.x,std=SD.x,ave.ref=fit.mean.0,std=SD.0)
  invisible(results)  
}


mapESDlocs <- function(get.data="data(esdsummary,envir=environment())") {
  if (is.character(get.data)) eval(parse(text=get.data)) else
                              esdsummary <- get.data
  plot(esdsummary$lons,esdsummary$lats,pch=19,cex=0.6,col="red")
  addland()
}

queryLocations <- function(nr=NULL,get.data="data(esdsummary,envir=environment())") {
  if (is.character(get.data)) eval(parse(text=get.data)) else
                              esdsummary <- get.data
  if (!is.null(nr)) esdsummary$locations <- esdsummary$locations[nr]
  esdsummary$locations
}


get5mintopo <- function(browser = "firefox", url = "http://marine.rutgers.edu/po/tools/gridpak/etopo5.nc") {
# get the binary from http://www.ngdc.noaa.gov/mgg/global/etopo5.HTML
 print(paste("Download netCDF file from: ",
             "'http://marine.rutgers.edu/po/tools/gridpak/etopo5.nc'",
             "before proceeding - hopefully a browser will start and download the file..."))
 system(paste(browser, url))
 readLine("press a key when 'etopo5.nc' is located in working directory")
 require(ncdf)
 ncid <- open.ncdf("etopo5.nc")
 lon.z <- get.var.ncdf(ncid,"topo_lon")
 lat.z <- get.var.ncdf(ncid,"topo_lat")
 topo <- get.var.ncdf(ncid,"topo")
 close.ncdf(ncid)
 srtx <- order(lon)
 lon <- lon[srtx]
 topo <- topo[srtx,]
 image(lon.z,lat.z,topo,col=topo.colors(100))
 save(file="data/etopo5.Rdata",lon.z,lat.z,topo)
}


fortegn <- function(a,b) {
  if (a < b) fortegn <- -1 else fortegn <- 1
  fortegn
}


# Fits an empirical model based on distance to coast, altitude, latitude, and longitude
# to geographically distributed parameters.
# R.E. Benestad, 17.10.02


KrigFields <- function(resid,lon.grd,lat.grd) {
  print("Kriging based on 'fields'")
  require(fields)
  x <- cbind(resid$x,resid$y)
  k <- Krig(x=x,Y=resid$z)
  print(summary(k))
  
  glist<- list( x= lon.grd, y=lat.grd)

  out<- predict.surface( k, grid.list = glist )
  invisible(out)
}


KrigSgeostat <- function(resid,lon.grd,lat.grd,do.km) {
  print("Kriging based on 'sgeodat'")
  regiure(sgeostat)
  print("Apply kriging to the residual...")

  z.point<-point(resid)

  if (do.km) maxdist <- 400 else
             maxdist <- 3
#z.pair<-pair(z.point,type='anisotropic', theta=60, dtheta=10, maxdist=maxdist) 
  z.pair<-pair(z.point,num.lags=50,maxdist=maxdist) 
  print('est.variogram:')

  print(summary(resid)); print(dim(z.point)); print(length(z.pair))
  save(file="test.krig.Rdata",resid,z.point,z.pair,lon.grd,lat.grd)
  print(resid)

  z.v<-est.variogram(z.point,z.pair,'z')
  print("krig.vmod.l")
  z.vmod<-fit.variogram("spherical",z.v,plot.it=plot)
  #z.vmod<-fit.variogram("linear",z.v,plot.it=plot)

  print(summary(z.vmod))

  grd<-list(x=lon.grd,
            y=lat.grd)

  grd$xr<-range(grd$x)
  grd$xs<-grd$xr[2] - grd$xr[1]
  grd$yr<-range(grd$y)
  grd$ys<-grd$yr[2] - grd$yr[1]
  grd$max<-max(grd$xs, grd$ys)
  grd$xy<-data.frame(cbind(c(matrix(grd$x,length(grd$x),length(grd$y))),
                           c(matrix(grd$y,length(grd$x),length(grd$y),
                           byrow=TRUE))))
  colnames(grd$xy)<-c("x","y")
  grd$point<-point(grd$xy)

  # Krige
  print("Krige...")

  print(summary(grd$x))
  print(summary(grd$y))
  print(summary(z.point))

  grd$krige<-krige(grd$point,z.point,'z',z.vmod,maxdist=100,
                 extrap=TRUE)

  map.residuals <- list(x=grd$x,y=grd$y,
                        z=matrix(grd$krige$zhat,length(grd$x),length(grd$y)))
  invisible(map.residuals)
}



geo.inf <- function(g.obj,do.km=TRUE,x.scale=1000,
                    predict=TRUE,krig=TRUE,krig.Nx=NULL,krig.Ny=NULL,
                    x.rng=c(-10,32),y.rng=c(44,70),plot=FALSE,
                    krig.package="fields",
                    use.previous.estimates=TRUE,linear.intp=TRUE) {
  
require(clim.pact)
require(akima)
require(ncdf)

print("geo.inf:")

# Distance from coast:

print("Load geo-data...")
data(addland,envir=environment())

print("etopo5...")
if (file.exists("data/etopo5.Rdata")) load("data/etopo5.Rdata") else
                                      get5mintopo()

print("focus on region...")
# eg image(lon.z,lat.z,t(topo))
x.keep <- (lon.z >=  min(x.rng)) & (lon.z <= max(x.rng))
y.keep <- (lat.z >=  min(y.rng)) & (lat.z <= max(y.rng))

lon.z <- lon.z[x.keep]
lat.z <- lat.z[y.keep]
topo <- topo[x.keep,y.keep]

print(paste("keep",sum(x.keep),"&",sum(y.keep)))
print(paste("Lon range=",range(lon.z)[1],"-",range(lon.z)[2]))
print(paste("Lat range=",range(lat.z)[1],"-",range(lat.z)[2]))

print("Dimensions of geo-date:")
print(c(dim(topo)))

land <- (topo > 0)
np <- sum(land)
nx <- length(lon.z)
ny <- length(lat.z)
nxy <- nx*ny

good <- is.finite(g.obj$y)   &  is.finite(g.obj$lon) &
        is.finite(g.obj$lat) &  is.finite(g.obj$alt) &
        is.finite(g.obj$dist)&  is.finite(log(g.obj$dist)) &
        is.finite(sqrt(g.obj$dist)) &  is.finite(sqrt(g.obj$alt))

lon <- g.obj$lon; lat <- g.obj$lat
lon2x <- lon; lat2x <- lat
mlon <- median(lon); mlat <- median(lat)
if (do.km) {
  for (i in 1:length(lon)) {
    lon2x[i] <- fortegn(lon[i],mlon)*distAB(mlon,lat[i],lon[i],lat[i])/x.scale
    lat2x[i] <- fortegn(lat[i],mlat)*distAB(lon[i],mlat,lon[i],lat[i])/x.scale
  }
  good <- good & is.finite(lon2x) & is.finite(lat2x)
}


if (!do.km) {
    print("distance in 'degrees'")
    geo.par <- data.frame(y=g.obj$y[good],alt=g.obj$alt[good],
                          dist=g.obj$dist[good],lat=g.obj$lat[good],lon=g.obj$lon[good],
                          log.dist=log(g.obj$dist[good]), sqrt.dist=sqrt(g.obj$dist[good]),
                          sqrt.alt=sqrt(g.obj$alt[good]))
  } else {
    print("semi-distance in 'km'")
    geo.par <- data.frame(y=g.obj$y[good],alt=g.obj$alt[good],
                          dist=g.obj$dist[good],lat=lat2x[good],lon=lon2x[good],
                          log.dist=log(g.obj$dist[good]), sqrt.dist=sqrt(g.obj$dist[good]),
                          sqrt.alt=sqrt(g.obj$alt[good]))
}

print(summary(geo.par))
geo.mod <- lm(y ~ alt + dist + lat + lon + log.dist + sqrt.dist + sqrt.alt, data=geo.par) 

print("geo.mod calibrated")
print(summary(geo.mod))

isvalid <- is.finite(geo.par$dist) & is.finite(geo.par$alt) & is.finite(geo.par$lat) &
           is.finite(geo.par$lon) & is.finite(log(geo.par$dist)) &
           is.finite(sqrt(geo.par$dist)) & is.finite(sqrt(geo.par$alt))

#print(cbind(dist,alt,lat,lon))
Anova <- summary(geo.mod)
print(Anova)
coefs <- Anova$coefficients
coefs[!is.finite(coefs)] <- 0

print(paste("#Valid data points = ",sum(isvalid)," do.km=",do.km,
            "sum(good)=",sum(good),"sum(is.finite(lon2x))=",sum(is.finite(lon2x)),
            "sum(is.finite(lat2x))=",sum(is.finite(lat2x))))


if (!predict) return()

map <- topo
map[,] <- NA

map[land] <- coefs[1]
  
if (plot) {
  par(mfrow=c(3,2))
  image(lon.z,lat.z,map,main=paste("Constant for ",attr(g.obj,'Description')),
        col=topo.colors(100))
  grid()
  addland()
  points(g.obj$lon,g.obj$lat,pch=20,cex=0.5,col="white")
  points(g.obj$lon,g.obj$lat,pch=20,cex=0.3,col="black")
}

map[land] <- map[land] + topo[land]*coefs[2]  + sqrt(topo[land])*coefs[8]

if (plot) {
  image(lon.z,lat.z,map,main="Constant+altitude",col=topo.colors(100))
  grid()
  addland()
  contour(lon.z,lat.z,map,add=TRUE)
  points(g.obj$lon,g.obj$lat,pch=20,cex=0.5,col="white")
  points(g.obj$lon,g.obj$lat,pch=20,cex=0.3,col="black")
}
alt.grd <- topo[land]
dist.grd <- rep(NA,np); xkm.grd <- dist.grd; ykm.grd <- dist.grd

lat.grd <- t(matrix(rep(lat.z,nx),ny,nx)); lat.grd <- lat.grd[land]
lon.grd <- matrix(rep(lon.z,ny),nx,ny); lon.grd <- lon.grd[land]

lons <- lon.grd[land]; lats <- lat.grd[land]; dist <- dist.grd[land]
print(paste("Distance # data points =",sum(land)))

# Predict from geo.mod

N.land <- sum(land)
print(paste("Predict the geographically dependent part np=",np))
use.previous <- (use.previous.estimates) & (file.exists("geo.grid_previous.est.Rdata"))
if (use.previous) {
  print("Uses previous values to save time")
  load("geo.grid_previous.est.Rdata")
  N.land <- length(dist.grd)
}

if ( ((N.land != sum(land)) & (use.previous)) | (!use.previous) ) {
  N.land <- sum(land)
  print(paste("Estimate distance from coast, eastings and westings for",N.land,"grid boxes"))
  for (i in 1:N.land) {
    dist.grd[i]<-min(distAB(lon.grd[i],lat.grd[i],lon.cont,lat.cont),na.rm=TRUE)/x.scale
    xkm.grd[i] <-fortegn(lon.grd[i],mlon)*distAB(mlon,lat.grd[i],lon.grd[i],lat.grd[i])/x.scale
    ykm.grd[i] <-fortegn(lat.grd[i],mlat)*distAB(lon.grd[i],mlat,lon.grd[i],lat.grd[i])/x.scale
    if (mod(i,500)==0) {
      print(paste(round(100*i/sum(land)),"% of dist, x, & y"))
    }
  }
  save(file="geo.grid_previous.est.Rdata",dist.grd,xkm.grd,ykm.grd)
}

if (do.km) {
  lat.grd <- ykm.grd
  lon.grd <- xkm.grd
}

print("dist.grd")
map[land] <- map[land] + dist.grd*coefs[3] + log(dist.grd)*coefs[6] + sqrt(dist.grd)*coefs[7]
if (plot) {
  image(lon.z,lat.z,map,main="Constant+alt+dist",col=topo.colors(100))
  addland()
  grid()
  contour(lon.z,lat.z,map,add=TRUE)
  points(g.obj$lon,g.obj$lat,pch=20,cex=1.2,col="white")
  points(g.obj$lon,g.obj$lat,pch=20,cex=0.9,col="black")
}
print("lat.grd")
map[land] <- map[land] + lat.grd*coefs[4]

if (plot) {
  image(lon.z,lat.z,map,main="Constant+alt+dist+lat",col=topo.colors(100))
  addland()
  grid()
  contour(lon.z,lat.z,map,add=TRUE)
  points(g.obj$lon,g.obj$lat,pch=20,cex=1.2,col="white")
  points(g.obj$lon,g.obj$lat,pch=20,cex=0.9,col="black")
}
print("lon.grd")
map[land] <- map[land] + lon.grd*coefs[5]

if (plot) {
  x11()
  image(lon.z,lat.z,map,main="Constant+alt+dist+lat+lon",col=topo.colors(100))
  addland()
  grid()
  contour(lon.z,lat.z,map,add=TRUE)
  points(g.obj$lon,g.obj$lat,pch=20,cex=1.2,col="white")
  points(g.obj$lon,g.obj$lat,pch=20,cex=0.9,col="black")

  dev.copy2eps(file="geo.inf_GRM-wo-krig.eps")
}

geo.dat <- data.frame(dist=dist.grd,lat=lat.grd,lon=lon.grd,
                      alt=alt.grd,log.dist=log(dist.grd),sqrt.dist=sqrt(dist.grd))

print(summary(geo.dat))

print(paste("Vector lengths: ",
            length(geo.mod$residual),length(geo.par$lon),length(geo.par$lat)))
resid<-data.frame(z=geo.mod$residual,x=geo.par$lon[isvalid],y=geo.par$lat[isvalid]) 

save(file="geo.inf_GRM.Rdata",geo.mod,map,geo.dat,geo.par,isvalid) 

if (is.null(krig.Nx)) lon.grd <- rep(lon.z,length(lat.z)) else
                      lon.grd <- rep(seq(min(lon.z),max(lon.z),length=krig.Nx),krig.Ny)
if (is.null(krig.Ny)) lat.grd <- sort(rep(lat.z,length(lon.z))) else
                      lat.grd <- sort(rep(seq(min(lat.z),max(lat.z),length=krig.Ny),krig.Nx))

#  xkm.grd[i] <-fortegn(lon.grd[i],mlon)*distAB(mlon,lat.grd[i],lon.grd[i],lat.grd[i])/x.scale
#  ykm.grd[i] <-fortegn(lat.grd[i],mlat)*distAB(lon.grd[i],mlat,lon.grd[i],lat.grd[i])/x.scale

print(paste("length(lon.z)=",length(lon.z),"length(lon.grd)=",length(lon.grd),
            "length(lat.z)=",length(lat.z),"length(lat.grd)=",length(lat.grd),"x.scale=",x.scale))
# print(lat.z)

if (do.km) {
  for (i in 1:length(lon.grd)) {
    lon.grd[i] <- fortegn(lon.grd[i],mlon)*distAB(mlon,lat.grd[i],lon.grd[i],lat.grd[i])/x.scale
  }
  for (i in 1:length(lat.grd)) {
    lat.grd[i] <- fortegn(lat.grd[i],mlat)*distAB(lon.grd[i],mlat,lon.grd[i],lat.grd[i])/x.scale
  }
} 

lon.grd <- seq(min(lon.grd),max(lon.grd),length=length(lon.z))
lat.grd <- seq(min(lat.grd),max(lat.grd),length=length(lat.z))

#---------------------------------------------------------------------------

if (krig) {
  print("Apply kriging to the residual using the fields package")
  if (krig.package=="fields") grd <- KrigFields(resid,lon.grd,lat.grd) else
                              grd <- KrigSgeostat(resid,lon.grd,lat.grd,do.km)
  map.residuals <- grd$z

  if (plot) {
    x11()
    image(grd$x,grd$y,map.residuals,main=paste("Residuals",attr(g.obj,'Description')))
    points(resid$x,resid$y,col="blue")
    text(resid$x,resid$y,round(resid$z,2))
  }

  print(summary(c(map.residuals)))
  good <- is.finite(map.residuals)
  print(paste("geo.inf: There are ",sum(good),"datapoints with valid data after kriging"))
  if (sum(good) < 10) {
   print("Kriging failed - use interpolation instead!")
   grd <- interp(resid$x,resid$y,resid$z,lon.grd,lat.grd,
                 linear=FALSE,duplicate="mean")
   map.residuals <- grd$z
  }
  if ( (length(lon.z)!= length(lon.grd)) |
       (length(lat.z)!= length(lat.grd)) )  {
    print("Interpolate to finer grid using interp:")
    grd.xy <- rep(lon.grd,length(lat.grd))
    grd.yx <- sort(rep(lat.grd,length(lon.grd)))
#    grd <- interp(resid$x,resid$y,resid$z,lon.grd,lat.grd,
#                 linear=linear.intp,duplicate="mean")
    grd <- interp(grd.xy[good],grd.yx[good],map.residuals[good],
                 lon.z,lat.z,linear=linear.intp,duplicate="mean")
   map.residuals <- grd$z
  }

# -----------------------------------------------------------------------------------

} else {

# Bi-linear interpolation ------------------------------------------------------------
  
  print("Check for duplicates:")
  for (ii in 1:length(resid$z)) {
    match.coord <- is.element(resid$x,resid$x[ii]) &
                   is.element(resid$y,resid$y[ii])
    if (sum(match.coord)!=1) {
      iii <- (1:length(resid$z))[match.coord]
      print(paste(ii,") Non-unique corrdinates: x=",resid$x[ii]," & y=",resid$y[ii],sep=""))
      print(resid$z[match.coord]); print(iii)
      resid$z[iii[1]] <- mean(resid$z[match.coord],na.rm=TRUE)
      resid$z[iii[2:sum(match.coord)]] <- NA
    }
  }
  good <- is.finite(resid$z) & is.finite(resid$x) & is.finite(resid$y)
  grd <- interp(resid$x[good],resid$y[good],resid$z[good],
                lon.grd,lat.grd,linear=linear.intp)
  map.residuals <- grd$z; D <- dim(map.residuals)
  for (i in 1:D[1]) map.residuals[i,] <- gauss.filt(map.residuals[i,],5)
  for (j in 1:D[2]) map.residuals[,j] <- gauss.filt(map.residuals[,j],5)
  print(dim(map.residuals))
  
  if (plot) {
    x11()
    image(grd$x,grd$y,map.residuals,main="bi-linear interpolation of residuals")
    points(resid$x,resid$y,col="blue")
    text(resid$x,resid$y,round(resid$z,2))
  }

# ------------------------------------------------------------------------------------
}
# Plot figures:

print(paste("geo.inf: Add gridded residuals: map.residuals has",
            100*sum(is.finite(map.residuals))/length(map.residuals),
            "% and map.residuals has[land]",
            100*sum(is.finite(map.residuals[land]))/sum(land),
            " of valid data for", attr(g.obj,'Description')))
z <- map
z[land] <- map[land] + map.residuals[land]
z <- map + map.residuals

if (plot) {
  x11()
  image(lon.z,lat.z,z,main="Geo dep + Res krig",col=topo.colors(100))
  addland()
  grid()
  contour(lon.z,lat.z,map,add=TRUE)
  points(g.obj$lon,g.obj$lat,pch=20,cex=1.2,col="white")
  points(g.obj$lon,g.obj$lat,pch=20,cex=0.9,col="black")
  dev.copy2eps(file="geo.inf_GRM+krig.eps")


  
  x11()
  par(mfrow=c(2,1))

  image(lon.z,lat.z,map.residuals,
        main="Residuals from kriging",
        sub="Model: dist. to coast, altitude, latitude & longitude",
        xlab="Longitude (deg E)",ylab="Latitude (deg N)")
  addland()
  contour(lon.z,lat.z,map.residuals,add=TRUE)
  points(g.obj$lon,g.obj$lat,pch=20,cex=1.2,col="white")
  points(g.obj$lon,g.obj$lat,pch=20,cex=0.9,col="black")
  grid()
                                        #  dev.off()
  plot(geo.par$y,predict(geo.mod),pch=20,col="darkblue",
       main="Obs vs. predicted values",xlab="Obs",ylab="Predicted")
  grid()

  x11()
  par(mfrow=c(1,1))
  image(grd$x,grd$y,map.residuals,main="Residuals")
  points(resid$x,resid$y,col="blue")
  text(resid$x,resid$y,round(resid$z,2),cex=0.5,pos=1)

}

geo.inf <- list(x=lon.z,y=lat.z,z=z,
                map.residuals=map.residuals,geo.dep=map,
                lon=g.obj$lon,lat=g.obj$lat,alt=g.obj$alt,dist=g.obj$dist,
                model=geo.mod,x.scale=x.scale,
                lon.reference=mlon,lat.reference=mlat,R2=summary(geo.mod)$r.squared)
print(summary(geo.inf))
rm(topo,lon.z,lat.z)

#print(paste("Save 'geo.inf' - grid size=",length(z)))
#if (krig) save(file="geo.inf_GRM+res-krig.Rdata",geo.inf) else
#          save(file="geo.inf_GRM+res-spline.Rdata",geo.inf)

invisible(geo.inf)

}


gridESD <- function(get.data="data(esdsummary,envir=environment())",plot=FALSE,
                    x.rng=c(-30,50),y.rng=c(40,72),x.scale=1000,do.km=TRUE,krig=TRUE,
                    new=TRUE,krig.Nx=30,krig.Ny=30,use.previous.estimates=TRUE,
                    linear.intp=TRUE,krig.package="fields",fname="gridded.c.rda") {
#  x.rng=c(-15,45),y.rng=c(33,72)
#  x.rng=c(-10,40),y.rng=c(40,70)
#  x.rng=c(5,12),y.rng=c(58,63)
#  x.rng=c(5,19),y.rng=c(53,65)
#  x.rng=c(5,19),y.rng=c(53,65),plot=TRUE
#  x.rng=c(5,25),y.rng=c(53,67),plot=TRUE
  if (krig) grid.res.method <- "kriging" else grid.res.method <- "bi-linear interp."
  if (is.character(get.data)) eval(parse(text=get.data)) else
                              esdsummary <- get.data
  print(summary(esdsummary))
  n.coefs <- length(esdsummary$coefs[1,,,])
  d.coefs <- dim(esdsummary$coefs[1,,,])

# Unrealistic results: Kaunas, Vilnius and Klaipeda: winter mean temperatures ~ 10-15C!!!
  
  inregion <- ( (esdsummary$lons>=min(x.rng)) & (esdsummary$lons<=max(x.rng)) &
                (esdsummary$lats>=min(y.rng)) & (esdsummary$lats<=max(y.rng)) &
                 !is.element(esdsummary$locations,c("KAUNAS","KLAIPEDA","VILNIUS")) )
  inregion[is.na(inregion)] <- FALSE
  print(summary(esdsummary$lons))
  ns <- sum(inregion)

  print(paste("gridESD: ",ns," locations within ",x.rng[1],"-",x.rng[2],"E and ",
              y.rng[1],"-",y.rng[2],"N - gridding coefficients. plot=",plot,sep=""))
  locs <- esdsummary$locations[inregion]
  lon <- esdsummary$lons[inregion]
  lat <- esdsummary$lats[inregion]
  alt <- esdsummary$alts[inregion]
  coefs <- esdsummary$coefs[inregion,,,]
  dim(coefs) <- c(ns,n.coefs)

  GRM.R2 <- rep(0,n.coefs)
  coef.tag <- rep("",n.coefs); dim(coef.tag) <- d.coefs
  coef.tag[,1,] <- "mean"
  coef.tag[,2,] <- "q05"
  coef.tag[,3,] <- "q95"
  coef.tag[,,1] <- paste(coef.tag[,,1],"DJF",sep=".")
  coef.tag[,,2] <- paste(coef.tag[,,2],"MAM",sep=".")
  coef.tag[,,3] <- paste(coef.tag[,,3],"JJA",sep=".")
  coef.tag[,,4] <- paste(coef.tag[,,4],"SON",sep=".")
  for (i in 1:d.coefs[1]) coef.tag[i,,] <- paste("c",i,".",coef.tag[i,,],sep="")
  coef.tag <- c(coef.tag)
  
  load("data/etopo5.Rdata")
  inlon.z <- (lon.z>=min(x.rng)) & (lon.z<=max(x.rng))
  inlat.z <- (lat.z>=min(y.rng)) & (lat.z<=max(y.rng)) 
  topo <- topo[inlon.z,inlat.z]
  lon.z <- lon.z[inlon.z]
  lat.z <- lat.z[inlat.z]
  data(addland,envir=environment())
  good <- is.finite(lon.cont) & is.finite(lat.cont)
  lon.cont <- lon.cont[good]
  lat.cont <- lat.cont[good]

  mfrow=c(3,n.coefs/3)
  z.levs <- seq(-0.5,0.5,by=0.05); nl <- length(z.levs)
  my.col <- rgb(seq(0,1,length=nl),
                seq(0,1,length=nl),
                seq(0,1,length=nl))

  las <- 1

  
  dist <- rep(NA,ns)
  print("Derive distance-to-coast:")
  if (plot) {
    x11()
    plot(x.rng,y.rng,type="n")
    grid()
    addland()
  }
  for (il in 1:ns) {
      dd <- distAB(lon[il],lat[il],lon.cont,lat.cont)
      dist[il] <-  min(dd,na.rm=T)/x.scale
      nearest.coast <- (dd==min(dd,na.rm=TRUE))
      if (plot) {
        points(lon[il],lat[il],pch=18,col="darkred")
        lines(c(lon[il],lon.cont[nearest.coast]),
              c(lat[il],lat.cont[nearest.coast]),lty=2,col="red")
        points(lon.cont[nearest.coast],lat.cont[nearest.coast],pch=".",col="red")
        print(paste(locs[il],"lon=",round(lon[il],2),"lat=",round(lat[il],2),
                    "dist=",round(dist[il]*10,2),"alt=",round(alt[il],2)))
      }
  }


  ic <- 1
  while (ic <= n.coefs) {
    print(paste("Coefficient ",ic,"of",n.coefs,"tag=",coef.tag[ic]))
    
    metric <- coefs[,ic]
 
    if (sum(is.finite(metric))>0) {
      geo.par <- data.frame(y=metric,dist=dist,lat=lat,
                            lon=lon,alt=alt,location=locs)
      attr(geo.par,'Description') <- coef.tag[ic]


      print("-----------------geo.inf-------------------------:")
      geo <- geo.inf(geo.par,do.km=do.km,x.scale=x.scale,krig=krig,
                     krig.Nx=krig.Nx,krig.Ny=krig.Ny,x.rng=x.rng,
                     y.rng=y.rng,plot=plot,
                     use.previous.estimates=use.previous.estimates)
      print("geo.inf completed")
      if (ic==1) {
        print("First loop - set up gridded.c")
        if ( (file.exists(fname)) & (!new) ) {
          load(fname)
          d.c <- dim(gridded.c); d.g <- dim(geo$z)
          if ( (d.c[1]==d.g[1]) & (d.c[1]==d.g[1]) ) {
            dim(gridded.c) <- c(d.c[1]*d.c[2],d.c[3])
            ic <- max( (1:n.coefs)[sum(is.finite(colMeans(gridded.c,na.rm=TRUE)))] )
            dim(gridded.c) <- c(d.c[1],d.c[2],d.c[3])                      
            print("###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###")
            print(paste("### Reading",fname,"with",ic,"maps already done ###"))
           }
        }
        if (new) {
          print(paste("create 'gridded.c'; Grid dim:",dim(geo$z)[1],dim(geo$z)[2]))
          gridded.c <- rep(geo$z,n.coefs); dim(gridded.c) <- c(dim(geo$z),n.coefs)
          d.c <- dim(gridded.c)
        }
        print(paste("dim(gridded.c):",d.c[1],d.c[2],d.c[3],
                    "length(coef.tag):",length(coef.tag),"n.coefs=",n.coefs,
                    "number of dimensions=",length(d.c),"coef.tag=",coef.tag[ic]))
        dim(gridded.c) <- c(d.c[1]*d.c[2],d.c[3])
        colnames(gridded.c) <- c(coef.tag)
        dim(gridded.c) <- c(d.c[1],d.c[2],d.c[3])
      } else print("Update the 3D matrix 'gridded.c[,,ic]'")
      GRM.R2[ic] <- geo$R2
      gridded.c[,,ic] <- geo$z[,]
      attr(gridded.c,'longitudes') <- lon.z
      attr(gridded.c,'latitudes') <- lat.z
      attr(gridded.c,'residual_gridding_method') <- grid.res.method
      attr(gridded.c,'type') <- coef.tag
      attr(gridded.c,'GRM-R2') <- round(100*GRM.R2) 
      print(attributes(gridded.c))
      print(colnames(gridded.c))
      save(file=fname,gridded.c)
   
      if (plot) { 
       dev.copy2eps(file=paste("gridESD_",ic,"a.eps",sep="")); dev.off()
       dev.copy2eps(file=paste("gridESD_",ic,"b.eps",sep="")); dev.off()
      }
    }
#    print("define the nc-variables")
#    varnm <- paste("coef.",ic,sep="")
#    eval(parse(text=paste("vardef",ic,
#                   " <- var.def.ncdf(name=varnm,units='not-specified',",
#                   "dim=list(dimlon,dimlat), missval=-999,",
#                   "longname='coefficient ",ic," of 5th-order polynomial fit',",
#                   "prec='short')",sep="")))
    ic <- ic+1
    rm(geo)
  }

  
  gridded.c[!is.finite(gridded.c)] <- -999
  print("Set up ncdf- dimensions")
  dimlon <- dim.def.ncdf( "lon", "deg E", lon.z)
  dimlat <- dim.def.ncdf( "lat", "deg N", lat.z)
  dimn <- dim.def.ncdf( "coefnr", "none", 1:d.c[3])

  vardef <- var.def.ncdf(name="coef",units="none",dim=list(dimlon,dimlat,dimn), missval=-999, 
        longname="coefficients for 5th-order polynomials", prec="single")
  metadef <- var.def.ncdf(name="meta",units="none",dim=dimn, missval="-999", 
        longname="coefficients for 5th-order polynomials", prec="char")
  print("Create netCDF-file")
  ncnew <- create.ncdf("Europe_E-SDS_t2m-trend_map.nc", list(vardef,metadef))
  put.var.ncdf( ncnew, vardef, round(gridded.c*100) )
#  print("Add metadata:")
#  put.var.ncdf( ncnew, metadef, coef.tag)

  print("add attributes")
  att.put.ncdf( ncnew, vardef, 'scale_factor', 0.01, prec="double")
  att.put.ncdf( ncnew, 0, 'scenario', 'IPCC SRES A1b')
  att.put.ncdf( ncnew, 0, 'author', 'Rasmus Benestad')
  att.put.ncdf( ncnew, 0, 'title', 'Climate change scenarios from multi-model IPCC AR4 climate simulations')
  att.put.ncdf( ncnew, 0, 'R-package', 'esd4all')
  att.put.ncdf( ncnew, 0, 'date', date())
  att.put.ncdf( ncnew, 0, 'key_words','statistical-empirical, downscaling, multi-model ensemble, gridding')
  att.put.ncdf( ncnew, 0, 'ensemble','CMIP3')
  grid.res.method <- "bi-polar interpolation of residuals"
  att.put.ncdf( ncnew, 0, 'residual_gridding_method',grid.res.method)
  att.put.ncdf( ncnew, 0, 'GCM_source','PCMDI - Earth System Grid: URL https://esg.llnl.gov:8443/index.jsp')
  close.ncdf(ncnew)
}



mapESDquants <- function(what="q95",season=3,year=2050,ref=NULL,
                         get.data1="data(gridded.c,envir=environment())",
                         get.data2="data(esdsummary,envir=environment())",
                         plot=TRUE) {
  seasons <- c("DJF","MAM","JJA","SON")

  if (is.null(ref)) {
    now <- date()
    ref <- as.integer(substr(now,21,24))
  }
  #print(c(what,seasons[season]))
  if (is.character(get.data1)) eval(parse(text=get.data1)) else
                              gridded.c <- get.data1
  attr(gridded.c,'type') -> coef.tag
  selection <- grep(seasons[season],coef.tag)
  iC <- selection[grep(what,coef.tag[selection])]
  tags <- coef.tag[iC]
  #print(tags); print(iC); print(selection)
  if (is.character(get.data2)) eval(parse(text=get.data2)) else
                              esdsummary <- get.data2
  data(addland)
  if (is.null(esdsummary$X.0)) esdsummary$X.0 <- 1900
  X <- year - esdsummary$X.0
  #print(summary(esdsummary))
  d <- dim(gridded.c)
    map <- matrix(rep(0,d[1]*d[2]),d[1],d[2])
  # print(d); print(iC); print(dim(map)); print(dim(gridded.c[,,iC]))
  map[,] <- gridded.c[,,iC[1]]
  #print("Add all terms")
  for (ic in 2:6) {
    #print(c(ic,iC[ic],X,X^(ic-1)))
    map[,] <- map[,] + gridded.c[,,iC[ic]]*X^(ic-1)
  }
  lon.z <- attr(gridded.c,'longitudes')
  lat.z <- attr(gridded.c,'latitudes')
  if (plot) {
    x11()
    X.ref <- ref - esdsummary$X.0
    map.ref <- matrix(rep(0,d[1]*d[2]),d[1],d[2])
    map.ref[,] <- gridded.c[,,iC[1]]
    for (ic in 2:6) {
       map.ref[,] <- map.ref[,] + gridded.c[,,iC[ic]]*X.ref^(ic-1)
    }
    par(mfrow=c(1,3))

    zlim <- range(c(map.ref,map),na.rm=TRUE)
    image(lon.z,lat.z,map.ref,main=paste(seasons[season],what,ref),
        col=cm.colors(21),xlab="",ylab="",zlim=zlim)
    addland()
    contour(lon.z,lat.z,map.ref,add=TRUE)
    points(esdsummary$lons,esdsummary$lats,pch=19,cex=0.5,col="grey30")
    
    image(lon.z,lat.z,map,main=paste(seasons[season],what,year),
        col=cm.colors(21),xlab="",ylab="",zlim=zlim)
    addland()
    contour(lon.z,lat.z,map,add=TRUE)
    points(esdsummary$lons,esdsummary$lats,pch=19,cex=0.5,col="grey30")

    z1 <- max(abs(c(map-map.ref)),na.rm=TRUE)
    image(lon.z,lat.z,map-map.ref,main=paste("Difference",what,year,"-",ref),
        col=cm.colors(21),xlab="",ylab="",zlim=c(-z1,z1))
    addland()
    contour(lon.z,lat.z,map-map.ref,add=TRUE)
    points(esdsummary$lons,esdsummary$lats,pch=19,cex=0.5,col="grey30")
    dev.copy2eps(file=paste("mapESDquants_",what,"-",seasons[season],
                   "-",year,"-",ref,".eps",sep=""))
  }
#  print("finished")
  attr(map,'longitudes') <- lon.z
  attr(map,'latitudes') <- lat.z
  attr(map,'season') <- seasons[season]
  attr(map,'what') <- what
  attr(map,'year') <- year
  invisible(map)
}

mapESDprobs <- function(thresh=0,season=1,year=2050,ref=NULL,
                        get.data="data(gridded.c,envir=environment())",plot=TRUE) {

  seasons <- c("DJF","MAM","JJA","SON")
  map.q95 <- mapESDquants(season=season,year=year,get.data1=get.data,plot=FALSE)
  map.ave <- mapESDquants(what="mean",season=season,year=year,
                          get.data1=get.data,plot=FALSE)
  map.std <- (map.q95 - map.ave)/1.64  # From Table B.1 in Wilks (1995) p. 434
  
  if (is.null(ref)) {
    now <- date()
    ref <- as.integer(substr(now,21,24))
  }

  map.q95.ref <- mapESDquants(season=season,year=ref,get.data1=get.data,plot=FALSE)
  map.ave.ref <- mapESDquants(what="mean",season=season,year=ref,
                          get.data1=get.data,plot=FALSE)
  map.std.ref <- (map.q95 - map.ave)/1.64  # From Table B.1 in Wilks (1995) p. 434

  lon.z <- attr(map.q95,'longitudes'); nx <- length(lon.z)
  lat.z <- attr(map.q95,'latitudes'); ny <- length(lat.z)

  map.P <- map.std*0; map.P.ref <- map.P
  
  for (i in 1:nx) {
    for (j in 1:ny) {
      map.P[i,j] <- 100*pnorm(thresh,mean=map.ave[i,j],sd=map.std[i,j])
      map.P.ref[i,j] <- 100*pnorm(thresh,mean=map.ave.ref[i,j],sd=map.std.ref[i,j])
    }
  }
  
  if (plot) {
    x11()
    par(mfrow=c(1,3))
    
    image(lon.z,lat.z,map.P.ref,
          main=paste("Prob( X <",thresh,") ", seasons[season],ref),
        col=heat.colors(21),xlab="",ylab="",zlim=c(0,150))
    addland()
    contour(lon.z,lat.z,map.P.ref,add=TRUE,levels=seq(0,100,by=5))
    
    image(lon.z,lat.z,map.P,main=paste("Prob( X <",thresh,") ",seasons[season],year),
        col=heat.colors(21),xlab="",ylab="",zlim=c(0,150))
    addland()
    contour(lon.z,lat.z,map.P,add=TRUE,levels=seq(0,100,by=10))

    z1 <- max(abs(c(map.P-map.P.ref)),na.rm=TRUE)
    image(lon.z,lat.z,map.P-map.P.ref,main=paste("P difference",year,"-",ref),
        col=heat.colors(21),xlab="",ylab="",zlim=c(-z1,z1))
    addland()
    contour(lon.z,lat.z,map.P-map.P.ref,add=TRUE)
    dev.copy2eps(file=paste("mapESDprobs_",seasons[season],
                   "-",year,"-",ref,".eps",sep=""))
  }
#  print("finished")
  attr(map.P,'longitudes') <- lon.z
  attr(map.P,'latitudes') <- lat.z
  attr(map.P,'season') <- seasons[season]
  attr(map.P,'what') <- "Probability"
  attr(map.P,'year') <- year
  invisible(map.P) 
}

ESDinGoogle <- function(browser = "firefox", url="http://eklima.met.no/metno/esd/esd.google.earthTemp.kmz") {
    print("View the ESD results in GoogleEarth - you may have to install a plug-in, or save the page as a KML-file and then open this with GoogleEarth.")
    system(paste(browser, url))

}

ESDdetails <- function(browser = "firefox", url="http://met.no/Forskning/Publikasjoner/") {
  print("Relevant reports are: met.no Note 03/2009 and met.no Note 15/2009")
  print("The report can be accessed through browser/PDF-viewer")
    system(paste(browser, url))
}



ESDreference <- function(browser = "firefox", url="http://www.agu.org/pubs/crossref/2005/2005GL023401.shtml") {
    print("Opens the compendium in browser/PDF-viewer")
    system(paste(browser, url))

}

