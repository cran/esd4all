reduce.rda.size <- function(get.data="data(gridded.c,envir=environment())",
                            reduce.res=TRUE,Nx=100,Ny=100) {
  eval(parse(text=get.data))
  nx <- Nx; ny <- Ny
  d <- dim(gridded.c)
  if (!reduce.res) {
    dim(gridded.c) <- c(d[1]*d[2],d[3])
    land <- is.finite(gridded.c[,1])
    gridded.C <- rep(NA,sum(land)*d[3])
    dim(gridded.C) <- c(sum(land),d[3])
    print(dim(gridded.C)); print(dim(gridded.c[land,]))
    for (i in 1:d[3]) gridded.C[,i] <- gridded.c[land,i]
    attr(gridded.C,"longitudes") <- attr(gridded.c,"longitudes")
    attr(gridded.C,"latitudes") <- attr(gridded.c,"latitudes")
    attr(gridded.C,"residual_gridding_method") <- attr(gridded.c,"residual_gridding_method")
    attr(gridded.C,"type") <- attr(gridded.c,"type")
    attr(gridded.C,"dimensions") <- d
    attr(gridded.C,"landmask") <- land
    print("Save reduced data size in 'gridded.c.Rdata'")
    save(file="gridded.c.Rdata",gridded.C,compress=TRUE)
  } else {
    lon <- as.numeric(attr(gridded.c,"longitudes"))
    lat <- as.numeric(attr(gridded.c,"latitudes"))
    lonxy <- rep(lon,length(lat))
    latxy <- sort(rep(lat,length(lon)))
    lono <- seq(min(lon),max(lon),length=nx)
    lato <- seq(min(lat),max(lat),length=ny)
    gridded.C <- rep(NA,nx*ny*d[3])
    dim(gridded.C) <- c(nx,ny,d[3])
    gridded.c[is.na(gridded.c)] <-  -9999
    print(dim(gridded.c[,,1])); print(length(lon)); print(length(lat))
    print(length(lono)); print(length(lato))
    for (i in 1:d[3]) {
      map <- interp(lonxy,latxy,gridded.c[,,i],lono,lato)$z
      gridded.C[,,i] <- map
    }
    
    gridded.C[gridded.C < -99] <- NA
    attr(gridded.C,"longitudes") <- seq(min(lon),max(lon),length=nx)
    attr(gridded.C,"latitudes") <- seq(min(lat),max(lat),length=ny)
    attr(gridded.C,"residual_gridding_method") <- attr(gridded.c,"residual_gridding_method")
    attr(gridded.C,"type") <- attr(gridded.c,"type")
    attr(gridded.C,"post-process") <- "reduction of resolution (reduce.rda.size)"
    gridded.c <- gridded.C; rm(gridded.C)
    save(file="gridded.c.rda",gridded.c,compress=TRUE)
}
}

rda2cdf <- function(get.data="data(gridded.c,envir=environment())") {
  eval(parse(text=get.data))
  lon.z <- attr(gridded.c,'longitudes')
  lat.z <- attr(gridded.c,'latitudes')
  gridded.c[!is.finite(gridded.c)] <- -99
  dimlon <- dim.def.ncdf( "lon", "degrees_east", lon.z)
  dimlat <- dim.def.ncdf( "lat", "degrees_north", lat.z)
  dim(gridded.c) -> d.c
  dimn <- dim.def.ncdf( "coefnr", "none", 1:d.c[3])
  vardef <- var.def.ncdf(name="coef",units="none",dim=list(dimlon,dimlat,dimn), missval=-99,
                       longname="coefficients for 5th-order polynomials", prec="short")
  metadef <- var.def.ncdf(name="meta",units="none",dim=dimn, missval="NA",
                        longname="coefficients for 5th-order polynomials", prec="char")
  print("Create netCDF-file")
  ncnew <- create.ncdf("Europe_E-SDS_t2m-trend_map.nc", list(vardef,metadef))
  put.var.ncdf( ncnew, vardef, round(gridded.c*10) )
  print("Add metadata")
  attr(gridded.c,'type') -> coef.tag
  put.var.ncdf( ncnew, metadef, coef.tag)
  print("add attributes")
  att.put.ncdf( ncnew, vardef, 'scale_factor', 0.1, prec="double")
  att.put.ncdf( ncnew, 0, 'scenario', 'IPCC SRES A1b')
  att.put.ncdf( ncnew, 0, 'author', 'Rasmus Benestad')
  att.put.ncdf( ncnew, 0, 'title', 'Climate change scenarios from multi-model IPCC AR4 climate simulations')
  att.put.ncdf( ncnew, 0, 'R-package', 'esd4all')
  att.put.ncdf( ncnew, 0, 'date', date())
  att.put.ncdf( ncnew, 0, 'key_words','statistical-empirical, downscaling, multi-model ensemble, gridding')
  att.put.ncdf( ncnew, 0, 'ensemble','CMIP3')
  grid.res.method<-'bi-linear residual interpolation'
  att.put.ncdf( ncnew, 0, 'residual_gridding_method',grid.res.method)
  att.put.ncdf( ncnew, 0, 'GCM_source','PCMDI - Earth System Grid: URL https://esg.llnl.gov:8443/index.jsp')
  close.ncdf(ncnew)
}

figures <- function(get.data="data(gridded.c,envir=environment())",
                    season.1=3,season.2=1,year=2050,thresh=0,what="q95") {

  seasons <- c("DJF","MAM","JJA","SON")
  print("Quantiles")
  mapESDquants(get.data1=get.data,season=season.1,year=year,what=what) -> map.q95
  print(dim(map.q95))
  dev2bitmap(file="JJA-q95.jpg",type="jpeg",res=200)
  longname=paste(what,"for",seasons[season.1],"mean T(2m)  in year",year)
  
  lon.z <- attr(map.q95,'longitudes')
  lat.z <- attr(map.q95,'latitudes')
  map.q95[!is.finite(map.q95)] <- -99
  dimlon <- dim.def.ncdf( "lon", "degrees_east", lon.z)
  dimlat <- dim.def.ncdf( "lat", "degrees_north", lat.z)
  vardef <- var.def.ncdf(name="q95",units="none",dim=list(dimlon,dimlat),
                         missval=-99,longname=longname, prec="short")
  ncnew <- create.ncdf("Gridded_ESD-q95map.nc", vardef)
  put.var.ncdf( ncnew, vardef, round(map.q95*10) )

  att.put.ncdf( ncnew, vardef, 'scale_factor', 0.1, prec="double")
  att.put.ncdf( ncnew, vardef, 'missing_value', -990.0, prec="double")
  att.put.ncdf( ncnew, 0, 'scenario', 'IPCC SRES A1b')
  att.put.ncdf( ncnew, 0, 'author', 'Rasmus Benestad')
  att.put.ncdf( ncnew, 0, 'title', 'Climate change scenarios from multi-model IPCC AR4 climate simulations')
  att.put.ncdf( ncnew, 0, 'R-script', 'figures.gridESD.R')
  att.put.ncdf( ncnew, 0, 'date', date())
  att.put.ncdf( ncnew, 0, 'key_words','statistical-empirical, downscaling, multi-model ensemble, gridding')
  att.put.ncdf( ncnew, 0, 'ensemble','CMIP3')
  grid.res.method <- "bi-polar interpolation of residuals"
  att.put.ncdf( ncnew, 0, 'residual_gridding_method',grid.res.method)
  att.put.ncdf( ncnew, 0, 'GCM_source','PCMDI - Earth System Grid: URL https://esg.llnl.gov:8443/index.jsp')
  close.ncdf(ncnew)

  print("Probabilities")
  
#  mapESDprobs(get.data="load('gridded.c.large.rda')") -> map.p0
  mapESDprobs(get.data=get.data,season=season.2,year=year,thresh=thresh) -> map.p0
  dev2bitmap(file="DJF-P-T-lt-0.jpg",type="jpeg",res=200) 
  longname=paste("probability for",seasons[season.2],
    "mean T(2m) < ",thresh,"in year",year)
  
  map.p0[!is.finite(map.p0)] <- -90
  print(summary(c(map.p0)))
  var2def <- var.def.ncdf(name="Pr",units="none",dim=list(dimlon,dimlat),
                          missval=-99,longname=longname,prec="short")
  ncnew2 <- create.ncdf("Gridded_ESD-p0map.nc", var2def)
  put.var.ncdf( ncnew2, var2def, round(map.p0) )

  att.put.ncdf( ncnew2, var2def, 'scale_factor', 1, prec="double")
  att.put.ncdf( ncnew2, var2def, 'missing_value', -90.0, prec="double")
  att.put.ncdf( ncnew2, 0, 'scenario', 'IPCC SRES A1b')
  att.put.ncdf( ncnew2, 0, 'author', 'Rasmus Benestad')
  att.put.ncdf( ncnew2, 0, 'title', 'Climate change scenarios from multi-model IPCC AR4 climate simulations')
  att.put.ncdf( ncnew2, 0, 'R-script', 'figures.gridESD.R')
  att.put.ncdf( ncnew2, 0, 'date', date())
  att.put.ncdf( ncnew2, 0, 'key_words','statistical-empirical, downscaling, multi-model ensemble, gridding')
  att.put.ncdf( ncnew2, 0, 'ensemble','CMIP3')
  grid.res.method <- "bi-polar interpolation of residuals"
  att.put.ncdf( ncnew, 0, 'residual_gridding_method',grid.res.method)
  att.put.ncdf( ncnew, 0, 'GCM_source','PCMDI - Earth System Grid: URL https://esg.llnl.gov:8443/index.jsp')
  close.ncdf(ncnew)
}
