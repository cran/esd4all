#get.num <- function(txt,index=1) {
#  #print(txt)
#  require(met.no.REB) # http://noserc.met.no/grtools/reb.html
#  txt <- replace.char("<",txt," ")
#  spc <- instring(" ",txt)
#  spc[c(diff(spc)==1,FALSE)] <- NA
#  spc <- c(spc[is.finite(spc)],nchar(txt))
#  ele <- rep(NA,length(spc))
#  for (ii in 1:(length(spc)-1)) {
#    ele[ii] <- substr(txt,spc[ii],spc[ii+1])
#  }
#  #print(ele)
#  as.numeric(ele[index])
#}
#
#get.c.anova <- function(fname="~/latex/publications/paper44/gridESD.out") {
#  op <- readLines(fname)
#  anova.head <- grep("Estimate Std. Error t value",op)
#  N <- length(anova.head)
#  c.est <- matrix(rep(NA,8*N),N,8*4)
#  dim(c.est) <- c(N,8,4)
#  for (i in 1:N) {
#    for (ii in 1:4) {
#      c.est[i,1,ii] <- get.num(op[anova.head[i]+1],ii)
#      c.est[i,2,ii] <- get.num(op[anova.head[i]+2],ii)
#      c.est[i,3,ii] <- get.num(op[anova.head[i]+3],ii)
#      c.est[i,4,ii] <- get.num(op[anova.head[i]+4],ii)
#      c.est[i,5,ii] <- get.num(op[anova.head[i]+5],ii)
#      c.est[i,6,ii] <- get.num(op[anova.head[i]+6],ii)
#      c.est[i,7,ii] <- get.num(op[anova.head[i]+7],ii)
#      c.est[i,8,ii] <- get.num(op[anova.head[i]+8],ii)
#    }
#  }
#  colnames(c.est) <- c("intercept","Z","d","Y","X",
#                       "ln(d)","sqrt(d)","sqrt(Z)")
#  c.est
#}

plotCoefs <- function(a1,a2,a3,types,year,it,icoef=2,pch=c(19,21,22),
                      col=c("blue","darkgreen","red","grey40"),x.scale=1000) {
data(esdsummary)
par(col.axis="white")
type <- paste("c",icoef,".",sep="")
iv <- grep(type,types)
#print(iv); print(type); print(types)
if (icoef>1) scaling <- c(0.3,1,0.1,0.1,3*log(x.scale),sqrt(x.scale),10) else
             scaling <- rep(1,7)
Scaling <- rep(scaling,length(iv)); dim(Scaling) <- c(7,length(iv)); Scaling <- t(Scaling)
#print(dim(Scaling)); print(dim(a1[iv,2:8,1])); print(Scaling)
ylim <- range(a1[iv,2:8,1]/Scaling-a1[iv,2:8,2]/Scaling,
              a2[iv,2:8,1]/Scaling-a2[iv,2:8,2]/Scaling,
              a3[iv,2:8,1]/Scaling-a3[iv,2:8,2]/Scaling,
              a1[iv,2:8,1]/Scaling+a1[iv,2:8,2]/Scaling,
              a2[iv,2:8,1]/Scaling+a2[iv,2:8,2]/Scaling,
              a3[iv,2:8,1]/Scaling+a3[iv,2:8,2]/Scaling)*
                esdsummary$X[it]^(icoef - 1)
ylim <- ylim + c(0,5)
plot(c(0,8),ylim,type="n",
     ylab=paste("arbitrary scaled associated with c_",icoef-1,"t^",icoef-1,sep=""),xlab="parameter",
     main=paste("Contribution",year,"polynomial coef",icoef-1),
     sub="Geographical regression")
par(col.axis="black")
axis(side=1,at=1:7,labels=colnames(a1)[2:8])
axis(side=2)
grid()
legend(1,max(ylim),
       c(attr(a1,"name"),attr(a2,"name"),attr(a3,"name")),
       col="black",pch=c(19,21,22),bg="grey95")
legend(3,max(ylim),
       c("DJF","MAM","JJA","SON"),col=col,
       lty=1,lwd=4,bg="grey95")

ii <- 1
  for (season in 1:4) {
    for (statistic in 1:3) {
      #print(type[ii])
      for (i in 1:6) {
        y1 <- a1[ii,2:8,1]
        yu1 <- y1 + a1[ii,2:8,2]
        yl1 <- y1 - a1[ii,2:8,2]
        y2 <- a2[ii,2:8,1]
        yu2 <- y2 + a2[ii,2:8,2]
        yl2 <- y2 - a2[ii,2:8,2]
        y3 <- a3[ii,2:8,1]
        yu3 <- y3 + a3[ii,2:8,2]
        yl3 <- y3 - a3[ii,2:8,2]

          y1 <- y1 * esdsummary$X[it]^(i - 1)/scaling
          yu1 <- yu1 * esdsummary$X[it]^(i - 1)/scaling
          yl1 <- yl1 * esdsummary$X[it]^(i - 1)/scaling
          y2 <- y2 * esdsummary$X[it]^(i - 1)/scaling
          yu2 <- yu2 * esdsummary$X[it]^(i - 1)/scaling
          yl2 <- yl2 * esdsummary$X[it]^(i - 1)/scaling
          y3 <- y3 * esdsummary$X[it]^(i - 1)/scaling
          yu3 <- yu3 * esdsummary$X[it]^(i - 1)/scaling
          yl3 <- yl3 * esdsummary$X[it]^(i - 1)/scaling

        if (i==icoef) {
          points(1:7 - 0.25,y1,pch=19,col=col[season])
          for (iii in 1:7) lines(rep(iii,2)-0.25,c(yu1[iii],yl1[iii]),
                                 lty=1,col=col[season])
          points(1:7,y2,pch=21,col=col[season])
          for (iii in 1:7) lines(rep(iii,2),c(yu2[iii],yl2[iii]),
                                 lty=1,col=col[season])
          points(1:7 + 0.25,y3,pch=22,col=col[season])
          for (iii in 1:7) lines(rep(iii,2)+0.25,c(yu3[iii],yl3[iii]),
                                 lty=1,col=col[season])
          #print(y1); print(yl1); print(yu1); 
        }
        ii <- ii+1
      }
    }
  }
}

library(esd4all)
data(gridded.c)
R2 <- attributes(gridded.c)$GRM.R2
dim(R2) <- c(6,12); R2 <- t(R2)
colnames(R2) <- paste("c",0:5,sep="_")
rownames(R2) <- paste(rep(c("mean","q05","q95"),4),
                      c(rep("DJF",3),rep("MAM",3),rep("JJA",3),rep("SON",3)))
print(R2)

types <- attributes(gridded.c)$type
rm("gridded.c"); gc(reset=TRUE)

data(gridded.africa.c)
R2 <- attributes(gridded.c)$GRM.R2
dim(R2) <- c(6,12); R2 <- t(R2)
colnames(R2) <- paste("c",0:5,sep="_")
rownames(R2) <- paste(rep(c("mean","q05","q95"),4),
                      rep(c("DJF","MAM","JJA","SON"),3))
print(R2)
rm("gridded.c"); gc(reset=TRUE)
data(gridded.ealat.c)
R2 <- attributes(gridded.c)$GRM.R2
dim(R2) <- c(6,12); R2 <- t(R2)
colnames(R2) <- paste("c",0:5,sep="_")
rownames(R2) <- paste(rep(c("mean","q05","q95"),4),
                      rep(c("DJF","MAM","JJA","SON"),3))
print(R2)
rm("gridded.c"); gc(reset=TRUE)

data(esdsummary)
it <- length(esdsummary$X)
year <- esdsummary$X[it] + esdsummary$X.0


#get.c.anova() -> grm.coef
#get.c.anova("~/latex/publications/paper44/gridESD-Africa.out") ->
#  grm.coef.africa
#get.c.anova("~/latex/publications/paper44/gridESD-EALAT.out") ->
#  grm.coef.ealat
#attr(grm.coef,"name") <- "Europe"
#attr(grm.coef.africa,"name") <- "Africa"
#attr(grm.coef.ealat,"name") <- "NW Russia"
#save(file="esd4all/data/grm.coef.rda",grm.coef)
#save(file="esd4all/data/grm.coef.africa.rda",grm.coef.africa)
#save(file="esd4all/data/grm.coef.ealat.rda",grm.coef.ealat)

data(grm.coef)
data(grm.coef.africa)
data(grm.coef.ealat)

plotCoefs(grm.coef,grm.coef.africa,grm.coef.ealat,types,year=year,it=it,icoef=1,x.scale)
dev.copy2eps(file="paper44-grm-c_0.eps")

x11()
plotCoefs(grm.coef,grm.coef.africa,grm.coef.ealat,types,year=year,it=it,icoef=2,x.scale)
dev.copy2eps(file="paper44-grm-c_1.eps")

x11()
plotCoefs(grm.coef,grm.coef.africa,grm.coef.ealat,types,year=year,it=it,icoef=3,x.scale)
dev.copy2eps(file="paper44-grm-c_2.eps")


x11()
plotCoefs(grm.coef,grm.coef.africa,grm.coef.ealat,types,year=year,it=it,icoef=4,x.scale)
dev.copy2eps(file="paper44-grm-c_3.eps")
  

x11()
plotCoefs(grm.coef,grm.coef.africa,grm.coef.ealat,types,year=year,it=it,icoef=5,x.scale)
dev.copy2eps(file="paper44-grm-c_4.eps")


x11()
plotCoefs(grm.coef,grm.coef.africa,grm.coef.ealat,types,year=year,it=it,icoef=6,x.scale)
dev.copy2eps(file="paper44-grm-c_5.eps")
  





