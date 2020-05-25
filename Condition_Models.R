#save.image("Condition_Models.Rdata")
#load("Condition_Models.Rdata")

d <- read.csv("elo cond repro.csv", header = TRUE, na.strings = "")

### repro state 
d$preg <- ifelse(d$newstate=="preg" , 1 , 0 )
d$lact <- ifelse(d$newstate=="lact" , 1 , 0 )
d$flat <- ifelse(d$newstate=="flat" , 1 , 0 )
d$swol <- ifelse(d$newstate=="swol" , 1 , 0 )

### groups
d$phg <- ifelse(d$group=="PHG" , 1 , 0 )
d$enk <- ifelse(d$group=="ENK" , 1 , 0 )
d$ynt <- ifelse(d$group=="YNT" , 1 , 0 )
d$chk <- ifelse(d$group=="CHK" , 1 , 0 )

#make character or factors into integers 
d$id <- as.factor(d$id)
d$id <- as.integer(as.factor(d$id))

dd <- d[complete.cases(d$id,
                       d$mean_elo,
                       d$phg,
                       d$enk,
                       d$ynt,
                       d$infage,
                       d$reverse_cond), ]

dd$s_focelo <- (dd$mean_elo - mean(dd$mean_elo))/sd(dd$mean_elo)
dd$s_foccondRev <- (dd$reverse_cond - mean(dd$reverse_cond))/sd(dd$reverse_cond)
dd$s_infage <- (dd$infage - mean(dd$infage))/sd(dd$infage)
dd$s_grpsize <- (dd$grpsize - mean(dd$grpsize))/sd(dd$grpsize)
dd<-dd[!(dd$chk==1),]

#redo bc NAs were removed
dd$id=as.factor(dd$id)
dd$id=as.integer(as.factor(dd$id))

condition_Model <- map2stan(
  alist(
    
    foccond ~ dnorm(mu , sigma ),
    
    mu <- ap + gp_id[id] + bp_enk*enk + bp_ynt*ynt + bp_elo1*focelo + 
      bp_preg*preg + bp_flat*flat + bp_swol*swol +
      bp_grpsize*grpsize,
    
    c(gp_id)[id] ~ dnorm(0,sigma_focal) ,
    c(ap) ~ dnorm(0,2) ,
    c(bp_enk,bp_ynt,bp_elo1,bp_grpsize,
      bp_preg,bp_flat,bp_swol) ~ dnorm(0,2) ,
    sigma_focal  ~ dexp(1)  ,
    sigma ~ dexp(1)
    
  ),
  data=list(
    foccond = dd$reverse_cond,
    id=dd$id,
    focelo = dd$s_focelo,
    phg = dd$phg,
    enk = dd$enk,
    ynt = dd$ynt,
    preg = dd$preg,
    lact = dd$lact,
    swol = dd$swol,
    flat = dd$flat,
    grpsize = dd$s_grpsize
  ),
  
  cores=3 , chains=2 , warmup=3000, iter=6000, WAIC=TRUE, types=list(adapt.delta=0.99)
)


#limit to pregnant
ddpreg<-dd[!(dd$preg==0),]
ddpreg$s_infage <- (ddpreg$infage - mean(ddpreg$infage))/sd(ddpreg$infage)

ddpreg$id=as.factor(ddpreg$id)
ddpreg$id=as.integer(as.factor(ddpreg$id))

condition_Pregnant <- map2stan(
  alist(
    
    foccond ~ dnorm(mu , sigma ),
    
    mu <- ap + gp_id[id] + bp_enk*enk + bp_ynt*ynt + bp_elo1*focelo + 
      bp_infage*infage + bp_grpsize*grpsize,
    
    c(gp_id)[id] ~ dnorm(0,sigma_focal) ,
    c(ap) ~ dnorm(0,3) ,
    c(bp_enk,bp_ynt,bp_elo1,bp_infage,bp_grpsize) ~ dnorm(0,2) ,
    sigma_focal  ~ dexp(1)  ,
    sigma ~ dexp(1)
    
  ),
  data=list(
    foccond = ddpreg$reverse_cond,
    id=ddpreg$id,
    focelo = ddpreg$s_focelo,
    phg = ddpreg$phg,
    enk = ddpreg$enk,
    ynt = ddpreg$ynt,
    infage = ddpreg$s_infage,
    grpsize = ddpreg$s_grpsize
  ),
  
  cores=3 , chains=3 , warmup=3000, iter=6000, WAIC=TRUE, types=list(adapt.delta=0.99)
)

#limit to lactating
ddlact<-dd[!(dd$lact==0),]
ddlact$s_infage <- (ddlact$infage - mean(ddlact$infage))/sd(ddlact$infage)

ddlact$id=as.factor(ddlact$id)
ddlact$id=as.integer(as.factor(ddlact$id))

condition_Lactating <- map2stan(
  alist(
    
    foccond ~ dnorm(mu , sigma ),
    
    mu <- ap + gp_id[id] + bp_enk*enk + bp_ynt*ynt + bp_elo1*focelo + 
      bp_infage*infage + bp_grpsize*grpsize,
    
    c(gp_id)[id] ~ dnorm(0,sigma_focal) ,
    c(ap) ~ dnorm(0,3) ,
    c(bp_enk,bp_ynt,bp_elo1,bp_infage,bp_grpsize) ~ dnorm(0,2) ,
    sigma_focal  ~ dexp(1)  ,
    sigma ~ dexp(1)
    
  ),
  data=list(
    foccond = ddlact$reverse_cond,
    id=ddlact$id,
    focelo = ddlact$s_focelo,
    phg = ddlact$phg,
    enk = ddlact$enk,
    ynt = ddlact$ynt,
    infage = ddlact$s_infage,
    grpsize = ddlact$s_grpsize
  ),
  
  cores=3 , chains=3 , warmup=3000, iter=6000, WAIC=TRUE, types=list(adapt.delta=0.99)
)




#Full model: elo
a_foc_z <- matrix(0,1000,length(unique(dd$id)))
elo.seq=seq(min(dd$s_focelo),max(dd$s_focelo),length=1000)

d.pred_elo <- list(
  id=rep(1,length(elo.seq)),
  grpsize=rep(mean(dd$s_grpsize),length(elo.seq)),
  enk=rep(mean(dd$enk),length(elo.seq)),
  ynt=rep(mean(dd$ynt),length(elo.seq)),
  preg=rep(mean(dd$preg),length(elo.seq)),
  flat=rep(mean(dd$flat),length(elo.seq)),
  swol=rep(mean(dd$swol),length(elo.seq)),
  focelo=elo.seq
)

pred <- link(condition_Model, n=1000 , data=d.pred_elo, replace=
               list(ap_id=a_foc_z , al_id=a_foc_z, WAIC=TRUE))
pred.median <- apply( pred , 2 , median )
pred.HPDI <- apply( pred , 2 , HPDI )

par(mfrow = c(1, 1), cex=1.1, mar=c(0,0,0,0), oma=c(4,3,1.1,4))
plot( reverse_cond ~ s_focelo, data=dd , col=alpha("#33CCFF",0.1),pch=19 ,xaxt='n',xlim = c(-1.5,2),xlab=NA,yaxt='n',ylab=NA,ylim = c(0,5),cex=1)

pred.lines = pred[sample(nrow(pred),100,replace=F),]
for (i in 1:100){
  lines( elo.seq , pred.lines[i,] ,lwd=2,col=alpha("#33CCFF",0.05))
}
lines( elo.seq , pred.median , lwd=1,col="black")
lines( elo.seq , pred.HPDI[1,],lty=2,lwd=1)
lines( elo.seq , pred.HPDI[2,],lty=2,lwd=1)
axis(1, at = seq(from=-2 , to=5, by = 1) ,tck=-0.02 , labels=T )
axis(2, at = seq(from=-4 , to=5, by = 1) ,tck=-0.02 , labels=T )

mtext(side=1,line=2,text="elo rank standardized",cex=1.5)
mtext(side=2,line=2,text="body condition",cex=1.5)

#Full model: group size
a_foc_z <- matrix(0,1000,length(unique(dd$id)))
grp.seq=seq(min(dd$s_grpsize),max(dd$s_grpsize),length=1000)

d.pred_grp <- list(
  id=rep(1,length(grp.seq)),
  enk=rep(mean(dd$enk),length(grp.seq)),
  ynt=rep(mean(dd$ynt),length(grp.seq)),
  focelo=rep(mean(dd$s_focelo),length(grp.seq)),
  preg=rep(mean(dd$preg),length(grp.seq)),
  flat=rep(mean(dd$flat),length(grp.seq)),
  swol=rep(mean(dd$swol),length(grp.seq)),
  grpsize=grp.seq
)

pred <- link(condition_Model, n=1000 , data=d.pred_grp, replace=
               list(ap_id=a_foc_z , al_id=a_foc_z, WAIC=TRUE))
pred.median <- apply( pred , 2 , median )
pred.HPDI <- apply( pred , 2 , HPDI )

par(mfrow = c(1, 1), cex=1.1, mar=c(0,0,0,0), oma=c(4,3,1.1,4))
plot( reverse_cond ~ s_grpsize, data=dd , col=alpha("#33CCFF",0.1),pch=19 ,xaxt='n',xlim = c(-4.7,3),xlab=NA,yaxt='n',ylab=NA,ylim = c(0,5),cex=1)

pred.lines = pred[sample(nrow(pred),100,replace=F),]
for (i in 1:100){
  lines( grp.seq , pred.lines[i,] ,lwd=2,col=alpha("#33CCFF",0.05))
}
lines( grp.seq , pred.median , lwd=1,col="black")
lines( grp.seq , pred.HPDI[1,],lty=2,lwd=1)
lines( grp.seq , pred.HPDI[2,],lty=2,lwd=1)
axis(1, at = seq(from=-5 , to=5, by = 1) ,tck=-0.02 , labels=T )
axis(2, at = seq(from=-4 , to=5, by = 1) ,tck=-0.02 , labels=T )

mtext(side=1,line=2,text="group size standardized",cex=1.5)
mtext(side=2,line=2,text="body condition",cex=1.5)


#Full model: repro state

a_foc_z <- matrix(0,1000,length(unique(dd$id)))

d.pred.flat <- list(
  id=c(1,1),
  enk=rep(mean(dd$enk),2),
  ynt=rep(mean(dd$ynt),2),
  focelo=rep(mean(dd$s_focelo),2),
  grpsize=rep(mean(dd$s_grpsize),2),
  flat=c(1,1),
  preg=c(0,0),
  swol=c(0,0)
)

d.pred.lact <- list(
  id=c(1,1),
  enk=rep(mean(dd$enk),2),
  ynt=rep(mean(dd$ynt),2),
  focelo=rep(mean(dd$s_focelo),2),
  grpsize=rep(mean(dd$s_grpsize),2),
  flat=c(0,0),
  preg=c(0,0),
  swol=c(0,0)
)

d.pred.preg <- list(
  id=c(1,1),
  enk=rep(mean(dd$enk),2),
  ynt=rep(mean(dd$ynt),2),
  focelo=rep(mean(dd$s_focelo),2),
  grpsize=rep(mean(dd$s_grpsize),2),
  flat=c(0,0),
  preg=c(1,1),
  swol=c(0,0)
)

d.pred.swol <- list(
  id=c(1,1),
  enk=rep(mean(dd$enk),2),
  ynt=rep(mean(dd$ynt),2),
  focelo=rep(mean(dd$s_focelo),2),
  grpsize=rep(mean(dd$s_grpsize),2),
  flat=c(0,0),
  preg=c(0,0),
  swol=c(1,1)
)

pred.flat <- link(condition_Model, n=1000 , data=d.pred.flat, replace=
                    list(ap_id=a_foc_z , al_id=a_foc_z), WAIC=TRUE)
pred.flat.median <- apply( pred.flat , 2 , median )
pred.flat.HPDI <- apply( pred.flat , 2 , HPDI )
median(pred.flat[,1])
HPDI(pred.flat[,1])

pred.lact <- link(condition_Model, n=1000 , data=d.pred.lact, replace=
                    list(ap_id=a_foc_z , al_id=a_foc_z), WAIC=TRUE)
pred.lact.median <- apply( pred.lact , 2 , median )
pred.lact.HPDI <- apply( pred.lact , 2 , HPDI )
median(pred.lact[,1])
HPDI(pred.lact[,1])

pred.preg <- link(condition_Model, n=1000 , data=d.pred.preg, replace=
                    list(ap_id=a_foc_z , al_id=a_foc_z), WAIC=TRUE)
pred.preg.median <- apply( pred.preg , 2 , median )
pred.preg.HPDI <- apply( pred.preg , 2 , HPDI )
median(pred.preg[,1])
HPDI(pred.preg[,1])

pred.swol <- link(condition_Model, n=1000 , data=d.pred.swol, replace=
                    list(ap_id=a_foc_z , al_id=a_foc_z), WAIC=TRUE)
pred.swol.median <- apply( pred.swol , 2 , median )
pred.swol.HPDI <- apply( pred.swol , 2 , HPDI )
median(pred.swol[,1])
HPDI(pred.swol[,1])

par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(1,1,1,1))

dens(pred.lact[,1], xlim=c(1,5) , ylim=c(-.035,3) , xlab="" , ylab="" , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE , cex.lab=1.5 , cex.axis=0.5 , adj=0.1)
ll <- dd$reverse_cond[dd$lact==1]
points(ll, rep(-.025,length(ll)), pch="|" , col=col.alpha("#33CCFF", alpha=0.5) , cex=0.5 )
shade( density(pred.lact[,1]) , lim= as.vector(HPDI(pred.lact[,1], prob=0.99999)) , col = col.alpha("#33CCFF", 0.5))
shade( density(pred.preg[,1]) , lim= as.vector(HPDI(pred.preg[,1], prob=0.99999)) , col = col.alpha("orange1", 0.5))
shade( density(pred.swol[,1]) , lim= as.vector(HPDI(pred.swol[,1], prob=0.99999)) , col = col.alpha("green", 0.5))
shade( density(pred.flat[,1]) , lim= as.vector(HPDI(pred.flat[,1], prob=0.99999)) , col = col.alpha("purple", 0.5))
ll <- dd$reverse_cond[dd$preg==1]
points(ll, rep(-.025,length(ll)), pch="|" , col=col.alpha("orange1", alpha=0.5) , cex=0.5 )
ll <- dd$reverse_cond[dd$swol==1]
points(ll, rep(-.025,length(ll)), pch="|" , col=col.alpha("green", alpha=0.5) , cex=0.5 )
ll <- dd$reverse_cond[dd$lat==1]
points(ll, rep(-.025,length(ll)), pch="|" , col=col.alpha("purple", alpha=0.5) , cex=0.5 )

abline(v=median(pred.lact[,1]) , lty=1)
abline(v=median(pred.preg[,1]) , lty=2)
abline(v=median(pred.swol[,1]) , lty=3)
abline(v=median(pred.flat[,1]) , lty=4)

legend(3.8,2.5, legend = c("", "", "", ""),
       col=c(1,1)  , lty=c(1,2,3,4),
       lw=1 , cex=.85, bty="n")

legend(4,2.5,, legend = c("lactating", "pregnant", "swollen", "flat"),
       col=c(col.alpha("#33CCFF", 0.5) , col.alpha("orange1", 0.5) ,  col.alpha("green", 0.5) , col.alpha("purple", 0.5)) , pch=c(15,17),
       pt.cex=2 , cex=.85, bty="n")

mtext("body condition" , side=1 , line=1 , cex=1.5)
#axis(1, at = seq(from=0 , to=1, by = .1) ,tck=-0.1 , labels=T )
axis(1, at = seq(from=0 , to=5, by = 1) ,tck=-0.01 , labels=T )



#Pregnant & Lactating Model Plots

a_foc_z <- matrix(0,1000,length(unique(ddpreg$id)))
timeP.seq=seq(min(ddpreg$s_infage),max(ddpreg$s_infage),length=1000)

d.pred_timeP <- list(
  id=rep(1,length(timeP.seq)),
  enk=rep(mean(ddpreg$enk),length(timeP.seq)),
  ynt=rep(mean(ddpreg$ynt),length(timeP.seq)),
  focelo=rep(mean(ddpreg$s_focelo),length(timeP.seq)),
  grpsize=rep(mean(ddpreg$s_grpsize),length(timeP.seq)),
  infage=timeP.seq
)

predP <- link(condition_Pregnant, n=1000 , data=d.pred_timeP, replace=
               list(ap_id=a_foc_z , al_id=a_foc_z, WAIC=TRUE))
predP.median <- apply( predP , 2 , median )
predP.HPDI <- apply( predP , 2 , HPDI )

par(mfrow = c(1, 1), cex=1.1, mar=c(0,0,0,0), oma=c(4,3,1.1,4))
plot( reverse_cond ~ s_infage, data=ddpreg , col=alpha("#33CCFF",0.1),pch=19 ,xaxt='n',xlim = c(-2,1.7),xlab=NA,yaxt='n',ylab=NA,ylim = c(0,5),cex=1)

predP.lines = predP[sample(nrow(predP),100,replace=F),]
for (i in 1:100){
  lines( timeP.seq , predP.lines[i,] ,lwd=2,col=alpha("#33CCFF",0.05))
}
lines( timeP.seq , predP.median , lwd=1,col="black")
lines( timeP.seq , predP.HPDI[1,],lty=2,lwd=1)
lines( timeP.seq , predP.HPDI[2,],lty=2,lwd=1)
axis(1, at = seq(from=-5 , to=5, by = 1) ,tck=-0.02 , labels=T )
axis(2, at = seq(from=-4 , to=5, by = 1) ,tck=-0.02 , labels=T )

mtext(side=1,line=2,text="time across pregnancy",cex=1.5)
mtext(side=2,line=2,text="body condition",cex=1.5)

a_foc_z <- matrix(0,1000,length(unique(ddlact$id)))
timeL.seq=seq(min(ddlact$s_infage),max(ddlact$s_infage),length=1000)

d.pred_timeL <- list(
  id=rep(1,length(timeL.seq)),
  enk=rep(mean(ddlact$enk),length(timeL.seq)),
  ynt=rep(mean(ddlact$ynt),length(timeL.seq)),
  focelo=rep(mean(ddlact$s_focelo),length(timeL.seq)),
  grpsize=rep(mean(ddlact$s_grpsize),length(timeL.seq)),
  infage=timeL.seq
)

predL <- link(condition_Lactating, n=1000 , data=d.pred_timeL, replace=
               list(ap_id=a_foc_z , al_id=a_foc_z, WAIC=TRUE))
predL.median <- apply( predL , 2 , median )
predL.HPDI <- apply( predL , 2 , HPDI )

par(mfrow = c(1, 1), cex=1.1, mar=c(0,0,0,0), oma=c(4,3,1.1,4))
plot( reverse_cond ~ s_infage, data=ddlact , col=alpha("#33CCFF",0.1),pch=19 ,xaxt='n',xlim = c(-1.3,4.2),xlab=NA,yaxt='n',ylab=NA,ylim = c(0,5),cex=1)

predL.lines = predL[sample(nrow(predL),100,replace=F),]
for (i in 1:100){
  lines( timeL.seq , predL.lines[i,] ,lwd=2,col=alpha("#33CCFF",0.05))
}
lines( timeL.seq , predL.median , lwd=1,col="black")
lines( timeL.seq , predL.HPDI[1,],lty=2,lwd=1)
lines( timeL.seq , predL.HPDI[2,],lty=2,lwd=1)
axis(1, at = seq(from=-5 , to=5, by = 1) ,tck=-0.02 , labels=T )
axis(2, at = seq(from=-4 , to=5, by = 1) ,tck=-0.02 , labels=T )

mtext(side=1,line=2,text="time across lactation",cex=1.5)
mtext(side=2,line=2,text="body condition",cex=1.5)


#same frame
par(mfrow = c(1, 2), cex=1.1, mar=c(0,0,0,0), oma=c(4,3,1.1,4))
plot( reverse_cond ~ s_infage, data=ddpreg , col=alpha("#33CCFF",0.1),pch=19 ,xaxt='n',xlim = c(-2,1.7),xlab=NA,yaxt='n',ylab=NA,ylim = c(0,5),cex=1)

predP.lines = predP[sample(nrow(predP),100,replace=F),]
for (i in 1:100){
  lines( timeP.seq , predP.lines[i,] ,lwd=2,col=alpha("#33CCFF",0.05))
}
lines( timeP.seq , predP.median , lwd=1,col="black")
lines( timeP.seq , predP.HPDI[1,],lty=2,lwd=1)
lines( timeP.seq , predP.HPDI[2,],lty=2,lwd=1)
axis(1, at = seq(from=-5 , to=5, by = 1) ,tck=-0.02 , labels=T )
axis(2, at = seq(from=-4 , to=5, by = 1) ,tck=-0.02 , labels=T )

mtext(side=1,line=2,text="time across pregnancy",cex=1.5)
mtext(side=2,line=2,text="body condition",cex=1.5)

plot( reverse_cond ~ s_infage, data=ddlact , col=alpha("#33CCFF",0.1),pch=19 ,xaxt='n',xlim = c(-1.3,4.2),xlab=NA,yaxt='n',ylab=NA,ylim = c(0,5),cex=1)

predL.lines = predL[sample(nrow(predL),100,replace=F),]
for (i in 1:100){
  lines( timeL.seq , predL.lines[i,] ,lwd=2,col=alpha("#33CCFF",0.05))
}
lines( timeL.seq , predL.median , lwd=1,col="black")
lines( timeL.seq , predL.HPDI[1,],lty=2,lwd=1)
lines( timeL.seq , predL.HPDI[2,],lty=2,lwd=1)
axis(1, at = seq(from=-5 , to=5, by = 1) ,tck=-0.02 , labels=T )
#axis(2, at = seq(from=-4 , to=5, by = 1) ,tck=-0.02 , labels=T )

mtext(side=1,line=2,text="time across lactation",cex=1.5)





d <- read.csv("MomConditXWean.csv", header = TRUE, na.strings = "")

library(rethinking)

### groups
d$phg <- ifelse(d$grp=="PHG" , 1 , 0 )
d$enk <- ifelse(d$grp=="ENK" , 1 , 0 )
d$ynt <- ifelse(d$grp=="YNT" , 1 , 0 )
d$chk <- ifelse(d$grp=="CHK" , 1 , 0 )

#make character or factors into integers 
d$mom <- as.factor(d$mom)
d$mom <- as.integer(as.factor(d$mom))

dd <- d[complete.cases(d$elo,
                       d$weanage,
                       d$ibi,
                       d$month4), ]
d<-dd

d$s_elo <- (d$momelo - mean(d$momelo))/sd(d$momelo)
d$s_cond <- (d$month4 - mean(d$month4))/sd(d$month4)
d$s_weanage <- (d$weanage - mean(d$weanage))/sd(d$weanage)
d$s_ibi <- (d$ibi - mean(d$ibi))/sd(d$ibi)

d$mom <- as.factor(d$mom)
d$mom <- as.integer(as.factor(d$mom))

condWean <- map2stan(
  alist(
    
    wean ~ dnorm(mu , sigma ),
    
    mu <- ap + gp_id[id] + bp_enk*enk + bp_ynt*ynt + bp_elo*elo + 
      bp_cond*cond + bp_elo_cond*elo*cond,
    
    c(gp_id)[id] ~ dnorm(0,sigma_focal) ,
    c(ap) ~ dnorm(0,2) ,
    c(bp_enk,bp_ynt,bp_elo,bp_cond,bp_elo_cond) ~ dnorm(0,2) ,
    sigma_focal  ~ dexp(1)  ,
    sigma ~ dexp(1)
    
  ),
  data=list(
    cond = d$s_cond,
    id=d$mom,
    elo = d$s_elo,
    phg = d$phg,
    enk = d$enk,
    ynt = d$ynt,
    wean = d$s_weanage
  ),
  
  cores=3 , chains=3 , warmup=3000, iter=7000, WAIC=TRUE, types=list(adapt.delta=0.99)
)


condIBI <- map2stan(
  alist(
    
    ibi ~ dnorm(mu , sigma ),
    
    mu <- ap + gp_id[id] + bp_enk*enk + bp_ynt*ynt + bp_elo*elo + 
      bp_cond*cond + bp_elo_cond*elo*cond,
    
    c(gp_id)[id] ~ dnorm(0,sigma_focal) ,
    c(ap) ~ dnorm(0,2) ,
    c(bp_enk,bp_ynt,bp_elo,bp_cond,bp_elo_cond) ~ dnorm(0,2) ,
    sigma_focal  ~ dexp(1)  ,
    sigma ~ dexp(1)
    
  ),
  data=list(
    cond = d$s_cond,
    id=d$mom,
    elo = d$s_elo,
    phg = d$phg,
    enk = d$enk,
    ynt = d$ynt,
    wean = d$s_weanage,
    ibi = d$s_ibi
  ),
  
  cores=3 , chains=3 , warmup=3000, iter=7000, WAIC=TRUE, types=list(adapt.delta=0.99)
)





a_foc_z <- matrix(0,1000,length(unique(d$mom)))
elo.seq=seq(min(d$s_cond),max(d$s_cond),length=1000)

d.pred_elo <- list(
  id=rep(1,length(elo.seq)),
  enk=rep(mean(d$enk),length(elo.seq)),
  ynt=rep(mean(d$ynt),length(elo.seq)),
  elo=rep(mean(d$s_elo),length(elo.seq)),
  cond=elo.seq
)

pred <- link(condWean, n=1000 , data=d.pred_elo, replace=
               list(ap_id=a_foc_z , al_id=a_foc_z, WAIC=TRUE))
pred.median <- apply( pred , 2 , median )
pred.HPDI <- apply( pred , 2 , HPDI )

par(mfrow = c(1, 2), cex=1.1, mar=c(0,0,0,0), oma=c(4,3,1.1,4))
plot( s_weanage ~ s_cond, data=d , col=alpha("#33CCFF",0.1),pch=19 ,xaxt='n',xlim = c(-2.2,1.7),xlab=NA,yaxt='n',ylab=NA,ylim = c(-1.3,3.7),cex=1)

pred.lines = pred[sample(nrow(pred),100,replace=F),]
for (i in 1:100){
  lines( elo.seq , pred.lines[i,] ,lwd=2,col=alpha("#33CCFF",0.05))
}
lines( elo.seq , pred.median , lwd=1,col="black")
lines( elo.seq , pred.HPDI[1,],lty=2,lwd=1)
lines( elo.seq , pred.HPDI[2,],lty=2,lwd=1)
axis(1, at = seq(from=-5 , to=5, by = 1) ,tck=-0.02 , labels=T )
axis(2, at = seq(from=-4 , to=5, by = 1) ,tck=-0.02 , labels=T )

mtext(side=1,line=2,text="body condition",cex=1.5)
mtext(side=2,line=2,text="weaning age",cex=1.5)


#ibi and wean plot
#wean cond
a_foc_z <- matrix(0,1000,length(unique(d$mom)))
elo.seq=seq(min(d$s_cond),max(d$s_cond),length=1000)

d.pred_elo1 <- list(
  id=rep(1,length(elo.seq)),
  enk=rep(mean(d$enk),length(elo.seq)),
  ynt=rep(mean(d$ynt),length(elo.seq)),
  elo=rep(-1,length(elo.seq)),
  cond=elo.seq
)

d.pred_elo3 <- list(
  id=rep(1,length(elo.seq)),
  enk=rep(mean(d$enk),length(elo.seq)),
  ynt=rep(mean(d$ynt),length(elo.seq)),
  elo=rep(1,length(elo.seq)),
  cond=elo.seq
)
pred1 <- link(condWean, n=1000 , data=d.pred_elo1, replace=
                list(ap_id=a_foc_z , al_id=a_foc_z, WAIC=TRUE))
pred1.median <- apply( pred1 , 2 , median )
pred1.HPDI <- apply( pred1 , 2 , HPDI )

pred3 <- link(condWean, n=1000 , data=d.pred_elo3, replace=
                list(ap_id=a_foc_z , al_id=a_foc_z, WAIC=TRUE))
pred3.median <- apply( pred3 , 2 , median )
pred3.HPDI <- apply( pred3 , 2 , HPDI )

par(mfrow = c(1, 2), cex=1.1, mar=c(0,0,0,0), oma=c(4,3,1.1,4))

plot( s_weanage[d$s_elo==-1] ~ s_cond[d$s_elo==-1], data=d , col=alpha("#33CCFF",0.1),pch=19 ,xaxt='n',xlim = c(-2.2,1.7),xlab=NA,yaxt='n',ylab=NA,ylim = c(-1.3,3.4),cex=1)
pred1.lines = pred1[sample(nrow(pred1),100,replace=F),]
for (i in 1:100){
  lines( elo.seq , pred1.lines[i,] ,lwd=2,col=alpha("#33CCFF",0.05))
}
lines( elo.seq , pred1.median , lwd=1,col="black")
lines( elo.seq , pred1.HPDI[1,],lty=2,lwd=1)
lines( elo.seq , pred1.HPDI[2,],lty=2,lwd=1)
axis(1, at = seq(from=-5 , to=5, by = 1) ,tck=-0.02 , labels=T )
axis(2, at = seq(from=-4 , to=5, by = 1) ,tck=-0.02 , labels=T )
mtext( "low rank" , side=3 , line=-1, outer=FALSE , cex=1)

plot( s_weanage[d$s_elo==1] ~ s_cond[d$s_elo==1], data=d , col=alpha("#33CCFF",0.1),pch=19 ,xaxt='n',xlim = c(-2.2,1.7),xlab=NA,yaxt='n',ylab=NA,ylim = c(-1.4,2.6),cex=1)
pred3.lines = pred3[sample(nrow(pred3),100,replace=F),]
for (i in 1:100){
  lines( elo.seq , pred3.lines[i,] ,lwd=2,col=alpha("#33CCFF",0.05))
}
lines( elo.seq , pred3.median , lwd=1,col="black")
lines( elo.seq , pred3.HPDI[1,],lty=2,lwd=1)
lines( elo.seq , pred3.HPDI[2,],lty=2,lwd=1)
axis(1, at = seq(from=-5 , to=5, by = 1) ,tck=-0.02 , labels=T )
mtext( "high rank" , side=3 , line=-1, outer=FALSE , cex=1)

mtext("condition at 4months standardized" , side=1 , line=2, outer=TRUE , cex=1.5)
mtext("weaning age standardized" , side=2 , line=2, outer=TRUE , cex=1.5)


#ibi cond
a_foc_z <- matrix(0,1000,length(unique(d$mom)))
elo.seq=seq(min(d$s_cond),max(d$s_cond),length=1000)

d.pred_elo <- list(
  id=rep(1,length(elo.seq)),
  enk=rep(mean(d$enk),length(elo.seq)),
  ynt=rep(mean(d$ynt),length(elo.seq)),
  elo=rep(mean(d$s_elo),length(elo.seq)),
  cond=elo.seq
)

pred <- link(condIBI, n=1000 , data=d.pred_elo, replace=
               list(ap_id=a_foc_z , al_id=a_foc_z, WAIC=TRUE))
pred.median <- apply( pred , 2 , median )
pred.HPDI <- apply( pred , 2 , HPDI )

#next plot
plot( s_ibi ~ s_cond, data=d , col=alpha("#33CCFF",0.1),pch=19 ,xaxt='n',xlim = c(-2.2,1.7),xlab=NA,yaxt='n',ylab=NA,ylim = c(-1.4,2.6),cex=1)

pred.lines = pred[sample(nrow(pred),100,replace=F),]
for (i in 1:100){
  lines( elo.seq , pred.lines[i,] ,lwd=2,col=alpha("#33CCFF",0.05))
}
lines( elo.seq , pred.median , lwd=1,col="black")
lines( elo.seq , pred.HPDI[1,],lty=2,lwd=1)
lines( elo.seq , pred.HPDI[2,],lty=2,lwd=1)
axis(1, at = seq(from=-5 , to=5, by = 1) ,tck=-0.02 , labels=T )
#axis(2, at = seq(from=-4 , to=5, by = 1) ,tck=-0.02 , labels=T )

mtext(side=1,line=2,text="body condition",cex=1.5)
mtext(side=4,line=1,text="ibi",cex=1.5)
