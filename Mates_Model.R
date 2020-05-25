### repro state 
d$focpreg <- ifelse(d$focstate=="preg" , 1 , 0 )
d$foclact <- ifelse(d$focstate=="lact" , 1 , 0 )
d$focflat <- ifelse(d$focstate=="flat" , 1 , 0 )
d$focswol <- ifelse(d$focstate=="swol" , 1 , 0 )

d$prtpreg <- ifelse(d$prtstate=="preg" , 1 , 0 )
d$prtlact <- ifelse(d$prtstate=="lact" , 1 , 0 )
d$prtflat <- ifelse(d$prtstate=="flat" , 1 , 0 )
d$prtswol <- ifelse(d$prtstate=="swol" , 1 , 0 )

### groups
d$phg <- ifelse(d$focgrp=="PHG" , 1 , 0 )
d$enk <- ifelse(d$focgrp=="ENK" , 1 , 0 )
d$ynt <- ifelse(d$focgrp=="YNT" , 1 , 0 )
d$chk <- ifelse(d$focgrp=="CHK" , 1 , 0 )

### kin categories
d$mom <- ifelse(d$kincat=="mom" , 1 , 0 )
d$sis <- ifelse(d$kincat=="sisters" , 1 , 0 )
d$gma <- ifelse(d$kincat=="grandma" , 1 , 0 )
d$aunt <- ifelse(d$kincat=="aunt" , 1 , 0 )
d$distant <- ifelse(d$kincat=="distant" , 1 , 0 )
d$nonkin <- ifelse(d$kincat=="nonkin" , 1 , 0 )

d$id1 <- as.factor(d$focal)
d$id2 <- as.factor(d$partner)
d$dyad <-as.factor(d$co_id)

d$id1 <- as.integer(as.factor(d$focal))
d$id2 <- as.integer(as.factor(d$partner))
d$dyad <- as.integer(as.factor(d$co_id))

#standardize
d$s_elodif <- d$rankdiff
d$s_hrs <- (d$coreshrs - mean(d$coreshrs))/sd(d$coreshrs)
d$s_osr <- (d$focosr - mean(d$focosr))/sd(d$focosr)

dd <- d[complete.cases(d$aggr1,d$id1,d$id2,d$dyad,d$focelo01,d$prtelo01,d$phg,d$enk,d$ynt,
                       d$focpreg,d$foclact,d$focswol,d$focflat,d$prtpreg,d$prtlact,d$prtswol,d$prtflat,
                       d$mom,d$sis,d$gma,d$aunt,d$distant,d$nonkin,
                       d$s_osr), ]

dd$s_focelo <- (dd$focelo01 - mean(dd$focelo01))/sd(dd$focelo01)
dd$s_prtelo <- (dd$prtelo01 - mean(dd$prtelo01))/sd(dd$prtelo01)

dd<-dd[!(dd$chk==1),]
dd<-dd[!(dd$focswol==0),]

dd$id1=as.factor(dd$id1)
dd$id2=as.factor(dd$id2)
dd$dyad=as.factor(dd$dyad)

unique(dd$id1[!(dd$id1 %in% dd$id2)]) 
unique(dd$id2[!(dd$id2 %in% dd$id1)]) 
dd$id1=droplevels(dd$id1)
dd$id2=droplevels(dd$id2)

dd$id1=as.integer(as.factor(dd$id1))
dd$id2=as.integer(as.factor(dd$id2))

Mates_Model <- map2stan(
  alist(
    aggr1 ~ dzipois( p , lambda ),
    
    logit(p) <- ap + s_hrs + gp_id[id1] + rp_id[id2] + ap_dyad[dyad] +
      bp_mom*mom + bp_sis*sis + bp_gma*gma + bp_aunt*aunt + bp_distant*distant +
      bp_enk*enk + bp_ynt*ynt + bp_elo1*focelo + bp_elo2*prtelo + bp_elointer*focelo*prtelo +
      bp_prtflat*prtflat + bp_prtpreg*prtpreg + bp_prtlact*prtlact + bp_osr*s_osr,
    
    log(lambda) <- al + s_hrs + gl_id[id1] + rl_id[id2] + al_dyad[dyad] +
      bl_mom*mom + bl_sis*sis + bl_gma*gma + bl_aunt*aunt + bl_distant*distant +
      bl_enk*enk + bl_ynt*ynt + bl_elo1*focelo + bl_elo2*prtelo + bl_elointer*focelo*prtelo +
      bl_prtflat*prtflat + bl_prtpreg*prtpreg + bl_prtlact*prtlact + bl_osr*s_osr,
    
    c(gp_id,gl_id)[id1] ~ dmvnormNC(sigma_id1, Rho_id1),
    c(rp_id,rl_id)[id2] ~ dmvnormNC(sigma_id2, Rho_id2),
    c(al_dyad,ap_dyad)[dyad] ~ dmvnormNC(sigma_dyad, Rho_dyad),
    c(ap,al) ~ dstudent(2,0,1),
    c(bp_mom,bp_sis,bp_gma,bp_aunt,bp_distant,
      bl_mom,bl_sis,bl_gma,bl_aunt,bl_distant,
      bp_enk,bp_ynt,bp_elo1,bp_elo2,bp_elointer,bp_prtflat,bp_prtpreg,bp_prtlact,
      bl_enk,bl_ynt,bl_elo1,bl_elo2,bl_elointer,bl_prtflat,bl_prtpreg,bl_prtlact,bp_osr,bl_osr) ~ dnorm(0,1),
    c(sigma_id1,sigma_id2,sigma_dyad) ~ dexp(1),
    c(Rho_id1,Rho_id2,Rho_dyad) ~ dlkjcorr(3)
    
  ),
  data=list(
    aggr1=dd$aggr1,
    s_hrs=dd$s_hrs,
    id1=dd$id1,
    id2=dd$id2,
    dyad = dd$dyad,
    mom = dd$mom,
    sis = dd$sis,
    gma = dd$gma,
    aunt = dd$aunt,
    distant = dd$distant,
    nonkin = dd$nonkin,
    focelo = dd$s_focelo,
    prtelo = dd$s_prtelo,
    phg = dd$phg,
    enk = dd$enk,
    ynt = dd$ynt,
    prtpreg = dd$prtpreg,
    prtlact = dd$prtlact,
    prtswol = dd$prtswol,
    prtflat = dd$prtflat,
    s_osr = dd$s_osr
  ),
  cores=3 , chains=2 , warmup=1000, iter=5000, WAIC=TRUE, types=list(adapt.delta=0.99)
)

Mates_Model_output=precis(Mates_Model)@output
write.csv(Mates_Model_output, file = "Mates_Model output.csv")


# limit analysis to swollen vs swollen 
ddsub <- subset(dd, prtswol==1)

ddsub$id1=as.factor(ddsub$id1)
ddsub$id2=as.factor(ddsub$id2)
ddsub$dyad=as.factor(ddsub$dyad)

unique(ddsub$id1[!(ddsub$id1 %in% ddsub$id2)]) 
unique(ddsub$id2[!(ddsub$id2 %in% ddsub$id1)]) 

ddsub$dyad=droplevels(ddsub$dyad)
ddsub$id1=droplevels(ddsub$id1)
ddsub$id2=droplevels(ddsub$id2)

ddsub$id1=as.integer(as.factor(ddsub$id1))
ddsub$id2=as.integer(as.factor(ddsub$id2))
ddsub$dyad=as.integer(as.factor(ddsub$dyad))

Mates_Model_swolvswol <- map2stan(
  alist(
    aggr1 ~ dzipois( p , lambda ),
    
    logit(p) <- ap + s_hrs + gp_id[id1] + rp_id[id2] + ap_dyad[dyad] +
      bp_mom*mom + bp_sis*sis + bp_gma*gma + bp_aunt*aunt + bp_distant*distant +
      bp_enk*enk + bp_ynt*ynt + bp_elo1*focelo + bp_elo2*prtelo + bp_elointer*focelo*prtelo +
      bp_osr*s_osr,
    
    log(lambda) <- al + s_hrs + gl_id[id1] + rl_id[id2] + al_dyad[dyad] +
      bl_mom*mom + bl_sis*sis + bl_gma*gma + bl_aunt*aunt + bl_distant*distant +
      bl_enk*enk + bl_ynt*ynt + bl_elo1*focelo + bl_elo2*prtelo + bl_elointer*focelo*prtelo +
      bl_osr*s_osr,
    
    c(gp_id,gl_id)[id1] ~ dmvnormNC(sigma_id1, Rho_id1),
    c(rp_id,rl_id)[id2] ~ dmvnormNC(sigma_id2, Rho_id2),
    c(al_dyad,ap_dyad)[dyad] ~ dmvnormNC(sigma_dyad, Rho_dyad),
    c(ap,al) ~ dstudent(2,0,1),
    c(bp_mom,bp_sis,bp_gma,bp_aunt,bp_distant,
      bl_mom,bl_sis,bl_gma,bl_aunt,bl_distant,
      bp_enk,bp_ynt,bp_elo1,bp_elo2,bp_elointer,
      bl_enk,bl_ynt,bl_elo1,bl_elo2,bl_elointer,bp_osr,bl_osr) ~ dnorm(0,1),
    c(sigma_id1,sigma_id2,sigma_dyad) ~ dexp(1),
    c(Rho_id1,Rho_id2,Rho_dyad) ~ dlkjcorr(3)
    
  ),
  data=list(
    aggr1=ddsub$aggr1,
    s_hrs=ddsub$s_hrs,
    id1=ddsub$id1,
    id2=ddsub$id2,
    dyad = ddsub$dyad,
    mom = ddsub$mom,
    sis = ddsub$sis,
    gma = ddsub$gma,
    aunt = ddsub$aunt,
    distant = ddsub$distant,
    nonkin = ddsub$nonkin,
    focelo = ddsub$s_focelo,
    prtelo = ddsub$s_prtelo,
    phg = ddsub$phg,
    enk = ddsub$enk,
    ynt = ddsub$ynt,
    prtpreg = ddsub$prtpreg,
    prtlact = ddsub$prtlact,
    prtswol = ddsub$prtswol,
    prtflat = ddsub$prtflat,
    s_osr = ddsub$s_osr
  ),
  cores=3 , chains=2 , warmup=2000, iter=5000, WAIC=TRUE, types=list(adapt.delta=0.99)
)
Mates_Model_swolvswol_output=precis(Mates_Model_swolvswol)@output
write.csv(Mates_Model_swolvswol_output, file = "Mates_Model_swolvswol output.csv")


# plot
a_foc_z <- matrix(0,1000,length(unique(dd$id1)))
a_prt_z <- matrix(0,1000,length(unique(dd$id2)))
a_dyad_z <- matrix(0,1000,length(unique(dd$dyad)))


d.pred.flat <- list(
  id1=c(1,1),
  id2=c(1,1),
  dyad=c(1,1),
  s_hrs=rep(mean(dd$s_hrs),2),
  mom=rep(mean(dd$mom),2),
  sis=rep(mean(dd$sis),2),
  gma=rep(mean(dd$gma),2),
  aunt=rep(mean(dd$aunt),2),
  distant=rep(mean(dd$distant),2),
  enk=rep(mean(dd$enk),2),
  ynt=rep(mean(dd$ynt),2),
  focelo=rep(mean(dd$s_focelo),2),
  prtelo=rep(mean(dd$s_prtelo),2),
  prtflat=c(1,1),
  prtpreg=c(0,0),
  prtlact=c(0,0),
  s_osr=rep(mean(dd$s_osr),2)
)

d.pred.lact <- list(
  id1=c(1,1),
  id2=c(1,1),
  dyad=c(1,1),
  s_hrs=rep(mean(dd$s_hrs),2),
  mom=rep(mean(dd$mom),2),
  sis=rep(mean(dd$sis),2),
  gma=rep(mean(dd$gma),2),
  aunt=rep(mean(dd$aunt),2),
  distant=rep(mean(dd$distant),2),
  enk=rep(mean(dd$enk),2),
  ynt=rep(mean(dd$ynt),2),
  focelo=rep(mean(dd$s_focelo),2),
  prtelo=rep(mean(dd$s_prtelo),2),
  prtflat=c(0,0),
  prtpreg=c(0,0),
  prtlact=c(1,1),
  s_osr=rep(mean(dd$s_osr),2)
)

d.pred.preg <- list(
  id1=c(1,1),
  id2=c(1,1),
  dyad=c(1,1),
  s_hrs=rep(mean(dd$s_hrs),2),
  mom=rep(mean(dd$mom),2),
  sis=rep(mean(dd$sis),2),
  gma=rep(mean(dd$gma),2),
  aunt=rep(mean(dd$aunt),2),
  distant=rep(mean(dd$distant),2),
  enk=rep(mean(dd$enk),2),
  ynt=rep(mean(dd$ynt),2),
  focelo=rep(mean(dd$s_focelo),2),
  prtelo=rep(mean(dd$s_prtelo),2),
  prtflat=c(0,0),
  prtpreg=c(1,1),
  prtlact=c(0,0),
  s_osr=rep(mean(dd$s_osr),2)
)

d.pred.swol <- list(
  id1=c(1,1),
  id2=c(1,1),
  dyad=c(1,1),
  s_hrs=rep(mean(dd$s_hrs),2),
  mom=rep(mean(dd$mom),2),
  sis=rep(mean(dd$sis),2),
  gma=rep(mean(dd$gma),2),
  aunt=rep(mean(dd$aunt),2),
  distant=rep(mean(dd$distant),2),
  enk=rep(mean(dd$enk),2),
  ynt=rep(mean(dd$ynt),2),
  focelo=rep(mean(dd$s_focelo),2),
  prtelo=rep(mean(dd$s_prtelo),2),
  prtflat=c(0,0),
  prtpreg=c(0,0),
  prtlact=c(0,0),
  s_osr=rep(mean(dd$s_osr),2)
)

link2flat <- link(Mates_Model , n=1000 , data=d.pred.flat, replace=
                    list(ap_id=a_foc_z , ap_dyad=a_dyad_z, al_id=a_foc_z , al_dyad=a_dyad_z), WAIC=TRUE)
pred.flat <- link2flat$lambda*(1-link2flat$p)
pred.flat.p <- (1-link2flat$p)
pred.flat.lambda <- link2flat$lambda
median(pred.flat[,1])
HPDI(pred.flat[,1])

link2lact <- link(Mates_Model , n=1000 , data=d.pred.lact, replace=
                    list(ap_id=a_foc_z , ap_dyad=a_dyad_z, al_id=a_foc_z , al_dyad=a_dyad_z), WAIC=TRUE)
pred.lact <- link2lact$lambda*(1-link2lact$p)
pred.lact.p <- (1-link2lact$p)
pred.lact.lambda <- link2lact$lambda
median(pred.lact[,1])
HPDI(pred.lact[,1])

link2preg <- link(Mates_Model , n=1000 , data=d.pred.preg, replace=
                    list(ap_id=a_foc_z , ap_dyad=a_dyad_z, al_id=a_foc_z , al_dyad=a_dyad_z), WAIC=TRUE)
pred.preg <- link2preg$lambda*(1-link2preg$p)
pred.preg.p <- (1-link2preg$p)
pred.preg.lambda <- link2preg$lambda
median(pred.preg[,1])
HPDI(pred.preg[,1])

link2swol <- link(Mates_Model , n=1000 , data=d.pred.swol, replace=
                    list(ap_id=a_foc_z , ap_dyad=a_dyad_z,al_id=a_foc_z , al_dyad=a_dyad_z), WAIC=TRUE)
pred.swol <- link2swol$lambda*(1-link2swol$p)
pred.swol.p <- (1-link2swol$p)
pred.swol.lambda <- link2swol$lambda
median(pred.swol[,1])
HPDI(pred.swol[,1])

par(mfrow=c(1,1), cex=1, mar=c(0,0,0,0), oma=c(2,2,2,2))

dens(pred.lact[,1], xlim=c(0,.1) , ylim=c(-.035,165) , xlab="" , ylab="" , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE , cex.lab=1.5 , cex.axis=0.5 , adj=0.1)
ll <- dd$aggr1[dd$prtlact==1]
points(ll, rep(-.025,length(ll)), pch="|" , col=col.alpha("#33CCFF", alpha=0.5) , cex=0.5 )
shade( density(pred.lact[,1]) , lim= as.vector(HPDI(pred.lact[,1], prob=0.99999)) , col = col.alpha("#33CCFF", 0.5))
shade( density(pred.preg[,1]) , lim= as.vector(HPDI(pred.preg[,1], prob=0.99999)) , col = col.alpha("orange1", 0.5))
shade( density(pred.swol[,1]) , lim= as.vector(HPDI(pred.swol[,1], prob=0.99999)) , col = col.alpha("green", 0.5))
shade( density(pred.flat[,1]) , lim= as.vector(HPDI(pred.flat[,1], prob=0.99999)) , col = col.alpha("purple", 0.5))
ll <- dd$aggr1[dd$prtpreg==1]
points(ll, rep(-.025,length(ll)), pch="|" , col=col.alpha("orange1", alpha=0.5) , cex=0.5 )
ll <- dd$aggr1[dd$prtswol==1]
points(ll, rep(-.025,length(ll)), pch="|" , col=col.alpha("green", alpha=0.5) , cex=0.5 )
ll <- dd$aggr1[dd$prtlat==1]
points(ll, rep(-.025,length(ll)), pch="|" , col=col.alpha("purple", alpha=0.5) , cex=0.5 )

abline(v=median(pred.lact[,1]) , lty=1)
abline(v=median(pred.preg[,1]) , lty=2)
abline(v=median(pred.swol[,1]) , lty=3)
abline(v=median(pred.flat[,1]) , lty=4)

legend(.07,160, legend = c("", "", "", ""),
       col=c(1,1)  , lty=c(1,2,3,4),
       lw=1 , cex=.85, bty="n")

legend(.08,160,, legend = c("lactating", "pregnant", "swollen", "flat"),
       col=c(col.alpha("#33CCFF", 0.5) , col.alpha("orange1", 0.5) ,  col.alpha("green", 0.5) , col.alpha("purple", 0.5)) , pch=c(15,17),
       pt.cex=2 , cex=.85, bty="n")

mtext("Rate of aggression from swollen females" , side=1 , line=.3, outer=TRUE , cex=1.5)
#axis(1, at = seq(from=0 , to=1, by = .1) ,tck=-0.1 , labels=T )
axis(1, at = seq(from=0 , to=1, by = .05) ,tck=-0.01 , labels=T )


#osr
a_foc_z <- matrix(0,1000,length(unique(dd$Focal)))
osr.seq=seq(min(dd$s_osr),max(dd$s_osr),length=1000)


d.pred_osr <- list(
  Focal=rep(1,length(osr.seq)),
  s_hrs=rep(mean(dd$s_hrs),length(osr.seq)),
  id1=rep(1,length(osr.seq)),
  id2=rep(1,length(osr.seq)),
  dyad=rep(1,length(osr.seq)),
  mom=rep(mean(dd$mom),length(osr.seq)),
  sis=rep(mean(dd$sis),length(osr.seq)),
  gma=rep(mean(dd$gma),length(osr.seq)),
  aunt=rep(mean(dd$aunt),length(osr.seq)),
  distant=rep(mean(dd$distant),length(osr.seq)),
  focelo=rep(mean(dd$s_focelo),length(osr.seq)),
  prtelo=rep(mean(dd$s_prtelo),length(osr.seq)),
  prtflat=rep(mean(dd$prtflat),length(osr.seq)),
  prtpreg=rep(mean(dd$prtpreg),length(osr.seq)),
  prtlact=rep(mean(dd$prtlact),length(osr.seq)),
  enk=rep(mean(dd$prtlact),length(osr.seq)),
  ynt=rep(mean(dd$prtlact),length(osr.seq)),
  s_osr=osr.seq
)


pred <- link(Mates_Model , n=1000 , data=d.pred_osr, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
pred.osr <- (1-pred$p)*pred$lambda
pred.osr.median <- apply( pred.osr , 2 , median )
pred.osr.HPDI <- apply( pred.osr , 2 , HPDI )

par(mar=c(0,0,0,0),oma=c(3,3,3,3))
plot( aggr1 ~ s_osr, data=dd , col=alpha("#33CCFF",0.2),pch=19 ,xaxt='n',xlim = c(-.7,3.4),xlab=NA,yaxt='n',ylab=NA,ylim = c(0,.1),cex=1)

pred.lines = pred.osr[sample(nrow(pred.osr),100,replace=F),]
for (i in 1:100){
  lines( osr.seq , pred.lines[i,] ,lwd=2,col=alpha("#33CCFF",0.05))
}
lines( osr.seq , pred.osr.median , lwd=1,col="black")
lines( osr.seq , pred.osr.HPDI[1,],lty=2,lwd=1)
lines( osr.seq , pred.osr.HPDI[2,],lty=2,lwd=1)

mtext(side=1,line=1.5,text="osr",cex=1.5)
mtext(side=2,line=1.5,text="Rate of aggression from swollen females",cex=1.3)
axis(1, at = seq(from=-1 , to=3, by = 1) ,tck=-0.02 , labels=T )
axis(2, at = seq(from=-1 , to=3, by = .05) ,tck=-0.02 , labels=T )



#osr limited to swollen vs swollen 
a_foc_z <- matrix(0,1000,length(unique(ddsub$Focal)))
osr.seq=seq(min(ddsub$s_osr),max(ddsub$s_osr),length=1000)

d.pred_osr <- list(
  Focal=rep(1,length(osr.seq)),
  s_hrs=rep(mean(ddsub$s_hrs),length(osr.seq)),
  id1=rep(1,length(osr.seq)),
  id2=rep(1,length(osr.seq)),
  dyad=rep(1,length(osr.seq)),
  mom=rep(mean(ddsub$mom),length(osr.seq)),
  sis=rep(mean(ddsub$sis),length(osr.seq)),
  gma=rep(mean(ddsub$gma),length(osr.seq)),
  aunt=rep(mean(ddsub$aunt),length(osr.seq)),
  distant=rep(mean(ddsub$distant),length(osr.seq)),
  focelo=rep(mean(ddsub$s_focelo),length(osr.seq)),
  prtelo=rep(mean(ddsub$s_prtelo),length(osr.seq)),
  prtflat=rep(mean(ddsub$prtflat),length(osr.seq)),
  prtpreg=rep(mean(ddsub$prtpreg),length(osr.seq)),
  prtlact=rep(mean(ddsub$prtlact),length(osr.seq)),
  enk=rep(mean(ddsub$prtlact),length(osr.seq)),
  ynt=rep(mean(ddsub$prtlact),length(osr.seq)),
  s_osr=osr.seq
)


pred <- link(Mates_Model_swolvswol , n=1000 , data=d.pred_osr, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
pred.osr <- (1-pred$p)*pred$lambda
pred.osr.median <- apply( pred.osr , 2 , median )
pred.osr.HPDI <- apply( pred.osr , 2 , HPDI )


#pdf(file="mates model osr.pdf" , width=7, height=7)

par(mar=c(0,0,0,0),oma=c(3,3,3,3))
plot( aggr1 ~ s_osr, data=ddsub , col=alpha("#33CCFF",0.2),pch=19 ,xaxt='n',xlim = c(-.7,3.4),xlab=NA,yaxt='n',ylab=NA,ylim = c(0,.1),cex=1)

pred.lines = pred.osr[sample(nrow(pred.osr),100,replace=F),]
for (i in 1:100){
  lines( osr.seq , pred.lines[i,] ,lwd=2,col=alpha("#33CCFF",0.05))
}
lines( osr.seq , pred.osr.median , lwd=1,col="black")
lines( osr.seq , pred.osr.HPDI[1,],lty=2,lwd=1)
lines( osr.seq , pred.osr.HPDI[2,],lty=2,lwd=1)

mtext(side=1,line=1.5,text="osr",cex=1.5)
mtext(side=2,line=1.5,text="Rate of aggression from swollen females",cex=1.3)
axis(1, at = seq(from=-1 , to=3, by = 1) ,tck=-0.02 , labels=T )
axis(2, at = seq(from=-1 , to=3, by = .05) ,tck=-0.02 , labels=T )

