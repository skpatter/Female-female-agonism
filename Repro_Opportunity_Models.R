


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

#make character or factors into integers 
d$id1 <- as.factor(d$focal)
d$id2 <- as.factor(d$partner)
d$dyad <-as.factor(d$co_id)

d$id1 <- as.integer(as.factor(d$focal))
d$id2 <- as.integer(as.factor(d$partner))
d$dyad <- as.integer(as.factor(d$co_id))

#pre-process data 
d$s_hrs <- (d$coreshrs - mean(d$coreshrs))/sd(d$coreshrs)
d$s_osr <- (d$focosr - mean(d$focosr))/sd(d$focosr)

dd <- d[complete.cases(d$aggr1rate,d$aggr2rate,d$id1,d$id2,d$dyad,d$focelo01,d$prtelo01,d$phg,d$enk,d$ynt,
                       d$focpreg,d$foclact,d$focswol,d$focflat,d$prtpreg,d$prtlact,d$prtswol,d$prtflat,
                       d$mom,d$sis,d$gma,d$aunt,d$distant,d$nonkin,
                       d$s_osr), ]
dd <- d[complete.cases(d$aggr1rate,d$aggr2rate,d$id1,d$id2,d$dyad,d$focelo01,d$prtelo01,d$phg,d$enk,d$ynt,
                       d$focpreg,d$foclact,d$focswol,d$focflat,d$prtpreg,d$prtlact,d$prtswol,d$prtflat,
                       d$mom,d$sis,d$gma,d$aunt,d$distant,d$nonkin,
                       d$s_osr,d$Panextsire,d$Sirenextsire,d$PaNextPa), ]
dd$s_focelo <- (dd$focelo01 - mean(dd$focelo01))/sd(dd$focelo01)
dd$s_prtelo <- (dd$prtelo01 - mean(dd$prtelo01))/sd(dd$prtelo01)

dd<-dd[!(dd$chk==1),]

#redo the as.integer bc NAs were removed
dd$id1=as.factor(dd$id1)
dd$id2=as.factor(dd$id2)
dd$dyad=as.factor(dd$dyad)

# id1 in id2 and id2 in id1 columns...
unique(dd$id1[!(dd$id1 %in% dd$id2)]) 
unique(dd$id2[!(dd$id2 %in% dd$id1)]) 

dd$dyad=droplevels(dd$dyad)
dd$id1=droplevels(dd$id1)
dd$id2=droplevels(dd$id2)

#factors back to integers
dd$id1=as.integer(as.factor(dd$id1))
dd$id2=as.integer(as.factor(dd$id2))
dd$dyad=as.integer(as.factor(dd$dyad))



Repro_Model <- map2stan(
  alist(
    aggr1 ~ dzipois( p , lambda ),
    
    logit(p) <- ap + s_hrs + gp_id[id1] + rp_id[id2] + ap_dyad[dyad] +
      bp_mom*mom + bp_sis*sis + bp_gma*gma + bp_aunt*aunt + bp_distant*distant +
      bp_enk*enk + bp_ynt*ynt + bp_elo1*focelo + bp_elo2*prtelo + bp_elointer*focelo*prtelo +
      bp_prtflat*prtflat + bp_prtpreg*prtpreg + bp_prtlact*prtlact ,
    
    log(lambda) <- al + s_hrs + gl_id[id1] + rl_id[id2] + al_dyad[dyad] +
      bl_mom*mom + bl_sis*sis + bl_gma*gma + bl_aunt*aunt + bl_distant*distant +
      bl_enk*enk + bl_ynt*ynt + bl_elo1*focelo + bl_elo2*prtelo + bl_elointer*focelo*prtelo +
      bl_prtflat*prtflat + bl_prtpreg*prtpreg + bl_prtlact*prtlact,
    
    c(gp_id,gl_id)[id1] ~ dmvnormNC(sigma_id1, Rho_id1),
    c(rp_id,rl_id)[id2] ~ dmvnormNC(sigma_id2, Rho_id2),
    c(al_dyad,ap_dyad)[dyad] ~ dmvnormNC(sigma_dyad, Rho_dyad),
    c(ap,al) ~ dstudent(2,0,1),
    c(bp_mom,bp_sis,bp_gma,bp_aunt,bp_distant,
      bl_mom,bl_sis,bl_gma,bl_aunt,bl_distant,
      bp_enk,bp_ynt,bp_elo1,bp_elo2,bp_elointer,bp_prtflat,bp_prtpreg,bp_prtlact,
      bl_enk,bl_ynt,bl_elo1,bl_elo2,bl_elointer,bl_prtflat,bl_prtpreg,bl_prtlact) ~ dnorm(0,1),
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
  cores=3 , chains=2 , warmup=1000, iter=2000, WAIC=TRUE, types=list(adapt.delta=0.99)
)
Repro_Model_output=precis(Repro_Model)@output
write.csv(Repro_Model_output, file = "Repro_Model output.csv")

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
  prtlact=c(0,0)
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
  prtlact=c(1,1)
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
  prtlact=c(0,0)
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
  prtlact=c(0,0)
)

link2flat <- link(Repro_Model , n=1000 , data=d.pred.flat, replace=
                    list(ap_id=a_foc_z , ap_dyad=a_dyad_z, am_id=a_foc_z , am_dyad=a_dyad_z), WAIC=TRUE)
pred.flat <- (1-link2flat$p)*link2flat$lambda
pred.flat.p <- (1-link2flat$p)
pred.flat.lambda <- link2flat$lambda
median(pred.flat[,1])
HPDI(pred.flat[,1])

link2lact <- link(Repro_Model , n=1000 , data=d.pred.lact, replace=
                    list(ap_id=a_foc_z , ap_dyad=a_dyad_z, am_id=a_foc_z , am_dyad=a_dyad_z), WAIC=TRUE)
pred.lact <- (1-link2lact$p)*link2lact$lambda
pred.lact.p <- (1-link2lact$p)
pred.lact.lambda <- link2lact$lambda
median(pred.lact[,1])
HPDI(pred.lact[,1])

link2preg <- link(Repro_Model , n=1000 , data=d.pred.preg, replace=
                    list(ap_id=a_foc_z , ap_dyad=a_dyad_z, am_id=a_foc_z , am_dyad=a_dyad_z), WAIC=TRUE)
pred.preg <- (1-link2preg$p)*link2preg$lambda
pred.preg.p <- (1-link2preg$p)
pred.preg.lambda <- link2preg$lambda
median(pred.preg[,1])
HPDI(pred.preg[,1])

link2swol <- link(Repro_Model , n=1000 , data=d.pred.swol, replace=
                    list(ap_id=a_foc_z , ap_dyad=a_dyad_z,am_id=a_foc_z , am_dyad=a_dyad_z), WAIC=TRUE)
pred.swol <- (1-link2swol$p)*link2swol$lambda
pred.swol.p <- (1-link2swol$p)
pred.swol.lambda <- link2swol$lambda
median(pred.swol[,1])
HPDI(pred.swol[,1])

#pdf(file="repro model joint_agg.pdf" , width=7, height=7)

par(mar=c(2,2,2,2))
par(mfrow=c(1,1), cex=1, mar=c(0,0,0,0), oma=c(2,2,2,2))
#dataplots.diff <- cbind(pred.lact[,1],pred.preg[,1],pred.swol[,1],pred.flat[,1])
#sublist <- cbind(dd$prtlact , dd$prtpreg , dd$prtswol , dd$prtflat)
#plot.titles <- c("a) Lactating","b) Pregnant","c) Swollen","d) Flat")

dens(pred.lact[,1], xlim=c(0,.5) , ylim=c(-.035,25) , xlab="" , ylab="" , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE , cex.lab=1.5 , cex.axis=0.5 , adj=0.1)
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

legend(.15,20, legend = c("", "", "", ""),
       col=c(1,1)  , lty=c(1,2,3),
       lw=1 , cex=.85, bty="n")

legend(.17,20,, legend = c("lactating", "pregnant", "swollen", "flat"),
       col=c(col.alpha("#33CCFF", 0.5) , col.alpha("orange1", 0.5) ,  col.alpha("green", 0.5) , col.alpha("purple", 0.5)) , pch=c(15,17),
       pt.cex=2 , cex=.85, bty="n")

mtext("Rate of aggression from pregnant and lactating females" , side=1 , line=.3, outer=TRUE , cex=1.5)
#axis(1, at = seq(from=0 , to=1, by = .1) ,tck=-0.1 , labels=T )
axis(1, at = seq(from=0 , to=1, by = .1) ,tck=-0.05 , labels=T )



## preg/lact vs swollen 

dd <- d[complete.cases(d$aggr1rate,d$aggr2rate,d$id1,d$id2,d$dyad,d$focelo01,d$prtelo01,d$phg,d$enk,d$ynt,
                       d$focpreg,d$foclact,d$focswol,d$focflat,d$prtpreg,d$prtlact,d$prtswol,d$prtflat,
                       d$mom,d$sis,d$gma,d$aunt,d$distant,d$nonkin,
                       d$s_osr,d$Panextsire,d$Sirenextsire,d$PaNextPa), ]

dd$s_focelo <- (dd$focelo01 - mean(dd$focelo01))/sd(dd$focelo01)
dd$s_prtelo <- (dd$prtelo01 - mean(dd$prtelo01))/sd(dd$prtelo01)

dd<-dd[!(dd$chk==1),]

dd$id1=as.factor(dd$id1)
dd$id2=as.factor(dd$id2)
dd$dyad=as.factor(dd$dyad)

unique(dd$id1[!(dd$id1 %in% dd$id2)]) 
unique(dd$id2[!(dd$id2 %in% dd$id1)]) 

dd$dyad=droplevels(dd$dyad)
dd$id1=droplevels(dd$id1)
dd$id2=droplevels(dd$id2)

dd$id1=as.integer(as.factor(dd$id1))
dd$id2=as.integer(as.factor(dd$id2))
dd$dyad=as.integer(as.factor(dd$dyad))



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

Repro_Model_PA_Sire <- map2stan(
  alist(
    aggr1 ~ dzipois( p , lambda ),
    
    logit(p) <- ap + s_hrs + gp_id[id1] + rp_id[id2] + ap_dyad[dyad] +
      bp_mom*mom + bp_sis*sis + bp_gma*gma + bp_aunt*aunt + bp_distant*distant +
      bp_enk*enk + bp_ynt*ynt + bp_elo1*focelo + bp_elo2*prtelo + bp_elointer*focelo*prtelo +
      bp_nextPA*nextPA + bp_nextSire*nextSire,
    
    log(lambda) <- al + s_hrs + gl_id[id1] + rl_id[id2] + al_dyad[dyad] +
      bl_mom*mom + bl_sis*sis + bl_gma*gma + bl_aunt*aunt + bl_distant*distant +
      bl_enk*enk + bl_ynt*ynt + bl_elo1*focelo + bl_elo2*prtelo + bl_elointer*focelo*prtelo +
      bl_nextPA*nextPA + bp_nextSire*nextSire,
    
    c(gp_id,gl_id)[id1] ~ dmvnormNC(sigma_id1, Rho_id1),
    c(rp_id,rl_id)[id2] ~ dmvnormNC(sigma_id2, Rho_id2),
    c(al_dyad,ap_dyad)[dyad] ~ dmvnormNC(sigma_dyad, Rho_dyad),
    c(ap,al) ~ dstudent(2,0,1),
    c(bp_mom,bp_sis,bp_gma,bp_aunt,bp_distant,
      bl_mom,bl_sis,bl_gma,bl_aunt,bl_distant,
      bp_enk,bp_ynt,bp_elo1,bp_elo2,bp_elointer,bp_nextPA,bp_nextSire,
      bl_enk,bl_ynt,bl_elo1,bl_elo2,bl_elointer,bl_nextPA,bl_nextSire) ~ dnorm(0,1),
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
    s_osr = ddsub$s_osr,
    nextPA = ddsub$PaNextPa,
    nextSire = ddsub$Panextsire
  ),
  cores=3 , chains=2 , warmup=2000, iter=7000, WAIC=TRUE, types=list(adapt.delta=0.99)
)

Repro_Model_PAnextPA.Sireoutput=precis(Repro_Model_PA_Sire)@output
write.csv(Repro_Model_PAnextPA.Sireoutput, file = "Repro_Model PA_Sire output.csv")


#Next PA
a_foc_z <- matrix(0,1000,length(unique(ddsub$id1)))
a_prt_z <- matrix(0,1000,length(unique(ddsub$id2)))
a_dyad_z <- matrix(0,1000,length(unique(ddsub$dyad)))

d.pred.nextPA <- list(
  id1=c(1,1),
  id2=c(1,1),
  dyad=c(1,1),
  s_hrs=rep(mean(ddsub$s_hrs),2),
  mom=rep(mean(ddsub$mom),2),
  sis=rep(mean(ddsub$sis),2),
  gma=rep(mean(ddsub$gma),2),
  aunt=rep(mean(ddsub$aunt),2),
  distant=rep(mean(ddsub$distant),2),
  enk=rep(mean(ddsub$enk),2),
  ynt=rep(mean(ddsub$ynt),2),
  focelo=rep(mean(ddsub$s_focelo),2),
  prtelo=rep(mean(ddsub$s_prtelo),2),
  nextSire=rep(mean(ddsub$Panextsire),2),
  nextPA=c(0,1)
)

link2nextPA <- link(Repro_Model_PA_Sire , n=1000 , data=d.pred.nextPA, replace=
                      list(ap_id=a_foc_z , ap_dyad=a_dyad_z, am_id=a_foc_z , am_dyad=a_dyad_z), WAIC=TRUE)
pred.nextPA <- (1-link2nextPA$p)*link2nextPA$lambda
median(pred.nextPA[,1])
HPDI(pred.nextPA[,1])
median(pred.nextPA[,2])
HPDI(pred.nextPA[,2])

par(mfrow = c(1, 2), cex=1.1, mar=c(3,0,0,0), oma=c(2,2,2,2))

dens(pred.nextPA[,1], xlim=c(0,.5), xlab="", ylim=c(0,10), ylab="",yaxt='n',col="white" , xaxt='n', zero.line=FALSE, cex.lab=1.5, cex.axis=0.5)
ll <- ddsub$aggr1[ddsub$PaNextPa==0]
points(ll, rep(-.01,length(ll)), pch=15 , col=col.alpha("orange1", alpha=0.1) , cex=0.75)

shade( density(pred.nextPA[,1]) , lim= as.vector(HPDI(pred.nextPA[,1], prob=0.9999)) , col = col.alpha("orange1", 0.5))
shade( density(pred.nextPA[,2]) , lim= as.vector(HPDI(pred.nextPA[,2], prob=0.9999)) , col = col.alpha("#33CCFF", 0.5))
ll <- ddsub$aggr1[ddsub$PaNextPa==1]

points(ll, rep(-.03, length(ll)), pch=17 , col=col.alpha("#33CCFF", alpha=0.1) , cex=0.75)
axis(1, at = seq(from=0 , to=1, by = .1) ,tck=-0.02 , labels=T )
#axis(1, at = seq(from=0 , to=3, by = 1) ,tck=-0.01 , labels=F)
abline(v=median(pred.nextPA[,1]) , lty=1)
abline(v=median(pred.nextPA[,2]) , lty=2)

legend(.15,8, legend = c("", ""),
       col=c(1,1)  , lty=c(1,2),
       lw=1 , cex=.85, bty="n")

legend(.17,8,, legend = c("Different PA", "Same PA"),
       col=c(col.alpha("orange1", 0.5) , col.alpha("#33CCFF", 0.5) ) , pch=c(15,17),
       pt.cex=2 , cex=.85, bty="n")

#next Sire
a_foc_z <- matrix(0,1000,length(unique(ddsub$id1)))
a_prt_z <- matrix(0,1000,length(unique(ddsub$id2)))
a_dyad_z <- matrix(0,1000,length(unique(ddsub$dyad)))

d.pred.nextSire <- list(
  id1=c(1,1),
  id2=c(1,1),
  dyad=c(1,1),
  s_hrs=rep(mean(ddsub$s_hrs),2),
  mom=rep(mean(ddsub$mom),2),
  sis=rep(mean(ddsub$sis),2),
  gma=rep(mean(ddsub$gma),2),
  aunt=rep(mean(ddsub$aunt),2),
  distant=rep(mean(ddsub$distant),2),
  enk=rep(mean(ddsub$enk),2),
  ynt=rep(mean(ddsub$ynt),2),
  focelo=rep(mean(ddsub$s_focelo),2),
  prtelo=rep(mean(ddsub$s_prtelo),2),
  nextPA=rep(mean(ddsub$PaNextPa),2),
  nextSire=c(0,1)
  
)

link2nextSire <- link(Repro_Model_PA_Sire , n=1000 , data=d.pred.nextSire, replace=
                        list(ap_id=a_foc_z , ap_dyad=a_dyad_z, am_id=a_foc_z , am_dyad=a_dyad_z), WAIC=TRUE)
pred.nextSire <- (1-link2nextSire$p)*link2nextSire$lambda
median(pred.nextSire[,1])
HPDI(pred.nextSire[,1])
median(pred.nextSire[,2])
HPDI(pred.nextSire[,2])


dens(pred.nextSire[,1], xlim=c(0,.5), xlab="", ylim=c(0,10), ylab="", yaxt='n',col="white" , xaxt='n', zero.line=FALSE, cex.lab=1.5, cex.axis=0.5)
ll <- ddsub$aggr1[ddsub$Panextsire==0]
points(ll, rep(-.01,length(ll)), pch=15 , col=col.alpha("orange1", alpha=0.1) , cex=0.75)

shade( density(pred.nextSire[,1]) , lim= as.vector(HPDI(pred.nextSire[,1], prob=0.9999)) , col = col.alpha("orange1", 0.5))
shade( density(pred.nextSire[,2]) , lim= as.vector(HPDI(pred.nextSire[,2], prob=0.9999)) , col = col.alpha("#33CCFF", 0.5))
ll <- ddsub$aggr1[ddsub$Panextsire==1]

points(ll, rep(-.03, length(ll)), pch=17 , col=col.alpha("#33CCFF", alpha=0.1) , cex=0.75)
#axis(1, at = seq(from=0 , to=3, by = .5) ,tck=-0.01 , labels=F)
abline(v=median(pred.nextSire[,1]) , lty=1)
abline(v=median(pred.nextSire[,2]) , lty=2)

legend(.15,8, legend = c("", ""),
       col=c(1,1)  , lty=c(1,2),
       lw=1 , cex=.85, bty="n")

legend(.18,8,, legend = c("PA is not next sire", "PA is next sire"),
       col=c(col.alpha("orange1", 0.5) , col.alpha("#33CCFF", 0.5) ) , pch=c(15,17),
       pt.cex=2 , cex=.85, bty="n")
axis(1, at = seq(from=0 , to=3, by = .1) ,tck=-0.02 , labels=T )

mtext("Rate of aggression from pregnant/lactating to swollen females", side = 1, line = .1, outer = TRUE , cex=1.5)


