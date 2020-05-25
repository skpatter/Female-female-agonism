### repro state 
#focal
d$focpreg <- ifelse(d$focstate=="preg" , 1 , 0 )
d$foclact <- ifelse(d$focstate=="lact" , 1 , 0 )
d$focflat <- ifelse(d$focstate=="flat" , 1 , 0 )
d$focswol <- ifelse(d$focstate=="swol" , 1 , 0 )

#partner
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

#standardize
d$s_hrs2 <- (d$coreshrs - mean(d$coreshrs))/sd(d$coreshrs)

dd <- d[complete.cases(d$aggr1,
                       d$s_hrs,
                       d$id1,
                       d$id2,
                       d$focelo01,
                       d$prtelo01,
                       d$mom,
                       d$sis,
                       d$gma,
                       d$aunt,
                       d$distant,
                       d$nonkin,
                       d$phg,
                       d$enk,
                       d$ynt,
                       d$focpreg,
                       d$foclact,
                       d$focswol,
                       d$focflat,
                       d$foc_condit_rev),]

dd$s_focelo <- (dd$focelo01 - mean(dd$focelo01))/sd(dd$focelo01)
dd$s_prtelo <- (dd$prtelo01 - mean(dd$prtelo01))/sd(dd$prtelo01)
dd$s_foccondRev <- (dd$foc_condit_rev - mean(dd$foc_condit_rev))/sd(dd$foc_condit_rev)
dd$s_grpsize <- (dd$focgrpsize - mean(dd$focgrpsize))/sd(dd$focgrpsize)
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

Food_Model_Full <- map2stan(
  alist(
    aggr1 ~ dzipois( p , lambda ),
    
    logit(p) <- ap + s_hrs + gp_id[id1] + rp_id[id2] + ap_dyad[dyad] +
      bp_mom*mom + bp_sis*sis + bp_gma*gma + bp_aunt*aunt + bp_distant*distant +
      bp_enk*enk + bp_ynt*ynt + bp_elo1*focelo + bp_elo2*prtelo + bp_elointer*focelo*prtelo +
      bp_focflat*focflat + bp_focpreg*focpreg + bp_focswol*focswol + bp_foccond*cond + bp_grpsize*grpsize,
    
    log(lambda) <- al + s_hrs + gl_id[id1] + rl_id[id2] + al_dyad[dyad] +
      bl_mom*mom + bl_sis*sis + bl_gma*gma + bl_aunt*aunt + bl_distant*distant +
      bl_enk*enk + bl_ynt*ynt + bl_elo1*focelo + bl_elo2*prtelo + bl_elointer*focelo*prtelo +
      bl_focflat*focflat + bl_focpreg*focpreg + bl_focswol*focswol + bl_foccond*cond + bl_grpsize*grpsize,
    
    c(gp_id,gl_id)[id1] ~ dmvnormNC(sigma_id1, Rho_id1),
    c(rp_id,rl_id)[id2] ~ dmvnormNC(sigma_id2, Rho_id2),
    c(al_dyad,ap_dyad)[dyad] ~ dmvnormNC(sigma_dyad, Rho_dyad),
    c(ap,al) ~ dstudent(2,0,1),
    c(bp_mom,bp_sis,bp_gma,bp_aunt,bp_distant,
      bl_mom,bl_sis,bl_gma,bl_aunt,bl_distant,
      bp_enk,bp_ynt,bp_elo1,bp_elo2,bp_elointer,bp_focflat,bp_focpreg,bp_focswol,bp_foccond,bp_grpsize,
      bl_enk,bl_ynt,bl_elo1,bl_elo2,bl_elointer,bl_focflat,bl_focpreg,bl_focswol,bl_foccond,bl_grpsize) ~ dnorm(0,1),
    c(sigma_id1,sigma_id2,sigma_dyad) ~ dexp(1),
    c(Rho_id1,Rho_id2,Rho_dyad) ~ dlkjcorr(3)
    
  ),
  data=list(
    aggr1=dd$aggr1,
    s_hrs=dd$s_hrs2,
    id1=dd$id1,
    id2=dd$id2,
    dyad = dd$dyad,
    mom = dd$mom,
    sis = dd$sis,
    gma = dd$gma,
    aunt = dd$aunt,
    distant = dd$distant,
    nonkin = dd$nonkin,
    elodif = dd$s_elodif,
    focelo = dd$s_focelo,
    prtelo = dd$s_prtelo,
    phg = dd$phg,
    enk = dd$enk,
    ynt = dd$ynt,
    focpreg = dd$focpreg,
    foclact = dd$foclact,
    focswol = dd$focswol,
    focflat = dd$focflat,
    cond = dd$s_foccondRev,
    grpsize = dd$s_grpsize
  ),
  cores=3 , chains=2 , warmup=3500, iter=7000, WAIC=TRUE, types=list(adapt.delta=0.99)
)



### plot
a_foc_z <- matrix(0,1000,length(unique(dd$id1)))
a_prt_z <- matrix(0,1000,length(unique(dd$id2)))
a_dyad_z <- matrix(0,1000,length(unique(dd$dyad)))


d.pred.flat <- list(
  id1=c(1,1),
  id2=c(1,1),
  dyad=c(1,1),
  s_hrs=rep(mean(dd$s_hrs2),2),
  mom=rep(mean(dd$mom),2),
  sis=rep(mean(dd$sis),2),
  gma=rep(mean(dd$gma),2),
  aunt=rep(mean(dd$aunt),2),
  distant=rep(mean(dd$distant),2),
  enk=rep(mean(dd$enk),2),
  ynt=rep(mean(dd$ynt),2),
  focelo=rep(mean(dd$s_focelo),2),
  prtelo=rep(mean(dd$s_prtelo),2),
  cond=rep(mean(dd$s_foccondRev),2),
  grpsize=rep(mean(dd$s_grpsize),2),
  focflat=c(1,1),
  focpreg=c(0,0),
  focswol=c(0,0)
)

d.pred.lact <- list(
  id1=c(1,1),
  id2=c(1,1),
  dyad=c(1,1),
  s_hrs=rep(mean(dd$s_hrs2),2),
  mom=rep(mean(dd$mom),2),
  sis=rep(mean(dd$sis),2),
  gma=rep(mean(dd$gma),2),
  aunt=rep(mean(dd$aunt),2),
  distant=rep(mean(dd$distant),2),
  enk=rep(mean(dd$enk),2),
  ynt=rep(mean(dd$ynt),2),
  focelo=rep(mean(dd$s_focelo),2),
  prtelo=rep(mean(dd$s_prtelo),2),
  cond=rep(mean(dd$s_foccondRev),2),
  grpsize=rep(mean(dd$s_grpsize),2),
  focflat=c(0,0),
  focpreg=c(0,0),
  focswol=c(0,0)
)

d.pred.preg <- list(
  id1=c(1,1),
  id2=c(1,1),
  dyad=c(1,1),
  s_hrs=rep(mean(dd$s_hrs2),2),
  mom=rep(mean(dd$mom),2),
  sis=rep(mean(dd$sis),2),
  gma=rep(mean(dd$gma),2),
  aunt=rep(mean(dd$aunt),2),
  distant=rep(mean(dd$distant),2),
  enk=rep(mean(dd$enk),2),
  ynt=rep(mean(dd$ynt),2),
  focelo=rep(mean(dd$s_focelo),2),
  prtelo=rep(mean(dd$s_prtelo),2),
  cond=rep(mean(dd$s_foccondRev),2),
  grpsize=rep(mean(dd$s_grpsize),2),
  focflat=c(0,0),
  focpreg=c(1,1),
  focswol=c(0,0)
)

d.pred.swol <- list(
  id1=c(1,1),
  id2=c(1,1),
  dyad=c(1,1),
  s_hrs=rep(mean(dd$s_hrs2),2),
  mom=rep(mean(dd$mom),2),
  sis=rep(mean(dd$sis),2),
  gma=rep(mean(dd$gma),2),
  aunt=rep(mean(dd$aunt),2),
  distant=rep(mean(dd$distant),2),
  enk=rep(mean(dd$enk),2),
  ynt=rep(mean(dd$ynt),2),
  focelo=rep(mean(dd$s_focelo),2),
  prtelo=rep(mean(dd$s_prtelo),2),
  cond=rep(mean(dd$s_foccondRev),2),
  grpsize=rep(mean(dd$s_grpsize),2),
  focflat=c(0,0),
  focpreg=c(0,0),
  focswol=c(1,1)
)

link2flat <- link(Food_Model_Full , n=1000 , data=d.pred.flat, replace=
                    list(ap_id=a_foc_z , ap_dyad=a_dyad_z, al_id=a_foc_z , al_dyad=a_dyad_z), WAIC=TRUE)
pred.flat <- link2flat$lambda*(1-link2flat$p)
pred.flat.p <- (1-link2flat$p)
pred.flat.lambda <- link2flat$lambda
median(pred.flat[,1])
HPDI(pred.flat[,1])

link2lact <- link(Food_Model_Full , n=1000 , data=d.pred.lact, replace=
                    list(ap_id=a_foc_z , ap_dyad=a_dyad_z, al_id=a_foc_z , al_dyad=a_dyad_z), WAIC=TRUE)
pred.lact <- link2lact$lambda*(1-link2lact$p)
pred.lact.p <- (1-link2lact$p)
pred.lact.lambda <- link2lact$lambda
median(pred.lact[,1])
HPDI(pred.lact[,1])

link2preg <- link(Food_Model_Full , n=1000 , data=d.pred.preg, replace=
                    list(ap_id=a_foc_z , ap_dyad=a_dyad_z, al_id=a_foc_z , al_dyad=a_dyad_z), WAIC=TRUE)
pred.preg <- link2preg$lambda*(1-link2preg$p)
pred.preg.p <- (1-link2preg$p)
pred.preg.lambda <- link2preg$lambda
median(pred.preg[,1])
HPDI(pred.preg[,1])

link2swol <- link(Food_Model_Full , n=1000 , data=d.pred.swol, replace=
                    list(ap_id=a_foc_z , ap_dyad=a_dyad_z,al_id=a_foc_z , al_dyad=a_dyad_z), WAIC=TRUE)
pred.swol <- link2swol$lambda*(1-link2swol$p)
pred.swol.p <- (1-link2swol$p)
pred.swol.lambda <- link2swol$lambda
median(pred.swol[,1])
HPDI(pred.swol[,1])

par(mfrow=c(1,1), cex=1, mar=c(0,0,0,0), oma=c(2,2,2,2))

dens(pred.lact[,1], xlim=c(0,.2) , ylim=c(-.035,40) , xlab="" , ylab="" , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE , cex.lab=1.5 , cex.axis=0.5 , adj=0.1)
ll <- dd$aggr1[dd$foclact==1]
points(ll, rep(-.025,length(ll)), pch="|" , col=col.alpha("#33CCFF", alpha=0.5) , cex=0.5 )
shade( density(pred.lact[,1]) , lim= as.vector(HPDI(pred.lact[,1], prob=0.99999)) , col = col.alpha("#33CCFF", 0.5))
shade( density(pred.preg[,1]) , lim= as.vector(HPDI(pred.preg[,1], prob=0.99999)) , col = col.alpha("orange1", 0.5))
shade( density(pred.swol[,1]) , lim= as.vector(HPDI(pred.swol[,1], prob=0.99999)) , col = col.alpha("green", 0.5))
shade( density(pred.flat[,1]) , lim= as.vector(HPDI(pred.flat[,1], prob=0.99999)) , col = col.alpha("purple", 0.5))
ll <- dd$aggr1[dd$focpreg==1]
points(ll, rep(-.025,length(ll)), pch="|" , col=col.alpha("orange1", alpha=0.5) , cex=0.5 )
ll <- dd$aggr1[dd$focswol==1]
points(ll, rep(-.025,length(ll)), pch="|" , col=col.alpha("green", alpha=0.5) , cex=0.5 )
ll <- dd$aggr1[dd$foclat==1]
points(ll, rep(-.025,length(ll)), pch="|" , col=col.alpha("purple", alpha=0.5) , cex=0.5 )

abline(v=median(pred.lact[,1]) , lty=1)
abline(v=median(pred.preg[,1]) , lty=2)
abline(v=median(pred.swol[,1]) , lty=3)
abline(v=median(pred.flat[,1]) , lty=4)

legend(.13,40, legend = c("", "", "", ""),
       col=c(1,1)  , lty=c(1,2,3,4),
       lw=1 , cex=.85, bty="n")

legend(.14,40,, legend = c("lactating", "pregnant", "swollen", "flat"),
       col=c(col.alpha("#33CCFF", 0.5) , col.alpha("orange1", 0.5) ,  col.alpha("green", 0.5) , col.alpha("purple", 0.5)) , pch=c(15,17), pt.cex=2 , cex=.85, bty="n")

mtext("Rate of aggression given" , side=1 , line=.3, outer=TRUE , cex=1.5)
#axis(1, at = seq(from=0 , to=1, by = .1) ,tck=-0.1 , labels=T )
axis(1, at = seq(from=0 , to=1, by = .1) ,tck=-0.01 , labels=T )

