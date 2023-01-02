###Mezclo 3 distribuciones, una gaussiana, una cauchy, 
rm(list=ls())
#######caso 1 Cambiando los grados de libertad
library(mvtnorm)
library(npmv)
library(energy)
library(SHT)
###########pvalor para el test normal 2010

pvalorn<-function(datos1,datos2){
  YY<-cbind(datos1,datos2)
  nn<-nrow(YY)
  d<-ncol(datos1)
  JJn<-matrix(1,ncol=nn,nrow=nn)
  IIn<-diag(nn)
  SSA<-t(YY)%*% (IIn-(1/nn)*JJn) %*% YY
  JJ2<-matrix(1,ncol=2,nrow=2)
  IIp<-diag(d)
  AA<-(1/2)*(JJ2%x% IIp)
  SSM<-t(YY-(1/nn)*JJn%*%YY%*%AA) %*% (YY-(1/nn)*JJn%*%YY%*%AA) 
  diag1<-SSM[1:d,1:d]
  diag2<-SSM[(d+1):(2*d),(d+1):(2*d)]
  diagprom<-(diag1+diag2)/2
  hSSM<- SSM
  hSSM[1:d,1:d]<-diagprom
  hSSM[(d+1):(2*d),(d+1):(2*d)]<-diagprom
  Wn<-nn*log(det(hSSM)/det(SSA))
  gr<-d*(d+3)/2
  print(1-pchisq(Wn,gr))
}



pvalorc<-function(datos1,datos2){  pp<-sim2.2018HN(datos1, datos2)

pp[[2]]
}

pvalore<-function(datos1,datos2){  
  pp<-eqdist.etest(rbind(datos1,datos2), sizes=c(nrow(datos1),nrow(datos2)),
                   method="original", R=1000)
  pp[[4]]
}


generadorDatos<-function(n,d,alfa,beta,mu, gr=grr){
  prob<-c(alfa, beta, 1-alfa-beta)
  ns<-rmultinom(1, n, prob)
  mus<-matrix(mu*rep(1,n*d),ncol=d,nrow=n)
  
if(ns[1]>0) {pp1<-matrix(mu*rep(1,ns[1]*d),ncol=d,nrow=ns[1])+ matrix(rnorm(d*ns[1]), ncol=d)} else  {pp1<-NULL}
  if(ns[2]>0) {pp2<-matrix(mu*rep(1,ns[2]*d),ncol=d,nrow=ns[2])+ rmvt(ns[2], sigma = diag(d), df = gr)} else  {pp2<-NULL}
  if(ns[3]>0) {
    nns1<-rbinom(1,ns[3], 0.5)
    nns2<-ns[3]-nns1
    nnnomal1<-matrix(1+mu+rnorm(d*nns1), ncol=d)
    nnnomal2<-matrix(-1-mu+rnorm(d*nns2), ncol=d)
    nonormal<-rbind(nnnomal1,nnnomal2)
    pp3<-matrix(mu*rep(1,ns[3]*d),ncol=d,nrow=ns[3])+ nonormal} else  {pp3<-NULL}
  
datos<-rbind(pp1,pp2,pp3)
print(datos)
}


n=100
#######proporcion de datos de al normal
alfa1<-(0:4)/4
beta1<-(0:4)/4
grilla<-expand.grid(alfa1, beta1)
grilla2<-round(grilla[apply(grilla,1,sum)<=1,],2)

potenciasF<-list()
for(ma in 1:nrow(grilla2)){
alfa=grilla2[ma,1]
#######proporcion de datos de la cauchy
beta=grilla2[ma,2]
grr=1
d<-5
datos1<-generadorDatos(n,d,alfa,beta,mu=0)
######
potencias<-list()
mmu2<-c(0,0.2,0.4,0.6)
for(ja in 1:length(mmu2)){
mu2<-mmu2[ja]
##############
datos2<-generadorDatos(n,d,alfa,beta,mu=mu2)
#plot(normal)
#lines(cauchy, type="p", col="red")
#lines(nonormal, type="p", col="blue")

##################Calculo del estadistico de RP  bajo H0 cierta
statdist<-c()
for(iis in 1:10000){
  datos<-rbind(datos1,datos2)
  datos1a<-datos[sample(1:nrow(datos), nrow(datos1), replace = T),]
  datos2b<-datos[sample(1:nrow(datos), nrow(datos2), replace = T),]
  direcciones<-matrix(rnorm(d*d), ncol=d)
  norma<-function(u){ sqrt(sum(u^2))}
  normaliza<-function(u){ u /norma(u)}
  direccionesn<-t(apply(direcciones,1,normaliza))
  #########vectores en los que tengo que proyectar (d^2+d)/2
  bb<-c()
  for(ii in 1:d){
    for(jj in 1:ii){ v<-direccionesn[ii,]+direccionesn[jj,]
    bb<-cbind(bb,normaliza(v))
    }}
  proy<-ncol(bb)
  datos1aProy<-datos1a%*%bb
  datos2bProy<-datos2b%*%bb
  dens<-c()
  for(gg in 1:proy){ a<-ks.test(datos1aProy[,gg],datos2bProy[,gg])
  dens[gg]<-a$statistic
  }
  statdist[iis]<-max(dens)
}


potencia<-c()
potencia_normal<-c()
potencia_cov<-c()
potencia_ene<-c()
for(ss in 1:1000){
  ########GENERO DATOS 
  datos1<-generadorDatos(n,d,alfa,beta,mu=0)
  datos2<-generadorDatos(n,d,alfa,beta,mu=mu2)
  #######################################################################
  #######generacion de los datos S.
  ###d vectores independientes unitarios
  #direcciones<-matrix(rnorm(d*d), ncol=d)
  #norma<-function(u){ sqrt(sum(u^2))}
  #normaliza<-function(u){ u /norma(u)}
  #direccionesn<-t(apply(direcciones,1,normaliza))
  #########vectores en los que tengo que proyectar (d^2+d)/2
  #bb<-c()
  #for(ii in 1:d){
  #  for(jj in 1:ii){ v<-direccionesn[ii,]+direccionesn[jj,]
  #  bb<-cbind(bb,normaliza(v))
  #  }}
  
  #proy<-ncol(bb)
  datos1Proy<-datos1%*%bb
  datos2Proy<-datos2%*%bb
  den<-c()
  for(gg in 1:proy){
    den[gg]<-ks.test(datos1Proy[,gg],datos2Proy[,gg])$statistic
  }
  stat<-max(den)
  
  
  pvalor_normal<-pvalorn(datos1,datos2)
  potencia_normal[ss]<-pvalor_normal<0.05
  
  potencia_co<-pvalorc(datos1,datos2)
  potencia_cov[ss]<- potencia_co<0.05
  
  potencia_en<-pvalore(datos1,datos2)
  potencia_ene[ss]<- potencia_en<0.05
  
  pvalor<-mean(statdist>stat)
  
  potencia[ss]<-pvalor<0.05
  
  #datos<-data.frame(rbind(cbind(datos1,1), cbind(datos2,2)))
  
  #datos$X11<-as.factor(datos$X11)
  
  
  #dd<-nonpartest(X1|X2|X3|X4|X5|X6|X7|X8|X9|X10~X11,permreps=1000, test=c(1,0,0,0),datos,alpha=.05, plots=F ,permtest=T)
  
  #potencia2[ss]<-dd[[1]][[4]][1]<0.05
  #potencia3[ss]<- eqdist.etest(datos[, 1:d], c(n,n) , R=999)$p.value
  #################################
  #########otro test
  
}

potencias[[ja]]<-c(mean(potencia),
mean(potencia_normal, na.rm=T),
mean(potencia_cov, na.rm=T),
mean( potencia_ene, na.rm=T)
)
}
potenciasF[[ma]]<-do.call(rbind,potencias)
}



potenciasF

jaa=list()
for(ja in 1:nrow(grilla2)){
pp<-grilla2[ja,]
nombre<-paste0("(", pp[1],",",pp[2], ",",1-pp[1]-pp[2], ")")
gr<-data.frame(potenciasF[[ja]])
names(gr)<-c("RPT","LRTN","L2Norm","eDistance")
gr$mus<-mmu2
gr$nombres<-nombre
jaa[[ja]]<-gr
}

grr<-do.call(rbind,jaa)

library(ggplot2)
library(reshape2)
base2<-melt(grr, id=c("mus","nombres"))
names(base2)<-c("mus", "nombres","test", "Power")

unique(base2$nombres)

base2<-subset(base2, nombres %in% c("(1,0,0)","(0,1,0)","(0,0,1)","(0.5,0.5,0)","(0.25,0.75,0)","(0,0.75,0.25)","(0.25,0.5,0.25)", "(0.75,0.25,0)", "(0,0.5,0.5)"     ))

library(tidyverse)
library(showtext)
library(thematic)
colors <- thematic::okabe_ito(4)
base2$test<-as.factor(base2$test)



# Line pl
p <- base2 %>% 
  ggplot(aes(x = mus , y= Power, col = test)) +facet_wrap(~ nombres, nrow = 3, labeller = label_bquote(alpha == .(nombres)))+
  geom_line(size = 0.5) +
  scale_color_manual(values = colors) +theme_bw()+labs(x = ~ mu[2], y = 'Power function', col = 'Test')


p

  save.image("potenciasF1.Rdata")

pdf("Pp.pdf", height=6, width = 7)
p
dev.off()








