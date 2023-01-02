rm(list=ls())
#######caso 1 Cambiando los grados de libertad
library(mvtnorm)
###dimension


MUESTRAS<-c(100,500,1000)
dimensiones<-c(5,10,50)
gr2<-c(1,2,3,4)
MM<-list()
for(tam in 1:length(MUESTRAS)){ print(tam)
  sim=MUESTRAS[tam]
  dimesioneslist<-list()
  for(lf in 1:length(dimensiones)){
    d<-dimensiones[lf]
    
    df1=1
    
    potenciass<-c()
    for(mw in 1 : length(gr2)){
      mu<-0
      df2=gr2[mw]
      mus<-matrix(mu*rep(1,sim*d),ncol=d,nrow=sim)
      
      datos1 <- rmvt(sim, sigma = diag(d), df = df1) # t_3(0, diag(2)) sample
      datos2 <- mus+ rmvt(sim, sigma = diag(d), df = df2) # t_3(0, diag(2)) sample
      
      #statiticsdistributionH0
      statdist<-c()
      for(iis in 1:1000){
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
      for(ss in 1:1000){ 
        datos1 <- rmvt(sim, sigma = diag(d), df = df1) # t_3(0, diag(2)) sample
        datos2 <- mus+ rmvt(sim, sigma = diag(d), df = df2) # t_3(0, diag(2)) sample
        #######################################################################
        #######generacion de los datos S.
        ###d vectores independientes unitarios
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
        datos1Proy<-datos1%*%bb
        datos2Proy<-datos2%*%bb
        den<-c()
        for(gg in 1:proy){
          den[gg]<-ks.test(datos1Proy[,gg],datos2Proy[,gg])$statistic
        }
        stat<-max(den)
        
        
        
        
        
        pvalor<-mean(statdist>stat)
        
        potencia[ss]<-pvalor<0.05
        
      }
      
      potenciass[mw]<-mean(potencia)
    }
    dimesioneslist[[lf]]<-potenciass
  }
  
  resultados<-data.frame(do.call(cbind,dimesioneslist ))
  resultados$N<-sim
  
  MM[[tam]]<-resultados
  
}

base<-do.call(rbind,MM)
setwd("~/Dropbox/2022_CramerWold_Elliptical/CodeR")
save.image("student_gr.Rdata")

names(base)<-c("dim=5","dim=10","dim=50","N")
base$gr<-rep(gr2,length(MUESTRAS))
library(reshape2)

  
base2<-melt(base, id=c("gr","N"))
library(tidyverse)
library(showtext)
library(thematic)
colors <- thematic::okabe_ito(8)[-6]
base2$Dimension<-as.factor(base2$variable)
base2$N<-as.factor(base2$N)
base2$Power<-base2$value
# Line plot
p <- base2 %>% 
  ggplot(aes(x = gr , y= Power, col = N)) +
  geom_line(size = 0.5) + facet_wrap(~ Dimension, nrow = 1)+
  scale_color_manual(values = colors) +theme_bw()+
  labs(x = ~ paste("Degrees of freedom (", nu[2], ")"), y = 'Power function', col = 'Sample size')

setwd("~/Dropbox/2022_CramerWold_Elliptical/CodeR")
pdf("Power_Gr.pdf", height=3, width = 9)
p
dev.off()


 



