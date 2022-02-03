SLModels<-function(data, algorithm="stepwise", scaling=FALSE){
  ##Scaling
  if(scaling==TRUE){
    MinMaxScaling <- function(x){return((x-min(x))/(max(x)-min(x)))}
    for(i in 1:(dim(data)[2]-1)){
      data[,i]<-MinMaxScaling(data[,i])
    }
  }
  if(algorithm=="stepwise" || algorithm=="minmax" || algorithm=="minmaxmedian" || algorithm=="minmaxiqr"){
    ##Inputs transformation
    if(algorithm=="minmax"){
      transf<-cbind(apply(data[,-c(dim(data)[2])],1,max),apply(data[,-c(dim(data)[2])],1,min))
      data<-cbind(transf,data[,dim(data)[2]])
      data<-as.data.frame(data)
    }else if(algorithm=="minmaxmedian"){
      transf<-cbind(apply(data[,-c(dim(data)[2])],1,max),apply(data[,-c(dim(data)[2])],1,min),apply(data[,-c(dim(data)[2])],1,median))
      data<-cbind(transf,data[,dim(data)[2]])
      data<-as.data.frame(data)
    }else if(algorithm=="minmaxiqr"){
      transf<-cbind(apply(data[,-c(dim(data)[2])],1,max),apply(data[,-c(dim(data)[2])],1,min),apply(data[,-c(dim(data)[2])],1,quantile)[3,]-apply(data[,-c(dim(data)[2])],1,quantile)[1,])
      data<-cbind(transf,data[,dim(data)[2]])
      data<-as.data.frame(data)
    }
    ##Algorithm
    indicesa<-expand.grid(seq(1,length(data)-1,by=1),seq(1,length(data)-1,by=1))
    indicesp<-subset(indicesa,indicesa$Var1!=indicesa$Var2)
    indicesn<-indicesp
    indicesn[,2]<- -indicesp[,2]
    indices<-rbind(indicesp,indicesn)
    lineal2varp<-function(x,y,k){data[,x]+k*data[,y]}
    SLMcombination2ap<-function(p,q) {mapply(lineal2varp,x=p,y=q,k=seq(-1,1,by=0.01))}
    SLMcombination2p<-matrix(mapply(SLMcombination2ap,p=indicesp[,2],q=indicesp[,1]),nrow=dim(data)[1])
    lineal2varn<-function(x,y,k){-data[,x]+k*data[,y]}
    SLMcombination2an<-function(p,q) {mapply(lineal2varn,x=p,y=q,k=seq(-1,1,by=0.01))}
    SLMcombination2n<-matrix(mapply(SLMcombination2an,p=-indicesn[,2],q=indicesn[,1]),nrow=dim(data)[1])
    SLMcombination2<-cbind(SLMcombination2p,SLMcombination2n)
    t<-matrix(cbind(rep(data[,length(data)],c(2*201*(length(data)-2)*(length(data)-1)))),nrow=dim(data)[1])
    youden2varf<-function(x){
      p<-prediction(SLMcombination2[,x],data[,length(data)])
      a<-attributes(performance(p,"sens","spec"))
      max(a$y.values[[1]]+a$x.values[[1]]-1)
    }
    youden2var<-sapply(seq(1,2*201*(length(data)-1)*(length(data)-2)),youden2varf)
    validation2var<-function(x){
      which(youden2var==max(youden2var))}
    Youdenposition2var<-validation2var(1)
    beta1<-((Youdenposition2var-(floor((Youdenposition2var-1)/201))*201)-101)/100
    v1<-indices[floor((Youdenposition2var-1)/201)+1,2]
    v2<-indices[floor((Youdenposition2var-1)/201)+1,1]
    datareducedf<-function(x){data[,c(subset(seq(1,length(data)-1,by=1),seq(1,length(data)-1,by=1)!=abs(v1[x])&seq(1,length(data)-1,by=1)!=v2[x]),length(data))]}
    datareduced<-lapply(seq(1,length(beta1),by=1),datareducedf)
    indicesselectedf<-function(x){subset(indices,indices$Var2==v1[x]&indices$Var1!=v2[x])}
    indicesselected<-lapply(seq(1,length(beta1),by=1),indicesselectedf)
    dataselectedf<-function(x){data.frame(SLMcombination2[,Youdenposition2var[x]],datareduced[[x]])}
    dataselected<-lapply(seq(1,length(beta1),by=1),dataselectedf)
    Youden<-youden2var[Youdenposition2var]
    informationstep1<-data.frame(beta1,v1,v2,Youden)
    information<-list(informationstep1)
    betanextstep2<-beta1
    marcador<-c(2)
    while(marcador<length(data)-1){
      combination<-function(h){
        linealvar<-function(y,k){dataselected[[h]][,1]+k*data[,y]}
        SLMcombinationa<-function(q) {mapply(linealvar,y=q,k=seq(-1,1,by=0.01))}
        SLMcombination<-matrix(mapply(SLMcombinationa,q=indicesselected[[h]][,1]),nrow=dim(dataselected[[h]])[1])
        tselected<-matrix(rep(data[,length(data)],dim(SLMcombination)[2]),nrow=dim(data)[1])
        youdenvarselectedf<-function(x){
          p<-prediction(SLMcombination[,x],data[,length(data)])
          a<-attributes(performance(p,"sens","spec"))
          max(a$y.values[[1]]+a$x.values[[1]]-1)
        }
        youdenvarselected<-sapply(seq(1,201*(dim(dataselected[[1]])[2]-2)),youdenvarselectedf)
        c(which(youdenvarselected==max(youdenvarselected)),max(youdenvarselected))
      }
      Youdenpositionvarselected<-lapply(seq(1,length(betanextstep2),by=1),combination)
      selectYoudenf<-function(x){Youdenpositionvarselected[[x]][length(Youdenpositionvarselected[[x]])]}
      Youdenvalueselected<-sapply(seq(1,length(Youdenpositionvarselected),by=1),selectYoudenf)
      selectYouden<-max(sapply(seq(1,length(Youdenpositionvarselected),by=1),selectYoudenf))
      positionlist<-which(Youdenvalueselected==selectYouden)
      Youdenpositionvarselectedstep1f<-function(x){Youdenpositionvarselected[[x]]}
      Youdenpositionvarselectedstep1<-lapply(positionlist,Youdenpositionvarselectedstep1f)
      indicesselectedstep1f<-function(x){indicesselected[[x]]}
      indicesselectedstep1<-lapply(positionlist,indicesselectedstep1f)
      Youdenpositionvarselectedstep2f<-function(x){Youdenpositionvarselectedstep1[[x]][1:(length(Youdenpositionvarselectedstep1[[x]])-1)]}
      Youdenpositionvarselectedstep2<-lapply(seq(1,length(Youdenpositionvarselectedstep1),by=1),Youdenpositionvarselectedstep2f)
      betanextstep1f<-function(x){((Youdenpositionvarselectedstep2[[x]]-(floor((Youdenpositionvarselectedstep2[[x]]-1)/201))*201)-101)/100}
      betanextstep1<-lapply(seq(1,length(Youdenpositionvarselectedstep2),by=1),betanextstep1f)
      vnextstep1f<-function(x){indicesselectedstep1[[x]][floor((Youdenpositionvarselectedstep2[[x]]-1)/201)+1,1]}
      vnextstep1<-lapply(seq(1,length(Youdenpositionvarselectedstep2),by=1),vnextstep1f)
      indicelistas1f<-function(x){length(betanextstep1[[x]])}
      indicelistas1<-sapply(seq(1,length(betanextstep1),by=1),indicelistas1f)
      indicelistas2<-rep(seq(1,length(betanextstep1),by=1),indicelistas1)
      indicelistas3fa<-function(x){seq(1,indicelistas1[x],by=1)}
      indicelistas3a<-sapply(seq(1,length(betanextstep1),by=1),indicelistas3fa)
      indicelistas3<-unlist(indicelistas3a)
      betanextstep2f<-function(x,y){as.list(betanextstep1[[x]][y])}
      betanextstep2<-mapply(betanextstep2f,x=indicelistas2,y=indicelistas3)
      vnextstep2f<-function(x,y){as.list(vnextstep1[[x]][y])}
      vnextstep2<-mapply(vnextstep2f,x=indicelistas2,y=indicelistas3)
      indicesselectedstep2f<-function(x){indicesselectedstep1[[x]]}
      indicesselectedstep2<-lapply(indicelistas2,indicesselectedstep2f)
      indicesselectedstep3f<-function(x){subset(indicesselectedstep2[[x]],indicesselectedstep2[[x]]$Var1!=vnextstep2[[x]])}
      indicesselectedstep3<-lapply(seq(1,length(betanextstep2),by=1),indicesselectedstep3f)
      dataselectedstep1f<-function(x){dataselected[[x]]}
      dataselectedstep1<-lapply(positionlist,dataselectedstep1f)
      dataselectedstepf<-function(x){dataselectedstep1[[x]]}
      dataselectedstep<-lapply(indicelistas2,dataselectedstepf)
      datareducedf<-function(x){data[c(indicesselectedstep3[[x]]$Var1,length(data))]}
      datareduced<-lapply(seq(1,length(betanextstep2),by=1),datareducedf)
      dataselectedf<-function(x){data.frame(dataselectedstep[[x]][,1]+betanextstep2[[x]]*data[,vnextstep2[[x]]],datareduced[[x]])}
      dataselected<-lapply(seq(1,length(betanextstep2),by=1),dataselectedf)
      indicesselected<-indicesselectedstep3
      Youdenselectedstepf<-function(x){
        p<-prediction(dataselectedstep[[x]][,1]+ betanextstep2[[x]]*data[,vnextstep2[[x]]],data[,length(data)])
        a<-attributes(performance(p,"sens","spec"))
        max(a$y.values[[1]]+a$x.values[[1]]-1)
      }
      Youdenselected<-lapply(seq(1,length(betanextstep2),by=1),Youdenselectedstepf)
      positionrepf<-function(x){rep(positionlist[x],length(vnextstep1[[x]]))}
      positionrep<-unlist(lapply(seq(1,length(vnextstep1),by=1),positionrepf))
      information[[length(information)+1]]<-data.frame(unlist(positionrep),unlist(betanextstep2),unlist(vnextstep2),unlist(Youdenselected))
      names(information[[length(information)]])<-c("last position","beta","var","Youden")
      marcador<-marcador+1
    }
    ##Final information
    finalinfo<-information
    marcadora<-c(length(information))
    while(marcadora>1){finalinfo[[marcadora-1]]<-information[[marcadora-1]][as.numeric(as.vector(as.data.frame(table(finalinfo[[marcadora]][,1]))[,1])),]
    marcadora<-marcadora-1}
    columnas<-c()
    filas<-c()
    for(i in 1:(dim(data)[2]-1)){
      if(i==1){
        columnas<-c(columnas,paste0("v",abs(as.numeric(finalinfo[[1]][1,][2]))))
        filas<-c(filas,sign(as.numeric(finalinfo[[1]][1,][2])))
      }else if(i==2){
        columnas<-c(columnas,paste0("v",abs(as.numeric(finalinfo[[1]][1,][3]))))
        filas<-c(filas,sign(as.numeric(finalinfo[[1]][1,][3]))*as.numeric(finalinfo[[1]][1,][1]))
      }else{
        columnas<-c(columnas,paste0("v",abs(as.numeric(finalinfo[[i-1]][1,][3]))))
        filas<-c(filas,sign(as.numeric(finalinfo[[i-1]][1,][3]))*as.numeric(finalinfo[[i-1]][1,][2]))
      }
    }
    columnas<-c(columnas,"Youden")
    filas<-c(filas,as.numeric(finalinfo[[dim(data)[2]-2]][1,][4]))
    filas<-filas[order(columnas)]
    columnas<-columnas[order(columnas)]
    comb<-0
    for(i in 1:(length(columnas)-1)){
      comb<-comb+filas[i]*data[,i]
    }
    data2<-matrix(comb)
    p<-prediction(data2,data[,dim(data)[2]])
    a<-attributes(performance(p,"sens","spec"))
    cutoff<-a$alpha.values[[1]][which.max(a$y.values[[1]]+a$x.values[[1]]-1)]
    columnas<-c(columnas,"Cutoff")
    filas<-c(filas,cutoff)
    if(algorithm=="minmax"){
      columnas<-replace(columnas, columnas=='v1','max')
      columnas<-replace(columnas, columnas=='v2','min')
    }else if(algorithm=="minmaxmedian"){
      columnas<-replace(columnas, columnas=='v1','max')
      columnas<-replace(columnas, columnas=='v2','min')
      columnas<-replace(columnas, columnas=='v3','median')
    }else if(algorithm=="minmaxiqr"){
      columnas<-replace(columnas, columnas=='v1','max')
      columnas<-replace(columnas, columnas=='v2','min')
      columnas<-replace(columnas, columnas=='v3','iqr')
    }
    output <- data.frame(matrix(nrow = 1,data = filas))
    colnames(output) <- columnas
    output

  }else{
    message('Argument \'algorithm\' is not correct. Available options: stepwise; minmax; minmaxmedian; minmaxiqr.')
  }
  }

