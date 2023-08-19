###################################################################
# This function is used to calculate gradient and heissan matrices
#        in the process of likelihood maximization
###################################################################
LL_ha <-function(X_sub,Y_s, betas, sigs, rhos, pairlist,outcome_type){
  n=N=nrow(X_sub)
  q=ncol(Y_s)
  p=ncol(X_sub)
  la=la2=la3=0

  for(j in 1:(q-1)){   #ncol(Ys): # of outcomes
    for(k in (j+1):q){
        data_1 <- cbind(1,X_sub,Y_s[,j]) #1 + X +outcome j
        data_2 <- cbind(1,X_sub,Y_s[,k]) #1 + X +outcome k
        beta_1 <- betas[j,]
        beta_2 <- betas[k,]
        sig_1 <- sigs[j]
        sig_2 <- sigs[k]
        rho_12 <- rhos[j,k]
        mup=pairfunlist[[pairlist[j,k]]](data_1=data_1, data_2=data_2, beta_1=beta_1, beta_2=beta_2, sig_1=sig_1, sig_2=sig_2, rho_12=rho_12) #binary first
        if (length(mup)==0) return()
        attr(la,paste("l_",j,k,"_grad",sep=""))=attr(mup,"gradient")
        attr(la2,paste("l_",j,k,"_hessian",sep=""))=attr(mup,"hessian")
        attr(la3,paste("l_",j,k,"_gradpoint",sep=""))=attr(mup,"grad_point")
      }
    }


  l_beta=vector(mode = "list", length = q)
  l_beta=lapply(l_beta, function(x) x=0)
  pl_beta=pl_sig=l_beta

  # l_beta=pl_beta=pl_sig=list()   #pl: point
  # for(j in 1:ncol(Y)){
  #   l_beta[[j]]=pl_beta[[j]]=pl_sig[[j]]=0
  # }
  l_sig=rep(0,q)
  l_rho=pl_rho=matrix(0,nrow=q,ncol=q)
  diag(l_rho)=1
  su=1+sqrt(1/N/2)
  for(j in 1:(q-1)){
    for(k in (j+1):q){
      if(pairlist[j,k]=="Case_2_discrete"){
        ca=attr(la,paste("l_",j,k,"_grad",sep=""))
        l_beta[[j]]=l_beta[[j]]+ca[1:(p+1)]
        l_beta[[k]]=l_beta[[k]]+ca[(p+2):(2*p+2)]
        l_rho[j,k]=l_rho[k,j]=ca[length(ca)]
        ca=attr(la3,paste("l_",j,k,"_gradpoint",sep=""))
        pl_beta[[j]]=pl_beta[[j]]+ca[,1:(p+1)]
        pl_beta[[k]]=pl_beta[[k]]+ca[,(p+2):(2*p+2)]
        attr(pl_rho,paste(j,k,sep=""))=ca[,dim(ca)[2]]
      }
      if(pairlist[j,k]=="Case_2_continuous"){
        ca=attr(la,paste("l_",j,k,"_grad",sep=""))
        l_beta[[j]]=l_beta[[j]]+ca[1:(p+1)]
        l_beta[[k]]=l_beta[[k]]+ca[(p+2):(2*p+2)]
        l_sig[j]=l_sig[j]+ca[2*p+3]
        l_sig[k]=l_sig[k]+ca[2*p+4]
        l_rho[j,k]=l_rho[k,j]=ca[length(ca)]

        ca=attr(la3,paste("l_",j,k,"_gradpoint",sep=""))
        pl_beta[[j]]=pl_beta[[j]]+ca[,1:(p+1)]
        pl_beta[[k]]=pl_beta[[k]]+ca[,(p+2):(2*p+2)]
        pl_sig[[j]]=pl_sig[[j]]+ca[,(2*p+3)]
        pl_sig[[k]]=pl_sig[[k]]+ca[,(2*p+4)]
        attr(pl_rho,paste(j,k,sep=""))=ca[,dim(ca)[2]]
      }
      if(pairlist[j,k]=="Case_1_discrete"){
        ca=attr(la,paste("l_",j,k,"_grad",sep=""))   #bin_cont:  input: bin, cont   output: cont,bin
        l_beta[[j]]=l_beta[[j]]+ca[(p+2):(2*p+2)]    #
        l_beta[[k]]=l_beta[[k]]+ca[1:(p+1)]    #
        l_sig[k]=l_sig[k]+ca[2*p+3]
        l_rho[j,k]=l_rho[k,j]=ca[length(ca)]

        ca=attr(la3,paste("l_",j,k,"_gradpoint",sep=""))
        pl_beta[[j]]=pl_beta[[j]]+ca[,(p+2):(2*p+2)]
        pl_beta[[k]]=pl_beta[[k]]+ca[,1:(p+1)]
        pl_sig[[k]]=pl_sig[[k]]+ca[,(2*p+3)]
        attr(pl_rho,paste(j,k,sep=""))=ca[,ncol(ca)]
      }
    }
  }

  l_gradient = list(l_beta,l_sig,l_rho)
  l_gradient_point = list(pl_beta,pl_sig,pl_rho)

  name_Y=colnames(Y_s)
  name_beta=c("Intercept",colnames(X_sub))
  v.name = unlist(lapply(name_Y,function(x) paste( rep(x,length(name_beta)), name_beta, sep = " " )))
  name_sig=paste("sig",1:q,sep="")
  name_rho=c()
  for(j in 1:(q-1)){
    for(k in (j+1):q){
      name_rho=c(name_rho,paste("rho",j,k,sep=""))
    }
  }
  name=c(v.name,name_sig,name_rho)


  #Hessian
  ll_hessian <- matrix(0,length(name),length(name))
  colnames(ll_hessian)=row.names(ll_hessian)=name

  for(j in 1:(q-1)){
    for(k in (j+1):q){
      yesindex=c()
      if(pairlist[j,k]=="Case_2_discrete"){
        ca=attr(la2,paste("l_",j,k,"_hessian",sep=""))
        yesnames=c(paste(rep(name_Y[j],length(name_beta)), name_beta, sep = " " ),paste(rep(name_Y[k],length(name_beta)), name_beta, sep = " "))
        yesnames=c(yesnames,paste("rho",j,k,sep=""))
        #yesnames=c(paste("intercept",j,sep=""),paste("beta",j,sep=""),paste("intercept",k,sep=""),paste("beta",k,sep=""),paste("rho",j,k,sep=""))
        for(m in 1:length(yesnames)){
          yesindex[m]=which(name==yesnames[m])
        }
        ll_hessian[yesindex,yesindex]=ll_hessian[yesindex,yesindex]+ca
      }
      if(pairlist[j,k]=="Case_2_continuous"){
        ca=attr(la2,paste("l_",j,k,"_hessian",sep=""))
        yesnames=c(paste(rep(name_Y[j],length(name_beta)), name_beta, sep = " " ),paste(rep(name_Y[k],length(name_beta)), name_beta, sep = " "))
        yesnames=c(yesnames,paste("sig",j,sep=""),paste("sig",k,sep=""),paste("rho",j,k,sep=""))
        #yesnames=c(paste("intercept",j,sep=""),paste("beta",j,sep=""),paste("intercept",k,sep=""),paste("beta",k,sep=""),paste("sig",j,sep=""),paste("sig",k,sep=""),paste("rho",j,k,sep=""))
        for(m in 1:length(yesnames)){
          yesindex[m]=which(name==yesnames[m])
        }
        ll_hessian[yesindex,yesindex]=ll_hessian[yesindex,yesindex]+ca
      }
      if(pairlist[j,k]=="Case_1_discrete"){
        ca=attr(la2,paste("l_",j,k,"_hessian",sep=""))

        yesnames=c(paste(rep(name_Y[k],length(name_beta)), name_beta, sep = " " ),paste(rep(name_Y[j],length(name_beta)), name_beta, sep = " "))
        yesnames=c(yesnames,paste("sig",k,sep=""),paste("rho",j,k,sep=""))
        #yesnames=c(paste("intercept",k,sep=""),paste("beta",k,sep=""),paste("intercept",j,sep=""),paste("beta",j,sep=""),paste("sig",k,sep=""),paste("rho",j,k,sep=""))
        for(m in 1:length(yesnames)){
          yesindex[m]=which(name==yesnames[m])
        }
        ll_hessian[yesindex,yesindex]=ll_hessian[yesindex,yesindex]+ca
      }
    }
  }



  l_gradient1=c()
  for(i in 1:q){
    l_gradient1=c(l_gradient1,l_gradient[[1]][[i]])
  }
  l_gradient1=c(l_gradient1,l_gradient[[2]])
  for(j in 1:(q-1)){
    for(k in (j+1):q){
      l_gradient1=c(l_gradient1,l_gradient[[3]][j,k])
    }
  }
  names(l_gradient1)=name


  pl_gradient1=matrix(nrow=nrow(X_sub),ncol=length(name))
  for(i in 1:q){
    pl_gradient1[,((p+1)*(i-1)+1):((p+1)*i)]=l_gradient_point[[1]][[i]]
    pl_gradient1[,(p+1)*q+i]=l_gradient_point[[2]][[i]]
  }
  for(j in 1:(q-1)){
    for(k in (j+1):q){
      lao=which(name==paste("rho",j,k,sep=""))
      pl_gradient1[,lao]=attr(l_gradient_point[[3]],paste(j,k,sep=""))
    }
  }
  colnames(pl_gradient1)=name


  ind=which(outcome_type=="binary")

  ll=list(l_gradient1,pl_gradient1,ll_hessian)
  if (any(is.na(ll_hessian))==T) return()

  #pl_nonzero<-ll[[2]][-((p+1)*q+ind),-((p+1)*q+ind)]
  #hessian_nonzero<-ll[[3]][-((p+1)*q+ind),-((p+1)*q+ind)]
  J_theta   <-  1/N * (t ( ll[[2]])  %*% ll[[2]] )
  H_theta   <- - ll[[3]]/N
  I_Godambe <- H_theta %*%  ginv(J_theta)  %*%  H_theta
  I_Godambe[is.nan(I_Godambe)] <- 0
  cova=su/N*ginv(I_Godambe)
  mama=diag2vec(ginv(I_Godambe))
  #if (sum(is.na(mama))>0) print(c(i,mama))
  mama[mama<0]=0
  sd<-sqrt( 1/N*mama)
  
  ll=list(l_gradient1,pl_gradient1,ll_hessian,sd,cova)

  return(ll)
}
