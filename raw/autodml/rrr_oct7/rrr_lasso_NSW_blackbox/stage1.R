arg_Forest<- list(clas_nodesize=1, reg_nodesize=5, ntree=1000, na.action=na.omit, replace=TRUE)
arg_Nnet3<- list(size=8,  maxit=1000, decay=0.01, MaxNWts=10000,  trace=FALSE)
arg_Nnet4<-  list(size=8,  maxit=1000, decay=0.01, MaxNWts=10000,  trace=FALSE, linout=TRUE)

get_stage1<-function(Y,T,X,p0,D_LB,D_add,max_iter,b,alpha_estimator,gamma_estimator,c_RR,c_CEF){
  
  p=length(b(T[1],X[1,]))
  n=length(T)
  MNG<-get_MNG(Y,T,X,b)
  M=MNG[[1]]
  N=MNG[[2]]
  G=MNG[[3]]
  B=MNG[[4]]
  
  ###########
  # alpha hat
  ###########
  if(alpha_estimator==0){ # dantzig with theoretical iteration
    rho_hat=RMD_stable(Y,T,X,p0,D_LB,D_add,max_iter,b,1,0,c_RR)
  } else if(alpha_estimator==1){ # lasso with theoretical iteration
    rho_hat=RMD_stable(Y,T,X,p0,D_LB,D_add,max_iter,b,1,1,c_RR)
  }
  
  alpha_hat<-function(d,z){
    return(b(d,z)%*%rho_hat)
  }
  
  ###########
  # gamma hat
  ###########
  if(gamma_estimator==0){ # dantzig with theoretical iteration
    
    beta_hat=RMD_stable(Y,T,X,p0,D_LB,D_add,max_iter,b,0,0,c_CEF)
    gamma_hat<-function(d,z){
      return(b(d,z)%*%beta_hat)
    }
    
  } else if(gamma_estimator==1){ # lasso with theoretical iteration
    
    beta_hat=RMD_stable(Y,T,X,p0,D_LB,D_add,max_iter,b,0,1,c_CEF)
    gamma_hat<-function(d,z){ 
      return(b(d,z)%*%beta_hat)
    }
    
  } else if(gamma_estimator==2){ # random forest
    
    forest<- do.call(randomForest, append(list(x=B,y=Y), arg_Forest))
    gamma_hat<-function(d,z){
      return(predict(forest,newdata=b(d,z), type="response"))
    }
    
  } else if(gamma_estimator==3 || gamma_estimator==4 || gamma_estimator==5){ # neural net
    
    # scale down, de-mean, run NN, scale up, remean so that NN works well
    # hack to ensure that constant covariates do not become NA in the scaling
    maxs_B <- apply(B, 2, max)
    mins_B <- apply(B, 2, min)
    maxs_Y<-max(Y)
    mins_Y<-min(Y)
    const=maxs_B==mins_B
    keep=(1-const)*1:length(const)
    NN_B<-B
    NN_B[,keep]<-scale(NN_B[,keep], center = mins_B[keep], scale = maxs_B[keep] - mins_B[keep])
    NN_Y<-scale(Y, center = mins_Y, scale = maxs_Y - mins_Y)
    
    # train
    if(gamma_estimator==3 || gamma_estimator==4){ # 1 layer NN (nnet)
      
      # choose architecture
      if(gamma_estimator==3){
        arg_Nnet=arg_Nnet3 # linout=FALSE
      } else {
        arg_Nnet=arg_Nnet4 # linout=TRUE
      }
      
      # use package
      nn<- do.call(nnet, append(list(x=NN_B,y=NN_Y), arg_Nnet))
      
    }else if(gamma_estimator==5){ # 2 layer NN (keras)
      
      # choose architecture
      build_model <- function() {
        model <- keras_model_sequential() %>% 
          layer_dense(units = 8, activation = "relu", 
                      input_shape = dim(NN_B)[[2]]) %>% 
          layer_dense(units = 8, activation = "relu") %>% 
          layer_dense(units = 1) 
        
        model %>% compile(
          optimizer = "rmsprop", 
          loss = "mse", 
          metrics = c("mae")
        )
      }
      
      # use package
      model <- build_model()
      num_epochs <- 100
      model %>% fit(NN_B, NN_Y,
                    epochs = num_epochs, batch_size = 1, verbose = 0)
    }
    
    # predict
    gamma_hat<-function(d,z){
      
      NN_b<-t(as.vector(b(d,z))) # test point
      NN_b[,keep]<-scale(t(NN_b[,keep]), 
                         center = mins_B[keep], 
                         scale = maxs_B[keep] - mins_B[keep])
      
      if(gamma_estimator==3 || gamma_estimator==4){ # 1 layer NN (nnet)
        NN_Y_hat<-predict(nn,newdata=NN_b)
      } else if(gamma_estimator==5){ # 2 layer NN (keras)
        NN_Y_hat<-model %>% predict(NN_b, verbose = 0)
      }
      Y_hat=NN_Y_hat*(maxs_Y-mins_Y)+mins_Y
      
      return(Y_hat)
    }
    
  }
  
  return(list(alpha_hat,gamma_hat))
  
}
