get_data<-function(df,spec){
  Y <- df[,"re78"]
  T <- df[,"treat"]
  
  ## Main effects and all squared terms
  X.1 <- cbind(df[,"age"], df[,"education"], df[,"married"], df[,"black"], df[,"hispanic"], df[,"re74"], df[,"re75"], df[,"age"]^2, df[,"education"]^2, df[,"re74"]^2, df[,"re75"]^2)
  
  ## Main effects and all squared terms and some indicators
  X.2 <- cbind(df[,"age"], df[,"education"], df[,"married"], df[,"black"], df[,"hispanic"], df[,"re74"], df[,"re75"], df[,"age"]^2, df[,"education"]^2, df[,"re74"]^2, df[,"re75"]^2,df[,"re74"]==0, df[,"re75"]==0,df[,"nodegree"])
  
  ## Many more terms for flexibility
  
  X.main <- cbind(df[,"age"], df[,"education"], df[,"married"], df[,"black"], df[,"hispanic"], df[,"re74"], df[,"re75"], df[,"nodegree"], df[,"re74"]==0, df[,"re75"]==0)
  
  X.cont.interactions <- cbind(
    df[,"age"]*df[,"married"],
    df[,"age"]*df[,"nodegree"],
    df[,"age"]*df[,"black"],
    df[,"age"]*df[,"hispanic"],
    df[,"age"]*(df[,"re74"]==0),
    df[,"age"]*(df[,"re75"]==0),
    df[,"education"]*df[,"married"],
    df[,"education"]*df[,"nodegree"],
    df[,"education"]*df[,"black"],
    df[,"education"]*df[,"hispanic"],
    df[,"education"]*(df[,"re74"]==0),
    df[,"education"]*(df[,"re75"]==0),
    df[,"re74"]*df[,"married"],
    df[,"re74"]*df[,"nodegree"],
    df[,"re74"]*df[,"black"],
    df[,"re74"]*df[,"hispanic"],
    df[,"re74"]*(df[,"re75"]==0),
    df[,"re75"]*df[,"married"],
    df[,"re75"]*df[,"nodegree"],
    df[,"re75"]*df[,"black"],
    df[,"re75"]*df[,"hispanic"],
    df[,"re75"]*(df[,"re74"]==0)
  )
  
  X.dummy.interactions <- cbind(
    df[,"married"]*df[,"nodegree"],
    df[,"married"]*df[,"black"],
    df[,"married"]*df[,"hispanic"],
    df[,"married"]*(df[,"re74"]==0),
    df[,"married"]*(df[,"re75"]==0),
    df[,"nodegree"]*df[,"black"],
    df[,"nodegree"]*df[,"hispanic"],
    df[,"nodegree"]*(df[,"re74"]==0),
    df[,"nodegree"]*(df[,"re75"]==0),
    df[,"black"]*(df[,"re74"]==0),
    df[,"black"]*(df[,"re75"]==0),
    df[,"hispanic"]*(df[,"re74"]==0),
    df[,"hispanic"]*(df[,"re75"]==0),
    (df[,"re74"]==0)*(df[,"re75"]==0)
  )
  
  X.poly <- poly(cbind(df[,"age"], df[,"education"], df[,"re74"], df[,"re75"]),degree=5)
  
  X.educ.dummies <- model.matrix(~as.factor(df[,"education"]) - 1)[,-1]
  X.3<-cbind(X.main,X.cont.interactions,X.dummy.interactions,X.poly)
  
  if (spec==1){
    X=X.1
  } else if (spec==2) {
    X=X.2
  } else if (spec==3) {
    X=X.3
  }
  
  X <- scale(X,center=TRUE,scale=TRUE) #can set center=TRUE
  
  #impose common support
  p.1 <- multinom(T~X-1, trace=FALSE)$fitted.values
  indexes.to.drop <- which(p.1 < min(p.1[T==1]) | max(p.1[T==1]) < p.1)
  if (length(indexes.to.drop)==0) {indexes.to.drop <- n+1}	#R throws a wobbly if [-indexes.to.drop] is negating an empty set. 
  n.per.treatment <- as.vector(table(T[-indexes.to.drop]))
  n.trim <- n.per.treatment[1]+n.per.treatment[2]
  
  Y.trimmed=Y[-indexes.to.drop]
  T.trimmed=T[-indexes.to.drop]
  X.trimmed=X[-indexes.to.drop,]
  
  return(list(Y.trimmed,T.trimmed,X.trimmed))
  
}


