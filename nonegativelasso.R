A <- read.csv('x.csv')[,2:146] #自变量
y <- read.csv('y.csv')$V1 #因变量

error = 1e-8
want = 30 #需要留下的变量数
lambda = 0.00314
lambda_list = seq(0.003,0.004,0.000001) #lambda参数的范围
iter = 1000 #迭代次数
######以上为参数########
######以下为函数定义######
nonnagative.lasso <- function(X,y,lambda=0.5,iter=10000,error=1e-8){
  
  X <- as.matrix(X)
  y <- as.matrix(y)
  
  m <- nrow(X) #X行数
  n <- ncol(X) #X列数
  #argmin beta'(X'X)beta+(lambda_nI - 2X'Y)'beta
  #初始化beta
  beta <- as.matrix(rep(1,n))
  A <- t(X)%*%X
  A1 <- A
  A2 <- A
  A1[A1<0] <- 0
  A2[A2<0] <- -1*A2[A2<0]
  A2[A2<0]
  A2[A2>=0] <- 0
  b <- lambda*as.matrix(rep(1,n)) - 2*t(X)%*%y
  #开始进行迭代
  for(each_iter in 1:iter){
    a <- A1%*%beta #n*1向量
    #将缺失值设置为0
    c <- A2%*%beta #n*1向量
    beta0 <- beta*(-b+sqrt(b^2 + 4*a*c))/a
    beta0[is.na(beta0)] <- 0
    if(sum(abs(beta0-beta))<1e-8){  #控制误差
      return(beta0)
    }
    beta <- beta0
  }
  return(beta)
}

justify.lambda <- function(A,y,lambda_list,want,error=1e-8){
  print(paste('需要留下的非0变量数为:',want))
  nlambda <- length(lambda_list)
  n <- ncol(A)
  beta_list <- as.data.frame(matrix(0,nrow=nlambda,ncol=n+1))
  #set beta name
  beta.name <- rep(0,ncol(A))
  for(j in 1:n){
    beta.name[j] <- paste('beta',j)  
  }
  names(beta_list)[1] <- c('lambda')
  names(beta_list)[c(2:(n+1))] <- beta.name
  i <- 1
  
  for(l in lambda_list){
    beta <- nonnagative.lasso(A,y,l)
    if(sum(beta>error)==want){
      beta <- nonnagative.lasso(A,y,l)
      beta_list[i,] <- c(l,beta)
      i <- i + 1
    }
  }
  return(beta_list[1:(i-1),])
}

plot.lambda <- function(A,y,lambda_list,error=1e-8){
  nlambda <- length(lambda_list)
  n <- ncol(A)
  beta_list <- as.data.frame(matrix(0,nrow=nlambda,ncol=n+2))
  #set beta name
  beta.name <- rep(0,ncol(A))
  for(j in 1:n){
    beta.name[j] <- paste('beta',j)  
  }
  names(beta_list)[c(1:2)] <- c('lambda','not zero num')
  names(beta_list)[c(3:(n+2))] <- beta.name
  colormap <- rainbow(n)
  i <- 1
  
  for(l in lambda_list){
    beta <- nonnagative.lasso(A,y,l)
    num <- sum(beta>error)
    beta_list[i,] <- c(l,num,beta)
    i <- i + 1
  }
  
  plot(0,0,ylim=c(-.0001,1),xlim=c(min(lambda_list),max(lambda_list)),type='n',xlab = 'lambda',ylab = 'beta',main='variable selection')
  abline(0,0,col='black')
  for(j in 1:n){
    lines(lambda_list,beta_list[,j],col=colormap[j],lty=2)
  }
  
  legend('topright',title='beta',beta.name,lty=rep(2,n),col=colormap)
  
  return(beta_list)
}

######以下为主函数######
## 单个lambda计算
beta <- nonnagative.lasso(A,y,lambda,iter,error)
sum(beta>1e-8)
## 画图
bb <- plot.lambda(A,y,lambda_list,error)
tail(bb)
## 输出使多少个变量为0时的lambda值
justify.lambda(A,y,lambda_list,want,error)