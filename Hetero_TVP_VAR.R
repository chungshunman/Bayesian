# Heteroscedasticity VAR
library("dplyr") 

path= "/Users/manchungshun/Documents/GitHub/Bayesian/"

Y = (read.csv(paste(path, 'ydata.csv', sep=""),
                       header = FALSE))
yearlab = read.csv(paste(path, "yearlab.csv", sep=""),
                  header = FALSE)

t = dim(y)[1]; M = dim(y)[2]
tau = 40
p = 2
numa = M*(M-1)/2

# create lagged variables
ylag = Y
for(i in p:1){
  yl=lag(y, n=i)
  ylag = cbind(yl, ylag)
}
ylag = ylag[p+tau+1:t, 1:(M*p)]; rownames(ylag)=NULL


K = M + p*(M**2)
Z = matrix(0, nrow=(t-tau-p)*M, ncol=K)

for(i in 1:(t-tau-p)){
  ztemp = diag(M)
  for(j in 1:p){
    xtemp = as.matrix(ylag[i, ((j-1)*M+1):(j*M)])
    xtemp = kronecker(diag(M), xtemp)
    ztemp = cbind(ztemp, xtemp)
  }
  Z[((i-1)*M+1):(i*M), ] = ztemp
}

y = Y[(tau+p+1):t,]
yearlab = yearlab[(tau+p+1):t,]
t = length(yearlab)


# Preliminaries
nrep = 50
nburn = 20
it_print = 100

output = ts_prior(Y, tau, M, p)






