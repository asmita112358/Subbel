##Choice of K
rm(list = ls())
sim_K = function(alpha, beta, gamma, other_covariates, coef_X_alpha, coef_X_beta, 
                 n, omega = 0.05)
{
  ##Generating Data
  Y<- G <- S <- c()
  
  S = runif(n, 1,2)
  G = S*alpha + rowSums(coef_X_alpha*other_covariates) + rnorm(n)
  Y = G*beta + S*gamma + rowSums(coef_X_beta*other_covariates) + rnorm(n)
  ##Subsampling
  
  K = 5*(1:5)
  niter = 100
  
  tstat = matrix(nrow = length(K), ncol = niter)
  pval = matrix(nrow = length(K), ncol = niter)
  CCT = c()
  p_cct = c()
  for(kk in 1:length(K))
  {
    for(iter in 1:niter)
    {
      #Sample partitioning
      all <- sample.int(n = n, size = n, replace = F)
      
      index <- lapply(1:K[kk], function(ii){
        m <- floor(n/K[kk])
        all[((ii-1)*m+1):(ii*m)]
      })
      
      #identical(sort(do.call('c', index)), 1:n) #check the sample partition
      
      S_K = c()
      for(i in 1:K[kk])
      {
        
        G_s <- G[index[[i]]]
        S_s <- S[index[[i]]]
        Y_s <- Y[index[[i]]]
        other_covariates_s <- other_covariates[index[[i]],]
        
        obj1 = lm(G_s ~ -1 + cbind(S_s, other_covariates_s))
        obj2 = lm(Y_s ~ -1 + cbind(G_s, S_s, other_covariates_s))
        t1 = summary(obj1)[["coefficients"]][1,3]
        t2 = summary(obj2)[["coefficients"]][1,3]
        
        S_K[i] = (t1*t2)/sqrt(t1^2 + t2^2)
      }
      
      tstat[kk,iter] = sqrt(K[kk])*mean(S_K)/sqrt(var(S_K))
      pval[kk,iter] = 2*pt(-abs(tstat[kk,iter]), df = K[kk]-1)
      ##print(tstat)
    }
    wts = runif(niter)
    wts = wts/sum(wts)
    CCT[kk] = sum(tan(3.141593*(0.5 - pval[kk,]))*wts) ##Cauchy combination statistic
    p_cct[kk] = 1- pcauchy(CCT[kk])
  }
  
  return(p_cct)
  
}

##parameters

n = 600
#Coefficients
alpha = 0.5
beta = 0.5
gamma = 0.1
omega = 0.05
coef_X_alpha = c(1,1,1)
coef_X_beta = c(1,1,1)


##Checking size, qqplot
prop = c(0.8, 0.1, 0.1)
N_rep = 200       ##Number of monte carlo repetitions
mc = matrix(nrow = N_rep, ncol = 5)
samp = sample(1:3, N_rep, replace = T, prob = prop)
for(j in 1: N_rep)
{
  other_covariates = cbind(rep(1,n), runif(n), runif(n))
  alphaj = 0 + alpha*(samp[j]== 2)
  betaj = 0 + beta*(samp[j] == 3)
  #print(alphaj*betaj)
  mc[j,] = sim_K(alphaj, betaj,  gamma = 0.1, other_covariates, coef_X_alpha, coef_X_beta, 
               n, omega = 0.05)
}
size = apply((mc<=omega), 2, mean)
se_size = sqrt(size*(1 - size)/ N_rep)
K = 5*(1:5)
data = data.frame(K, size, se_size)
sizep <- data.table(data) %>%
  ggplot(aes(x = factor(K), y = size)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = size - 1.96*se_size, ymax = size + 1.96*se_size), width = 0.2,position = position_dodge(width = 0.9)) + 
  geom_hline(yintercept = 0.05) +
  labs(title = "Empirical Size") +
  xlab("K")+
  theme(legend.position = "none")



##K vs power
n = 600
#Coefficients
alpha = 0.5
beta = 0.5
gamma = 0.1
omega = 0.05
coef_X_alpha = c(1,1,1)
coef_X_beta = c(1,1,1)


##Checking size, qqplot
prop = c(0.8, 0.1, 0.1)
N_rep = 200       ##Number of monte carlo repetitions
mc = matrix(nrow = N_rep, ncol = 5)

for(j in 1: N_rep)
{
  other_covariates = cbind(rep(1,n), runif(n), runif(n))
  
  #print(alphaj*betaj)
  mc[j,] = sim_K(0.5, 0.5,  gamma = 0.1, other_covariates, coef_X_alpha, coef_X_beta, 
                 n, omega = 0.05)
}
power = apply((mc<=omega), 2, mean)
#se_size = sqrt(size*(1 - size)/ N_rep)
data = data.frame(K, power)
powerp <- data.table(data) %>%
  ggplot(aes(x = factor(K), y = power, group = 1))+
  geom_line(linetype = "dashed")+ geom_point()+
  labs(title = "Power") +
  xlab("K")+
  theme(legend.position = "none")

allp <- ggarrange(sizep, powerp)
ggsave("choice_K.png", allp, device = "png",
       width = 8, height = 4.5, units = "in", 
       dpi = 600)
