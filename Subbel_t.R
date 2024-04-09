##sobel -> sample partiioning -> tstatistic p value -> heavy tailed combination test
rm(list = ls())
source("~/Documents/OneDrive - Texas A&M University/Project3/sobel/HT_combtest.R")

library(data.table)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)

sim = function(alpha, beta, gamma, other_covariates, coef_X_alpha, coef_X_beta, 
               n, K, omega = 0.05)
{
  
  ##Generating Data
  Y<- G <- S <- c()
  
  S = runif(n, 1,2)
  G = S*alpha + rowSums(coef_X_alpha*other_covariates) + rnorm(n)
  Y = G*beta + S*gamma + rowSums(coef_X_beta*other_covariates) + rnorm(n)
  
  
  
  ##Subsampling
  
  #K = 10
  niter = 100
  
  tstat = c()
  pval = c()
  for(iter in 1:niter)
  {
    #Sample partitioning
    all <- sample.int(n = n, size = n, replace = F)
    
    index <- lapply(1:K, function(ii){
      m <- floor(n/K)
      all[((ii-1)*m+1):(ii*m)]
    })
    
    #identical(sort(do.call('c', index)), 1:n) #check the sample partition
    
    S_K = c()
    for(i in 1:K)
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
    
    tstat[iter] = sqrt(K)*mean(S_K)/sqrt(var(S_K))
    pval[iter] = 2*pt(-abs(tstat[iter]), df = K-1)
    ##print(tstat)
  }
  
  #p_global = combination.test(pval, method ="Levy" ,truncate.threshold = 0.97)
  wts = runif(niter)
  wts = wts/sum(wts)
  #hist(wts)
  CCT = sum(tan(3.141593*(0.5 - pval))*wts) ##Cauchy combination statistic
  p_cct = 1- pcauchy(CCT)
  #ABtest
  
  mat = DBmypackage6two::One_MC_DB_X_noProj_C(cbind(S,G,Y), other_covariates,
                                              n, omega, num_med = 1, lambda1all = 2, lambda2all = 2, N_boot_out = 500, N_boot_in = 0)
  pvalue = mat[1,1]
  
  #fitting full model for remaining methods
  object1 = lm(G ~ -1 + cbind(S, other_covariates))
  object2 = lm(Y ~ -1 + cbind(G, S, other_covariates))
  t1 = summary(object1)[["coefficients"]][1,3]
  t2 = summary(object2)[["coefficients"]][1,3]
  p1 = summary(object1)[["coefficients"]][1,4]
  p2 = summary(object2)[["coefficients"]][1,4]
  ##fix this part. wtf is happening
  
  
  ##sobel
  T.sobel = t2/sqrt(1+(t2/t1)^2)
  p.sobel = 2*pnorm(-abs(T.sobel),mean = 0, sd=1, lower.tail = TRUE)
  
  return(c(p_cct, pvalue, max(p1,p2), p.sobel))
  
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
prop = c(0.4, 0.3, 0.3)
N_rep = 200       ##Number of monte carlo repetitions
mc = matrix(nrow = N_rep, ncol = 4)
samp = sample(1:3, N_rep, replace = T, prob = prop)
for(j in 1: N_rep)
{
  other_covariates = cbind(rep(1,n), runif(n), runif(n))
  alphaj = 0 + alpha*(samp[j]== 2)
  betaj = 0 + beta*(samp[j] == 3)
  #print(alphaj*betaj)
  mc[j,] = sim(alphaj, betaj,  gamma = 0.1, other_covariates, coef_X_alpha, coef_X_beta, 
               n, K = 12, omega = 0.05)
}
size = apply((mc<=omega), 2, mean)
se_size = sqrt(size*(1 - size)/ N_rep)
##qqplot

theo_q = seq(0.01,1, 0.01)
samp_q = matrix(nrow = 100, ncol = 4)
for(i in 1:4)
{
  samp_q[,i] = quantile(mc[,i], probs = theo_q)
}
data = data.frame(theo_q, samp_q)
colnames(data) <- c("theo_q", "Subbel", "ABtest", "MaxP", "Sobel")


data.plot <- data.table(data) 


qqp <- data.plot %>%
  melt(id.vars = 1,
       variable.name = "Method") %>%
  ggplot(aes(x = theo_q, y = value,
             group = Method, color = Method,
             shape = Method)) +
  geom_point() +
  geom_abline() +
  labs(title = "QQplot for p-values",
       x = "Theoretical quantile",
       y = "Sample quantile") +
  theme(legend.position = "bottom")

barp <- data.table(Method = factor(colnames(data.plot)[2:5],
                                   levels = c("Subbel", "ABtest",
                                              "MaxP", "Sobel")),
                   Size = size) %>%
  ggplot(aes(x = Method, y = Size, fill = Method)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = size - 1.96*se_size, ymax = size + 1.96*se_size), width = 0.2,position = position_dodge(width = 0.9)) + 
  geom_hline(yintercept = 0.05) +
  labs(title = "Empirical Size") +
  theme(legend.position = "none")

allp <- ggarrange(qqp, barp, common.legend = T)

ggsave("size_cct_case_4.png", allp, device = "png",
       width = 8, height = 4.5, units = "in", 
       dpi = 600)

rm(mc)
allp


##Power check
n = 200
#n = 500
intensity = seq(0, 0.5, length.out = 10)

power = matrix(nrow = length(intensity), ncol = 4)
N_rep = 100
for(count in 1:length(intensity))
{
  alphaj = intensity[count]
  betaj = intensity[count]
  mc = matrix(nrow = N_rep, ncol = 4)
  for(j in 1:N_rep)
  {
    other_covariates = cbind(rep(1,n), runif(n), runif(n))
    
    mc[j,] = sim(alphaj, betaj,  gamma = 0.1, other_covariates, coef_X_alpha, coef_X_beta, 
                 n, K = 5, omega = 0.05)
    
  }
  power[count,] = apply((mc <= omega), 2, mean)
  print(c(alphaj, betaj))
}

colnames(power) = c("Subbel", "ABtest", "MaxP", "Sobel")
data = data.frame(intensity, power)
write.csv(power, "power200.csv")
data.plot = data.table(data)

power <- data.plot %>% 
  melt(id.vars = 1,
       variable.name = "Method") %>%
  ggplot(aes(x = intensity, y = value,
             group = Method, color = Method,
             shape = Method))  + geom_line() + labs(title = "Empirical Power", x = "Signal strength", y = "power")
ggsave("power_cct_200.png", power, device = "png",
       width = 4, height = 3, units = "in", 
       dpi = 600)
