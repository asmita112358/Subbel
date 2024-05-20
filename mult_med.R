##Multivariate mediators

rm(list = ls())
library(tidyverse)
library(ggplot2)
sim.multmed = function(alpha, beta, gamma, other_covariates, coef_X_alpha, coef_X_beta, 
                       n, K, num.med, omega = 0.05)
{
  ##genrate data
  Y <- S <- c() ##outcome, exposure
  G = matrix(nrow = n, ncol = num.med)
  S = runif(n, 1,2)
  for(j in 1:num.med)
  {
    G[,j] = S*alpha[j] + rowSums(coef_X_alpha*other_covariates) + rnorm(n)
    
  }
  Y = G%*%beta + S*gamma + rowSums(coef_X_beta*other_covariates) + rnorm(n)

 #print(mean(Y))
  
  niter = 100
  tstat <- pval <- matrix(nrow = num.med, ncol = niter)
  CCT <- p_cct <- c()
  for(iter in 1:niter)
  {
    
    #sample_partitioning
    all <- sample.int(n = n, size = n, replace = F)
    
    index <- lapply(1:K, function(ii){
      m <- floor(n/K)
      all[((ii-1)*m+1):(ii*m)]
    })
    S_K = matrix(nrow = num.med, ncol = K)
    for(i in 1:K)
    {
      S_s = S[index[[i]]]
      G_s = G[index[[i]],]
      Y_s = Y[index[[i]],]
      if(is.null(dim(G_s))){G_s = matrix(G_s, ncol = num.med)}
      other_covariates_s <- other_covariates[index[[i]],]
      t1 <- t2<- c()
      
      for(med in 1:num.med)
      {
        obj1 = lm(G_s[,med] ~ -1 + cbind(S_s, other_covariates_s)) 
        t1[med] = summary(obj1)[["coefficients"]][1,3]
        
      }
      
      obj2 = lm(Y_s ~ -1 + cbind(G_s, S_s, other_covariates_s))
      
      t2 = summary(obj2)[["coefficients"]][(1:num.med),3]
      
      #if(i == K)print(t1*t2)
      S_K[,i] = (t1*t2)/sqrt(t1^2 + t2^2)
      
      #print(S_K[,i])
    }
    #test_stat
    #print(rowMeans(S_K))
   
    tstat[,iter] = apply(S_K,1, function(dummy){sqrt(length(dummy))*mean(dummy)/sqrt(var(dummy))})
    
    pval[,iter] = 2*pt(-abs(tstat[,iter]), df = K-1)
    rm(S_K)
    rm(obj1)
    rm(obj2)
    #print(pval[,iter])
  }
  for(med in 1:num.med)
  {
    wts = runif(niter)
    wts = wts/sum(wts)
    CCT[med] = sum(tan(3.141593*(0.5 - pval[med,]))*wts) 
    p_cct[med] = 1- pcauchy(CCT[med])
  }
  #print(p_cct)
  ##Fitting full model
  t1 <- t2 <- c()
  p1 <- p2 <- c()
  for(med in 1:num.med)
  {
    obj1 = lm(G[,med] ~ -1 + cbind(S, other_covariates))
    t1[med] = summary(obj1)[["coefficients"]][1,3]
    p1[med] = summary(obj1)[["coefficients"]][1,4]
              
  }
  obj2 = lm(Y ~ -1 + cbind(G, S, other_covariates))
  t2 = summary(obj2)[["coefficients"]][(1:num.med),3]
  p2 = summary(obj2)[["coefficients"]][(1:num.med),4]
  ##Sobel's test
  T.sobel = t2/sqrt(1+(t2/t1)^2)
  p_sobel = 2*pnorm(-abs(T.sobel),mean = 0, sd=1, lower.tail = TRUE)
  
  
  #Pmax test
 
  return(cbind(p_cct, p_max = pmax(p1, p2), p_sobel))
}

##parameters

n = 600
#num.med = 1
K = 12

#Coefficients

alpha = c(1,1,0,0,1)
beta = c(0,0,1,1,0)
gamma = 0.1
omega = 0.05
coef_X_alpha = c(1,1,1)
coef_X_beta = c(1,1,1)
N_rep = 100       ##Number of monte carlo repetitions
mc = matrix(nrow = length(alpha), ncol = 3)
t = matrix(0, nrow = length(alpha), ncol = 3)
#samp = sample(1:3, N_rep, replace = T, prob = prop)
for(j in 1: N_rep)
{
  other_covariates = cbind(rep(1,n), runif(n), runif(n))
 
  #print(alphaj*betaj)
  mc = sim.multmed(alpha, beta,  gamma = 0.1, other_covariates, coef_X_alpha, coef_X_beta, 
                 n, K = 12, num.med = length(alpha), omega = 0.05)
  print(mc)
  
  mc = 1* (mc <= 0.05)
  
  mc = t + mc
  t = mc
  
}
size = mc/N_rep
colnames(size) = c("Subbel", "MaxP", "Sobel")
se_size = sqrt(size*(1 - size)/ N_rep)
lci = size - 1.96*se_size
uci = size + 1.96*se_size
mediator = rep(c("G1","G2", "G3", "G4", "G5"),3)
size = gather(data.frame(size), key = "method")
se_size = gather(data.frame(se_size))[,2]
size = data.frame(mediator, size, se_size)
colnames(size) <- c("mediator","method","size", "se_size")
size$method <- factor(size$method, levels = c("Subbel", "MaxP", "Sobel"))
#size Graph
#F8766D
#7CAE00
#00BFC4
#C77CFF

barp <- ggplot(size, aes(x = method,y = size, fill = method))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#F8766D","#00BFC4" , "#C77CFF" ))+
  geom_errorbar(aes(ymin = size - 1.96*se_size, ymax = size + 1.96*se_size), width = 0.2,position = position_dodge(width = 0.9)) +
  facet_grid(.~mediator)+
  geom_hline(yintercept = 0.05)+
  labs(title = "Empirical Size for multiple mediators")
ggsave("size_mult.png", barp, device = "png",
       width = 8, height = 4.5, units = "in", 
       dpi = 600)
  
  ##power
  intensity = seq(0, 0.5, length.out = 10)
alpha = intensity%o%rep(1,5)
beta = intensity%o%rep(1,5)
gamma = 0.1
omega = 0.05
coef_X_alpha = c(1,1,1)
coef_X_beta = c(1,1,1)
N_rep = 100
power = array(dim = c(length(intensity), 5, 3))
for(count in 1:length(intensity))
{
  alphaj = alpha[count,]
  betaj = beta[count,]
  mc = matrix(nrow = 5, ncol = 3)
  t = matrix(0, nrow = 5, ncol = 3)
  for(j in 1:N_rep)
  {
    other_covariates = cbind(rep(1,n), runif(n), runif(n))
    
    mc = sim.multmed(alphaj, betaj,  gamma = 0.1, other_covariates, coef_X_alpha, coef_X_beta, 
                         n, K = 12,num.med = length(alphaj), omega = 0.05)
    mc = 1* (mc <= 0.05)
    
    mc = t + mc
    t = mc
  }
  power[count,,] = mc/N_rep
  print(c(alphaj, betaj))
}

write.csv(power, "powermultmed.csv")

power_tab = read.csv("powermultmed.csv")
power_tab = power_tab[,-1]
colnames(power_tab) = rep(c("Subbel", "MaxP", "Sobel"), each = 5)
power_tab = gather(data.frame(power_tab), key = "method")
mediator = rep(c("G1", "G2", "G3", "G4", "G5"), each = 10, times = 3)
method = rep(c("Subbel", "MaxP", "Sobel"), each = 50)
intensity = rep(seq(0, 0.5, length.out = 10), 15)
pow_tab = data.frame(method, mediator,intensity, power_tab$value)
colnames(pow_tab) = c("method", "mediator", "intensity", "power")
pow_tab$method <- factor(pow_tab$method, levels = c("Subbel", "MaxP", "Sobel"))
write.csv(pow_tab, "power_tab_multmed.csv")
pow_plot = ggplot(pow_tab, aes(x = intensity, y = power, group = method, color = method, shape = method))+
  scale_color_manual(values = c("#F8766D","#00BFC4" , "#C77CFF" ))+
  geom_line() + facet_grid(.~mediator) + labs(x = "Signal strength", y = "Empirical power")+
  theme(legend.position = "bottom")

ggsave("power_mult.png", pow_plot, device = "png",
       width = 10, height = 3, units = "in", 
       dpi = 600)
