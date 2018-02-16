setwd("~/Desktop/ENSAE_MS/S2/PAC_bayesian_online_clustering")
data = read.csv('data')
data[1] <- NULL 

mydata = data
N_iterations = 500
coeff = 2
R = max(sqrt(rowSums(mydata^2)))
K_max = 50

?PACBO
install.packages('PACBO')
library('PACBO')

# Exemple: generating 4 clusters of 100 points
set.seed(100)
Nb <- 4
d <- 5
T <- 10
proportion = rep(1/Nb, Nb)
Mean_vectors <- matrix(runif(d*Nb,min=-10, max=10),nrow=Nb,ncol=d, byrow=TRUE)
mydata <- matrix(replicate(T, rmnorm(1, mean= Mean_vectors[sample(1:Nb, 1, prob = proportion),],
varcov = diag(1,d))), nrow = T, byrow=T)
R <- max(sqrt(rowSums(mydata^2)))
##run the algorithm.
result <- PACBO(mydata, R,var_ind = TRUE, plot_ind = FALSE)
result

mydata[1:6,]

#######Modification of the package to see the evolution of the clustering with online clustering

instantaneous_loss = function(centers, instant_observation){
  if (class(centers) != 'matrix'){
    return (sum((centers-instant_observation)^2))
  }else{
    replicate_instant_observation = t(replicate(length(centers[,1]), instant_observation))
    distances = apply((centers - replicate_instant_observation)^2, 1, sum)
    return (min(distances))
  }
}

cumulative_loss = function(centers, observations){
  if (class(observations) != 'matrix'){
    return (instantaneous_loss(centers, observations))
  }else{
    if (class(centers) != 'matrix'){
      number_of_obs = length(observations[,1])
      centers_transform = t(replicate(number_of_obs, centers))
      distances = apply((centers_transform - observations)^2, 1, sum)
      return (sum(distances))
      
    }else{
      
      d = length(observations[1,])
      lth_centers=length(centers[,1])
      size_obs=length(observations[,1])
      C_mtx=replicate(size_obs,centers)
      Data_mtx=array(apply(observations,1,function(x) matrix(rep(x,lth_centers),nrow                           =lth_centers, byrow=T)),dim=c(lth_centers,d,size_obs))
      Diff=apply((C_mtx-Data_mtx)^2, c(1,3), sum)
      return(sum(apply(Diff,2,min)))
    }
  }
}

return(apply(Diff, 2, function(x) which(x == min(x))))


labels_function = function(centers, observations){
  if (class(observations) != 'matrix'){
    if (class(centers) != 'matrix'){
      return (c(1))
    }else{
      replicate_instant_observation = t(replicate(length(centers[,1]), observations))
      distances = apply((centers - replicate_instant_observation)^2, 1, sum)
      return (which(distances == min(distances)))
    }
  }else{
    if (class(centers) != 'matrix'){
      number_of_obs = length(observations[,1])
      centers_transform = t(replicate(number_of_obs, centers))
      distances = apply((centers_transform - observations)^2, 1, sum)
      return (rep(1, number_of_obs))
      
    }else{
      d = length(observations[1,])
      lth_centers=length(centers[,1])
      size_obs=length(observations[,1])
      C_mtx=replicate(size_obs,centers)
      Data_mtx=array(apply(observations,1,function(x) matrix(rep(x,lth_centers),nrow=lth_centers, byrow=T)),dim=c(lth_centers,d,size_obs))
      Diff=apply((C_mtx-Data_mtx)^2, c(1,3), sum)
      return(apply(Diff, 2, function(x) which(x == min(x))))
    }
  }
}

PACBO = function(mydata, R, coeff = 2, K_max = 50, scaling = FALSE, var_ind = FALSE, N_iterations = 500, plot_ind = FALSE, axis_ind = c(1,2)){
  
  if (class(mydata) != 'matrix'){
    mydata = matrix(mydata, nrow =1)
  }
  
  if (R < max(sqrt(rowSums(mydata^2)))){
    print(c('R should be bigger than the maximum Euclidean distance of observations'))
  }else{
    
    if (scaling){
      mydata = scale(mydata)
    }
    
    d = length(mydata[1,])
    T = length(mydata[,1])
    multiplier_R = 1.5
    
    Nclusters = rep(1, N_iterations)
    Niter_centers = list()
    parameter_means_proposal = list()
    sum_loss = rep(0, N_iterations)
    labels_inter = list()
    
    if (var_ind){
      nb_of_clusters = rep(1, T)
      pred_centers = list()
      predicted_loss = rep(0, T)
      proposal_loss = rep(0, T)
      labels_inter = list()
      
      lambda_1 = rep(1, T)
      lambda_1[2:T] = sapply(2:T, function(t) 2.5*coeff *(d+2)*(t-1)^(-0.5)/R)
      lambda_2 = sapply(1:T, function(t) (d+2)*(t)^(-0.5)/(multiplier_R * R)^4)
      
      
      
      c_1 = runiform_ball(5000, d, multiplier_R*R)
      index_1 = which(apply((c_1-t(replicate(5000, mydata[1,])))^2, 1, sum) == instantaneous_loss(c_1, mydata[1,]))
      pred_centers[[1]] = c_1[index_1,]
      predicted_loss[1] = instantaneous_loss(pred_centers[[1]], mydata[1,])
      labels_inter[[1]] = labels_function(pred_centers[[1]], mydata[1,])
      
      pb <- txtProgressBar(0, 1, char = '=')
      
      
      for (t in 2:T){
        
        tau_proposal = (K_max*t*d)^(-0.5)
        
        Nclusters[1] = nb_of_clusters[t-1]
        parameter_means_proposal[[1]] = kmeans(matrix(mydata[1:(t-1),], ncol=d), Nclusters[1], nstart=2, iter.max=10)$centers
        
        Niter_centers[[1]] = t(apply(parameter_means_proposal[[1]],1, function(x) rmt(1, mean=x, S = diag(tau_proposal, d), df =3)))
        
        
        while ( sum(sqrt(rowSums((Niter_centers[[1]]^2))) > multiplier_R * R) > 0 ){
          Niter_centers[[1]] = t(apply(parameter_means_proposal[[1]],1, function(x) rmt(1, mean=x, S = diag(tau_proposal, d), df =3)))
        }
        
        
        proposal_loss[1:(t-1)] = apply(matrix(mydata[1:(t-1),], nrow = t-1), 1, function (x) instantaneous_loss(Niter_centers[[1]], x))
        
        sum_loss[1] = sum(proposal_loss[1:(t-1)]) + 0.5 * sum(lambda_2 * (proposal_loss - predicted_loss)^2)
        
        for (n in 2:N_iterations){
          proposal_loss_temp = rep(0, T)
          transition_prob = c(transition_probability(Nclusters[n-1],Nclusters[n-1]-1,K_max), transition_probability(Nclusters[n-1],Nclusters[n-1],K_max), transition_probability(Nclusters[n-1],Nclusters[n-1]+1,K_max))
          transition_prob = transition_prob/sum(transition_prob)
          new_k = sample(c(Nclusters[n-1]-1,Nclusters[n-1],Nclusters[n-1]+1),1, prob = transition_prob)
          if (new_k == Nclusters[n-1]){
            m_t = parameter_means_proposal[[n-1]]
          }else{
            if(new_k >= t-1){
              new_k = t-1
              m_t = matrix(mydata[1:(t-1),], ncol=d, byrow=F)
            }else{
              m_t = kmeans(matrix(mydata[1:(t-1),],ncol=d), new_k, nstart=2, iter.max=10)$centers
            }
          }
          c_k_prime = t(apply(m_t,1, function(x) rmt(1, mean=x, S = diag(tau_proposal, d), df =3)))
          
          if (sum( sqrt(rowSums(c_k_prime^2)) < multiplier_R * R ) == new_k) {
            
            log_numerator_prop = log(apply(matrix(1:Nclusters[n-1], ncol=1),1,function(x) dmt(Niter_centers[[n-1]][x,], mean = parameter_means_proposal[[n-1]][x,], S = diag(tau_proposal,d), df = 3)))
            
            log_denominator_prop = log(apply(matrix(1:new_k, ncol=1),1,function(x) dmt(c_k_prime[x,], mean=m_t[x,], S = diag(tau_proposal,d), df = 3)))
            
            log_numerator_prior = rep((log(gamma(d/2+1)) - (d/2)*log(pi) - d*log(multiplier_R*R)), new_k)
            log_denominator_prior = rep((log(gamma(d/2+1)) - (d/2)*log(pi) - d*log(multiplier_R*R)), Nclusters[n-1])
            
            
            ln_division = sum(c(log_numerator_prior, log_numerator_prop)-c(log_denominator_prior, log_denominator_prop))
            
            
            proposal_loss_temp[1:(t-1)] = apply(matrix(mydata[1:(t-1),], nrow = t-1), 1, function (x) instantaneous_loss(c_k_prime, x))
            
            s_loss_prime = sum(proposal_loss_temp) + 0.5 * sum(lambda_2 * (proposal_loss_temp - predicted_loss)^2)
            
            ln_accept_ratio=(-lambda_1[t-1]*(s_loss_prime-sum_loss[n-1]))+ln_division+log(transition_probability(new_k, Nclusters[n-1],K_max))-log(transition_probability(Nclusters[n-1], new_k, K_max))
          }else{
            ln_accept_ratio = log(0)
          }
          
          
          bool= (log(runif(1)) < ln_accept_ratio)
          
          if (bool){
            Niter_centers[[n]] = c_k_prime
            parameter_means_proposal[[n]] = m_t
            Nclusters[n] = new_k
            sum_loss[n] = s_loss_prime
          }else{
            Niter_centers[[n]] = Niter_centers[[n-1]]
            parameter_means_proposal[[n]] = parameter_means_proposal[[n-1]]
            Nclusters[n] = Nclusters[n-1]
            sum_loss[n] = sum_loss[n-1]
          }
          
        }
        
        nb_of_clusters[t] = Nclusters[N_iterations]
        pred_centers[[t]]= Niter_centers[[N_iterations]]
        predicted_loss[t] = instantaneous_loss(pred_centers[[t]], mydata[t,])
        labels_inter[[t]] = labels_function(pred_centers[[t]], mydata[1:t,])
        setTxtProgressBar(pb, t/T)
        
      }
      
      labels = labels_function(pred_centers[[T]], mydata)
      if (plot_ind){
        
        plot(mydata[, axis_ind[1]], mydata[, axis_ind[2]], xlim =c(min(mydata[, axis_ind[1]])-5, max(mydata[, axis_ind[1]])+5), ylim = c(min(mydata[, axis_ind[2]])-5, max(mydata[, axis_ind[2]])+5) , col = 1+labels, xlab = paste(c('axis_'), axis_ind[1], sep = ''), ylab = paste(c('axis_'), axis_ind[2], sep = ''))
        points(pred_centers[[T]][, axis_ind[1]], pred_centers[[T]][, axis_ind[2]], pch = 17, col = 1)
        lgd = c('pred_centers       ', sapply(1:length(unique(labels)), function(x) paste(c('cluster'), as.character(x), sep = ' '), simplify = 'vector'))
        legend('topright', legend = lgd, col = seq(length(unique(labels))+1), pch = c(17, rep(1, length(unique(labels)))), cex = 0.7, xjust = 1)
        
      }
      return (list('centers' = pred_centers, 'nb_of_clusters' = nb_of_clusters, 'labels' = labels_inter))
      #return (list('centers' = pred_centers[[T]], 'nb_of_clusters' = nb_of_clusters[T], 'labels' = labels))
    }else{
      
      lambda_1 = coeff *(d+2)*(T-1)^(-0.5)/R
      tau_proposal = (K_max*T*d)^(-0.5)
      Nclusters[1] = 1
      parameter_means_proposal[[1]] = kmeans(matrix(mydata, ncol=d), Nclusters[1], nstart=2, iter.max=10)$centers
      
      Niter_centers[[1]] = matrix(rmt(1, mean=as.vector(parameter_means_proposal[[1]]),  S = diag(tau_proposal, d), df =3), nrow = 1)
      
      while ( sqrt(sum((Niter_centers[[1]]^2))) > multiplier_R * R ){
        Niter_centers[[1]] = matrix(rmt(1, mean=as.vector(parameter_means_proposal[[1]]), S = diag(tau_proposal, d), df =3),nrow = 1)
      }
      
      
      sum_loss[1] = cumulative_loss(Niter_centers[[1]], mydata)
      
      
      for (n in 2:N_iterations){
        transition_prob = c(transition_probability(Nclusters[n-1],Nclusters[n-1]-1,K_max), transition_probability(Nclusters[n-1],Nclusters[n-1],K_max), transition_probability(Nclusters[n-1],Nclusters[n-1]+1,K_max))
        transition_prob = transition_prob/sum(transition_prob)
        
        new_k = sample(c(Nclusters[n-1]-1,Nclusters[n-1],Nclusters[n-1]+1),1, prob = transition_prob)
        if (new_k == Nclusters[n-1]){
          m_t = parameter_means_proposal[[n-1]]
        }else{
          if(new_k >= T){
            new_k = T
            m_t = matrix(mydata, ncol=d, byrow=F)
          }else{
            m_t = kmeans(matrix(mydata,ncol=d), new_k, nstart=2, iter.max=50)$centers
          }
        }
        
        
        c_k_prime = t(apply(m_t,1, function(x) rmt(1, mean=x, S = diag(tau_proposal, d), df =3)))
        
        if (sum( sqrt(rowSums(c_k_prime^2)) < multiplier_R * R ) == new_k) {
          
          log_numerator_prop = log(apply(matrix(1:Nclusters[n-1],ncol=1),1,function(x) dmt(Niter_centers[[n-1]][x,], mean=parameter_means_proposal[[n-1]][x,], S=diag(tau_proposal,d), df=3)))
          
          log_denominator_prop = log(apply(matrix(1:new_k,ncol=1),1,function(x) dmt(c_k_prime[x,], mean=m_t[x,], S=diag(tau_proposal,d), df=3)))
          
          log_numerator_prior = rep((log(gamma(d/2+1)) - (d/2)*log(pi) - d*log(multiplier_R*R)), new_k)
          log_denominator_prior = rep((log(gamma(d/2+1)) - (d/2)*log(pi) - d*log(multiplier_R*R)), Nclusters[n-1])
          
          
          ln_division = sum(c(log_numerator_prior, log_numerator_prop)-c(log_denominator_prior, log_denominator_prop))
          
          s_loss_prime = cumulative_loss(c_k_prime, mydata)
          
          
          ln_accept_ratio=(-lambda_1*(s_loss_prime-sum_loss[n-1]))+ln_division+log(transition_probability(new_k, Nclusters[n-1],K_max))-log(transition_probability(Nclusters[n-1], new_k, K_max))
        }else{
          ln_accept_ratio = log(0)
        }
        bool=(log(runif(1)) < ln_accept_ratio)
        if (bool){
          Niter_centers[[n]] = c_k_prime
          parameter_means_proposal[[n]] = m_t
          Nclusters[n] = new_k
          sum_loss[n] = s_loss_prime
        }else{
          Niter_centers[[n]] = Niter_centers[[n-1]]
          parameter_means_proposal[[n]] = parameter_means_proposal[[n-1]]
          Nclusters[n] = Nclusters[n-1]
          sum_loss[n] = sum_loss[n-1]
        }
        
        
      }
      pred_centers = Niter_centers[[N_iterations]]
      nb_clusters = Nclusters[N_iterations]
      labels = labels_function(pred_centers, mydata)
      
      if (plot_ind){
        
        plot(mydata[, axis_ind[1]], mydata[, axis_ind[2]], xlim =c(min(mydata[, axis_ind[1]])-5, max(mydata[, axis_ind[1]])+5), ylim = c(min(mydata[, axis_ind[2]])-5, max(mydata[, axis_ind[2]])+5), col = 1+labels, xlab = paste(c('axis_'), axis_ind[1], sep = ''), ylab = paste(c('axis_'), axis_ind[2], sep = ''))
        points(pred_centers[, axis_ind[1]], pred_centers[, axis_ind[2]], col= 1, pch=17)
        lgd = c('pred_centers       ', sapply(1:length(unique(labels)), function(x) paste(c('cluster'), as.character(x), sep = ' '), simplify = 'vector'))
        legend('topright', legend = lgd, col = seq(length(unique(labels))+1), pch = c(17, rep(1, length(unique(labels)))), cex = 0.7, xjust = 1)
        
      }
      
      return (list('predicted_centers' = pred_centers, 'nb_of_clusters' = nb_clusters, 'labels' = labels))
    }
  }
}

runiform_ball = function(n, d, R){
  simulations = matrix(NA, nrow = n, ncol = d)
  for (i in 1:n){
    simul = runif(d, min = -R, max = R)
    while(sqrt(sum(simul^2)) > R){
      simul = runif(d, min = -R, max = R)
    }
    simulations[i,] = simul
  }
  return (simulations)
}


transition_probability = function(x, z, K_max, prob = 1/3){
  if((x > K_max)|(z > K_max)|(x < 1)|(z < 1)){
    q = 0
  }else{ if(x == 1){
    q = 0.5*((z == 1)|(z == 2))
  }else{ if(x == K_max){
    q = 0.5*((z == K_max)|(z == (K_max-1)))
  }else{if(x == z){
    q = prob
  }else{q = (1-prob)/2 * (( x < z & z < (x+2))|(z < x & (z > (x-2))))}
    
  }
  }
  }
  return(q)
}




