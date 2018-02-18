library('PACBO')
setwd("~/Desktop/ENSAE_MS/S2/PAC_bayesian_online_clustering")

coeff = 2
K_max = 50
N_iterations = 500
mydata <- read.csv('test.csv')
mydata <- as.matrix(sapply(mydata, as.numeric)) 
mydata <- mydata[,2:length(mydata[1,])]
R <- max(sqrt(rowSums(mydata^2)))


PACBO(mydata, R, coeff = 2, K_max = 50,var_ind = TRUE, N_iterations = 500, plot_ind=TRUE)



#####
library('PACBO')
setwd("~/Desktop/ENSAE_MS/S2/PAC_bayesian_online_clustering")

coeff = 2
K_max = 50
N_iterations = 500
mydata <- read.csv('test.csv')
mydata <- as.matrix(sapply(mydata, as.numeric)) 
mydata <- mydata[,2:length(mydata[1,])]
R <- max(sqrt(rowSums(mydata^2)))


d = length(mydata[1,])
T = length(mydata[,1])
multiplier_R = 1.5
Nclusters = rep(1, N_iterations)
Niter_centers = list()
parameter_means_proposal = list()
sum_loss = rep(0, N_iterations)

# IF

nb_of_clusters = rep(1, T)
pred_centers = list()
predicted_loss = rep(0, T)
proposal_loss = rep(0, T)

lambda_1 = rep(1, T)
lambda_1[2:T] = sapply(2:T, function(t) 2.5*coeff *(d+2)*(t-1)^(-0.5)/R)
lambda_2 = sapply(1:T, function(t) (d+2)*(t)^(-0.5)/(multiplier_R * R)^4)



c_1 = runiform_ball(5000, d, multiplier_R*R)
index_1 = which(apply((c_1-t(replicate(5000, mydata[1,])))^2, 1, sum) == instantaneous_loss(c_1, mydata[1,]))
pred_centers[[1]] = c_1[index_1,]
predicted_loss[1] = instantaneous_loss(pred_centers[[1]], mydata[1,])




t=4




tau_proposal = (K_max*t*d)^(-0.5)

Nclusters[1] = nb_of_clusters[t-1]
parameter_means_proposal[[1]] = kmeans(matrix(mydata[1:(t-1),], ncol=d), Nclusters[1], nstart=2, iter.max=10)$centers

Niter_centers[[1]] = t(apply(parameter_means_proposal[[1]],1, function(x) rmt(1, mean=x, S = diag(tau_proposal, d), df =3)))


while ( sum(sqrt(rowSums((Niter_centers[[1]]^2))) > multiplier_R * R) > 0 ){
  Niter_centers[[1]] = t(apply(parameter_means_proposal[[1]],1, function(x) rmt(1, mean=x, S = diag(tau_proposal, d), df =3)))
}

proposal_loss[1:(t-1)] = apply(matrix(mydata[1:(t-1),], nrow = t-1), 1, function (x) instantaneous_loss(Niter_centers[[1]], x))

sum_loss[1] = sum(proposal_loss[1:(t-1)]) + 0.5 * sum(lambda_2 * (proposal_loss - predicted_loss)^2)

sqrt(rowSums(Niter_centers[[1]]^2))


# Check


tau_proposal
Nclusters
parameter_means_proposal
Niter_centers
proposal_loss
sum_loss



#n = 4
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
### Check
proposal_loss_temp
transition_prob
new_k
m_t
c_k_prime
ln_division
s_loss_prime
ln_accept_ratio
bool
Niter_centers
parameter_means_proposal
Nclusters
sum_loss



# 

nb_of_clusters[t] = Nclusters[N_iterations]
pred_centers[[t]]= Niter_centers[[N_iterations]]
predicted_loss[t] = instantaneous_loss(pred_centers[[t]], mydata[t,])


##Check
nb_of_clusters
pred_centers
predicted_loss


#

labels = labels_function(pred_centers[[T]], mydata)

