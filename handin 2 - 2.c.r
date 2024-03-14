r0 <- 0.02
thetaQ <- 0.05
kappa <- 0.1
sigma <- 0.015
timetomat <- 1
N <- 5
alphas <- rep(1, N)  # Since alpha_i = 1 for all i
Ti <- seq(1, N) + 1  # Ti = i + 1
K <- 4.5


ZCBvasicek<-function(r,tau,Pparam,riskpremium=0)
{ thetaQ<-Pparam[1]+riskpremium
kap<-Pparam[2]
sig<-Pparam[3]

Btau<-(1-exp(-kap*tau))/kap
Atau<-((Btau-tau)*(kap^2*thetaQ-0.5*sig^2)/kap^2 - sig^2*Btau^2/(4*kap))

ZCBvasicek<-exp(-r*Btau+Atau) 

}

PUTvasicek<-function(r,tauPUT,tauZCB,strike,Pparam,riskpremium=0){ 
  kap<-Pparam[2]
  sig<-Pparam[3]
  
  PtauPUT<-ZCBvasicek(r,tauPUT,Pparam,riskpremium)
  PtauZCB<-ZCBvasicek(r,tauZCB,Pparam,riskpremium)
  
  sigZCB<-sig*(1-exp(-kap*(tauZCB-tauPUT)))*sqrt((1-exp(-2*kap*(tauPUT)))/(2*kap))/kap
  
  h<-log(PtauZCB/(PtauPUT*strike))/sigZCB+0.5*sigZCB
  
  PUTvasicek<-strike*PtauPUT*pnorm(-h+sigZCB)-PtauZCB*pnorm(-h)
  
}

CALLvasicek<-function(r,tauCALL,tauZCB,strike,Pparam,riskpremium=0){
  PtauCALL<-ZCBvasicek(r,tauCALL,Pparam,riskpremium)
  PtauZCB<-ZCBvasicek(r,tauZCB,Pparam,riskpremium)
  CALLvasicek<-PUTvasicek(r,tauCALL,tauZCB,strike,Pparam,riskpremium)+PtauZCB-strike*PtauCALL
}



# Vasicek model functions using our notation
B <- function(t, T) {
  return((1 - exp(-kappa * (T - t))) / kappa)
}

A <- function(t, T) {
  return((thetaQ - 0.5 * (sigma / kappa) ^ 2) * (B(t, T) - (T - t)) - (sigma ^ 2) / (4 * kappa) * B(t, T) ^ 2)
}

# Calculate B(t, Ti), A(t, Ti), and P(t, Ti) for each Ti
B_values <- B(timetomat, Ti)
A_values <- A(timetomat, Ti)
P_values <- exp(A_values - B_values * r0)

# Define the price of the bond at T = 1
bond_price <- function(T, Ti, r_star) {
  B_values <- B(T, Ti)
  A_values <- A(T, Ti)
  return(sum(exp(A_values - B_values * r_star)))
}


# Solve for r_star
r_star_equation <- function(r_star, K) {
  return(bond_price(timetomat, Ti, r_star) - K)
}


# Solving for r_star
r_star <- uniroot(r_star_equation, interval = c(-0.1, 0.2), K = K)$root


AdjustedStrikes <- rep(NA, N)
for (i in 1:N){
  T[i] <- i+1
  AdjustedStrikes[i] <- bond_price(timetomat, T[i], r_star)
}

# Calculate adjusted strikes Ki using the solved r_star
Ki <- exp(A(timetomat, Ti) - B(timetomat, Ti) * r_star)


Ti <- seq(1, N, 1)
call <- rep(NA, length(Ti))
for (i in 1:length(Ti)){
  call[i] <- CALLvasicek(r0, tauCALL = timetomat, tauZCB = timetomat + Ti[i], 
                         AdjustedStrikes[i], c(thetaQ, kappa, sigma), riskpremium=0)
}

CallPrice <- sum(alphas * call)


#Results
Ki
sum(Ki)
r_star
CallPrice
call


