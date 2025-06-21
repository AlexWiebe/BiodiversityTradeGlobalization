
### Alex Wiebe
## Princeton University
## Description: Various helper functions to compute
## dynamics using systems of ODEs.

library(deSolve)
library(tidyverse)

faup.model = function (t, x, params){
  # t is vector of times
  # x is vector of initial values
  # params is vector of independent variables
  
  F1 = x[1]
  A11 = x[2]
  A12 = x[3]
  U1 = x[4]
  P1 = x[5]
  F2 = x[6]
  A22 = x[7]
  A21 = x[8]
  U2 = x[9]
  P2 = x[10]
  
  s1 = params["s1"]
  d1 = params["d1"]
  a1 = params["a1"]
  e = params["e"]
  h = params["h"]
  r = params["r"]
  m1 = params["m1"]
  
  s2 = params["s2"]
  d2 = params["d2"]
  a2 = params["a2"]
  m2 = params["m2"]
  
  dF1dt = s1*U1 - d1*F1*P1 - m1*d1*F1*P1*P2
  dA11dt = d1*F1*P1 - a1*A11 - e*(m2*d2*F2*P2*P1) + e*a2*A21
  dA12dt = m1*d1*F1*P1*P2 - a1*A12
  dU1dt = a1*A11 + a1*A12 + e*(m2*d2*F2*P2*P1) - e*a2*A21 - s1*U1
  dP1dt = r*P1*(A11 + A21 - h*P1)/(A11 + A21)
  
  dF2dt = s2*U2 - d2*F2*P2 - m2*d2*F2*P2*P1
  dA22dt = d2*F2*P2 - a2*A22 - e*(m1*d1*F1*P1*P2) + e*a1*A12
  dA21dt = m2*d2*F2*P2*P1 - a2*A21
  dU2dt = a2*A22 + a2*A21 + e*(m1*d1*F1*P1*P2) - e*a1*A12 - s2*U2
  dP2dt = r*P2*(A22 + A12 - h*P2)/(A22 + A12)
  
  dxdt = c(dF1dt, dA11dt, dA12dt, dU1dt, dP1dt, dF2dt, dA22dt, dA21dt, dU2dt, dP2dt)
  
  list(dxdt)
  
}


### Single country FAUP model
single.model = function (t, x, parms){
  # t is vector of times
  # x is vector of initial values
  # params is vector of independent variables
  
  F = x[1]
  A = x[2]
  U = x[3]
  P = x[4]
  
  s = params["s"]
  d = params["d"]
  a = params["a"]
  h = params["h"]
  r = params["r"]
  
  dFdt = s*U - d*F*P
  dAdt = d*F*P - a*A
  dUdt = a*A - s*U
  dPdt = r*P*(A - h*P)/(A)
  
  dxdt = c(dFdt, dAdt, dUdt, dPdt)
  
  list(dxdt)
  
}




















