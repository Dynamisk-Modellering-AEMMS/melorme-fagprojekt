# SIR model in R
#
# Demonstration of the model
# Example: Influenza in a boarding school 
# (British Med Journal 4th March 1978)
#
rm(list=ls())
# get the ODE solver
library(deSolve)

# make the function of the system
sir_model = function (current_timepoint, state_values, parameters)
{
  # create state variables (local variables)
  sb = state_values [1]        # susceptibles
  eb = state_values [2]        # infectious
  ab = state_values [3]        # recovered
  ib = state_values [4]
  ibm = state_values [5]
  rb = state_values [6]
  sw = state_values [7]
  ew = state_values [8]
  aw = state_values [9]
  iw = state_values [10]
  iwm = state_values [11]
  rw = state_values [12]
  sv = state_values [13]
  ev = state_values [14]
  iv = state_values [15]
  
  
  with ( 
    as.list (parameters),     # variable names within parameters can be used 
    {
      # compute derivatives
      pa   = (t*b*tarm*zigma)
      fa   = t*b*tarm*zigma
      ka   = t*b*tarm*zigma
      fptz = p*t*zigma
      g    = t*P/zigma
      h    = t*f/zigma*k/zigma*v/zigma
      i    = t*k/zigma*v/zigma
      j    = t*H*v/zigma
      

      
      zigma = p+f+k+v
      dPx = pi-pa
      dFx = fi-fa
      dKx = ki-ka
      
      pp = pi/(pi+fi+ki)
      fp = fi/(pi+fi+ki)
      kp = ki/(pi+fi+ki)
      
      #1 = pp+fp+kp
      dtarm = fptz -(pi+fi+ki)
      #dtarm = 0
      dP = pa*g
      dUx = g
      dF = fa-h
      dk = ka+h-i
      dv = i-j
      dvx = j
      
      
      
      # combine results
      results = c (dsb,deb,dab,dib,dibm,rb,sw,ew,aw,iw,iwm,rw,sv,ev,iv)
      list (results)
    }
  )
}
# parameters
alpha = 1/(16*365)
Rhob=0.5
Rhow=0.5
p=0.5
r=0.5
qa=0.5
qi=0.5
qr=0.5
betaw=0.33
rhow=1/7.5
gammaw=1/8.5
myw=1/(70*365)
pib=1/(15*365)
betab=0.33
rhob=1/7.5
gammab=1/8.5
myb=1/(18.60*365)
piv=500
betav=0.33
bv=0.5
rhov=1/3.5
myv=1/21
ny=0.5
parameter.list <- c(alpha,Rhob,Rhow,p,r,qa,qi,qr,betaw,rhow,gammaw,myw,pib,betab,rhob,gammab,myb,piv,betav,bv,rhov,myv,ny)


# initial values
sw0 = 5000
ew0 =  200
aw0 =  200
iw0 =  50
iwm0 =  20
rw0 = 10
sb0 = 10000
eb0 =  220
ab0 =  400
ib0 =  100
ibm0 = 120
rb0 = 200 
sv0 =  500
ev0 =  100
iv0 =  100

initial.values <- c(sw=sw0,ew=ew0,aw=aw0,iw=iw0,iwm=iwm0,rw=rw0,sb=sb0,eb=eb0,ab=ab0,ib=ib0,ibm=ibm0,rb=rb0,sv=sv0,ev=ev0,iv=iv0)

# Output timepoints
time.points <- seq(0,15,by=1)

# simulate the epidemic
output <- ode(y=initial.values,times = time.points,func= sir_model, parms = parameter.list)
plot(output)
# Plot the result 
plot(sb~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,8000),ylab='# of boys',xlab='Time (days)')
# Add S
lines(eb~time,data=output,type='l',lwd=3,lty=1,col='blue')
# And R
lines(ab~time,data=output,type='l',lwd=3,lty=3,col='red')
lines(ib~time,data=output,type='l',lwd=3,lty=3,col='cyan')
lines(ibm~time,data=output,type='l',lwd=3,lty=3,col='black')
lines(rb~time,data=output,type='l',lwd=3,lty=3,col='brown')
lines(sw~time,data=output,type='l',lwd=3,lty=3,col='purple')
lines(ew~time,data=output,type='l',lwd=3,lty=3,col='grey')
lines(aw~time,data=output,type='l',lwd=3,lty=3,col='hotpink')
lines(iw~time,data=output,type='l',lwd=3,lty=3,col='green')
lines(iwm~time,data=output,type='l',lwd=3,lty=3,col='green')
lines(rw~time,data=output,type='l',lwd=3,lty=3,col='green')
lines(sv~time,data=output,type='l',lwd=3,lty=3,col='green')
lines(ev~time,data=output,type='l',lwd=3,lty=3,col='green')
lines(iv~time,data=output,type='l',lwd=3,lty=3,col='green')

# # The data of infected boys
flulist <- c(1,3,7,25,72,222,282,256,233,189,123,70,25,11,4)
time.flu <-c(1:15) 
# # add to the plot
points(flulist~time.flu,pch=16)


# fitting SIR model to data

# slightly more robust SIR model
fit.sir_model = function (current_timepoint, state_values, parameters)
{
  # create state variables (local variables)
  S = state_values [1]        # susceptibles
  I = state_values [2]        # infectious
  R = state_values [3]        # recovered
  beta=parameters[1]
  gamma=parameters[2]
  # compute derivatives
  dS = (-beta * S * I)
  dI = ( beta * S * I) - (gamma * I)
  dR = (gamma * I)
  
  # combine results
  results = c (dS, dI, dR)
  list (results)
}



sse.sir <- function(start.parm, t, sir.data, sir.start) {
  beta = start.parm[1]
  gamma = start.parm[2]
  S0 = sir.start[1]
  I0 = sir.start[2]
  R0 = sir.start[3]
  inits <- c(S=S0, I=I0, R=R0)
  output <- as.data.frame(ode(y=inits,times = t,func= fit.sir_model, parms = c(beta=beta,gamma=gamma)))
  sse <- sum((output$I-sir.data)^2)
}

# fix gamma to make the point
sse.example <- data.frame(beta=1*10^-seq(2,3,length.out=20),gamma=rep(0.5,20),sse=NA)
for (i in 1:dim(sse.example)[1]){
  sse.example$sse[i] <- sse.sir(as.numeric(sse.example[i,]),time.flu,flulist,c(762,1,9))
}
par(mfrow=c(1,1))
plot(sse~beta,data=sse.example, type = 'b',log='xy')

# solve for beta and gamma
start.guess <- c(beta=0.002,gamma=0.5)
(fit.sir <- optim(start.guess,sse.sir,t=time.flu,sir.data=flulist,sir.start=c(762,1,0)))

# rerun with best guess
(refit.sir <- optim(fit.sir$par,sse.sir,t=time.flu,sir.data=flulist,sir.start=c(762,1,0)))

# Plot
plot(flulist~time.flu,type = 'p', pch=15,col='red')
t <- seq(1,15,by=0.05)
fit.pred <- as.data.frame(ode(c(S=762,I=1,R=0),times=t,fit.sir_model,refit.sir$par))   
lines(fit.pred$I~t)


