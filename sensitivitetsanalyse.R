#Sensitivitets analyse
# SIR model in R
#

# get the ODE solver
library(deSolve)

# make the function of the system

growth = function(t, y, params) {
  with(as.list(c(params,y)), {
    ## function values
    zigma = P+FF+K+V
    f = (zigma/80)*(637.27 - 0.305*density)/12/7 # [density] = 1/(dm^2), (zigma/80) er et gæt for hver meget mere større larver spiser
    g = log(exp((P-0.46*zigma))+1) # "gæt", sørger for maks 46% proteinindhold
    h = 0
    i = -0.3485+0.033054*temp #For temp 10-30 grader. Ellers død?
    radius = (zigma/(10*pi))^(1/3) # approx radius (as if cylinder of water 5x as long as the diameter)
    j = 3.4 * #(0.22 * pi * radius^2) * # surface area
      1440 * # minutes/day
      248.33265 * # constant
      (V/zigma - H) * # humidity difference, diffusion boundary assumed linear, replaces (1-H) in eq
      temp^(-1.4) * sqrt(vAir) * Pwater(temp)

    pi = pp*f
    fi = fp*f
    ki = kp*f

    pa = pi
    fa = fi
    ka = ki

    ## derivatives
    dPx = pi-pa
    dFx = fi-fa
    dKx = ki-ka

    dP = pa-g
    dUx = 0.5*g
    dF = fa-h
    dK = ka+h-i
    dV = 0.5*g + i - j
    dVx = j

    return(list(c(dPx, dFx, dKx, dUx, dVx,
                  dP, dF, dK, dV)))
  })
}

Pwater = function(temp) { # Antoine ligningen
  10^(
    2.124903 # for at omregne til Pascal
    + 8.07131 # A
    - (
      1730.63/ # B
        (233.426 # C
         + temp) # [temp] = °C
    )
  )
}

params = c(
  temp = 20,
  b = 0,
  density = 50,
  H = 0.2,
  pressure = 101325, # Pa
  vAir = 0.15, # m/s

  pp = 0.2,
  fp = 0.4,
  kp = 0.4
)

initials = c(
  Px = 0, Fx = 0, Kx = 0, Ux = 0, Vx = 0,
  P = 2, FF = 3, K = 5, V = 1
)

sols = ode(initials,c(1:84),growth,params)
plot(sols)

# LOkal sensitivitet
# juster hver parameter op eller ned fra mean
params <- c(temp.value=temp*1.25,
            b = 0,
            density.value = density*1.25,
            H = 0.2,
            pressure = 101325, # Pa
            vAir = 0.15,

            pp = 0.2,
            fp = 0.4,
            kp = 0.4
) # m/s)

output <- ode(y=initials,times = c(1:84),func= growth, parms = params)
# antal syge på samme tid
(FF.max <- max(output[,'FF']))
# t.max
#t.max <- output[which(output[,'t'] == i.max),'time'])
# totale antal syge
#(r.max <- max(output[,'R']))




# Scatterplot
iter = 1000
temp.list <- runif(1000,min=10,max=30)
density.list <- runif(1000,min=10,max=100)
FF.max <- numeric(1000)
#i.max <- numeric(iter)
#r.max <- numeric(iter)

for ( i in 1:1000){
  # simulate the epidemic
  params <- c(temp=temp.list[i],
              b = 0,
              denisity = density.list[i],
              H = 0.2,
              pressure = 101325,
              vAir = 0.15,
              pp = 0.2,
              fp = 0.4,
              kp = 0.4
  )
  output <- ode(y=initials,times = c(1:84),func= growth, parms = params)
  # antal syge på samme tid
  FF.max[i] <- max(output[,'FF'])
  # t.max
  #t.max[i] <- output[which(output[,'I'] == i.max[i]),'time']
  # totale antal syge
  #r.max[i] <- max(output[,'R'])
}

GSA.res <- data.frame(temp=temp.list,denistet=density.list,FF.max)

# Scatterplots
pairs(GSA.res)


# PRCC

library(epiR)

epi.prcc(GSA.res[,c(1,2,3)]) # beta/gamma er positiv/negativ korreleret med i.max

epi.prcc(GSA.res[,c(1,2,4)]) # beta/gamma er negativ/svagt positiv korreleret med t.max

epi.prcc(GSA.res[,c(1,2,5)]) # beta/gamma er positiv/negativ korreleret med r.max

