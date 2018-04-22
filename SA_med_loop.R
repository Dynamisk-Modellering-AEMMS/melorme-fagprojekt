library("deSolve")

growth = function(t, y, params) {
  with(as.list(c(params,y)), {
    ## function values
    zigma = P+FF+K+V
    f = (zigma/80)*(637.27 - 0.305*params[2])/12/7 # [density] = 1/(dm^2), (zigma/80) er et gæt for hver meget mere større larver spiser
    g = max((P-0.46*zigma),0) # "gæt", sørger for maks 46% proteinindhold
    h = 0
    i = -0.3485+0.033054*params[1] #For temp 10-30 grader. Ellers død?
    radius = (zigma/(10*pi))^(1/3) # approx radius (as if cylinder of water 5x as long as the diameter)
    j = 3.4 * #(0.22 * pi * radius^2) * # surface area
      1440 * # minutes/day
      248.33265 * # constant
      (V/zigma - params[3]) * # humidity difference, diffusion boundary assumed linear, replaces (1-H) in eq
      params[1]^(-1.4) * sqrt(params[5]) * Pwater(params[1])

    pi = params[8]*f
    fi = params[9]*f
    ki = params[10]*f

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

Pwater = function(params) { # Antoine ligningen
  10^(
    2.124903 # for at omregne til Pascal
    + 8.07131 # A
    - (
      1730.63/ # B
        (233.426 # C
         + params[1]) # [temp] = °C
    )
  )
}

params = c(
  temp = 20,
  density = 50,
  H = 0.3,

  pressure = 101325, # Pa
  vAir = 0.15, # m/s
  Dp = 1,
  Df = 1,

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

# Scatterplot
mi = c(10,10,0,10,10,10,10,10,0,0,0)
ma = c(100,100,1,100,100,100,100,100,1,1,1)
attribute.names = c('temp','density','H','pressure','vAir','Dp','Df','Dk','pp','fp','kp')
N.iter <- 12
P.max <- matrix(data = NA,nrow = length(params),ncol = N.iter)
FF.max <- matrix(data = NA,nrow = length(params),ncol = N.iter)
K.max <- matrix(data = NA,nrow = length(params),ncol = N.iter)
param.list = matrix(data = NA,nrow = length(params),ncol = N.iter)
#parms = matrix(data = NA,nrow = length(params),ncol = N.iter)
parms = vector("numeric", length = length(params))
GSA.res = matrix(data = NA,nrow = length(params),ncol = N.iter)
#output = matrix(data = NA,nrow = length(params),ncol = N.iter)
j = 1
#loop param
set.seed(3)
for (j in 1:length(params)){
  param.list[j,] = runif(N.iter,min=mi[j],max=ma[j])

  for ( i in 1:N.iter){
    # simulate the epidemic
    parms[j] = param.list[j,i]
    output <- ode(initials,times = c(1:84), growth, parms)
    # antal syge p� samme tid
    P.max[j,i] <- max(output[,'P'])
    # t.max
    FF.max[j,i] <- max(output[,'FF'])
    # totale antal syge
    K.max[j,i] <- max(output[,'K'])
  }
  GSA.res <- data.frame(parm1=param.list[j,],P.max[j,],FF.max[j,],K.max[j,])

  pairs(GSA.res)
}
# Andre metoder f�lger...

