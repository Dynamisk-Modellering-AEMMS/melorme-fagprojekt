library("deSolve")

growth = function(t, y, params) {
  with(as.list(c(params,y)), {
    ## function values
    zigma = P+F+K+V
#    temp  = -0.3485+0.033054*t #For temp 10-30 grader. Ellers død?
    rho   = 143.23-.144*t      # food consumed vs.Density pr dm^2
temp = 20
    f     = density+temp+zigma
    g     = temp+P/zigma
    h     = temp+F/zigma+K/zigma+V/zigma
    i     = temp+K/zigma+V/zigma
    j     = temp+H+V/zigma

    pi = pp*f
    fi = fp*f
    ki = kp*f

    pa    = temp+b+zigma
    fa    = temp+b+zigma
    ka    = temp+b+zigma

    ## derivatives
    dPx   = pi-pa
    dFx   = fi-fa
    dKx   = ki-ka

    dP    = pa+g
    dUx   = g
    dF    = fa-h
    dK    = ka+h-i
    dV    = i-j
    dVx   = j

    return(list(c(dPx, dFx, dKx, dUx, dVx,
                  dP, dF, dK, dV)))
  })
}

params = c(
  temp = 20,
  b = 0,
  density = 10,
  H = 0.4,

  pp = 0.2,
  fp = 0.3,
  kp = 0.5
)

initials = c(
  Px = 0, Fx = 0, Kx = 0, Ux = 0, Vx = 0,
  P = 1, FF = 1, K = 1, V = 1
)

sols = ode(initials,c(1:30),growth,params)
plot(sols)
