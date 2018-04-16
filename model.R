library("deSolve")

growth = function(t, y, params) {
  with(as.list(c(params,y)), {
    ## function values
    zigma = P+FF+K+V
    f     = K/1000*(143.23-.144*density)      # food consumed vs.Density pr dm^2
    g     = log(exp(P-200)+1)
    h     = 0
    i     = -0.3485+0.033054*temp #For temp 10-30 grader. Ellers død?
    j     = V/(zigma*H)

    pi = pp*f
    fi = fp*f
    ki = kp*f

    pa    = pi
    fa    = fi
    ka    = ki

    ## derivatives
    dPx   = pi-pa
    dFx   = fi-fa
    dKx   = ki-ka

    dP    = pa-g
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
