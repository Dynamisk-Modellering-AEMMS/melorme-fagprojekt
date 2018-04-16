library("deSolve")

growth = function(t, y, params) {
  with(as.list(c(params,y)), {
    ## function values
    zigma = P+FF+K+V
    f = (zigma/80)*(637.27 - 0.305*density)/12/7 # [density] = 1/(dm^2), (zigma/80) er et gæt for hver meget mere større larver spiser
    print(c(P,FF,K,V))
    print(c(P*100/zigma,FF*100/zigma,K*100/zigma,V*100/zigma))
    #f = max(zigma,0)+1
    g = log(exp((P-0.46*zigma))+1) # "gæt", sørger for maks 46% proteinindhold
    print(g)
    h = 0
    i = -0.3485+0.033054*temp #For temp 10-30 grader. Ellers død?
    j = V/(zigma*H)

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
    dUx = g
    dF = fa-h
    dK = ka+h-i
    dV = i-j
    dVx = j

    return(list(c(dPx, dFx, dKx, dUx, dVx,
                  dP, dF, dK, dV)))
  })
}

params = c(
  temp = 20,
  b = 0,
  density = 50,
  H = 0.4,

  pp = 0.2,
  fp = 0.3,
  kp = 0.5
)

initials = c(
  Px = 0, Fx = 0, Kx = 0, Ux = 0, Vx = 0,
  P = 2, FF = 3, K = 5, V = 1
)

sols = ode(initials,c(1:84),growth,params)
plot(sols)
