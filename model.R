library("deSolve")

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
