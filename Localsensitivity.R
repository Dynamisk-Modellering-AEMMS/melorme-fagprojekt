library("deSolve")

growth = function(t, y, params) {
  with(as.list(c(params,y)), {
    ## function values
    zigma = P+FF+K+V
    f = (zigma/80)*(637.27 - 0.305*density)/12/7 # [density] = 1/(dm^2), (zigma/80) er et gæt for hver meget mere større larver spiser
    g = max((P-0.46*zigma),0) # "gæt", sørger for maks 46% proteinindhold
    h = 0
    kCombust = -0.3485+0.033054*temp #For temp 10-30 grader. Ellers død?
    fCombust = 0
    radius = (zigma/(10*pi))^(1/3) # approx radius (as if cylinder of water 5x as long as the diameter)
    j = (0.22 * pi * radius^2) * # surface area
      1440 * # minutes/day
      248.33265 * # constant
      (V/zigma - H) * # humidity difference, diffusion boundary assumed linear, replaces (1-H) in eq
      temp^(-1.4) * sqrt(vAir) * Pwater(temp)

    lumenSize = 0.15*zigma

    ppl = Pl/lumenSize
    fpl = Fl/lumenSize
    kpl = Kl/lumenSize

    pa = Dp * Pl/lumenSize
    fa = Df * Fl/lumenSize
    ka = (36/5) * Kl / lumenSize # (0.015umol/min)/(3 mol/L) ⋅ Kl [mg] ⋅ 1 kg/L / lumenSize [mg]

    ## derivatives
    dPx = ppl*f
    dFx = fpl*f
    dKx = kpl*f

    dP = pa-g
    dUx = 0.5*g
    dF = fa + h - fCombust
    dK = ka - h - kCombust
    dV = 0.5*g + kCombust + fCombust - j
    dVx = j

    dPl = pp*f - ppl*f - pa
    dFl = fp*f - fpl*f - fa
    dKl = kp*f - kpl*f - ka
    #print(j)

    return(list(c(dPx, dFx, dKx, dUx, dVx,
                  dP, dF, dK, dV,
                  dPl, dFl, dKl),
                c(f = f, pa = pa, fa = fa, ka = ka, zigma = zigma, kCombust = kCombust, fCombust = fCombust, j = j, pp = P/zigma, kp = K/zigma, fp = FF/zigma, vp = V/zigma)))
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
  P = 2, FF = 3, K = 5, V = 5,
  Pl = 0.2, Fl = 0.3, Kl = 0.5
)

sols = ode(initials,c(1:84),growth,params)
plot(sols)

Pmax = c(NULL,length = length(params)-1)
FFmax = c(NULL,length = length(params)-1)
Kmax = c(NULL,length = length(params)-1)

# LOkal sensitivitet
# juster hver parameter op eller ned fra mean

x = 1.25
n=1
for (n in range(1:length(params))){
  params[n]=params[n]*x
  output <- ode(y=initials,times = c(1:84),func= growth, parms = params)
  # antal syge p� samme tid
  (Pmax[n] <- max(output[,'P']))
  # t.max
  (FFmax[n] <- max(output[,'FF']))
  # totale antal syge
  (Kmax[n]<- max(output[,'K']))

}
