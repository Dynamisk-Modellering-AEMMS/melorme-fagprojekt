library("deSolve")

growth = function(t, y, params) {
  with(as.list(c(params,y)), {
    ## function values
    zigma = P+FF+K+V
    f = (zigma/80)*(637.27 - 0.305*density)/12/7 # [density] = 1/(dm^2), (zigma/80) er et gÃ¦t for hver meget mere stÃ¸rre larver spiser
    g = 0.125 * P # "gÃ¦t", sÃ¸rger for maks 46% proteinindhold
    kCombustC = 0.08473955701/0.5
    h = kCombustC/2*K
    #kCombust = -0.3485+0.033054*temp #For temp 10-30 grader. Ellers dÃ¸d?
    kCombust = kCombustC/2*K
    fCombust = 0.05298229217*FF
    radius = (zigma/(10*pi))^(1/3) # approx radius (as if cylinder of water 5x as long as the diameter)
    j = (0.22 * pi * radius^2) * # surface area
      1440 * # minutes/day
      0.89750342 *
      (V/zigma - H) * # humidity difference, diffusion boundary assumed linear, replaces (1-H) in eq
      (temp + 273.15)^(-1.4) * sqrt(vAir) * Pwater(temp)

    lumenSize = 0.15*zigma

    ppl = Pl/lumenSize
    fpl = Fl/lumenSize
    kpl = Kl/lumenSize

    pa = Dp * Pl/lumenSize
    fa = Df * Fl/lumenSize
    ka = Dk * Kl/lumenSize # / lumenSize [mg]

    ## derivatives
    dPx = ppl*f
    dFx = fpl*f
    dKx = kpl*f

    dP = pa-g
    dUx = 0.5*g
    dF = fa + h - fCombust
    dK = ka - h - kCombust
    #dV = 0.5*g + kCombust + fCombust - j
    dV = 0.5*g + (0.02108 * zigma) - j # (kCombust + fCombust) = 0.021808 uden skelnen, desvÃ¦rre
    dVx = j

    dPl = pp*f - ppl*f - pa
    dFl = fp*f - fpl*f - fa
    dKl = kp*f - kpl*f - ka
    #print(j)

    return(list(c(dPx, dFx, dKx, dUx, dVx,
                  dP, dF, dK, dV,
                  dPl, dFl, dKl, f, pa+fa+ka),
                c(f = f, pa = pa, fa = fa, ka = ka, zigma = zigma, kCombust = kCombust, fCombust = fCombust, j = j, pp = 100*P/zigma, kp = 100*K/zigma, fp = 100*FF/zigma, vp = 100*V/zigma, ppt = 100*P/(P+K+FF), kpt = 100*K/(P+K+FF), fpt = 100*FF/(P+K+FF), urp = Ux/(Px+Fx+Kx+Ux), growthrate=dK+dF+dK)))
  })
}

Pwater = function(temp) { # Antoine ligningen
  10^(
    2.124903 # for at omregne til Pascal
    + 8.07131 # A
    - (
      1730.63/ # B
        (233.426 # C
         + temp) # [temp] = Â°C
    )
  )
}

params = c(
  temp = 28,
  density = 50,
  H = 0.65,

  pressure = 101325, # Pa
  vAir = 0.15, # m/s
  Dp = 100, # gÃ¦t
  Df = 100, # gÃ¦t
  Dk = 100, # (0.015umol/min)/(3 mol/L) * 1 kg/L omregnet til mg/d

  pp = 0.33,
  fp = 0.06,
  kp = 0.61
)

initials = c(
  Px = 0, Fx = 0, Kx = 0, Ux = 0, Vx = 0,
  P = 2, FF = 3, K = 5, V = 5,
  Pl = 0.2, Fl = 0.3, Kl = 0.5,
  FoodConsumed = 0, FoodDigested = 0
)

sols = ode(initials,c(1:84),growth,params)
plot(sols)
print(sols[84,])

# Scatterplot
mi <- as.vector(c(0,10,0,101325/2,0.001,1,1,1,0,0,0),"numeric")
ma <- as.vector(c(100,100,1,101325*2,10,100,100,100,1,1,1),"numeric")
attribute.names = c('temp','density','H','pressure','vAir','Dp','Df','Dk','pp','fp','kp')
N.iter <- 500
P.end <- vector("numeric",length = N.iter)
FF.end <- vector("numeric",N.iter)
K.end <- vector("numeric",N.iter)
param.list = matrix(data = NA,nrow = length(params),ncol = N.iter)
#parms = matrix(data = NA,nrow = length(params),ncol = N.iter)
parms = vector("numeric", length = length(params))
GSA.res = matrix(data = NA,nrow = length(params),ncol = N.iter)
#output = matrix(data = NA,nrow = length(params),ncol = N.iter)
#j = 1
#loop param
set.seed(3)
for (j in 1:length(params)){
  param.list[j,] = runif(N.iter,min=mi[j],max=ma[j])
}
for ( i in 1:N.iter){
  # simulate the epidemic
  parms = param.list[,i]
  names(parms)<-attribute.names
  output <- ode(initials,times = c(1:84), growth, parms)
  # antal syge på samme tid
  P.end[i] <- output[84,'ppt']
  # t.max
  FF.end[i] <- output[84,'fpt']
  # totale antal syge
  K.end[i] <- output[84,'kpt']
}
GSA.res <- data.frame(id=1:N.iter,P.end,FF.end,K.end)
names(GSA.res) <- c("id","P end","F end", "K end")

pairs(GSA.res)

# Andre metoder følger...

library(ggplot2)
library(reshape2)
#head(GSA.res)
#GSA.res[1] <- as.factor(P.max)
nyGSA = melt(GSA.res,id.vars=c("id"))
nyGSA$variable <- as.factor(nyGSA$variable)
ggplot(nyGSA, aes(x=variable,y=value)) +
  geom_violin(trim=T, fill='pink', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle("Violinplot")+xlab("Slut næringsindhold")+ylab("Værdier")

