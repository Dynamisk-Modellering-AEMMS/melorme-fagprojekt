library("deSolve")




=======
growth = function(t, y, params) {
  with(as.list(c(params,y)), {
    ## function values
    pa    = (temp*b*tarm*zigma)
    fa    = temp*b*tarm*zigma
    ka    = temp*b*tarm*zigma
    fptz  = p*temp*zigma
    g     = temp*P/zigma
    h     = temp*F/zigma*K/zigma*V/zigma
    i     = temp*K/zigma*V/zigma
    j     = temp*H*V/zigma
    zigma = P+F+K+V

    pp    = pi/(pi+fi+ki)
    fp    = fi/(pi+fi+ki)
    kp    = ki/(pi+fi+ki)

    ## derivatives
    #1    = pp+fp+kp
    dtarm = fptz -(pi+fi+ki)
    #dtarm  = 0
    dPx   = pi-pa
    dFx   = fi-fa
    dKx   = ki-ka

    dP    = pa*g
    dUx   = g
    dF    = fa-h
    dK    = ka+h-i
    dV    = i-j
    dVx   = j

    return(list(c(dPx, dFx, dKx, dUx, dVx,
                  dtarm,
                  dP, dF, dK, dV)))
  })
}

params = c(
  temp = 20,
  b = 0,

  pp = 0.2,
  fp = 0.3,
  fk = 0.5,
)
