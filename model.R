library("deSolve")




=======
growth = function(t, y, params) {
  with(as.list(c(params,y)), {
    ## function values
    
    # compute derivatives
    pa    = (t*b*tarm*zigma)
    fa    = t*b*tarm*zigma
    ka    = t*b*tarm*zigma
    fptz  = p*t*zigma
    g     = t*P/zigma
    h     = t*f/zigma*k/zigma*v/zigma
    i     = t*k/zigma*v/zigma
    j     = t*H*v/zigma
    zigma = p+f+k+v

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
    dk    = ka+h-i
    dv    = i-j
    dvx   = j
    
    return(list(c(dPx, dFx, dKx, dUx, dVx,
                  dtarm,
                  dP, dF, dK, dV)))
  })
}
