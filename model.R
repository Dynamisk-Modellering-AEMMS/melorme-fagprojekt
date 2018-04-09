library("deSolve")

growth = function(t, y, params) {
  with(as.list(c(params,y)), {
    ## function values

    ## derivatives

    return(list(c(dPx, dFx, dKx, dUx, dVx,
                  dtarm,
                  dP, dF, dK, dV)))
  })
}
