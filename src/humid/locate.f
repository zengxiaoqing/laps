cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps
cdis
cdis    this software and its documentation are in the public domain and
cdis    are furnished "as is."  the united states government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  they assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis    technical support to users.
cdis
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  if significant modifications or enhancements
cdis    are made to this software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis


      subroutine locate(xx,n,x,j)


c       adapted for the laps moisture routine june 1996
c       author dan birkenheuer
c
c
c
c       theory:
c
c
c       given input array xx, dimensioned by n, and value x
c       the routine returns the index j such that x exists
c       between xx(j) and xx(j+1).
c
c       the routine returns 0 if x does not exist inside the array
c
c       $log: locate.for,v $
c revision 1.2  1996/06/25  19:43:47  birk
c added fsl disclaimer.
c




      dimension xx(n)
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      go to 10
      endif
      j=jl
      return
      end
