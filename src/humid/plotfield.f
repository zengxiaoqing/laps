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
cdis     technical support to users.
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


        subroutine plotfield (data_anal,igrid,jgrid)


c       $log: plotfield.for,v $
c revision 1.1  1996/08/30  20:48:28  birk
c initial revision
c

c ported to the unix environment 1 oct 1993

        implicit none

        integer igrid,jgrid


        real data_anal(igrid,jgrid)
        integer i,j
        real finc
        integer pn



        integer ioffp,ioffm
        real spval,epsval,cntmin,cntmax,cntint

        common /conre1/ioffp,spval,epsval,cntmin,cntmax,cntint,ioffm

        ioffp = 1   ! says that the missing value exists.
        spval = 1.0e+37  ! sets the missing value flag for plots
        ioffm = 0 !omit the conrec message on the plot -- doesn't work





c        call opngks


        ioffm = 1




        call conrec (data_anal,igrid,igrid,jgrid,0.,0.,finc,-1,0,0)
        call frame

c        call clsgks


        end
