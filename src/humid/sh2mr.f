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



        subroutine sh2mr (sh,mr)


c       $log: sh2mr.for,v $
c revision 1.1  1996/08/30  20:57:36  birk
c initial revision
c

c       given the specific humidity (g/g)
c       this routine computes the mixing ratio (g/g)

        real sh,mr
        real epsilon

        epsilon = 0.622 !(ratio of molecular weight of water ~ 18,
c                       over dry air~29)

        if (sh/epsilon .eq. 1.) then
                write(6,*) 'error in input to routine sh2mr'
                write(6,*) 'divide by zero has occurred'
                write(6,*) 'exit with return value (1)'
                call exit (1)
        endif


        if (sh/epsilon .gt. 1.) then
                write(6,*) 'error in input to routine sh2mr'
                write(6,*) 'impossible condition has occurred'
                write(6,*) 'exit with return value (1)'
                call exit (1)
        endif

        mr = sh / (1.-sh/epsilon)

        return
        end





