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


       subroutine  check_nan2 (var2,n,m,istatus)

c  this routine checks for Not a Number conditions (NaN).  It returns an
c  istatus of 0 if NaNs are detected. Otherwise it returns a value of 1.

c  author: Dan Birkenheuer
c  date:   12/11/96
c  date:    1/15/97  modified into 4 routines Brandy can't handle 
c                    entry points.
       integer istatus
       real var2(n,m)
       integer i,j
       integer m,n



c  double dimension

        do j = 1,m
        do i = 1,n
        if(var2(i,j).ne.var2(i,j)) then ! NaN detected
         istatus = 0
         return
        endif
        enddo
        enddo


        istatus = 1
        return

        end


