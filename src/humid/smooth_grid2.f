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



        subroutine smooth_grid2(m,n,data,range)


c       $log: smooth_grid2.for,v $
c revision 1.1  1996/08/30  20:59:01  birk
c initial revision
c

        implicit none

        integer m,n,range
        real data(m,n)
        real sum


        real counter
        integer i,j
        integer k,l



        do i = 1,m
        do j = 1,n

        sum = 0.0


        counter= 0.0

        do k=i-range,i+range
        do l=j-range,j+range
        if ( (k.gt.0  .and. l.gt.0) .and. (k.le.m .and. l.le.n) ) then
                counter = counter+1.
                sum = sum + data(k,l)
        endif
        enddo
        enddo

        data(i,j) = (sum + data(i,j))/(counter+1.)

        enddo
        enddo

        return
        end

