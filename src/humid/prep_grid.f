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


        subroutine prep_grid(m,n,data,pn_max,points,pn)

c       $log: prep_grid.for,v $
c revision 1.1  1996/08/30  20:48:53  birk
c initial revision
c

        implicit none

        integer m,n,pn,pn_max

        real points(pn_max,3),data(m,n),weight_t,dist,weight


        integer i,j,k


        do i = 1,m,m-1
        do j = 1,n


        weight_t = 0.0
        data(i,j) = 0.0


        do k = 1,pn

                if (points(k,2).eq.i  .and. points(k,3) .eq.j) then
                        data(i,j) = points(k,1)
                        go to 22
                elseif (points(k,2) .ne. 0) then

                dist = sqrt( (i-points(k,2))**2+(j-points(k,3))**2)
                weight = 1./dist
                data(i,j) = data(i,j) + weight*points(k,1)
                weight_t = weight_t + weight


                endif

        enddo !k

        if(weight_t.ne.0.) then
        data(i,j) = data(i,j) / weight_t
        else
        write(6,*) 'error in divide'
        endif

22      continue

        enddo !j
        enddo !k



        do i = 1,m
        do j = 1,n,n-1


        weight_t = 0.0
        data(i,j) = 0.0


        do k = 1,pn

                if (points(k,2).eq.i  .and. points(k,3) .eq.j) then
                        data(i,j) = points(k,1)
                        go to 23
                elseif (points(k,2).ne.0) then

                dist = sqrt( (i-points(k,2))**2+(j-points(k,3))**2)
                weight = 1./dist
                data(i,j) = data(i,j) + weight*points(k,1)
                weight_t = weight_t + weight


                endif

        enddo !k

        if(weight_t.ne.0.) then
        data(i,j) = data(i,j) / weight_t
        else
        write (6,*) 'error in divide'
        endif

23      continue

        enddo !j
        enddo !k





        return
        end

