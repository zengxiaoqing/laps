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



        subroutine slv_laplc (data,mask, nx, ny)

c       $log: slv_laplc.for,v $
c revision 1.1  1996/08/30  20:57:55  birk
c initial revision
c

        implicit none

        integer nx,ny
        integer mask(nx,ny)
        real data(nx,ny),error
        real maxerror
        real typical_data

        integer i,j,k



        do k = 1,6000

        maxerror =0.0


        do j = 2,ny -1
           do i = 2,nx -1

              if(mask(i,j) .eq. 0) then

                 error = 0.25 * ( data(i+1,j) + data(i-1,j) +
     1                data(i,j+1) + data (i,j-1) ) - data (i,j)

                 data (i,j) = error + data(i,j)

                 maxerror = max(maxerror,abs(error))

              else 
                 typical_data = data(i,j)


              endif

           enddo
        enddo

c       print*, maxerror

        if(typical_data.ne.0.0) then
           if (maxerror/typical_data.le. 1.e-3) go to 22
        else
           if (maxerror.le. 1.e-3) then
             write(6,*) 'typical_data = 0.0, divide avoided'
             go to 22
           endif
        endif


        enddo

        write(6,*) 'diriclet terminated on iterations ', k

22      write (6,*) 'max error solving dirichlet problem ', 
     1       maxerror/typical_data,k

        return
        end

