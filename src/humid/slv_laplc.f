cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis cdis
cdis
cdis
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

        if(typical_data > 1.e-9) then
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
     1       maxerror,'/',typical_data,k

c     fill boarders (normally zero) with nearest neighbor values

        do j = 1,ny
           do i = 1,nx
              if(i .eq. 1 ) then !boarder
                 data(i,j) = data(i+1, j)
              elseif(i .eq. nx) then
                 data(i,j) = data(i-1,j)
              endif
           enddo
        enddo
     
        do j = 1,ny
           do i = 1,nx
              if(j .eq. 1 ) then !boarder
                 data(i,j) = data(i, j+1)
              elseif(j .eq. ny) then
                 data(i,j) = data(i,j-1)
              endif
           enddo
        enddo
                   

        return
        end

