cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 

        subroutine bilinear_laps(ri,rj,imax,jmax,array_2d,result)

        real*4 array_2d(imax,jmax)

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error in bilinear_laps: STOP'
            stop
        endif

        i = int(ri)
        if(i .eq. imax)i=i-1

        j = int(rj)
        if(j .eq. jmax)j=j-1

        if(i .ge. 1 .and. i .le. imax .and.
     1   j .ge. 1 .and. j .le. jmax) then

            fraci = ri - i
            fracj = rj - j

            Z1=array_2d(i  , j  )
            Z2=array_2d(i+1, j  )
            Z3=array_2d(i+1, j+1)
            Z4=array_2d(i  , j+1)

            if(  z1 .ne. r_missing_data
     1     .and. z2 .ne. r_missing_data
     1     .and. z3 .ne. r_missing_data
     1     .and. z4 .ne. r_missing_data)then

                result =  Z1+(Z2-Z1)*fraci+(Z4-Z1)*fracj
     1                - (Z2+Z4-Z3-Z1)*fraci*fracj

            else
                result = r_missing_data

            endif

        else
            result = r_missing_data

        endif

        return
        end

        subroutine bilinear_interp_extrap(ri,rj,imax,jmax
     1                                   ,array_2d,result,istatus)

        real*4 array_2d(imax,jmax)

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error in bilinear_laps_extrap'
            return
        endif

        j = nint(rj)
        i = nint(ri)
        if(j .lt. jmax .and. j .gt. 1 .and.
     &     i .lt. imax .and. i .gt. 1)then
c standard bilinear interpolation
              fraci = ri - i
              fracj = rj - j

              Z1=array_2d(i  , j  )
              Z2=array_2d(i+1, j  )
              Z3=array_2d(i+1, j+1)
              Z4=array_2d(i  , j+1)

              if(    z1 .ne. r_missing_data
     1         .and. z2 .ne. r_missing_data
     1         .and. z3 .ne. r_missing_data
     1         .and. z4 .ne. r_missing_data)then

                  result =  Z1+(Z2-Z1)*fraci+(Z4-Z1)*fracj
     1                   - (Z2+Z4-Z3-Z1)*fraci*fracj
              else
                  result = r_missing_data
              endif
        elseif(j .gt. jmax .or. j .lt. 1 .or.
     &            i .gt. imax .or. i .lt. 1)then
                  result = r_missing_data
        elseif(j .eq. jmax .or. j .eq. 1)then
              if(i .eq. imax .or. i .eq. 1)then
                 result = array_2d(i,j)
              else
                 frac1 = 1-(ri-int(ri))
                 frac2 = 1-frac1
                 Z1=array_2d(i  , j  )
                 Z2=array_2d(i+1, j  )
                 result = Z1*frac1+Z2*frac2
              endif
        elseif(i .eq. imax .or. i .eq. 1)then
              if(j .eq. jmax .or. j .eq. 1)then
                 result = array_2d(i,j)
              else
                 frac1 = 1-(rj-int(rj))
                 frac2 = 1-frac1
                 Z1=array_2d(i  , j  )
                 Z2=array_2d(i, j+1  )
                 result = Z1*frac1+Z2*frac2
              endif
        endif

        istatus = 1

        return
        end

