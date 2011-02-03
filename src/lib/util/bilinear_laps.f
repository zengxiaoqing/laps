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
cdis
cdis
cdis   
cdis

        subroutine bilinear_laps(ri,rj,imax,jmax,array_2d,result)

cdoc    Interpolate 2-d array to find the field value at a fractional grid
cdoc    point.

        real array_2d(imax,jmax)

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
     1     j .ge. 1 .and. j .le. jmax) then

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

cdoc    Interpolate 2-d array to find the field value at a fractional grid
cdoc    point. This one allows you to extrapolate very slightly outside the 
cdoc    grid.

        real array_2d(imax,jmax)

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error in bilinear_interp_extrap'
            return
        endif

        j = nint(rj)
        i = nint(ri)
        if(j .lt. jmax .and. j .gt. 1 .and.
     &     i .lt. imax .and. i .gt. 1)then

           fraci = ri - i
           fracj = rj - j
           if(j.gt.rj)fracj=j-rj
           if(i.gt.ri)fraci=i-ri

c standard bilinear interpolation

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
     &         i .gt. imax .or. i .lt. 1)then
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

        subroutine bilinear_laps_3do(ri_a,rj_a,imax,jmax,nx_laps,ny_laps
     1                             ,nz_laps,array_3d,result)

cdoc    Interpolate 3-d array to find the field values at fractional grid
cdoc    points.

        real array_3d(imax,jmax,nz_laps)                 ! I
        real ri_a(nx_laps,ny_laps),rj_a(nx_laps,ny_laps) ! I
        real result(nx_laps,ny_laps,nz_laps)             ! O

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error in bilinear_laps: STOP'
            stop
        endif

        do k  = 1,nz_laps
        do il = 1,nx_laps
        do jl = 1,ny_laps

          ri = ri_a(il,jl)

          i = int(ri)
          if(i .eq. imax)i=i-1

          rj = rj_a(il,jl)

          j = int(rj)
          if(j .eq. jmax)j=j-1

          if(i .ge. 1 .and. i .le. imax .and.
     1       j .ge. 1 .and. j .le. jmax) then

            fraci = ri - i
            fracj = rj - j

            Z1=array_3d(i  , j  ,k)
            Z2=array_3d(i+1, j  ,k)
            Z3=array_3d(i+1, j+1,k)
            Z4=array_3d(i  , j+1,k)

            if(  z1 .ne. r_missing_data
     1     .and. z2 .ne. r_missing_data
     1     .and. z3 .ne. r_missing_data
     1     .and. z4 .ne. r_missing_data)then

                result(il,jl,k) =  Z1+(Z2-Z1)*fraci+(Z4-Z1)*fracj
     1                          - (Z2+Z4-Z3-Z1)*fraci*fracj

            else
                result(il,jl,k) = r_missing_data

            endif

          else
            result(il,jl,k) = r_missing_data

          endif

        enddo ! j
        enddo ! i  
        enddo ! k

        return
        end
 
        subroutine bilinear_laps_3df(ri_a,rj_a,imax,jmax,nx_laps,ny_laps
     1                             ,nz_laps,array_3d,result)

cdoc    Interpolate 3-d array to find the field values at fractional grid
cdoc    points.

        real array_3d(imax,jmax,nz_laps)                 ! I
        real ri_a(nx_laps,ny_laps),rj_a(nx_laps,ny_laps) ! I
        real result(nx_laps,ny_laps,nz_laps)             ! O

        real z1(nz_laps),z2(nz_laps),z3(nz_laps),z4(nz_laps)

        write(6,*)' Subroutine bilinear_laps_3df...'

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error in bilinear_laps: STOP'
            stop
        endif

        do jl = 1,ny_laps
        do il = 1,nx_laps

          ri = ri_a(il,jl)

          i = int(ri)
          if(i .eq. imax)i=i-1

          rj = rj_a(il,jl)

          j = int(rj)
          if(j .eq. jmax)j=j-1

          if(i .ge. 1 .and. i .le. imax .and.
     1       j .ge. 1 .and. j .le. jmax) then

            fraci = ri - i
            fracj = rj - j
            fracij = fraci*fracj

            Z1(:)=array_3d(i  , j  ,:)
            Z2(:)=array_3d(i+1, j  ,:)
            Z3(:)=array_3d(i+1, j+1,:)
            Z4(:)=array_3d(i  , j+1,:)

!           result(il,jl,:) 
!    1                  =  Z1(:)+(Z2(:)-Z1(:))*fraci+(Z4(:)-Z1(:))*fracj
!    1                  - (Z2(:)+Z4(:)-Z3(:)-Z1(:))*fracij

            result(il,jl,:) 
     1                  =  Z1(:) * (1.0 - fraci - fracj + fracij)             
     1                  +  Z2(:) * (      fraci         - fracij) 
     1                  +  Z3(:) * (                      fracij) 
     1                  +  Z4(:) * (              fracj - fracij) 

            where(abs(result(il,jl,:)) .gt. 1e10)
     1          result(il,jl,:) = r_missing_data
  
          else
            result(il,jl,:) = r_missing_data

          endif

        enddo ! i  
        enddo ! j

        return
        end

        subroutine bilinear_laps_3d(ri_a,rj_a,imax,jmax,nx_laps,ny_laps
     1                             ,nz_laps,array_3d,result)

cdoc    Interpolate 3-d array to find the field values at fractional grid
cdoc    points. Hopefully it's more efficient to call bilinear_laps_2d.

        real array_3d(imax,jmax,nz_laps)                 ! I
        real ri_a(nx_laps,ny_laps),rj_a(nx_laps,ny_laps) ! I
        real result(nx_laps,ny_laps,nz_laps)             ! O

        do k = 1,nz_laps
           call bilinear_laps_2d(ri_a,rj_a,imax,jmax,nx_laps,ny_laps
     1                          ,array_3d(1,1,k),result(1,1,k))
        enddo ! k

        return
        end

        subroutine bilinear_laps_2d(ri_a,rj_a,imax,jmax,nx_laps,ny_laps
     1                             ,array_2d,result)

cdoc    Interpolate 2-d array to find the field values at fractional grid
cdoc    points.

        real array_2d(imax,jmax)                         ! I
        real ri_a(nx_laps,ny_laps),rj_a(nx_laps,ny_laps) ! I
        real result(nx_laps,ny_laps)                     ! O

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error in bilinear_laps_2d: STOP'
            stop
        endif

        do il = 1,nx_laps
        do jl = 1,ny_laps

          ri = ri_a(il,jl)

          i = int(ri)
          if(i .eq. imax)i=i-1

          rj = rj_a(il,jl)

          j = int(rj)
          if(j .eq. jmax)j=j-1

          if(i .ge. 1 .and. i .le. imax .and.
     1       j .ge. 1 .and. j .le. jmax) then

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

                result(il,jl) =  Z1+(Z2-Z1)*fraci+(Z4-Z1)*fracj
     1                          - (Z2+Z4-Z3-Z1)*fraci*fracj

            else
                result(il,jl) = r_missing_data

            endif

          else
            result(il,jl) = r_missing_data

          endif

        enddo ! j
        enddo ! i  

        return
        end
