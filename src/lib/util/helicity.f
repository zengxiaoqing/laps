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

        subroutine helicity_laps(uanl,vanl,ustorm,vstorm,topo
     1  ,u_rel_last ! Local
     1  ,v_rel_last ! Local
     1  ,area_sum    ! Local
     1        ,klow        ! Local
     1  ,imax,jmax,kmax,helicity,istatus)

        include 'lapsparms.inc' ! for vert grid, r_missing_data

        real*4 ustorm(imax,jmax),vstorm(imax,jmax)
        real*4 uanl(imax,jmax,kmax),vanl(imax,jmax,kmax)
        real*4 helicity(imax,jmax)
        real*4 topo(imax,jmax)

!       Local
        real*4 u_rel_last(imax,jmax)
        real*4 v_rel_last(imax,jmax)
        real*4 area_sum(imax,jmax)
        integer*4 klow(imax,jmax)

        icount_write = 0

        write(6,*)' Computing Helicity'
        hel_max = +1000.

        if    (vertical_grid .eq. 'HEIGHT')then
            khigh = nint(height_to_zcoord(5000.,istatus))

        elseif(vertical_grid .eq. 'PRESSURE')then
            pres_mb = 500.
            pres_pa = pres_mb * 100.
            khigh = nint(zcoord_of_pressure(pres_pa))

        endif

        write(6,*)' Top level of helicity computation = ',khigh

        do j = 1,jmax
          do i = 1,imax
             area_sum(i,j) = 0.
             klow(i,j) = max(nint(height_to_zcoord(topo(i,j),istatus)),1
     1)
          enddo ! i
        enddo ! j

        do k = 1,khigh
          do j = 1,jmax
            do i = 1,imax
              if(uanl(i,j,k) .ne. r_missing_data .and.
     1    vanl(i,j,k) .ne. r_missing_data      )then

!               Determine u & v relative to storm mean motion vector
                u_rel = uanl(i,j,k) - ustorm(i,j)
                v_rel = vanl(i,j,k) - vstorm(i,j)

                if(k .gt. klow(i,j))then

!                 Cross product of wind vectors on top and bottom of layer
                  xprod = u_rel * v_rel_last(i,j) - v_rel * u_rel_last(i
     1,j)

!                 Incremental area of hodograph
                  area = .5 * xprod

!                 Total area of hodograph
                  area_sum(i,j) = area_sum(i,j) + area

                endif

                u_rel_last(i,j) = u_rel
                v_rel_last(i,j) = v_rel

              endif ! Missing Data

            enddo ! j

          enddo ! i

        enddo ! k

        do j = 1,jmax
          do i = 1,imax
             helicity(i,j) = area_sum(i,j) /
     1      (zcoord_of_level(khigh) - zcoord_of_level(klow(i,j)))

             if(helicity(i,j) .lt. hel_max)then
                 hel_max = helicity(i,j)
                 if(icount_write .eq. icount_write / 10 * 10)then
                     write(6,101)i,j,klow(i,j),khigh,
     1          helicity(i,j),ustorm(i,j),(uanl(i,j,k),k=1,min(kmax,21))
101                  format(/4i4,f10.5,f8.3/7e11.3/7e11.3/7e11.3)
                 endif
                 icount_write = icount_write + 1
             endif

          enddo ! i
        enddo ! j


        return
        end
