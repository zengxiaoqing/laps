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

        subroutine mean_wind(uanl,vanl,topo,imax,jmax,kmax
     1                                  ,umean,vmean
     1                                  ,ustorm,vstorm,istatus)

        logical ltest_vertical_grid

        real*4 umean(imax,jmax),vmean(imax,jmax)                  ! Output
        real*4 ustorm(imax,jmax),vstorm(imax,jmax)                ! Output
        real*4 uanl(imax,jmax,kmax),vanl(imax,jmax,kmax)          ! Input

        real*4 topo(imax,jmax)                                    ! Input

        real*4 sum(imax,jmax)                                     ! Local
        real*4 usum(imax,jmax)                                    ! Local
        real*4 vsum(imax,jmax)                                    ! Local
        integer*4 klow(imax,jmax)                                 ! Local

        if(ltest_vertical_grid('HEIGHT'))then
            khigh = nint(height_to_zcoord(5000.,istatus))
        elseif(ltest_vertical_grid('PRESSURE'))then
            pres_mb = 300.
            pres_pa = pres_mb * 100.
            khigh = nint(zcoord_of_pressure(pres_pa))
        else
            write(6,*)' mean_wind: unknown vertical grid'
            istatus = 0
            return
        endif

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)return

        write(6,*)' Top level of mean wind computation = ',khigh

        do j = 1,jmax
          do i = 1,imax
             klow(i,j) =
     1            max(nint(height_to_zcoord(topo(i,j),istatus)),1)
             if(istatus .ne. 1)then
                 write(6,*)' mean_wind: ERROR in height_to_zcoord'
                 return
             endif
             sum(i,j) = 0.
             usum(i,j) = 0.
             vsum(i,j) = 0.
          enddo ! j
        enddo ! i

        do k = 1,khigh
          do j = 1,jmax
            do i = 1,imax
              if(uanl(i,j,k) .ne. r_missing_data .and.
     1           vanl(i,j,k) .ne. r_missing_data
     1                        .and. k .ge. klow(i,j))then
                sum(i,j) = sum(i,j) + 1.
                usum(i,j) = usum(i,j) + uanl(i,j,k)
                vsum(i,j) = vsum(i,j) + vanl(i,j,k)
              endif
             enddo
          enddo
        enddo

        do j = 1,jmax
          do i = 1,imax

!            Mean wind through the layer
             umean(i,j) = usum(i,j) / sum(i,j)
             vmean(i,j) = vsum(i,j) / sum(i,j)

!            Shear Vector through the layer
             ushear = uanl(i,j,khigh) - uanl(i,j,klow(i,j))
             vshear = vanl(i,j,khigh) - vanl(i,j,klow(i,j))

!            Estimate storm motion of a right moving storm
!            Rotate shear vector by 90 deg, multiply by .15, add to mean wind
             ustorm(i,j) = umean(i,j) ! + 0.15 * vshear
             vstorm(i,j) = vmean(i,j) ! - 0.15 * ushear

          enddo ! i
        enddo ! j

        return
        end
