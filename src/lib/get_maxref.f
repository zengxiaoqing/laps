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

        subroutine get_max_ref(grid_ra_ref,imax,jmax,kmax,radar_array)

!       This routine converts 3D reflectivity volume data to a 2D reflectivity
!       field. This is done by vertically projecting the maximum reflectivity
!       onto a horizontal surface.

!       Steve Albers            1990
!       Steve Albers            1994     Test for missing data values

        real*4 grid_ra_ref(imax,jmax,kmax)      ! Input 3D Array
        real*4 radar_array(imax,jmax)           ! Output 2D Array

        common /laps_diag/ no_laps_diag

        if(no_laps_diag .eq. 0)
     1   write(6,*)' Projecting maximum reflectivity onto horizontal sur
     1face'

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error in get_max_ref, STOP'
            stop
        endif

        call get_ref_base(ref_base,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error in get_max_ref, STOP'
            stop
        endif

!       Initialize Radar Array
        do j = 1,jmax
        do i = 1,imax
            radar_array(i,j) = ref_base
        enddo
        enddo

        do k = 1,kmax
        do j = 1,jmax
        do i = 1,imax
            if(grid_ra_ref(i,j,k) .ne. r_missing_data)
     1       radar_array(i,j) = max(radar_array(i,j),grid_ra_ref(i,j,k))
        enddo
        enddo
        enddo

        return
        end

