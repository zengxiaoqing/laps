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
        subroutine plot_winds_2d(u,v,interval,size_in
     1          ,imax,jmax,lat,lon,r_missing_data)

!       include 'lapsparms.for'
!       include 'lapsgrid.cmmn'

        logical l_atms

        common /atms/ l_atms

        real*4 u(imax,jmax),v(imax,jmax)
        real*4 lat(imax,jmax),lon(imax,jmax)
        real*4 mspkt
        data mspkt/.518/

!       This variable keeps the barbs away from the boundary
        isize = 0 ! interval + 1

        relsize = 61. / 200.

        do j = 1+isize,jmax-isize,interval
        do i = 1+isize,imax-isize,interval

            alat = lat(i,j)
            alon = lon(i,j)

            if( u(i,j) .ne. r_missing_data
     1    .and. v(i,j) .ne. r_missing_data
     1    .and. abs(u(i,j)) .lt. 1e6               ! Deals with old data
     1    .and. abs(v(i,j)) .lt. 1e6
     1                                                  )then

                call         uv_to_disp(u(i,j),
     1                          v(i,j),
     1                          dir,
     1                          speed)
                spd_kt = speed / mspkt
                call latlon_to_rlapsgrid(alat,alon,lat,lon,imax,jmax
     1                                                  ,ri,rj,istatus)
                call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax
     1                          ,relsize)

            endif


        enddo ! i
        enddo ! j

        if(l_atms)then
            call setusv_dum(2HSR,0) ! Turn off relative addressing for Apollo
        endif

        return
        end
