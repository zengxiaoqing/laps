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
        subroutine plot_winds_2d(u,v,interval,size_in,zoom
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

        relsize = size_in

        write(6,*)
        write(6,*) ' Plot_winds_2d: interval/size=',interval,relsize
        write(6,*)
        write(6,*) ' winds are assumed to be GRID north at this point'       

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
     1                          ,relsize,'grid')

            endif


        enddo ! i
        enddo ! j

        if(l_atms)then
            call setusv_dum(2HSR,0) ! Turn off relative addressing for Apollo
        endif

        return
        end
