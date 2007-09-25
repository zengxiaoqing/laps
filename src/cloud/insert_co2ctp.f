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
c
c

        subroutine insert_co2ctp(i4time,cld_hts,heights_3d                ! I
     1            ,imax,jmax,kmax,kcloud,r_missing_data                   ! I
     1            ,l_use_co2,latency_co2                                  ! I 
     1            ,default_clear_cover                                    ! I
     1            ,lat,lon,ix_low,ix_high,iy_low,iy_high                  ! I
     1            ,cld_snd,wt_snd,i_snd,j_snd,n_cld_snd,max_cld_snd       ! I/O
     1            ,istatus)                                               ! O

!       Input
        real lat(imax,jmax),lon(imax,jmax)
        real cld_hts(kcloud)
        real heights_3d(imax,jmax,kmax)

!       Arrays for cloud soundings
        real cld_snd(max_cld_snd,kcloud)
        real wt_snd(max_cld_snd,kcloud)
        integer i_snd(max_cld_snd)
        integer j_snd(max_cld_snd)

        logical lstat_co2_a(imax,jmax)
        logical l_use_co2

!       Local
        real pct_pa(imax,jmax)
        real lca(imax,jmax)
        real cldtop_co2_pa_a(imax,jmax)
        real cloud_frac_co2_a(imax,jmax)

        character*31 ext
        character var*3,units*10
        character*125 comment

        write(6,*)' Subroutine insert_c02ctp...'

        write(6,*)
     1    ' Number of cumulative cloud soundings before CO2-slice = '
     1                                                     ,n_cld_snd

!       Obtain NESDIS Cloud-top pressure
        if(l_use_co2)then
            i4_co2_window = latency_co2

            write(6,*)' Getting NESDIS Cloud-top pressure'
            ext = 'ctp'
            var = 'PCT'
            ilevel = 0
            call get_laps_2dgrid(i4time,i4_co2_window,i4time_nearest
     1                      ,EXT,var       
     1                      ,units,comment,imax,jmax,pct_pa,ilevel
     1                      ,istat_pct)       
            if(abs(istat_pct) .ne. 1)then
                write(6,*)' Note: cannot read NESDIS Cloud-top pressure'       
            endif

!           Obtain NESDIS Cloud-fraction
            write(6,*)' Getting NESDIS Cloud-fraction'
            ext = 'ctp'
            var = 'LCA'
            ilevel = 0
            call get_laps_2dgrid(i4time,i4_co2_window,i4time_nearest
     1                          ,EXT,var
     1                          ,units,comment,imax,jmax,lca,ilevel
     1                          ,istat_lca)       
            if(abs(istat_lca) .ne. 1)then
                write(6,*)' Note: cannot read NESDIS Cloud-fraction'
            endif

        endif ! l_use_co2

!       Place CO2-Slicing Cloud-top pressure in sparse array

        cloud_frac_co2_a = r_missing_data
        cldtop_co2_pa_a  = r_missing_data
        lstat_co2_a = .false.
        icount = 0

        if(abs(istat_pct) .eq. 1 .and. abs(istat_lca) .eq. 1 
     1                           .and. l_use_co2                )then
            write(6,*)' Extracting CO2-Slicing info from NESDIS data'
            do j = 1,jmax
            do i = 1,imax
!               Test for partial cloudiness
                if(lca(i,j) .gt. 0. .and. lca(i,j) .lt. 1.0)then ! use co2 data
                    if(pct_pa(i,j) .lt. 100000.)then
                        icount = icount + 1
                        cloud_frac_co2_a(i,j) = lca(i,j)
                        cldtop_co2_pa_a(i,j) = pct_pa(i,j) 
                        lstat_co2_a(i,j) = .true.
                    endif
                endif
            enddo ! i
            enddo ! j

            istat_co2 = 0 ! for now
        else
            istat_co2 = 0
        endif

!       Add to cloud sounding arrays

        n_cloud = 0 ! # of CO2 slicing soundings

        do j = 1,jmax
        do i = 1,imax
            if(lstat_co2_a(i,j))then ! Create this cloud sounding
                call pressure_to_height(cldtop_co2_pa_a(i,j),heights_3d       
     1                                 ,imax,jmax,kmax,i,j,ctop_m
     1                                 ,istatus)       

                cbase_m = ctop_m - cld_thk(ctop_m)
                cbuf_low = cbase_m - 500.

                cover_rpt = cloud_frac_co2_a(i,j)

                if(n_cloud .le. 10)then
                    write(6,*)' Working with CO2 cloud top:'
     1                       ,nint(cldtop_co2_pa_a(i,j))
     1                       ,nint(cbase_m),nint(ctop_m),cover_rpt
     1                       ,' at ',i,j
                endif

                n_cloud = n_cloud + 1

                n_cld_snd = n_cld_snd + 1 ! # of cumulative cloud soundings

                do k=1,kcloud
                    cover = cover_rpt

!                   Fill in Cloud Layer
                    if(cld_hts(k) .ge. cbase_m .and.
     1                 cld_hts(k) .lt. ctop_m                  )then
                        call spread2(cld_snd,wt_snd,i_snd,j_snd
     1                              ,n_cld_snd,max_cld_snd
     1                              ,kcloud,i,j,k,cover,1.)
                        if(n_cloud .le. 10)then
                            write(6,*)' Fill in k = ',k,cover,cld_hts(k)
                        endif
                    endif

                    cover = default_clear_cover

!                   Fill in clear buffer below cloud layer
                    if(cld_hts(k) .lt. cbase_m .and.
     1                 cld_hts(k) .ge. cbuf_low               )then      
                        call spread2(cld_snd,wt_snd,i_snd,j_snd
     1                              ,n_cld_snd,max_cld_snd
     1                              ,kcloud,i,j,k,cover,1.)
                        if(n_cloud .le. 10)then
                            write(6,*)' Fill in k = ',k,cover,cld_hts(k)       
                        endif
                    endif

!                   Fill in clear sky above cloud layer
                    if(cld_hts(k) .ge. ctop_m                 )then      
                        call spread2(cld_snd,wt_snd,i_snd,j_snd
     1                              ,n_cld_snd,max_cld_snd
     1                              ,kcloud,i,j,k,cover,1.)
                        if(n_cloud .le. 10)then
                            write(6,*)' Fill in k = ',k,cover,cld_hts(k)       
                        endif
                    endif

                enddo ! K cld_hts

            endif ! CO2 measurement at this grid point
        enddo ! i
        enddo ! j

        write(6,*)' Number of valid CO2-Slicing soundings = ',n_cloud

        write(6,*)
     1    ' Number of cumulative cloud soundings after CO2-slice = '
     1                                                     ,n_cld_snd

        return
        end
