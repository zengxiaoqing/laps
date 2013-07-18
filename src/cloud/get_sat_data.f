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

        subroutine get_sat_data(i4time,
     1  i4_sat_window,i4_sat_window_offset,                              ! I
     1  imax,jmax,r_missing_data,                                        ! I
     1  l_use_39,l_use_co2,latency_co2,                                  ! I
     1  lat,lon,                                                         ! I
     1  subpoint_lat_clo_s8a,subpoint_lon_clo_s8a,                       ! O
     1  s8a_k,istat_s8a,comment_s8a,                                     ! O
     1  s3a_k,istat_s3a,comment_s3a,                                     ! O
     1  sst_k,istat_sst,comment_sst,                                     ! O
     1  cldtop_co2_pa_a,cloud_frac_co2_a,istat_co2,lstat_co2_a)          ! O

!       Input
        real lat(imax,jmax)
        real lon(imax,jmax)

!       Output
        real s8a_k(imax,jmax)
        real s3a_k(imax,jmax)
        real sst_k(imax,jmax)
        real cldtop_co2_pa_a(imax,jmax)
        real cloud_frac_co2_a(imax,jmax)

        logical lstat_co2_a(imax,jmax)
        logical l_use_39, l_use_co2

!       Local
        real pct_pa(imax,jmax)
        real lca(imax,jmax)
        real    subpoint_lat_clo_s8a(imax,jmax)
        real    subpoint_lon_clo_s8a(imax,jmax)
        real    subpoint_lat_clo_s3a(imax,jmax)
        real    subpoint_lon_clo_s3a(imax,jmax)

        character*3 lvd_ext
        data lvd_ext /'lvd'/

        character*31 ext
        character var*3,units*10
        character*125 comment_s8a,comment_s3a,comment_sst,comment

        write(6,*)' Subroutine get_sat_data...'

        write(6,*)' Getting IR satellite data from LVD file'
        ext = lvd_ext
        var = 'S8A'
        ilevel = 0
        call get_laps_2dvar(i4time+i4_sat_window_offset,i4_sat_window       
     1                     ,i4time_s8a,lat,lon
     1                     ,subpoint_lat_clo_s8a,subpoint_lon_clo_s8a     ! O 
     1                     ,EXT,var,units
     1                     ,comment_s8a,imax,jmax,s8a_k,ilevel
     1                     ,istat_s8a)      
        if(istat_s8a .ne. 1)then
            write(6,*)' No S8A data available, i4_sat_window is: '
     1               ,i4_sat_window
        endif

!       Final QC check on band 8 (11.2 mm) brightness temps
!       Hopefully, bad/missing  values were filtered out in the creation of LVD
!       Any remaining bad/missing pixels will be evaluated in this QC check
        icount = 0
        do j = 1,jmax
        do i = 1,imax
            if(s8a_k(i,j) .lt. 173. .or. s8a_k(i,j) .gt. 350.
     1                              .or. istat_s8a  .ne. 1      )then
                if(icount .le. 100 .and. istat_s8a .eq. 1)then
                    write(6,*)' Bad LVD/S8A Satellite Brightness '       
     1                       ,'Temperature of'
     1                       ,s8a_k(i,j),' at',i,j
                endif
                icount = icount + 1
                s8a_k(i,j) = r_missing_data
            endif
        enddo
        enddo


!       Obtain Sea Surface Temperature
        write(6,*)' Getting Sea Sfc Temps data from SST file'
        ext = 'sst'
        var = 'SST'
        ilevel = 0
        call get_laps_2dgrid(i4time,3600,i4time_nearest,EXT,var
     1                ,units,comment_sst,imax,jmax,sst_k,ilevel
     1                ,istat_sst)       
        if(istat_sst .ne. 1)then
            write(6,*)' Note: cannot read sst_k'
        endif

        if(l_use_39)then
            write(6,*)' Getting 3.9 micron satellite data from LVD file'
            ext = lvd_ext
            var = 'S3A'
            ilevel = 0
            call get_laps_2dvar(i4time_s8a,0       
     1                         ,i4time_nearest,lat,lon
     1                         ,subpoint_lat_clo_s3a           ! O
     1                         ,subpoint_lon_clo_s3a           ! O 
     1                         ,EXT,var,units
     1                         ,comment_s3a,imax,jmax,s3a_k,ilevel
     1                         ,istat_s3a)
            if(istat_s3a .ne. 1)then
                write(6,*)' No S3A data available'
                s3a_k = r_missing_data
            endif

        else
            write(6,*)' Namelist flag set for not using 3.9u (S3A) data'
            istat_s3a = 0
            s3a_k = r_missing_data

        endif

!       Obtain NESDIS Cloud-top pressure
        if(l_use_co2)then ! This is the "mode 2" type of CO2 usage
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

!       Calculate CO2-Slicing Cloud-top pressure

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

        write(6,*)' Number of utilized CO2-Slicing image points = '
     1           ,icount

        if(l_use_co2)then
            percent_co2_pot = float(icount) / float(imax*jmax) * 100.

            write(6,101)percent_co2_pot
101         format(' CO2-Slicing imagery potentially used over ',f6.2
     1            ,'% of domain')
        endif

        return
        end
