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
     1  s8a_k,istat_s8a,                                                 ! O
     1  s3a_k,istat_s3a,                                                 ! O
     1  sst_k,istat_sst)                                                 ! O

!       Output
        real*4 s8a_k(imax,jmax)
        real*4 s3a_k(imax,jmax)
        real*4 sst_k(imax,jmax)

        character*3 lvd_ext
        data lvd_ext /'lvd'/

        character*31 ext
        character var*3,comment*125,units*10

        write(6,*)' Subroutine get_sat_data...'

        write(6,*)' Getting IR satellite data from LVD file'
        ext = lvd_ext
        var = 'S8A'
        ilevel = 0
        call get_laps_2dvar(i4time+i4_sat_window_offset,i4_sat_window       
     1                     ,i4time_s8a,EXT,var,units
     1                     ,comment,imax,jmax,s8a_k,ilevel,istat_s8a)
        if(istat_s8a .ne. 1)then
            write(6,*)' No S8A data available'
        endif

!       Final QC check on band 8 (11.2 mm) brightness temps
!       Hopefully, bad/missing  values were filtered out in the creation of LVD
!       Any remaining bad/missing pixels will be evaluated in this QC check
        icount = 0
        do j = 1,jmax
        do i = 1,imax
            if(s8a_k(i,j) .lt. 173. .or. s8a_k(i,j) .gt. 350.)then
                if(icount .le. 100)then
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
     1                ,units,comment,imax,jmax,sst_k,ilevel,istat_sst)       
        if(istat_sst .ne. 1)then
            write(6,*)' Note: cannot read sst_k'
        endif

        write(6,*)' Getting 3.9 micron satellite data from LVD file'
        ext = lvd_ext
        var = 'S3A'
        ilevel = 0
        call get_laps_2dvar(i4time_s8a,0       
     1                     ,i4time_nearest,EXT,var,units
     1                     ,comment,imax,jmax,s3a_k,ilevel,istat_s3a)
        if(istat_s3a .ne. 1)then
            write(6,*)' No S3A data available'
        endif

        return
        end
