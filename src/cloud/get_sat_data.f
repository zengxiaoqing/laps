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
     1                     ,i4time_nearest,EXT,var,units
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
        call get_laps_2dvar(i4time+i4_sat_window_offset,i4_sat_window       
     1                     ,i4time_nearest,EXT,var,units
     1                     ,comment,imax,jmax,s3a_k,ilevel,istat_s3a)
        if(istat_s3a .ne. 1)then
            write(6,*)' No S3A data available'
        endif

        return
        end
