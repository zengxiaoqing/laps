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

        subroutine get_radar_max_pd(i4time_beg,i4time_end,imax,jmax,kmax
     1                             ,heights_3d,ext_radar
     1                             ,max_radar_files                      ! I
     1                             ,lat,lon,topo
     1                             ,radar_max,frac_sum,istatus)

!       1992         Steve Albers
!       1996 Feb     Steve Albers  Call read_radar_2dref for radar data

!       Input
        real lat(imax,jmax)
        real lon(imax,jmax)
        real topo(imax,jmax)
        real heights_3d(imax,jmax,kmax)

!       Output
        real radar_max(imax,jmax) ! M

        real dbz_2d(imax,jmax)

        character*9 asc_tim_9,asc_tim_9_beg,asc_tim_9_end
        integer i4time_file(MAX_RADAR_FILES)
        real frac(MAX_RADAR_FILES)
        character c_fnames(MAX_RADAR_FILES)*80

        character*255 c255_radar_filename
        character*255 c_filespec

        real grid_ra_ref(imax,jmax,kmax)

        character*3 ext_radar
        character*150  directory

        character*4 radar_name

        integer max_radar_gap
        parameter (max_radar_gap = 4800)

        data mode_radar/1/

        istatus = 1

c       rmax_so_far = 0.

        call make_fnam_lp(i4time_beg,asc_tim_9_beg,istatus)
        call make_fnam_lp(i4time_end,asc_tim_9_end,istatus)

        write(6,*)' Max Radar Reflectivity from ',asc_tim_9_beg,
     1                         ' to ',asc_tim_9_end

!       Get File Times
        call get_directory(ext_radar,directory,len_dir)
        c_filespec = directory(1:len_dir)//'*.'//ext_radar

        call Get_file_names(c_filespec,
     1                      i_nbr_files_ret,
     1                      c_fnames,
     1                      max_radar_files,
     1                      i_status )

        min_diff = 1999999999

        if(i_nbr_files_ret .gt. 0)then
            call get_directory_length(c_fnames(1),lenf)
        else ! Error Condition
            write(6,*)' No Radar Data Available'
            istatus = 0
            return
        endif

10      do i=1,i_nbr_files_ret
            asc_tim_9 = c_fnames(i)(lenf+1:lenf+9)
            call i4time_fname_lp(asc_tim_9,I4time_file(i),istatus)
        enddo

        call get_fracs_radar(i4time_beg,i4time_end,max_radar_gap
     1                      ,MAX_RADAR_FILES
     1                      ,i_nbr_files_ret
     1                      ,i4time_file,frac,frac_sum,istatus)

        if(istatus .ne. 1)then
            return
        endif

!       Initialize
        do j = 1,jmax
        do i = 1,imax
            radar_max(i,j) = 0.
        enddo ! i
        enddo ! j

        i4_interval = (i4time_end - i4time_beg)
        i4time_mid = i4time_beg + i4_interval / 2
        i4time_temp = 0

        do ifile = 1,i_nbr_files_ret
          if(frac(ifile) .gt. 0.)then

            i4time_radar = i4time_file(ifile)

            call make_fnam_lp(i4time_radar,asc_tim_9,istatus)

!           write(6,101)asc_tim_9,frac(ifile)
!101        format(' Time, Frac = ',a9,2x,f6.3)

            c255_radar_filename
     1      = c_fnames(ifile)(1:lenf)//asc_tim_9//'.'//ext_radar

            write(6,*)

            if(ext_radar .ne. 'vrc')then

              call get_ref_base(ref_base,istatus)
              if(istatus .ne. 1)return

              i4_tol = 1200

              call read_radar_3dref(i4time_radar,
     1                 i4_tol,i4_ret,                                   ! I/O
     1                 .true.,ref_base,
     1                 imax,jmax,kmax,ext_radar,
     1                 lat,lon,topo,.false.,.false.,
     1                 heights_3d,
     1                 grid_ra_ref,
     1                 rlat_radar,rlon_radar,rheight_radar,radar_name,     
     1                 n_ref_grids,istatus_2dref,istatus_3dref)       


              if(istatus_2dref .eq. 0)then
                write(6,*)' Error in reading radar data'
                istatus = 0
                return
              endif

              call get_max_ref(grid_ra_ref,imax,jmax,kmax,dbz_2d)

            else

              call read_radar_2dref(i4time_radar,radar_name,
     1                  imax,jmax,
     1                  dbz_2d,istatus_2dref)

              if(istatus_2dref .eq. 0)then
                write(6,*)' Error in reading radar data'
                istatus = 0
                return
              endif

            endif

            write(6,*)' Increment Radar Data'

            do j = 1,jmax
            do i = 1,imax
                radar_max(i,j) = max(radar_max(i,j),dbz_2d(i,j))
c               rmax_so_far = max(rmax_so_far,dbz_2d(i,j))
            enddo ! i
            enddo ! j

          endif ! Frac > 0

c         write(6,*)' Cycle to next file',rmax_so_far

        enddo ! i

        istatus = 1

        return
        end

        subroutine get_fracs_radar(i4time_beg,i4time_end,max_radar_gap
     1          ,MAX_RADAR_FILES
     1          ,i_nbr_files_ret,i4time_file,frac,frac_sum,istatus)

        integer MAX_RADAR_FILES

        character*9 asc_tim_9
        real frac(MAX_RADAR_FILES)
        integer i4time_file(MAX_RADAR_FILES)

        i4_interval = i4time_end - i4time_beg

        do i = 1,i_nbr_files_ret
!           write(6,301)i4time_beg,i4time_file(i),i4time_end
!301        format(1x,3i12)
            frac(i) = 0.
        enddo

        do i = 1,i_nbr_files_ret
            if(i4time_file(i) .ge. i4time_beg .and.
     1       i4time_file(i) .le. i4time_end        )then ! Within Pd
                frac(i) = 1.0

            endif
        enddo ! i

        frac_sum = 0.
        do ifile = 1,i_nbr_files_ret
          if(frac(ifile) .gt. 0.)then
            call make_fnam_lp(i4time_file(ifile),asc_tim_9,istat_fnam)
            write(6,101)asc_tim_9,frac(ifile)
 101        format('        FileTime, Frac = ',a9,2x,f8.5)
            frac_sum = frac_sum + frac(ifile)
          endif
        enddo

        if(frac_sum .eq. 0.)then
            write(6,*)' No radar data in the time window'
            istatus = 0
        endif

999     if(istatus .ne. 1)then
            write(6,*)' Insufficient files within time window'
            return
        endif

        return
        end
