

      subroutine map_aoml_sub(nilaps,njlaps,aoml_path_in,vrc_outdir    ! I
     1                       ,i4time_sys,laps_cycle_time               ! I
     1                       ,istatus)                                 ! O

      real lat(nilaps,njlaps)
      real lon(nilaps,njlaps)
      real topo(nilaps,njlaps)
      real rlaps_land_frac(nilaps,njlaps)

      real dbz(nilaps,njlaps)

      integer cvt_wfo_fname13_i4time

      character fname_in*140, path*140, basename*40
      character*13 wfo_fname13
      character*200 aoml_path_in
      character*7 vrc_outdir

      character*255 c_filespec_in,vrc_full_path
      integer max_files
      parameter(max_files = 3000)
      character*255 c_fnames_in(max_files),c_fnames_out(max_files)
!     integer i4times_in(max_files)
      integer i4times_out(max_files)

      logical l_match

      call get_r_missing_data(r_missing_data,istatus)

c read in laps lat/lon and topo
      call get_laps_domain_95(nilaps,njlaps,lat,lon,topo
     1                       ,rlaps_land_frac
     1                       ,grid_spacing_cen_m,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Error getting LAPS domain'
          istatus = 0
          return
      endif

      i4time_proc_end  = i4time_sys + laps_cycle_time/2 + 1800
      i4time_proc_strt = i4time_sys - laps_cycle_time/2 - 1800

!     Get times of output files
      i_vrc = 2
      call get_vrc_full_path(i_vrc,vrc_outdir
     1                        ,vrc_full_path,lenv,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Return from map_aoml_sub without AOML processing'
          return
      endif

      call get_file_times(vrc_full_path,max_files,c_fnames_out
     1                   ,i4times_out,i_nbr_files_out,istatus)


!     Get names and times of input files
      call s_len(aoml_path_in,lenp)
      if(lenp .eq. 0)then
          write(6,*)'map_aoml_sub: aoml_path_in has zero length'
          istatus = 0
          return
      endif

      c_filespec_in = aoml_path_in
      call get_file_names(c_filespec_in,i_nbr_files_in,c_fnames_in
     1                   ,max_files,istatus)

      if(i_nbr_files_in .eq. 0)then
          write(6,*)'No AOML files available'
          istatus = 0
          return
      endif

      do ifile_in = 1,i_nbr_files_in

          write(6,*)

!         basename = 'lf050705I1_stm_230058.swp'
!         call s_len(basename,lenb)
!         fname_in = path(1:lenp)//'/'//basename(1:lenb)

!         Split out basename from full file name
          fname_in = c_fnames_in(ifile_in)
          call s_len(fname_in,lenf)

          if(fname_in(lenp:lenp) .eq. '/')then
              lenf_start = lenp+1
          else
              lenf_start = lenp+2
          endif

          if(lenf .ge. lenf_start)then
              basename = fname_in(lenf_start:lenf)
          else
              write(6,*)' lenp/lenf mismatch ',lenp,lenf_start,lenf
              istatus = 0
              return
          endif

!         Obtain i4time corresponding to the basename, we handle the hours
!         and minutes separately since AOML can go to more than 24 hours
          wfo_fname13 = '20'//basename(3:8)//'_'//basename(16:19)
          write(6,*)'wfo_fname13 = ',wfo_fname13 

          wfo_fname13 = '20'//basename(3:8)//'_'//'0000'

          read(basename(16:17),*)ihour
          read(basename(18:19),*)imin

          i4time_file_in = cvt_wfo_fname13_i4time(wfo_fname13)

          i4time_file_in = i4time_file_in + ihour*3600 + imin*60

          write(6,*)' ihour/imin/i4time_file_in = '
     1               ,ihour,imin,i4time_file_in

          if(i4time_file_in .ge. i4time_proc_strt .AND.
     1       i4time_file_in .le. i4time_proc_end        )then 

!             Check for a match with output file times
              l_match = .false.
              do ifile_out = 1,i_nbr_files_out
                  if(i4time_file_in .eq. i4times_out(ifile_out))then
                    l_match = .true.
                  endif
              enddo ! ifile_out

              if(.not. l_match)then
                  write(6,*)' No output match - processing this time'
                  call map_aoml_sweep(nilaps,njlaps,lat,lon,dbz,fname_in
     1                           ,r_missing_data
     1                           ,rlat_radar,rlon_radar)

                  write(6,*)' calling put_vrc ',vrc_outdir

                  call put_vrc(i4time_file_in,comment_2d 
     1                        ,rlat_radar,rlon_radar,rheight_radar
     1                        ,lat,lon,topo
     1                        ,dbz,nilaps,njlaps,i_vrc
     1                        ,vrc_outdir,r_missing_data,istatus)

              else
                  write(6,*)' Output already exists - skip this time'

              endif

          else
              write(6,*)' Outside Time Window'

          endif

      enddo ! ifile_in

      istatus = 1

      return
      end

      subroutine map_aoml_sweep(nilaps,njlaps,lat,lon,dbz,name
     1                         ,r_missing_data
     1                         ,olat,olon)

      real lat(nilaps,njlaps)
      real lon(nilaps,njlaps)
      real dbz(nilaps,njlaps)

      real sum(nilaps,njlaps)
      integer npt(nilaps,njlaps)

      CHARACTER      keyword*4, fltname*8, stmname*12, radarid*4,
     +               trackname*32, createtime*32, ramname*28,
     +               name*140, jfile*140

      parameter (imax=240)
      parameter (jmax=240)
      parameter (max_x=imax)
      parameter (max_y=jmax)
      INTEGER*2    zarray(0:max_x-1,0:max_y-1) !ROW OF DBZ VALUES TO BE 
      INTEGER*2    kpac(32767)
      lunit = 1
      luout = 2

      write(6,*)' filename = ',name

      ierr = 0

      CALL READHEAD(luout, name, lunit, ierr, keyword, fltname,
     +  stmname, radarid, trackname, createtime, ramname, imax, jmax,
     +  kmax, nsweeps, nmosmflag, iunfoldflag, intattflag, ieditflag,
     +  iextra2, iextra3, stime, etime, olat, olon, sx, sy, sz, xdis,
     +  ydis, z0, rot, radaralt, calibco1, calibco2, azmcor, elcor,
     +  thresh, dbzrs, pcor, dcor, rcor, starthorelev, htbb, dbb)
C     OPENS RADAR DATA FILE AND READS RADAR HEADER DATA VALUES.
C     PAUL A. LEIGHTON, USDC/NOAA/AOML/HRD, 4 JUN 1991

      if(ierr .ne. 0)then
          write(6,*)' Note: ierr from AOML READHEAD = ',ierr
          return
      endif

      write(6,*)' stmname = ',stmname
      write(6,*)' olat/olon = ',olat,olon
      write(6,*)' imax,jmax = ',imax,jmax
      write(6,*)' xdis,ydis,rot = ',xdis,ydis,rot
      write(6,*)' sx,sy = ',sx,sy
      write(6,*)' radaralt = ',radaralt

      CALL PASSXCMP(lunit, imax, jmax, max_x, max_y, zarray)

!     write(10,*)zarray

      icen = imax/2
      jcen = jmax/2

      elev = 0.
      rheight_radar = 3000.

!     Initialize arrays
      sum = 0.
      npt = 0

!     Use pixel averaging onto laps grid
      do i = 0,imax-1
      do j = 0,jmax-1

        if(zarray(i,j) .gt. -999)then
          riradar = float(i-icen)
          rjradar = float(j-jcen)
          range = sqrt(riradar**2 + rjradar**2) * sx * 1000.
          azimuth = atan3d(riradar,rjradar) + rot
          call radar_to_latlon(rlat,rlon,height_grid
     1                  ,azimuth,range,elev
     1                  ,olat,olon,rheight_radar)

          call latlon_to_rlapsgrid(rlat,rlon,lat,lon,nilaps,njlaps
     1                            ,rilaps,rjlaps,istatus)     

          ilaps = nint(rilaps)
          jlaps = nint(rjlaps)

          if(ilaps .ge. 1 .and. ilaps .le. nilaps .and.
     1       jlaps .ge. 1 .and. jlaps .le. njlaps       )then
              sum(ilaps,jlaps) = sum(ilaps,jlaps) + zarray(i,j)
              npt(ilaps,jlaps) = npt(ilaps,jlaps) + 1
!             write(6,*)i,j,ilaps,jlaps,npt(ilaps,jlaps),zarray(i,j)
          endif
        endif

      enddo ! j
      enddo ! i

!     Divide arrays
      do il = 1,nilaps
      do jl = 1,njlaps
          if(npt(il,jl) .gt. 0)then
              dbz(il,jl) = sum(il,jl) / float(npt(il,jl))
          else
              dbz(il,jl) = r_missing_data
          endif   
      enddo ! jl
      enddo ! il

      return
      END
