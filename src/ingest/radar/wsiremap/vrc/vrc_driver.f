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
       program vrc_driver
       include 'lapsparms.cmn'
       
       Call get_laps_config('nest7grid',IStatus)
       if(istatus.eq.1)then
         write(*,*)'LAPS Parameters obtained'
       else
          write(*,*)'IStatus = ',IStatus,'Error - Get_LAPS_Config'
          write(*,*)'Terminating LAPS-VRC. WSI remapping'
          stop
       end if
       call vrc_driver_sub(NX_L_CMN,NY_L_CMN,c_raddat_type)
       stop
       end

       subroutine vrc_driver_sub(nx_l,ny_l,c_raddat_type)
c
c Program drives transformation of NOWRAD high density (hd) radar to LAPS
c domain (subroutine NOWRAD_to_LAPS). 'hd' files are assumed to reside in
c /public/data/radar/wsi/nowrad/netcdf.  The pathway to these data can be
c changed pretty easily.
c
ccc       include 'lapsparms.for'

       integer extreme_thrsh_70
       integer extreme_thrsh_47

       parameter (max_files=100,
     &            extreme_thrsh_47=0.30,
     &            extreme_thrsh_70=0.10)

       character*50 dir_vrc
       character*31 ext_vrc
       character*200 dir_static

       character*125 comment_ll(2),comment_vrc
       character*10 units_ll(2),units_vrc
       character*3 var_ll(2),var_vrc

       character*11 laps_dom_file
       character*4 lvl_coord_2d

       real*4 lat(nx_l,ny_l),lon(nx_l,ny_l),rdbz(nx_l,ny_L)
       real*4 grid_spacing
       real*4 data(nx_l,ny_l,2)

       real*4 percent_extreme_47
       real*4 percent_extreme_70

       integer i4time_cur
       integer i4time_latest_diff
       integer i4time_data
       integer i4time_latest_vrc
       integer i4time_latest_wsi
       integer i4time_now_gg
       integer i4_validTime
       integer istatus
       integer lvl_2d
       integer n,nn,nd, len
       integer n_vars_req

       character*100 c_values_req
       character*40  c_vars_req

       character*13 c_fname_cur
       character*13 c_fname_cur_temp
       character*13 cvt_i4time_wfo_fname13
       character*14 c_filetime
       character*9 wfo_fname13_to_fname9
       character*9 c_filename
       character*100 wsi_dir_path
       character*255 c_filespec, fname
       character*200 c_filenames_proc(max_files)
       character*3   c_raddat_type
     
       data lvl_2d/0/

c
c set filename. wfo data is 13 character. HOWEVER, if reading from /public
c then we must fool the wfo filename to still be 9 characters.
c
       i4time_cur = i4time_now_gg()
       c_fname_cur_temp = cvt_i4time_wfo_fname13(i4time_cur)

       if(c_raddat_type.eq.'wfo'.and.wsi_dir_path(2:7).ne.
     1'public')then
          c_fname_cur = c_fname_cur_temp
       else
          c_fname_cur = wfo_fname13_to_fname9(c_fname_cur_temp)
       endif
c
c get systime.
c
       write(6,*)'Current (nominal) file time: ',c_fname_cur
c
c build directory path and filename for expected nowrad data file.
c
       write(6,*)'Read Static info - path_to_wsi_2d_radar'
       n_vars_req = 1
       c_vars_req='path_to_wsi_2d_radar'
       call get_static_info(c_vars_req,c_values_req,n_vars_req
     1                                               ,istatus)
       if(istatus.eq.1)then

          write(6,*)'Successful getting path_to_wsi_2d_radar'
          write(6,*)'c_values_req = ',c_values_req
          wsi_dir_path=c_values_req

       else

          write(6,*)'Error getting path_to_wsi_2d_radar'
          write(6,*)'Terminating'
          goto 898

       endif

       n=index(wsi_dir_path,' ')
       write(6,*)'wsi_dir_path = ',wsi_dir_path(1:n-1)
c
c new code to accomodate the submission of this code by "at"
c
       if(c_raddat_type.eq.'wfo')then
          c_filespec=wsi_dir_path(1:n-1)//'*'
       else
          c_filespec=wsi_dir_path(1:n-1)//'*_hd'
       endif
       write(6,*)
       write(6,*)'Latest time for wsi data'

       call get_latest_file_time(c_filespec,i4time_latest_wsi)
c
c convert to fname9 and determine if this time has already been processed
c into vrc.
c
       call get_directory('vrc',dir_vrc,nd)
c       dir_vrc = '../lapsprd/vrc/'    !this also use below for output
c       nd = index(dir_vrc,' ')-1
       c_filespec = dir_vrc(1:nd)//'*'
       write(6,*)
       write(6,*)'Latest time for vrc files'

       call get_latest_file_time(c_filespec,i4time_latest_vrc)
c
c get wait-for-data parameters
c
       call get_directory('static',fname,len)
       
       open(15,
     &      file=fname(1:len)//'vrc/vrc_wait.parms',
     &      form='formatted',err=995)
       read(15,87)i4_check_interval
       read(15,87)i4_total_wait
       read(15,87)i4_thresh_age
       close(15)
87     format(i4)
c
       i4time_latest_diff = i4time_latest_vrc-i4time_latest_wsi

       if(i4time_latest_diff .eq. 0)then

          i4time_latest_vrc = i4time_latest_vrc + 900

          n=index(wsi_dir_path,' ')
          if(c_raddat_type.eq.'wfo')then
             c_filespec=wsi_dir_path(1:n-1)//'*'
          else
             c_filespec=wsi_dir_path(1:n-1)//'*_hd'
          endif
          n=index(c_filespec,' ')

          write(6,*)'Directory_wait: ',c_filespec(1:n)
          write(6,*)
          write(6,*)'Calling wait_for_data'
          write(6,*)'check-interval, total_wait, thresh age'
          write(6,*)i4_check_interval,i4_total_wait,i4_thresh_age

          call wait_for_data(c_filespec,i4time_latest_vrc
     1               ,i4_check_interval,i4_total_wait
     1               ,i4_thresh_age
     1               ,istatus)
          if(istatus .ne.1)then
             write(6,*)'wait for data did not find the data'
             goto 998
          endif

       elseif(i4time_latest_diff.gt.i4_thresh_age)then

          write(6,*)'Data too old to wait'
          goto 998

       endif
c
c since we found the data, set nfiles to 1 and set the filenames for input.
c in theory we could find more than one new file to process. In such a case
c nfiles > 1.
c
       nfiles=1
       n=index(wsi_dir_path,' ')
       if(c_raddat_type.eq.'wfo')then
          c_filespec=wsi_dir_path(1:n-1)//'*'
c          c_raddat_type=c_raddat_types(1)
       else
          c_filespec=wsi_dir_path(1:n-1)//'*_hd'
c          c_raddat_type=c_raddat_types(1)
       endif

       call get_latest_file_time(c_filespec,i4time_latest_wsi)

       if(c_raddat_type.eq.'wfo'.and.wsi_dir_path(2:7).ne.'public'
     1)then
          c_filetime=cvt_i4time_wfo_fname13(i4time_latest_wsi)
          nn=index(c_filetime,' ')-1
          c_filenames_proc(nfiles)=wsi_dir_path(1:n-1)//c_filetime(1:nn)
       else
          call make_fnam_lp (i4time_latest_wsi, c_filename, ISTATUS)
          c_filetime=c_filename
          nn=index(c_filetime,' ')-1
          c_filenames_proc(nfiles)=wsi_dir_path(1:n-1)//c_filetime(1:nn)
     1//'_hd'
       endif
c
c This for output.  LAPS VRC files as indicated.
c
       ext_vrc = 'VRC'
       var_vrc = 'REF'
       units_vrc = 'DBZ'
       comment_vrc = 'WSI Nowrad data remapped to LAPS domain'
c
c Definitions needed for acquiring LAPS latitude and longitude arrays.
c
c       dir_static='../static/'
       call get_directory('static',dir_static,len)
ccc       laps_dom_file=laps_domain_file
       var_ll(1)='LAT'
       var_ll(2)='LON'

       call rd_laps_static(dir_static, 'nest7grid', nx_l, ny_l, 2,
     &     var_ll, units_ll, comment_ll, data, grid_spacing,
     &     istatus)

       if(istatus.eq.1)then
          write(*,*)'LAPS lat/lon grid obtained'
          write(*,*)
          do j=1,ny_l
             do i=1,nx_l
                lat(i,j)=data(i,j,1)
                lon(i,j)=data(i,j,2)
             end do
          end do
       else
          write(*,*)'Unable to get lat/lon data'
          write(*,*)'NOWRAD-VRC process terminating'
          stop
       end if
c
c *****************************************************************************
c process the nowrad high density data. Remap to LAPS domain.  Generate output.
c *****************************************************************************
c
       do k=1,nfiles
       call NOWRAD_to_LAPS(c_raddat_type,
     &                 c_filenames_proc(k),
     &                 nx_l,ny_l,
     &                 lat,
     &                 lon,
     &                 i4_validTime,
     &                 rdbz,
     &                 istatus)
       if(istatus .eq. 1)then
          write(6,*)'WSI data properly remapped'
       else
          goto 18
       endif
c
c quick QC check
c
       n_extreme=0
       do j=1,ny_l
       do i=1,nx_l
          if(rdbz(i,j) .ge. 70.)then
             n_extreme=n_extreme+1
          end if
       end do
       end do
       percent_extreme_70 = n_extreme/(nx_l*ny_l)
       n_extreme=0
       do j=1,ny_l
       do i=1,nx_l
          if(rdbz(i,j) .ge. 47.)then
             n_extreme=n_extreme+1
          end if
       end do
       end do
       percent_extreme_47 = n_extreme/(nx_l*ny_l)


       n=index(c_filenames_proc(k),' ')

       if(percent_extreme_70 .le. extreme_thrsh_70 .and.
     &    percent_extreme_47 .le. extreme_thrsh_47)then
c
c   Output...i4time is i4 time 1960.
c
          i4time_data=i4_validTime+315619200

          call write_laps_data(i4time_data,
     &                         dir_vrc,
     &                         ext_vrc,
     &                         nx_l,ny_l,1,1,
     &                         var_vrc,
     &                         lvl_2d,
     &                         lvl_coord_2d,
     &                         units_vrc,
     &                         comment_vrc,
     &                         rdbz,
     &                         istatus)

          if(istatus.eq.1)then
             write(*,*)'VRC file successfully written'
             write(*,*)'for: ',c_filenames_proc(k)(1:n-1)
             write(*,*)'i4 time: ',i4time_data
          else
             goto 14
          end if

       else

          write(6,*)'Either'
          write(6,*)'More than 10% of dBZ >= 70 OR'
          write(6,*)'More than 25% of dBZ >= 47; thus,'
          write(6,*)'Not writing this vrc'
          write(*,*)'i4Time: ',i4time_data,
     &              ' fTime: ',c_filenames_proc(k)(1:n-1)

       end if

       goto 25

18     write(6,*)'Error with the wsi data. Not writing a'
       write(6,*)'vrc file for this time'
       write(6,*)'Filename: ',c_filenames_proc(k)(1:n-1)

       goto 25

14     write(*,*)'Error writing VRC file - Terminating'
       write(*,*)'i4Time: ',i4time_data,
     &           ' fTime: ',c_filenames_proc(k)(1:n-1)

25     enddo

990    write(6,*)'Normal completion of vrc_driver'
       goto 16

898    write(6,*)'Error getting raw data path'
       goto 16

995    write(6,*)'Error opening wait.parms'
       goto 16

999    write(*,*)'Error reading systime.dat - terminating'
       write(*,*)
       goto 16
998    write(6,*)'Stopping because no data was found'

16     stop
       end
