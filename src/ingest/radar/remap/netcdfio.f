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

 
       subroutine radar_init(i_radar,path_to_radar,path_to_vrc,itimes   ! I
     1                      ,b_missing_data                             ! I
     1                      ,i_tilt_proc                                ! I/O
     1                      ,i_last_scan,istatus)                       ! O
 
!      Open/Read Polar NetCDF file for the proper time
       integer max_files
       parameter(max_files=1000)

       character*150 path_to_radar,c_filespec,filename,directory
     1              ,c_fnames(max_files)
       character*15 path_to_vrc
       character*9 a9_time
       integer*4 i4times_raw(max_files),i4times_lapsprd(max_files)
       character*2 c2_tilt
       character*3 laps_radar_ext, c3_radar_subdir
       character*8 radar_subdir, c8_fname_format
       logical l_multi_tilt,l_exist,l_output,l_realtime
       character*13 a9_to_rsa13

       save a9_time, i_nbr_files_raw, i_nbr_files_2nd
       save i4times_raw

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
       include 'remap_constants.dat' ! for debugging only? (also for structure)
!      include 'remap.cmn' ! for debugging only

!      This call is still needed for return of 'laps_radar_ext/c3_radar_subdir'
!      We could change this to pass these in through the 'radar_init' call
       call get_remap_parms(i_radar,n_radars_remap                        ! I/O
     1                    ,max_times,path_to_radar                        ! O
     1                    ,laps_radar_ext,c3_radar_subdir                 ! O
     1                    ,path_to_vrc                                    ! O
     1                    ,ref_min,min_ref_samples,min_vel_samples,dgr    ! O
     1                    ,namelist_parms,istatus)                        ! O
       if(istatus .ne. 1)then
           write(6,*)'Warning: bad status return from get_remap_parms'       
           return
       endif

       call get_systime_i4(i4time_sys,istatus)
       if(istatus .ne. 1)then ! use wall clock time if systime is not available
           i4time_sys = i4time_now_gg()
           l_realtime = .true.

       else ! success reading 'systime.dat'
           i4_age = i4time_now_gg() - i4time_sys
           if(i4_age .gt. 43200)then ! archive case assumed
               l_realtime = .false.
           else
               l_realtime = .true.
           endif
       endif

       call get_laps_cycle_time(laps_cycle_time,istatus)   
       if(istatus .ne. 1)return   

       c8_fname_format = 'UNKNOWN'
c
c      Determine filename extension
       write(6,*)' radar_init: laps_ext = ',laps_radar_ext
       if(laps_radar_ext(1:1) .ne. 'v')then ! Sanity check
           laps_radar_ext = 'v01'
       endif

       if(laps_radar_ext .eq. 'vrc')then 
           radar_subdir = c3_radar_subdir
           write(6,*)' radar_init: radar_subdir = ',radar_subdir
       endif

       i_last_scan = 0

       if(i_tilt_proc .lt. 10)then
           write(c2_tilt,101)i_tilt_proc
 101       format('0',i1)
       else
           write(c2_tilt,102)i_tilt_proc
 102       format(i2)
       endif

       if(i_tilt_proc .eq. 1)then
           call s_len(path_to_radar,len_path)
 
!          Get filecount of 02 elevation raw files
           c2_tilt = '02'
           c_filespec = path_to_radar(1:len_path)//'/*elev'//c2_tilt       

           if(itimes .eq. 1)then
               call get_file_names(c_filespec,i_nbr_files_2nd,c_fnames
     1                            ,max_files,istatus)
               if(istatus .ne. 1)then
                   return
               endif
           endif

           write(6,*)' # of 2nd tilt raw files = ',i_nbr_files_2nd

           I4_elapsed = ishow_timer()

           if(i_nbr_files_2nd .gt. 0)then
               l_multi_tilt = .true.
               write(6,*)' We have multiple tilt data'
           else
               l_multi_tilt = .false.
               write(6,*)' We have single tilt data'
           endif
 
!          Get i4times of 01 elevation raw files
           c2_tilt = '01'
           c_filespec = path_to_radar(1:len_path)//'/*elev'//c2_tilt       

           if(itimes .eq. 1)then
               call get_file_times(c_filespec,max_files,c_fnames
     1                            ,i4times_raw,i_nbr_files_raw,istatus)
           endif

           write(6,*)' # of 1st tilt raw files = ',i_nbr_files_raw

           I4_elapsed = ishow_timer()

           if(istatus .ne. 1 .or. i_nbr_files_raw .eq. 0)then
               istatus = 0
               return
           endif

!          Get output filespec
           if(laps_radar_ext .ne. 'vrc')then
               call get_filespec(laps_radar_ext,1,c_filespec,istatus)

           else ! laps_radar_ext = 'vrc', now check path_to_vrc
               if(path_to_vrc .eq. 'rdr')then
                   call get_directory('rdr',directory,len_dir)
                   c_filespec = directory(1:len_dir)//radar_subdir(1:3)
     1                          //'/vrc'     
               else ! path_to_vrc = 'lapsprd'
                   call get_filespec(laps_radar_ext,1,c_filespec
     1                              ,istatus)      
               endif

           endif

           call s_len(c_filespec,lenspec)
           write(6,*)' Output filespec = ',c_filespec(1:lenspec)

!          Get i4times of output files 
           call get_file_times(c_filespec,max_files,c_fnames
     1                   ,i4times_lapsprd,i_nbr_lapsprd_files,istatus)
           write(6,*)' # of output files = ',i_nbr_lapsprd_files
           
           i4time_process = 0

           if(l_realtime)then
               needed_raw_files = 2
               latest_raw_file = i_nbr_files_raw-1
           else
               needed_raw_files = 1
               latest_raw_file = i_nbr_files_raw
           endif

           write(6,*)' l_realtime/needed/latest'
     1                ,l_realtime,needed_raw_files,latest_raw_file

!          Get input filetime to process
           if(i_nbr_files_raw .ge. needed_raw_files)then ! use earliest time
               i4_earliest_window = i4time_sys - laps_cycle_time - 1800       
               call make_fnam_lp(i4_earliest_window,a9_time,istatus)

               write(6,*)
     1           ' Looking for earliest unprocessed input file back to '       
     1           ,a9_time

               do i = latest_raw_file,1,-1
                   if(i .eq. latest_raw_file)then ! Write latest raw filetime
                       call make_fnam_lp(i4times_raw(i),a9_time,istatus)
                       write(6,*)' Latest raw filetime = ',a9_time
                   endif

                   l_output = .false.
                   do j = 1,i_nbr_lapsprd_files
                       if(i4times_raw(i) .eq. i4times_lapsprd(j))then
                           l_output = .true.
                       endif
                   enddo ! j

                   if( (.not. l_output)                      .AND. 
     1                i4times_raw(i) .gt. i4_earliest_window       )then
                       i4time_process = i4times_raw(i)
                   endif

               enddo

           endif

           if(i4time_process .ne. 0)then
               call make_fnam_lp(i4time_process,a9_time,istatus)
               write(6,*)' Processing filetime ',a9_time
           else
               write(6,*)' No new filetimes to process'
               istatus = 0
               return
           endif

       endif ! i_tilt_proc = 1

!      Pull in housekeeping data from 1st tilt
       write(6,*)' radar_init: looking for file... '

       i_skip = 0

!      Test existence of raw 'yyjjjhhmm_elevtt / yyyymmdd_hhmm.elevtt' input
 200   call check_input_file(path_to_radar,a9_time,i_tilt_proc      ! I
     1                      ,c8_fname_format                        ! I
     1                      ,filename,l_exist)                      ! O     

       if(l_exist)then
           call get_tilt_netcdf_data(filename
     1                               ,radarName
     1                               ,siteLat                        
     1                               ,siteLon                        
     1                               ,siteAlt                        
     1                               ,elevationAngle
     1                               ,numRadials
     1                               ,elevationNumber
     1                               ,VCP
     1                               ,r_nyquist
     1                               ,radialAzim
     1                               ,Z
     1                               ,V
     1                               ,resolutionV
     1                               ,gateSizeV,gateSizeZ
     1                               ,firstGateRangeV,firstGateRangeZ
     1                               ,MAX_VEL_GATES, MAX_REF_GATES
     1                               ,MAX_RAY_TILT
     1                               ,istatus)

           if(istatus .eq. 1 .and. namelist_parms%l_line_ref_qc)then
               write(6,*)' radar_init: calling rayqckz...'
               iscale = 2
               miss = nint(b_missing_data)
               call rayqckz(MAX_REF_GATES,numRadials,iscale
     1                     ,miss,Z,radialAzim)
           endif

       elseif(i_tilt_proc .le. 20)then
           i_tilt_proc = i_tilt_proc + 1
           i_skip = 1
           goto 200

       else
           istatus = 0

       endif ! file exists

       if(istatus .eq. 1
     1             .AND.
     1     .not. (laps_radar_ext .eq. 'vrc' .and. i_tilt_proc .gt. 1)
     1                                                             )then       
           if(i_tilt_proc .eq. 1)then
               write(6,201)elevationNumber, i_tilt_proc
 201           format(' elevationNumber, i_tilt_proc',2i4)

           else
               write(6,202)elevationNumber, i_tilt_proc
 202           format(' elevationNumber, i_tilt_proc',2i4
     1               ,' (upcoming tilt)')

               if(i_skip .eq. 1)then
                   write(6,*)
     1               ' WARNING: We had to skip past some missing tilts'       
               endif

           endif

       else
           write(6,*)' Could not read tilt # ',i_tilt_proc
           i_last_scan = 1

       endif

       write(6,*)

       istatus = 1
       return
       end
 
       function get_altitude()
       integer get_altitude          ! Site altitude (meters)

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
 
       get_altitude = nint(siteAlt)

       return
       end
 
 
       function get_latitude()
       integer get_latitude          ! Site latitude (degrees * 100000)

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
 
       get_latitude = nint(siteLat*100000)
       return
       end
 
 
       function get_longitude()
       integer get_longitude         ! Site longitude (degrees * 100000)

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'

       character*8 c8_project

!      call get_c8_project(c8_project,istatus)
!      if(istatus .ne. 1)then
!          write(6,*)' Error, no c8_project'
!          stop
!      endif

!      if(c8_project(1:3) .ne. 'CWB')then

       if(.false.)then
           get_longitude = -abs(nint(siteLon*100000))
       else
           get_longitude =      nint(siteLon*100000)
       endif

       return
       end

       subroutine get_radarname(c4_radarname,istatus)

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
       character*4 c4_radarname

       c4_radarname = radarName
       call upcase(c4_radarname,c4_radarname)
       write(6,*)' c4_radarname = ',c4_radarname

       istatus = 1
       return
       end
       
 
       function get_field_num(c3_field)
       integer get_field_num
       character*3 c3_field
 
       if(c3_field .eq. 'DBZ')get_field_num = 1
       if(c3_field .eq. 'VEL')get_field_num = 2
 
       return
       end
 
 
       function read_radial()
 
       read_radial = 0
       return
       end
 
 
       function get_status()
       integer get_status
 
       get_status = 0
       return
       end
 
 
       function get_fixed_angle()
       integer get_fixed_angle     ! Beam tilt angle (degrees * 100)

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'

       if( abs(elevationAngle) .le. 1e10 )then
           get_fixed_angle = nint(elevationAngle * 100.)
       else
           write(6,*)' warning in get_fixed_angle, invalid value'
           get_fixed_angle = -999 ! i_missing_data
       endif

 
       return
       end
 
 
       function get_scan()
       integer get_scan            ! Scan #

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
 
       get_scan = elevationNumber
       return
       end
 
 
       function get_tilt()
       integer get_tilt            ! Tilt #

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
 
       get_tilt = elevationNumber
       return
       end
 
 
       function get_num_rays()
       integer get_num_rays

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
 
       get_num_rays = numRadials
       return
       end
 
 
       subroutine get_volume_time(i4time_process_ret)

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
 
       i4time_process_ret = i4time_process
       return
       end
 
 
       function get_vcp()
       integer get_vcp

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
 
       get_vcp = VCP
       return
       end
 
 
       function get_azi(iray) ! azimuth * 100.
       integer get_azi

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
 
       if( abs(radialAzim(iray)) .le. 1e10 )then
           get_azi = nint(radialAzim(iray)*100.)
       else
           write(6,*)' warning in get_azi, azimuth = ',iray
     1                                   , radialAzim(iray)
           get_azi = -999 ! i_missing_data
       endif

       return
       end
 
 
       function get_nyquist()
       real get_nyquist                 ! Nyquist velocity of the radial (M/S)

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
 
       get_nyquist = r_nyquist
       return
       end
 
 
       function get_number_of_gates(index)
       integer get_number_of_gates
 
       get_number_of_gates = 0
       return
       end
 
 
       subroutine get_first_gate(index,first_gate_m,gate_spacing_m)

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'

       call get_r_missing_data(r_missing_data,istatus)
       if(istatus .ne. 1)then
           write(6,*)' Error calling get_r_missing_data'
           stop
       endif
 
       if(index .eq. 1)then
           if(firstGateRangeZ .ge. -10. .and. 
     1        firstGateRangeZ .lt. 1000.)then     
               first_gate_m = firstGateRangeZ * 1000.
           else
               write(6,*)' Warning: firstGateRangeZ is outside range'
     1                  ,firstGateRangeZ
               first_gate_m = r_missing_data
           endif

           if(gateSizeZ .ge. 0. .and. gateSizeZ .lt. 1000.)then     
               gate_spacing_m = gateSizeZ * 1000.
           else
               write(6,*)' Warning: gateSizeZ is outside range'
     1                  ,gateSizeZ
               gate_spacing_m = r_missing_data
           endif

       elseif(index .eq. 2)then
           if(firstGateRangeV .ge. -10. .and. 
     1        firstGateRangeV .lt. 1000.)then     
               first_gate_m = firstGateRangeV * 1000.
           else
               write(6,*)' Warning: firstGateRangeV is outside range'
     1                  ,firstGateRangeV
               first_gate_m = r_missing_data
           endif

           if(gateSizeV .ge. 0. .and. gateSizeV .lt. 1000.)then     
               gate_spacing_m = gateSizeV * 1000.
           else
               write(6,*)' Warning: gateSizeV is outside range'
     1                  ,gateSizeV      
               gate_spacing_m = r_missing_data
           endif

       endif

       return
       end
 
 
       function get_data_field(index, data, n_ptr, n_gates
     1                                           , b_missing_data)
       integer get_data_field

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'

!      Note that Z and V are integer arrays that have been read in from the
!      NetCDF files as byte values. These are here converted to dBZ and 
!      meters per second, respectively. This is done for one radial for
!      each call to this routine.
 
       real*4 data(n_gates)

       if(index .eq. 1)then ! reflectivity
           do i = 1,n_gates
               data(i) = Z(n_ptr + (i-1))

!              Convert from signed to unsigned
               if(data(i) .gt. 127.) then
                   print *, 'error in Reflectivity: ',i,data(i)
                   stop
               endif
               if(data(i) .lt. 0.) then
                   data(i) = 256. + data(i)
               endif

               if(data(i) .ne. b_missing_data)then ! Scale
                   data(i) = (data(i) - 2.)/2.0 - 32.
               endif

           enddo

       elseif(index .eq. 2)then ! velocity
           do i = 1,n_gates
               data(i) = V(n_ptr + (i-1))

!              Convert from signed to unsigned
               if(data(i) .gt. 127.) then
                   print *, 'error in Velocity: ',i,data(i)
                   stop
               endif
               if(data(i) .lt. 0.) then
                   data(i) = 256. + data(i)
               endif

               if(data(i) .eq. 1. .or. data(i) .eq. 0.)then 
                   data(i) = b_missing_data  ! Invalid Measurement
               endif

               if(resolutionV .eq. 0.)then ! QC Check
                   data(i) = b_missing_data
               endif

               if(data(i) .ne. b_missing_data)then ! Scale valid V
                   data(i) = (data(i) - 129.) * resolutionV
               endif

           enddo

       endif

       get_data_field = 1
       return
       end
 
 
       function cvt_fname_data()
 
       cvt_fname_data = 0
       return
       end
 
 
       subroutine check_input_file(path_to_radar,a9_time,i_tilt           ! I
     1                            ,c8_fname_format                        ! I/O
     1                            ,filename,l_exist)                      ! O

       logical l_exist

       character*(*) path_to_radar
       character*(*) filename
       character*2 c2_tilt
       character*9 a9_time
       character*8 c8_fname_format
       character*13 a13_time, fname9_to_wfo_fname13

!      Determine type of filename (if not yet known) and test for existence
!      of radar file for given time/tilt

       if(i_tilt .lt. 10)then
           write(c2_tilt,101)i_tilt
 101       format('0',i1)
       else
           write(c2_tilt,102)i_tilt
 102       format(i2)
       endif

       call s_len(path_to_radar,len_path)

       if(c8_fname_format .eq. 'UNKNOWN')then
           write(6,*)'check_input_file: resolving filename format'
       endif

       if(c8_fname_format .eq. 'NIMBUS' .or. 
     1    c8_fname_format .eq. 'UNKNOWN')then
           filename = path_to_radar(1:len_path)//'/'//a9_time//'_elev'
     1                //c2_tilt

!          Test existence of radar file
           inquire(file=filename,exist=l_exist)

           if(l_exist)c8_fname_format = 'NIMBUS'
       endif

       if(c8_fname_format .eq. 'WFO' .or. 
     1    c8_fname_format .eq. 'UNKNOWN')then
           a13_time = fname9_to_wfo_fname13(a9_time)
           filename = path_to_radar(1:len_path)//'/'//a13_time//'.elev'
     1                //c2_tilt

!          Test existence of radar file
           inquire(file=filename,exist=l_exist)

           if(l_exist)c8_fname_format = 'WFO'
       endif


       call s_len(filename,len_file)
       write(6,*)' check_input_file: ',filename(1:len_file),' ',l_exist       

       return
       end
