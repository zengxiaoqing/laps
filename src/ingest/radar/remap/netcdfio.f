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

 
       subroutine radar_init(i_radar,path_to_radar,path_to_vrc,itimes  ! I
     1                      ,b_missing_data                            ! I
     1                      ,i_tilt_proc                               ! I/O
     1                      ,l_realtime                                ! O
     1                      ,i_last_scan,istatus)                      ! O

       use mem_vol
 
!      Open/Read Polar NetCDF file for the proper time
       integer max_files
       parameter(max_files=20000)  ! max_radar_files

       character*150 path_to_radar,c_filespec,filename,directory
     1              ,c_fnames(max_files),c_fnames_in(max_files)
       character*15 path_to_vrc
       character*9 a9_time
       integer i4times_raw(max_files),i4times_lapsprd(max_files)
       character*2 c2_tilt
       character*4 laps_radar_ext, c3_radar_subdir
       character*8 radar_subdir                             
       logical l_multi_tilt,l_exist,l_output,l_realtime
       character*13 a13_time
       integer cvt_wfo_fname13_i4time
       integer Z_bin, V_bin, radial

       character*31 station

!      integer gateR, gateR_HI, gateV, gateV_HI, radialR, radialR_HI,
!    +     radialV, radialV_HI, scanR, scanR_HI, scanV,
!    +     scanV_HI,nf_fid, nf_vid, nf_status

       save a9_time, i_nbr_files_raw, i_nbr_files_2nd, i_nbr_files_vol
       save i4times_raw, i4times_lapsprd, i_nbr_lapsprd_files

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
       include 'remap_constants.dat' ! for debugging only? (also for structure)
!      include 'remap.cmn' ! for debugging only

!      Make this allocatable?
!      integer   MAX_TILTS   
!      parameter (MAX_TILTS=10) 
!      integer   Z_vol(MAX_REF_GATES*MAX_RAY_TILT*MAX_TILTS)
!      integer   V_vol(MAX_VEL_GATES*MAX_RAY_TILT*MAX_TILTS)

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

       call get_r_missing_data(r_missing_data,istatus)
       if(istatus .ne. 1)then
           write(6,*)' Error calling get_r_missing_data'
           return
       endif

!      c8_fname_format = 'UNKNOWN'
c
c      Determine output filename extension
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
           if(len_path .le. 0)then
               write(6,*)' len_path <= 0, return from radar_init'
               istatus = 0
               return
           endif
 
!          Get filecount of 02 elevation raw files
           c2_tilt = '02'
           c_filespec = path_to_radar(1:len_path)//'/*elev'//c2_tilt       

           write(6,*)' itimes = ',itimes

           if(itimes .eq. 1)then
               i_nbr_files_vol = 0

               call get_file_names(c_filespec,i_nbr_files_2nd,c_fnames
     1                            ,max_files,istatus)
               if(istatus .ne. 1)then
                   write(6,*)' istatus returned from get_file_names ='        
     1                      ,istatus
                   write(6,*)' assume no directory present (tilt)'
                   return
               endif

               if(i_nbr_files_2nd .eq. 0)then ! try for volume files
                   c_filespec = path_to_radar(1:len_path)//'/*.nc'
                   write(6,*)' Volume filespec = ',trim(c_filespec)
                   call get_file_names(c_filespec,i_nbr_files_vol
     1                                ,c_fnames_in,max_files,istatus)
                   if(istatus .ne. 1)then
                       write(6,*)
     1                       ' istatus returned from get_file_names ='  
     1                        ,istatus
                       write(6,*)' assume no directory present (vol)'
                       return
                   else
                       write(6,*)' i_nbr_files_vol = ',i_nbr_files_vol
                   endif
               endif

           endif

           call s_len(c_filespec,lenspec)
           write(6,*)' Input filespec = ',c_filespec(1:lenspec)

           write(6,*)' # of 2nd tilt raw files = ',i_nbr_files_2nd

           I4_elapsed = ishow_timer()

           if(i_nbr_files_vol .gt. 0)then
               write(6,*)' We have volume data, # files ='
     1                                        ,i_nbr_files_vol
               c8_fname_format = 'VOLUME'
               c_filespec = path_to_radar(1:len_path)//'/*.nc'

               if(itimes .eq. 1)then ! determine file times
                   do i = 1,i_nbr_files_vol
                       if(i .eq. 1)then
                           call get_directory_length(c_fnames_in(i)
     1                                              ,lend)
                           call s_len(c_fnames_in(i),len_fname)
                           lenf = len_fname - lend
                       endif
                       a13_time = c_fnames_in(i)(lend+5:lend+17)
                       write(6,*)' vol a13_time = ',a13_time
                       i4times_raw(i) = cvt_wfo_fname13_i4time(a13_time)
                   enddo ! i
               endif

               i_nbr_files_raw = i_nbr_files_vol

               write(6,*)' # of (volume) raw files = ',i_nbr_files_raw

           else
               if(i_nbr_files_2nd .gt. 0)then
                   l_multi_tilt = .true.
                   write(6,*)' We have multiple tilt data'
               else
                   l_multi_tilt = .false.
                   write(6,*)' We have single tilt data'
               endif
 
!              Get i4times of 01 elevation raw files
               c2_tilt = '01'
               c_filespec = path_to_radar(1:len_path)//'/*elev'//c2_tilt       

               if(itimes .eq. 1)then
                   call get_file_times(c_filespec,max_files,c_fnames_in
     1                                ,i4times_raw,i_nbr_files_raw
     1                                ,istatus)
               endif

               write(6,*)' # of (1st tilt) raw files = ',i_nbr_files_raw

           endif

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
           if(itimes .eq. 1)then
               call get_file_times(c_filespec,max_files,c_fnames
     1                            ,i4times_lapsprd,i_nbr_lapsprd_files
     1                            ,istatus)

           else ! more efficiency for subsequent radars
               i_nbr_lapsprd_files = i_nbr_lapsprd_files + 1
               i4times_lapsprd(i_nbr_lapsprd_files) = i4time_process
               if(i4time_process .eq. 0)then
                   write(6,*)' WARNING: i4time_process = 0'
               else
                   write(6,*)' adding to i4times_lapsprd:'
     1                      ,i4times_lapsprd(i_nbr_lapsprd_files)
               endif

           endif
               
           write(6,*)' # of output files = ',i_nbr_lapsprd_files

           if(i_nbr_lapsprd_files .ge. 1)then ! Write latest output filetime
               call make_fnam_lp(i4times_lapsprd(i_nbr_lapsprd_files)
     1                          ,a9_time,istatus)
               write(6,*)' Latest output filetime = ',a9_time
           endif
           
           I4_elapsed = ishow_timer()

           i4time_process = 0

           if(l_realtime .and. i_nbr_files_vol .eq. 0)then
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

               if(l_realtime)then
                   i4_latest_window = i4time_now_gg()
               else
                   i4_latest_window = i4time_sys + laps_cycle_time + 1800
               endif
               call make_fnam_lp(i4_latest_window,a9_time,istatus)

               write(6,*)' Files accepted up to ',a9_time

               i_process = 0
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
     1                i4times_raw(i) .ge. i4_earliest_window       
     1                              .AND.
     1                i4times_raw(i) .le. i4_latest_window
     1                                                           )then
                       i4time_process = i4times_raw(i)
                       i_process = i
                   endif

               enddo

           endif

           if(i4time_process .ne. 0)then
               call make_fnam_lp(i4time_process,a9_time,istatus)
               write(6,*)' Processing file/a9time ',i_process,a9_time
           else
               write(6,*)' No new filetimes to process'
               istatus = 0
               return
           endif

       endif ! i_tilt_proc = 1

!      Pull in housekeeping data from 1st tilt
       write(6,*)' radar_init: looking for file for tilt... '
     1          ,i_tilt_proc

       i_skip = 0

!      Test existence of raw 'yyjjjhhmm_elevtt / yyyymmdd_hhmm.elevtt' input
 200   if(i_nbr_files_vol .eq. 0)then ! tilt data
           c8_fname_format = 'UNKNOWN'
           call check_input_file(path_to_radar,a9_time,i_tilt_proc      ! I
     1                          ,c8_fname_format                        ! I/O
     1                          ,filename,l_exist)                      ! O     
       else                           ! volume data
         l_exist = .true.
         if(i_tilt_proc .eq. 1)then   ! first tilt
!          filename = path_to_radar(1:len_path)//'/'//a9_time//'.nc'
           write(6,*)' processing volume directory file # ',i_process
           if(i_process .lt. 0 .or. i_process .gt. max_files)then
               write(6,*)' ERROR: i_process out of bounds in radar_init'
               stop
           endif
           filename = c_fnames_in(i_process)
           write(6,*)' c_fnames array',(c_fnames_in(i),i=1,i_process)
         endif ! first tilt
       endif

       if(l_exist)then ! these calls will fill the variables in 
                       ! 'netcdfio_radar_common.inc'

         write(6,*)' c8_fname_format = ',trim(c8_fname_format)

         if(c8_fname_format .ne. 'VOLUME')then

           write(6,*)filename
           call get_tilt_netcdf_hdr  (filename,nf_fid
     1                               ,radarName  ! Is this returned?
     1                               ,siteLat                        
     1                               ,siteLon                        
     1                               ,siteAlt                        
     1                               ,elevationAngle
     1                               ,numRadials 
     1                               ,elevationNumber
     1                               ,VCP
     1                               ,radialAzim
     1                               ,resolutionV
     1                               ,gateSizeV,gateSizeZ
     1                               ,firstGateRangeV,firstGateRangeZ
     1                               ,MAX_VEL_GATES, MAX_REF_GATES ! I
     1                               ,MAX_RAY_TILT                 ! I
     1                               ,V_bin,     Z_bin,     radial ! O
     1                               ,istatus)

!          Equivalent names, particularly in 'netcdfio_radar_common.inc'
           numRadials = radial
           ngates_ref_cdf = Z_bin
           ngates_vel_cdf = V_bin

!          Initialize in case they aren't present in the NetCDF file
           Z_scale  = r_missing_data
           Z_offset = r_missing_data
           V_scale  = r_missing_data
           V_offset = r_missing_data

!          Note that file remains open from call to 'get_tilt_netcdf_hdr'
           call get_tilt_netcdf_data(filename,nf_fid
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
     1                               ,Z_scale, Z_offset
     1                               ,V_scale, V_offset
     1                               ,V_bin,     Z_bin,     radial    ! I
     1                               ,istatus)

!          Use Default values if needed
           if(abs(Z_scale ) .gt. 1e10)Z_scale  = 2.0
           if(abs(Z_offset) .gt. 1e10)Z_offset = 66.
           if(abs(V_scale ) .gt. 1e10)V_scale  = 2.0
           if(abs(V_offset) .gt. 1e10)V_offset = 129.

           write(6,*)' Z/V scale/offset = ',Z_scale,Z_offset
     1                                     ,V_scale,V_offset

         else ! c8_fname_format = 'VOLUME'

           if(i_tilt_proc .eq. 1)then
               write(6,*)' call get_vol_netcdf_hdr'
               write(6,*)trim(filename)

               call get_vol_netcdf_hdr(filename,
     +            gateR, gateR_HI, gateV, gateV_HI, radialR, radialR_HI,
     +            radialV, radialV_HI, scanR, scanR_HI, scanV,
     +            scanV_HI,nf_fid, nf_vid, nf_status)

               write(6,*)' call get_attribute_vol'  

               call get_attribute_vol(nf_fid,siteLat,siteLon,siteAlt
     +            ,station
     +            ,istatus)

               radarName = station

               write(6,*)' lat/lon/alt/name ',siteLat,siteLon,siteAlt
     1                                       ,trim(radarName)

!              Deallocate old volume scan (if needed)
               if(allocated(Reflectivity))deallocate(Reflectivity)
               if(allocated(Reflectivity_HI))deallocate(Reflectivity_HI)
               if(allocated(RadialVelocity))deallocate(RadialVelocity)
               if(allocated(RadialVelocity_HI)) 
     1                                     deallocate(RadialVelocity_HI)

               if(allocated(elevationR))deallocate(elevationR)
               if(allocated(elevationR_HI))deallocate(elevationR_HI)
               if(allocated(elevationV))deallocate(elevationV)
               if(allocated(elevationV_HI))deallocate(elevationV_HI)

               if(allocated(azimuthR))deallocate(azimuthR)
               if(allocated(azimuthR_HI))deallocate(azimuthR_HI)
               if(allocated(azimuthV))deallocate(azimuthV)
               if(allocated(azimuthV_HI))deallocate(azimuthV_HI)

               if(allocated(distanceR))deallocate(distanceR)
               if(allocated(distanceR_HI))deallocate(distanceR_HI)
               if(allocated(distanceV))deallocate(distanceV)
               if(allocated(distanceV_HI))deallocate(distanceV_HI)
               if(allocated(nyquistVelocityV))
     1           deallocate(nyquistVelocityV)
               if(allocated(nyquistVelocityV_HI))
     1           deallocate(nyquistVelocityV_HI)

!              Allocate new volume scan
               write(6,*)'gateR,radialR,scanR = ',gateR,radialR,scanR
               allocate(Reflectivity(gateR,radialR,scanR))
               allocate(elevationR(radialR,scanR))
               allocate(azimuthR(radialR,scanR))
               allocate(distanceR(gateR))

               write(6,*)'gateR_HI,radialR_HI,scanR_HI = ',
     1                    gateR_HI,radialR_HI,scanR_HI
               allocate(Reflectivity_HI(gateR_HI,radialR_HI,scanR_HI))
               allocate(elevationR_HI(radialR_HI,scanR_HI))
               allocate(azimuthR_HI(radialR_HI,scanR_HI))
               allocate(distanceR_HI(gateR_HI))

               write(6,*)'gateV,radialV,scanV = ',gateV,radialV,scanV
               allocate(RadialVelocity(gateV,radialV,scanV))
               allocate(elevationV(radialV,scanV))
               allocate(azimuthV(radialV,scanV))
               allocate(distanceV(gateV))              
               allocate(nyquistVelocityV(gateV))              

               write(6,*)'gateV_HI,radialV_HI,scanV_HI = ',
     1                    gateV_HI,radialV_HI,scanV_HI
               allocate(RadialVelocity_HI(gateV_HI,radialV_HI,scanV_HI))
               allocate(elevationV_HI(radialV_HI,scanV_HI))
               allocate(azimuthV_HI(radialV_HI,scanV_HI))
               allocate(distanceV_HI(gateV_HI))                   
               allocate(nyquistVelocityV_HI(gateV_HI))                   

               write(6,*)' call get_vol_netcdf_data'
               call get_vol_netcdf_data(nf_fid, gateR, gateR_HI, 
     +              gateV, gateV_HI,
     +              radialR, radialR_HI, radialV, radialV_HI, 
     +              scanR, scanR_HI,
     +              scanV, scanV_HI, 
     +              Reflectivity, Reflectivity_HI,
     +              RadialVelocity, RadialVelocity_HI,
     +              elevationR, elevationR_HI,
     +              elevationV, elevationV_HI,
     +              azimuthR, azimuthR_HI,
     +              azimuthV, azimuthV_HI, 
     +              distanceR, distanceR_HI,
     +              distanceV, distanceV_HI,
     +              nyquistVelocityV, nyquistVelocityV_HI)

!              Use Default values
               Z_scale  = 2.0
               Z_offset = 66.
               V_scale  = 2.0
               V_offset = 129.
               resolutionV = 1. / V_scale

               write(6,*)' Z/V scale/offset = ',Z_scale,Z_offset
     1                                         ,V_scale,V_offset

               write(6,*)' elevationR = '   ,elevationR(1,:)
               write(6,*)' elevationR_HI = ',elevationR_HI(1,:)
               write(6,*)' elevationV = '   ,elevationV(1,:)
               write(6,*)' elevationV_HI = ',elevationV_HI(1,:)

           endif ! i_tilt_proc = 1

!          Transfer tilt to tilt arrays
           write(6,*)' Transfer tilt to tilt arrays'

!          radarName = 'KFTG'
!          c4_radarname = 'KFTG'

           elevationNumber = i_tilt_proc

!          Note the _HI scans are lowest elevations

!          High elevation (R)
           if(i_tilt_proc .gt. scanR_HI .and. 
     1        i_tilt_proc .le. scanR_HI + scanR)then
               i_array = i_tilt_proc - scanR_HI
               iscr = 0   
!              irmax = 0                                ! debug
               do j = 1,radialR
                 do i = 1,gateR
                   iscr = iscr + 1
                   Z(iscr) = Reflectivity(i,j,i_array)
!                  if(Z(iscr) .gt. 0)then               ! debug
!                    irmax = i                          ! debug
!                  endif                                ! debug
                 enddo ! i
!                if(irmax .gt. 0)then                   ! debug
!                 write(6,*)j,azimuthR(j,i_array),irmax ! debug
!    1                     ,distanceR(irmax)/1000.      ! debug
!                endif                                  ! debug
                 radialAzim(j) = azimuthR(j,i_array)
               enddo ! j
               elevationAngle = elevationR(1,i_array)
               numRadials = radialR
               ngates_ref_cdf = gateR
               firstGateRangeZ = distanceR(1) / 1000.
               gateSizeZ = (distanceR(2) - distanceR(1)) / 1000.
               write(6,*)' R i_tilt_proc/i_array/elev/gsp = ',
     1                     i_tilt_proc,i_array,elevationAngle,gateSizeZ
           endif

!          Low elevation (R_HI)
           if(i_tilt_proc .le. scanR_HI)then
               i_array = i_tilt_proc                   
               iscr = 0   
!              irmax = 0                                   ! debug
               if(radialR_HI .gt. MAX_RAY_TILT)then
                 write(6,*)' Error: radialR_HI > MAX_RAY_TILT',
     1                     radialR_HI,MAX_RAY_TILT
                 istatus = 0
                 return
               endif
               do j = 1,radialR_HI
                 do i = 1,gateR_HI
                   iscr = iscr + 1
                   Z(iscr) = Reflectivity_HI(i,j,i_array)
!                  if(Z(iscr) .gt. 0)then                  ! debug
!                    irmax = i                             ! debug
!                  endif                                   ! debug
                 enddo ! i
!                if(irmax .gt. 0)then                      ! debug
!                 write(6,*)j,azimuthR_HI(j,i_array),irmax ! debug
!    1                     ,distanceR_HI(irmax)/1000.      ! debug
!                endif                                     ! debug
                 radialAzim(j) = azimuthR_HI(j,i_array)
               enddo ! j
               elevationAngle = elevationR_HI(1,i_array)
               numRadials = radialR_HI
               ngates_ref_cdf = gateR_HI
               firstGateRangeZ = distanceR_HI(1) / 1000.
               gateSizeZ = (distanceR_HI(2) -  distanceR_HI(1)) / 1000.
               write(6,*)' R_HI i_tilt_proc/i_array/elev/gsp = ',
     1                     i_tilt_proc,i_array,elevationAngle,gateSizeZ
           endif
 
!          High elevation (V)
           i_v_match = 0
           do i_v = 1,scanV    
               if(elevationV(1,i_v) .eq. elevationAngle)then
                   i_v_match = i_v
                   write(6,*)' V elevation match ',i_v
               endif
           enddo

           if(i_v_match .gt. 0)then                              
               i_array = i_v_match                     
               iscr = 0   
               do j = 1,radialV
                 do i = 1,gateV
                   iscr = iscr + 1
                   V(iscr) = RadialVelocity(i,j,i_array)
                 enddo ! i
                 radialAzim(j) = azimuthV(j,i_array)         
               enddo ! j
               elevationAngle = elevationV(1,i_array)
               numRadials = radialV
               ngates_vel_cdf = gateV
               firstGateRangeV = distanceV(1) / 1000.
               gateSizeV = (distanceV(2) - distanceV(1)) / 1000.
               r_nyquist = nyquistVelocityV(i_array)
               write(6,*)' V i_tilt_proc/i_array/elev/nyq = ',
     1                     i_tilt_proc,i_array,elevationAngle,r_nyquist       
           else
               write(6,*)' No V match at elev ',elevationAngle
           endif

!          Low elevation (V_HI)
           i_v_match = 0
           do i_v = 1,scanV_HI    
               if(elevationV_HI(1,i_v) .eq. elevationAngle)then
                   i_v_match = i_v
                   write(6,*)' V_HI elevation match ',i_v
               endif
           enddo

           if(i_v_match .gt. 0)then                              
               i_array = i_v_match                     
               iscr = 0   
               do j = 1,radialV_HI
                 do i = 1,gateV_HI
                   iscr = iscr + 1
                   V(iscr) = RadialVelocity_HI(i,j,i_array)
                 enddo ! i
                 radialAzim(j) = azimuthV_HI(j,i_array)         
               enddo ! j
               elevationAngle = elevationV_HI(1,i_array)
               numRadials = radialV_HI
               ngates_vel_cdf = gateV_HI
               firstGateRangeV = distanceV_HI(1) / 1000.
               gateSizeV = (distanceV_HI(2) -  distanceV_HI(1)) / 1000. 
               r_nyquist = nyquistVelocityV_HI(i_array)
               write(6,*)' V_HI i_tilt_proc/i_array/elev/nyq = ',
     1                     i_tilt_proc,i_array,elevationAngle,r_nyquist       
           else
               write(6,*)' No V_HI match at elev ',elevationAngle
           endif

           if(i_tilt_proc .eq. scanR + scanR_HI)then
               i_last_scan = 1
           endif

!          if(i_tilt_proc .eq. scanV + scanV_HI)then
!              i_last_scan = 1
!          endif

           istatus = 1

           write(6,*)' End of volume processing for tilt ',i_tilt_proc       
     1                                                    ,i_last_scan

         endif ! end of volume processing for this tilt

       elseif(i_tilt_proc .le. 20)then ! l_exist is .false.
           i_tilt_proc = i_tilt_proc + 1
           i_skip = 1
           goto 200 ! check for existence of next tilt

       else
           istatus = 0

       endif ! file exists

       if(i_tilt_proc .gt. 50)then
           write(6,*)' ERROR in radar_init, i_tilt_proc = ',i_tilt_proc
           istatus = 0
           return
       endif

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

               if(i_skip .eq. 1)then ! tilt is found and prior one didn't exist
                   if(.false.)then
                       write(6,*)
     1               ' WARNING: We had to skip past some missing tilts'       
                   else
                       call s_len(radarName,lenr)
                       write(6,*)' ERROR: Missing tilts present for '
     1                          ,radarName(1:lenr),' at ',a9_time
     1                          ,' ',trim(laps_radar_ext)
                       istatus = 0
                       return
                   endif
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

       include 'remap_dims.inc'
       include 'netcdfio_radar_common.inc'
 
       get_number_of_gates = 0

       if(index .eq. 1)then
           get_number_of_gates = ngates_ref_cdf
       endif

       if(index .eq. 2)then
           get_number_of_gates = ngates_vel_cdf
       endif

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
!              write(6,*)' Warning: firstGateRangeV is outside range'
!    1                  ,firstGateRangeV
               first_gate_m = r_missing_data
           endif

           if(gateSizeV .ge. 0. .and. gateSizeV .lt. 1000.)then     
               gate_spacing_m = gateSizeV * 1000.
           else
!              write(6,*)' Warning: gateSizeV is outside range'
!    1                  ,gateSizeV      
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
 
       real data(n_gates)

       if(index .eq. 1)then ! reflectivity
           do i = 1,n_gates
               call counts_to_dbz(Z(n_ptr + (i-1))                       ! I
     1                           ,Z_scale,Z_offset,b_missing_data        ! I
     1                           ,data(i),istatus)                       ! O
           enddo

           if(istatus .ne. 1)then
               get_data_field = 0
               return
           endif

       elseif(index .eq. 2)then ! velocity
           do i = 1,n_gates
               call counts_to_vel(V(n_ptr + (i-1))                       ! I
     1                           ,b_missing_data,V_scale,V_offset        ! I
     1                           ,resolutionV,data(i))                   ! O
           enddo

       endif

       get_data_field = 1
       return
       end
 
       subroutine counts_to_dbz(zcounts,Z_scale,Z_offset,b_missing_data  ! I
     1                         ,dbz,istatus)                             ! O

!      Convert integer Z count value to dbz

!      From the NetCDF header
!      Z:valid_range = 2b, -2b ;     (2 through 254)
!      Z:below_threshold = 0b ;      (0)
!      Z:range_ambiguous = 1b ;      (1)
!      Z:_FillValue = -1b ;          (255)

       integer zcounts
       real dbz,dbz_hold,b_missing_data                         

       dbz_hold = zcounts

!      Convert from signed to unsigned
       if(dbz_hold .gt. 127.) then
           print *, 'error in Reflectivity: ',dbz_hold,zcounts
           istatus = 0
       endif

       if(dbz_hold .lt. 0.) then
           dbz_hold = 256. + dbz_hold
       endif

!      if(dbz_hold .eq. 1.)then 
!          dbz_hold = b_missing_data  ! Range Ambiguous
!      endif

       if(dbz_hold .ne. b_missing_data)then ! Scale
           dbz_hold = (dbz_hold - Z_offset) / Z_scale
       endif

       dbz = dbz_hold

       istatus = 1
       return
       end


       subroutine counts_to_vel(vcounts,b_missing_data,V_scale          ! I
     1                         ,V_offset,resolutionV,vel_ms)            ! O

!      Convert integer V count value to radial velocity (meters/sec)

       integer vcounts
       real vel_ms,vel_hold,b_missing_data                         

       vel_hold = vcounts

!      Convert from signed to unsigned
       if(vel_hold .gt. 127.) then
           print *, 'error in Velocity: ',vel_hold
           stop
       endif

       if(vel_hold .lt. 0.) then
           vel_hold = 256. + vel_hold
       endif

       if(vel_hold .eq. 1. .or. vel_hold .eq. 0.)then 
           vel_hold = b_missing_data  ! Invalid Measurement
       endif

       if(resolutionV .eq. 0.)then ! QC Check
           vel_hold = b_missing_data
       endif

       if(vel_hold .ne. b_missing_data)then ! Scale valid V
           vel_hold = (vel_hold - V_offset) / V_scale
       endif

       vel_ms = vel_hold
   
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

!      if(c8_fname_format .eq. 'VOLUME')then
!          a15_time = fname9_to_wfo_fname15(a9_time)
!          filename = path_to_radar(1:len_path)//'/KFTG'//a15_time//'_V03.nc'
!    1                //c2_tilt
!      endif

       call s_len(filename,len_file)
       write(6,*)' check_input_file: ',filename(1:len_file),' ',l_exist       

       return
       end
