 
       subroutine radar_init()       ! Open Polar NetCDF file for the proper time
 
       integer max_files
       parameter(max_files=1000)

       character*150 path_to_wideband,c_filespec,filename
     1              ,c_fnames(max_files)
       character*9 a9_time
       integer*4 i4times(max_files),i4times_lapsprd(max_files)

       integer*2 VCP, elevationNumber
       common/radar_housekeeping/
     1                                siteLat                        
     1                               ,siteLon                        
     1                               ,siteAlt                        
     1                               ,elevationAngle
     1                               ,elevationNumber
     1                               ,VCP
     1                               ,i4time_process

 !     Get path to file
       call get_remap_parms(path_to_wideband,istatus)
       call s_len(path_to_wideband,len_path)
 
 !     Get i4time of 01 elevation file nearest to 15 minutes ago
       i4time_now = i4time_now_gg() 
       c_filespec = path_to_wideband//'/*_elev01'

       call get_file_times(c_filespec,max_files,c_fnames
     1                    ,i4times,i_nbr_files_out,istatus)
       if(istatus .ne. 1)then
           stop
       endif

       call get_filespec('v01',1,c_filespec,istatus)
       call get_file_times(c_filespec,max_files,c_fnames
     1                   ,i4times_lapsprd,i_nbr_lapsprd_files,istatus)

       if(i_nbr_files_out .ge. 2)then
           i4time_process = i4times(i_nbr_files_out-1)
           call make_fnam_lp(i4time_process,a9_time,istatus)
           do i = 1,i_nbr_lapsprd_files
               if(i4time_process .eq. i4times(i))then
                   write(6,*)' Product file already exists ',a9time
               endif
           enddo ! i
       else
           write(6,*)' # of files = ',i_nbr_files_out
       endif

!      Pull in housekeeping data from 1st tilt
       filename = path_to_wideband(1:len_path)//'/'//a9_time//'_elev01'
       write(6,*)' We will read this file for housekeeping: '
       write(6,*)filename

       call get_tilt_netcdf_data(filename
     1                               ,siteLat                        
     1                               ,siteLon                        
     1                               ,siteAlt                        
     1                               ,elevationAngle
     1                               ,elevationNumber
     1                               ,VCP
     1                               ,istatus)
      
       return
       end
 
       subroutine get_remap_parms(path_to_wideband,istatus)

       character*150 path_to_wideband,filename
       namelist /remap_nl/ path_to_wideband
 
       character*150 static_dir,filename
 
       call get_directory('nest7grid',static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/nest7grid.parms'
 
       open(1,file=filename,status='old',err=900)
       read(1,remap_nl,err=901)
       close(1)

       istatus = 1
       return

  900  print*,'error opening file ',filename
       istatus = 0
       return

  901  print*,'error reading remap_nl in ',filename
       istatus = 0
       return

       end

 
       function get_altitude()
       integer get_altitude          ! Site altitude (meters)

       integer*2 VCP, elevationNumber
       common/radar_housekeeping/
     1                                siteLat                        
     1                               ,siteLon                        
     1                               ,siteAlt                        
     1                               ,elevationAngle
     1                               ,elevationNumber
     1                               ,VCP
     1                               ,i4time_process
 
       get_altitude = nint(siteAlt)

       return
       end
 
 
       function get_latitude()
       integer get_latitude          ! Site latitude (degrees * 100000)

       integer*2 VCP, elevationNumber
       common/radar_housekeeping/
     1                                siteLat                        
     1                               ,siteLon                        
     1                               ,siteAlt                        
     1                               ,elevationAngle
     1                               ,elevationNumber
     1                               ,VCP
     1                               ,i4time_process
 
       get_latitude = nint(siteLat*100000)
       return
       end
 
 
       function get_longitude()
       integer get_longitude         ! Site longitude (degrees * 100000)

       integer*2 VCP, elevationNumber
       common/radar_housekeeping/
     1                                siteLat                        
     1                               ,siteLon                        
     1                               ,siteAlt                        
     1                               ,elevationAngle
     1                               ,elevationNumber
     1                               ,VCP
     1                               ,i4time_process
 
       get_longitude = nint(siteLon*100000)
       return
       end
 
 
       function get_field_num(ext)
       integer get_field_num
       character*3 ext
 
       if(ext .eq. 'DBZ')get_field_num = 1
       if(ext .eq. 'VEL')get_field_num = 2
 
       return
       end
 
 
       function read_radial()
 
       read_radial = 0
       return
       end
 
 
       function get_status()
 
       get_status = 0
       return
       end
 
 
       function get_fixed_angle()
       integer get_fixed_angle     ! Beam tilt angle (degrees * 100)

       integer*2 VCP, elevationNumber
       common/radar_housekeeping/
     1                                siteLat                        
     1                               ,siteLon                        
     1                               ,siteAlt                        
     1                               ,elevationAngle
     1                               ,elevationNumber
     1                               ,VCP
     1                               ,i4time_process
 
       get_fixed_angle = nint(elevationAngle * 100.)
       return
       end
 
 
       function get_scan()
       integer get_scan            ! Scan #

       integer*2 VCP, elevationNumber
       common/radar_housekeeping/
     1                                siteLat                        
     1                               ,siteLon                        
     1                               ,siteAlt                        
     1                               ,elevationAngle
     1                               ,elevationNumber
     1                               ,VCP
     1                               ,i4time_process
 
       get_scan = elevationNumber
       return
       end
 
 
       function get_tilt()
       integer get_tilt            ! Tilt #

       integer*2 VCP, elevationNumber
       common/radar_housekeeping/
     1                                siteLat                        
     1                               ,siteLon                        
     1                               ,siteAlt                        
     1                               ,elevationAngle
     1                               ,elevationNumber
     1                               ,VCP
     1                               ,i4time_process
 
       get_tilt = elevationNumber
       return
       end
 
 
       function get_volume_time()
       integer get_volume_time

       integer*2 VCP, elevationNumber
       common/radar_housekeeping/
     1                                siteLat                        
     1                               ,siteLon                        
     1                               ,siteAlt                        
     1                               ,elevationAngle
     1                               ,elevationNumber
     1                               ,VCP
     1                               ,i4time_process
 
       get_volume_time = i4time_process
       return
       end
 
 
       function get_vcp()
       integer get_vcp

       integer*2 VCP, elevationNumber
       common/radar_housekeeping/
     1                                siteLat                        
     1                               ,siteLon                        
     1                               ,siteAlt                        
     1                               ,elevationAngle
     1                               ,elevationNumber
     1                               ,VCP
     1                               ,i4time_process
 
       get_vcp = VCP
       return
       end
 
 
       function get_azi()
 
       get_azi = 0
       return
       end
 
 
       function get_nyquist()
       integer get_nyquist        ! Nyquist velocity of the radial (M/S*100)
 
       get_nyquist = 0
       return
       end
 
 
       function get_number_of_gates(index)
       integer get_number_of_gates
 
       get_number_of_gates = 0
       return
       end
 
 
       function get_first_gate()
 
       get_first_gate = 0
       return
       end
 
 
       function get_data_field(index, n_gates)
 
       get_data_field = 0
       return
       end
 
 
       function cvt_fname_data()
 
       cvt_fname_data = 0
       return
       end
 
 
