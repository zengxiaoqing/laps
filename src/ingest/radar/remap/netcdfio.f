 
       subroutine radar_init()       ! Open Polar NetCDF file for the proper time
 
       character*150 path_to_wideband,filename
       character*9 a9_time

 !     Get path to file
       call get_remap_parms(path_to_wideband,istatus)
       call s_len(path_to_wideband,len_path)
 
 !     Get i4time of file
       i4time_now = i4time_now_gg()
       call get_file_time(path_to_wideband,i4time_now,i4time_nearest)
       call make_fnam_lp(i4time_nearest,a9_time,istatus)
       if(istatus .ne. 1)then
           stop
       endif

       filename = path_to_wideband(1:len_path)//'/'//a9_time
       write(6,*)' We should here open this file: ',filename
      
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
       integer get_altitude          ! Site altitude (meters * 100000)
 
       get_altitude = 0
       return
       end
 
 
       function get_latitude()
       integer get_latitude          ! Site latitude (meters * 100000)
 
       get_latitude = 0
       return
       end
 
 
       function get_longitude()
       integer get_longitude         ! Site longitude (meters * 100000)
 
       get_longitude = 0
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
 
       get_fixed_angle = 0
       return
       end
 
 
       function get_scan()
       integer get_scan            ! Scan #
 
       get_scan = 0
       return
       end
 
 
       function get_tilt()
       integer get_tilt            ! Tilt #
 
       get_tilt = 0
       return
       end
 
 
       function get_year()
 
       get_year = 0
       return
       end
 
 
       function get_month()
 
       get_month = 0
       return
       end
 
       function get_day()
 
       get_day = 0
       return
       end
 
 
       function get_hour()
 
       get_hour = 0
       return
       end
 
 
       function get_min()
 
       get_min = 0
       return
       end
 
 
       function get_sec()
 
       get_sec = 0
       return
       end
 
 
       function get_vcp()
       integer get_vcp
 
       get_vcp = 0
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
 
 
