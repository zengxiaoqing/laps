      subroutine get_systime(i4time_sys,a9_time,istatus)

      integer i4time_sys
      character*9 a9_time,a9_time_save
      integer istatus

      character*100 dir
      integer length

      integer init,i4time_sys_save
      data init/0/
      save init,i4time_sys_save,a9_time_save

      if(init .eq. 0)then
          call GETENV('LAPS_A9TIME',a9_time)
          call s_len(a9_time,ilen)

          if(ilen .eq. 9)then ! override file using environment variable
              write(6,*)' systime (from env) = ',a9_time
              call i4time_fname_lp(a9_time,i4time_sys,istatus)

          else ! read systime from file
              call get_directory('time',dir,length)

              open(11,file=dir(1:length)//'systime.dat',status='old')

              read(11,*,err=999)i4time_sys
              read(11,2,err=999)a9_time
2             format(1x,a9)
              close(11)
              write(6,*)' systime (from systime.dat) = ',a9_time
          endif

          init = 1
          i4time_sys_save = i4time_sys
          a9_time_save = a9_time

      else ! used saved variables for efficiency
          i4time_sys = i4time_sys_save
          a9_time = a9_time_save
          
      endif

      istatus = 1
      return

 999  print*,'Error reading systime file'
      istatus = 0
      return

      end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      subroutine get_systime_i4(i4time_sys,istatus)

      integer i4time_sys
      integer istatus

      character*9 a9_time

      call get_systime(i4time_sys,a9_time,istatus)

      return
      end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      subroutine get_systime_all(i4time_sys,a9_time,analysis_hr,
     1             analysis_min,asctim_str,yr_jday,istatus)

      integer i4time_sys
      character*9 a9_time
      character*(*) analysis_hr, analysis_min
      character*(*) asctim_str
      character*(*) yr_jday
      integer istatus

      character*30 asc_str
      character*10 yjd_str
      character*100 dir
      integer length, loc_len, loc_asc_len, found, i

      call get_directory('time',dir,length)

      open(11,file=dir(1:length)//'systime.dat',status='old')

      read(11,*,err=999)i4time_sys
      read(11,2,err=999)a9_time
      read(11,4,err=999)analysis_hr
      read(11,4,err=999)analysis_min
      read(11,4,err=999)asc_str
      read(11,4,err=999)yjd_str

2     format(1x,a9)
4     format(a)
      close(11)

      length = len(asctim_str)
      loc_asc_len = len(asc_str)
      loc_len = 0
      found = 0
      i = 1
      do while ((i .lt. loc_asc_len) .and. (found .eq. 0))
        if (asc_str(i:i+1) .eq. '  ') then
          found = 1
          loc_len = i-1
        endif 
        i = i + 1
      enddo
      if (loc_len .eq. 0) loc_len = loc_asc_len
      if (loc_len .gt. length) then 
        print*,'Error reading systime file - asctime string truncated'
        asctim_str = asc_str(1:length)
      else
        asctim_str = asc_str(1:loc_len)
      endif

      length = len(yr_jday)
      loc_asc_len = len(yjd_str)
      loc_len = 0
      found = 0
      i = 1
      do while ((i .lt. loc_asc_len) .and. (found .eq. 0))
        if (yjd_str(i:i+1) .eq. '  ') then
          found = 1
          loc_len = i-1
        endif 
        i = i + 1
      enddo
      if (loc_len .eq. 0) loc_len = loc_asc_len
      if (loc_len .gt. length) then
        print*,'Error reading systime file - YYJJJ string truncated'
        yr_jday = yjd_str(1:length)
      else
        yr_jday = yjd_str(1:loc_len)
      endif

      istatus = 1
      return

 999  print*,'Error reading systime file'
      istatus = 0
      return
      end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine get_modeltime(i4time_mdl,a9_time,istatus)

      implicit none 
      integer i4time_mdl
      character*9 a9_time,a9_time_save
      integer            :: istatus
      integer, parameter :: lun=120
      logical            :: file_exists
      character*300 laps_data_root,dir_t,filenamet
      integer len_dir_t

      call get_directory('time',dir_t,istatus)
      call s_len(dir_t,len_dir_t)
      filenamet = dir_t(1:len_dir_t)//'/modeltime.dat'
      inquire(file=filenamet,exist=file_exists)
      if (file_exists) then
         open(lun,file=filenamet,status='old')
         read(lun,*) a9_time
         close(lun)
         !print *, ' modeltime (from modeltime.dat) = ',a9_time
         call i4time_fname_lp(a9_time,i4time_mdl,istatus)
         if (istatus .ne. 1) then
           print *, "-- Error getting LAPS modeltime.dat!"
           stop "STOP in get_laps_modeltime"
         endif
      else
        print *, "-- No modeltime.dat file found: ",trim(filenamet)
        stop "STOP in get_laps_modeltime"
      endif
      end subroutine get_modeltime
