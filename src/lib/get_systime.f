      subroutine get_systime(i4time_sys,a9_time,istatus)

      integer*4 i4time_sys
      character*9 a9_time
      integer *4 istatus

      character*100 dir
      integer *4 length

      call get_directory('time',dir,length)

      open(11,file=dir(1:length)//'systime.dat',status='old')

      read(11,*,err=999)i4time_sys
      read(11,2,err=999)a9_time
2     format(1x,a9)
      close(11)
      istatus = 1
      return

 999  print*,'Error reading systime file'
      istatus = 0
      return
      end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      subroutine get_systime_i4(i4time_sys,istatus)

      integer*4 i4time_sys
      integer*4 istatus

      character*100 dir
      integer*4 length

      call get_directory('time',dir,length)

      open(11,file=dir(1:length)//'systime.dat',status='old',err=999)

      read(11,*,err=999)i4time_sys
      close(11)
      istatus = 1
      return

 999  print*,'Error reading systime file'
      istatus = 0
      return
      end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      subroutine get_systime_all(i4time_sys,a9_time,analysis_hr,
     1             analysis_min,asctim_str,yr_jday,istatus)

      integer*4 i4time_sys
      character*9 a9_time
      character*(*) analysis_hr, analysis_min
      character*(*) asctim_str
      character*(*) yr_jday
      integer*4 istatus

      character*30 asc_str
      character*10 yjd_str
      character*100 dir
      integer*4 length, loc_len, loc_asc_len, found, i

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

