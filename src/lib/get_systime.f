      subroutine get_systime(i4time_sys,a9_time,istatus)
      character*9 a9_time
      character*100 dir
      integer i4time_sys, len, istatus

      call get_directory('etc',dir,len)
c      print *, dir(1:len), len      

      open(11,file=dir(1:len)//'systime.dat',status='old')

      read(11,*,err=999)i4time_sys

c      print *,i4time_sys
      read(11,2,err=999)a9_time
c      print *, a9_time
2     format(1x,a9)

      close(11)

      istatus = 1

      return

 999  print*,'Error reading systime file'
      istatus = 0
      return
      end
