      subroutine get_balance_nl(lrunbal,lstagger,icon,gamo,delo,tau,
     .                          lnon_linear,istatus)
c
      implicit none

      integer    icon
      integer    istatus
      integer    len_dir
      logical    lrunbal
      logical    lstagger
      logical    lnon_linear
      real*4     gamo,delo,tau
      character  nest7grid*150

      namelist /balance_nl/lrunbal,lstagger,icon,gamo,delo,tau,
     1lnon_linear

      istatus = 0

      call get_directory('nest7grid',nest7grid,len_dir)
      if(nest7grid(len_dir:len_dir).ne.'/') then
        len_dir=len_dir+1
        nest7grid(len_dir:len_dir)='/'
      endif

      nest7grid = nest7grid(1:len_dir)//'balance.nl'

      open(1,file=nest7grid,status='old',err=900)
      read(1,balance_nl,err=901)
      close(1)
      return
900   print*,'error opening file ',nest7grid
      stop
901   print*,'error reading balance.nl in ',nest7grid
      write(*,balance_nl)
      stop
      end

