      subroutine get_ln3_parameters(msng_radar,ickint,
     +itotwait,iageth,istatus)
c
      implicit none

      integer    istatus
      integer    len_dir
      integer    ickint,itotwait,iageth
      integer    msng_radar

      character  nest7grid*150

      namelist /ln3_nl/msng_radar,ickint,itotwait,iageth

      istatus = 1

      call get_directory('nest7grid',nest7grid,len_dir)
      if(nest7grid(len_dir:len_dir).ne.'/') then
        len_dir=len_dir+1
        nest7grid(len_dir:len_dir)='/'
      endif

      nest7grid = nest7grid(1:len_dir)//'ln3.nl'

      open(1,file=nest7grid,status='old',err=900)
      read(1,ln3_nl,err=901)
      close(1)
      istatus=0
      return
900   print*,'error opening file ',nest7grid
      stop
901   print*,'error reading ln3.nl in ',nest7grid
      write(*,ln3_nl)
      stop
      end
c
c------------------------------------------------------
c
      subroutine read_vrc_nl(iwsimsng,icheckint,iwaittime,
     +iagethresh,istatus)

      integer imsngrad
      integer icheckint
      integer iwaittime
      integer iagethresh
      integer istatus

      character nest7grid*150

      namelist /vrc_nl/iwsimsng,icheckint,iwaittime,iagethresh

      istatus = 1

      call get_directory('nest7grid',nest7grid,len_dir)
      if(nest7grid(len_dir:len_dir).ne.'/') then
        len_dir=len_dir+1
        nest7grid(len_dir:len_dir)='/'
      endif

      nest7grid = nest7grid(1:len_dir)//'vrc.nl'

      open(1,file=nest7grid,status='old',err=900)
      read(1,vrc_nl,err=901)

      close(1)
      istatus = 1
      return

 900  print*,'error opening file ',nest7grid
      return
 901  print*,'error reading vrc_nl in ',nest7grid
      write(*,vrc_nl)
      return
      end
