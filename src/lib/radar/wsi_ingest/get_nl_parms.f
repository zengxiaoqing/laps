      subroutine get_ln3_parameters(msng_radar,ickint,
     +itotwait,iageth,istart,jstart,iend,jend,istatus)
c
      implicit none

      integer    istatus
      integer    len_dir
      integer    istart,jstart
      integer    iend,jend
      integer    ickint,itotwait,iageth
      integer    msng_radar

      character  nest7grid*150

      namelist /ln3_nl/ msng_radar,ickint,itotwait,
     + iageth,istart,jstart,iend,jend

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
      subroutine get_wsi_parms_vrc(irad,lines,elems,
     +dx,dy,rla1,rlo1,rla2,rlo2,rlat,rlon,rlatin,istatus)

      integer nrad_types
      parameter (nrad_types=2)

      integer nlines(nrad_types)
      integer nelems(nrad_types)
      integer irad
      integer lines,elems
      integer istatus

      real*4  resx(nrad_types)
      real*4  resy(nrad_types)
      real*4  radla1(nrad_types)
      real*4  radlo1(nrad_types)
      real*4  radla2(nrad_types)
      real*4  radlo2(nrad_types)
      real*4  radlat(nrad_types)
      real*4  radlon(nrad_types)
      real*4  radlatin(nrad_types)

      real*4  rla1,rlo1,rlat,rlon,rlatin
      real*4  dx,dy

      character nest7grid*150

      namelist /vrc_nl/nlines,nelems,resx,resy,radla1,
     +radlo1,radla2,radlo2,radlat,radlon,radlatin

      istatus = 0

      call get_directory('nest7grid',nest7grid,len_dir)

      nest7grid = nest7grid(1:len_dir)//'/vrc.nl'

      open(1,file=nest7grid,status='old',err=900)
      read(1,vrc_nl,err=901)

      lines=nlines(irad)
      elems=nelems(irad)
      dx=resx(irad)
      dy=resy(irad)
      rla1=radla1(irad)
      rlo1=radlo1(irad)
      rla2=radla2(irad)
      rlo2=radlo2(irad)
      rlat=radlat(irad)
      rlon=radlon(irad)
      rlatin=radlatin(irad)

      close(1)
      istatus = 1
      return

 900  print*,'error opening file ',nest7grid
      return
 901  print*,'error reading vrc_nl in ',nest7grid
      write(*,vrc_nl)
      return
      end
