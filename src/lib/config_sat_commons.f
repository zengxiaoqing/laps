      subroutine config_satellite_lvd(istatus)

      character nest7grid*150
      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'
      include 'satellite_namelist_lvd.cmn'
      include 'satdata_lvd.for'

      istatus = 0

      call get_directory('nest7grid',nest7grid,len_dir)

      nest7grid = nest7grid(1:len_dir)//'/satellite_lvd.nl'

      open(1,file=nest7grid,status='old',err=900)
      read(1,satellite_lvd_nl,err=901)
      close(1)
      istatus = 1
      return

 900  print*,'error opening file ',nest7grid
      return
 901  print*,'error reading satellite_nl in ',nest7grid
      write(*,satellite_lvd_nl)
      return
      end
c
c =-=-===-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c
      subroutine get_sat_sounder_info(n_sndr,c_sndr_id,
     +n_sndr_channels,path_to_sat_sounder,n_elems,n_lines,
     +channel_wavelength_u,imsng_sndr_pix,istatus)

      integer maxsndr
      integer maxch
      parameter (maxsndr=4,maxch=19)

      integer len_dir
      integer n_sndr
      integer n_sndr_channels
      integer n_elems(maxsndr)
      integer n_lines(maxsndr)
      integer imsng_sndr_pix
      integer istatus
      character*6   c_sndr_id(maxsndr)
      character*150 nest7grid
      character*200 path_to_sat_sounder(maxsndr)

      real*4 channel_wavelength_u(maxch,maxsndr)

      namelist /satellite_sounder_nl/ n_sndr,c_sndr_id,path_to_sat_sound
     +er,n_elems,n_lines,n_sndr_channels,channel_wavelength_u,imsng_sndr
     +_pix

      call get_directory('nest7grid',nest7grid,len_dir)

      nest7grid = nest7grid(1:len_dir)//'/sat_sounder.nl'

      open(1,file=nest7grid,status='old',err=900)
      read(1,satellite_sounder_nl,err=901)
      close(1)

      istatus = 1
      return
 900  print*,'error opening file ',nest7grid
      stop
 901  print*,'error reading satellite_sounder_nl in ',nest7grid
      write(*,satellite_sounder_nl)
      stop
      end
