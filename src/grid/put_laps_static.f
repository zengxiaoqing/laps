
c     Program Gridmap

c      include 'lapsparms.for'
c     character*125 comment(4)
c     character*131 model
c     real*4 data(NX_L,NY_L,4)

c     grid_spacing = 10000.0

c     call put_laps_static(grid_spacing,model,comment,data)

c     end


      subroutine put_laps_static(grid_spacing,model,comment,data,nx,ny)       

      implicit none

!     include 'lapsparms.for'
               
      integer*4 NX,NY,mkmax
! make sure these match the values in laps_grid_def.h
!     parameter (NX = NX_L)
!     parameter (NY = NY_L)
      parameter (mkmax = 7)

      INTEGER*4	IMAX,		!I4time of data
     1		JMAX,KMAX,	!# cols, # rows, # fields
     1		KDIM,i,j,k,		!K dimension of DATA array
     1		ISTATUS
C
      REAL*4	DATA(NX,NY,mkmax),	!Raw data to be written
     1          zin(NX,NY),
     1		grid_spacing,val
C
      CHARACTER*50	DIR_in		!Directory to be written to
      CHARACTER*50	DIR_out		!Directory to be written to
      CHARACTER*31	EXT		!File name ext (up to 31 chars)
      CHARACTER*3	VAR(mkmax)	        !3 letter ID of each field
      CHARACTER*10	UNITS(mkmax)	!units of each field
      CHARACTER*125	COMMENT(mkmax)	!Comments for each field
      character*131     model
      character*9       laps_dom_file
      integer len      
      imax = NX
      jmax = NY 

      var(1) = 'LAT'
      var(2) = 'LON'
      var(3) = 'AVG'
      var(4) = 'LDF'

      call get_directory('static',dir_out,len) 
      ext = 'nest7grid'
      kmax = mkmax-3
      kdim = mkmax-3

c      dir_out = './model/'
c     dir_out = '/home/navaho/wharton/laps/'
      laps_dom_file = 'nest7grid'


c     write(6,*) var(1),' ', var(2),' ', var(3),' ', var(4)
c     call rd_laps_static(dir_in,laps_dom_file,imax,jmax,kdim,
c    1                    var, units, comment, data, 
c    1                    grid_spacing,istatus)
c     write(6,*) var(1),' ', var(2),' ', var(3),' ', var(4)
c     write (6,*),'rd_laps_static:status = ',istatus

c     val = 0.0
c     do k = 1,mkmax
c       val = val + 1.0
c       do i = 1, NX
c         do j = 1, NY
c           data(i,j,k) = val
c         enddo
c       enddo
c     enddo

      model = 'RAMS 4 delta x smoothed filter'
c     write(6,*) var(1),' ', var(2),' ', var(3),' ', var(4)
      call wrt_laps_static (dir_out(1:len),laps_dom_file,imax,jmax,kmax,
     1                      var,comment,data,zin,model,
     1                      grid_spacing,istatus)
      
      write (6,*)'wrt_laps_static:status = ',istatus

      call rd_laps_static(dir_out(1:len),laps_dom_file,imax,jmax,kdim,
     1                    var,units,comment,data,grid_spacing,istatus)
      
      write (6,*)'rd_laps_static:status = ',istatus

      return
      end
