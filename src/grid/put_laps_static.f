
      subroutine put_laps_static(grid_spacing,model,comment,data
     1                          ,imax,jmax, mkmax
     1                          ,std_lat,std_lat2,std_lon      
     1                          ,c6_maproj,deltax,deltay)
 
      integer*4 mkmax

      INTEGER*4	IMAX,		        !I4time of data
     1		JMAX,KMAX,	        !# cols, # rows, # fields
     1		KDIM,i,j,k,		!K dimension of DATA array
     1		ISTATUS
C
      REAL*4	DATA(imax,jmax,mkmax),	!Raw data to be written
     1		grid_spacing,val
C
      CHARACTER*50	DIR_in		!Directory to be written to
      CHARACTER*50	DIR_out		!Directory to be written to
      CHARACTER*31	EXT		!File name ext (up to 31 chars)
      CHARACTER*3	VAR(mkmax)	!3 letter ID of each field
      CHARACTER*10	UNITS(mkmax)	!units of each field
      CHARACTER*(*)	COMMENT(mkmax)	!Comments for each field
      character*(*)     model
      character*80      origin          !Run time parameter - c80_description
      character*9       laps_dom_file
      character*6       c6_maproj       !Map projection
      integer len      

      call get_directory('static',dir_out,len) 
      ext = 'nest7grid'
      kmax = mkmax-3
      kdim = mkmax-3

      laps_dom_file = 'nest7grid'

      call get_c80_description(origin,istatus)

      var(1) = 'LAT'
      var(2) = 'LON'
      var(3) = 'AVG'
      var(4) = 'LDF'
      var(5) = 'ZIN'

!     Do zin calc (note this is 5th element in data array)
      do i = 1,imax
      do j = 1,jmax
          psa = ztopsa(data(i,j,3)) ! This is the AVG data (3rd element)
          data(i,j,5) = (20.0 - ((psa - 100.0) * 0.02))
      enddo ! j
      enddo ! i
      
      

      write(6,*) dir_out(1:len),len
      call wrt_laps_static (dir_out(1:len),laps_dom_file,imax,jmax,
     1                      kmax,deltax,deltay,std_lon,std_lat,
     1                      std_lat2,origin,var,comment,
     1                      data,model,grid_spacing,
     1                      c6_maproj,istatus)

      write (6,*)'wrt_laps_static:status = ',istatus

      call rd_laps_static(dir_out(1:len),laps_dom_file,imax,jmax,kdim,
     1                    var,units,comment,data,grid_spacing,istatus)
      
      write (6,*)'rd_laps_static:status = ',istatus

      return
      end

