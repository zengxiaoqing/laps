      subroutine getradarinfo(nradars,max_radars_natl,
     &radar_id,radar_name,radar_lat,radar_lon,radar_elev)
c
c routine reads static/radarinfo.dat file and returns to the calling
c program the 4 character radar_id, the 8 character radar_name,
c the real radar_lat and lon and the integer radar_elev in meters
c
      implicit none

      integer    i
      integer    lat_deg,lat_min,lat_sec
      integer    lon_deg,lon_min,lon_sec
      integer    max_radars_natl

      real       rlat_deg,rlat_min,rlat_sec
      real       rlon_deg,rlon_min,rlon_sec

      real       radar_lat(max_radars_natl)
      real       radar_lon(max_radars_natl)

      integer    radar_elev(max_radars_natl)

      logical      end_of_file

      character    radar_id(max_radars_natl)*4
      character    radar_name(max_radars_natl)*8

c ======================================================================
      open(88,file=../../../../../static/radarinfo.dat,form='formatted',
     &     status='old',err=999)

      read(88,*)     !header line

      end_of_file=.false.
      i=1
      do while (.not.end_of_file)

         read(88,90,end=99)radar_id(i),radar_name(i),lat_deg,
     &lat_min,lat_sec,lon_deg,lon_min,lon_sec,radar_elev(i)
c
c convert to fraction
c
         if(lat_deg.ne.0 .and. lon_deg.ne.0)then

            rlat_min=float(lat_min)/60.
            rlat_sec=float(lat_sec)/3600.
            rlon_min=float(lon_min)/60.
            rlon_sec=float(lon_sec)/3600.

            radar_lat(i)=float(lat_deg)+rlat_min+rlat_sec
            radar_lon(i)=float(lon_deg)+rlon_min+rlon_sec
            i=i+1

         endif

         goto 100

99       end_of_file=.true.

100   enddo
      nradars=i

90    format(a4,8x,a8,3(6x,i2),7x,i3,5x,i2,6x,i2,6x,i4)
      return
      end
