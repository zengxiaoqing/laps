cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
      subroutine radar_init(drive)
      implicit none
      character*80 drive
      integer iyear,imonth,iday,ihour,imin,isec
      common /timeinfo/  iyear,imonth,iday,ihour,imin,isec
      integer iazi,iangle,itilt,iscan
      common /counters/ iazi,iangle,itilt,iscan
      print *, ' '
      print *, ' Initializing a2io_fake for file testing only'
      print *, ' '
      iazi=0
      iangle=100
      itilt=1
      iscan=1
      return
      end
c
      subroutine set_radar_name(site)
      implicit none
      character*80 site
      character*80 csite
      csite=site
      return
      end
c
      function read_radial()
c
      implicit none
      integer read_radial
      integer itime,istatus
      integer iyear,imonth,iday,ihour,imin,isec
      integer time
      common /timeinfo/  iyear,imonth,iday,ihour,imin,isec
      integer iazi,iangle,itilt,iscan
      common /counters/ iazi,iangle,itilt,iscan
      integer maxtilt
      parameter (maxtilt=14)
      integer kelev(maxtilt)
      data kelev /50,150,240,340,430,530,620,750,
     :            870,1000,1200,1400,1670,1950/
      iazi=iazi+100
      IF(iazi.ge.360 00) THEN
        iazi=iazi-360 00
        itilt=itilt+1
c       call system('sleep 1')
        IF(itilt.gt.maxtilt) THEN
          itilt=1
          iscan=iscan+1
        END IF
      END IF
      iangle=kelev(itilt)
      itime= time()
c      print *, ' itime = ',itime
      call cv_i4tim_int_lp(itime,iyear,imonth,iday,ihour,
     :                     imin,isec,istatus)
      iyear=iyear+10
c      print *, ' day,mon,year : ',iday,imonth,iyear
c      print *, ' hour,min,sec : ',ihour,imin,isec
      read_radial=0
      return
      end
c
      function get_status(index)
      implicit none
      integer get_status
      integer index
      get_status=1
      return
      end
c
      subroutine get_data(ref,vel,spw,snr)
      implicit none
      integer nsize
      parameter (nsize=920)
      real ref(*),vel(*)
      real spw(*),snr(*)
      integer i
      DO 100 i=1,nsize
        ref(i)=35.
        vel(i)=8.
        spw(i)=3.
        snr(i)=10.
  100 CONTINUE
      return
      end
c
      function get_data_field(ifield,field,nsize)
      implicit none
      integer get_data_field
      integer nsize
      integer ifield
      real field(*)
      integer iazi,iangle,itilt,iscan
      common /counters/ iazi,iangle,itilt,iscan
      integer i
      IF(ifield.eq.0) THEN
        DO 100 i=1,460
          field(i)=0.1*float(i)
  100   CONTINUE
      ELSE
        DO 200 i=1,920
          field(i)=0.025*float(i)
  200   CONTINUE
      END IF
      get_data_field=0
      return
      end
c
      function get_rt_mode()
      implicit none
      integer get_rt_mode
      get_rt_mode=1
      return
      end
c
      function get_field_num(fchar)
      implicit none
      integer get_field_num
      character*3 fchar
      IF(fchar.eq.'DBZ') THEN
        get_field_num=0
      ELSE IF(fchar.eq.'VEL') THEN
        get_field_num=1
      ELSE IF(fchar.eq.'SPW') THEN
        get_field_num=2
      ELSE
        get_field_num=3
      END IF
      return
      end
c
      function get_azi()
      implicit none
      integer iazi,iangle,itilt,iscan
      common /counters/ iazi,iangle,itilt,iscan
      integer get_azi
      get_azi=iazi
      return
      end
c
      function get_fixed_angle()
      implicit none
      integer iazi,iangle,itilt,iscan
      common /counters/ iazi,iangle,itilt,iscan
      integer get_fixed_angle
      get_fixed_angle=iangle
      return
      end
c
      function get_scan()
      implicit none
      integer get_scan
      integer iazi,iangle,itilt,iscan
      common /counters/ iazi,iangle,itilt,iscan
      get_scan=iscan
      return
      end
c
      function get_tilt()
      implicit none
      integer get_tilt
      integer iazi,iangle,itilt,iscan
      common /counters/ iazi,iangle,itilt,iscan
      get_tilt=itilt
      return
      end
c
      function get_nyquist()
      implicit none
      integer get_nyquist
      get_nyquist=35 00
      return
      end
c
      function get_year()
      implicit none
      integer get_year
      integer iyear,imonth,iday,ihour,imin,isec
      common /timeinfo/  iyear,imonth,iday,ihour,imin,isec
      get_year=iyear 
      return
      end
c
      function get_month()
      implicit none
      integer get_month
      integer iyear,imonth,iday,ihour,imin,isec
      common /timeinfo/  iyear,imonth,iday,ihour,imin,isec
      get_month=imonth
      return
      end
c
      function get_day()
      implicit none
      integer get_day
      integer iyear,imonth,iday,ihour,imin,isec
      common /timeinfo/  iyear,imonth,iday,ihour,imin,isec
      get_day=iday
      return
      end
c
      function get_hour()
      implicit none
      integer get_hour
      integer iyear,imonth,iday,ihour,imin,isec
      common /timeinfo/  iyear,imonth,iday,ihour,imin,isec
      get_hour=ihour
      return
      end
c
      function get_min()
      implicit none
      integer get_min
      integer iyear,imonth,iday,ihour,imin,isec
      common /timeinfo/  iyear,imonth,iday,ihour,imin,isec
      get_min=imin
      return
      end
c
      function get_sec()
      implicit none
      integer get_sec
      integer iyear,imonth,iday,ihour,imin,isec
      common /timeinfo/  iyear,imonth,iday,ihour,imin,isec
      get_sec=isec
      return
      end
c
      function get_vcp()
      implicit none
      integer get_vcp
      get_vcp=41
      return
      end
c
      function get_latitude()
      implicit none
      integer get_latitude
      get_latitude=3523000
      return
      end
c
      function get_longitude()
      implicit none
      integer get_longitude
      get_longitude=9746000
      return
      end
c
      function get_altitude()
      implicit none
      integer get_altitude
      get_altitude=381
      return
      end
