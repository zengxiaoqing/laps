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
      program test
      implicit none
      real elev,slantkm,slant_range
      real height,range
      real height_grid,r_range
      real rheight_radar
c
c     Doviak Zrnic stuff
c
      real re,forthre,fthrsqd
      parameter (re=6371 000.,
     :           forthre=(4.*re/3.),
     :           fthrsqd=(forthre*forthre))
c
c     Albers' stuff
c
      real rpd,mpd,radius_earth,radius_earth_8_thirds
      real hor_dist,curvature,height_factor
  10  print *, '  Enter elev (deg), slant range (km): '
      read(5,*) elev,slantkm
      if(slantkm.eq.0.) STOP
      slant_range=1000.*slantkm
      rheight_radar=0.
c
c     Doviak Zrnic Calculation
c
      height=sqrt(slant_range*slant_range + fthrsqd +
     :           2.*forthre*slant_range*sind(elev)) -
     :           forthre
c
      range=forthre*asin((slant_range*cosd(elev)) /
     :                   (forthre + height))
      range=0.001*range
c
      height=height+rheight_radar
c
      print *, ' Height, Range (D & Z ): ',height,range
c
c     Albers' Calculation
c
       rpd = 3.141592653589/180.
       mpd = 111194.
       radius_earth = 6371.e3
       radius_earth_8_thirds = 6371.e3 * 2.6666666

       hor_dist = slant_range * cosd(elev)

       curvature = hor_dist **2 / radius_earth_8_thirds
       height_grid =
     1  slant_range * sind(elev) + curvature + rheight_radar

       height_factor =
     1   (radius_earth + 0.5 * (rheight_radar + height_grid))
     1  /radius_earth

      r_range = hor_dist / height_factor * .001

      print *, ' Height, Range (Albers): ',height_grid,r_range

      GO TO 10
      END
