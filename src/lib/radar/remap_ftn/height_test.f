cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
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
