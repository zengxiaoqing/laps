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


      subroutine gen_btemps (kk,temp,pres,sfc_temp,sfc_pres,q,
     1     jday,i4time,lat,lon,sat,tbest,radest)

c     this module is to be used either by analysis or model input to generate
c     synthetic btemps for model comparison to satellite imagery.


      implicit none

      integer kk                !number of vertical levels in the call 
                                !including surface.  
                                !note the organization of profiles is 1=surface
                                !kk = top
      real temp(kk)             !level temp in K
      real pres(kk)             !level pressure in hPa
      real sfc_temp, sfc_pres   !surface temp and pressure (K, hPa)
      real q (kk)               !level moisture in specific humidity g/g
      integer jday              !julian day for climo data fill
      integer i4time            !current time for solar position comp
      real lat                  !latitude (degrees) for climio fill
      real lon                  !longitude (east, degrees) for zenith ang.
      integer  sat              !goes satellite number (10,9,12,8, etc)
      character*256  c_path,t_path
                                !coefficient path and transmittance path
      real tbest (18)           !goes sounder channel brightness temps (return)
      real radest(18)           !goes sounder channel radiance (return)

c     internal variables

      real :: za,sec_za         !zenith angle degrees
      real :: sfc_emis = 0.98   !surface emissivity
      real :: sfc_refl          !surface reflectance
      real :: sec_solar         !local secant of the suns position
      real :: alt               !altitude (unused output of solar routine)
      real :: sol_sub_lat, sol_sub_lon !solar position
      real :: hour_angle         !solar hour angle
      real :: sat_sub           !satellite subpoint longitude
      real :: sza               !solar zenith angle
      character*2 csat          !character version of sat
      real zenith               !external function
      external zenith

c     satellite number

      integer goes_number
      common /sat_id/ goes_number

c     coefficient path data

      character*256 sndr_coeff, sndr_trans
      integer sndr_coeff_len, sndr_trans_len
      common /optn_coef/ sndr_coeff, sndr_trans, sndr_coeff_len,
     1     sndr_trans_len


      real :: pi, d2r


code--------------------------------------------------

      pi = acos(-1.0)
      d2r = pi/180.

      goes_number = sat  ! put in common block for ofm
      write (csat,23) sat
 23   format (i2)

c     construct the full coefficient and transmittance file names

      sndr_coeff =  'sndr_g'//csat//'_spectral_coefficients_big'
      sndr_trans =  'sndr_g'//csat//'_transmittance_coefficients_big'

      sndr_coeff_len = len_trim (sndr_coeff)
      sndr_trans_len = len_trim (sndr_trans)

c     compute the zenith angle from lat lon and also sec zen angle


      if (sat .eq. 8 .or. sat .eq.12) sat_sub = -75.
      if (sat .eq. 10 ) sat_sub = -135
      if (sat .eq. 9)  sat_sub = +155. !
      za = zenith(lat*d2r,lon*d2r,0.0,sat_sub*d2r) !zenith deg
      sec_za = 1./cos(za*d2r)

c     compute the solar zenith secant

      call solar_position (lat,lon,i4time,alt,sol_sub_lat,hour_angle)

      sol_sub_lon = lon-hour_angle

      sza = zenith(lat*d2r,lon*d2r,sol_sub_lat*d2r,sol_sub_lon*d2r)
      sza = abs(sza)
      if (sza .gt. 85.) sza  = 85.

      sec_solar = 1./cos(sza*d2r)



c     set reflectivity based on emissivity (isotropic)

      sfc_refl = (1.0-sfc_emis)/pi
      
      call ofm (kk,pres,temp,q,sfc_temp,sfc_pres,jday,lat,za,
     1     tbest,radest,sec_za,sfc_emis,sfc_refl,sec_solar)

      return
      end
