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
cdis
cdis
cdis
cdis
cdis
cdis
c
c
c
      subroutine hypso (
     1     z,                   !height meters n levels     
     1     t,                   !temp kelvin, nlevels
     1     rmd,                 !missing data flag for real numbers
     1     nn,                  !levels nn
     1     p,                   !pressure hpa, n, levels
     1     abort                !abort flag, 0= abort, 1= good
     1     )

c     NOTES:
c
c     This routine assumes dry air but could use virtual temperature
c     as input if so desired.
c
c     Even though this routine is tailored for radiometer use, 
c     the 2km cutoff in height (according to Stick Ware) is not
c     performed here but in the routine above this one.
c
c     AUTHOR: Dan Birkenheuer 12/14/2010
c     Purpose: add Radiometrics (brand) of radiometer data to LAPS
c


      implicit none

c     parameters for radiometer hypsometric eqn to determine pressure

      integer nn                !levels passed
      real z(nn)                !height (m)
      real t(nn)                !temperature (c)
      real rmd                  !missing data flag
      real p(nn)                !pressure (hPa) (variable to be "descovered"
      integer abort             !0=fail, 1= success

c     internal variables and constants

      real r                    !gas constant
      real g                    !gravity
      real c2k                  !conversion to kelvin
      real p2p                  !converson Pa to hPa
      real delz                 !delta z (m)
      real tbar                 !average temperature (k)
      integer n                 !index 

c     code

c     assign constants
      r = 287.04                !gas constant for dry air
      g = 9.80665
      c2k = 273.15
      p2p = 0.01 

c     run hypsometric equation to determine p from sfc values an hts

      if ( (nn .eq. 1) .or. (z(1) .eq. rmd) ) then ! bad p value abt
         abort = 0 
         return
      endif


      do n = 1, nn-1            !stop shy of top level

         delz = z(n+1) - z(n)   !z2-z1 (m)
         tbar = (t(n+1)+t(n)+2*c2k)/2.0 ! (k)
         p(n+1) = p(n)*p2p/(exp (g/(r*tbar)*(delz))) ! in Pa
         p(n+1) = p(n+1)/p2p    !in hPa (as returned desire)


c     this is coded for clarity.  improvement would be to remove the
c     p2p operationa in mult and later divide to improve speed.  some
c     complilers might do this.

      enddo

      abort = 1

      return
      end


