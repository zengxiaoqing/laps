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
      subroutine make_optran_layers (a,b,c,aa,bb,cc,p,pp,kk,istatus)
c
c     Author: D. Birkenheuer
c     Date: 9/4/2002
c     Project: JCDSA
c     Function: 
c     This routine is designed to generate density weighted 
c     layer averages for use in OPTRAN90, the FORTRAN90 
c     version of OPTRAN initially written by Tom Kleespies 
c     in F77 and later converted to F90 by Paul vanDelst.  
c     OPTRAN90 requires layer-density weighted values 
c     for temperature, water vapor and ozone mixing ratios along 
c     with pressure.  These are accomplished by calling 
c     this routine with these independent variables as the first 
c     4 arguments, and the p variable is pressure.  The 
c     routine is set to process profiles independently thus allowing 
c     for the LAPS domain to possess vertical profiles 
c     with different pressure levels.


C     PARAMETER VARIABLES
      integer, intent (in) :: kk
      integer, intent (in out) ::  istatus
      real, dimension (0:kk), intent (in) :: a,b,c,p ! a-c are arbitrary
                                !  but abc are for temp, h2o and ozone
                                !                      p is pressure
      real, dimension (kk), intent (out) :: aa,bb,cc,pp ! arbitrary outputs

C     INTERNAL VARIABLES
      integer :: k              ! index to KK

      istatus = 0               ! set fault

c     compute layer pressure (density weighted) per Van Delst
      do k = 1, kk
         pp(k) = ( p(k-1)-p(k))/log( p(k-1)/p(k) )
      enddo                     !k

c     in certain cases where the surface is very close to the lowest level
c     a roundoff error will cause the pp(kk) value to be larger than p(kk)! 
c     this in turn causes an out-of bounds extrapolation in the following code
c     plus may damage the integration in the optran code.  To get around this,
c     here a test is done to compare the values and if the pp(kk) is larger
c     than p(kk), and average approximation is used for pp(kk)

      if (pp(kk) .gt. p(kk) .or. pp(kk) .lt. p(kk-1)) then ! error condition
c         write (6,*) 'ppkk too close', pp(kk), p(kk-1), p(kk)
         pp(kk) = (p(kk)+p(kk-1)) /2.
c         write (6,*) 'ppkk changed to: ', pp(kk), p(kk-1), p(kk)
      endif

c     use the above pressure to derive the temp and other parameter values
c     at the density weighted pressure by linear interpolation in log p
c     space

      do k = 1,kk
         call interp (log(pp(k)), log(p(k-1)), log (p(k)), 
     1        a(k-1), a(k), aa(k) )
         call interp (log(pp(k)), log(p(k-1)), log (p(k)), 
     1        b(k-1), b(k), bb(k) )
         call interp (log(pp(k)), log(p(k-1)), log (p(k)), 
     1        c(k-1), c(k), cc(k) )
      enddo ! k

c     remove numerical artifacts
      where (aa < 0.0) aa = 0.0
      where (bb < 0.0) bb = 0.0
      where (cc < 0.0) cc = 0.0

c     test for fault condition -- may now be superfluous

      if ( any ( aa < 0.0 ) .or.
     1     any ( bb < 0.0 ) .or.
     1     any ( cc < 0.0 ) )  then
         write (6,*) 'error in make_optran_layers.f, out of range error'
         stop
         istatus = 0 
      endif


      istatus = 1               ! set success
      return

      end subroutine make_optran_layers
