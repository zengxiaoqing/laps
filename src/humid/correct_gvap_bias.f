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
      subroutine correct_gvap_bias (wv,wv1,wv2,wv3,n)
c
c     routine to use the bias information gleened comparing the IHOP 3x3 
c     gvap/gps comparison statistics with the goes data.  Note the correction 
c     is applied first to the total water as this is what is compared to GPS.
c     then the difference factor is computed (scale amount) and this is then 
c     applied to the layer pw values that sum to wv

c     note that this routine modifies the values of WV such that the returned
c     data are not the same as those input thus minimizing the impact on
c     other parts of the main routine.

      implicit none
      integer n
      real, dimension (n) :: wv, wv1,wv2, wv3
c     internal
      real a,b
      parameter (a=0.83885946)
      parameter (b=0.19417955)
      real scale, adjusted_wv, sum_scale, counter


      integer i !loop veriable

c     code

      sum_scale = 0.0
      counter = 0.0

      
c     start a large loop recall that correction in in cm, data in mm
      do i = 1,n                ! for all values of wv
         if (wv(i)/10. .gt. b/(1.-a)) then ! adjust it is over 1.17cm
            adjusted_wv =  a*wv(i)/10.+b ! adjusted in cm
            scale = adjusted_wv*10. / wv(i) ! scale from mm/mm
            wv(i) = wv(i)*scale
            wv1(i) = wv1(i)*scale
            wv2(i) = wv2(i)*scale
            wv3(i) = wv3(i)*scale
            sum_scale = sum_scale + scale
            counter= counter+1.0
         endif
      enddo                     !all data

      write (6,*) 'Bias Corrected GOES IPW data, factor ~ ',
     1     sum_scale/counter

      return
      end

