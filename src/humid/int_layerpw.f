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
      subroutine  int_layerpw(x,data,kstart,qs,ps,p_1d,
     1     p1,p2,p3,lpw1,lpw2,lpw3,kk,mdf)

c     subroutine to integrate the layer precipitable water from the
c     specific humidity data.  consolodating this process, since it is
c     done numerous times in the code
      
c     Birkenheuer  01/09/2001
      
      implicit none
      
c     include 'lapsparms.for'
c     include 'parmtrs.inc'
      
c     parameter variables
      
      integer kk    
      real data(kk)
      real x(3)
      integer kstart
      real qs
      real ps
      real p_1d (kk)
      real lpw1,lpw2,lpw3
      real mdf
      real p1,p2,p3
      
      
c     variables requireing dynamic allocation
      
c     NONE
      
c     internal variables
      
      integer i,j,k

c     check bounds

      if (p1.gt.ps) then ! below ground for this level
         lpw1 = mdf
         return
      endif
      
c     integrate the tpw field
      

c     integrate q for the pressure layers
      lpw1 = 0.0
      lpw2 = 0.0
      lpw3 = 0.0
      do k = kstart,kk-1
         if (p1.le.p_1d(k) .and. data(k) .ne. mdf) then
            if(p_1d(k).ge.700.) then
               lpw1 = lpw1 + abs(x(1))*( data(k)+data(k+1) )/2. 
     1              * (p_1d(k)-p_1d(k+1))
            elseif(p_1d(k).ge.500.) then
               lpw1 = lpw1 + abs(x(2))*( data(k)+data(k+1) )/2. 
     1              * (p_1d(k)-p_1d(k+1))
            elseif(p_1d(k).lt.500.) then
               lpw1 = lpw1 + abs(x(3))*( data(k)+data(k+1) )/2. 
     1              * (p_1d(k)-p_1d(k+1))
            endif
               
         elseif(p2.le.p_1d(k) .and. data(k) .ne. mdf) then
            if(p_1d(k).ge.700.) then
               lpw2 = lpw2 + abs(x(1))*( data(k)+data(k+1) )/2. 
     1              * (p_1d(k)-p_1d(k+1))
            elseif(p_1d(k).ge.500.) then
               lpw2 = lpw2 + abs(x(2))*( data(k)+data(k+1) )/2. 
     1              * (p_1d(k)-p_1d(k+1))
            elseif(p_1d(k).lt.500.) then
               lpw2 = lpw2 + abs(x(3))*( data(k)+data(k+1) )/2. 
     1              * (p_1d(k)-p_1d(k+1))
            endif
               
         elseif(p3.le.p_1d(k) .and. data(k) .ne. mdf) then
            if(p_1d(k).ge.700.) then
               lpw3 = lpw3 + abs(x(1))*( data(k)+data(k+1) )/2. 
     1              * (p_1d(k)-p_1d(k+1))
            elseif(p_1d(k).ge.500.) then
               lpw3 = lpw3 + abs(x(2))*( data(k)+data(k+1) )/2. 
     1              * (p_1d(k)-p_1d(k+1))
            elseif(p_1d(k).lt.500.) then
               lpw3 = lpw3 + abs(x(3))*( data(k)+data(k+1) )/2. 
     1              * (p_1d(k)-p_1d(k+1))
            endif

         endif
      enddo
c     change units of integrated(q) to  g/kg     
      lpw1 = lpw1 *1.e3 
      lpw2 = lpw2 *1.e3 
      lpw3 = lpw3 *1.e3 
c     add surface layer qs (already in g/kg) to first level
      lpw1 = lpw1 
     1     + ( qs + data(kstart) *1.e3 ) /2.
     1     * (ps-p_1d(kstart) )
c     comvert g/kg to mm <-- millimeters, not cm!
      lpw1 = lpw1 / 10. / 9.8
      lpw2 = lpw2 / 10. / 9.8
      lpw3 = lpw3 / 10. / 9.8
      
      return
      
      end
