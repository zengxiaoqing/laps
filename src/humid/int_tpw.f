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
cdis cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
      subroutine  int_tpw(data,kstart,qs,ps,plevel,tpw,ii,jj,kk)

c     $log: int_tpw.for,v $
c     revision 1.3  1995/09/13  21:35:23  birk
c     added disclaimer to files
c     
c     revision 1.2  1994/07/22  22:02:22  birk
c     made code compatible with lapsparms.inc
c     
c     revision 1.1  1994/04/25  15:01:49  birk
c     initial revision
c     
      
c     subroutine to integrate the total precipitable water from the
c     specific humidity data.  consolodating this process, since it is
c     done numerous times in the code
      
c     birkenheuer  feb 9 1993
      
      implicit none
      
c     include 'lapsparms.for'
c     include 'parmtrs.inc'
      
c     parameter variables
      
      integer ii,jj,kk
      real data(ii,jj,kk)
      integer kstart (ii,jj)
      real qs(ii,jj)
      real ps(ii,jj)
      real plevel (kk)
      real tpw (ii,jj)
      
      
c     variables requireing dynamic allocation
      
c     NONE
      
c     internal variables
      
      integer i,j,k
      
c     integrate the tpw field
      
      do j = 1,jj
         do i = 1,ii
c     integrate q
            tpw(i,j) = 0.0
            do k = kstart(i,j),kk-1
               tpw(i,j) = tpw(i,j) + ( data(i,j,k)+data(i,j,k+1) )/2. 
     1              * (plevel(k)-plevel(k+1))
            enddo
c     change units of q to g/kg
            tpw(i,j) = tpw(i,j) *1.e3 
c     add surface layer (already in g/kg)
            tpw(i,j) = tpw(i,j)
     1           + ( qs(i,j) + data(i,j,kstart(i,j)) *1.e3 ) /2.
     1           * (ps(i,j)-plevel(kstart(i,j)) )
c     comvert g/kg to cm
            tpw(i,j) = tpw(i,j) / 100. / 9.8
            
         enddo
      enddo
      
      return
      
      end
