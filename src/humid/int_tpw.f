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
      subroutine  int_tpw(data,kstart,qs,ps,p_3d,tpw,mdf,ii,jj,kk)

c     subroutine to integrate the total precipitable water from the
c     specific humidity data.  consolodating this process, since it is
c     done numerous times in the code
      
c     birkenheuer  feb 9 1993
      
      implicit none
      
c     include 'lapsparms.for'
c     include 'parmtrs.inc'
      
c     parameter variables
      
      integer ::  ii,jj,kk
      real, dimension(ii,jj,kk) :: data,p_3d
      integer kstart (ii,jj)
      real, dimension(ii,jj) :: qs, ps, tpw
      real :: mdf               !laps missing data flag
      real, dimension (500) :: p_vert,data_vert
      real, dimension(3) :: xdum = (/1.,1.,1./) !dummy variable for call
      
      
c     variables requiring dynamic allocation
      
c     NONE
      
c     internal variables
      
      integer i,j,k
      
c     integrate the tpw field
      
      do j = 1,jj
         do i = 1,ii
            
            tpw(i,j) = 0.0

c     fill vertical specific arrays
            do k = 1,kk
               p_vert(k) = p_3d(i,j,k)
               data_vert(k) = data(i,j,k)
            enddo               !end k

            call int_ipw (xdum,p_vert,data_vert,kstart(i,j),
     1           qs(i,j),ps(i,j),tpw(i,j),mdf,kk)

         enddo                  !i
      enddo                     !j

      return
      
      end
