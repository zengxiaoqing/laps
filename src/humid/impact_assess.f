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
      subroutine impact_assess (data_back, data_final, tpw1,tpw2, 
     1     ii,jj,kk, 
     1     gps_data,gps_w, p, mdf)

c     this routine compares the original and final analysis gridpoints to the
c     gps value for tpw at gps sites.  

c     the purpose of this routine is to assess the change in tpw over gps areas
c     for an assessment of improvement to the analysis.

c     parameter variables

      implicit none

      integer :: ii,jj,kk       !passed in dimensions
      real, dimension (ii,jj,kk) :: data_back, data_final
      real :: p                 !location pressure
      real, dimension (ii,jj) :: gps_data, gps_w,tpw1,tpw2
      real :: mdf               !laps missing data flag

      integer :: i,j,k

c     compute the ipw data for each place where we have gps data

      do i = 1,ii
         do j = 1,jj

c     call to ipw routine at gps point

            if(gps_w(i,j) == 1.0) then
               write(6,*) 
     1              'TASSESS',tpw1(i,j),tpw2(i,j),gps_data(i,j),i,j
            endif

c     ipw_data_back = f(i,j)
c     call int_ipw (x,p,data_back(i,j),ipw_data_back(i,j),mdf,kk)
c     call int_ipw (x,p,dadta_final(i,j),ipw_data_final(i,j),mdf,kk)
c     ipw_data_final = f(i,j)

         enddo                  !j
      enddo                     !i

c     print out comparisons of co-located points only
 

c     do n = 1,nn
c     write(6,*) ipw_data_back(point(2,n),point(3,n), ipw.....

      return

      end subroutine impact_assess
 
