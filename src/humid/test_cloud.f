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
c FORTRAN 90
      subroutine test_cloud (t1,t2,q1,q2,p,qsat, qadjust)

      implicit none

      real :: tbar, ubar, qbar, delu, delq, delt
      real :: terror, qerror
      real :: t1,t2, q1,q2, p, tref = -132.0
      real :: qsat
      real :: qadjust,uadj
      real :: temp_uncertainty = 0.005 ! 0.5% precision error in temp
      real :: q_uncertainty = 0.01 ! 1% precision error in SH

      real, external :: ssh2    !function for saturation specific 
                                ! humidity

c     tref = temperature for ice phase reference
c     p = pressure in mb (something close is all that is needed)
c     level 1 is low, level 2 is higher, level 1 is (k)
c     bar values are mean for the range
c     del values are allowed deltas
c     u = relative humidity (fraction)
c     t = temperature degK
c     terror, qerror are the approximate error in the analyzed values of
                                !these components, the natural error
                                !in these variables (precision) will 
                                !have a direct bearing on the probabily 
                                !of cloud.   
c     q = specific humidity 
c     qbarsat = the computed saturated q based on tbar
c     qadjust = the min amount that qbar can be increased to make a cloud
c     routine based on climo scheme developed by Peixoto and Oort 1997
c     refer to NOAA Tech Note for reference

      tbar = (t1+t2)/2.         ! average temperature in K
      qbar = (q1+q2)/2.         ! average specific humidity

      terror = temp_uncertainty* tbar !approximate temperature error
      
      qerror = q_uncertainty* qbar !approximate SH error
      
      delq = abs(q2-q1)         !box moisture variabiliy
      delt = abs(t2-t1)         !box  temperature variabilty
      
      ubar = q1/ qsat           !use gridpoint value here 
      
      delu = ubar * ( (delq+qerror)/qbar + 20.*(delt+terror)/tbar )
      
c     this is the accepted variation of Relhumidity in the "box"
c     immediately adjacent to a a level  we now apply this variation
c     to that level (k) being evaluated.
     
c     if this range extends RH to saturation, then the box is
c     assumed to support cloud.  If clouds are associated with this
c     box, then no adjustment to moisture is needed if the cloud
c     analysis spcifies clouds

c     test the box for cloudy conditions.

      if (ubar+delu/2. >= 1.) then ! clouds are seen as a possibility
         qadjust = 0.0  ! no need to enhance q if clouds are present
      else                      !clouds are not a possiblity compute
                                !the minimum amount of extra q needed 
                                ! to make a cloud somewhere in the box

         uadj = 1. - (ubar+delu/2.)

         qadjust = uadj * qsat

         if (qadjust < 0 ) qadjust = 0.0 !slight roundoff correction

c     this is the approximate increment in q needed to achieve saturation

      endif

      return

      end subroutine test_cloud
