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
      real function func (x)
      
c     $log: func.for,v $
c     revision 1.1  1996/08/30  20:44:48  birk
c     initial revision
c     
      
c     basically a simple cost function for the moisture function giving
c     goes-8 radiances for the moisture channels only  image case,
c     the passed parameter is the variable changed in the iteration
c     other variable (invariant for this problem) are passed via common
c     note: w_cost is same as w, but w has to be in a common for the
c     forwardmodel call reqd by wisconsin code.
      
c     the function blows up whenever w_cost is negative. therefore we test
c     to see if any w_cost element is negative if it is we assign a high
c     cost and return.
      
      implicit none
      save
      real x(3)
      integer lsfc
      real psfc
      real emiss
      real gimrad
      integer ngoes,kan(18)
      integer ngoes_cost,isnd_cost,chan_used
      real radiance (18)
      real theta,tau(40),tskin
      
      integer lsfc_cost
      real radiance_ob(18),p_cost(40),t_cost(40),
     1     ozo_cost(40),tskin_cost,psfc_cost,
     1     theta_cost,w_cost(40)
      real p(40),t(40),w(40),ozo(40)
      common/atmos/p,t,w,ozo
      common/cost_var/radiance_ob, p_cost, t_cost,ozo_cost,
     1     tskin_cost, lsfc_cost,psfc_cost,theta_cost,w_cost,
     1     ngoes_cost,isnd_cost,chan_used
      common /cost_gvap/cost_w1,cost_w2,cost_w3,cost_gvap_p,cost_weight,
     1     cost_gvap_istatus
      real cost_w1,cost_w2,cost_w3,cost_gvap_p,cost_weight
      integer cost_gvap_istatus
      
      
      integer i,j
      
      
      ngoes = ngoes_cost
      
c     set up for sounder if needed instead of imager
      
      if (isnd_cost.eq.1) then  ! sounder radiances used
         
         
         
         kan(1) = 10            ! 7.4 like ch3
         kan(2) = 8             ! 11.02 like ch4
         kan(3) = 7             ! 12.02 like ch5
         kan(4) = 11
         kan(5) = 16
         kan(6) = 6
         kan(7) = 13
         
      else                      ! imager radiances used
         
         
         kan(1) = 23            !  7.4  imager ch3
         kan(2) = 24            ! 11.02  imager ch4
         kan(3) = 25            ! 12.02  imager ch5
         
         
      endif
      
      
      emiss = .99
      
      
c     set up wisconsin variables in common
      
      do i = 1,40
         
         if(i.le.40 .and. i.gt. 35) then ! sfc to 780
            w(i) = abs(x(1)) * w_cost(i)
         elseif (i.le.35 .and. i.gt. 30) then ! 700 to 500
            w(i) = abs(x(2)) * w_cost(i)
         elseif (i.gt. 20) then ! between 475 and 100
            w(i) = abs(x(3)) * w_cost(i)
         else
            w(i) =  w_cost(i)
         endif
         
         tskin = tskin_cost
         
         
         t(i) = t_cost(i)
         p(i) = p_cost(i)
         ozo(i) = ozo_cost(i)
      enddo
      lsfc = lsfc_cost
      psfc = psfc_cost
      theta = theta_cost
      
      
c     perform forward model computation for radiance
      
      
      do j=1,chan_used
         
         call taugim(t,w,ozo,theta,ngoes,kan(j),tau)
         radiance(j) = gimrad(tau,t,tskin,kan(j),lsfc,psfc,emiss)
         
      enddo                     ! j
      
c     compute cost function
      
      func = 0.0
      
      do j = 1,chan_used
         func = func + (radiance_ob(j)-radiance(j))**2/2.
      enddo
      do j = 1,3
         func = func + ((x(j) - 1.)**2 )
      enddo

      if (cost_gvap_istatus .eq. 1) then
c     integrate q for gvap layers

c     generate modfied cost function based on these layers
         func = func
      endif

      
c     print *, (radiance(i),i=1,3),x
      
      
      return
      end
      
      
