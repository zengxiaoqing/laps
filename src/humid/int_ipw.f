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
      subroutine int_ipw (x,p,data,kstart,qs,ps,ipw,mdf,kk)


c     subroutine int_ipw generates single column total integrated water for
c     comparison with GPS data for the variational minimizaiton technique.
c     basically the forward model for this data source.

      implicit none

c     input variables
      real ::  x(3)             !adjustment vector
      integer ::  kk            !vertical dimension
      real ::  p(500)           !pressure (hPa) through the column
      real :: data(500)         !data vector (actually sh)
      integer :: kstart         !first k in vertical with real data 
      real :: qs                !surface q
      real :: ps                !surface pressure
      real :: mdf               !missing data flag
      real :: ipw               !returned computed ipw

c     internal variables
      integer k
      real d(kk)                !modified data (modified by x below 100hpa)
      

c     code
c     modify data vector
      do k = 1,kk
        d(k) = 0.0
      enddo

      do k = kstart,kk
         if(data(k) < 1000. ) then
            if(data(k) >= 0.0) then

               if(p(k) > 700.) then
                  d(k) = data(k) * x(1)
               elseif(p(k) < 700. .and. p(k) >  500) then
                  d(k) = data(k) * x(2)
               elseif(p(k) < 500.and. p(k) > 100. ) then
                  d(k) = data(k) * x(3)
               else
                  d(k) = data(k) ! above 100 hpa no augmentation
               endif

            else
              d(k) =    0.0  !transfers negative missing data flag
            endif
         else
            d(k) = 0.0 ! data equals missing data flag transfer sum zero 
         endif
      enddo
      

c     integrate

      ipw = 0.0

      do k = kstart,kk-1
         if( abs(d(k)) < 1000. ) then   ! skip if missing data flag
            if(d(k) >= 0.0) then !skip if negative
               ipw = ipw + (d(k)+d(k+1))/2.*(p(k)-p(k+1))
            endif
         endif
      enddo

c     change units of q ti g/kg

      ipw = ipw * 1.e3

c     add surface layer (already in g/kg)

      if (abs(qs) < 1000 .and. abs(ps) < 1500.) then ! assume sfc vals good
         ipw = ipw + (qs + d(kstart))/2. * (ps - p(kstart))
      endif

c     convert from g/kg to cm

      ipw = ipw/100./9.8

      return
      end

