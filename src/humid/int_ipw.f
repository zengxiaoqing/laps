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
      subroutine int_ipw (x,p,data,ipw,mdf,kk)
c     subroutine int_ipw generates single column total integrated water for
c     comparison with GPS data for the variational minimizaiton technique.
c     basically the forward model for this data source.

c     input variables
      real x(3)                 !adjustment vector
      integer kk                !vertical dimension
      real p(500)                !pressure (hPa) through the column
      real data(500)             !data vector (actually sh)
      real mdf                  !missing data flag
      real ipw                  !returned computed ipw

c     internal variables
      integer k
      real d(kk)                !modified data (modified by x)
      

c     code


c     modify data vector

      do k = 1,kk
         if(data .ne. mdf) then
            if(data(k) .ge. 0.0) then

               if(p(k).ge.700.) then
                  d(k) = data(k) * x(1)
               elseif(p(k).lt. 700. .and. p(k) .ge. 500) then
                  d(k) = data(k) * x(2)
               elseif(p(k) .lt.500. .and. p(k).ge.100.) then
                  d(k) = data(k) * x(3)
               else
                  d(k) = data(k)
               endif

            else
               d(k) = data(k)
            endif
         else
            d(k) = data (k)
         endif
      enddo
      

c     integrate

      ipw = 0.0

      do k = 1,kk-1
         if(d(k) .ne. mdf) then
            if(d(k) .ge. 0.0) then
               ipw = ipw + (d(k)+d(k+1))/2.*(p(k)-p(k+1))
            endif
         endif
      enddo

c     change units of q ti g/kg

      ipw = ipw * 1.e3

c     convert from g/kg to cm

      ipw = ipw/100/9.8

      return
      end

