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
      subroutine analq
     1     (i4time,p_3d,ps,t,ph,td,data,cg,tpw,bias_one,
     1     kstart,qs,glat,glon,mdf,ii,jj,kk)

c     variables
c     
c     p_3d(ii,jj,kk)  the laps pressure levels
c     ps (ii,jj)  is the surface pressure field hpa
c     t (ii,jj)  surfact temp (c)
c     ph (ii,jj)  is the top of boundary layer pressure hpa
c     td (ii,jj)  is the surface dew point
c     data (ii,jj,kk) is the laps sh field (input and output here)
c     cg is the cloud grid.  if radiometer column is cloudy, no scaling is
c     performed
c     
c     other variables
c     
c     qs (ii,jj) surface q g/kg
c     qk  (ii,jj) the k index of the top of the b.l. at each point
c     tpw_point is the radiometer's total precipitable water
c     tpw is the field tpw
c     irad  is the i index of the radiometer position in the laps domain
c     jrad is the j index of the radiometer position in laps domain
c     i4time is the i4time for the call to get_radiometer data (get_rad)
c     pw is a get_rad input
c     plat is a get_rad output
c     plon is a get_rad output
c     npts is a getrad output indicating the nprofilepts in calling routine
c     istatus is normal status indicator
c     ix,jy are indexes of laps gridpoints nearest the radiometer.
      
c     preliminary computations
      
      implicit none
      
c     input variables
      
      integer i4time,ii,jj,kk
      real, dimension(ii,jj,kk) :: p_3d, data,cg
      real, dimension(ii,jj) :: ps,t,ph,td,tpw,qs,glat,glon
      real bias_one
      integer kstart (ii,jj)
      real mdf                  ! missing data flag
      
c     internal variables requiring dynamic dimensions
      
      integer qk(ii,jj)
      
c     regular internal variables
      
      real tpw_point
      real bias, bias_correction ! bias one used for validation
      integer irad,jrad
      integer i,j,k
      real cgsum                ! cloud grid sum for vertical cloud check
      real frac                 ! fraction used in linear 
                                ! interpolation BL moisture
      integer
     1     npts,                !npfilepts in calling routine
     1     istatus
      integer loop_counter
      integer qk_max
      
c     internal variables with fixed dimensions
      
      real pw(4)
      real plat(4)
      real plon(4)
      integer ix(4),jy(4)

c     --- code ---
      
c     preliminary computations
      
      loop_counter = 0
      bias_one = -500.
     
c     compute top of b.l. q from data and ph
c     compute level of bl (qk)
      
      do j = 1,jj
         do i = 1,ii
            
c     check integrity
            if (ph(i,j) .le. 0.0)then
               write(6,*) 'SEVERE error in routine analq.f '
               write(6,*) 'Aborting analq.f'
               return           ! this cannot be 0.0
            endif
            
            k = kstart (i,j)    ! k is "surface"
 10         if (p_3d(i,j,k) .gt. ph(i,j)) then !  haven't reached the bl top
               k = k+1
               if (k.gt.kk) then
                  write(6,*)  'error in analq, top of column reached'
                  return
               endif
               goto 10
            endif
            
            ph(i,j) = p_3d(i,j,k) ! this is the top of bl in reg grid
            qk(i,j) = k
            
         enddo
      enddo

      do j = 1,jj
         do i = 1,ii
            qk_max = max (qk(i,j), qk_max)
         enddo
      enddo
      write(6,*) 'QK max of field is', qk_max

      call check_nan2(ph,ii,jj,istatus)

      if(istatus.ne.1) then
         write(6,*) 'NaN detected in var:ph  routine:analq.f'
         return
      endif

      
c     ----------  mixing step
      
      do j = 1,jj
         do i = 1,ii
            

               
            do k = kstart(i,j),qk(i,j)
               if (kstart (i,j) .ge. qk(i,j) ) goto 11 ! skip loop
                  
c     new method is simple linear approximation in pressure space.
c     this method is independent of grid spacing, the old method was
c     not going to be consistent if vertical grid were to vary.  However, 
c     the old method didn't overestimate Q so badily.  So to emulate that
c     aspect, the new analysis will average the "background" with the
c     linear approximation at all levels.
                  
c                  frac = (float(k) - float(qk(i,j)) )/
c     1                 (float(kstart(i,j) ) -  float(qk(i,j)) )
c                  
c                  frac = abs (frac)
c                  
c                  data(i,j,k) = data(i,j,kstart(i,j)) * frac +
c     1                 (1. - frac) * data(i,j,qk(i,j)) 
c     1                 + data(i,j,k)
c                  
c                  data (i,j,k) = data (i,j,k) /2.

c     new experimental method 10/29/99 should eliminate problems with
c     over saturating aloft as above case does.  In effect the mixing
c     is bottom up and is limited by mixing length theory.  50% possible
c     at 50 hpa, and exponentially down from there,  using pressure as
c     independent variable.
c     0.01386294 <--> 50 hPa 
c     0.006931471 <--> 100
c     0.027725887 <-->  25 hPa
c     0.092103403 <-->  10% at 25 hpa
c     0.11982929 <-->  5% at 25 hpa


               frac = p_3d(i,j,k)-p_3d(i,j,kstart(i,j))
               frac = exp(frac*0.01386294) ! 
               data(i,j,k) = data(i,j,kstart(i,j))*frac +
     1              (1.-frac)*data(i,j,k)


               
            enddo

 11         continue            ! loop skipped  no boundary layer, no mixing
            
         enddo
      enddo
      
      call check_nan3(data,ii,jj,kk,istatus)
      if(istatus.ne.1) then
         write(6,*) 'NaN detected in var:data  routine:analq.f'
         write(6,*) 'check after frac computation'
         return
      endif
c     integrate the tpw field
         
      call int_tpw(data,kstart,qs,ps,p_3d,tpw,mdf,ii,jj,kk)

      call check_nan2(tpw,ii,jj,istatus)
      if(istatus.ne.1) then
         write(6,*) 'NaN detected in var:tpw  routine:analq.f'
         return
      endif
      
      
c     ---------- end mixing step at this point return and ignore 
c     radiometer code

      return

c--------------- this code is now defunct and not used.......

c     note that qs  has units of g/kg

c     now obtain the radiometer tpw units of cm or gm/cm**2

      call get_rad(i4time,pw,plat,plon,npts,istatus)
      
      if (istatus.eq.1) then
         
         do i = 1,npts
            plon(i) = -plon(i)
         enddo
         
         call get_laps_gp (npts,plat,plon,ix,jy,glat,glon,ii,jj)
         
c     now determine the tpw at the selected gridpoints and decide on
c     tpw_point, irad and jrad
         
c     integrate the tpw field
         
         call int_tpw(data,kstart,qs,ps,p_3d,tpw,mdf,ii,jj,kk)
         
c     determine the radiometer with the lowest bias correction and assume
c     that it is "true"
         
         do i = 1,npts
            
            if (i .eq. 1) then
               tpw_point = pw(i)
               irad = ix(i)
               jrad = jy(i)
               bias = tpw_point - tpw(irad,jrad)
               
            else
               
               if (abs(pw(i) - tpw(ix(i),jy(i)) ) .lt.abs(bias)) then
                  
                  tpw_point = pw(i)
                  irad = ix(i)
                  jrad = jy(i)
                  bias = tpw_point - tpw(irad,jrad)
                  
               endif
               
            endif
            
         enddo
         
c     now begins a loop to converge the bias correction and integrate the
c     moisture.  the surface
c     sh from td is not changed since we have faith in this value.
c     the values aloft are however changed interatively to agree with
c     an integral of tpw.  this loop also improves the tpw.
c     note: this loop does not mix water from the surface up or down
         
c     mod for version 6.0  2.2.93  db
c     check vertical column over radiometer and do not scale if cloudy
         
         cgsum = 0.
         k = 1
         
         do while (p_3d(irad,jrad,k).gt.600) ! integrate cloud in lower trop
            cgsum=cgsum+cg(irad,jrad,k)
            k = k+1
         enddo
         
         if(cgsum .gt. 0.6) then ! cloudy over radiometer
            write(6,*) 'clds over radiometer scaling bypassed'
            bias = 0.005        ! bias is assigned this value for tracking
         endif
         
c     the value of bias will cause the code to branch at this point
         
         do while (abs(bias) .gt. 0.01 )
            
c     integrate the tpw field
            
            call int_tpw (data,kstart,qs,ps,p_3d,tpw,mdf,ii,jj,kk)
            
c     compute the bias at the radiometer site
            
            bias = tpw_point - tpw(irad,jrad)
            
c     record first iteration bias (bias_one) to monitor the process
            
            if (bias_one.eq.-500.) bias_one = bias
            bias_correction = tpw_point/tpw(irad,jrad)
            print*, bias_correction, 'factor', bias , 'bias'
            
c     apply the bias change at upper levels only
            
            do k = 1,kk
               do j = 1,jj
                  do i = 1,ii
                     if(data(i,j,k) .gt. 0.0 )
     1                    data(i,j,k) = data(i,j,k) * bias_correction
                  enddo
               enddo
            enddo
            
c     repeat the integration of tpw
            
c     increment loop-counter to prevent run-away situation experienced
c     1 0/11/94 db
            
            loop_counter = loop_counter+1
            if (loop_counter.gt.15) go to 123
            
         enddo                  !  (while)
         
 123     continue               ! bailout for loop counter
         
c     at this point tpw has been integrated
c     at this point data array has been modified.
         
      else
         
c     the routine goes here if there are no radiometer data avail
c     integrate the tpw field one pass and beleive it to be good
         
         write(6,*) 'radiometer un-avail, no bias correction to tpw'
         write (6,*) 'computing ipw from laps field for later comps'
         
         call int_tpw(data,kstart,qs,ps,p_3d,tpw,mdf,ii,jj,kk)

         call check_nan2(tpw,ii,jj,istatus)
         if(istatus.ne.1) then
            write(6,*) 'NaN detected in var:tpw  routine:analq.f'
            return
         endif
      
         
      endif
      
      return
      end











