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

      subroutine lsin (i4time,p_3d,sfc_data,lt1dat,data,cg,tpw,
     1     bias_one,t,td,
     1     kstart,qs,p,mdf,ii,jj,kk,istatus)

c     this routine is the laps surface interface for water vapor
c     its function is to get the relevant boundary layer moisture
c     and insert it properly into the data array.
c     one other function of this routine is to establish the bounds
c     of the bottom of the column
      
c     updated analysis september 21 1993 -- modified the treatment of the
c     surface level.  moved the surface moisture to the immediate laps level
c     above the surface.   this was done because workstation plotting
c     software would not recognize the extra moisture below ground for
c     asthetic reasons.  also the workstation total precipitable water comp-
c     utation did not match the one made in this code for that saem reason
c     
c     the new arrangement of code does the following.
      
c     basically we prepare all for a call to int_tpw
c     we leave the boundary layer mixing process to analq along with
c     radiometer adjustment.
      
c     1 ) puts qs (surface q) at the p_3d above ps (surface pressure)
c     2) does not write any data below ground
c     3) maintains an integral from the true surface ps for tpw
c     this might cause some difference between this output and
c     the workstation integrated soundings, but not the error of
c     putting the data below ground.

      use module_sfc_structure

      implicit none

c     input variables


      integer i4time,istatus,ii,jj,kk
      type(lbsi), intent(in out), dimension (ii,jj) :: sfc_data
      real p_3d (ii,jj,kk)
      real lt1dat (ii,jj,kk)    ! laps 3-d temp field
      real data (ii,jj,kk)
      real cg (ii,jj,kk)
      real tpw (ii,jj)
      real bias_one
      integer kstart (ii,jj)
      real qs (ii,jj)
      real p (ii,jj)            !surface pressure (topo level)
      real mdf
      
c     internal variables with lapsparms.inc dependence
      
      real
     1     t(ii,jj),            !surface temperature k
     1     pu(ii,jj),           !pressure if top of boundary layer
     1     td(ii,jj),           !dew point temperature of surf. k -> c
     1     blsh(ii,jj)          !boundary layer specific humidity
      
c     normal internal variables
      
      integer
     1     i,j,k
      
c     constants
      save r, bad, g
      real
     1     r,                   !the gas constant for dry air
     1     bad,                 !bad data flag
     1     g                    !the acceleration of gravity
      
      data r /287.04/
      data g /9.80665/
      data bad/-1e30/
      
c     function names
      
      real ssh2                 ! type the funciton ssh2
      
c     special notes:   td will be converted to c for call to ssh
c     p will be converted to mb for comparison to vert coord
      
c-------------------------------code-----------------------------
      
c     get required field variables
      

c     fill sfc data structure with temperature (K)
c     t is now passed into routine from above

      sfc_data%sfc_temp = t



      

c     p is now passed into routine from above
c     convert p to mb

      p = p*0.01
      sfc_data%sfc_pres = p




c     td is now passed into routine from above










c     compute boundry mixing
      
      call ghbry (i4time,p_3d,p,t,lt1dat,pu,ii,jj,kk,
     1     istatus)
      if(istatus.ne.1) return

      call check_nan2(pu,ii,jj,istatus)

      if(istatus.ne.1) then
         write(6,*) 'NaN detected in var:pu  routine:lsin.f'
         return
      endif
      
      istatus = 0               ! begin with bad istatus
      
c     convert td and t to c then compute surface specific h.

      td = td - 273.15
      t  =  t - 273.15
      
      do j = 1,jj
         do i = 1,ii
            qs (i,j) = ssh2 (p(i,j),t(i,j),
     1           td(i,j),0.0)   ! qs is gm/kg
         enddo
      enddo

      call check_nan2 (qs,ii,jj,istatus)

      if(istatus.ne.1) then
         write(6,*) 'NaN detected in var:qs  routine:lsin.f'
         return
      endif

      blsh = qs * 1.e-3 ! blsh is gm/gm
      
      
c     write surface data in at the bottom of the column
c     define kstart (k index of bottom of the column)
      
      do j = 1,jj
         do i = 1,ii
            do k = 1,kk
               
               if (p(i,j).lt. p_3d(i,j,k)) then ! p_3d is underground
                  data(i,j,k) = bad
                  cg(i,j,k)   = 0.0 ! no clouds under ground
               else
                  data(i,j,k) = blsh(i,j) ! assign boundary layer sh (gm/gm)
c     to the bottom level of the column
c     define kstart
                  kstart(i,j) = k
c     jump out of loop
                  go to 2001
               endif
               
            enddo
 2001       continue
         enddo
      enddo

      call check_nan3 (data,ii,jj,kk,istatus)

      if(istatus.ne.1) then
         write(6,*) 'NaN detected in var:data  routine:lsin.f'
         return
      endif

      call check_nan3 (cg,ii,jj,kk,istatus)

      if(istatus.ne.1) then
         write(6,*) 'NaN detected in var:cg  routine:lsin.f'
         write(6,*) 'detected after boundary layer adjust'
         return
      endif

c     compute the total precipitable water and bias correct total 3-d field to
c     radiometer data
      
      print*, 'call routine analq'
      call analq(i4time,p_3d,p,t,pu,td,data,cg,tpw,bias_one,kstart,
     1     qs,sfc_data%lat,sfc_data%lon,mdf,ii,jj,kk)
      
      print*, 'done with routine analq'
      
      istatus = 1
      
      return
      
      end
