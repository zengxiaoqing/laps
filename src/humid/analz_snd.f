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
      subroutine analz_snd (rh,ll_n,lat_s,lon_s,p_s,ll,nnn,nn,
     1     raob_radius,
     1     lat,lon,p_3d,
     1     RH_fields,W_field,ii,jj,kk)

c     another legacy in software

      implicit none

c     input variables

      integer ll,nnn            ! dimensions of snd arrays
      integer nn                ! number of soundings
      real rh(ll,nnn)           ! snd rh by level and nnn soundings
      integer ll_n(nnn)         ! number of levels at each sounding
      real raob_radius          ! radius of influence (meters)
      real lat_s(ll,nnn)        ! latitude of SND lvl data
      real lon_s(ll,nnn)        ! longitude of SND lvl data
      real p_s(ll,nnn)          ! pressure of SND data
      integer ii,jj,kk          ! dimensions of LAPS data fields
      real lat(ii,jj)           ! laps grid lats
      real lon(ii,jj)           ! laps grid lons
      real p_3d(ii,jj,kk)       ! laps 3D pressure field
      real RH_fields(ii,jj,kk)  ! laps analyzed RH field
      real W_field(ii,jj,kk)    ! SND horizontal weights for variational step

c     internal variables

      real rh_k(kk)             ! interpolated RH at k level
      real lat_k(kk)            ! interpolated lat at k level 
      real lon_k(kk)            ! interpolate lon at k level
      integer i,j,k,l,n         ! indexes
      real ri,rj,rk             ! real counterparts of indexes
      integer mask(ii,jj,kk)    ! mask field used for weights
      real r50                  ! 50% weight distance (weighting)
      integer istatus           ! used in function calls
      real pres(kk)             ! local pressure array for interpolation
      real points(nnn,3,kk)     ! array of points for prep_grid.f
      integer pn(kk)            ! number of points in each level of points

c     interp each sounding in three dimensions and place the result value
c     into the 3-D array.

c     first step is to interp in the horizontal, this is because the
c     pressure grid is not uniform in the horizontal and to get a good
c     vertical interpolation, i and j must be known first

c     initialzed summed arrays and variables


      do k = 1, kk              !all levels in laps

         pn(k) = 0
         rh_k(k) = 0.0
         lat_k(k) = 0.0
         lon_k(k) = 0.0
         pres(k) = 0.0
      enddo                     !k laps levels

      do n = 1,nnn
         do i = 1,3
            do k = 1,kk
               points(n,i,k) = 0.0
            enddo
         enddo
      enddo

c     code     


      do n = 1,nn  ! for all soundings
         do l = 1,ll_n(n) ! for all levels in a particular sounding

            call latlon_to_rlapsgrid(lat_s(l,n),lon_s(l,n), lat,lon, 
     1           ii,jj,  ri,rj, istatus)

            if (istatus .eq. 1) then ! continue to modify location

               i = int (ri)
               j = int (rj)
               if (i.le.0) i = 1
               if (i.gt.ii) i = ii
               if (j.le.0) j = 1
               if (j.gt.jj) j = jj

c     next step is to interp in the vertical (since pressure can be
c     horizontally dependent given the possibility of p_3d).

c     write the 3d pressure vector into pres

               do k = 1,kk
                  pres(k) = p_3d(i,j,k)
               enddo            !end k vertical variable
               
               call locate (pres,kk,p_s(l,n),k) ! returns k level

c     now i,j,k is known for the data element RH_s and can be transferred
c     to the correct 3-d location.

c     transfer data to 3-d array and 3-d vsn of points
               if (k.gt.0 .and.  k.lt.kk) then ! valid k

                  pn(k) = pn(k) +1
                  points(pn(k),1,k) = rh (l,n)
                  points(pn(k),2,k) = i
                  points(pn(k),3,k) = j

                  mask(i,j,k) = 1
                  rh_fields(i,j,k) = rh(l,n)

               elseif (k.eq.kk.and.p_s(l,n).eq.pres(k)) then 
c     ! special case (top limit)

                  pn(k) = pn(k) +1
                  points(pn(k),1,k) = rh (l,n)
                  points(pn(k),2,k) = i
                  points(pn(k),3,k) = j

                  mask(i,j,k) = 1
                  rh_fields(i,j,k) = rh(l,n)

               endif            !valid k test
            endif               !test for valid ij
         enddo                  !l vertical dimension of sounding data
      enddo                     !n sounding number

c     all data now in the 2 arrays in question.

c     next step is to distribute the data using existing routines




      do k = 1,kk

         if(pn(k).ne.0) then

            call prep_grid (ii,jj,rh_fields(1,1,k),nnn,
     1           points(1,1,k),pn(k),istatus)
            call slv_laplc (rh_fields(1,1,k),mask(1,1,k),ii,jj)
            r50 = raob_radius       ! spatial influence from namelist
            call weight_field (w_field(1,1,k), mask(1,1,k), ii,jj,
     1           r50, istatus)

         endif

      enddo                     !k for all levles

      return

      end
      
      
