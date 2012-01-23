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
      subroutine interp_to_sfc(sfc_2d,field_3d,heights_3d,ni,nj,nk,
     &                         badflag,interp_2d)
c
c=================================================================
c
c     Routine to interpolate a 3-d field to a 2-d surface.  The
c     2-d surface can be the actual "surface" (the topography),
c     or any other 2-d surface within the grid.
c
c     Based on code in the LAPS wind analysis, Steve Albers, FSL.
c
c     Original: 04-09-99  Peter A. Stamus, NOAA/FSL
c     Changes:
c
c=================================================================
c
      real sfc_2d(ni,nj), field_3d(ni,nj,nk), heights_3d(ni,nj,nk)
      real interp_2d(ni,nj)
      
      logical ltest_vertical_grid
c
      write(6,*)' Interpolate 3-d field to 2-d surface (interp_to_sfc)'

c
c..... Interpolate from the 3-d grid to the 2-d surface at each point.
c
      do j=1,nj
      do i=1,ni
c
         if(.not. ltest_vertical_grid('SIGMA_HT'))then
            zlow = height_to_zcoord2(sfc_2d(i,j),heights_3d,ni,nj,nk,       
     &                                                  i,j,istatus)
            if(istatus .ne. 1)then
               write(6,*) 
     &             ' Error in height_to_zcoord2 in interp_to_sfc',
     &                 istatus
               write(6,*) i,j,zlow,sfc_2d(i,j),(heights_3d(i,j,k)
     1                                          ,k=1,nk)
               return
            endif

         elseif(.true.)then ! SIGMA_HT grid
            zlow = 1.       ! Surface defined as lowest SIGMA_HT level

         else               ! .FALSE. 
            zlow = rlevel_of_field(sfc_2d(i,j),heights_3d,ni,nj,nk,       
     &                                                  i,j,istatus)
            if(istatus .ne. 1)then
               write(6,*) 
     &             ' Error in rlevel_of_field in interp_to_sfc',
     &                 istatus
               write(6,*) i,j,zlow,sfc_2d(i,j),(heights_3d(i,j,k)
     1                                          ,k=1,nk)
               return
            endif

         endif
c
         klow = max(zlow, 1.)
         khigh = klow + 1
         fraclow = float(khigh) - zlow
         frachigh = 1.0 - fraclow
c
         if( field_3d(i,j,klow)  .eq. badflag .or.
     &       field_3d(i,j,khigh) .eq. badflag) then

            write(6,3333)i,j
 3333       format(' Warning: cannot interpolate to sfc at ',2i5)
            interp_2d(i,j) = badflag

         else

            interp_2d(i,j) = field_3d(i,j,klow ) * fraclow  +
     &                       field_3d(i,j,khigh) * frachigh

         endif
c
      enddo !i
      enddo !j
c
c..... That's all.
c
      return
      end
