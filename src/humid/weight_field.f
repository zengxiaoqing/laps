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


      subroutine weight_field (field, mask, ii,jj, r50, istatus)

      implicit none
      integer ii,jj             !field dimension
      real field (ii,jj)        !weight field to be determined
      integer mask (ii,jj)      !mask (indicating where stations are)
      real r50                  !range where 50% wieght is applied
      integer istatus           !istatus

      real template (0:ii,0:jj)     !internal template array
      real grid_spacing         !spacing between gridpoints
      integer i,j               !indexes (generic)
      integer ix,jx             !indexes for template assignment
      real frac                 !fraction in template exponent

c     determine TEMPLATE array (uses r50, gridspacing) exponential function

c     zero out TEMPLATE it may have garbage values in it

      do i = 0, ii
         do j = 0,jj
            template(i,j) = 0.0
         enddo
      enddo

c     fill spacing between gridpoints

      call get_grid_spacing(grid_spacing, istatus)
      if (istatus.ne.1) then 
         write(6,*) 'grid spacing not available'
         write(6,*) 'in routine weight_field.f, aborting now'
         return
      endif


c     using r50 and the spacing, determine function for filling template

      r50 = r50 / grid_spacing

      frac = alog (0.5) /r50

c     note frac is negative

      do j  = 0, jj
         do i = 0, ii

            template (i,j) = exp(sqrt(float(i**2+j**2))*frac)

         enddo
      enddo

c     apply template array to weight grid for each data point in MASK

      do j = 1, jj
         do i = 1,ii
            if(mask (i,j) .eq. 1) then ! at a point to fill

               do jx = 1,jj
                  do ix = 1,ii

                     if(field(ix,jx).lt.template(abs(ix-i),abs(jx-j)))
     1                    field(ix,jx) = template(abs(ix-i),abs(jx-j))
                  enddo
               enddo
            endif
         enddo
      enddo

      istatus = 1

      return
      end
