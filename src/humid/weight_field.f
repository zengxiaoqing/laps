cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps
cdis
cdis    this software and its documentation are in the public domain and
cdis    are furnished "as is."  the united states government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  they assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  if significant modifications or enhancements
cdis    are made to this software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
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
