cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
       SUBROUTINE rrvdat2laps_v01(imax,jmax,
     &                  r_grid_ratio,
     &                  r_missing_data,
     &                  image_v01,
     &                  r_llij_lut_ri,
     &                  r_llij_lut_rj,
     &                  line_dim,elem_dim, ! input array dimensions
     &                  sa,sc,st,
     &                  istatus)

c.....  This routine can be called for any GOES Sounder channel data.
c.....  This is McAlbers version returing warmest and coldest pixels
c.....  The warmest go into ST, the "extremum" to SC, the mean to SA
c       SC is the result of an edge enhancing filter that tries to filter
c       out inbetween pixels in favor of either the warmest or the coldest.
C
c       J. Smart          Jul 1995          Orginal Subroutine for GOES 8
c.....          Changes:  02-OCT-1990       Set up for prodgen.
c.....          Changes:     SEP-1993       Add average (SA)
c.....  J. Smart          Dec 1996          Modified the satdat2laps_ir for sounder data
c.......J. Smart          Jul 1997          Allows for partial domain coverage mapping. Checks
c                                           for 0 values in lookuptables.
c
c.......J. Smart          Nov 1997          Modification for rrv to v01 remapping.
c
	implicit none

	integer*4 max_elem,max_line
        integer*4 imax,jmax
	parameter (max_elem = 15)
	parameter (max_line = 15)
	integer*4 line_dim,elem_dim
        real*4    r_grid_ratio

	real image_v01(elem_dim,line_dim)
        real*4 st(imax,jmax)
        real*4 sc(imax,jmax)
	real*4 sa(imax,jmax)
        real*4 r_llij_lut_ri(imax,jmax)
        real*4 r_llij_lut_rj(imax,jmax)

	real*4 line_mx,line_mn,elem_mx,elem_mn
	real*4 t_array(max_elem*max_line)
	real*4 wm, wc, btemp, tmean
        real*4 frac
        real*4 fraci,fracj
        real*4 f_inv_ratio_h
        real*4 f_inv_ratio_l
        real*4 pixsum 
        real*4 r_missing_data
        real*4 result
        real*4 percent_msng
        real*4 percent_zero
        real*4 percent_gt_zero
        real*4 percent_lt_zero

        integer*4 icnt_msng
        integer*4 icnt_zero
        integer*4 icnt_gt_zero
        integer*4 icnt_lt_zero
        integer*4 itot
        integer*4 i,j,ii,jj
        integer*4 istart,jstart
        integer*4 iend,jend
	integer*4 npix, nwarm
        integer*4 maxpix
        integer*4 ipix
        integer*4 istatus
        integer*4 qcstatus
        integer*4 fcount
        integer*4 ipt,jpt
        integer*4 nipt,njpt
        integer*4 icnt
        integer*4 iwrite
        save
        data iwrite/1/

        CALL ZERO(ST,IMAX,JMAX)
        istatus = -1
        qcstatus=0
        fcount=0
c
c       write(6,*)'   I   J   WarmPix  ColdPix  NPix Nwarm  CldTemp'
c
c The "10" loop represents input image resolution < output grid resolution such
c that there are enough pixels from the input image to get a representative
c mean value for the remapped output grid value
c
c rrv data is 5km resolution and it is sparse data so use this switch. 
c
        itot=0
        icnt_msng=0
        icnt_zero=0
        icnt_gt_zero=0
        icnt_lt_zero=0
        fcount=0
        icnt=0

        if(r_grid_ratio .lt. 1.0)then  !0.75)then

          write(6,*)'Grid ratio .lt. 1.0'   !0.75'
          write(6,*)'Use pixel avg to get dbz from rrv'
          maxpix=0
          DO 10 J=1,JMAX
          DO 10 I=1,IMAX

          IF(ST(I,J).NE.0.) GO TO 10
c
c line/elem are floating point i/j positions in ISPAN grid for input lat/lon
c also, use the lat/lon to real i/j look up table (r_llij_lut) to map out points
c needed for satellite pixels.
c****************************************************************************
             elem_mx = r_llij_lut_ri(i,j) + ((1./r_grid_ratio) * 0.5)
             elem_mn = r_llij_lut_ri(i,j) - ((1./r_grid_ratio) * 0.5)
             line_mx = r_llij_lut_rj(i,j) + ((1./r_grid_ratio) * 0.5)
             line_mn = r_llij_lut_rj(i,j) - ((1./r_grid_ratio) * 0.5)
             jstart = nint(line_mn+0.5)
             jend   = int(line_mx)
             istart = nint(elem_mn+0.5)
             iend   = int(elem_mx)

             if(istart.le.0 .or. jstart.le.0 .or.
     &iend.gt.elem_dim .or. jend.gt.line_dim)then

c               write(*,*)'insufficient data for lat/lon sector'
c               write(*,1020)i,j
1020	        format(1x,'LAPS grid (i,j) = ',i3,1x,i3)
c               write(6,1021)elem_mx,elem_mn,line_mx,line_mn
1021            format(1x,'elem mx/mn  line mx/mn ',4f7.1)

                icnt = icnt + 1

             else
c
                npix = 0
                pixsum = 0.
                itot=itot+1

                DO 3 JJ=jstart,jend
                DO 3 II=istart,iend

                 if(image_v01(II,JJ).ne.r_missing_data)then
                    npix=npix+1
                    t_array(npix) = image_v01(II,JJ)
                 endif
    3           continue  
c
                if(npix.gt.0)then
                   Do ii=1,npix
                      pixsum = pixsum + t_array(ii)
                   enddo
                else   
                   pixsum=r_missing_data
                   fcount=fcount+1
                endif

                if(npix .ge. maxpix)maxpix=npix

                if(npix .gt. 0)then
                   sa(i,j) = pixsum / float(npix)
                else
                   sa(i,j)=pixsum
                   sc(i,j)=pixsum
                   sT(i,j)=pixsum
                endif ! npix .gt. 1

                if(pixsum.eq.r_missing_data)then
                   icnt_msng=icnt_msng+1
                elseif(pixsum.eq.0.0)then
                   icnt_zero=icnt_zero+1
                elseif(pixsum.gt.0.0)then
                   icnt_gt_zero=icnt_gt_zero+1
                elseif(pixsum.lt.0.0)then
                   icnt_lt_zero=icnt_lt_zero+1
                endif

             end if  ! Enough data for num_lines .gt. 0

cd           if(i .eq. i/10*10 .and. j .eq. j/10*10)then
cd              write(6,5555)i,j,wm,wc,npix,nwarm,sc(i,j)
cd5555         format(1x,2i4,2f10.2,2i5,f10.2)
cd           endif

   10     CONTINUE ! I,J

          write(6,*)'Number of pixels with no radar ',icnt
          write(6,*)'Max num sndr pix for avg: ',maxpix
          write(6,*)'Number of gridpts missing',fcount

c ---------------------------------------------------------------------
        elseif(r_grid_ratio .le. 1.1)then
c
c input image resolution is large relative to output grid spacing
c this section uses bilinear interpolation to map
c the four surrounding input pixels to the output grid.
c
          if(iwrite.eq.1)then
             write(6,*)'Image res .ge. output grid spacing'
             write(6,*)'Using bilinear interp for sndr Rad'
             iwrite=0
          endif

          DO 20 J=1,JMAX
          DO 20 I=1,IMAX
c            if(i_found_it .ne. 1)then
c              write(6,*)'Enter i and j point'
c              read(5,*)ipt,jpt
c              i_found_it = 1
c            end if
            IF(ST(I,J).NE.0.) GO TO 20
c
c line/elem are floating point i/j positions in ISPAN grid for input lat/lon
c
c bilinear_interp_extrap checks on boundary conditions and
c uses r_missing_data if out of bounds.
c
            if( (r_llij_lut_ri(i,j).gt.0.0) .and.
     &          (r_llij_lut_rj(i,j).gt.0.0) .and.
     &          (r_llij_lut_ri(i,j).le.elem_dim) .and.
     &          (r_llij_lut_rj(i,j).le.line_dim) )then

c               call bilinear_interp_extrap(r_llij_lut_ri(i,j),
c    &                                      r_llij_lut_rj(i,j),
c    &                                      elem_dim,line_dim,
c    &                                      image_v01,
c    &                                      result,
c    &                                      istatus)


               call bilinear_laps(r_llij_lut_ri(i,j),
     &                            r_llij_lut_rj(i,j),
     &                            elem_dim,line_dim,
     &                            image_v01,
     &                            result,
     &                            istatus)

               sa(i,j)=result
               sc(i,j)=sa(i,j)
               sT(i,j)=sa(i,j)
               itot=itot+1

               if(result.eq.r_missing_data)then
                  icnt_msng=icnt_msng+1
               elseif(result.eq.0)then
                  icnt_zero=icnt_zero+1
               elseif(result.gt.0.0)then
                  icnt_gt_zero=icnt_gt_zero+1
               elseif(result.lt.0.0)then
                  icnt_lt_zero=icnt_lt_zero+1
               endif

            endif

c            if(i.eq.ipt .and. j.eq.jpt)then
c               write(29,39)ipt,jpt
c               write(29,49)image_v01(ii,jj)
c               i_found_it = -1
c 39          format(1x,'ipt = ',i3,2x,'jpt = ',i3)
c 49      format(1x,'ir_count = ',1x,f6.1,/,'--------------------',//)
c            end if

   20     CONTINUE ! I,J
c          close(29)
c --------------------------------------------------------------------------
        else       !this for large resolutions in which we take nearest pixel
                   !the nearest pixel method is constrained by using the inverse
c                   of r_grid_ratio which may lead to more points = r_missing_data.

           write(6,*)'r_grid_ratio > 1.1  Use nearest point'
           f_inv_ratio_l=1.0/r_grid_ratio
           f_inv_ratio_h=1.0-f_inv_ratio_l
           do j=1,jmax
           do i=1,imax
              if(r_llij_lut_ri(i,j).ne.r_missing_data.and.
     &           r_llij_lut_rj(i,j).ne.r_missing_data)then
c
c constrain the points that get to the laps grid using r_grid_ratio
c
                 nipt=int(r_llij_lut_ri(i,j))
                 njpt=int(r_llij_lut_rj(i,j))
                 fraci=r_llij_lut_ri(i,j)-nipt
                 fracj=r_llij_lut_rj(i,j)-njpt
c                if((fraci.le.f_inv_ratio_l.and.fracj.le.f_inv_ratio_l)
c    &.or.          (fraci.ge.f_inv_ratio_h.and.fracj.ge.f_inv_ratio_h)
c    &.or.          (fraci.le.f_inv_ratio_l.and.fracj.ge.f_inv_ratio_h)
c    &.or.          (fraci.ge.f_inv_ratio_h.and.fracj.le.f_inv_ratio_l))
c    &then
                    ipt=nint(r_llij_lut_ri(i,j))
                    jpt=nint(r_llij_lut_rj(i,j))
                    if((ipt.gt.0.and.jpt.gt.0)        .and.
     &              (ipt.le.elem_dim.and.jpt.le.line_dim) )then
                       sa(i,j)=image_v01(ipt,jpt)
                    else
                       sa(i,j)=r_missing_data
                    endif
c                else
c                   sa(i,j)=r_missing_data
c                endif
              else
                 sa(i,j)=r_missing_data
              endif
           enddo
           enddo

        end if     ! r_grid_ratio

c        WRITE(6,1234) IB,I4VTIME,ICT
c1234       FORMAT(1X,'BAND ',I4,' COUNT FOR I4TIME ',I10,' IS ',I8)

        percent_msng=float(icnt_msng)/float(itot)
        percent_zero=float(icnt_zero)/float(itot)
        percent_gt_zero=float(icnt_gt_zero)/float(itot)
        percent_lt_zero=float(icnt_lt_zero)/float(itot)

        if(icnt_msng .ne. itot)then
           write(*,*)'------------------------------'
           write(*,*)'Total potential useable points: ',itot
           write(*,*)'N missing/% ',icnt_msng,percent_msng
           write(*,*)'N =  zero/% ',icnt_zero,percent_zero
           write(*,*)'N > zero/% ',icnt_gt_zero,percent_gt_zero
           write(*,*)'N < zero/% ',icnt_lt_zero,percent_lt_zero
           write(*,*)'------------------------------'
           write(*,*)
        else
           write(*,*)'All points = missing'
        endif

        istatus = 1
c
        return
        end
