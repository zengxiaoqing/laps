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
      SUBROUTINE  process_nowrad_z(imax,jmax,
     &     r_grid_ratio,
     &     image_to_dbz,
     &     image_data,
     &     r_llij_lut_ri,
     &     r_llij_lut_rj,
     &     nline,nelem,   ! input array dimensions
     &     laps_dbz,
     &     istatus)

C
c       J. Smart          Jul 1995          Orginal Subroutine for GOES 8
c.....          Changes:  02-OCT-1990       Set up for prodgen.
c.....          Changes:     SEP-1993       Add average (SA)
c
c       J. Smart          Sept 1996         Modified this satellite processing code
c                                           to do similar function with the wsi (nowrad)
c      				    reflectivity data.
c
      implicit none

      integer max_elem,max_line
      integer imax,jmax
      parameter (max_elem = 20)
      parameter (max_line = 20)
      integer   nline,nelem
      real    r_grid_ratio
      real    r_missing_data
      real    ref_base

      integer image_data(nelem,nline)
      integer image_to_dbz(0:15)

      real laps_dbz(imax,jmax)
      real wsi_dbz_data(nelem,nline)
      real r_llij_lut_ri(imax,jmax)
      real r_llij_lut_rj(imax,jmax)

      real line_mx,line_mn,elem_mx,elem_mn
      real z_array(max_elem*max_line)
      real zmax, zmin, zmean
      real rdbz
      real pixsum 
      real result
      real bad_data_flag_dbz

      integer i,j,ii,jj
      integer istart,jstart
      integer iend,jend
      integer npix, nwarm
      integer maxpix
      integer ipix
      integer istatus
      integer qcstatus
      integer fcount,icnt_out
      integer bad_data_flag_cnt

      istatus = -1
      bad_data_flag_cnt = 16
      bad_data_flag_dbz = 100
      qcstatus=0
      fcount=0

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)then
         print*,'Error getting r_missing_data'
         goto 1000
      endif
      call get_ref_base(ref_base,istatus)
      if(istatus.ne.1)then
         print*,'Error getting ref_base'
         goto 1000
      endif
c
c initialize laps dbz field to ref_base
c
      do j=1,jmax
      do i=1,imax
         laps_dbz(i,j)=ref_base
      enddo
      enddo
c
      if(r_grid_ratio .le. 0.5)then

         write(6,*)'Grid ratio .le. 0.5' !0.75'
         write(6,*)'Use pixel avg to get dbz'

         icnt_out = 0

         DO J=1,JMAX
         DO I=1,IMAX

            if(r_llij_lut_ri(i,j).ne.r_missing_data.and.
     &         r_llij_lut_rj(i,j).ne.r_missing_data)then

c
c line/elem are floating point i/j positions in ISPAN grid for input lat/lon
c also, use the lat/lon to real i/j look up table (r_llij_lut) to map out points
c needed for satellite pixels.
c****************************************************************************
               elem_mx = r_llij_lut_ri(i,j)+((1./r_grid_ratio) * 0.5)
               elem_mn = r_llij_lut_ri(i,j)-((1./r_grid_ratio) * 0.5)
               line_mx = r_llij_lut_rj(i,j)+((1./r_grid_ratio) * 0.5)
               line_mn = r_llij_lut_rj(i,j)-((1./r_grid_ratio) * 0.5)
               jstart = nint(line_mn+0.5)
               jend   = int(line_mx)
               istart = nint(elem_mn+0.5)
               iend   = int(elem_mx)

               if(istart.le.0 .or. jstart.le.0 .or.
     &            iend.gt.nelem .or. jend.gt.nline)then
                  icnt_out=icnt_out+1
c                 write(*,*)'insufficient data for lat/lon sector'
c                 write(*,1020)i,j
c1020             format(1x,'LAPS grid (i,j) = ',i3,1x,i3)
c                 write(6,1021)elem_mx,elem_mn,line_mx,line_mn
c1021             format(1x,'elem mx/mn  line mx/mn ',4f7.1)
               else
c
c **** FIND PIXELS AROUND GRID POINT
c
                  zmax=-1.e15
                  zmin=1.e15
                  npix = 0
                  pixsum = 0.
                  maxpix = 0

                  DO JJ=jstart,jend
                  DO II=istart,iend
                        
                     if(image_data(II,JJ).lt.bad_data_flag_cnt.and.
     &                  image_data(II,JJ).ge.0)then
                        npix=npix+1
                        z_array(npix)=
     +                       image_to_dbz(image_data(II,JJ))
                     endif
                           
                  enddo
                  enddo
 
                  if(npix.gt.1)then
c
c...  this section finds the highest, lowest, and mean pixel.
c
                     do ii=1,npix
                        rdbz  = z_array(ii)
                        pixsum = pixsum + rdbz
                        if(rdbz.gt.zmax)zmax=rdbz
                        if(rdbz.lt.zmin)zmin=rdbz
                     enddo

                  elseif(npix.eq.1)then
                     rdbz = z_array(npix)
                  else   
                     rdbz=ref_base
                     fcount=fcount+1
                  endif

                  if(npix .ge. maxpix)maxpix=npix

                  if(npix .gt. 1)then
                     if(pixsum.gt.0)then
                        laps_dbz(i,j) = pixsum / float(npix)
                     else
                        laps_dbz(i,j) = ref_base
                     endif

                  elseif(rdbz.gt.0.0)then

                     laps_dbz(i,j) = rdbz

                  else

                     laps_dbz(i,j) = ref_base

                  endif      ! npix .gt. 1

               end if        ! Enough data for num_lines .gt. 0

cccd           if(i .eq. i/10*10 .and. j .eq. j/10*10)then
cccd              write(6,5555)i,j,wm,wc,npix,nwarm,sc(i,j)
cccd5555         format(1x,2i4,2f10.2,2i5,f10.2)
cccd           endif

            else !laps grid point outside wsi nowrad domain

               laps_dbz(i,j) = r_missing_data

            endif            !r_llij's = r_missing_data

         enddo
         enddo

         write(6,*)'Max num WSI pix for avg: ',maxpix
         write(6,*)'Number of LAPS gridpoints missing',
     &              fcount
         write(6,*)'Number of grid points without data for ',
     &'this domain: ',icnt_out

      else
c       ---------------------------------------
c input image resolution is large relative to output grid spacing
c this section uses bilinear interpolation to map
c the four surrounding input pixels to the output grid.
c
         write(6,*)'Image res .ge. output grid spacing'
         write(6,*)'Using bilinear interp to get dBZ '

         do j=1,nline
            do i=1,nelem
               if(image_data(i,j).lt.bad_data_flag_cnt.and.
     &              image_data(i,j).ge.0)then
                  wsi_dbz_data(i,j)=image_to_dbz(image_data(i,j))
               else
                  wsi_dbz_data(i,j)=r_missing_data
               endif
            enddo
         enddo

         DO J=1,JMAX
         DO I=1,IMAX

            if(r_llij_lut_ri(i,j).ne.r_missing_data.or.
     &         r_llij_lut_rj(i,j).ne.r_missing_data)then
c
c bilinear_interp_extrap checks on boundary conditions and
c uses ref_base if out of bounds. NO; out-of-bounds now = r_missing_data
c
               call  bilinear_interp_extrap(
     &               r_llij_lut_ri(i,j),
     &               r_llij_lut_rj(i,j),
     &               nelem,nline,wsi_dbz_data,
     &               result,istatus)

               if(result .lt. bad_data_flag_dbz .and.
     &            result .gt. 0.0)then

                  laps_dbz(i,j) = result

               endif

            else

               laps_dbz(i,j) = r_missing_data

            endif

         enddo
         enddo

      end if                    ! r_image_ratio

      istatus = 1
c
1000  return
      end



