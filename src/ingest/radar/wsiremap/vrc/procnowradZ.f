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
     &                  r_grid_ratio,
     &                  image_to_dbz,
     &                  image_data,
     &                  r_llij_lut_ri,
     &                  r_llij_lut_rj,
     &                  line_dim,elem_dim, ! input array dimensions
     &                  laps_dbz,
     &                  istatus)

C
c       J. Smart          Jul 1995          Orginal Subroutine for GOES 8
c.....          Changes:  02-OCT-1990       Set up for prodgen.
c.....          Changes:     SEP-1993       Add average (SA)
c
c       J. Smart          Sept 1996         Modified this satellite processing code
c                                           to do similar function with the wsi (nowrad)
c					    reflectivity data.
c
	implicit none

        include 'lapsparms.cmn'
        include 'vrc.inc'

	integer*4 max_elem,max_line
        integer*4 imax,jmax
	parameter (max_elem = 15)
	parameter (max_line = 15)
	integer*4 line_dim,elem_dim
        real*4    r_grid_ratio

	byte image_data(elem_dim,line_dim)
        integer*4 image_to_dbz(0:15)

        real*4 laps_dbz(imax,jmax)
        real*4 wsi_dbz_data(nelems,nlines)
        real*4 r_llij_lut_ri(imax,jmax)
        real*4 r_llij_lut_rj(imax,jmax)

	real*4 line_mx,line_mn,elem_mx,elem_mn
	real*4 z_array(max_elem*max_line)
	real*4 zmax, zmin, zmean
        real*4 rdbz
        real*4 pixsum 
        real*4 result

        integer*4 i,j,ii,jj
        integer*4 istart,jstart
        integer*4 iend,jend
	integer*4 npix, nwarm
        integer*4 maxpix
        integer*4 ipix
        integer*4 istatus
        integer*4 qcstatus
        integer*4 fcount
        integer*4 bad_data_flag

        CALL ZERO(laps_dbz,IMAX,JMAX)
        istatus = -1
        bad_data_flag = 16
        qcstatus=0
        fcount=0
c
c The "10" loop represents input image resolution < output grid resolution such
c that there are enough pixels from the input image to get a representative
c mean value for the remapped output grid value
c
        if(r_grid_ratio .le. 0.5)then

          write(6,*)'Grid ratio .le. 0.5'   !0.75'
          write(6,*)'Use pixel avg to get dbz'
          DO 10 J=1,JMAX
          DO 10 I=1,IMAX

          IF(laps_dbz(I,J).NE.0.) GO TO 10
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
             write(*,*)'insufficient data for lat/lon sector'
                write(*,1020)i,j
1020	        format(1x,'LAPS grid (i,j) = ',i3,1x,i3)
                write(6,1021)elem_mx,elem_mn,line_mx,line_mn
1021            format(1x,'elem mx/mn  line mx/mn ',4f7.1)
             else
c
c **** FIND PIXELS AROUND GRID POINT
c
                zmax=-1.e15
                zmin=1.e15
                npix = 0
                pixsum = 0.
                maxpix = 0

                DO 3 JJ=jstart,jend
                DO 3 II=istart,iend

                  if(image_data(II,JJ).le.bad_data_flag.and.
     &image_data(II,JJ).ge.0)then
                      npix=npix+1
                      z_array(npix) = image_to_dbz(image_data(II,JJ))
                  endif

    3           continue  
 
                if(npix.gt.1)then
c
c...  this section finds the warmest pixel, coldest pixel, and mean pixel.
c
                   Do ii=1,npix

                      rdbz  = z_array(ii)
                      pixsum = pixsum + rdbz

                      if(rdbz.gt.zmax) then
                         zmax=rdbz
                      end if

                      if(rdbz.lt.zmin) then
                         zmin=rdbz
                      end if

                   enddo

                elseif(npix.eq.1)then

                   rdbz = z_array(npix)

                else   

                   rdbz=ref_base_cmn
                   fcount=fcount+1

                endif

                if(npix .ge. maxpix)maxpix=npix

                if(npix .gt. 1)then

                   if(pixsum.gt.0)then
                      laps_dbz(i,j) = pixsum / float(npix)
                   else
                      laps_dbz(i,j) = ref_base_cmn
                   endif

                elseif(rdbz.gt.0.0)then

                   laps_dbz(i,j) = rdbz

                else

                   laps_dbz(i,j) = ref_base_cmn

                endif ! npix .gt. 1

             end if  ! Enough data for num_lines .gt. 0

d           if(i .eq. i/10*10 .and. j .eq. j/10*10)then
d              write(6,5555)i,j,wm,wc,npix,nwarm,sc(i,j)
d5555         format(1x,2i4,2f10.2,2i5,f10.2)
d           endif

   10     CONTINUE ! I,J
          write(6,*)'Max num WSI pix for avg: ',maxpix
          write(6,*)'Number of LAPS gridpoints missing',
     &fcount

        else
c       ---------------------------------------
c input image resolution is large relative to output grid spacing
c this section uses bilinear interpolation to map
c the four surrounding input pixels to the output grid.
c
          write(6,*)'Image res .ge. output grid spacing'
          write(6,*)'Using bilinear interp to get dBZ '

          do j=1,nlines
          do i=1,nelems
             if(image_data(i,j).le.bad_data_flag.and.image_data(i,j)
     &.ge.0)then
                wsi_dbz_data(i,j)=image_to_dbz(image_data(i,j))
             else
                wsi_dbz_data(i,j)=r_missing_data_cmn
             endif
          enddo
          enddo

          DO 20 J=1,JMAX
          DO 20 I=1,IMAX

            IF(laps_dbz(I,J).NE.0.) GO TO 20
c
c bilinear_interp_extrap checks on boundary conditions and
c uses ref_base_cmn if out of bounds.
c
            call bilinear_laps(
     &           r_llij_lut_ri(i,j),
     &           r_llij_lut_rj(i,j),
     &           elem_dim,line_dim,wsi_dbz_data,
     &           result)

	    if(result .ne. r_missing_data_cmn .and.
     &         result .gt. 0.0)then

	        laps_dbz(i,j) = result

            else

                laps_dbz(i,j) = ref_base_cmn

            endif

   20     CONTINUE ! I,J

        end if     ! r_image_ratio

        istatus = 1
c
        return
        end



