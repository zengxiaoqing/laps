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
       SUBROUTINE satdat2laps_vis(
     &                  r_grid_ratio,
     &                  r_llij_lut_ri,
     &                  r_llij_lut_rj,
     &                  imax,jmax,
     &                  line_dim,elem_dim, ! image_vis array dimensions
     &                  image_vis,
     &                  sv,
     &                  istatus)
c
c.....       This is the Chandran version returing visible values
c.....       Modified by Steve Albers to average a window the size of a LAPS
c.....       Grid box, giving it better resolution and faster run times.
c           Changes:  P. Stamus  11-13-92  Install in operational sfc code.
c
c                     J. Smart   3-1-94    Install for using ISPAN feed for
c                                          satellite data.
c                                          Allow analysis technique to vary
c                                          as function of input/output ratio
       Implicit none

       Integer max_elem,max_line
       parameter (max_elem = 15)
       parameter (max_line = 15)

       Integer imax, jmax
       Integer line_dim, elem_dim
       Real*4 sv(IMAX,JMAX)
       Real*4 t_array(max_line*max_elem)
       Real*4 image_vis(elem_dim,line_dim)
       Real*4 r_llij_lut_ri(imax,jmax)
       Real*4 r_llij_lut_rj(imax,jmax)
       Real*4 elem_mn,elem_mx
       Real*4 line_mn,line_mx
       Integer npix
       Integer maxpix
       Integer i,j,ii,jj
       Integer istart
       Integer iend
       Integer jstart
       Integer jend
       Integer istatus
       Integer qcstatus

       Real*4 r_missing_data
       Real*4 Temp
       Real*4 r_grid_ratio
       Real*4 result
c
c -----------------------------begin--------------------------------
c
       call get_r_missing_data(r_missing_data, istatus)
       CALL ZERO(SV,IMAX,JMAX)
       istatus = 0
       qcstatus =0
c
c       grid_spacing_deg = sqrt( 
c    1       (  xlat(1,2) - xlat(1,1)                   )**2
c    1     + ( (xlon(1,2) - xlon(1,1))*cosd(xlat(1,1))  )**2 
c    1                         )   

c.....       Define half of the window dimensions

c      wdw_lat =  grid_spacing_deg / 2.
c      wdw_lon = (grid_spacing_deg / 2.) / cosd(xlat(1,1))
c      write(6,*)' GET VIS: wdw_lat, wdw_lon = ',wdw_lat,wdw_lon

       if(r_grid_ratio .lt. 0.5)then  !0.75)then
c      ------------------------------
c In this block the average pixel value is used for remapping the visible
c satellite to the output LAPS grid.
c
          write(6,*)'Image ratio .lt. 0.5'   !0.75'
          write(6,*)'Use pixel avg for VIS count'
          DO 10 J=1,JMAX
          DO 10 I=1,IMAX
            IF(SV(I,J).NE.0.) GO TO 10
c
c compute the line and element for window surrounding LAPS grid point.
c
            elem_mx = r_llij_lut_ri(i,j) + ((1./r_grid_ratio) * 0.5)
            elem_mn = r_llij_lut_ri(i,j) - ((1./r_grid_ratio) * 0.5)
            line_mx = r_llij_lut_rj(i,j) + ((1./r_grid_ratio) * 0.5)
            line_mn = r_llij_lut_rj(i,j) - ((1./r_grid_ratio) * 0.5)
            istart = nint(elem_mn+0.5)
            iend   = int(elem_mx)
            jstart = nint(line_mn+0.5)
            jend   = int(line_mx)

            if(istart .le. 0 .or. jstart .le. 0 .or.
     &iend .gt. elem_dim .or. jend .gt. line_dim)then
               write(*,*)'insufficient data for visible lat/lon sector'
               write(*,1020)i,j
1020	       format(1x,'LAPS grid (i,j) = ',i3,1x,i3)
               sv(i,j)=r_missing_data
 	    else
 
c **** FIND THE AVERAGE VISIBLE PIXELS AROUND GRID POINT
c
               Temp = 0.0
               npix = 0
               maxpix = 0

               DO 3 JJ=jstart,jend
                  DO 3 II=istart,iend

                  IF(image_vis(II,JJ) .eq. r_missing_data) GO TO 3

                     npix = npix + 1
                     t_array(npix) = image_vis(II,JJ)

    3          continue

               if(npix .ge. maxpix)maxpix=npix
c
c Added 11-6-95... A quality control test on the group of pixels used to
c                  derive the laps average pixel value.
c
               if(npix.ge.2)then

                  do ii=1,npix
                     temp = temp + t_array(ii)
                  enddo

                  sv(i,j) = Temp/ npix

               elseif(npix.eq.1)then

                  sv(i,j) = t_array(npix)

               else

                  sv(i,j) = r_missing_data

               endif

Cd             if(i .eq. i/10*10 .and. j .eq. j/10*10)then
Cd                write(6,5555)i,j,wm,wc,npix,nwarm,sc(i,j)
Cd5555            format(1x,2i4,2f10.2,2i5,f10.2)
Cd             endif

            end if

   10     CONTINUE ! I,J

          write(6,*)'Max num vis pixels for avg: ',maxpix
          write(6,*)'Number of vis satellite pixels modified'
          write(6,*)'           by statistical qc: ',qcstatus

       else
c      ----------------------------
c This block bilinearly interpolates four surrounding grid points to the
c output grid point.
c
          write(6,*)'grid ratio .ge. 0.5'  !0.75'
          write(6,*)'use bilinear interp for VIS count'
          DO 20 J=1,JMAX
          DO 20 I=1,IMAX

            IF(SV(I,J).NE.0.) GO TO 20
c
c line/elem are floating point i/j positions in ISPAN grid for input lat/lon
c 
c bilinear_interp_extrap checks on boundary conditions and
c uses r_missing_data if out of bounds.
c
            call bilinear_laps(
     &           r_llij_lut_ri(i,j),
     &           r_llij_lut_rj(i,j),
     &           elem_dim,line_dim,image_vis,
     &           result,istatus)

            if(result .lt. r_missing_data .and.
     &         result .gt. 0.0)then

                sv(i,j) = result

            else

                sv(i,j) = r_missing_data

            endif

   20     CONTINUE ! I,J

        end if

C       WRITE(6,1234) IB,I4VTIME,ICT
1234       FORMAT(1X,'BAND ',I4,' COUNT FOR I4TIME ',I10,' IS ',I8)

        istatus = 1
c
        return
        end
