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
       subroutine NOWRADWSI_to_LAPS(ctype,
     &                    filename,
     &                    nlines,nelems,
     &                    imax,jmax,
     &                    lat,lon,
     &                    validTime,
     &                    rdbz,
     &                    istatus)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Routine reads netCDF nowrad (WSI) data using subroutine read_wsi_cdf.
c Nowrad data is then remapped to LAPS domain given lat/lon of domain.
c Routine automatically moves the boundary for any domain with the nowrad
c confines.

       real*4 lat(imax,jmax)
       real*4 lon(imax,jmax)
       real*4 rdbz(imax,jmax)
       real*4 ri(imax,jmax)
       real*4 rj(imax,jmax)
       real*4 grid_spacing_m
       real*4 dlat,dlon
       real*4 LoV, Latin, La1,La2,Lo1,Lo2
       real*4 centerLon,topLat,dx,dy
       real*4 nw(2),se(2)
       real*4 pi,rad2dg

       integer image_to_dbz(0:15)
       integer validTime
       integer istatus
       integer status
       integer i_base_value
       integer increment

       character filename*200
       character ctype*3

       integer image(nelems,nlines)

       common /cegrid/nx,ny,nz,nw,se,rlatc,rlonc

c ----------------------------------------------------------------
c Read data and get navigation info to build lat/lon to i/j table
c ----------------------------------------------------------------
      istatus=1

      call read_nowrad_cdf(ctype,filename,nlines,nelems,
     + dlat,dlon,La1,Lo1,La2,Lo2,centerlon,topLat,validTime,
     + dx,dy,Lov, Latin, image, istatus)
c
c this routine converts the raw bytes to integers scaled 0 - 64
c
      call cvt_wsi_nowrad(ctype,nelems,nlines,image,istatus)

      pi=acos(-1.)
      rad2dg=180.0/pi

      if(ctype.eq.'wfo')then

         Lo1=-Lo1
         Lo2=-Lo2
         if(Lov.lt.-180.0)Lov=Lov+360.
         if(Lo1.lt.-180.0)Lo1=Lo1+360.
         if(Lo2.lt.-180.0)Lo2=Lo2+360.

         if(Lov.gt.180.0)Lov=Lov-360.
         if(Lo1.gt.180.0)Lo1=Lo1-360.
         if(Lo2.gt.180.0)Lo2=Lo2-360.

         call gen_rirj_lam(imax,jmax,lat,lon,nelems,nlines
     &        ,La1,Lo1,La2,Lo2,Latin,Lov,dx,dy,ri,rj
     &        ,istatus)
c
 
      elseif(ctype.eq.'wsi')then
c
c load common block for cyclindrical equidistant
c
           nx=nelems
           ny=nlines
           nz=1
           nw(1)=la1
           nw(2)=lo1
           se(1)=la2
           se(2)=lo2
           rlatc=(topLat - (dlat*(nlines-1)*0.5))*rad2dg
           rlonc=centerLon*rad2dg
           
           call latlon_2_ceij(imax*jmax,lat,lon,ri,rj)

           call check_domain_vrc(imax,jmax,ri,rj,nx,ny)

      endif

      write(*,*)' NOWRAD data prepared for requested time '
      write(*,*)
      write(*,*)'Valid Time ',validTime

c ---------------------------------------------------
c    Build radar count level to dbz look up
c ---------------------------------------------------
      increment=5
      i_base_value=7
      do i = 1,15
         image_to_dbz(i)=(i-1)*increment+i_base_value
      end do
      image_to_dbz(0)=0
c ---------------------------------------------------
c compute grid ratio 
c ---------------------------------------------------
      call get_grid_spacing(grid_spacing_m,istatus)
      r_grid_ratio=sqrt(Dx*Dx + Dy*Dy)/grid_spacing_m
c ---------------------------------------------------
c  Remap NOWRAD-WSI data to LAPS grid.
c ---------------------------------------------------
       Call process_nowrad_z(imax,jmax,
     &                  r_grid_ratio,
     &                  image_to_dbz,
     &                  image,
     &                  ri,
     &                  rj,
     &                  nlines,nelems, ! input array dimensions
     &                  rdbz,
     &                  istatus)

       goto 16
19     write(6,*)'Error in nowrad_2_laps, terminating'
       goto 16
14     write(*,*)'NOWRAD data not found for given time or'
       write(6,*)'Data could be bad. No vrc produced'

16     write(6,*)'Finished in nowrad_2_laps'
       return
       end
