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
     &                    maxradars,
     &                    radar_dist_min,
     &                    istatus)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Routine reads netCDF nowrad (WSI) data using subroutine read_wsi_cdf.
c Nowrad data is then remapped to LAPS domain given lat/lon of domain.
c Routine automatically moves the boundary for any domain with the nowrad
c confines.

      integer maxradars
      integer nsites_present
      integer nsites_absent
      real    present_site_loc_i(maxradars)
      real    present_site_loc_j(maxradars)
      real    absent_site_loc_i(maxradars)
      real    absent_site_loc_j(maxradars)
      real    present_site_lat(maxradars)
      real    present_site_lon(maxradars)
      real    absent_site_lat(maxradars)
      real    absent_site_lon(maxradars)
      real    radar_lat(maxradars)
      real    radar_lon(maxradars)
      real    radar_elev(maxradars)
      real    radar_ri(maxradars)
      real    radar_rj(maxradars)

      real lat(imax,jmax)
      real lon(imax,jmax)
      real height_grid(imax,jmax)
      real rdbz(imax,jmax)
      real radar_dist_min(imax,jmax)
      real ri(imax,jmax)
      real rj(imax,jmax)
      real distmin,dist
      real rdistmax
      real xgrddismx
      real ygrddismx
      real grid_spacing_m
      real dlat,dlon
      real LoV, Latin, La1,La2,Lo1,Lo2
      real centerLon,topLat,dx,dy
      real nw(2),se(2)
      real pi,rad2dg,rii,rjj
      real r_missing_data
      real heigth_mean_grid
      real height_mean
      real height_sum
      real ridis,rjdis

      integer image_to_dbz(0:15)
      integer validTime
      integer istatus
      integer jstatus
      integer status
      integer i_base_value
      integer increment

      character filename*200
      character ctype*3

      integer image(nelems,nlines)
      integer indxofthosein(maxradars)

      common /cegrid/nx,ny,nz,nw,se,rlatc,rlonc

c ----------------------------------------------------------------
c Read data and get navigation info to build lat/lon to i/j table
c ----------------------------------------------------------------
      istatus=1

      nradars_in = 0 ! initialize this

      call read_nowrad_cdf(ctype,filename,nlines,nelems,
     + dlat,dlon,La1,Lo1,La2,Lo2,centerlon,topLat,validTime,
     + dx,dy,Lov, Latin, image, istatus)
c
c this routine converts the raw bytes to integers scaled 0 - 64
c it also determines the wsi grid i/j locations for absent
c or presently reporting radars.

      call cvt_wsi_nowrad(ctype,nelems,nlines,image
     +,nsites_present,present_site_loc_i,present_site_loc_j
     +,nsites_absent,absent_site_loc_i,absent_site_loc_j
     +,maxradars,istatus)
      if(istatus .ne. 1)then
          write(6,*)'Bad status returned from cvt_wsi_nowrad'
          return
      endif

c     call cvt_wsi_nowrad(ctype,nelems,nlines,image,istatus)

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

c--------------------------------------------------
c compute the wsi grid i/j values for domain lat/lon
c--------------------------------------------------           
           call latlon_2_ceij(imax*jmax,lat,lon,ri,rj)

           call check_domain_vrc(imax,jmax,ri,rj,nx,ny)

c--------------------------------------------------
c compute radar site lat/lon locations both present
c and missing for the wsi mosaic.
c--------------------------------------------------
           call ceij_2_latlon(nsites_present,present_site_loc_i
     .,present_site_loc_j,present_site_lat,present_site_lon)

           call ceij_2_latlon(nsites_absent,absent_site_loc_i
     .,absent_site_loc_j,absent_site_lat,absent_site_lon)

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

c----------------------------------------------------
c Determine which radars in wsi grid are in this domain
c and compute/save the minimum distance to radar.
c----------------------------------------------------
      call get_r_missing_data(r_missing_data,istatus)
      call  read_static_grid(imax,jmax,'AVG',height_grid,istatus)

      if(ctype.eq.'wsi')then

         height_sum=0
         do j=1,jmax
         do i=1,imax
            height_sum=height_sum+height_grid(i,j)
         enddo
         enddo
         height_mean_grid=height_sum/(imax*jmax)
         nradars_tot=0
         rdistmax=480000.
         xgrddismx=imax+rdistmax/grid_spacing_m
         ygrddismx=jmax+rdistmax/grid_spacing_m
         xgrddismn=-rdistmax/grid_spacing_m
         ygrddismn=-rdistmax/grid_spacing_m
         height_sum=0.0

         print*,'mx num grid points for search x: ',xgrddismx
         print*,'mx num grid points for search y: ',ygrddismx
         print*,'mn num grid points for search x: ',xgrddismn
         print*,'mn num grid points for search y: ',ygrddismn

         do i=1,nsites_present
c
c fudge factor on site lat/lon due to unknown systematic
c offset.
           present_site_lat(i)=present_site_lat(i)-0.07
           present_site_lon(i)=present_site_lon(i)+0.02

           call latlon_to_rlapsgrid(present_site_lat(i),
     &                              present_site_lon(i),
     &                              lat,lon,            !LAPS lat/lon arrays
     &                              imax,jmax,          !LAPS horiz domain
     &                              rii,rjj,            !Output: real i,j, scalars
     &                              jstatus)

           ridis=abs(rii)
           rjdis=abs(rjj)

c          if(ridis.le.xgrddismx .and. rjdis.le.ygrddismx)then

           if( (rii.le.xgrddismx .and. rii.ge.xgrddismn) .and.
     &         (rjj.le.ygrddismx .and. rjj.ge.ygrddismn)  )then

              nradars_tot = nradars_tot + 1
              radar_lat(nradars_tot)=present_site_lat(i)
              radar_lon(nradars_tot)=present_site_lon(i)
              radar_ri(nradars_tot)=rii
              radar_rj(nradars_tot)=rjj
              radar_elev(nradars_tot)=r_missing_data

              if( (rii.gt.0.0  .and.  rii.lt.imax) .and.
     &            (rjj.gt.0.0  .and.  rjj.lt.jmax)  )then

                 call bilinear_interp_extrap(rii,rjj,imax,jmax
     &                    ,height_grid,result,istatus)

                 height_sum=height_sum+result
                 radar_elev(nradars_tot)=result
                 nradars_in=nradars_in+1
                 indxofthosein(nradars_in)=nradars_tot
              endif

           endif

        enddo

        if(nradars_in.gt.0)then
           print*
           print*,'Found ',nradars_in,' in this domain'
           height_mean=height_sum/nradars_in
        else
           print*
           print*,'Interesting! No radars in domain'
           print*,'set mean radar height to mean grid elev'
           print*,' = ',height_mean_grid
           height_mean = height_mean_grid
        endif
        print*
        print*,'Found ',nradars_tot,' to consider in min dist array'
        print*,'Mean height of radars = ',height_mean
        print*

        do i=1,nradars_tot
        do j=1,nradars_in
         ii=indxofthosein(j)
         if(ii.eq.i)then
          print*,'#/ri/rj/lat/lon/elev: ',ii,radar_ri(ii),radar_rj(ii)
     .,radar_lat(ii),radar_lon(ii),radar_elev(ii)
         endif
        enddo
        enddo

        where(radar_elev.eq.r_missing_data)radar_elev=height_mean_grid
        
        print*,'Determine nearest radar distance array'
        print*

        do j=1,jmax
        do i=1,imax

           distmin=r_missing_data
           do k=1,nradars_tot
              if(radar_elev(k).lt.10000.)then
                call latlon_to_radar(lat(i,j),lon(i,j),height_grid(i,j)
     1,azimuth,slant_range,elev,radar_lat(k),radar_lon(k),radar_elev(k))     
                if(slant_range.lt.distmin.and.slant_range.lt.rdistmax)
     &             distmin=slant_range
              endif
           enddo
           radar_dist_min(i,j)=distmin
       enddo
       enddo 

      endif

      goto 16
19    print*,'Error in nowrad_2_laps, terminating'
      goto 16
14    print*,'NOWRAD data not found for given time or'
      print*,'Data could be bad. No vrc produced'

16    print*,'Finished in nowrad_2_laps'
      return
      end
