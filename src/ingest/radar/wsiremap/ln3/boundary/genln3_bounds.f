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
        Program gen_3dradar_bounds
cc
c
c program builds static file for WSI radar processor. Requires lat/lon static
c file for laps. Output is i and j starting and ending points for processing the
c WSI 3d radar data.  This program is run when localizing laps since ln3_driver.exe
c needs the correct istart/jstart/iend/jend for the domain settings.
c
c J Smart  9/95 original
c J Smart  9/97 newlaps/dynamic memory mods
c    "     9/98 removed routine readwrite_ln3_parms and added
c               routines to get array dimensions and other file
c               header stuff directly from NIMBUS data file.

      character cfilespec*100
      character cftime*9
      integer   i4time_latest
      integer   max_files
      integer   dataLevel,elems,lines
      parameter (max_files=1000)

      character c_filenames(max_files)*200

      include 'lapsparms.cmn'
       
      Call get_laps_config('nest7grid',IStatus)
      if(istatus.eq.1)then
         write(*,*)'LAPS Parameters obtained'
      else
         write(*,*)'IStatus = ',IStatus,'Error - Get_LAPS_Config'
         write(*,*)'Terminating LAPS-VRC. WSI remapping'
         stop
      end if

      call s_len(path_to_wsi_3d_radar_cmn,lp)

      cfilespec=path_to_wsi_3d_radar_cmn(1:lp)//'*_ll'

      call get_file_names(cfilespec,numoffiles,c_filenames
     1        ,max_files,istatus)
      if(istatus.ne.1.or.numoffiles.eq.0)then
         call s_len(cfilespec,nn)
         print*,'No files found for cfilespec:'
         print*,cfilespec(1:nn)
         print*,'boundary info not updated - Terminating'
         stop
      endif

      call get_wsi_3drad_dims(c_filenames(numoffiles),
     &    dataLevel,elems,lines,istatus)

      call genradarbounds_sub(c_filenames(numoffiles),
     &      NX_L_CMN,NY_L_CMN,dataLevel,elems,lines)

      stop
      end

c ----------------------------------------------------

      subroutine genradarbounds_sub(cfname,nx_l,ny_l,
     & dataLevel,elems,lines)

      Implicit NONE

      integer  nx_l,ny_l

      real*4 lat(nx_l,ny_l)
      real*4 lon(nx_l,ny_l)
      real*4 data(nx_l,ny_l,2)

      integer istart,iend,jstart,jend
      integer istart1,iend1,jstart1,jend1
      integer istart2,iend2,jstart2,jend2
c
      real*4 rdtodg, dlat_deg, dlon_deg
      real*4 lon1_deg, lat2_deg
      real*4 rlat1,rlat2,rlon1,rlon2,rdlat,rdlon
      real*4 wsi_lat, wsi_lon
      real*4 rstartlat1,rstartlat2
      real*4 rstartlon1,rstartlon2
      real*4 rendlat1,rendlat2
      real*4 rendlon1,rendlon2
      real*4 grid_spacing

      integer i,j
      integer istatus, istatus_io
      integer lend
      integer lines,elems
      character*150 dir_static
      character*125 comment_ll(2)
      character*10 units_ll(2)
      character*3 var_ll(2)

      integer dataLevel

      logical       lnorth_flag
      logical       lsouth_flag
      logical       lwest_flag
      logical       least_flag

      character*50  grid_name
      character*50  grid_type
      character*50  level_prefix(dataLevel)
      character*20  product_units
      character*50  x_dim
      character*50  y_dim
      character*(*) cfname
      integer Nx, Ny, data_levels(dataLevel), 
     +   num_levels, validTime
      real Dx,Dy,centerLon,diffLon,topLat

c
c  BEGIN
c
c
      call get_directory('static',dir_static,lend)
      var_ll(1)='LAT'
      var_ll(2)='LON'

      call rd_laps_static(dir_static, 'nest7grid', nx_l, ny_l, 2,
     &     var_ll, units_ll, comment_ll, data, grid_spacing,
     &     istatus)

      if(istatus.eq.1)then
          write(*,*)'LAPS lat/lon grid obtained'
          write(*,*)
          do j=1,ny_l
             do i=1,nx_l
                lat(i,j)=data(i,j,1)
                lon(i,j)=data(i,j,2)
             end do
          end do
      else
          write(*,*)'Unable to get lat/lon data'
          write(*,*)'generate_radar_bounds process terminating'
          stop
      end if

      rdtodg=180.0/3.1415926
c
c wsi data starts (ie, [1,1]) in NW corner  
c
c      call readwrite_ln3_parms('noinstall',lines,elements,
c     &rdlat,rdlon,rlat2,rlon1,istart,iend,jstart,jend,istatus_io)

      call read_cdf_wsi3d_head(cfname,dataLevel, Dx, Dy,
     +    rlat1, rlat2, rlon1, rlon2, Nx, Ny, centerLon, data_levels,
     +    diffLon, grid_name, grid_type, level_prefix,
     +    num_levels, product_units, rdlon, rdlat,
     +    topLat, validTime, x_dim, y_dim,istatus)


      dlat_deg=rdlat*rdtodg
      dlon_deg=rdlon*rdtodg

      istart=999999
      jstart=999999
      iend  =0
      jend  =0
      istart1=999999
      jstart1=999999
      iend1  =0
      jend1  =0
      istart2=999999
      jstart2=999999
      iend2  =0
      jend2  =0
      lnorth_flag=.false.
      lsouth_flag=.false.
      lwest_flag=.false.
      least_flag=.false.


      do i=1,nx_l
         if(lat(i,ny_l).gt.rlat1)lnorth_flag=.true.
         if(lat(i,1).lt.rlat2)lsouth_flag=.true.
      enddo
      do j=1,ny_l
         if(lon(1,j).lt.rlon1)lwest_flag=.true.
         if(lon(nx_l,j).gt.rlon2)least_flag=.true.
      enddo

      do j =1,lines

         wsi_lat = rlat1 - (dlat_deg * (j-1))

         do i = 1,elems
c
c determine starting/ending i and j points for the remapping
c process.
c dateline adjustment= if(wsi_lon.gt.180.0)wsi_lon=-180.0+(wsi_lon-180.0)
c

            wsi_lon = rlon1 + (dlon_deg * (i-1))
c
c find istarting position; either SW or NW corner, furthest western point.
c
            if( (wsi_lat.le.lat(1,2)) .and.
     &          (wsi_lat.ge.lat(1,1)) .and.
     &          (wsi_lon.le.lon(2,1)) .and.
     &          (wsi_lon.ge.lon(1,1)) .and.
     &         istart1 .eq. 999999)then
               istart1 = i-15
               rstartlon1=wsi_lon
            endif
            if( (wsi_lat.le.lat(1,ny_l)) .and.
     &          (wsi_lat.ge.lat(1,ny_l-1)) .and.
     &          (wsi_lon.ge.lon(1,ny_l)) .and.
     &          (wsi_lon.le.lon(2,ny_l)) .and.
     &         istart2 .eq. 999999)then
               istart2=i-15
               rstartlon2=wsi_lon
            endif
c
c find jstarting position; either NW or NE corner, furthest northern point.
c
            if( (wsi_lat.le.lat(1,ny_l))   .and.
     &          (wsi_lat.ge.lat(1,ny_l-1)) .and.
     &          (wsi_lon.ge.lon(1,ny_l))   .and.
     &          (wsi_lon.le.lon(2,ny_l))   .and.
     &         jstart1 .eq. 999999)then
               jstart1 = j-15
               rstartlat1=wsi_lat
            endif
            if( (wsi_lat.le.lat(nx_l,ny_l)) .and.
     &          (wsi_lat.ge.lat(nx_l,ny_l-1)).and.
     &          (wsi_lon.le.lon(nx_l,ny_l)) .and.
     &          (wsi_lon.ge.lon(nx_l-1,ny_l)).and.
     &         jstart2 .eq.999999)then
               jstart2=j-15
               rstartlat2=wsi_lat
            endif
c
c find iending position; either NE or SE corner, furtherest eastern point.
c
            if( (wsi_lat.ge.lat(nx_l,1)) .and.
     &          (wsi_lat.le.lat(nx_l,2)) .and.
     &          (wsi_lon.le.lon(nx_l,1)) .and.
     &          (wsi_lon.ge.lon(nx_l-1,1)).and.
     &         iend1 .eq. 0)then
               iend1 = i+15
               rendlon1=wsi_lon
            endif
            if( (wsi_lat.le.lat(nx_l,ny_l)) .and.
     &          (wsi_lat.gt.lat(nx_l,ny_l-1)) .and.
     &          (wsi_lon.le.lon(nx_l,ny_l)) .and.
     &          (wsi_lon.ge.lon(nx_l-1,ny_l)) .and.
     &         iend2 .eq. 0)then
               iend2 = i+15
               rendlon2=wsi_lon
            endif
c
c find jending position; either SE or SW corner, furtherest southern point.
c
            if( (wsi_lat.ge.lat(nx_l,1)) .and.
     &          (wsi_lat.le.lat(nx_l,2)) .and.
     &          (wsi_lon.ge.lon(nx_l-1,1)) .and.
     &          (wsi_lon.le.lon(nx_l,1)) .and.
     &         jend1 .eq. 0)then
               jend1 = j+15
               rendlat1=wsi_lat
            endif
            if( (wsi_lat.ge.lat(1,1))  .and.
     &          (wsi_lat.le.lat(1,2))  .and.
     &          (wsi_lon.ge.lon(1,1)) .and.
     &          (wsi_lon.le.lon(2,1)) .and.
     &         jend2 .eq. 0)then
               jend2 = j+15
               rendlat2=wsi_lat
            endif

34       end do    !all elems
35    end do      !all lines

      if(lnorth_flag)then
         jstart1=1
      endif
      if(lsouth_flag)then
         jend1=ny
      endif
      if(lwest_flag)then
         istart1=1
      endif
      if(least_flag)then
         iend1=nx
      endif

      write(6,*)'istart1,jstart1,iend1,jend1:',
     &istart1,jstart1,iend1,jend1
      write(6,*)'istart2,jstart2,iend2,jend2:',
     &istart2,jstart2,iend2,jend2

      write(6,*)'rstartlat1/rstartlat2 ',rstartlat1,rstartlat2
      write(6,*)'rstartlon1/rstartlon2 ',rstartlon1,rstartlon2
      write(6,*)'rendlat1/rendlat2 ',rendlat1,rendlat2
      write(6,*)'rendlon1/rendlon2 ',rendlon1,rendlon2

      jstart=min(jstart1,jstart2)
      jend=max(jend1,jend2)
      istart=min(istart1,istart2)
      iend=max(iend1,iend2)

      if(jstart .le. 1)then
         write(6,*)'Nrthn LAPS Bndry may be out of radar Domain'
      endif
      if(jend .ge. ny)then
         write(6,*)'Sthrn LAPS Bndry may be out of radar Domain'
      endif
      if(istart .le. 1)then
         write(6,*)'Wstrn LAPS Bndry may be out of radar Domain'
      endif
      if(iend .ge. nx)then
         write(6,*)'Estrn LAPS Bndry may be out of radar Domain'
      endif

c     call readwrite_ln3_parms('install',nlines,nelements,
c    &rdlat,rdlon,rlat2,rlon1,istart,iend,jstart,jend,istatus_io)
c
c rewrite the i/j start/end values to ln3 namelist
c
      call rewrite_ijse(istart,jstart,iend,jend,istatus)
      if(istatus.ne.0)then
         print*,'error returned from rewrite_ijse'
         goto 1000
      endif

      goto 999

997   write(6,*)'istart/iend and jstart/jend not updated'
      return

999   write(6,*)'istart/jstart and iend/jend written'
      return
1000  end
