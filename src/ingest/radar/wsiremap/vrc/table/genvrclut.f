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
      program gen_vrc_llij_lut
      implicit none
      include 'lapsparms.cmn'      
      integer istatus
      Call get_laps_config('nest7grid',IStatus)
      if(IStatus.eq.1)then
         write(6,*)'LAPS Parameters obtained'
      else
          write(6,*)'IStatus = ',IStatus,'Error - Get_LAPS_Config'
          write(6,*)'Terminating LAPS-LVD. ISPAN Satellite remapping'
          stop
      end if
c

      call gen_vrc_llij_lut_sub(nx_l_cmn,ny_l_cmn,c_raddat_type)

      stop
      end
c
c-------------------------------------------------------------
      subroutine gen_vrc_llij_lut_sub(nx_l,ny_l,c_raddat_type)
c
      implicit none

      integer nx_l, ny_l
      character*3 c_raddat_type
      real*4    lat(nx_l,ny_l)
      real*4    lon(nx_l,ny_l)
      real*4    grid_spacing
      real*4    data(nx_l,ny_l,2)
      real*4    rdum
      integer   i,j
      integer   istatus
      integer   len
      integer   nlines,nelems
      logical   test/.false./
c
c dimensions for lat/lon
c
      character*125 comment_ll(2)
      character*10 units_ll(2)
      character*3 var_ll(2)
      character*200 dir_static
c
c
c Definitions needed for acquiring LAPS latitude and longitude arrays.
c -------------------------------------------------------------------
      call get_directory('static',dir_static,len)

      var_ll(1) = 'LAT'
      var_ll(2) = 'LON'

      write(6,*)'Get LAPS lat/lon grid'
      call rd_laps_static(dir_static,'nest7grid',nx_l,ny_l, 2,
     &     var_ll, units_ll, comment_ll, data, grid_spacing,
     &     istatus)

      if(istatus.eq.1)then
         write(6,*)'LAPS lat/lon grid obtained'
         write(6,*)
         do j=1,ny_l
            do i=1,nx_l
               lat(i,j)=data(i,j,1)
               lon(i,j)=data(i,j,2)
            end do
         end do
      else
         write(6,*)'Unable to get lat/lon data'
         write(6,*)'lvd process terminating'
         stop
      end if
c
c
c      do i=1,n_radar_types
      
      if(c_raddat_type.eq.'wfo')then

         write(6,*)'Gen LUT for Conus (lambert) wsi (WFO) data '
         call gen_vrc_wfo_cdf_lut(2,nx_l,ny_l,lat,lon,istatus)

         if(istatus.eq.1)then
            write(6,*)'CDF look-up table generated'
         else
            write(6,*)'Error generating look-up-table'
            goto 901
         endif

      elseif(c_raddat_type.eq.'wsi')then

         write(6,*)'Gen LUT for wsi - Cylindrical Equidistant'

         if(test)then

            print*,'Using latlon_to_ceij routine'
            call gen_llij_lut_wsi(1,nx_l,ny_l,lat,lon,c_raddat_type
     +,istatus)

         else

            print*,'Using inefficient latlon_to_ll routine'
            call get_wsi_parms_vrc(1,nlines,nelems,
     +rdum,rdum,rdum,rdum,rdum,rdum,rdum,rdum,rdum,
     +istatus) 

            call gen_llij_lut_ll(1,nx_l,ny_l,lat,lon,c_raddat_type
     +,nlines,nelems,istatus)

         endif

         if(istatus.eq.1)then
            write(6,*)'CDF look-up table generated'
         else
            write(6,*)'Error generating look-up-table'
            goto 901
         endif
      else
        print*, 'Error in nest7grid.parms c_raddat_type:',c_raddat_type
        goto 901
      endif

c      enddo
c
      goto 900

901   write(6,*)'Error - genvrclut.f'

900   stop
      end
