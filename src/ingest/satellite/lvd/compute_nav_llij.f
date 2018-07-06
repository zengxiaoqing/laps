      subroutine compute_nav_llij(nx_l,ny_l,maxchannels,nchannels,
     &csatid,csattype,chtype,isat,jtype,cfname_cur,gri,grj,nstatus)

      implicit none

      integer  nx_l,ny_l
      integer  maxchannels
      integer  nchannels
      integer  isat,jtype
      integer  istatus
      integer  gstatus
      integer  nstatus

      real     lat(nx_l,ny_l)
      real     lon(nx_l,ny_l)
      real     topo(nx_l,ny_l)
      real     rdumij(nx_l,ny_l,maxchannels)
      real     gri(nx_l,ny_l,maxchannels)
      real     grj(nx_l,ny_l,maxchannels)
      real     grid_spacing_laps_m

      character csatid*6
      character csattype*3
      character cfname_cur*9
      character chtype(maxchannels)*3
      character cgenericdataroot*255
      character c_gridfname*50

      logical   lut_flag

      nstatus=-1

c ---------------------------------------------
c acquiring LAPS latitude and longitude arrays.
c ---------------------------------------------
      call find_domain_name(cgenericdataroot,c_gridfname,istatus)
      call get_domain_laps(nx_l,ny_l,c_gridfname,lat,lon,topo,
     &grid_spacing_laps_m,istatus)
      if(istatus.eq.1)then
         write(6,*)'LAPS lat/lon/grid_spacing obtained'
         write(6,*)
      else
         write(6,*)'Error getting LAPS lat/lon data'
         stop
      end if
c
c try obtaining look-up-table. failure indicates no lut; generation required. 
c ---------------------------------------------------------------------------
c     call readlut(csatid,csattype,maxchannels,nchannels,
c    &chtype,nx_l,ny_l,rdumij,rdumij,istatus)

      if(istatus.eq.1)then

c        write(6,*)'LUT not obtained: ',csatid,'/',csattype
         print*,'-------------------------------------------'
         write(6,*)'Computing gri and grj with genlvdlut_lvd'
         print*,'-------------------------------------------'
c        call genlvdlut_sub(nx_l,ny_l,gstatus)
         call genlvdlut_lvd(nx_l,ny_l,maxchannels,lat,lon
     +,jtype,isat,gri,grj,gstatus)
         if(gstatus.lt.0)then
            write(6,*)'Error generating mapping arrays - terminating'
            goto 910
c        endif
c        else
c           print*,'Still in test mode and rewriting nav info'
c           print*,'Current Version: 9-29-06'
c           print*,'rewrite satellite_lvd.nl'
c           print*
c           call rewrite_satellite_lvd_nl(istatus)
         else
            write(6,*)' Returned from genlvdlut_lvd'
         endif

         nstatus=1

      endif

c     else
c        write(6,*)'Compute ll/ij mapping arrays '
c        write(6,*)'Check if luts are up-to-date'
c        call check_luts(cfname_cur,isat,jtype,
c    &chtype,maxchannels,nchannels,lut_flag,istatus)
c        if(lut_flag.and.istatus.eq.0)then
c           print*,'************************************************'
c           write(6,*)'Found difference in nav parms',
c    +' - rebuild the lut'
c           print*,'************************************************'
c           call genlvdlut_lvd(nx_l,ny_l,lat,lon,jtype,isat,
c    +gstatus)
c           if(gstatus.lt.0)then
c              write(6,*)'Error generating LUT - terminating'
c              goto 910
c           else
c              write(6,*)'**********************************'
c              write(6,*)
c              call rewrite_satellite_lvd_nl(istatus)
c           endif
c           nstatus=1
c        elseif(istatus.eq.0)then
c           write(6,*)'Lut checked out ok'
c           write(6,*)
c           nstatus=0
c        else
c           write(6,*)'Error status returned from check_lut'
c           goto 910
c        endif
c     endif

      print*,' Returning from compute_nav_llij '

      goto 1000

910   print*,' ERROR: compute_nav_llij '

1000  return
      end
