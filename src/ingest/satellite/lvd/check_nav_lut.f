      subroutine check_nav_lut(nx_l,ny_l,maxchannels,nchannels,
     &csatid,csattype,chtype,isat,jtype,cfname_cur,nstatus)

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
      real     grid_spacing_laps_m

      character csatid*6
      character csattype*3
      character cfname_cur*9
      character chtype(maxchannels)*3

      logical   lut_flag

      nstatus=-1

c ---------------------------------------------
c acquiring LAPS latitude and longitude arrays.
c ---------------------------------------------
      call get_domain_laps(nx_l,ny_l,'nest7grid',lat,lon,topo,
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
      call readlut(csatid,csattype,maxchannels,nchannels,
     &chtype,nx_l,ny_l,rdumij,rdumij,istatus)

      if(istatus.eq.1)then

         write(6,*)'LUT not obtained: ',csatid,'/',csattype
         write(6,*)'Computing lut using genlvdlut_lvd'
c        call genlvdlut_sub(nx_l,ny_l,gstatus)
         call genlvdlut_lvd(nx_l,ny_l,lat,lon,jtype,isat,
     +gstatus)
         if(gstatus.lt.0)then
            write(6,*)'Error generating LUT - terminating'
            goto 910
         else
            write(6,*)'rewrite satellite_lvd.nl'
            write(6,*)
            call rewrite_satellite_lvd_nl(istatus)
         endif
         nstatus=1
      else

         write(6,*)'Got the mapping look-up-tables '
         write(6,*)'Check if luts are up-to-date'

         call check_luts(cfname_cur,isat,jtype,
     &chtype,maxchannels,nchannels,lut_flag,istatus)

         if(lut_flag.and.istatus.eq.0)then
            print*,'************************************************'
            write(6,*)'Found difference in nav parms',
     +' - rebuild the lut'
            print*,'************************************************'
            call genlvdlut_lvd(nx_l,ny_l,lat,lon,jtype,isat,
     +gstatus)
            if(gstatus.lt.0)then
               write(6,*)'Error generating LUT - terminating'
               goto 910
            else
               write(6,*)'**********************************'
               write(6,*)
               call rewrite_satellite_lvd_nl(istatus)
            endif
            nstatus=1
         elseif(istatus.eq.0)then
            write(6,*)'Lut checked out ok'
            write(6,*)
            nstatus=0
         else
            write(6,*)'Error status returned from check_lut'
            goto 910
         endif

      endif

      goto 1000

910   print*,' error in check_nav_lut '

1000  return
      end
