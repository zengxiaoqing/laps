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
      subroutine genlvdlut_sub(nx_l,ny_l,istatus)
c
      implicit none

      integer nx_l,ny_l

      real    lat(nx_l,ny_l)
      real    lon(nx_l,ny_l)
      real    topo(nx_l,ny_l)    !is not used.

      integer istatus
      integer it,js,lc
c
c dimensions for lat/lon
c
      character cgrid_fname*10

      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'
c
c ======================== START ==============================
c Acquiring LAPS latitude and longitude arrays.
c -------------------------------------------------------------------
      call find_domain_name(cgrid_fname,istatus)
      if(istatus.ne.1)then
         print*,'Error returned: find_domain_name'
         goto 900
      endif
      write(6,*)'genlvdlut_sub: Get lat/lon grid'
      call get_laps_domain(nx_l,ny_l,cgrid_fname,lat,lon,topo,istatus)
      if(istatus.ne.1)Then
         write(6,*)'Error - Unable to get lat/lon data'
         goto 900 
      end if
c-----------------------------------------------------------------
c
c Detemine which look-up-tables need generating based upon the specific
c lvd setup defined in satellite_master.nl. Specifically the namelist values
c isats, itypes, and ichannels indicate which sats, types, and channels within
c the master namelist get processed. Thus it, js, and lc are indices to these
c values.
c
      istatus = -1

      do js=1,maxsat

       if(isats(js).eq.1)then

        do it = 1,maxtype

         if(itypes(it,js).eq.1)then

          if(c_sat_types(it,js).eq.'cdf'.or.
     &       c_sat_types(it,js).eq.'wfo')then

           do lc=1,maxchannel

            if(ichannels(lc,it,js).eq.1)then

             write(6,59)c_sat_id(js),c_sat_types(it,js),
     &c_channel_types(lc,it,js)
59          format(1x,'Generating lookup table: ',a6,"/",a3,"/",a3)

             call gen_lut_lambert(js,it,lc,nx_l,ny_l,
     &                     lat,lon,istatus)

             if(istatus.eq.1)then
              write(6,*)'LUT generated'
             elseif(istatus.eq.0)then
              write(*,*)'ir LUT already generated'
             else
              write(6,*)'Error! LUT not generated ',
     &c_sat_id(js),'/',c_sat_types(it,js),'/',c_channel_types(lc,it,js)
             endif

            else
             write(6,49)c_sat_id(js),c_sat_types(it,js),
     &c_channel_types(lc,it,js)
            endif

           enddo

          elseif(c_sat_types(it,js).eq.'gvr'.or.
     &           c_sat_types(it,js).eq.'gwc')then

           do lc=1,maxchannel

            if(ichannels(lc,it,js).eq.1)then

             write(6,59)c_sat_id(js),c_sat_types(it,js),
     &c_channel_types(lc,it,js)
             call gen_gvarimage_lut(js,it,lc,nx_l,ny_l,lat,lon,
     &istatus)
             if(istatus.eq.1)then
              write(6,*)'LUT generated'
             elseif(istatus.eq.0)then
              write(*,*)'ir LUT already generated'
             else
              write(6,*)'LUT not generated ',c_channel_types(lc,it,js)
             endif
            else
             write(6,49)c_sat_id(js),c_sat_types(it,js),
     &c_channel_types(lc,it,js)
49           format(3x,a6,"/",a3,"/",a3,1x
     1                ,'not on in satellite namelist (via ichannels-2)')
            endif
           enddo

c         elseif(c_sat_types(it,js).eq.'asc')then
c
c          write(6,*)'Generate LUTs for ascii data'

c           do lc=1,maxchannel
c              lc=ichannels(lc,it,js)
c              call get_ascii_dimensions(path_to_data_cdf,
c    &c_channel_types(lch,i,j),nelem,nlines,istatus)

c              call gen_ascii_lut(path_to_data_cdf,
c    &csatid,c_sat_types(i,j),c_channel_types(lch,i,j),
c    &nelem,nlines,nx_l,ny_l,lat,lon,istatus)

c              if(istatus.eq.1)then
c                 write(6,*)'LUT generated'
c              else
c                 write(6,*)'LUT not generated ',c_channel_types(lch,i)
c              endif

c           enddo

          elseif(c_sat_types(it,js).ne.'     ')then

           write(6,*)'Unknown satellite data type! '
           write(6,*)'Check static/satellite_master.nl'

          endif

         elseif(c_sat_types(it,js).ne.'   ')then
          write(6,47)c_sat_id(js), c_sat_types(it,js)
47        format(2x,a6,"/",a3
     1                ,' not on in satellite namelist (via itypes)')

         endif ! if(itypes(it,js).eq.1)

        enddo

       else
        write(*,*)c_sat_id(js)
     1           ,' not on in satellite namelist (via isats)'
       endif

      enddo

900   return 
      end
