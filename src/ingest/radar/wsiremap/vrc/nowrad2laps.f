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
       subroutine NOWRAD_to_LAPS(c_radtype,
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
       real*4 r_llij_lut_ri(imax,jmax)
       real*4 r_llij_lut_rj(imax,jmax)
       real*4 grid_spacing_m
       real*8 valTime

       integer image_to_dbz(0:15)
       integer validTime
       integer istatus
       integer status
       integer i_base_value,increment

cccccccccccccccccccccccccccccccccccccccccccccccccc

       character filename*200
       character c_radtype*3

       integer image(nelems,nlines)
c
c ***************************************************************************
c ************************** GET NOWRAD DATA ********************************
c
      istatus=1

      if(c_radtype.eq.'wfo')then
         call read_wsi_cdf_wfo(filename,nlines,nelems,
     1Dx,Dy,valTime,image,istatus)
         if(istatus.ne.0)then
            write(6,*)'Error reading nowrad data'
            goto 14
         end if
         validTime = int(valTime)
         write(*,*)'Valid Time ',validTime

      elseif(c_radtype.eq.'wsi')then

           call read_wsi_cdf_wsi(filename,nlines,nelems,
     +dlat,dlon,lat2,lon1,validTime,Dx,Dy,image,status)
           if(status.ne.0)then
              write(6,*)'Bad data detected'
              write(6,*)'Returning without data'
              goto 14
           endif

      endif

      write(*,*)' Found 5-min NOWrad data for requested time '
      write(*,*)
      write(*,*)'Valid Time ',validTime

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c    Build radar count level to dbz look up
c
      increment=5
      i_base_value=7
      do i = 1,15
         image_to_dbz(i)=(i-1)*increment+i_base_value
      end do
      image_to_dbz(0)=0
c
c get remapping look-up-table
c
      
      call readvrclut(c_radtype,imax,jmax,
     &     r_llij_lut_ri,r_llij_lut_rj,istatus)
c
c compute grid ratio 
c
      call get_grid_spacing(grid_spacing_m,istatus)
      r_grid_ratio=sqrt(Dx*Dx + Dy*Dy)/grid_spacing_m
c
c  Remap NOWrad onto LAPS grid.
c
       Call process_nowrad_z(imax,jmax,
     &                  r_grid_ratio,
     &                  image_to_dbz,
     &                  image,
     &                  r_llij_lut_ri,
     &                  r_llij_lut_rj,
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
