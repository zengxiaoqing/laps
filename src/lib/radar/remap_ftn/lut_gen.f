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
      subroutine lut_gen(c4_radarname,rlat_radar,rlon_radar
     :                  ,rheight_radar,NX_L,NY_L,NZ_L)
c
c     PURPOSE:
c        Generate look-up tables for radar remapping.
c
      implicit none
c
      include 'remap_constants.dat'
      include 'remap.cmn'
c
c     Variables from LAPS domain file
c
      integer NX_L,NY_L,NZ_L
      real*4 lat(NX_L,NY_L)
      real*4 lon(NX_L,NY_L)
      real*4 topo(NX_L,NY_L)
c
c     Functions
c
      real height_to_zcoord
c
c     Misc interval variables
c
      integer i,j,igate_lut,iaz,ielev,iran,iz_grid,istatus,len_dir
      real rlat_grid,rlon_grid,height_grid
      real rlat_radar,rlon_radar,rheight_radar
      real elev,elev_deg,coselev,azimuth,azi_deg
      real slant_range,sl_range_m,ri,rj,dbz,z
      real cosd
      character*4 c4_radarname
      character*150 static_dir,filename
      character*3 ext
c
c     Fill arrays with initial values
c
      DO 50 i = 1,10
        i4time_old(i) = 0
        n_ref_obs_old(i) = 99999
   50 CONTINUE
c
      write(6,801)NX_L,NY_L
  801 format('REMAP > Getting LAPS Domain: ',2i5)
c
      call get_laps_domain(NX_L,NY_L,'nest7grid',
     :                      lat,lon,topo,istatus)
      IF (istatus .eq. 0) THEN
        write(6,*)' Error getting LAPS domain'
        stop
!       RETURN
      END IF

      write(6,*)' Corners of domain:'
      write(6,*)lat(1,1),lon(1,1),lat(NX_L,NY_L),lon(NX_L,NY_L)

c
      write(6,810) range_interval
  810 format(' Range interval for LUTs is ',F10.2)
      write(6,812) gate_spacing_m
  812 format(' Effective gate spacing for radar ',F10.2)

      write(6,*)' Radar Name  ',c4_radarname
      write(6,*)' Radar Coords',rlat_radar,rlon_radar,rheight_radar
      
      ext = 'dat'
      call get_directory(ext,static_dir,len_dir)

c     Put the coords into common so remap_process can access them
      rlat_radar_cmn = rlat_radar
      rlon_radar_cmn = rlon_radar
      rheight_radar_cmn = rheight_radar
      c4_radarname_cmn = c4_radarname      
c
c     Generate Gate/Elev to Projran lut
c
!     Try to read lut
      filename = static_dir(1:len_dir)//'vxx/'
     1         //'gate_elev_to_projran_lut.'//c4_radarname
      write(6,*)' Reading file: ',filename
      open(11,file=filename,form='unformatted',status='old',err=90)
      read(11,err=90)gate_elev_to_projran_lut
      close(11)
      goto 110
90    write(6,*)' Generating LUT - no valid file exists'

!     Calculate lut
      write(6,820) lut_gates
  820 format(' REMAP > Building Gate/Elev to Projran lut, gates =',
     :        I12)

      DO ielev = 0,lut_elevs

        elev_deg = ielev * elev_interval
        coselev = cosd(elev_deg)

        DO igate_lut = 1,lut_gates

          sl_range_m = igate_lut * gate_spacing_m * gate_interval
          iran = sl_range_m * coselev / range_interval
          gate_elev_to_projran_lut(igate_lut,ielev) = iran

        ENDDO
      ENDDO

!     Write lut
      open(12,file=filename,form='unformatted',status='new')
      write(12)gate_elev_to_projran_lut
      close(12)
  110 continue 

c
c     Generate Gate/Elev to Z lut
c
!     Try to read lut
      filename = static_dir(1:len_dir)//'vxx/'
     1         //'gate_elev_to_z_lut.'//c4_radarname
      write(6,*)' Reading file: ',filename
      open(11,file=filename,form='unformatted',status='old',err=190)
      read(11,err=190)gate_elev_to_z_lut
      close(11)
      goto 210
190   write(6,*)' Generating LUT - no valid file exists'

!     Calculate lut
      write(6,830)
  830 format(' REMAP > Building Gate/Elev to Z lut')
c
      DO 200 ielev = 0,lut_elevs
 
        elev_deg = ielev * elev_interval
 
        DO 180 igate_lut = 1,lut_gates

          sl_range_m = igate_lut * gate_spacing_m * GATE_INTERVAL
          azi_deg=0.
 
          call radar_to_latlon
     :         (rlat_grid,rlon_grid,height_grid
     :                   ,azi_deg,sl_range_m,elev_deg
     :                  ,rlat_radar,rlon_radar,rheight_radar)

          iz_grid = nint(height_to_zcoord(height_grid,istatus))
          gate_elev_to_z_lut(igate_lut,ielev) = min(iz_grid,NZ_L)

  180   CONTINUE
  200 CONTINUE

!     Write lut
      open(12,file=filename,form='unformatted',status='new')
      write(12)gate_elev_to_z_lut
      close(12)
  210 continue 

c
c     Generate Az/Ran to i,j lut
c
!     Try to read lut
      filename = static_dir(1:len_dir)//'vxx/'
     1         //'azran_to_ijgrid_lut.'//c4_radarname
      write(6,*)' Reading file: ',filename
      open(11,file=filename,form='unformatted',status='old',err=290)
      read(11,err=290)azran_to_igrid_lut,azran_to_jgrid_lut
      close(11)
      goto 310
290   write(6,*)' Generating LUT - no valid file exists'

!     Calculate lut
      write(6,840)
  840 format(' REMAP > Building Az/Ran to i,j lut')
c
      DO 300 iran = 0,lut_ranges

        slant_range = iran*range_interval

        DO 280 iaz = 0,lut_azimuths

          azimuth = float(iaz)
          elev = 0.

          call radar_to_latlon(rlat_grid,rlon_grid,height_grid
     :                  ,azimuth,slant_range,elev
     :                  ,rlat_radar,rlon_radar,rheight_radar)

          call latlon_to_rlapsgrid(rlat_grid,rlon_grid,lat,lon,
     :                   NX_L,NY_L,ri,rj,istatus)
          i = nint(ri)
          j = nint(rj)

          IF (i.le.0.or.i.gt.NX_L.or.j.le.0.or.j.gt.NY_L) THEN

            azran_to_igrid_lut(iaz,iran) = 0
            azran_to_jgrid_lut(iaz,iran) = 0

          ELSE

            azran_to_igrid_lut(iaz,iran) = i
            azran_to_jgrid_lut(iaz,iran) = j

          END IF

  280   CONTINUE
  300 CONTINUE

!     Write lut
      open(12,file=filename,form='unformatted',status='new')
      write(12)azran_to_igrid_lut,azran_to_jgrid_lut
      close(12)
  310 continue 

c     Generate DbZ Lookup Table (graduated for each tenth of a dbz)
      do i = -1000,+1000
          dbz = float(i) / 10.
          z = 10**(dbz / 10.)
          dbz_to_z_lut(i) = z
      enddo ! i

c     These lookup tables flag which gates actually need processing
      do i = 1,MAX_GATES
           
          if(i .le. 920       .and. i .ge. INITIAL_VEL_GATE)then
              lgate_vel_lut(i) = .true.
          else
              lgate_vel_lut(i) = .false.
          endif

          if(i .eq. (i/4) * 4 .and. i .ge. INITIAL_REF_GATE)then
              lgate_ref_lut(i) = .true.
          else
              lgate_ref_lut(i) = .false.
          endif

          if(lgate_vel_lut(i) .or. lgate_ref_lut(i))then
              lgate_lut(i) = .true.
          else
              lgate_lut(i) = .false.
          endif

      enddo ! i


      write(6,850)
  850 format(' REMAP > Lookup Tables Complete')

      RETURN
      END

   


      subroutine read_radar_info(c4_radarname,rlat_radar,rlon_radar
     :                                       ,rheight_radar,istatus)

      character*4 c4_radarname
      character*80 c80_line

      call getenv('RADARNAME',c4_radarname)

      write(6,*)'read_radar_info: RADARNAME = ',c4_radarname

      open(11,file='radarinfo.dat',status='old')

 10   read(11,1,err=900)c80_line
 1    format(a)

      if(c80_line(1:4) .eq. c4_radarname)then ! We've got the right radar
          read(c80_line,2)ideg_lat,imin_lat,isec_lat
     1                   ,ideg_lon,imin_lon,isec_lon
     1                   ,iheight
 2        format(26x,i2,6x,i2,6x,i2,7x,i3,5x,i2,6x,i2,6x,i4)
          rlat_radar =  float(ideg_lat) + float(imin_lat)/60.
     1                                  + float(isec_lat)/3600.
          rlon_radar = -float(ideg_lon) - float(imin_lon)/60.
     1                                  - float(isec_lon)/3600.
          rheight_radar = iheight
          write(6,*)' rlat_radar,rlon_radar,rheight_radar '
     1               ,rlat_radar,rlon_radar,rheight_radar
          istatus = 1
          close(11)
          return
      endif       

      goto 10

 900  istatus = 0
      close(11)
      write(6,*)'read_radar_info: Radar data not found'
      return

 999  istatus = 1

      return
      end
