cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   
cdis
      subroutine lut_gen(c4_radarname,rlat_radar,rlon_radar ! I
     :                  ,rheight_radar,ioffset,joffset      ! I
     :                  ,NX_L,NY_L,NZ_L)                    ! I
c
c     PURPOSE:
c        Generate look-up tables for radar remapping.
c
      include 'trigd.inc'
      implicit none
c
      include 'remap_constants.dat'
      include 'remap.cmn'
c
c     Variables from LAPS domain file
c

      integer NX_L,NY_L,NZ_L
      real lat(NX_L,NY_L)
      real lon(NX_L,NY_L)
      real topo(NX_L,NY_L)
c
c     Functions
c
      real height_to_zcoord
      integer ishow_timer
c
c     Misc interval variables
c
      integer i,j,igate_lut,iaz,ielev,iran,iz_grid,istatus,len_dir
      integer I4_elapsed, lenr, irecl_lutge, irecl_lutar
      integer ioffset,joffset,io,jo
      real rlat_grid,rlon_grid,height_grid
      real rlat_radar,rlon_radar,rheight_radar
      real elev,elev_deg,coselev,azimuth,azi_deg,azimuth_interval
      real slant_range,sl_range_m,ri,rj,dbz,z,grid_spacing_cen_m
      real difflat,difflon
      character*4 c4_radarname
      character*150 static_dir,filename
      character*3 ext
      logical l_readwrite_lut 
      data l_readwrite_lut /.true./
c
c     Fill arrays with initial values
c
      DO 50 i = 1,10
        i4time_old(i) = 0
        n_ref_obs_old(i) = 99999
   50 CONTINUE
c
      write(6,801)NX_L,NY_L
  801 format('REMAP > LUT_GEN - Getting LAPS Domain: ',2i5)
c
      call get_domain_laps(NX_L,NY_L,'nest7grid',
     :                     lat,lon,topo,grid_spacing_cen_m,istatus)    
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

      if(range_interval .gt. grid_spacing_cen_m)then
          write(6,*)' WARNING: Range interval is > grid spacing of'
     1             ,grid_spacing_cen_m, 'leading to inaccurate remap'
      endif

      write(6,812) gate_spacing_m
  812 format(' Effective gate spacing for radar ',F10.2)

      write(6,*)' Radar Name  ',c4_radarname
      write(6,*)' Radar Coords',rlat_radar,rlon_radar,rheight_radar
      call latlon_to_rlapsgrid(rlat_radar,rlon_radar,lat,lon,
     1                         NX_L,NY_L,ri,rj,istatus)
      write(6,*)' LAPS Grid Coords of radar ',ri,rj


      call s_len(c4_radarname,lenr)
      if(lenr .eq. 0)then
          write(6,*)' Error in lut_gen - no radarname '
          stop
      endif

      if(abs(rlat_radar) .gt. 90.)then
          write(6,*)' Error in lut_gen - invalid radar lat '
          stop
      endif
      
      ext = 'dat'
      call get_directory(ext,static_dir,len_dir)

c     Put the coords into common so remap_process can access them
      rlat_radar_cmn = rlat_radar
      rlon_radar_cmn = rlon_radar
      rheight_radar_cmn = rheight_radar
      c4_radarname_cmn = c4_radarname      

      irecl_lutge = LUT_GATES * (LUT_ELEVS+1) * 4
      irecl_lutar = (LUT_AZIMUTHS+1) * (LUT_RANGES+1) * 2 * 4
c
c     Generate Gate/Elev to Projran lut
c
      if(l_readwrite_lut)then ! Try to read lut
          filename = static_dir(1:len_dir)//'vxx/'
     1             //'gate_elev_to_projran_lut.'//c4_radarname
          write(6,*)' Reading file: ',filename
          open(11,file=filename,form='unformatted'
     1         ,access='direct',recl=irecl_lutge,status='old',err=90)       
          read(11,err=90)gate_elev_to_projran_lut
          close(11)
          goto 110
90        write(6,*)' Generating LUT - no valid file exists'
          close(11)
      endif

!     Calculate lut
      write(6,820) lut_gates
  820 format(' REMAP > Building Gate/Elev to Projran lut, gates =',
     :        I12)

      DO ielev = 0,lut_elevs

        elev_deg = min_elev + (ielev * elev_interval)
        coselev = cosd(elev_deg)

        DO igate_lut = 1,lut_gates

          sl_range_m = igate_lut * gate_spacing_m * gate_interval
          iran = sl_range_m * coselev / range_interval
          gate_elev_to_projran_lut(igate_lut,ielev) = iran

        ENDDO
      ENDDO

      if(l_readwrite_lut)then ! Write lut
          open(12,file=filename,form='unformatted',access='direct'
     1               ,recl=irecl_lutge,status='unknown')
          write(12)gate_elev_to_projran_lut
          close(12)
      endif
  110 continue 

      I4_elapsed = ishow_timer()
c
c     Generate Gate/Elev to Z lut
c
      if(l_readwrite_lut)then ! Try to read lut
          filename = static_dir(1:len_dir)//'vxx/'
     1             //'gate_elev_to_z_lut.'//c4_radarname
          write(6,*)' Reading file: ',filename
          open(11,file=filename,form='unformatted',access='direct'
     1                           ,recl=irecl_lutge,status='old',err=190)
          read(11,err=190)gate_elev_to_z_lut
          close(11)
          goto 210
190       write(6,*)' Generating LUT - no valid file exists'
          close(11)
      endif

!     Calculate lut
      write(6,830)
  830 format(' REMAP > Building Gate/Elev to Z lut')
c
      DO 200 ielev = 0,lut_elevs
 
        elev_deg = min_elev + (ielev * elev_interval)
 
        DO 180 igate_lut = 1,lut_gates

          sl_range_m = igate_lut * gate_spacing_m * GATE_INTERVAL
          azi_deg=0.
 
          call radar_to_latlon
     :         (rlat_grid,rlon_grid,height_grid
     :                   ,azi_deg,sl_range_m,elev_deg
     :                  ,rlat_radar,rlon_radar,rheight_radar)

!         This implies use of the standard atmosphere
          iz_grid = nint(height_to_zcoord(height_grid,istatus))
          gate_elev_to_z_lut(igate_lut,ielev) = min(iz_grid,NZ_L)

  180   CONTINUE
  200 CONTINUE

      if(l_readwrite_lut)then ! Write lut
          open(12,file=filename,form='unformatted',access='direct'
     1               ,recl=irecl_lutge,status='unknown')
          write(12)gate_elev_to_z_lut
          close(12)
      endif
  210 continue 

      I4_elapsed = ishow_timer()
c
c     Generate Az/Ran to i,j lut
c
      if(l_readwrite_lut)then ! Try to read lut
          filename = static_dir(1:len_dir)//'vxx/'
     1             //'azran_to_ijgrid_lut.'//c4_radarname
          write(6,*)' Reading file: ',filename
          open(11,file=filename,form='unformatted',access='direct'
     1                           ,recl=irecl_lutar,status='old',err=290)
          read(11,err=290)azran_to_igrid_lut,azran_to_jgrid_lut
          close(11)
          goto 310
290       write(6,*)' Generating LUT - no valid file exists'
          close(11)
      endif

!     Calculate lut
      write(6,840)
  840 format(' REMAP > Building Az/Ran to i,j lut')
c
      if(range_interval .le. 0.)then ! QC check
          write(6,*)' ERROR: range_interval = ',range_interval
          stop
      else
          write(6,*)' range_interval = ',range_interval
          azimuth_interval = 360. / float(lut_azimuths)
          write(6,*)' azimuth_interval = ',azimuth_interval
      endif

      DO 300 iran = 0,lut_ranges

        slant_range = iran*range_interval

        DO 280 iaz = 0,lut_azimuths

          azimuth = float(iaz) * azimuth_interval
          elev = 0.

          call radar_to_latlon(rlat_grid,rlon_grid,height_grid
     :                  ,azimuth,slant_range,elev
     :                  ,rlat_radar,rlon_radar,rheight_radar)

!         QC check
          difflat = rlat_grid - rlat_radar
          difflon = rlon_grid - rlon_radar
          if(slant_range .gt. 10000. .and. abs(difflat) .lt. .01
     1                               .and. abs(difflon) .lt. .01)then
              write(6,*)
     1           ' ERROR: QC check failed after radar_to_latlon call'
              write(6,*)'difflat,difflon,slant_range'
     1                  ,difflat,difflon,slant_range
              write(6,*)'azimuth,elev,rlat_radar,rlon_radar'
     1                  ,azimuth,elev,rlat_radar,rlon_radar
              
              stop
          endif

          call latlon_to_rlapsgrid(rlat_grid,rlon_grid,lat,lon,
     :                   NX_L,NY_L,ri,rj,istatus)

          i = nint(ri)
          j = nint(rj)

          io = i - ioffset
          jo = i - joffset

          IF (i.le.0.or.i.gt.NX_L.or.j.le.0.or.j.gt.NY_L) THEN

            azran_to_igrid_lut(iaz,iran) = 0
            azran_to_jgrid_lut(iaz,iran) = 0

          ELSE

            azran_to_igrid_lut(iaz,iran) = i ! io
            azran_to_jgrid_lut(iaz,iran) = j ! jo

          END IF

  280   CONTINUE
  300 CONTINUE

      if(l_readwrite_lut)then ! Write lut
          open(12,file=filename,form='unformatted',access='direct'
     1               ,recl=irecl_lutar,status='unknown')
          write(12)azran_to_igrid_lut,azran_to_jgrid_lut
          close(12)
      endif
  310 continue 

!     Debugging output
      write(6,*)' azran_to_[ij]grid_lut...'
      do iran = 0,lut_ranges,10
        slant_range = iran*range_interval
        iaz = 0
        write(6,*)iran,slant_range,azran_to_igrid_lut(iaz,iran),
     1                             azran_to_jgrid_lut(iaz,iran)
      enddo ! iran

      I4_elapsed = ishow_timer()

c     Generate DbZ Lookup Table (graduated for each tenth of a dbz)
      do i = -1000,+1000
          dbz = float(i) / 10.
          z = 10**(dbz / 10.)
          dbz_to_z_lut(i) = z
      enddo ! i

!     call lgate_lut_gen()

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


      subroutine lgate_lut_gen(gsp_ref_m_cmn,gsp_vel_m_cmn
     1                        ,n_ref_gates_cmn,n_vel_gates_cmn) 

!     This routine is currently called from 'lut_gen'. We may want to move
!     the call to 'ld_ray' (closer to where these luts are used) to help 
!     make things more dynamic and understandable.

      include 'remap_constants.dat'
      include 'remap.cmn'

      ratio_vel = gsp_vel_m_cmn / GATE_SPACING_M
      ratio_ref = gsp_ref_m_cmn / GATE_SPACING_M

c     These lookup tables flag which gates actually need processing

      lgate_ref_lut = .false.
      lgate_vel_lut = .false.

      DO igate_88d = 1,n_vel_gates_cmn
          i = nint(float(igate_88d) * ratio_vel)

          if(i .le. MAX_GATES .and. i .ge. INITIAL_VEL_GATE
     1 .and. i .le. 920                                    )then
              lgate_vel_lut(i) = .true.
          endif ! i
      ENDDO

      DO igate_88d = 1,n_ref_gates_cmn
          i = nint(float(igate_88d) * ratio_ref)

          if(i .le. MAX_GATES .and. i .ge. INITIAL_REF_GATE)then
              lgate_ref_lut(i) = .true.
          endif ! i
      ENDDO

      do i = 1,MAX_GATES
           
!         if(i .le. 920       .and. i .ge. INITIAL_VEL_GATE)then
!             lgate_vel_lut(i) = .true.
!         else
!             lgate_vel_lut(i) = .false.
!         endif

!         if(i .eq. (i/4) * 4 .and. i .ge. INITIAL_REF_GATE)then
!             lgate_ref_lut(i) = .true.
!         else
!             lgate_ref_lut(i) = .false.
!         endif

          if(lgate_vel_lut(i) .or. lgate_ref_lut(i))then
              lgate_lut(i) = .true.
          else
              lgate_lut(i) = .false.
          endif

      enddo ! i

      return
      end
