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

        subroutine get_laps_domain(ni,nj,grid_fnam,lat,lon,topo,istatus)

!       1992 Steve Albers
c
c  get maps grid information
c
c   lat() contains grid latitudes (degrees north)
c   lon() contains grid longitudes (degrees east; negative=west long.)
c   topo() contains grid elevations (m)
c
        integer*4 ni,nj              ! Input

        real*4 lat(ni,nj)            ! Output
        real*4 lon(ni,nj)            ! Output
        real*4 topo(ni,nj)           ! Output

        character*(*) grid_fnam      ! Input

        character*3 var
        character*50  directory
        character*31  ext
        character*10  units
        character*125 comment

        character*80 grid_fnam_common
        common / grid_fnam_cmn / grid_fnam_common
c
        write(6,*)'    Reading in lat/lon/topo ',grid_fnam

        grid_fnam_common = grid_fnam  ! Used in get_directory to modify
                                      ! extension based on the grid domain

        ext = grid_fnam

!       Get the location of the static grid directory
        call get_directory(ext,directory,len_dir)

!       directory = ''

        var = 'LAT'
        call rd_laps_static(directory,ext,ni,nj,1,var,units,comment
     1                                 ,lat,grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LAT field'
            return
        endif

        var = 'LON'
        call rd_laps_static(directory,ext,ni,nj,1,var,units,comment
     1                                  ,lon,grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LON field'
            return
        endif

        var = 'AVG'
        call rd_laps_static(directory,ext,ni,nj,1,var,units,comment
     1                                  ,topo,grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading AVG (topo) field'
            return
        endif

c       write(6,*)' LAT/LON Corner > ',lat(   1,   1),lon(   1,   1)
c       write(6,*)' LAT/LON Corner > ',lat(   1,nj),lon(   1,nj)
c       write(6,*)' LAT/LON Corner > ',lat(ni,   1),lon(ni,   1)
c       write(6,*)' LAT/LON Corner > ',lat(ni,nj),lon(ni,nj)

        call get_laps_config(grid_fnam,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error in get_laps_config'
            return
        endif

        call check_domain(lat,lon,ni,nj,grid_spacing_m,5,istat_chk)
        if(istat_chk .ne. 1)then
            write(6,*)' Warning or Error in check_domain'
        endif

        return

        end


        subroutine get_domain_laps(ni,nj,grid_fnam,lat,lon,topo
     1                             ,grid_spacing_m,istatus)

!       1994 Steve Albers
c
c  get maps grid information
c
c   lat() contains grid latitudes (degrees north)
c   lon() contains grid longitudes (degrees east; negative=west long.)
c   topo() contains grid elevations (m)
c
        integer*4 ni,nj              ! Input

        real*4 lat(ni,nj)            ! Output
        real*4 lon(ni,nj)            ! Output
        real*4 topo(ni,nj)           ! Output

        character*(*) grid_fnam      ! Input

        character*3 var
        character*50  directory
        character*31  ext
        character*10  units
        character*125 comment

        character*80 grid_fnam_common
        common / grid_fnam_cmn / grid_fnam_common

c
        write(6,*)'    Reading in lat/lon/topo ',grid_fnam

        grid_fnam_common = grid_fnam  ! Used in get_directory to modify
                                      ! extension based on the grid domain

        ext = grid_fnam

!       Get the location of the static grid directory
        call get_directory(ext,directory,len_dir)

!       directory = ''

        var = 'LAT'
        call rd_laps_static(directory,ext,ni,nj,1,var,units,comment
     1                                 ,lat,grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LAT field'
            return
        endif

        var = 'LON'
        call rd_laps_static(directory,ext,ni,nj,1,var,units,comment
     1                                  ,lon,grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LON field'
            return
        endif

        var = 'AVG'
        call rd_laps_static(directory,ext,ni,nj,1,var,units,comment
     1                                  ,topo,grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading AVG (topo) field'
            return
        endif

c       write(6,*)' LAT/LON Corner > ',lat(   1,   1),lon(   1,   1)
c       write(6,*)' LAT/LON Corner > ',lat(   1,nj),lon(   1,nj)
c       write(6,*)' LAT/LON Corner > ',lat(ni,   1),lon(ni,   1)
c       write(6,*)' LAT/LON Corner > ',lat(ni,nj),lon(ni,nj)

        call get_laps_config(grid_fnam,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error in get_laps_config'
            return
        endif

        call check_domain(lat,lon,ni,nj,grid_spacing_m,5,istat_chk)
        if(istat_chk .ne. 1)then
            write(6,*)' Warning or Error in check_domain'
        endif

        return

        end



        subroutine get_laps_domain_95(ni,nj,grid_fnam,lat,lon,topo
     1            ,rlaps_land_frac,grid_spacing_m,istatus)

!       1994 Steve Albers
c
c  get maps grid information
c
c   lat() contains grid latitudes (degrees north)
c   lon() contains grid longitudes (degrees east; negative=west long.)
c   topo() contains grid elevations (m)
c
        integer*4 ni,nj              ! Input

        real*4 lat(ni,nj)            ! Output
        real*4 lon(ni,nj)            ! Output
        real*4 topo(ni,nj)           ! Output
        real*4 rlaps_land_frac(ni,nj) ! Output

        character*(*) grid_fnam      ! Input

        character*3 var
        character*50  directory
        character*31  ext
        character*10  units
        character*125 comment

        character*80 grid_fnam_common
        common / grid_fnam_cmn / grid_fnam_common
c
        write(6,*)'    Reading in lat/lon/topo/land frac ',grid_fnam

        grid_fnam_common = grid_fnam  ! Used in get_directory to modify
                                      ! extension based on the grid domain

        ext = grid_fnam

!       Get the location of the static grid directory
        call get_directory(ext,directory,len_dir)

!       directory = ''

        var = 'LAT'
        call rd_laps_static(directory,ext,ni,nj,1,var,units,comment
     1                                 ,lat,grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LAT field'
            return
        endif

        var = 'LON'
        call rd_laps_static(directory,ext,ni,nj,1,var,units,comment
     1                                  ,lon,grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LON field'
            return
        endif

        var = 'AVG'
        call rd_laps_static(directory,ext,ni,nj,1,var,units,comment
     1                                  ,topo,grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading AVG (topo) field'
            return
        endif

        var = 'LDF'
        call rd_laps_static(directory,ext,ni,nj,1,var,units,comment
     1                 ,rlaps_land_frac,grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LDF (land fraction) field'
            return
        endif

c       write(6,*)' LAT/LON Corner > ',lat(   1,   1),lon(   1,   1)
c       write(6,*)' LAT/LON Corner > ',lat(   1,nj),lon(   1,nj)
c       write(6,*)' LAT/LON Corner > ',lat(ni,   1),lon(ni,   1)
c       write(6,*)' LAT/LON Corner > ',lat(ni,nj),lon(ni,nj)


        do i = 1,ni
        do j = 1,nj
            if(rlaps_land_frac(i,j) .le. 0.5)then
                rlaps_land_frac(i,j) = 0.           ! Water
            else
                rlaps_land_frac(i,j) = 1.           ! Land
            endif
        enddo
        enddo

        call get_laps_config(grid_fnam,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error in get_laps_config'
            return
        endif

        call check_domain(lat,lon,ni,nj,grid_spacing_m,5,istat_chk)
        if(istat_chk .ne. 1)then
            write(6,*)' Warning or Error in check_domain'
        endif

        return
        end


        subroutine get_laps_config(grid_fnam,istatus)

!       1992 Steve Albers

!       Read in parameters from parameter file. The first time the routine
!       is called, parameters are read in from the namelist. Subsequent
!       calls will return these parameters (now stored in common), without
!       rereading the namelist. This is more efficient if the routine is called
!       many times.

        character*(*) grid_fnam   ! Input (Warning: trailing blanks won't work)
        character*150  directory
        character*31  ext
        character*8 a8
        character*200 tempchar

        integer ipass
        data ipass/0/
        save ipass

        character*80 grid_fnam_common
        common / grid_fnam_cmn / grid_fnam_common

        include 'lapsparms.cmn'

        NAMELIST /lapsparms_NL/ iflag_lapsparms_cmn
     1  ,PRESSURE_BOTTOM_L,PRESSURE_INTERVAL_L,PRESSURE_0_L
     1  ,vertical_grid,nk_laps,standard_latitude,standard_latitude2
     1  ,standard_longitude,NX_L_CMN, NY_L_CMN, I_PERIMETER_CMN
     1  ,c50_lowres_directory,c6_maproj
     1  ,l_highres,l_pad1,l_pad2,l_pad3
     1  ,grid_spacing_m_cmn,grid_cen_lat_cmn,grid_cen_lon_cmn
     1  ,laps_cycle_time_cmn
     1  ,radarext_3d_cmn,radarext_3d_accum_cmn
     1  ,path_to_raw_pirep_cmn
     1  ,path_to_raw_rass_cmn,path_to_raw_profiler_cmn
     1  ,path_to_raw_blprass_cmn,path_to_raw_blpprofiler_cmn
     1  ,path_to_wsi_2d_radar_cmn,path_to_wsi_3d_radar_cmn
     1  ,path_to_qc_acars_cmn,path_to_raw_raob_cmn
     1  ,path_to_metar_data_cmn,path_to_local_data_cmn
     1  ,i2_missing_data_cmn, r_missing_data_cmn, MAX_RADARS_CMN
     1  ,ref_base_cmn,ref_base_useable_cmn,maxstns_cmn,N_PIREP_CMN
     1  ,vert_rad_meso_cmn,vert_rad_sao_cmn
     1  ,vert_rad_pirep_cmn,vert_rad_prof_cmn     
     1  ,silavwt_parm_cmn,toptwvl_parm_cmn,c8_project_common
     1  ,maxstations_cmn,maxobs_cmn
     1  ,c_raddat_type, c80_description, path_to_topt30s
     1  ,path_to_topt10m, path_to_pctl10m


        if(ipass.eq.1 .and. iflag_lapsparms_cmn .eq. 1) then ! Data already read in
!          print *, 'It works'
           goto999                                           ! Normal Return
        endif

!       While we are here, let's put the grid name into the common area
        grid_fnam_common = grid_fnam  ! Used in get_directory to modify
                                      ! extension based on the grid domain

!       Get the location of the parameter directory
        ext = grid_fnam
        call get_directory(ext,directory,len_dir)
        if(directory(len_dir:len_dir).ne.'/') then
          tempchar = directory(1:len_dir)//'/'//grid_fnam//'.parms'
        else
          tempchar = directory(1:len_dir)//grid_fnam//'.parms'
        endif
        call s_len(tempchar,len_dir)
 
        open(92,file=tempchar(1:len_dir),status='old',err=900)

        read(92,lapsparms_nl,err=910)

        if(iflag_lapsparms_cmn .ne. 1)then                      ! Error Return
            goto910                                    

        else                                                    ! Normal Return
            PRESSURE_0_L = PRESSURE_BOTTOM_L + PRESSURE_INTERVAL_L
            write(6,*)' get_laps_config - parameters read in OK'
            goto999                                     

        endif


900     write(6,*)                                              ! Error Return
     1  ' Open error in get_laps_config, parameter file not found'
        write(6,*)tempchar
        iflag_lapsparms_cmn = 0
        istatus = 0
        close(92)
        return

910     write(6,*)' Read error in get_laps_config'              ! Error Return
        write(6,*)' Check runtime parameter file ',tempchar
        write(6,lapsparms_nl)
        close(92)
        iflag_lapsparms_cmn = 0
        istatus = 0
        return

920     write(6,*)' Read error in get_laps_config'              ! Error Return
        write(6,*)' Truncated runtime parameter file ',tempchar
        close(92)
        iflag_lapsparms_cmn = 0
        istatus = 0
        return

999     close(92)                                               ! Normal Return

        ipass = 1
        istatus = 1
        return

        end


      subroutine get_standard_longitude(std_lon,istatus)

      include 'lapsparms.cmn' ! standard_longitude

!     This routine accesses the standard_longitude variable from the
!     .parms file via the common block. Note the variable name in the
!     argument list may be different in the calling routine

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' ERROR, get_laps_config not called'
          istatus = 0
          return
!         stop
      endif

      std_lon = standard_longitude

      istatus = 1
      return
      end


      subroutine get_grid_spacing(grid_spacing_m,istatus)

      include 'lapsparms.cmn' ! grid_spacing_m_cmn

!     This routine accesses the standard_longitude variable from the
!     .parms file via the common block.

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' ERROR, get_laps_config not called'
          istatus = 0
          return
!         stop
      endif

      grid_spacing_m = grid_spacing_m_cmn

      istatus = 1
      return
      end


      subroutine get_grid_center(grid_cen_lat,grid_cen_lon,istatus)

      include 'lapsparms.cmn' ! grid_spacing_m_cmn

!     This routine accesses the standard_longitude variable from the
!     .parms file via the common block.

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' ERROR, get_laps_config not called'
          istatus = 0
          return
!         stop
      endif

      grid_cen_lat = grid_cen_lat_cmn
      grid_cen_lon = grid_cen_lon_cmn

      istatus = 1
      return
      end

      subroutine get_standard_latitude(std_lat,istatus)

      include 'lapsparms.cmn' ! standard_latitude

!     This routine accesses the standard_latitude variable from the
!     .parms file via the common block. Note the variable name in the
!     argument list may be different in the calling routine

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' ERROR, get_laps_config not called'
          istatus = 0
          return
!         stop
      endif

      std_lat = standard_latitude

      istatus = 1
      return
      end


      subroutine get_standard_latitudes(std_lat1,std_lat2,istatus)

      include 'lapsparms.cmn' ! standard_latitude, standard_latitude2

!     This routine accesses the standard_latitude variables from the
!     .parms file via the common block. Note the variable name in the
!     argument list may be different in the calling routine

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' ERROR, get_laps_config not called'
          istatus = 0
          return
!         stop
      endif

      std_lat1 = standard_latitude
      std_lat2 = standard_latitude2

      istatus = 1
      return
      end


      subroutine get_maxstns(maxstns,istatus)

      include 'lapsparms.cmn' ! maxstns_cmn

!     This routine accesses the maxstns_cmn variable from the
!     .parms file via the common block. Note the variable name in the
!     argument list may be different in the calling routine

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' ERROR, get_laps_config not called'
          istatus = 0
          return
!         stop
      endif

      maxstns = maxstns_cmn

      istatus = 1
      return
      end


      subroutine get_c8_project(c8_project,istatus)

      include 'lapsparms.cmn' ! c8_project

      character*8 c8_project

!     This routine accesses the c8_project_common variable from the
!     .parms file via the common block. Note the variable name in the
!     argument list may be different in the calling routine

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' ERROR, get_laps_config not called'
          istatus = 0
          return
!         stop
      endif

      c8_project = c8_project_common

      istatus = 1
      return
      end


      subroutine get_c6_maproj(c6_maproj_ret,istatus)

      include 'lapsparms.cmn' ! c6_maproj

      character*6 c6_maproj_ret

!     This routine accesses the c6_maproj variable from the
!     .parms file via the common block. Note the variable name in the
!     argument list may be different in the calling routine

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' ERROR, get_laps_config not called'
          istatus = 0
          return
!         stop
      endif

      c6_maproj_ret = c6_maproj

      istatus = 1
      return
      end


      subroutine get_c80_description(c80_description_ret,istatus)

      include 'lapsparms.cmn' ! c80_description

      character*80 c80_description_ret
!     This routine accesses the c80_description variable from the
!     .parms file via the common block. Note the variable name in the
!     argument list may be different in the calling routine

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' ERROR, get_laps_config not called'
          istatus = 0
          return
!         stop
      endif

      c80_description_ret = c80_description

      istatus = 1
      return
      end


      subroutine get_laps_dimensions(nk,istatus)

      include 'lapsparms.cmn'              ! nk_laps

!     This routine accesses the nk variable from the
!     .parms file via the common block. Note the variable name in the
!     argument list may be different in the calling routine

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' get_laps_dimensions: calling get_laps_config'

          call get_laps_config('nest7grid',istatus)
          if(istatus .ne. 1 .or. iflag_lapsparms_cmn .ne. 1)then
              write(6,*)' Error detected in calling get_laps_config'
              istatus = 0
              return
          else
              write(6,*)' Success in calling get_laps_config'
          endif

      endif

      nk = nk_laps

      istatus = 1
      return
      end


      subroutine get_laps_cycle_time(laps_cycle_time,istatus)

      include 'lapsparms.cmn'              ! laps_cycle_time_cmn

!     This routine accesses the laps_cycle_time variable from the
!     .parms file via the common block. Note the variable name in the
!     argument list may be different in the calling routine

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' get_laps_cycle_time: calling get_laps_config'

          call get_laps_config('nest7grid',istatus)
          if(istatus .ne. 1 .or. iflag_lapsparms_cmn .ne. 1)then
              write(6,*)' Error detected in calling get_laps_config'
              istatus = 0
              return
          else
              write(6,*)' Success in calling get_laps_config'
          endif

      endif

      laps_cycle_time = laps_cycle_time_cmn

      istatus = 1
      return
      end


      subroutine get_grid_dim_xy(NX_L,NY_L,istatus)

      include 'lapsparms.cmn'              ! NX_L_CMN, NY_L_CMN

!     This routine accesses the NX_L and NY_L variables from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' get_grid_dim_xy: calling get_laps_config'

          call get_laps_config('nest7grid',istatus)
          if(istatus .ne. 1 .or. iflag_lapsparms_cmn .ne. 1)then
              write(6,*)' Error detected in calling get_laps_config'
              istatus = 0
              return
          else
              write(6,*)' Success in calling get_laps_config'
          endif

      endif

      NX_L = NX_L_CMN
      NY_L = NY_L_CMN

      istatus = 1
      return
      end

      subroutine get_topo_parms(silavwt_parm,toptwvl_parm,istatus)

      include 'lapsparms.cmn' ! silavwt_cmn, toptwvl_cmn

!     This routine accesses the silavwt_parm and toptwvl_parm
!     variables from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' get_topo_parms: calling get_laps_config'

          call get_laps_config('nest7grid',istatus)
          if(istatus .ne. 1 .or. iflag_lapsparms_cmn .ne. 1)then
              write(6,*)' Error detected in calling get_laps_config'
              istatus = 0
              return
          else
              write(6,*)' Success in calling get_laps_config'
          endif

      endif

      silavwt_parm = silavwt_parm_cmn
      toptwvl_parm = toptwvl_parm_cmn

      istatus = 1
      return
      end

      subroutine get_meso_sao_pirep(n_meso,n_sao,n_pirep,istatus)

      include 'lapsparms.cmn' ! maxstns_cmn, n_pirep_cmn

!     This routine accesses the maxstns and n_pirep
!     variables from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' ERROR, get_laps_config not called'
          istatus = 0
          return
!         stop
      endif

      n_meso  = maxstns_cmn
      n_sao   = maxstns_cmn
      n_pirep = n_pirep_cmn

      istatus = 1
      return
      end

      subroutine get_max_radars (max_radars, istatus)

      include 'lapsparms.cmn' ! max_radars_cmn

!     This routine accesses the max_radars variable from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' ERROR, get_laps_config not called'
          istatus = 0
          return
!         stop
      endif

      max_radars = max_radars_cmn

      istatus = 1
      return
      end

      subroutine get_max_stations (maxstns, istatus)

      include 'lapsparms.cmn' ! maxstns_cmn

!     This routine accesses the maxstns variable from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' ERROR, get_laps_config not called'
          istatus = 0
          return
!         stop
      endif

      maxstns = maxstns_cmn

      istatus = 1
      return
      end

      subroutine get_vert_rads (vert_rad_pirep,
     1                          vert_rad_sao,
     1                          vert_rad_meso,
     1                          vert_rad_prof,
     1                          istatus)

      include 'lapsparms.cmn' ! vert_rad_pirep_cmn, etc.

      integer*4 vert_rad_pirep
      integer*4 vert_rad_sao
      integer*4 vert_rad_meso
      integer*4 vert_rad_prof

!     This routine accesses the vert_rad_pirep, etc., variables from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' ERROR, get_laps_config not called'
          istatus = 0
          return
!         stop
      endif

      vert_rad_pirep = vert_rad_pirep_cmn
      vert_rad_sao = vert_rad_sao_cmn
      vert_rad_meso = vert_rad_meso_cmn
      vert_rad_prof = vert_rad_prof_cmn

      istatus = 1
      return
      end

      subroutine get_r_missing_data(r_missing_data, istatus)

      include 'lapsparms.cmn' ! r_missing_data_cmn

!     This routine accesses the r_missing_data variable from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' ERROR, get_laps_config not called'
          istatus = 0
          return
!         stop
      endif

      r_missing_data = r_missing_data_cmn

      istatus = 1
      return
      end

      subroutine get_i2_missing_data(i2_missing_data, istatus)

      include 'lapsparms.cmn' ! i2_missing_data_cmn

!     This routine accesses the i2_missing_data variable from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' ERROR, get_laps_config not called'
          istatus = 0
          return
!         stop
      endif

      i2_missing_data = i2_missing_data_cmn

      istatus = 1
      return
      end


      subroutine get_i_perimeter(i_perimeter, istatus)

      include 'lapsparms.cmn' ! i_perimeter_cmn

!     This routine accesses the 'i_perimeter' variable from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' ERROR, get_laps_config not called'
          istatus = 0
          return
!         stop
      endif

      i_perimeter = i_perimeter_cmn

      istatus = 1
      return
      end


      subroutine get_ref_base(ref_base, istatus)

      include 'lapsparms.cmn' ! ref_base_cmn

!     This routine accesses the 'ref_base' variable from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' ERROR, get_laps_config not called'
          istatus = 0
          return
!         stop
      endif

      ref_base = ref_base_cmn

      istatus = 1
      return
      end

      subroutine get_ref_base_useable(ref_base_useable, istatus)

      include 'lapsparms.cmn' ! ref_base_useable

!     This routine accesses the 'ref_base_useable' variable from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' ERROR, get_laps_config not called'
          istatus = 0
          return
!         stop
      endif

      ref_base_useable = ref_base_useable_cmn

      istatus = 1
      return
      end

      subroutine get_background_info(len,bgpaths,bgmodels
     +     ,oldest_forecast,max_forecast_delta,use_analysis)
      implicit none
      integer maxbgmodels,len
      parameter (maxbgmodels=10)
      character*150 nest7grid
      character*150 bgpaths(maxbgmodels)
      integer bgmodels(maxbgmodels), len_dir
      integer oldest_forecast, max_forecast_delta
      logical use_analysis
      namelist /background_nl/bgpaths,bgmodels
     +         ,oldest_forecast,max_forecast_delta,use_analysis   

      max_forecast_delta=6
      oldest_forecast=18
      use_analysis=.false.
      call get_directory('nest7grid',nest7grid,len_dir)
      if(nest7grid(len_dir:len_dir).ne.'/') then
        len_dir=len_dir+1
        nest7grid(len_dir:len_dir)='/'
      endif
      nest7grid = nest7grid(1:len_dir)//'background.nl'
      
      open(1,file=nest7grid(1:len_dir+13),status='old',err=900)
      read(1,background_nl,err=901)
      close(1)
      return
 900  print*,'error opening file ',nest7grid
      stop
 901  print*,'error reading background_nl in ',nest7grid
      write(*,background_nl)
      stop
      end
c
      subroutine get_sat_sounder_info(n_sndr,c_sndr_id,
     +n_sndr_channels,path_to_sat_sounder,n_elems,n_lines,
     +channel_wavelength_u,imsng_sndr_pix,istatus)

      integer maxsndr
      integer maxch
      parameter (maxsndr=4,maxch=19)

      integer len_dir
      integer n_sndr
      integer n_sndr_channels
      integer n_elems(maxsndr)
      integer n_lines(maxsndr)
      integer imsng_sndr_pix
      integer istatus
      character*6   c_sndr_id(maxsndr)
      character*150 nest7grid
      character*200 path_to_sat_sounder(maxsndr)

      real*4 channel_wavelength_u(maxch,maxsndr)

c-----------------------------------------------------------------------
      namelist /satellite_sounder_nl/ n_sndr,c_sndr_id,path_to_sat_sound
     +er,n_elems,n_lines,n_sndr_channels,channel_wavelength_u,imsng_sndr
     +_pix

      call get_directory('nest7grid',nest7grid,len_dir)

      nest7grid = nest7grid(1:len_dir)//'/sat_sounder.nl'

      open(1,file=nest7grid,status='old',err=900)
      read(1,satellite_sounder_nl,err=901)
      close(1)

      istatus = 1
      return
 900  print*,'error opening file ',nest7grid
      stop
 901  print*,'error reading satellite_sounder_nl in ',nest7grid
      write(*,satellite_sounder_nl)
      stop 
      end 
c
c-----------------------------------------------------------------------
c
      subroutine config_satellite_lvd(istatus)

      character nest7grid*150
      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'
      include 'satellite_namelist_lvd.cmn'

      istatus = 0

      call get_directory('nest7grid',nest7grid,len_dir)

      nest7grid = nest7grid(1:len_dir)//'/satellite_lvd.nl'

      open(1,file=nest7grid,status='old',err=900)
      read(1,satellite_lvd_nl,err=901)
      close(1)
      iflag_lvd_common=1
      istatus = 1
      return

 900  print*,'error opening file ',nest7grid
      return
 901  print*,'error reading satellite_nl in ',nest7grid
      write(*,satellite_lvd_nl)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine config_ingest_rrv_parms(nxv01,nyv01,nzv01,istatus)

      character nest7grid*150
      include 'ingest_rrv_dims.inc'
      include 'ingest_rrv_common.inc'

      namelist /ingest_rrv_nl/ path_to_raw_rrv,nxv01,nyv01,nzv01,
     +dxv01,dyv01

      istatus = 0

      call get_directory('nest7grid',nest7grid,len_dir)

      nest7grid = nest7grid(1:len_dir)//'/ingest_rrv.nl'

      open(1,file=nest7grid,status='old',err=900)
      read(1,ingest_rrv_nl,err=901)
      close(1)
      istatus = 1
      return

 900  print*,'error opening file ',nest7grid
      return
 901  print*,'error reading ingest_rrv_nl in ',nest7grid
      write(*,ingest_rrv_nl)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine get_wsi_parms_vrc(irad,lines,elems,
     +dx,dy,rla1,rlo1,rla2,rlo2,rlat,rlon,rlatin,istatus)

      integer nrad_types
      parameter (nrad_types=2)

      integer nlines(nrad_types)
      integer nelems(nrad_types)
      integer irad
      integer lines,elems
      integer istatus

      real*4  resx(nrad_types)
      real*4  resy(nrad_types)
      real*4  radla1(nrad_types)
      real*4  radlo1(nrad_types)
      real*4  radla2(nrad_types)
      real*4  radlo2(nrad_types)
      real*4  radlat(nrad_types)
      real*4  radlon(nrad_types)
      real*4  radlatin(nrad_types)

      real*4  rla1,rlo1,rlat,rlon,rlatin
      real*4  dx,dy

      character nest7grid*150

      namelist /vrc_nl/nlines,nelems,resx,resy,radla1,
     +radlo1,radla2,radlo2,radlat,radlon,radlatin

      istatus = 0

      call get_directory('nest7grid',nest7grid,len_dir)

      nest7grid = nest7grid(1:len_dir)//'/vrc.nl'

      open(1,file=nest7grid,status='old',err=900)
      read(1,vrc_nl,err=901)

      lines=nlines(irad)
      elems=nelems(irad)
      dx=resx(irad)
      dy=resy(irad)
      rla1=radla1(irad)
      rlo1=radlo1(irad)
      rla2=radla2(irad)
      rlo2=radlo2(irad)
      rlat=radlat(irad)
      rlon=radlon(irad)
      rlatin=radlatin(irad)

      close(1)
      istatus = 1
      return

 900  print*,'error opening file ',nest7grid
      return
 901  print*,'error reading vrc_nl in ',nest7grid
      write(*,vrc_nl)
      return
      end

      function ltest_vertical_grid(c_vertical_grid)

!     The input c_vertical_grid is a string you are testing
!     The output (.true. OR .false.) states whether the input vertical_grid 
!     name equals the vertical_grid in the Namelist

      include 'lapsparms.cmn' ! vertical_grid

      logical ltest_vertical_grid
      character*(*) c_vertical_grid

      call get_laps_config('nest7grid',istatus)

      if(istatus .ne. 1)then
          write(6,*)' ltest_vertical_grid: Error detected in '
     1             ,'calling get_laps_config'
          ltest_vertical_grid = .false.
          return
      endif

      call s_len(vertical_grid,len_cmn)
      call s_len(c_vertical_grid,len_in)

      if(vertical_grid(1:len_cmn) .eq. c_vertical_grid(1:len_in))then
          ltest_vertical_grid = .true.
      else
          ltest_vertical_grid = .false.
      endif

      return
      end
