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
c
cdoc    Reads in lat/lon/topo/landfrac fields
c
c       1994 Steve Albers

        integer*4 ni,nj               ! Input

        real*4 lat(ni,nj)             ! Output
        real*4 lon(ni,nj)             ! Output
        real*4 topo(ni,nj)            ! Output
        real*4 rlaps_land_frac(ni,nj) ! Local

        character*(*) grid_fnam       ! Dummy

        call get_laps_domain_95(ni,nj,lat,lon,topo
     1            ,rlaps_land_frac,grid_spacing_cen_m,istatus)

        return
        end


        subroutine get_domain_laps(ni,nj,grid_fnam,lat,lon,topo
     1                             ,grid_spacing_cen_m,istatus)
c
cdoc    Reads in lat/lon/topo fields with grid spacing
c
c       1994 Steve Albers
c
        integer*4 ni,nj               ! Input

        real*4 lat(ni,nj)             ! Output
        real*4 lon(ni,nj)             ! Output
        real*4 topo(ni,nj)            ! Output
        real*4 rlaps_land_frac(ni,nj) ! Local

        character*(*) grid_fnam       ! Dummy

        call get_laps_domain_95(ni,nj,lat,lon,topo
     1            ,rlaps_land_frac,grid_spacing_cen_m,istatus)

        return
        end



        subroutine get_laps_domain_95(ni,nj,lat,lon,topo
     1            ,rlaps_land_frac,grid_spacing_cen_m,istatus)
c
cdoc    Reads in lat/lon/topo/landfrac fields with grid spacing
c
c       1994 Steve Albers

        integer*4 ni,nj               ! Input

        real*4 lat(ni,nj)             ! Output
        real*4 lon(ni,nj)             ! Output
        real*4 topo(ni,nj)            ! Output
        real*4 rlaps_land_frac(ni,nj) ! Output

        character*3   var
        character*150 directory
        character*31  ext
        character*10  units
        character*125 comment

        include 'grid_fname.cmn'

        write(6,*)'    Reading in lat/lon/topo/land frac '

        ext = grid_fnam_common

!       Get the location of the static grid directory
        call get_directory(ext,directory,len_dir)

!       directory = ''

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)return

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
       
        call array_minmax(topo,ni,nj,rmin,rmax,r_missing_data)
        if(rmin .lt. -1000. .or. rmax .gt. 9000.)then
            write(6,*)' Error, topo range out of bounds',rmin,max
            istatus = 0
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

        call check_domain(lat,lon,ni,nj,grid_spacing_m,5,istat_chk)
        if(istat_chk .ne. 1)then
            write(6,*)' Warning or Error in check_domain'
        endif

        icen = ni/2 + 1
        jcen = nj/2 + 1
        call get_grid_spacing_actual(lat(icen,jcen),lon(icen,jcen)
     1                              ,grid_spacing_cen_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error return from get_grid_spacing_actual'       
            return
        endif

        return
        end

        subroutine read_static_grid(ni,nj,var,static_grid,istatus)

cdoc    Reads an arbitrary static field given an input 'var' string

c       2000    Steve Albers

        integer*4 ni,nj                       ! Input
        character*(*) var                     ! Input
        real*4 static_grid(ni,nj)             ! Output

        character*150 directory
        character*31  ext
        character*10  units
        character*125 comment

        include 'grid_fname.cmn'

        write(6,*)'    Reading in static grid: ',var

        ext = grid_fnam_common

!       Get the location of the static grid directory
        call get_directory(ext,directory,len_dir)

        call rd_laps_static(directory,ext,ni,nj,1,var,units,comment
     1                     ,static_grid,grid_spacing_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading field ',var
            return
        endif

        return
        end

      subroutine force_get_laps_config(grid_fnam,istatus)
      include 'lapsparms.cmn'
      character*(*) grid_fnam
      integer istatus
      iflag_lapsparms_cmn=0
      call get_laps_config(grid_fnam,istatus)
      return
      end

      subroutine get_laps_config(grid_fnam,istatus)
c
c higher level wrapper to allow get_laps_config to use new WRF SI
c namelist with hardwire grid_fnam (as is the case in most of the repository).
c
      implicit none

      integer istatus
      character*(*) grid_fnam

      include 'grid_fname.cmn'

      grid_fnam_common = grid_fnam

      call get_config(istatus)
      if(istatus .ne. 1)then
         print*,'Error returned from get_config'
         return
      endif
      
      return
      end
c --------------------------------------------------------------------

      subroutine get_config(istatus)

        implicit none

        integer max_files
        parameter (max_files = 500)

        character*80  domain_name
        integer       len_grid_fnam

        include   'grid_fname.cmn'

        character c_filenames(max_files)*255
        character cfilespec*255
        character directory*200
        integer   i,j,jj
        integer   lend,lens,len_root
        integer   numoffiles
        integer   narg,iargc
        integer   istatus

        integer   idir_len
        data      idir_len/0/
        save      idir_len
c
c LAPS_DATA_ROOT is either environment variable or command line input.
c command line overrides env variable if it exists.
c
        istatus = 0
        if(idir_len.eq.0) then
           narg = iargc()
           if(narg.gt.0)then
              call getarg(narg,generic_data_root)
              call s_len(generic_data_root,len_root)
              if(generic_data_root(len_root-1:len_root).eq.'.x'.or.
     1generic_data_root(len_root-3:len_root).eq.'.exe')then 
                 print*,'Not a typical dataroot on command line'
                 print*,'Trying LAPS_DATA_ROOT environment variable'
                 call GETENV('LAPS_DATA_ROOT',generic_data_root)
              endif
           else
              call GETENV('LAPS_DATA_ROOT',generic_data_root)
           endif
           call get_directory_length(generic_data_root,len_root)
           if(len_root.eq.0)then
              call s_len(generic_data_root,len_root)
              if(len_root.eq.0)then
                 print*,'Use either command line or ENV variable for ',
     +'system DATA_ROOT'
                 stop
              endif
           endif

           call find_domain_name(generic_data_root,domain_name,istatus)
           if(istatus.ne.1)then
              print*,'Error returned from find_domain_name'
              return
           endif
           call s_len(grid_fnam_common,len_grid_fnam)
           idir_len = len_grid_fnam

        endif


c use the "grid namelist" to load lapsparms.cmn with appropriate values.
        if(grid_fnam_common(1:len_grid_fnam).eq.'nest7grid')then
           call get_laps_config_sub(grid_fnam_common,istatus)
           if(istatus.ne.1)then
              print*,'Error returned from get_laps_config_sub'
              return
           endif
        elseif(grid_fnam_common(1:len_grid_fnam).eq.'wrfsi')then
           call get_wrfsi_config(istatus)
           if(istatus.ne.1)then
              print*,'Error returned from get_wrfsi_config'
              return
           endif
        else
           print*,'Unknown grid filename spec'
           return
        endif

        istatus = 1

        return
        end

      subroutine get_laps_config_sub(grid_fnam,istatus)

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

        include 'lapsparms.cmn'
        include 'grid_fname.cmn'

        NAMELIST /lapsparms_NL/ iflag_lapsparms_cmn
     1  ,PRESSURE_BOTTOM_L,PRESSURE_INTERVAL_L
     1  ,nk_laps,standard_latitude,standard_latitude2       
     1  ,standard_longitude,NX_L_CMN, NY_L_CMN, I_PERIMETER_CMN
     1  ,l_highres
     1  ,grid_spacing_m_cmn,grid_cen_lat_cmn,grid_cen_lon_cmn
     1  ,laps_cycle_time_cmn, min_to_wait_for_metars_cmn
     1  ,i2_missing_data_cmn, r_missing_data_cmn, MAX_RADARS_CMN
     1  ,ref_base_cmn,ref_base_useable_cmn,maxstns_cmn,N_PIREP_CMN
     1  ,vert_rad_meso_cmn,vert_rad_sao_cmn
     1  ,vert_rad_pirep_cmn,vert_rad_prof_cmn     
     1  ,silavwt_parm_cmn,toptwvl_parm_cmn
     1  ,vertical_grid,c50_lowres_directory,c6_maproj
     1  ,radarext_3d_cmn,radarext_3d_accum_cmn
     1  ,path_to_raw_pirep_cmn
     1  ,path_to_raw_rass_cmn,path_to_raw_profiler_cmn
     1  ,path_to_raw_blprass_cmn,path_to_raw_blpprofiler_cmn
     1  ,path_to_wsi_2d_radar_cmn,path_to_wsi_3d_radar_cmn
     1  ,path_to_qc_acars_cmn
     1  ,c8_project_common
     1  ,c_raddat_type, c80_description, path_to_topt30s
     1  ,path_to_topt10m, path_to_pctl10m, path_to_soil2m
     1  ,fdda_model_source_cmn


        if(ipass.eq.1 .and. iflag_lapsparms_cmn .eq. 1) then ! Data already read in
!          print *, 'It works'
           goto999                                           ! Normal Return
        endif

!       While we are here, let's put the grid name into the common area
        grid_fnam_common = grid_fnam  ! Used in get_directory to modify
                                      ! extension based on the grid domain

!       Get the location of the parameter directory
        call s_len(grid_fnam,len_grid_fnam)
c        len_dir = index(grid_fnam_common,'/',.true.)
        call get_directory_length(grid_fnam_common,len_dir)

        if(len_dir.gt.0) then
           ext='nest7grid'
        else
           ext=grid_fnam
        endif

        call get_directory(ext,directory,len_dir)

! this is laps specific for nest7grid.parms ... the laps namelist file.
        call s_len(ext,len_ext)
        if(directory(len_dir:len_dir).ne.'/') then
          tempchar = directory(1:len_dir)//'/'//ext(1:len_ext)//'.parms'
        else
          tempchar = directory(1:len_dir)//ext(1:len_ext)//'.parms'
        endif
        call s_len(tempchar,len_dir)
 

        min_to_wait_for_metars_cmn=10
        open(92,file=tempchar(1:len_dir),status='old',err=900)

        read(92,lapsparms_nl,err=910)

        if(iflag_lapsparms_cmn .ne. 1)then                      ! Error Return
            goto910                                    

        else                                                    ! Normal Return
            PRESSURE_0_L = PRESSURE_BOTTOM_L + PRESSURE_INTERVAL_L
            write(6,*)tempchar(1:len_dir),
     +            ' get_laps_config - parameters read in OK'
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
        print*,'Here: dumping lapsparms_nl'
        print*
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
c ----------------------------------------------------------------------------
      subroutine find_domain_name(data_root, domain_name,
     +                            istatus)

        implicit none

        integer max_files
        parameter (max_files = 500)

        character*(*)  domain_name
        character*(*)  data_root
        integer   istatus

        character grid_fnam*80
        integer   len_grid_fnam, len_gfc
        integer   lend,lens,len_root
        character c_filenames(max_files)*255
        character cfilespec*255
        integer   numoffiles, idir_len, lgfc
        integer   i, j, jj
        logical   lfound_name

        include 'grid_fname.cmn'

        istatus = 0

        call s_len(grid_fnam_common,lgfc)
        if(lgfc.gt.0)then
           domain_name = grid_fnam_common
           istatus = 1
           return
        endif

c  determine the grid_fnam by searching the directory for
c  known namelists files used by this system
c (currently nest7grid.parms and wrfsi.nl)
c  assume the namelist file is always in DATA_ROOT/static.

        lfound_name = .false.
        i=1
        call s_len(data_root,len_root)
        cfilespec=data_root(1:len_root)//'/static/*'
        call get_file_names(cfilespec,numoffiles,c_filenames,
     +max_files,istatus)
        do while (i.le.numoffiles.and.(.not.lfound_name))
          call get_directory_length(c_filenames(i),lend)
          call s_len(c_filenames(i),lens)
          if(c_filenames(i)(lend+1:lens).eq.'nest7grid.parms'
     +  .or. c_filenames(i)(lend+1:lens).eq.'wrfsi.nl')then
              jj = 0
              do j=lend,lens
                 jj=jj+1
                 if(c_filenames(i)(j:j).eq.'.')then
                    grid_fnam = c_filenames(i)(lend+1:lend+jj-2)
                 endif
              enddo
              lfound_name = .true.
           endif
           i=i+1
        enddo

        call s_len(grid_fnam,len_grid_fnam)
        if (len_grid_fnam .gt. 0) then

          if ((grid_fnam(1:len_grid_fnam) .eq. 'nest7grid') .or.
     +        (grid_fnam(1:len_grid_fnam) .eq. 'wrfsi')) then
            istatus = 1
c     set the grid fname common block value
            len_gfc = len(grid_fnam_common)
            if (len_grid_fnam .le. len_gfc) then
              grid_fnam_common = grid_fnam
              domain_name = grid_fnam
            else
              print*, 'ERROR length of domain name exceeds size of'
     +               ,' grid_fnam_common: '
              print*, 'name= ', grid_fnam
              istatus = 0
            endif
          else
            print*, 'ERROR in getting domain name from ',
     +             data_root(1:len_root),'/static/'
            istatus = 0
          endif
        else
          print*, 'ERROR in getting domain name from ',
     +           data_root(1:len_root),'/static/'
          istatus = 0
        endif

        return
        end
c ---------------------------------------------------------------------
      subroutine get_standard_longitude(std_lon,istatus)

      include 'lapsparms.cmn' ! standard_longitude
      include 'grid_fname.cmn'! grid_fnam_common

!     This routine accesses the standard_longitude variable from the
!     .parms file via the common block. Note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
      endif

      std_lon = standard_longitude

      istatus = 1
      return
      end


      subroutine get_grid_spacing(grid_spacing_m,istatus)

      include 'lapsparms.cmn' ! grid_spacing_m_cmn
      include 'grid_fname.cmn'! grid_fnam_common

!     This routine accesses the standard_longitude variable from the
!     .parms file via the common block.

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
      endif

      grid_spacing_m = grid_spacing_m_cmn

      istatus = 1
      return
      end


      subroutine get_grid_center(grid_cen_lat,grid_cen_lon,istatus)

      include 'lapsparms.cmn' ! grid_spacing_m_cmn
      include 'grid_fname.cmn'! grid_fnam_common

!     This routine accesses the standard_longitude variable from the
!     .parms file via the common block.

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
      endif

      grid_cen_lat = grid_cen_lat_cmn
      grid_cen_lon = grid_cen_lon_cmn

      istatus = 1
      return
      end

      subroutine get_standard_latitude(std_lat,istatus)

      include 'lapsparms.cmn' ! standard_latitude
      include 'grid_fname.cmn'! grid_fnam_common

!     This routine accesses the standard_latitude variable from the
!     .parms file via the common block. Note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
      endif

      std_lat = standard_latitude

      istatus = 1
      return
      end


      subroutine get_standard_latitudes(std_lat1,std_lat2,istatus)

      include 'lapsparms.cmn' ! standard_latitude, standard_latitude2
      include 'grid_fname.cmn'! grid_fnam_common

!     This routine accesses the standard_latitude variables from the
!     .parms file via the common block. Note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
      endif

      std_lat1 = standard_latitude
      std_lat2 = standard_latitude2

      istatus = 1
      return
      end

      subroutine get_vertical_grid(vert_grid,istatus)

      include 'lapsparms.cmn' ! standard_latitude, standard_latitude2
      include 'grid_fname.cmn'! grid_fnam_common

      character*40 vert_grid

!     This routine accesses the vertical_grid variables from the
!     .parms file via the common block. Note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'
          return
      endif

      vert_grid = vertical_grid

      istatus = 1
      return
      end


      subroutine get_maxstns(maxstns,istatus)

      include 'lapsparms.cmn' ! maxstns_cmn
      include 'grid_fname.cmn'! grid_fnam_common

!     This routine accesses the maxstns_cmn variable from the
!     .parms file via the common block. Note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
      endif

      maxstns = maxstns_cmn

      istatus = 1
      return
      end


      subroutine get_c8_project(c8_project,istatus)

      include 'lapsparms.cmn' ! c8_project
      include 'grid_fname.cmn'! grid_fnam_common

      character*8 c8_project

!     This routine accesses the c8_project_common variable from the
!     .parms file via the common block. Note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' WARNING, get_laps_config not successfully called'       
          return
      endif

      c8_project = c8_project_common

      istatus = 1
      return
      end

c------------------------------------------------------------
      subroutine get_raddat_type(c_raddat_type_ret,istatus)
      include 'lapsparms.cmn' ! c_raddat_type
      include 'grid_fname.cmn'! grid_fnam_common

      character*3 c_raddat_type_ret

!     This routine accesses the c8_project_common variable from the
!     .parms file via the common block. Note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' WARNING, get_laps_config not successfully called'
          return
      endif

      c_raddat_type_ret = c_raddat_type

      istatus = 1
      return
      end
c----------------------------------------------------------

      subroutine get_c6_maproj(c6_maproj_ret,istatus)

      include 'lapsparms.cmn' ! c6_maproj
      include 'grid_fname.cmn'! grid_fnam_common

      character*6 c6_maproj_ret

!     This routine accesses the c6_maproj variable from the
!     .parms file via the common block. Note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
      endif

      c6_maproj_ret = c6_maproj

      istatus = 1
      return
      end


      subroutine get_c80_description(c80_description_ret,istatus)

      include 'lapsparms.cmn' ! c80_description
      include 'grid_fname.cmn'! grid_fnam_common

      character*80 c80_description_ret
!     This routine accesses the c80_description variable from the
!     .parms file via the common block. Note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
      endif

      c80_description_ret = c80_description

      istatus = 1
      return
      end


      subroutine get_laps_dimensions(nk,istatus)

      include 'lapsparms.cmn'              ! nk_laps
      include 'grid_fname.cmn'             ! grid_fnam_common

!     This routine accesses the nk variable from the
!     .parms file via the common block. Note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
      endif

      nk = nk_laps

      istatus = 1
      return
      end


      subroutine get_laps_cycle_time(laps_cycle_time,istatus)

      include 'lapsparms.cmn'              ! laps_cycle_time_cmn
      include 'grid_fname.cmn'             ! grid_fnam_common

!     This routine accesses the laps_cycle_time variable from the
!     .parms file via the common block. Note the variable name in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
      endif

      laps_cycle_time = laps_cycle_time_cmn

      istatus = 1
      return
      end


      subroutine get_grid_dim_xy(NX_L,NY_L,istatus)

      include 'lapsparms.cmn'              ! NX_L_CMN, NY_L_CMN
      include 'grid_fname.cmn'             ! grid_fnam_common

!     character*80 grid_fnam

!     This routine accesses the NX_L and NY_L variables from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

!     grid_fnam = grid_fnam_common
      print*,'Here: get_grid_dim_xy'
      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
      endif

      NX_L = NX_L_CMN
      NY_L = NY_L_CMN

      istatus = 1
      return
      end

      subroutine get_topo_parms(silavwt_parm,toptwvl_parm,istatus)

      include 'lapsparms.cmn' ! silavwt_cmn, toptwvl_cmn
      include 'grid_fname.cmn'! grid_fnam_common

!     This routine accesses the silavwt_parm and toptwvl_parm
!     variables from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
      endif

      silavwt_parm = silavwt_parm_cmn
      toptwvl_parm = toptwvl_parm_cmn

      istatus = 1
      return
      end

      subroutine get_meso_sao_pirep(n_meso,n_sao,n_pirep,istatus)

      include 'lapsparms.cmn' ! maxstns_cmn, n_pirep_cmn
      include 'grid_fname.cmn'! grid_fnam_common

!     This routine accesses the maxstns and n_pirep
!     variables from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
      endif

      n_meso  = maxstns_cmn
      n_sao   = maxstns_cmn
      n_pirep = n_pirep_cmn

      istatus = 1
      return
      end

      subroutine get_max_radars (max_radars, istatus)

      include 'lapsparms.cmn' ! max_radars_cmn
      include 'grid_fname.cmn'! grid_fnam_common

!     This routine accesses the max_radars variable from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
      endif

      max_radars = max_radars_cmn

      istatus = 1
      return
      end

      subroutine get_max_stations (maxstns, istatus)

      include 'lapsparms.cmn' ! maxstns_cmn
      include 'grid_fname.cmn'! grid_fnam_common

!     This routine accesses the maxstns variable from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
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
      include 'grid_fname.cmn'! grid_fnam_common

      integer*4 vert_rad_pirep
      integer*4 vert_rad_sao
      integer*4 vert_rad_meso
      integer*4 vert_rad_prof

!     This routine accesses the vert_rad_pirep, etc., variables from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
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
      include 'grid_fname.cmn'! grid_fnam_common

!     This routine accesses the r_missing_data variable from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
      endif

      r_missing_data = r_missing_data_cmn

      istatus = 1
      return
      end

      subroutine get_i2_missing_data(i2_missing_data, istatus)

      include 'lapsparms.cmn' ! i2_missing_data_cmn
      include 'grid_fname.cmn'! grid_fnam_common

!     This routine accesses the i2_missing_data variable from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
      endif

      i2_missing_data = i2_missing_data_cmn

      istatus = 1
      return
      end


      subroutine get_i_perimeter(i_perimeter, istatus)

      include 'lapsparms.cmn' ! i_perimeter_cmn
      include 'grid_fname.cmn'! grid_fnam_common

!     This routine accesses the 'i_perimeter' variable from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
      endif

      i_perimeter = i_perimeter_cmn

      istatus = 1
      return
      end


      subroutine get_ref_base(ref_base, istatus)

      include 'lapsparms.cmn' ! ref_base_cmn
      include 'grid_fname.cmn'! grid_fnam_common

!     This routine accesses the 'ref_base' variable from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
      endif

      ref_base = ref_base_cmn

      istatus = 1
      return
      end

      subroutine get_ref_base_useable(ref_base_useable, istatus)

      include 'lapsparms.cmn' ! ref_base_useable
      include 'grid_fname.cmn'! grid_fnam_common

!     This routine accesses the 'ref_base_useable' variable from the
!     .parms file via the common block. Note the variable names in the
!     argument list may be different in the calling routine

      call get_laps_config(grid_fnam_common,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
      endif

      ref_base_useable = ref_base_useable_cmn

      istatus = 1
      return
      end
c
c-----------------------------------------------------------------------
c
      function ltest_vertical_grid(c_vertical_grid)

!     The input c_vertical_grid is a string you are testing
!     The output (.true. OR .false.) states whether the input vertical_grid 
!     name equals the vertical_grid in the Namelist

      include 'lapsparms.cmn' ! vertical_grid
      include 'grid_fname.cmn'! grid_fnam_common

      logical ltest_vertical_grid
      character*(*) c_vertical_grid
      character*40  cvgrid_test
      character*40  cvgrid_laps

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
          write(6,*)' ltest_vertical_grid: Error detected in '
     1             ,'calling get_laps_config'
          ltest_vertical_grid = .false.
          return
      endif

      call downcase(c_vertical_grid,cvgrid_test)
      call downcase(vertical_grid,cvgrid_laps)

      call s_len(cvgrid_laps,len_cmn)
      call s_len(cvgrid_test,len_in)


      if(cvgrid_laps(1:len_cmn).eq.cvgrid_test(1:len_in))then
          ltest_vertical_grid = .true.
      else
          ltest_vertical_grid = .false.
      endif

      return
      end
c


      subroutine get_pressure_interval(pressure_interval,istatus)

      include 'lapsparms.cmn' ! PRESSURE_INTERVAL_L
      include 'grid_fname.cmn'! grid_fnam_common

!     This routine accesses the pressure_interval variable from the
!     .parms file via the common block. Note the variable name in the
!     argument list is different in the calling routine

!     The result 'pressure_interval' has real meaning only when we are
!     using a uniform 'PRESSURE' grid.

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'       
          return
      endif

      pressure_interval = PRESSURE_INTERVAL_L

      return
      end

      subroutine get_pressure_bottom(pressure_bottom,istatus)

      include 'lapsparms.cmn' ! PRESSURE_BOTTOM_L
      include 'grid_fname.cmn'! grid_fnam_common

!     This routine accesses the pressure_bottom variable from the
!     .parms file via the common block. Note the variable name in the
!     argument list is different in the calling routine

!     The result 'pressure_bottom' has real meaning only when we are
!     using a uniform 'PRESSURE' grid.

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'
          return
      endif

      pressure_bottom = PRESSURE_BOTTOM_L

      return
      end


      subroutine get_earth_radius(erad,istatus)

      erad = 6367000.

      istatus = 1

      return
      end

      subroutine get_fdda_model_source(fdda_model_source
     1,n_fdda_models,istatus)

      include 'lapsparms.cmn' ! FDDA_MODEL_SOURCE_CMN
      include 'grid_fname.cmn'! grid_fnam_common

      character*(*) fdda_model_source(maxbgmodels)
      integer n_fdda_models,istatus

!     This routine accesses the fdda_model_source variable from the
!     .parms file via the common block. Note the variable name in the
!     argument list is different in the calling routine

      do i=1,maxbgmodels
         fdda_model_source_cmn(i) = ' '
      enddo
      iflag_lapsparms_cmn=0
      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'
          return
      endif

      n_fdda_models = 0
      do i=1,maxbgmodels
         if(fdda_model_source_cmn(i).ne. ' ')then
            n_fdda_models=n_fdda_models+1
            fdda_model_source(n_fdda_models)=fdda_model_source_cmn(i)
         endif
      enddo

      return
      end

      subroutine get_path_to_topo_10m(path_to_topo_10m,istatus)

      include 'lapsparms.cmn' ! path_to_topt_10m
      include 'grid_fname.cmn'

!     This routine accesses the path_to_topt_10m  variable from the
!     .parms file via the common block.

      character*200 path_to_topo_10m

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'
          return
      endif

      path_to_topo_10m =  path_to_topt10m

      return
      end

      subroutine get_path_to_topo_30s(path_to_topo_30s,istatus)

      include 'lapsparms.cmn' ! path_to_topt_30s
      include 'grid_fname.cmn'! grid_fnam_common

!     This routine accesses the path_to_topt_10m  variable from the
!     .parms file via the common block. 

      character*200 path_to_topo_30s

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'
          return
      endif

      path_to_topo_30s =  path_to_topt30s

      return
      end

      subroutine get_path_to_pctl_10m(path_to_pctl_10m,istatus)

      include 'lapsparms.cmn' ! path_to_pctl_10m
      include 'grid_fname.cmn'! grid_fnam_common

!     This routine accesses the path_to_pctl_10m  variable from the
!     .parms file via the common block. 

      character*200 path_to_pctl_10m

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'
          return
      endif

      path_to_pctl_10m =  path_to_pctl10m

      return
      end

      subroutine get_path_to_soil_2m(path_to_soil_2m,istatus)

      include 'lapsparms.cmn' ! path_to_soil_2m
      include 'grid_fname.cmn'! grid_fnam_common

!     This routine accesses the path_to_soil_2m  variable from the
!     .parms file via the common block.

      character*200 path_to_soil_2m

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
          write(6,*)' ERROR, get_laps_config not successfully called'
          return
      endif

      path_to_soil_2m =  path_to_soil2m

      return
      end


       subroutine array_minmax(a,ni,nj,rmin,rmax,r_missing_data)

       real*4 a(ni,nj)

       rmin =  abs(r_missing_data)
       rmax = -abs(r_missing_data)

       do i = 1,ni
       do j = 1,nj
           rmin = min(rmin,a(i,j))
           rmax = max(rmax,a(i,j))
       enddo ! j
       enddo ! i

       return
       end
