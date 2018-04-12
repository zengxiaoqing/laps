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

       subroutine laps_anl(uobs,vobs
     1     ,obs_point,max_obs,nobs_point                               ! I
     1     ,n_radars,istat_radar_vel                                   ! I
     1     ,vr_obs_unfltrd,vr_nyq,v_nyquist_in,idx_radar_a             ! I
     1     ,n_var                                                      ! I
     1     ,uanl,vanl                                                  ! O
     1     ,wt_p,weight_bkg_const,rms_thresh_wind                      ! I
     1     ,r0_barnes_max_m,brns_conv_rate_wind                        ! I
     1     ,max_radars                                                 ! I
     1     ,n_radarobs_tot_unfltrd,rlat_radar,rlon_radar,rheight_radar ! I
     1     ,thresh_2_radarobs_lvl_unfltrd                              ! I
     1     ,thresh_4_radarobs_lvl_unfltrd                              ! I
     1     ,thresh_9_radarobs_lvl_unfltrd                              ! I
     1     ,thresh_25_radarobs_lvl_unfltrd                             ! I
     1     ,u_laps_bkg,v_laps_bkg                                      ! I/L
     1     ,imax,jmax,kmax,lat,lon                                     ! I
     1     ,nx_r,ny_r,ioffset,joffset                                  ! I
     1     ,i4time,grid_spacing_m                                      ! I
     1     ,r_missing_data                                             ! I
     1     ,heights_3d                                                 ! I
     1     ,i4_loop_total                                              ! O
     1     ,l_derived_output,l_grid_north,l_3pass,l_correct_unfolding  ! I
     1     ,weight_cdw,weight_sfc,weight_pirep,weight_prof,weight_radar! I
     1     ,istatus)

!     This routine uses the inputted wind data and actually does the analysis

!     1992          Steve Albers
!     1992 Apr      Steve Albers      Subroutine to make u/v radar obs
!     1992 Apr      Steve Albers      Looping for multiple radars
!     1992 Aug      LincLb/S. Albers  More efficient, combining u+v analyses
!     1992 Aug      LincLb/S. Albers  Barnes_multivariate changed to be more
!                                     parallelizable
!     1992 Oct      S. Albers         Add intermediate pass with multi-Doppler
!                                     obs only
!     1994 May      S. Albers         Add Local wt_p_spread array so that
!                                     wt_p array is not modified by analysis
!     1994 Nov 01   S. Albers         Remove lapsprms.for include, add arguments
!     1995 Aug      S. Albers         Some U/V arrays now equivalenced
!     1995 Sep      K. Brewster       Added nyquist as input grid array

!********************ARGUMENT LIST********************************************

      integer max_obs
!     parameter (max_obs = 40000)       
      include 'barnesob.inc'
      type (barnesob) :: obs_point(max_obs)   ! Full Wind Obs  - Non-radar data
      type (barnesob) :: obs_point_qced(max_obs) ! QC'd Obs    - Non-radar data
      type (barnesob) :: obs_radar(max_obs)   ! Full Wind Obs  - Radar data
      type (barnesob) :: obs_radar_multi(max_obs) ! Full Wind Obs  - Multi-Radar data
      type (barnesob) :: obs_barnes(max_obs)  ! QC'd obs       - All Data

      integer n_var                                                ! Input
      integer imax,jmax,kmax        ! 3D array dimensions        ! Input

!     3D arrays of u/v observations, all data sources, one datum per gridpoint.
      real uobs(imax,jmax,kmax),vobs(imax,jmax,kmax)             ! Input

      integer n_radars   ! Actual number of radars having data     Input

!     Final pass analyzed winds
      real uanl(imax,jmax,kmax),vanl(imax,jmax,kmax)             ! Output

      real, allocatable, dimension(:,:,:,:) :: varbuff           ! Local
!     real varbuff(imax,jmax,kmax,n_var)                         ! Local

!     3D array of observation weights, depends on data type
!     The choices are outlined below
      real    wt_p(imax,jmax,kmax)                               ! Input
!     real    wt_p_radar(imax,jmax,kmax)                         ! Local
      real, allocatable, dimension(:,:,:) :: wt_p_radar          ! Local

!     Model background field
      real u_laps_bkg(imax,jmax,kmax),v_laps_bkg(imax,jmax,kmax) ! Input

!     Arrays of lat and lon for each gridpoint
      real lat(imax,jmax),lon(imax,jmax)                         ! Input

      integer i4time                                             ! Input

      real grid_spacing_m                                        ! Input
      real r_missing_data                  ! missing data value    Input

!-----Radar Data -----------------------------------------------------------

      integer max_radars            ! max possible # of radars         Input

!     4D Radial velocity array (all radars)
      real vr_obs_unfltrd(nx_r,ny_r,kmax,max_radars)                 ! Input
      real vr_nyq(nx_r,ny_r,kmax,max_radars)                         ! Input
      integer idx_radar_a(max_radars)                                ! Input

!     Nyquist velocity (if known and constant) for each radar
      real v_nyquist_in(max_radars)                                  ! Input

!     Location of each radar
      real rlat_radar(max_radars),rlon_radar(max_radars)             ! Input
     1                     ,rheight_radar(max_radars)

      integer ioffset(max_radars),joffset(max_radars)

      integer thresh_2_radarobs_lvl_unfltrd                          ! Input
     1         ,thresh_4_radarobs_lvl_unfltrd
     1         ,thresh_9_radarobs_lvl_unfltrd
     1         ,thresh_25_radarobs_lvl_unfltrd

!     # of radar obs before filtering for each radar (modified by QC)
      integer n_radarobs_tot_unfltrd(max_radars)                     ! Input/Modified
      real   heights_3d(imax,jmax,kmax)                              ! Input

      real vr_obs_unfltrd_3d(imax,jmax,kmax)                         ! Local
      real vr_nyq_3d(imax,jmax,kmax)                                 ! Local

!--------------------------------------------------------------------------------

!     First pass analyzed winds (innovation analysis with non-radar data)
      real, allocatable, dimension(:,:,:) :: upass1                  ! Local
      real, allocatable, dimension(:,:,:) :: vpass1                  ! Local

      real, allocatable, dimension(:,:,:) :: pres_3d                 ! Local

!     Note that aerr is only a dummy ATTM 
!     real aerr(imax,jmax,kmax,n_var)                                ! Local
      real, allocatable, dimension(:,:,:,:) :: aerr                  ! Local
      real, allocatable, dimension(:,:,:,:) :: varobs_diff_spread    ! Local

      integer n_obs_lvl(kmax)                                        ! Local
      logical  l_analyze(kmax) ! This depends on presence of radar obs ! Local
      logical  l_derived_output ! Flag for producing derived output    ! Input
      logical  l_grid_north     ! Flag for grid north or true north    ! Input
      logical  l_3pass          ! Flag for doing 3 pass analysis       ! Input
      logical  l_correct_unfolding ! Flag for dealiasing               ! Input
      logical  l_point_struct, l_variational, l_withheld_only

      real   rms_thresh                                              ! Input
      real   weight_bkg_const                                        ! Input

!     These are the weights of the various data types (filling the 3D array)
      real weight_cdw,weight_sfc,weight_pirep,weight_prof
     1      ,weight_radar                                              ! Input

      integer istatus         ! (1 is good)                          ! Output

      integer  n_fnorm_dum

!****************END DECLARATIONS *********************************************

      write(6,*)' Subroutine laps_anl...'

      l_point_struct = .true. 
      l_variational = .false.
      l_withheld_only = .false.

      single_ratio = 0.4

      ialloc_varobs_diff_spread = 0

csms$serial(default=ignore)  begin              

!     Compare background to obs
      call compare_wind(
     1            u_laps_bkg,v_laps_bkg,' FG ',
     1            istat_radar_vel,max_radars,vr_obs_unfltrd,n_radars,
     1            nx_r,ny_r,ioffset,joffset,
     1            rlat_radar,rlon_radar,rheight_radar,
     1            lat,lon,
     1            imax,jmax,kmax,r_missing_data,
     1            obs_point,max_obs,nobs_point,l_point_struct,
     1            l_withheld_only,
     1            weight_pirep,weight_prof,weight_sfc,weight_cdw,
     1            uobs,vobs,wt_p,istatus)
      rms_thresh_norm = rms_thresh_wind          

!     Subtract the background from the non-radar obs, then apply QC thresholds
!     and spread the obs vertically.

      write(6,91)
91    format(1x,' Subtracting the background from the obs'
     1         ,' then spreading the obs vertically.'
     1         /'       i    j    k  kk   udf   vdf     '
     1         ,'uob   vob     ubg   vbg vcdf  wt')

      allocate( pres_3d(imax,jmax,kmax), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate pres_3d'
          stop
      endif

      call get_pres_3d(i4time,imax,jmax,kmax,pres_3d,istatus)
      if(istatus .ne. 1)return

      call get_rep_pres_intvl(pres_3d,imax,jmax,kmax,rep_pres_intvl
     1                       ,istatus)

      if(allocated(pres_3d))deallocate(pres_3d)

!     QC the obs and place in different data structure
      call calc_qced_obs( 
     &  max_obs,nobs_point,u_laps_bkg,v_laps_bkg,
     &  imax,jmax,kmax,
     &  obs_point(:)%i,obs_point(:)%j,obs_point(:)%k,
     &  obs_point(:)%ri,obs_point(:)%rj,obs_point(:)%rk,
     &  obs_point(:)%valuef(1),obs_point(:)%valuef(2),
     &  obs_point(:)%value(1), obs_point(:)%value(2),
     &  obs_point(:)%weight,
     &  obs_point(:)%vert_rad_rat,
     &  obs_point(:)%elev,
     &  obs_point(:)%ldf,
     &  obs_point(:)%mask_sea,
     &  obs_point(:)%i4time,
     &  obs_point(:)%l_withhold,
     &  obs_point(:)%type,
     &  obs_point(:)%file,
     &  obs_point_qced(:)%i,obs_point_qced(:)%j,obs_point_qced(:)%k,
     &  obs_point_qced(:)%ri,obs_point_qced(:)%rj,obs_point_qced(:)%rk,
     &  obs_point_qced(:)%valuef(1),obs_point_qced(:)%valuef(2),
     &  obs_point_qced(:)%value(1),obs_point_qced(:)%value(2),
     &  obs_point_qced(:)%weight,
     &  obs_point_qced(:)%vert_rad_rat,
     &  obs_point_qced(:)%elev,
     &  obs_point_qced(:)%ldf,
     &  obs_point_qced(:)%mask_sea,
     &  obs_point_qced(:)%i4time,
     &  obs_point_qced(:)%type,
     &  obs_point_qced(:)%file, 
     &  n_qc_total_good)

!     Debugging only
      write(6,*)'obs_point(1)'
      write(6,*)obs_point(1)
      write(6,*)'obs_point_qced(1)'
      write(6,*)obs_point_qced(1)
csms$serial end

      ncnt_total = n_qc_total_good

      if(l_variational)then ! call variational routine
          return
      endif

      call get_inst_err2(r_missing_data                               ! I
     1                  ,obs_point_qced,max_obs,n_qc_total_good       ! I
     1                  ,rms_thresh_norm                              ! I
     1                  ,i4time                                       ! I
     1                  ,rms_inst,rms_thresh)                         ! O

      allocate( aerr(imax,jmax,kmax,n_var), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate aerr'
          stop
      endif

      aerr = weight_bkg_const ! set array to constant weight for now

      allocate( varbuff(imax,jmax,kmax,n_var), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate varbuff'
          stop
      endif

!     No Radar Data
      call barnes_multivariate(varbuff                                ! O
     1        ,n_var,ncnt_total,obs_point_qced                        ! I
     1        ,imax,jmax,kmax,grid_spacing_m,rep_pres_intvl           ! I
     1        ,aerr,i4_loop_total                                     ! I/O 
     1        ,wt_p,fnorm_dum,n_fnorm_dum                             ! I
     1        ,l_analyze_dum,.false.,rms_thresh                       ! I
     1        ,r0_barnes_max_m,brns_conv_rate_wind                    ! I
     1        ,topo_dum,rland_frac_dum,1,1                            ! I
     1        ,n_obs_lvl,istatus)                                     ! O
      if(istatus .ne. 1)return

      do ivar=1,n_var
          call qc_field_3d('U3',varbuff(1,1,1,ivar),imax,jmax,kmax
     1                                             ,istatus)       
          if(istatus .ne. 1)then
              write(6,*)
     1        ' ERROR: QC flag in barnes_multivariate output ',ivar       
              return
          endif
      enddo ! ivar

      allocate( wt_p_radar(imax,jmax,kmax), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate wt_p_radar'
          stop
      endif

      wt_p_radar = wt_p

csms$print_mode(async) begin
csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 1 processor=',me

csms$serial(default=ignore)  begin              

      write(6,*)' Allocating upass1,vpass1'

      allocate( upass1(imax,jmax,kmax), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate upass1'
     1             ,istat_alloc,imax,jmax,kmax
          stop
      endif
!     call maxminavIJK(upass1,imax,jmax,kmax)

      allocate( vpass1(imax,jmax,kmax), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate vpass1'
     1             ,istat_alloc,imax,jmax,kmax
          stop
      endif
!     call maxminavIJK(vpass1,imax,jmax,kmax)

      call move_3d(varbuff(1,1,1,1),upass1,imax,jmax,kmax)
      call move_3d(varbuff(1,1,1,2),vpass1,imax,jmax,kmax)

      I4_elapsed = ishow_timer()

!     Perform radar QC by differencing radial velocities and first pass analysis
      do i_radar = 1,n_radars

          v_nyquist_global = v_nyquist_in(i_radar)
          
          write(6,*)' Radar QC for radar #/v_nyq/lat/lon ',i_radar
     1             ,v_nyquist_global
     1             ,rlat_radar(i_radar),rlon_radar(i_radar)
          call qc_radar_obs(
     1           imax,jmax,kmax                             ! Input
     1          ,r_missing_data                             ! Input
     1          ,nx_r,ny_r,ioffset(i_radar),joffset(i_radar)! Input
     1          ,vr_obs_unfltrd(1,1,1,i_radar)              ! Input/Output
     1          ,vr_nyq(1,1,1,i_radar)                      ! Input
     1          ,n_radarobs_tot_unfltrd(i_radar)            ! Input/Output
     1          ,lat,lon                                    ! Input
     1          ,rlat_radar(i_radar),rlon_radar(i_radar)    ! Input
     1          ,rheight_radar(i_radar)                     ! Input
     1          ,upass1,vpass1  ! 1st pass anal             ! Input
     1          ,u_laps_bkg,v_laps_bkg                      ! Input
     1          ,v_nyquist_global                           ! Input
     1          ,l_correct_unfolding,l_grid_north           ! Input
     1          ,istatus                                    ! Input/Output
     1                                                          )
      enddo ! i_radar

      I4_elapsed = ishow_timer()

!     Perform analysis with radar data added in
      do k = 1,kmax
          l_analyze(k) = .false.
      enddo ! k

csms$serial end

!     Fill 'obs_barnes' at this point in case there is no radar data
      obs_barnes = obs_point_qced
      ncnt_total = n_qc_total_good

      if(n_radars .le. 1 .or. .not. l_3pass)then ! Single Doppler (or no radar) Option

csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 2 processor=',me

csms$serial(<wt_p_radar , 
csms$>       icount_radar_total, out>:default=ignore)  begin

          mode = 1 ! All radar obs (in this case single Doppler)

!         Take the data from all the radars and add the derived radar obs into
!         uobs_diff_spread and vobs_diff_spread (varobs_diff_spread)
          if(l_point_struct)then
              wt_p_radar = r_missing_data                 ! Initialize

              if(ialloc_varobs_diff_spread .eq. 0)then
                  allocate(varobs_diff_spread(imax,jmax,kmax,n_var)
     1                    ,STAT=istat_alloc)      
                  if(istat_alloc .ne. 0)then
                      write(6,*)
     1                   ' ERROR: Could not allocate varobs_diff_spread'      
                      stop
                  else
                      ialloc_varobs_diff_spread = 1
                  endif
              endif

              varobs_diff_spread = r_missing_data         ! Initialize
          endif

          call insert_derived_radar_obs(
     1         mode                                       ! Input
     1        ,n_radars,max_radars,idx_radar_a            ! Input
     1        ,imax,jmax,kmax                             ! Input
     1        ,nx_r,ny_r,ioffset,joffset                  ! Input
     1        ,r_missing_data                             ! Input
     1        ,heights_3d                                 ! Input
     1        ,vr_obs_unfltrd                             ! Input
     1        ,i4time                                     ! Input
     1        ,lat,lon                                    ! Input
     1        ,rlat_radar,rlon_radar                      ! Input
     1        ,rheight_radar                              ! Input
     1        ,upass1,vpass1                              ! Input
     1        ,u_laps_bkg,v_laps_bkg                      ! Input
!    1        ,weight_radar                               ! Input
     1        ,l_derived_output,l_grid_north              ! Input
     1        ,wt_p_radar                                 ! Input/Output
     1        ,varobs_diff_spread(1,1,1,1),varobs_diff_spread(1,1,1,2) ! I/O
     1        ,l_analyze,icount_radar_total               ! Output
     1        ,n_radarobs_tot_unfltrd                     ! Input
     1        ,istatus                                    ! Input/Output
     1                                                          )

csms$serial end

          if(icount_radar_total .gt. 0)then 

csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 3 processor=',me

csms$serial(<rms_thresh, out>:default=ignore)  begin              
              I4_elapsed = ishow_timer()

              write(6,*)' Calling barnes with single radar obs added'       

csms$serial end

              if(l_point_struct)then
                  call arrays_to_barnesobs(imax,jmax,kmax             ! I
     1                              ,r_missing_data                   ! I
     1                              ,varobs_diff_spread,wt_p_radar    ! I
     1                              ,n_var,max_obs                    ! I
     1                              ,obs_radar                        ! O
     1                              ,ncnt_radar,weight_radar_total    ! O
     1                              ,istatus)                         ! O
                  if(istatus .ne. 1)return

                  if(allocated(varobs_diff_spread))
     1              deallocate(varobs_diff_spread)
                  ialloc_varobs_diff_spread = 0

!                 Combine radar (obs_radar) and non-radar (obs_point_qced) 
!                 data structures into new structure (obs_barnes)
                  obs_barnes = obs_point_qced
                  ncnt_total = n_qc_total_good
                  do i = 1,ncnt_radar
                      ncnt_total = ncnt_total + 1
                      obs_radar(i)%i4time = i4time
                      obs_barnes(ncnt_total) = obs_radar(i)
                      obs_barnes(ncnt_total)%type = 'radar'
                  enddo ! i
              endif

              call get_inst_err2(r_missing_data                       ! I
     1                  ,obs_barnes,max_obs,ncnt_total                ! I
     1                  ,rms_thresh_norm                              ! I
     1                  ,i4time                                       ! I
     1                  ,rms_inst,rms_thresh)                         ! O

!             Single Doppler Added
              call barnes_multivariate(varbuff                           ! O
     1           ,n_var,ncnt_total,obs_barnes                            ! I
     1           ,imax,jmax,kmax,grid_spacing_m,rep_pres_intvl           ! I
     1           ,aerr,i4_loop_total                                     ! I/O 
     1           ,wt_p_radar,fnorm_dum,n_fnorm_dum                       ! I
     1           ,l_analyze_dum,.false.,rms_thresh                       ! I
     1           ,r0_barnes_max_m*single_ratio,brns_conv_rate_wind       ! I
     1           ,topo_dum,rland_frac_dum,1,1                            ! I
     1           ,n_obs_lvl,istatus)                                     ! O

csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 4 processor=',me

csms$serial(default=ignore)  begin              

              call move_3d(varbuff(1,1,1,1),uanl,imax,jmax,kmax)
              call move_3d(varbuff(1,1,1,2),vanl,imax,jmax,kmax)

              if(istatus .ne. 1)return

              I4_elapsed = ishow_timer()

csms$serial end

          endif ! There is any radar data

      else ! n_radars .gt. 1

csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 5 processor=',me

csms$serial(<wt_p_radar, icount_radar_total, 
csms$>                    out>:default=ignore) begin       

          mode = 2 ! Only multi-Doppler obs

!         Take the data from all the radars and add the derived radar obs into
!         uobs_diff_spread and vobs_diff_spread (varobs_diff_spread)
          if(l_point_struct)then
              wt_p_radar = r_missing_data                 ! Initialize

              if(ialloc_varobs_diff_spread .eq. 0)then
                  allocate(varobs_diff_spread(imax,jmax,kmax,n_var)
     1                    ,STAT=istat_alloc)      
                  if(istat_alloc .ne. 0)then
                      write(6,*)
     1                   ' ERROR: Could not allocate varobs_diff_spread'   
                      stop
                  else
                      ialloc_varobs_diff_spread = 1
                  endif
              endif

              varobs_diff_spread = r_missing_data         ! Initialize
          endif

          call insert_derived_radar_obs(
     1         mode                                       ! Input
     1        ,n_radars,max_radars,idx_radar_a            ! Input
     1        ,imax,jmax,kmax                             ! Input
     1        ,nx_r,ny_r,ioffset,joffset                  ! Input
     1        ,r_missing_data                             ! Input
     1        ,heights_3d                                 ! Input
     1        ,vr_obs_unfltrd                             ! Input
!    1        ,thresh_2_radarobs_lvl_unfltrd              ! Input
!    1        ,thresh_4_radarobs_lvl_unfltrd              ! Input
!    1        ,thresh_9_radarobs_lvl_unfltrd              ! Input
!    1        ,thresh_25_radarobs_lvl_unfltrd             ! Input
     1        ,i4time                                     ! Input
     1        ,lat,lon                                    ! Input
     1        ,rlat_radar,rlon_radar                      ! Input
     1        ,rheight_radar                              ! Input
     1        ,upass1,vpass1                              ! Input
     1        ,u_laps_bkg,v_laps_bkg                      ! Input
!    1        ,weight_radar                               ! Input
     1        ,l_derived_output,l_grid_north              ! Input
     1        ,wt_p_radar                                 ! Input/Output
     1        ,varobs_diff_spread(1,1,1,1),varobs_diff_spread(1,1,1,2) ! I/O
     1        ,l_analyze,icount_radar_total               ! Output
     1        ,n_radarobs_tot_unfltrd                     ! Input
     1        ,istatus                                    ! Input/Output
     1                                                          )

csms$serial end

          if(icount_radar_total .gt. 0)then 

csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 6 processor=',me

csms$serial(<rms_thresh, out>:default=ignore)  begin              

              I4_elapsed = ishow_timer()

              write(6,401)icount_radar_total
 401          format(1x,' Analyzing with ',i5
     1                 ,' multi-doppler grid points')     

csms$serial end

              if(l_point_struct)then
                  call arrays_to_barnesobs(imax,jmax,kmax             ! I
     1                              ,r_missing_data                   ! I
     1                              ,varobs_diff_spread,wt_p_radar    ! I
     1                              ,n_var,max_obs                    ! I
     1                              ,obs_radar_multi                  ! O
     1                              ,ncnt_radar,weight_radar_total    ! O
     1                              ,istatus)                         ! O
                  if(istatus .ne. 1)return

                  if(allocated(varobs_diff_spread))
     1              deallocate(varobs_diff_spread)
                  ialloc_varobs_diff_spread = 0

!                 Combine radar (obs_radar_multi) and non-radar (obs_point_qced) 
!                 data structures into new structure (obs_barnes)
                  obs_barnes = obs_point_qced
                  ncnt_total = n_qc_total_good
                  do i = 1,ncnt_radar
                      ncnt_total = ncnt_total + 1
                      obs_radar_multi(i)%i4time = i4time
                      obs_barnes(ncnt_total) = obs_radar_multi(i)
                      obs_barnes(ncnt_total)%type = 'radar'
                  enddo ! i

              endif

              call get_inst_err2(r_missing_data                         ! I
     1                  ,obs_barnes,max_obs,ncnt_total                  ! I
     1                  ,rms_thresh_norm                                ! I
     1                  ,i4time                                         ! I
     1                  ,rms_inst,rms_thresh)                           ! O

!             Multi-Doppler Added
              call barnes_multivariate(varbuff                          ! O
     1          ,n_var,ncnt_total                                       ! I
     1          ,obs_barnes,imax,jmax,kmax,grid_spacing_m,rep_pres_intvl! I   
     1          ,aerr,i4_loop_total                                     ! I/O 
     1          ,wt_p_radar,fnorm_dum,n_fnorm_dum                       ! I
     1          ,l_analyze_dum,.false.,rms_thresh                       ! I
     1          ,r0_barnes_max_m*single_ratio,brns_conv_rate_wind       ! I
     1          ,topo_dum,rland_frac_dum,1,1                            ! I
     1          ,n_obs_lvl,istatus)                                     ! O

csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 7 processor=',me

csms$serial(default=ignore)  begin              

              call move_3d(varbuff(1,1,1,1),uanl,imax,jmax,kmax)
              call move_3d(varbuff(1,1,1,2),vanl,imax,jmax,kmax)

              if(istatus .ne. 1)return

csms$serial end

          endif

csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 8 processor=',me

csms$serial(<wt_p_radar , rms_thresh, out>
csms$>                                     :default=ignore)  begin

!         Make sure each level of uanl and vanl is initialized in the event it
!         was not analyzed.
          do k = 1,kmax
              if(.not. l_analyze(k)
     1               .or. icount_radar_total .eq. 0)then
                   write(6,412)k
412                format(' No multi-radar obs at lvl',i3,
     1                          ' Insert 1st pass into uanl,vanl')
                   do j=1,jmax
                   do i=1,imax
                       uanl(i,j,k) = upass1(i,j,k)
                       vanl(i,j,k) = vpass1(i,j,k)

                   enddo ! i
                   enddo ! j
              endif

          enddo ! k

          I4_elapsed = ishow_timer()

          mode = 1 ! Single and multi-Doppler obs

!         Take the data from all the radars and add the derived radar obs into
!         uobs_diff_spread and vobs_diff_spread
          if(l_point_struct)then
              wt_p_radar = r_missing_data                 ! Initialize

              if(ialloc_varobs_diff_spread .eq. 0)then
                  allocate(varobs_diff_spread(imax,jmax,kmax,n_var)
     1                    ,STAT=istat_alloc)      
                  if(istat_alloc .ne. 0)then
                      write(6,*)
     1                   ' ERROR: Could not allocate varobs_diff_spread'      
                      stop
                  else
                      ialloc_varobs_diff_spread = 1
                  endif
              endif

              varobs_diff_spread = r_missing_data         ! Initialize
          endif

          call insert_derived_radar_obs(
     1         mode                                       ! Input
     1        ,n_radars,max_radars,idx_radar_a            ! Input
     1        ,imax,jmax,kmax                             ! Input
     1        ,nx_r,ny_r,ioffset,joffset                  ! Input
     1        ,r_missing_data                             ! Input
     1        ,heights_3d                                 ! Input
     1        ,vr_obs_unfltrd                             ! Input
!    1        ,thresh_2_radarobs_lvl_unfltrd              ! Input
!    1        ,thresh_4_radarobs_lvl_unfltrd              ! Input
!    1        ,thresh_9_radarobs_lvl_unfltrd              ! Input
!    1        ,thresh_25_radarobs_lvl_unfltrd             ! Input
     1        ,i4time                                     ! Input
     1        ,lat,lon                                    ! Input
     1        ,rlat_radar,rlon_radar                      ! Input
     1        ,rheight_radar                              ! Input
     1        ,uanl,vanl                                  ! Input
     1        ,u_laps_bkg,v_laps_bkg                      ! Input
!    1        ,weight_radar                               ! Input
     1        ,l_derived_output,l_grid_north              ! Input
     1        ,wt_p_radar                                 ! Input/Output
     1        ,varobs_diff_spread(1,1,1,1),varobs_diff_spread(1,1,1,2)  ! I/O
     1        ,l_analyze,icount_radar_total               ! Output
     1        ,n_radarobs_tot_unfltrd                     ! Input
     1        ,istatus                                    ! Input/Output
     1                                                          )

          I4_elapsed = ishow_timer()

          write(6,*)' Calling barnes with single+multi radar obs added'       

csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 9 processor=',me
csms$serial end
csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 10 processor=',me

          if(l_point_struct)then
              call arrays_to_barnesobs(imax,jmax,kmax             ! I
     1                              ,r_missing_data                   ! I
     1                              ,varobs_diff_spread,wt_p_radar    ! I
     1                              ,n_var,max_obs                    ! I
     1                              ,obs_radar                        ! O
     1                              ,ncnt_radar,weight_radar_total    ! O
     1                              ,istatus)                         ! O
              if(istatus .ne. 1)return

              if(allocated(varobs_diff_spread))
     1          deallocate(varobs_diff_spread)
              ialloc_varobs_diff_spread = 0

!             Combine radar (obs_radar) and non-radar (obs_point_qced) 
!             data structures into new structure (obs_barnes)
              obs_barnes = obs_point_qced
              ncnt_total = n_qc_total_good

              do i = 1,ncnt_radar
                  ncnt_total = ncnt_total + 1
                  obs_radar(i)%i4time = i4time
                  obs_barnes(ncnt_total) = obs_radar(i)
                  obs_barnes(ncnt_total)%type = 'radar'
              enddo ! i

          endif

          call get_inst_err2(r_missing_data                           ! I
     1                  ,obs_barnes,max_obs,ncnt_total                ! I
     1                  ,rms_thresh_norm                              ! I
     1                  ,i4time                                       ! I
     1                  ,rms_inst,rms_thresh)                         ! O

!         Single and Multi-Doppler radar added
          call barnes_multivariate(varbuff                            ! O
     1       ,n_var,ncnt_total,obs_barnes                             ! I
     1       ,imax,jmax,kmax                                          ! I
     1       ,grid_spacing_m,rep_pres_intvl                           ! I
     1       ,aerr,i4_loop_total                                      ! I/O 
     1       ,wt_p_radar,fnorm_dum,n_fnorm_dum                        ! I
     1       ,l_analyze_dum,.false.,rms_thresh                        ! I
     1       ,r0_barnes_max_m*single_ratio,brns_conv_rate_wind        ! I
     1       ,topo_dum,rland_frac_dum,1,1                             ! I
     1       ,n_obs_lvl,istatus)                                      ! O

csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 11 processor=',me

csms$serial(default=ignore)  begin              

          call move_3d(varbuff(1,1,1,1),uanl,imax,jmax,kmax)
          call move_3d(varbuff(1,1,1,2),vanl,imax,jmax,kmax)

          if(istatus .ne. 1)return

          I4_elapsed = ishow_timer()

csms$serial end

      endif ! n_radars

      if(allocated(aerr))deallocate(aerr)
      if(allocated(varbuff))deallocate(varbuff)
      if(allocated(wt_p_radar))deallocate(wt_p_radar)

csms$serial(default=ignore)  begin              

      write(6,*)' Adding analyzed differences to the background '
     1         ,'to reconstruct full analyses'

      do k=1,kmax ! Add back differences for first pass

          do j=1,jmax
          do i=1,imax
              if(upass1(i,j,k) .ne. r_missing_data)then
                  upass1(i,j,k) = upass1(i,j,k) + u_laps_bkg(i,j,k)
                  vpass1(i,j,k) = vpass1(i,j,k) + v_laps_bkg(i,j,k)
              else
                  write(6,*)
     1            ' ERROR: Missing data value(s) detected in first'
     1           ,' pass at lvl',k
                  istatus = 0
                  return
              endif
          enddo ! i
          enddo ! j


          if(l_analyze(k) .OR. (.true. .and. icount_radar_total .gt. 0) ! l_3d
     1                         )then ! This depends on the presence of radar obs
              write(6,511)k
511           format(' Use 2nd Pass for lvl',i3)

              do j=1,jmax
              do i=1,imax
                  uanl(i,j,k) = uanl(i,j,k) + u_laps_bkg(i,j,k)
                  vanl(i,j,k) = vanl(i,j,k) + v_laps_bkg(i,j,k)

              enddo ! i
              enddo ! j

          else
              write(6,512)k
512           format(' Use 1st Pass for lvl',i3)
              do j=1,jmax
              do i=1,imax
                  uanl(i,j,k) = upass1(i,j,k)
                  vanl(i,j,k) = vpass1(i,j,k)

              enddo ! i
              enddo ! j

          endif
      enddo ! k

      if(.true.)then
!         Compare 1st pass analysis to obs
          call compare_wind(
     1            upass1,vpass1,'PS1 ',
     1            istat_radar_vel,max_radars,vr_obs_unfltrd,n_radars,
     1            nx_r,ny_r,ioffset,joffset, 
     1            rlat_radar,rlon_radar,rheight_radar,
     1            lat,lon,
     1            imax,jmax,kmax,r_missing_data,
     1            obs_barnes,max_obs,ncnt_total,l_point_struct,
     1            l_withheld_only,
     1            weight_pirep,weight_prof,weight_sfc,weight_cdw,
     1            uobs,vobs,wt_p,istatus)

      endif ! Last iteration

      write(6,*)' Deallocate upass1, vpass1'
      deallocate(upass1)
      deallocate(vpass1)

!     Compare final analysis to QC'd obs
      call compare_wind(
     1            uanl,vanl,'LAPS',
     1            istat_radar_vel,max_radars,vr_obs_unfltrd,n_radars,
     1            nx_r,ny_r,ioffset,joffset,
     1            rlat_radar,rlon_radar,rheight_radar,
     1            lat,lon,
     1            imax,jmax,kmax,r_missing_data,
     1            obs_barnes,max_obs,ncnt_total,l_point_struct,
     1            l_withheld_only,
     1            weight_pirep,weight_prof,weight_sfc,weight_cdw,
     1            uobs,vobs,wt_p,istatus)

!     Compare final analysis to withheld (Non-QC'd) obs
      l_withheld_only = .true.
      call compare_wind(
     1            uanl,vanl,'LAPS',
     1            istat_radar_vel,max_radars,vr_obs_unfltrd,n_radars,
     1            nx_r,ny_r,ioffset,joffset,
     1            rlat_radar,rlon_radar,rheight_radar,
     1            lat,lon,
     1            imax,jmax,kmax,r_missing_data,
     1            obs_point,max_obs,nobs_point,l_point_struct,
     1            l_withheld_only,
     1            weight_pirep,weight_prof,weight_sfc,weight_cdw,
     1            uobs,vobs,wt_p,istatus)

      istatus = 1

      write(6,*)' End of subroutine laps_anl'

csms$serial end

      return
      end


      subroutine get_inst_err(imax,jmax,kmax,r_missing_data        ! I
     1                       ,wt_p,rms_thresh_norm                 ! I
     1                       ,rms_inst,rms_thresh)                 ! O

      real    wt_p(imax,jmax,kmax)                               ! Input

csms$ignore begin

      write(6,*)
      write(6,*)' subroutine get_inst_err...'

      n_obs_total = 0
      wt_p_inv_total = 0.

      do i = 1,imax
      do j = 1,jmax
      do k = 1,kmax
          if(wt_p(i,j,k) .ne. r_missing_data)then
              n_obs_total = n_obs_total + 1
              wt_p_inv_total = wt_p_inv_total + 1.0 / wt_p(i,j,k)
          endif
      enddo ! k
      enddo ! j
      enddo ! i

      if(n_obs_total .gt. 0)then
          wt_p_inv_ave = wt_p_inv_total / float(n_obs_total)
          rms_inst = sqrt(wt_p_inv_ave)
      else
          wt_p_inv_ave = 0.
          rms_inst = 0.
      endif

      rms_thresh = rms_inst * rms_thresh_norm

      write(6,*)' n_obs_total = ',n_obs_total
      write(6,*)' wt_p_inv_total,wt_p_inv_ave = '
     1           ,wt_p_inv_total,wt_p_inv_ave
      write(6,*)' rms_inst, rms_thresh = ',rms_inst,rms_thresh
csms$ignore end

      return
      end


      subroutine get_inst_err2(r_missing_data                       ! I
     1                        ,obs_barnes,max_obs,nobs_barnes       ! I
     1                        ,rms_thresh_norm                      ! I
     1                        ,i4time_sys                           ! I
     1                        ,rms_inst,rms_thresh)                 ! O

      integer max_obs
      include 'barnesob.inc'
      type (barnesob) :: obs_barnes(max_obs)                           

csms$ignore begin

      write(6,*)
      write(6,*)' subroutine get_inst_err2...'

      n_obs_total = 0
      wt_p_inv_total = 0.

      do i = 1,nobs_barnes
          n_obs_total = n_obs_total + 1
          call get_time_wt(i4time_sys,obs_barnes(i)%i4time,time_wt
     1                    ,istatus)
          wt_p_inv_total = wt_p_inv_total + 1.0 
     1                   / (obs_barnes(i)%weight * time_wt)
      enddo ! i

      if(n_obs_total .gt. 0)then
          wt_p_inv_ave = wt_p_inv_total / float(n_obs_total)
          rms_inst = sqrt(wt_p_inv_ave)
      else
          wt_p_inv_ave = 0.
          rms_inst = 0.
      endif

      rms_thresh = rms_inst * rms_thresh_norm

      write(6,*)' n_obs_total = ',n_obs_total
      write(6,*)' wt_p_inv_total,wt_p_inv_ave = '
     1           ,wt_p_inv_total,wt_p_inv_ave
      write(6,*)' rms_inst, rms_thresh = ',rms_inst,rms_thresh
csms$ignore end

      return
      end



      subroutine calc_qced_obs (
     &  max_obs,nobs_point,u_laps_bkg,v_laps_bkg,
     &  imax,jmax,kmax,
     &  obs_i,obs_j,obs_k,
     &  obs_ri,obs_rj,obs_rk,
     &  obs_fu,obs_fv,
     &  obs_u,obs_v,
     &  obs_weight,
     &  obs_vert_rad_rat,
     &  obs_elev,
     &  obs_ldf,
     &  obs_mask_sea,
     &  obs_i4time,
     &  obs_l_withhold,
     &  obs_type,
     &  obs_file,
     &  qced_i,qced_j,qced_k,
     &  qced_ri,qced_rj,qced_rk,
     &  qced_fu,qced_fv,
     &  qced_u,qced_v,
     &  qced_weight,
     &  qced_vert_rad_rat,
     &  qced_elev,
     &  qced_ldf,
     &  qced_mask_sea,
     &  qced_i4time,
     &  qced_type,
     &  qced_file, 
     &  n_qc_total_good)

       use mem_namelist, only : qc_thresh_wind_def, qc_thresh_wind_pin,
     &                          qc_thresh_wind_cdw, qc_thresh_wind_pro 

!      QC the obs and place in different arrays 

       real        , intent(in ) :: u_laps_bkg(imax,jmax,kmax)
       real        , intent(in ) :: v_laps_bkg(imax,jmax,kmax)
       integer     , intent(in ) :: obs_i(max_obs)
       integer     , intent(in ) :: obs_j(max_obs)
       integer     , intent(in ) :: obs_k(max_obs)
       real        , intent(in ) :: obs_ri(max_obs)
       real        , intent(in ) :: obs_rj(max_obs)
       real        , intent(in ) :: obs_rk(max_obs)
       real        , intent(in ) :: obs_fu(max_obs)
       real        , intent(in ) :: obs_fv(max_obs)
       real        , intent(in ) :: obs_u(max_obs)
       real        , intent(in ) :: obs_v(max_obs)
       real        , intent(in ) :: obs_weight(max_obs)
       real        , intent(in ) :: obs_vert_rad_rat(max_obs)
       real        , intent(in ) :: obs_elev(max_obs)
       real        , intent(in ) :: obs_ldf(max_obs)
       integer     , intent(in ) :: obs_mask_sea(max_obs)
       integer     , intent(in ) :: obs_i4time(max_obs)
       logical     , intent(in ) :: obs_l_withhold(max_obs)
       character*12, intent(in ) :: obs_type(max_obs)
       character*12, intent(in ) :: obs_file(max_obs)
       integer     , intent(out) :: qced_i(max_obs)
       integer     , intent(out) :: qced_j(max_obs)
       integer     , intent(out) :: qced_k(max_obs)
       real        , intent(out) :: qced_ri(max_obs)
       real        , intent(out) :: qced_rj(max_obs)
       real        , intent(out) :: qced_rk(max_obs)
       real        , intent(out) :: qced_fu(max_obs)
       real        , intent(out) :: qced_fv(max_obs)
       real        , intent(out) :: qced_u(max_obs)
       real        , intent(out) :: qced_v(max_obs)
       real        , intent(out) :: qced_weight(max_obs)
       real        , intent(out) :: qced_vert_rad_rat(max_obs)
       real        , intent(out) :: qced_elev(max_obs)
       real        , intent(out) :: qced_ldf(max_obs)
       integer     , intent(out) :: qced_mask_sea(max_obs)
       integer     , intent(out) :: qced_i4time(max_obs)
       character*12, intent(out) :: qced_type(max_obs)
       character*12, intent(out) :: qced_file(max_obs)
       integer     , intent(out) :: n_qc_total_good

       character*3 c3_string ! local
 
csms$serial(< qced_i, qced_k, qced_ri, qced_rj, qced_rk, qced_u, qced_v,
csms$>       qced_weight, qced_elev, qced_mask_sea, qced_i4time,
csms$>       qced_type, qced_file, n_qc_total_good, out>: default=ignore) begin
!
      write(6,*)' Subroutine calc_qced_obs...'

      n_qc_acrft_bad = 0
      n_qc_cdw_bad = 0
      n_qc_sfc_bad = 0
      n_qc_prof_bad = 0
      n_qc_total_bad = 0

      n_qc_acrft_good = 0
      n_qc_cdw_good = 0
      n_qc_sfc_good = 0
      n_qc_prof_good = 0
      n_qc_total_good = 0

      qc_thresh = qc_thresh_wind_def ! Threshold speed for throwing out the ob

      if(.true.)then ! l_point_struct

          do i_ob = 1,nobs_point
              i = obs_i(i_ob)
              j = obs_j(i_ob)
              k = obs_k(i_ob)

              u = obs_fu(i_ob)
              v = obs_fv(i_ob)

              speed_bkg  = sqrt(u_laps_bkg(i,j,k)**2
     1                        + v_laps_bkg(i,j,k)**2)

              u_diff = u - u_laps_bkg(i,j,k)
              v_diff = v - v_laps_bkg(i,j,k)
              speed_diff = sqrt(u_diff**2 + v_diff**2)

!             speed_thresh = max(10.,0.2 * speed_bkg)

!             Apply QC check to the OB against the background analysis
              if(
!                Make sure we actually have a real reference background
     1           (speed_bkg .gt. 0.) .and.

!                General QC check
     1           (speed_diff .gt. qc_thresh

!              Stricter QC check for pireps
     1                       .OR. 
     1         (speed_diff .gt. qc_thresh_wind_pin .and. 
     1                              obs_file(i_ob) .eq. 'pin')        

!              Stricter QC check for Cloud Drift Winds
     1                       .OR. 
     1         (speed_diff .gt. qc_thresh_wind_cdw .and. 
     1                              obs_file(i_ob) .eq. 'cdw')        

!              Stricter QC check for profilers
     1                       .OR. 
     1         (speed_diff .gt. qc_thresh_wind_pro .and. 
     1                              obs_file(i_ob) .eq. 'pro')        

!              Withhold observations for independent verification
     1                       .OR.
     1                obs_l_withhold(i_ob)
     1                                                                 )

     1                                                          )then

!                 Throw out the ob
                  if(obs_file(i_ob) .eq. 'pin')then
                      n_qc_acrft_bad = n_qc_acrft_bad + 1
                  elseif(obs_file(i_ob) .eq. 'cdw')then
                      n_qc_cdw_bad = n_qc_cdw_bad + 1
                  elseif(obs_file(i_ob) .eq. 'lso')then
                      n_qc_sfc_bad = n_qc_sfc_bad + 1
                  elseif(obs_file(i_ob) .eq. 'pro')then
                      n_qc_prof_bad = n_qc_prof_bad + 1
                  endif

                  write(6,231,err=232) obs_type(i_ob)(1:5)
     1                  ,i,j,k
     1                  ,u_diff
     1                  ,v_diff
     1                  ,u
     1                  ,v
     1                  ,u_laps_bkg(i,j,k)
     1                  ,v_laps_bkg(i,j,k)
     1                  ,speed_diff
     1                  ,obs_weight(i_ob)
231               format(a5,' QCed out - ',2i5,i4,1x,3(2x,2f5.0)
     1                                        ,f5.0,f5.2)
232               continue

                  n_qc_total_bad = n_qc_total_bad + 1

              else ! keep and write out the good OB
                  if( obs_file(i_ob) .eq. 'pin')then
                      n_qc_acrft_good = n_qc_acrft_good + 1
                      c3_string = 'Prp'
                  endif

                  if( obs_file(i_ob) .eq. 'cdw')then
                      n_qc_cdw_good  = n_qc_cdw_good + 1
                      c3_string = 'Cdw'
                  endif

                  if( obs_file(i_ob) .eq. 'lso')then
                      n_qc_sfc_good   = n_qc_sfc_good + 1
                      c3_string = 'Sfc'
                  endif

                  if( obs_file(i_ob) .eq. 'pro')then
                      n_qc_prof_good  = n_qc_prof_good + 1
                      c3_string = 'Prf'
                  endif

                  n_qc_total_good = n_qc_total_good + 1

!                 Assign array [data structure] element (using difference ob)
!                 obs_point_qced(n_qc_total_good) = obs_point(i_ob) [old way]
                  qced_i(n_qc_total_good) = obs_i(i_ob)
                  qced_j(n_qc_total_good) = obs_j(i_ob)
                  qced_k(n_qc_total_good) = obs_k(i_ob)
                  qced_ri(n_qc_total_good) = obs_ri(i_ob)
                  qced_rj(n_qc_total_good) = obs_rj(i_ob)
                  qced_rk(n_qc_total_good) = obs_rk(i_ob)
                  qced_fu(n_qc_total_good) = obs_fu(i_ob)
                  qced_fv(n_qc_total_good) = obs_fv(i_ob)
                  qced_weight(n_qc_total_good) = obs_weight(i_ob)
                  qced_vert_rad_rat(n_qc_total_good) = 
     1                                         obs_vert_rad_rat(i_ob)
                  qced_elev(n_qc_total_good) = obs_elev(i_ob)
                  qced_ldf(n_qc_total_good) = obs_ldf(i_ob)
                  qced_mask_sea(n_qc_total_good) = obs_mask_sea(i_ob)
                  qced_i4time(n_qc_total_good) = obs_i4time(i_ob)
                  qced_type(n_qc_total_good) = obs_type(i_ob)
                  qced_file(n_qc_total_good) = obs_file(i_ob)

                  qced_u(n_qc_total_good) = u_diff
                  qced_v(n_qc_total_good) = v_diff

                  if(n_qc_total_good .le. 500 .OR. 
     1               n_qc_total_good .eq. (n_qc_total_good/10)*10)then
                      iwrite = 1
                  else
                      iwrite = 0
                  endif

                  if(iwrite .eq. 1)then
                      write(6,201,err=302)c3_string,i,j,k
     1                  ,u_diff
     1                  ,v_diff
     1                  ,u
     1                  ,v
     1                  ,u_laps_bkg(i,j,k)
     1                  ,v_laps_bkg(i,j,k)
     1                  ,speed_diff
     1                  ,qced_weight(n_qc_total_good)
201                   format(1x,a3,2i5,i4,4x,f6.1,f6.1,2(2x,2f6.1)
     1                         ,f5.1,f5.2)
302                   continue
                  endif

              endif ! passed the QC test

          enddo ! i_ob

      endif

      write(6,*)
      write(6,*)' QC info for non-radar data (after remapping to grid)'
      write(6,601)n_qc_acrft_good,n_qc_acrft_bad
     1           ,pct_rejected(n_qc_acrft_good,n_qc_acrft_bad)
 601  format(' # of Aircraft   GOOD/BAD QC = ',2i6,7x
     1      ,'% rejected = ',f6.1)

      write(6,602)n_qc_cdw_good,n_qc_cdw_bad
     1           ,pct_rejected(n_qc_cdw_good,n_qc_cdw_bad)
 602  format(' # of CDWs       GOOD/BAD QC = ',2i6,7x
     1      ,'% rejected = ',f6.1)

      write(6,603)n_qc_sfc_good,n_qc_sfc_bad
     1           ,pct_rejected(n_qc_sfc_good,n_qc_sfc_bad)
 603  format(' # of SFC        GOOD/BAD QC = ',2i6,7x
     1      ,'% rejected = ',f6.1)

      write(6,604)n_qc_prof_good,n_qc_prof_bad
     1           ,pct_rejected(n_qc_prof_good,n_qc_prof_bad)
 604  format(' # of PROFs      GOOD/BAD QC = ',2i6,7x
     1      ,'% rejected = ',f6.1)

      write(6,605)n_qc_total_good,n_qc_total_bad
     1           ,pct_rejected(n_qc_total_good,n_qc_total_bad)
 605  format(/' # of Non-Radar  GOOD/BAD QC = ',2i6,7x
     1       ,'% rejected = ',f6.1)

      I4_elapsed = ishow_timer()

csms$serial end

      return
      end
