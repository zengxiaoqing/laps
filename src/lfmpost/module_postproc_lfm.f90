!dis   
!dis    Open Source License/Disclaimer, Forecast Systems Laboratory
!dis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
!dis    
!dis    This software is distributed under the Open Source Definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis    
!dis    In particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis    
!dis    - Redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis    
!dis    - Redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis    
!dis    - All modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis    
!dis    - If significant modifications or enhancements are made to this
!dis    software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis    
!dis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
!dis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
!dis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
!dis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
!dis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
!dis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
!dis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
!dis   
!dis 

MODULE postproc_lfm

  ! Contains routines necessary to post-process MM5v3 data file for LAPS.

  ! Brent Shaw, NOAA/OAR/FSL/FRD,  Dec 2000

  USE mm5v3_io
  USE setup
  USE constants
  USE vinterp_utils
  USE fire

  IMPLICIT NONE

  PRIVATE

  INTEGER                    :: current_lun
  CHARACTER(LEN=24)          :: time_to_proc
  INTEGER                    :: i,j,k
 
  ! Some variables needed on sigma for various derivations

  REAL, ALLOCATABLE          :: psig               ( : , : , : )
  REAL, ALLOCATABLE          :: tsig               ( : , : , : )
  REAL, ALLOCATABLE          :: tvsig              ( : , : , : )
  REAL, ALLOCATABLE          :: thetasig           ( : , : , : )
  REAL, ALLOCATABLE          :: mrsig              ( : , : , : )
  REAL, ALLOCATABLE          :: zsig               ( : , : , : )
  REAL, ALLOCATABLE          :: rhsig              ( : , : , : )
  REAL, ALLOCATABLE          :: rhodrysig          ( : , : , : )
  REAL, ALLOCATABLE          :: rhomoistsig        ( : , : , : )
  REAL, ALLOCATABLE          :: tvprs              ( : , : , : )
  REAL, ALLOCATABLE          :: mrprs              ( : , : , : )
  ! Some variables to hold interpolation coefficients
  INTEGER, ALLOCATABLE       :: trap_bot_ind       ( : , : , : )
  INTEGER, ALLOCATABLE       :: trap_top_ind       ( : , : , : )
  REAL, ALLOCATABLE          :: weight_top_lin     ( : , : , : )
  REAL, ALLOCATABLE          :: weight_top_log     ( : , : , : )

  ! Output variables are on pressure or are at the surface
  REAL, ALLOCATABLE, PUBLIC  :: tprs               ( : , : , : )
  REAL, ALLOCATABLE, PUBLIC  :: tdprs              ( : , : , : )
  REAL, ALLOCATABLE, PUBLIC  :: thetaprs           ( : , : , : )
  REAL, ALLOCATABLE, PUBLIC  :: uprs               ( : , : , : )
  REAL, ALLOCATABLE, PUBLIC  :: vprs               ( : , : , : )
  REAL, ALLOCATABLE, PUBLIC  :: wprs               ( : , : , : )
  REAL, ALLOCATABLE, PUBLIC  :: omprs              ( : , : , : )
  REAL, ALLOCATABLE, PUBLIC  :: shprs              ( : , : , : )
  REAL, ALLOCATABLE, PUBLIC  :: rhprs              ( : , : , : )
  REAL, ALLOCATABLE, PUBLIC  :: zprs               ( : , : , : )
  REAL, ALLOCATABLE, PUBLIC  :: cldliqmr_prs       ( : , : , : )
  REAL, ALLOCATABLE, PUBLIC  :: cldicemr_prs       ( : , : , : )
  REAL, ALLOCATABLE, PUBLIC  :: snowmr_prs         ( : , : , : )
  REAL, ALLOCATABLE, PUBLIC  :: rainmr_prs         ( : , : , : )
  REAL, ALLOCATABLE, PUBLIC  :: graupelmr_prs      ( : , : , : )
  REAL, ALLOCATABLE, PUBLIC  :: refl_prs           ( : , : , : )
  REAL, ALLOCATABLE, PUBLIC  :: pcptype_prs        ( : , : , : )
  REAL, ALLOCATABLE, PUBLIC  :: abs_vort           ( : , : , : )
  REAL, ALLOCATABLE, PUBLIC  :: tkeprs             ( : , : , : )
  REAL, ALLOCATABLE, PUBLIC  :: psfc               ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: pmsl               ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: redp               ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: tsfc               ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: tdsfc              ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: usfc               ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: vsfc               ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: upbl               ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: vpbl               ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: wsfc               ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: rhsfc              ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: cldbase            ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: cldtop             ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: cldamt             ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: ceiling            ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: heatind            ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: intliqwater        ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: totpcpwater        ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: pcp_inc            ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: pcp_init           ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: con_pcp_inc        ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: con_pcp_init       ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: snow_inc           ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: snow_init          ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: pcp_tot            ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: con_pcp_tot        ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: snow_tot           ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: pcptype_sfc        ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: thetasfc           ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: thetaesfc          ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: cape               ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: cin                ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: liftedind          ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: srhel              ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: max_refl           ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: echo_tops          ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: refl_sfc           ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: visibility         ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: thick_10_5         ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: snowcover          ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: lwout              ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: swout              ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: lwdown             ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: swdown             ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: albedo             ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: shflux             ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: lhflux             ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: pblhgt             ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: ground_t           ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: clwmrsfc           ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: icemrsfc           ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: rainmrsfc          ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: snowmrsfc          ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: graupmrsfc         ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: vnt_index          ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: ham_index          ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: hah_index          ( : , : )
  REAL, ALLOCATABLE, PUBLIC  :: fwi_index          ( : , : )



  REAL, PARAMETER            :: nonecode = 0.
  REAL, PARAMETER            :: raincode = 1.
  REAL, PARAMETER            :: snowcode = 2.
  REAL, PARAMETER            :: zraincode =3.
  REAL, PARAMETER            :: sleetcode = 4.
  REAL, PARAMETER            :: hailcode = 5.
  REAL, PARAMETER            :: drizzlecode = 6.
  REAL, PARAMETER            :: zdrizzlecode = 7.
  REAL, PARAMETER            :: rainsnowcode = 8.
  REAL, PARAMETER            :: rainicecode = 9.

  INTEGER                    :: k300, k500, k700, k850, k1000

  REAL, EXTERNAL             :: dewpt
  REAL, EXTERNAL             :: relhum
  REAL, EXTERNAL             :: potential_temp
  REAL, EXTERNAL             :: eq_potential_temp
  LOGICAL                    :: initialize 
  REAL                       :: smth
  PUBLIC process_one_lfm
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE process_one_lfm (lun,time)
  
  ! Main driver subroutine for processing one time of MM5v3
    IMPLICIT NONE
    INTEGER, INTENT(IN)           :: lun
    CHARACTER(LEN=24),INTENT(IN)  :: time

    current_lun  = lun
   
    ! Depending on time step used in the model, there can be
    ! some rounding errors in the computation of the time
    ! string in the subheaders of the model output (i.e., it can
    ! be off by 1 second so it does not exactly match.  Since the
    ! I/O routine can be called with a time string less than
    ! 1960-01-01 to tell it to use whatever it finds, we can
    ! excercise this option when using split output, since we
    ! know by the file number which output time is contained 
    ! within.
    IF (split_output) THEN
      time_to_proc = '0000-00-00_00:00:00.0000'
    ELSE
      time_to_proc = time
    ENDIF
    
    ! Find some key indices of pressure levels for use in some routines
    DO k = 1 , kprs
      IF (prslvl(k) .EQ. 100000) k1000 = k
      IF (prslvl(k) .EQ. 85000) k850 = k
      IF (prslvl(k) .EQ. 70000) k700 = k
      IF (prslvl(k) .EQ. 50000) k500 = k
      IF (prslvl(k) .EQ. 30000) k300 = k
    ENDDO
 
    PRINT '(A)', 'PROCESS_ONE_LFM: Allocating internal arrays...'
    CALL allocate_internal
    PRINT '(A)', 'PROCESS_ONE_LFM: Getting state variables on sigma...'
    CALL get_sigma_vars
    PRINT '(A)', 'PROCESS_ONE_LFM: Interp thermo variables to pressure...'
    CALL interp_thermo_3d
    PRINT '(A)', 'PROCESS_ONE_LFM: Interp momentum variables to pressure...'
    CALL interp_winds
    PRINT '(A)', 'PROCESS_ONE_LFM: Processing clouds and reflectivity...'
    CALL get_clouds_reflectivity
    PRINT '(A)', 'PROCESS_ONE_LFM: Getting precipitation amounts...'
    CALL get_precip
    PRINT '(A)', 'PROCESS_ONE_LFM: Computing stability indices...'
    CALL make_stability
    PRINT '(A)', 'PROCESS_ONE_LFM: Miscellaneous derivations...'
    CALL make_misc
    PRINT '(A)', 'PROCESS_ONE_LFM: Deallocating internal arrays...'
    CALL deallocate_internal
  
  END SUBROUTINE process_one_lfm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE allocate_internal

    ! This routine allocates all of the arrays internal to this module that
    ! may be shared among more than one routine in this module.  It also
    ! allocates the public (output) arrays if they are not already allocated
    IMPLICIT NONE
    
    ALLOCATE ( psig            ( nx , ny , ksigh ) )
    ALLOCATE ( tsig            ( nx , ny , ksigh ) )
    ALLOCATE ( tvsig           ( nx , ny , ksigh ) )
    ALLOCATE ( thetasig        ( nx , ny , ksigh ) )
    ALLOCATE ( mrsig           ( nx , ny , ksigh ) )
    ALLOCATE ( rhsig           ( nx , ny , ksigh ) )
    ALLOCATE ( rhodrysig       ( nx , ny , ksigh ) ) 
    ALLOCATE ( rhomoistsig     ( nx , ny , ksigh ) )
    ALLOCATE ( mrprs           ( nx , ny , kprs ) )
    ALLOCATE ( tvprs           ( nx , ny,  kprs ) )

    IF (.NOT. ALLOCATED(zsig) ) THEN

      ! The non-hydrostatic version of MM5 (only version now
      ! supported since version 3) has constant height values
      ! for each sigma based on the reference state.  Thus, they 
      ! only need to be allocated/computed the first time through.

      initialize = .true.
      ALLOCATE ( zsig          ( nx , ny , ksigh ) )
    ELSE
      initialize = .false.
    ENDIF
    ALLOCATE ( trap_top_ind    ( nx , ny , kprs ) )
    ALLOCATE ( trap_bot_ind    ( nx , ny , kprs ) )
    ALLOCATE ( weight_top_lin  ( nx , ny , kprs ) )
    ALLOCATE ( weight_top_log  ( nx , ny , kprs ) )

    ! Allocation of public variables
    IF (.NOT.ALLOCATED(psfc))          ALLOCATE ( psfc         ( nx , ny ) )
    IF (.NOT.ALLOCATED(tsfc))          ALLOCATE ( tsfc         ( nx , ny ) )
    IF (.NOT.ALLOCATED(thetasfc))      ALLOCATE ( thetasfc     ( nx , ny ) )
    IF (.NOT.ALLOCATED(thetaesfc))     ALLOCATE ( thetaesfc    ( nx , ny ) )
    IF (.NOT.ALLOCATED(rhsfc))         ALLOCATE ( rhsfc        ( nx , ny ) )
    IF (.NOT.ALLOCATED(tdsfc))         ALLOCATE ( tdsfc        ( nx , ny ) )
    IF (.NOT.ALLOCATED(thetaprs))      ALLOCATE ( thetaprs     ( nx , ny , kprs ) )
    IF (.NOT.ALLOCATED(zprs))          ALLOCATE ( zprs         ( nx , ny , kprs ) )
    IF (.NOT.ALLOCATED(rhprs))         ALLOCATE ( rhprs        ( nx , ny , kprs ) )
    IF (.NOT.ALLOCATED(tprs))          ALLOCATE ( tprs         ( nx , ny , kprs ) )
    IF (.NOT.ALLOCATED(tdprs))         ALLOCATE ( tdprs        ( nx , ny , kprs ) )
    IF (.NOT.ALLOCATED(shprs))         ALLOCATE ( shprs        ( nx , ny , kprs ) )
    IF (.NOT.ALLOCATED(redp))        ALLOCATE ( redp         ( nx , ny ) )
    IF (.NOT.ALLOCATED(pmsl))          ALLOCATE ( pmsl         ( nx , ny ) )
    IF (.NOT.ALLOCATED(usfc))          ALLOCATE ( usfc         ( nx , ny ) )
    IF (.NOT.ALLOCATED(uprs))          ALLOCATE ( uprs         ( nx , ny , kprs ) )
    IF (.NOT.ALLOCATED(vsfc))          ALLOCATE ( vsfc         ( nx , ny ) )
    IF (.NOT.ALLOCATED(vprs))          ALLOCATE ( vprs         ( nx , ny , kprs ) )
    IF (.NOT.ALLOCATED(wsfc))          ALLOCATE ( wsfc         ( nx , ny ) )
    IF (.NOT.ALLOCATED(wprs))          ALLOCATE ( wprs         ( nx , ny , kprs ) )
    IF (.NOT.ALLOCATED(omprs))         ALLOCATE ( omprs        ( nx , ny , kprs ) )
    IF (.NOT.ALLOCATED(cldbase))       ALLOCATE ( cldbase      ( nx , ny ) )
    IF (.NOT.ALLOCATED(cldtop))        ALLOCATE ( cldtop       ( nx , ny ) )   
    IF (.NOT.ALLOCATED(cldamt))        ALLOCATE ( cldamt       ( nx , ny ) )
    IF (.NOT.ALLOCATED(ceiling))       ALLOCATE ( ceiling      ( nx , ny ) )
    IF (.NOT.ALLOCATED(intliqwater))   ALLOCATE ( intliqwater  ( nx , ny ) )
    IF (.NOT.ALLOCATED(totpcpwater))   ALLOCATE ( totpcpwater  ( nx , ny ) )
    IF (.NOT.ALLOCATED(max_refl))      ALLOCATE ( max_refl     ( nx , ny ) )
    IF (.NOT.ALLOCATED(echo_tops))     ALLOCATE ( echo_tops    ( nx , ny ) )
    IF (.NOT.ALLOCATED(cldliqmr_prs))  ALLOCATE ( cldliqmr_prs ( nx , ny , kprs ) )
    IF (.NOT.ALLOCATED(cldicemr_prs))  ALLOCATE ( cldicemr_prs ( nx , ny , kprs ) )
    IF (.NOT.ALLOCATED(rainmr_prs))    ALLOCATE ( rainmr_prs   ( nx , ny , kprs ) )
    IF (.NOT.ALLOCATED(snowmr_prs))    ALLOCATE ( snowmr_prs   ( nx , ny , kprs ) )
    IF (.NOT.ALLOCATED(graupelmr_prs)) ALLOCATE ( graupelmr_prs( nx , ny , kprs ) )
    IF (.NOT.ALLOCATED(refl_prs))      ALLOCATE ( refl_prs     ( nx , ny , kprs ) )
    IF (.NOT.ALLOCATED(refl_sfc))      ALLOCATE ( refl_sfc     ( nx , ny ) )
    IF (.NOT.ALLOCATED(pcptype_sfc))   ALLOCATE ( pcptype_sfc  ( nx , ny ) )
    IF (.NOT.ALLOCATED(pcptype_prs))   ALLOCATE ( pcptype_prs  ( nx , ny , kprs) )
    IF (.NOT.ALLOCATED(pcp_init))      ALLOCATE ( pcp_init     ( nx , ny ) )
    IF (.NOT.ALLOCATED(pcp_inc))       ALLOCATE ( pcp_inc      ( nx , ny ) )
    IF (.NOT.ALLOCATED(pcp_tot))       ALLOCATE ( pcp_tot      ( nx , ny ) )
    IF (.NOT.ALLOCATED(con_pcp_init))  ALLOCATE ( pcp_init     ( nx , ny ) )
    IF (.NOT.ALLOCATED(con_pcp_inc))   ALLOCATE ( pcp_inc      ( nx , ny ) )
    IF (.NOT.ALLOCATED(con_pcp_tot))   ALLOCATE ( pcp_tot      ( nx , ny ) )
    IF (.NOT.ALLOCATED(snow_init))     ALLOCATE ( snow_init    ( nx , ny ) )
    IF (.NOT.ALLOCATED(snow_inc))      ALLOCATE ( snow_inc     ( nx , ny ) )
    IF (.NOT.ALLOCATED(snow_tot))      ALLOCATE ( snow_tot     ( nx , ny ) )
    IF (.NOT.ALLOCATED(srhel))         ALLOCATE ( srhel        ( nx , ny ) )
    IF (.NOT.ALLOCATED(cape))          ALLOCATE ( cape         ( nx , ny ) )
    IF (.NOT.ALLOCATED(cin))           ALLOCATE ( cin          ( nx , ny ) )
    IF (.NOT.ALLOCATED(liftedind) )    ALLOCATE ( liftedind    ( nx , ny ) )
    IF (.NOT.ALLOCATED(visibility) )   ALLOCATE ( visibility   ( nx , ny ) )
    IF (.NOT.ALLOCATED(heatind) )      ALLOCATE ( heatind      ( nx , ny ) )
    IF (.NOT.ALLOCATED(lwout) )        ALLOCATE (lwout      ( nx , ny ) )
    IF (.NOT.ALLOCATED(swout) )        ALLOCATE (swout      ( nx , ny ) )
    IF (.NOT.ALLOCATED(lwdown) )        ALLOCATE (lwdown      ( nx , ny ) )
    IF (.NOT.ALLOCATED(swdown) )        ALLOCATE (swdown      ( nx , ny ) )
    IF (.NOT.ALLOCATED(albedo) )        ALLOCATE (albedo      ( nx , ny ) )
    IF (.NOT.ALLOCATED(shflux))        ALLOCATE (shflux  ( nx , ny ) )
    IF (.NOT.ALLOCATED(lhflux) )       ALLOCATE (lhflux      ( nx , ny ) )
    IF (.NOT.ALLOCATED(pblhgt) )       ALLOCATE (pblhgt      ( nx , ny ) )
    IF (.NOT.ALLOCATED(ground_t) )     ALLOCATE (ground_t    ( nx , ny ) )
    IF (.NOT.ALLOCATED(clwmrsfc) )     ALLOCATE (clwmrsfc    ( nx , ny ) )
    IF (.NOT.ALLOCATED(icemrsfc) )     ALLOCATE (icemrsfc    ( nx , ny ) )
    IF (.NOT.ALLOCATED(rainmrsfc) )    ALLOCATE (rainmrsfc   ( nx , ny ) )
    IF (.NOT.ALLOCATED(snowmrsfc) )    ALLOCATE (snowmrsfc   ( nx , ny ) )
    IF (.NOT.ALLOCATED(graupmrsfc) )   ALLOCATE (graupmrsfc  ( nx , ny ) )
    IF ((.NOT.ALLOCATED(abs_vort)).AND.(make_v5d(domain_num))) &
         ALLOCATE ( abs_vort (nx, ny, kprs) )
    IF ((.NOT.ALLOCATED(thick_10_5)).AND.(make_v5d(domain_num))) &
         ALLOCATE ( thick_10_5 (nx,ny))
    IF ((.NOT.ALLOCATED(snowcover)).AND.(make_v5d(domain_num))) &
         ALLOCATE ( snowcover (nx,ny))
    RETURN

  END SUBROUTINE allocate_internal    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_sigma_vars

    ! This subroutine populates some of the state variables needed on 
    ! sigma by reading them in and/or deriving them
   
    IMPLICIT NONE
    REAL, ALLOCATABLE             :: pstar ( : , : )
    REAL, ALLOCATABLE             :: ppsig ( : , : , : )
    REAL, ALLOCATABLE             :: pbase ( : , : , : )
    REAL, ALLOCATABLE             :: phb   ( : , : , : )
    REAL, ALLOCATABLE             :: ph    ( : , : , : )
    REAL, ALLOCATABLE             :: zsigf ( : , : , : )
    REAL, ALLOCATABLE             :: mu    ( : , : )
    REAL, ALLOCATABLE             :: mub   ( : , : )
    INTEGER                       :: status
    REAL                          :: tvbar
    REAL, EXTERNAL                :: mixsat 
    REAL                          :: dz, dtdz, pbot,tvbot,zbot     
    REAL, ALLOCATABLE             :: mrsfc(:,:)
    REAL, EXTERNAL                :: tcvp
    LOGICAL                       :: made_pbl

    made_pbl = .false.
    ! Get the 3D pressure
 
    IF (mtype .EQ. 'mm5') THEN
      ! Compute pressure on sigma and sfc from the perturbation pressure and
      ! pstar

      ALLOCATE ( pstar ( nx , ny ) )
      ALLOCATE ( ppsig ( nx , ny , ksigh ) )

      CALL get_mm5_3d(current_lun, 'PP       ', time_to_proc, ppsig, &
                        'D    ', status)
      IF (status.NE.0) THEN
        PRINT '(A)', 'Problem getting peturbation pressure!  Aborting...'
        CALL abort
      ENDIF
      CALL get_mm5_2d(current_lun, 'PSTARCRS ', time_to_proc, pstar, &
                        'D   ', status)
      IF (status.NE.0) THEN
        PRINT '(A)', 'Problem getting pstar!  Aborting...'
        CALL abort
      ENDIF

      IF (do_smoothing) THEN
        DO smth = 0.5, -0.5, -1
          CALL SMOOTH(ppsig,nx,ny,ksigh,smth)
          CALL SMOOTH(pstar,nx,ny,1,smth)
        ENDDO
      ENDIF

      DO k = 1, ksigh
        psig(:,:,k) = pstar(:,:)*sigmah(k) + Ptop + ppsig(:,:,k)
      ENDDO
      psfc = pstar + ppsig(:,:,1) + Ptop
    
     ! Free up memory by deallocating local var ppsig and pstar

      DEALLOCATE (ppsig)
      DEALLOCATE (pstar)
    
      ! Get the temperature data
      CALL get_mm5_3d(current_lun, 'T        ', time_to_proc, tsig, &
                      'D    ', status)
      IF (status.NE.0) THEN
        PRINT '(A)', 'Problem getting temperature!  Aborting...'
        CALL abort
      ENDIF
    
      ! Get the mixing ratio data
      CALL get_mm5_3d(current_lun, 'Q        ', time_to_proc, mrsig, &
                      'D    ', status)

      IF (status.NE.0) THEN
        PRINT '(A)', 'Problem getting mixing ratio!  Aborting...'
        CALL abort
      ENDIF

      IF (do_smoothing) THEN
        DO smth = 0.5, -0.5, -1
          CALL SMOOTH(tsig,nx,ny,ksigh,smth)
          CALL SMOOTH(mrsig,nx,ny,ksigh,smth)
        ENDDO
      ENDIF
      WHERE(mrsig .LE. 0.) mrsig = 0.000001
    
      ! Compute Tv on sigma
      PRINT *, 'Computing virtual temperature on sigma layers...'
      tvsig = tsig * ( 1. + 0.61 * mrsig ) 
      PRINT '(A,2F7.1)', '   MIN/MAX = ', minval(tvsig),maxval(tvsig)
      ! Compute theta on sigma
      PRINT *, 'Computing theta on sigma layers...'

      thetasig = tsig * (100000./psig)**kappa
      PRINT '(A,2F7.1)', '   MIN/MAX = ', minval(thetasig),maxval(thetasig)
      rhodrysig = psig / (R * tsig)
      rhomoistsig = psig / (R * tvsig)

      ! Heights are computed using the hypsometric equation, starting with
      ! the surface pressure and terrain height and integrating upward
      ! one level at a time.  The earlier versions of this program computed
      ! the heights as a static field based on the MM5 base state per the
      ! documentation.  However, sometimes the heights at the first sigma
      ! level were slightly below ground when computed this way.  This new
      ! way gives nearly identical results, but the heights are never below
      ! ground on any sigma surface.

      PRINT *, 'Computing heights on sigma layers...'
      DO j = 1, ny
        DO i = 1, nx
          ! Initialize the values for the bottom of the layer we 
          ! are going to apply the hypsometric equation to.
          zbot = terdot(i,j)
          pbot = psfc(i,j)
          tvbot = tvsig(i,j,1)
          DO k = 1,ksigh
            ! Compute mean virtual temp for this layer
            tvbar = 0.5*(tvsig(i,j,k) + tvbot)
            ! Use hypsometric equation to compute thickness
            dz = r*tvbar*ALOG(pbot/psig(i,j,k))/grav
            ! Add thickness to height of bottom of layer
            zsig(i,j,k) = zbot + dz

            ! Set the bottom level to the current level for next go-around.
            zbot = zsig(i,j,k)
            pbot = psig(i,j,k)
            tvbot = tvsig(i,j,k)
          ENDDO
        ENDDO
      ENDDO
 
    ELSEIF (mtype(1:3).EQ.'wrf') THEN
      ! Get pressure on sigma...we do this by getting the base state
      ! and perturbation pressures and adding them together
      ALLOCATE(pbase(nx,ny,ksigh))
      ALLOCATE(ppsig(nx,ny,ksigh))
      CALL get_wrfnc_3d(current_lun, "PB","A",nx,ny,ksigh,1,pbase,status)
      IF (status.NE.0) THEN
        PRINT *, 'Could not properly obtain WRF base state pressure.'
        CALL ABORT
      ENDIF
      CALL get_wrfnc_3d(current_lun, "P", "A",nx,ny,ksigh,1,ppsig,status)
      IF (status.NE.0) THEN
        PRINT *, 'Could not properly obtain WRF perturbation pressure.'
        CALL ABORT
      ENDIF
      psig = pbase+ppsig
      DEALLOCATE(pbase)
      DEALLOCATE(ppsig) 
      IF (do_smoothing) THEN
        DO smth = 0.5, -0.5, -1
          CALL SMOOTH(psig,nx,ny,ksigh,smth)
        ENDDO
      ENDIF
      ! Get theta on sigma
      CALL get_wrfnc_3d(current_lun, "T","A",nx,ny,ksigh,1,thetasig,status)
      IF (status.NE.0) THEN
        PRINT *, 'Could not properly obtain WRF perturbation theta.'
        CALL ABORT
      ENDIF 
      thetasig = thetasig + 300.
      IF (do_smoothing) THEN
        DO smth = 0.5, -0.5, -1
          CALL SMOOTH(thetasig,nx,ny,ksigh,smth)
        ENDDO
      ENDIF

      ! Get Q on sigma
      CALL get_wrfnc_3d(current_lun, "QVAPOR","A",nx,ny,ksigh,1,mrsig,status)
      IF (status.NE.0) THEN
        PRINT *, 'Could not properly obtain WRF mixing ratio.'
        CALL ABORT
      ENDIF
      IF (do_smoothing) THEN
        DO smth = 0.5, -0.5, -1
          CALL SMOOTH(mrsig,nx,ny,ksigh,smth)
        ENDDO
      ENDIF
      ! Compute temperature on sigma
      tsig = thetasig/ ( (100000./psig)**kappa) 
 
      ! Compute dry density
      rhodrysig = psig / (R * tsig)
  
      ! Compute virtual temperature on sigma
      tvsig = tsig * (1. + 0.61 * mrsig)
  
      ! Get GPH and destagger vertically
      ALLOCATE( ph (nx,ny,ksigf) )
      ALLOCATE( phb(nx,ny,ksigf) )
      ALLOCATE( zsigf(nx,ny,ksigf) )
      CALL get_wrfnc_3d(current_lun, "PH","A",nx,ny,ksigf,1,ph,status)
      IF (status.NE.0) THEN
        PRINT *, 'Could not properly obtain WRF perturbation geopotential.'
        CALL ABORT
      ENDIF
      CALL get_wrfnc_3d(current_lun, "PHB","A",nx,ny,ksigf,1,phb,status)
      IF (status.NE.0) THEN
        PRINT *, 'Could not properly obtain WRF base-state geopotential.'
        CALL ABORT
      ENDIF

      zsigf = (ph + phb)/grav
      DO k=1,ksigh
        zsig(:,:,k) =  0.5*(zsigf(:,:,k)+zsigf(:,:,k+1))
      ENDDO
      
      !  Check lowest level heights
      do j=1,ny
        do i=1,nx
          if (zsig(i,j,1).LE.terdot(i,j))then
            print*, 'zsig1 < ter',i, j,zsig(i,j,1),terdot(i,j)
            stop
          endif
        enddo
      enddo
      DEALLOCATE(ph)
      DEALLOCATE(phb)
      DEALLOCATE(zsigf)
      IF (do_smoothing) THEN
        DO smth = 0.5, -0.5, -1
          CALL SMOOTH(zsig,nx,ny,ksigh,smth)
        ENDDO
      ENDIF

      ! Compute moist density on sigma
      rhomoistsig = psig / (R * tvsig)

      ! Compute psfc
         ! Get sfc dry pressure (mu+mub+ptop)
      ALLOCATE(mu(nx,ny))
      ALLOCATE(mub(nx,ny))
      CALL get_wrfnc_2d(current_lun, "MUB","A",nx,ny,1,mub,status) 
      IF (status.NE.0) THEN
        PRINT *, 'Could not properly obtain WRF base-state mu.'
        CALL ABORT
      ENDIF
      CALL get_wrfnc_2d(current_lun, "MU","A",nx,ny,1,mu,status)
      IF (status.NE.0) THEN
        PRINT *, 'Could not properly obtain WRF perturbation mu.'
        CALL ABORT
      ENDIF
      psfc = mu + mub 
      psfc = psfc + ptop
      DEALLOCATE(mu)
      DEALLOCATE(mub)
         ! Compute integrated column vapor pressure
      DO j = 1,ny
        DO i = 1,nx
          psfc(i,j) = psfc(i,j) + tcvp(psig(i,j,:), &
                                 mrsig(i,j,:),zsig(i,j,:), &
                                 rhomoistsig(i,j,:),ksigh)
        ENDDO
      ENDDO
      IF (do_smoothing) THEN
        DO smth = 0.5, -0.5, -1
          CALL SMOOTH(psfc,nx,ny,1,smth)
        ENDDO
      ENDIF

    ENDIF


    ! Compute RH on sigma levels, as it will be better to interpolate
    ! RH vertically than mixing ratio

    PRINT *, 'Computing RH on sigma...'
    DO k = 1, ksigh
      DO j = 1,ny
        DO i=1,nx
           rhsig(i,j,k) = relhum(tsig(i,j,k),mrsig(i,j,k),psig(i,j,k))*100.
        ENDDO
      ENDDO
    ENDDO
    PRINT *, '    MIN/MAX = ', minval(rhsig), maxval(rhsig)
    WHERE (rhsig .gt. 100.) rhsig = 100.
    WHERE (rhsig .lt. 1.) rhsig = 1.   
   ! Print diagnostics from center column
    PRINT '(A)','Diagnostics from domain center:'
    PRINT '(A,F7.0)', 'Terrain height at center = ', terdot(nx/2,ny/2)
    PRINT '(A)', &
       '----------------------------------------------------------------------'
    PRINT '(A)', &
       'LEVEL   SIGMA  PRESSURE(Pa)  HEIGHT     T      Tv     THETA   RH   QV'
    PRINT '(A)', &
       '----------------------------------------------------------------------'
    DO k = 1, ksigh
      PRINT '(I5,2x,F6.4,2x,F12.1,2x,F6.0,2x,F6.2,2x,F6.2,2x,F6.2,2x,F4.0,2x,F8.6)', &
               k, sigmah(k),psig(nx/2,ny/2,k),zsig(nx/2,ny/2,k), &
               tsig(nx/2,ny/2,k),tvsig(nx/2,ny/2,k),thetasig(nx/2,ny/2,k), &
               rhsig(nx/2,ny/2,k),mrsig(nx/2,ny/2,k)
   
    ENDDO 

    ! Get surface temperature and moisture.  MM5v3 may have T2 and Q2
    ! variables from similarity theory.  WRFv1 has TH2 and Q2 variables.

    ALLOCATE(mrsfc(nx,ny))
    tsfc(:,:) = 0.
    mrsfc(:,:)= 0.

    ! 2m temperature
    IF (mtype.EQ. 'mm5') THEN
      CALL get_mm5_2d(current_lun, 'T2       ', time_to_proc, tsfc, &
                      'D   ', status)

    ELSEIF(mtype(1:3).EQ.'wrf')THEN
   
       ! Get TH2 and convert to temperature if non-zero
       CALL get_wrfnc_2d(current_lun,'TH2','A',nx,ny,1,tsfc,status)
       thetasfc = tsfc
       IF ((MINVAL(tsfc) .GT. 100.) .AND. (MAXVAL(tsfc) .LT. 1000.)) THEN
         tsfc = tsfc/((100000./psfc)**kappa)
       ELSE 
         tsfc = 0. 
         thetasfc = 0.
       ENDIF
    ENDIF
    
    ! 2m qvapor
    IF (mtype.EQ. 'mm5') THEN
      CALL get_mm5_2d(current_lun, 'Q2       ', time_to_proc, mrsfc, &
                      'D   ', status)

    ELSEIF(mtype(1:3).EQ.'wrf')THEN
   
       ! Get Q2
      CALL get_wrfnc_2d(current_lun,'Q2','A',nx,ny,1,mrsfc,status)

    ENDIF 

    ! Make sure we had Q2.  If not, use lowest sigma value

    IF (MAXVAL(mrsfc).LT. 0.000001)THEN
      print *, 'Using lowest model level mixing ratio for surface.'
      mrsfc(:,:) = mrsig(:,:,1)
    ELSE 
      IF (do_smoothing) THEN
        DO smth = 0.5, -0.5, -1
          CALL SMOOTH(mrsfc,nx,ny,1,smth)
        ENDDO
      ENDIF 
    ENDIF

!DIAGNOSTIC
dz = zsig(nx/2,ny/2,1) - (terdot(nx/2,ny/2) + 2.)
dtdz = ( tsig(nx/2,ny/2,2) - tsig(nx/2,ny/2,1) ) / &
       ( zsig(nx/2,ny/2,2) - zsig(nx/2,ny/2,1) )
print '(A,4F6.1,F10.5)','SFCTEMPTEST:T1 Tsim Texp DZ DTDZ =',tsig(nx/2,ny/2,1),&
  tsfc(nx/2,ny/2),tsig(nx/2,ny/2,1)-dtdz*dz,dz,dtdz

    ! Make sure we have tsfc
    IF (MAXVAL(tsfc) .LT. 150.) THEN
      PRINT '(A)', 'T2 not available...will extrapolate from lowest level'

      ! Assume a constant mixing ratio
      ! Loop over all horizontal points for extrapolation

      DO j = 1 , ny
        DO i = 1 , nx

          ! Compute dz between lowest sigma layer and 2m level
          ! Added parenthesis around 2nd Term -- BLS 27 Sep 01
          dz = zsig(i,j,1) - (terdot(i,j) + 2.) 
          ! Compute the lapse rate of temp for the two lowest
          ! sigma levels.  Use of the hypsometric equation
          ! did not work well for these thin layers.

          dtdz = ( tsig(i,j,2) - tsig(i,j,1) ) / &
                 ( zsig(i,j,2) - zsig(i,j,1) ) 
          tsfc(i,j) = tsig(i,j,1) - dtdz*dz
          IF (mtype(1:3) .EQ. 'wrf') THEN
            thetasfc(i,j) = potential_temp(tsfc(i,j),psfc(i,j))
          ENDIF
        ENDDO
      ENDDO
    ELSE 
      IF (do_smoothing) THEN
        DO smth = 0.5, -0.5, -1
          CALL SMOOTH(tsfc,nx,ny,1,smth)
        ENDDO
      ENDIF
    ENDIF
    ! Compute some things that are derived from T and q
    print *, 'Computing RH...'
    DO j = 1 , ny

      DO i = 1 , nx     
     
       ! Compute sfc relative humidity
       rhsfc(i,j) = MIN(relhum(tsfc(i,j),mrsfc(i,j),psfc(i,j)),1.)*100.
       

       ! Compute sfc dewpoint
       tdsfc(i,j) = dewpt(tsfc(i,j),rhsfc(i,j)*0.01)

       IF (mtype(1:3).NE.'wrf') THEN
         ! Compute theta at the surface
         thetasfc(i,j) =  potential_temp(tsfc(i,j),psfc(i,j))
       ENDIF
       ! Compute thetae at the surface
       thetaesfc(i,j) = eq_potential_temp(tsfc(i,j),psfc(i,j),mrsfc(i,j), &
                        rhsfc(i,j)*0.01)

     ENDDO
   ENDDO
   PRINT *, '   MIN/MAX tsfc      = ', minval(tsfc),maxval(tsfc) 
   PRINT *, '   MIN/MAX rhsfc     = ', minval(rhsfc),maxval(rhsfc)   
   PRINT *, '   MIN/MAX tdsfc     = ', minval(tdsfc),maxval(tdsfc) 
   PRINT *, '   MIN/MAX thetasfc  = ', minval(thetasfc),maxval(thetasfc)
   PRINT *, '   MIN/MAX thetaesfc = ', minval(thetaesfc),maxval(thetaesfc)

   PRINT *, 'Diagnostic from lowest sigma layer at domain center:'
   PRINT '("P:",F8.1," T:",F7.1," MR:",F6.4," RH:",F5.1)', &
     psig(nx/2,ny/2,1),tsig(nx/2,ny/2,1),mrsig(nx/2,ny/2,1),rhsig(nx/2,ny/2,1)
   PRINT *, 'Diagnostic from surface at domain center:'
   PRINT '("P:",F8.1," T:",F7.1," MR:",F6.4," RH:",F5.1," TD:",F7.1)', &
     psfc(nx/2,ny/2),tsfc(nx/2,ny/2),mrsfc(nx/2,ny/2),rhsfc(nx/2,ny/2),tdsfc(nx/2,ny/2) 
   DEALLOCATE(mrsfc)
 
   ! Get a few other miscellaneous variables that were added for diagnostics

   IF (mtype.EQ.'mm5') THEN
     CALL get_mm5_2d(current_lun, 'LWOUT    ', time_to_proc, lwout, &
                    'D   ', status)
     IF (status .NE. 0) lwout(:,:) = 1.e37
   
     CALL get_mm5_2d(current_lun, 'SWOUT    ', time_to_proc, swout, &
                    'D   ', status)
     IF (status .NE. 0) swout(:,:) = 1.e37

     CALL get_mm5_2d(current_lun, 'LWDOWN   ', time_to_proc, lwdown, &
                    'D   ', status)
     IF (status .NE. 0) lwdown(:,:) = 1.e37

     CALL get_mm5_2d(current_lun, 'SWDOWN   ', time_to_proc, swdown, &
                    'D   ', status)
     IF (status .NE. 0) swdown(:,:) = 1.e37

     CALL get_mm5_2d(current_lun, 'SHFLUX   ', time_to_proc, shflux, &
                    'D   ', status)
     IF (status .NE. 0) shflux(:,:) = 1.e37

     CALL get_mm5_2d(current_lun, 'LHFLUX   ', time_to_proc, lhflux, &
                    'D   ', status)
     IF (status .NE. 0) lhflux(:,:) = 1.e37

     IF (use_model_pbl) THEN
       PRINT *, 'Trying to obtain MM5 PBLHGT'
       CALL get_mm5_2d(current_lun, 'PBL HGT  ', time_to_proc, pblhgt, &
                    'D   ', status)
       IF ((status .NE. 0) .or. (maxval(pblhgt) .LE. 0))THEN
         PRINT *, '  MM5 PBLHGT not available, generating PBLHGT with LAPS algorithm'
         CALL model_pblhgt(thetasig,thetasfc,psig,zsig,terdot,nx,ny,ksigh,pblhgt)
         made_pbl = .true.
       ELSE 
         PRINT *, ' PBLHGT found and used'
         made_pbl = .false.
       ENDIF
     ELSE
       PRINT *, 'Using internally generated PBLHGT based on lfmpost.nl settings'
       CALL model_pblhgt(thetasig,thetasfc,psig,zsig,terdot,nx,ny,ksigh,pblhgt)
       made_pbl = .true.
     ENDIF
     IF (minval(pblhgt) .LE. 0.) THEN
       print *, 'Correcting PBL returned'
       do j=1,ny
         do i=1,nx
           if (pblhgt(i,j) .le. 0) then
             print *, 'PBL < 0 at i/j/val',i,j,pblhgt(i,j)
             pblhgt(i,j) = zsig(i,j,1)
           endif
         enddo
       enddo
     ENDIF
     print *, 'Min/Max PBL Height: ', minval(pblhgt),maxval(pblhgt)
     print *, 'Min/Max zsig(:,:,1): ',minval(zsig(:,:,1)),maxval(zsig(:,:,1)) 
     CALL get_mm5_2d(current_lun, 'GROUND T ', time_to_proc, ground_t, &
                    'D   ', status)
     IF (status .NE. 0) ground_t(:,:) = 1.e37
   
   ELSEIF(mtype(1:3).EQ.'wrf') THEN
     lwout(:,:) = 1.e37
     swout(:,:) = 1.e37
     lwdown(:,:) = 1.e37
     swdown(:,:) = 1.e37
     albedo(:,:) = 1.e37
     CALL get_wrfnc_2d(current_lun,'HFX','A',nx,ny,1,shflux,status)
     CALL get_wrfnc_2d(current_lun,'QFX','A',nx,ny,1,lhflux,status)
     CALL get_wrfnc_2d(current_lun,'GSW','A',nx,ny,1,swdown,status)

     ! Unlike mm5, wrf downward shortwave is net...i.e., albedo
     ! effect has been removed.  Lets get the WRF albedo (required
     ! making albedo an output variable in the WRF registry) and
     ! divide it back out.
     CALL get_wrfnc_2d(current_lun,'ALBEDO','A',nx,ny,1,albedo,status)
     IF (status .EQ. 0) THEN
       ! make sure destaggering did not create values outside the range
       WHERE (albedo .LT. 0.) albedo = 0.
       WHERE (albedo .GT. 1.0) albedo = 1.0
       PRINT '(A)', '  SWDOWN being changed from net to total'
       PRINT '(A,2F10.1)', 'Min/Max net SWDOWN:',MINVAL(swdown),MAXVAL(swdown)
       swdown = swdown / (1. - albedo)
       PRINT '(A,2F10.1)', 'Min/Max tot SWDOWN:',MINVAL(swdown),MAXVAL(swdown)
     ELSE
       PRINT '(A)', 'SWDOWN is still net value...ALBEDO not obtained!'
     ENDIF

     CALL get_wrfnc_2d(current_lun,'GLW','A',nx,ny,1,lwdown,status) 
     CALL model_pblhgt(thetasig,thetasfc,psig,zsig,terdot,nx,ny,ksigh,pblhgt)
     made_pbl = .true.
     CALL get_wrfnc_2d(current_lun,'TSK','A',nx,ny,1,ground_t,status)
   ENDIF

!  Commented out smoothing of these fields. 01/21/2004 BLS
!  IF (do_smoothing) THEN
!    DO smth = 0.5, -0.5, -1
!       CALL SMOOTH(lwout,nx,ny,1,smth)
!       CALL SMOOTH(swout,nx,ny,1,smth)
!       CALL SMOOTH(lwdown,nx,ny,1,smth)
!       CALL SMOOTH(swdown,nx,ny,1,smth)
!       CALL SMOOTH(shflux,nx,ny,1,smth)
!       CALL SMOOTH(lhflux,nx,ny,1,smth)
!       IF (.NOT. made_pbl) CALL SMOOTH(pblhgt,nx,ny,1,smth)
!       CALL SMOOTH(ground_t,nx,ny,1,smth)
!    ENDDO
!  ENDIF

  END SUBROUTINE get_sigma_vars   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  SUBROUTINE interp_thermo_3d

    ! Subroutine to compute two 3-D arrays of weighting coefficients
    ! and an integer array of trapping indices that are used for 
    ! the vertical interpolation from sigma to pressure 

    IMPLICIT NONE

    REAL                    :: deltalnp,deltap
    REAL                    :: dTvdlnPs,dqdp
    INTEGER                 :: status
    INTEGER                 :: ks, kp
    REAL                    :: p_lower_mb, p_upper_mb, p_mid_mb
    LOGICAL                 :: found_trap
    REAL                    :: tvbar,tvtop,tvbot,dz
    REAL, EXTERNAL          :: mixsat
    REAL                    :: weight_bot
 
    ! Interpolate some of the basic variables to pressure levels.  While
    ! doing this lets save the trapping indices and the weight assigned
    ! to the top trapping value for use in later interpolations
 
    trap_bot_ind = -1
    trap_top_ind = -1
    weight_top_lin = 0.0
    weight_top_log = 0.0
    pressure_loop: DO kp = 1, kprs
      ns_loop: DO j = 1, ny
        ew_loop: DO i = 1, nx
          
          ! Case 1: Is the pressure level below ground?
          IF (prslvl(kp) .GT. psig(i,j,1)) THEN
            trap_bot_ind(i,j,kp) = 0
            trap_top_ind(i,j,kp) = 1
            weight_top_lin(i,j,kp) = 1.0
            weight_top_log(i,j,kp) = 1.0
            
            ! Estimate the height of this level that is below
            ! ground.  First, estimate mean virtual temperature by
            ! using the MM5 base lapse rate to extrapolate the lowest
            ! virtual temp (above the boundary layer) downward. 

            deltalnp = ALOG(prslvl(kp))-ALOG(psig(i,j,1))
            tvbot = tvsig(i,j,1) + deltalnp*dTdlnPBase
            tvbar = (tvsig(i,j,1)+tvbot)*0.5
            tprs(i,j,kp) = tvbot/(1.+0.61*mrsig(i,j,1))
            
            ! Derive relative humidity assuming constant mixing ratio
            mrprs(i,j,kp) = mrsig(i,j,1)
            rhprs(i,j,kp) = MIN(relhum(tprs(i,j,kp),mrprs(i,j,kp), &
                            prslvl(kp))*100.,100.)
            thetaprs(i,j,kp) = potential_temp(tprs(i,j,kp),prslvl(kp))
            IF (ABS(prslvl(kp)-psig(i,j,1) ) .LT. 0.1) THEN
               zprs(i,j,kp) = zsig(i,j,1)
               dz = 0.0
            ELSE
              dz = tvbar * rog * ALOG(prslvl(kp)/psig(i,j,1))
              zprs(i,j,kp) = zsig(i,j,1) - dz
            ENDIF
          
          ! Case 2: Is the pressure level above the model top?
          ELSE IF (prslvl(kp) .LT. psig(i,j,ksigh)) THEN
            trap_bot_ind(i,j,kp) = ksigh
            trap_top_ind(i,j,kp) = 0
            weight_top_lin(i,j,kp) = 0.0
            weight_top_log(i,j,kp) = 0.0
  
            ! Now, we simply assume we are high enough up that the
            ! atmosphere is isothermal (above the tropopause).  For
            ! mixing ratio, we will use half the value of either the
            ! model top sigma layer or the next lowest pressure level,
            ! whichever is physically higher in the atmosphere.  The
            ! idea is to reduce the moisture toward extremely dry as
            ! you approach the outer edges of the atmosphere.

            tprs(i,j,kp) = tsig(i,j,ksigh)  
            IF (prslvl(kp-1) .GT. psig(i,j,ksigh)) THEN
              mrprs(i,j,kp) = mrsig(i,j,ksigh) * 0.5
            ELSE
              mrprs(i,j,kp) = mrprs(i,j,kp-1) * 0.5
            ENDIF
            tvtop = tprs(i,j,kp)*(1.+0.61*mrprs(i,j,kp) )
            tvbar = (tvsig(i,j,ksigh)+tvtop)*0.5 

            ! Derive relative humidity from the new mixing ratio   
            rhprs(i,j,kp) = MAX(relhum(tprs(i,j,kp),mrprs(i,j,kp), &
                               prslvl(kp))*100.,1.)
            rhprs(i,j,kp) = MIN(rhprs(i,j,kp),100.)
            thetaprs(i,j,kp) = potential_temp(tprs(i,j,kp),prslvl(kp))
            IF (ABS(prslvl(kp)-psig(i,j,ksigh) ) .LT. 0.1) THEN
               zprs(i,j,kp) = zsig(i,j,ksigh)
               dz = 0.0
            ELSE
              dz = tvbar * rog * ALOG(psig(i,j,ksigh)/prslvl(kp))
              zprs(i,j,kp) = zsig(i,j,ksigh) + dz
            ENDIF
        
          ! Case 3: We can trap this level between 2 valid sigma levels
          ELSE 
            sigma_loop: DO ks = 1, ksigh - 1
              IF ( (prslvl(kp) .LE. psig(i,j,ks) ) .AND. &
                   (prslvl(kp) .GE. psig(i,j,ks+1)) ) THEN

                 p_lower_mb = psig(i,j,ks)/100.
                 p_upper_mb = psig(i,j,ks+1)/100.
                 p_mid_mb   = prslvl(kp)/100.
                 CALL compute_lin_weights(p_mid_mb,p_lower_mb,p_upper_mb,&
                                          weight_bot, weight_top_lin(i,j,kp))
                 CALL compute_log_weights(p_mid_mb,p_lower_mb,p_upper_mb, &
                                          weight_bot,weight_top_log(i,j,kp)) 
                 trap_bot_ind(i,j,kp) = ks
                 trap_top_ind(i,j,kp) = ks + 1
                 tprs(i,j,kp) = weight_top_log(i,j,kp)*tsig(i,j,ks+1)+&
                            (1.0-weight_top_log(i,j,kp))*tsig(i,j,ks)
                 rhprs(i,j,kp) = weight_top_lin(i,j,kp)*rhsig(i,j,ks+1) + &
                           (1.-weight_top_lin(i,j,kp))*rhsig(i,j,ks)
                 zprs(i,j,kp) =weight_top_log(i,j,kp)*zsig(i,j,ks+1)+&
                           (1.-weight_top_log(i,j,kp))*zsig(i,j,ks)
                 
                 ! Diagnose theta and mixing ratio
                 thetaprs(i,j,kp) =  potential_temp(tprs(i,j,kp),prslvl(kp))
                 mrprs(i,j,kp) = mixsat(tprs(i,j,kp),prslvl(kp)) * &
                                 rhprs(i,j,kp)*0.01

                 EXIT sigma_loop 
               END IF
             ENDDO sigma_loop
           ENDIF
         ENDDO ew_loop
       ENDDO ns_loop
     ENDDO pressure_loop
     
     ! Compute 1000-500mb thickness if make_v5d is set.

     IF (make_v5d(domain_num)) thick_10_5 = zprs(:,:,k500)-zprs(:,:,k1000)

     ! Convert THETA into temperature
     !print *, 'Converting interpolated theta to temp..'
     !DO k = 1,kprs
     !  tprs(:,:,k) = thetaprs(:,:,k) / (p0/prslvl(k))**kappa
     !ENDDO
     print *, 'Computing dewpoint on pressure...'
     tdprs = tprs/((-rvolv * ALOG(rhprs*.01)*tprs) + 1.0) 
     
     ! Compute specific humidity and Tv pressure
     print *, 'Computing specific humidity on pressure...'
     DO k=1,kprs
       DO j=1,ny
         DO i=1,nx
           shprs(i,j,k) = mrprs(i,j,k) / ( 1. + mrprs(i,j,k) )
           tvprs(i,j,k) = tprs(i,j,k) * (1. + 0.61*mrprs(i,j,k))
         ENDDO
       ENDDO
       print *, 'Min/Max SH at ',prslvl(k),' = ',&
           minval(shprs(:,:,k)),maxval(shprs(:,:,k))
     ENDDO  

     ! Get reduced pressure for LAPS usage
     print *, 'Reducing pressure to ',redp_lvl, ' meters'
     CALL interp_press_to_z(prslvl, zprs, redp_lvl, redp,nx,ny,kprs)


     ! Use same routine to interpolate sea-level pressure.  This method
     ! is used in lieu of original reduction routine, because it will
     ! keep the MSL field consistent with the height field, which has been
     ! properly reduced.  It also produces a bit smoother field over the mountains.

     print *, 'Reducing pressure to MSL...'
     CALL interp_press_to_z(prslvl, zprs, 0., pmsl, nx,ny,kprs)

     PRINT *, 'Diagnostic print of interpolated values at center'
     PRINT '(A)', &
       '--------------------------------------------------------------------'
     PRINT '(A)', &
       'LEVEL  PRESSURE(Pa)  HEIGHT  THETA   RH  TEMP   DEWPT  MR      Tv'
     PRINT '(A)', &
       '--------------------------------------------------------------------'
     DO k = 1, kprs
      PRINT '(I5,2x,F12.1,2x,F6.0,2x,F6.2,2x,F4.0,2x,F5.1,2x,F5.1,2x,F8.6,2x,F5.1)', &
               k, prslvl(k),zprs(nx/2,ny/2,k),  &
               thetaprs(nx/2,ny/2,k), &
               rhprs(nx/2,ny/2,k),tprs(nx/2,ny/2,k),tdprs(nx/2,ny/2,k), &
               mrprs(nx/2,ny/2,k),tvprs(nx/2,ny/2,k)

     ENDDO


     PRINT *, '      Min/Max Sea Level Press (Pa): ' ,minval(pmsl), maxval(pmsl)
     PRINT *, '      Min/Max Reduced Press (Pa):     ' ,minval(redp),&
                                                      maxval(redp)     
     RETURN
     
  END SUBROUTINE interp_thermo_3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE interp_winds
    
    ! Interpolates the momentum variables to pressure and surface

    IMPLICIT NONE

    REAL, ALLOCATABLE       :: usig ( : , : , : )
    REAL, ALLOCATABLE       :: vsig ( : , : , : )
    REAL, ALLOCATABLE       :: wsig ( : , : , : )
    REAL, ALLOCATABLE       :: wsigf( : , : , : )
    REAL, ALLOCATABLE       :: below_ground ( : , : )
    REAL, ALLOCATABLE       :: tkesig(:,:,:)
    REAL, ALLOCATABLE       :: tkesigcomp(:,:,:)
    INTEGER                 :: status
    REAL                    :: recipdx, dudy, dvdx

    ALLOCATE (below_ground (nx,ny))
    ! Get u-component of wind
    ALLOCATE (usig ( nx , ny , ksigh ) )
    IF (mtype.EQ.'mm5') THEN
      CALL get_mm5_3d(current_lun, 'U        ', time_to_proc, usig, &
                      'D    ', status)
    ELSEIF(mtype(1:3).EQ.'wrf')THEN
      ! Get/destagger WRF u wind
     CALL get_wrfnc_3d(current_lun,'U','A',nx,ny,ksigh,1,usig,status)

    ENDIF

    IF (status.NE.0) THEN
      PRINT '(A)', 'Problem getting u-wind!  Aborting...'
      CALL abort
    ENDIF
    IF (do_smoothing) THEN
      DO smth = 0.5,-0.5,-1
        CALL SMOOTH(usig,nx,ny,ksigh,smth)
      ENDDO
    ENDIF
    ! See if U10 is available (should be if using MRF PBL
    ! under MM5v3 R4)
    usfc(:,:) = 0.
    IF (mtype.EQ.'mm5')THEN
      CALL get_mm5_2d(current_lun, 'U10      ', time_to_proc, usfc, &
                      'D   ', status)
    ELSEIF(mtype(1:3).EQ.'wrf')THEN
      ! Get WRF U10 
      CALL get_wrfnc_2d(current_lun,'U10','A',nx,ny,1,usfc,status)
    ENDIF
    IF ((status .ne. 0.).OR. &
        ((maxval(usfc).eq.0.).AND.(minval(usfc).eq.0))) THEN
      PRINT *, 'U10 not available, using lowest sigma layer for usfc'
      usfc(:,:) = usig(:,:,1)
    ELSE
      IF (do_smoothing) THEN
        DO smth = 0.5,-0.5,-1
          CALL SMOOTH(usfc,nx,ny,1,smth)
        ENDDO
      ENDIF
    ENDIF
    below_ground = usfc
    CALL vinterp_3d(usig, trap_bot_ind, trap_top_ind, &
                     weight_top_lin, below_ground, uprs, &
                     nx, ny, ksigh, kprs)
    
    ! Get v-component of wind
    ALLOCATE (vsig ( nx , ny , ksigh ) )

    IF (mtype .EQ. 'mm5') THEN
      CALL get_mm5_3d(current_lun, 'V        ', time_to_proc, vsig, &
                      'D    ', status)
    ELSEIF(mtype(1:3) .EQ. 'wrf') THEN

      ! Get WRF V
      CALL get_wrfnc_3d(current_lun,'V','A',nx,ny,ksigh,1,vsig,status)
 
    ENDIF
    IF (status.NE.0) THEN
      PRINT '(A)', 'Problem getting v-wind!  Aborting...'
      CALL abort
    ENDIF
    IF (do_smoothing) THEN
      DO smth = 0.5,-0.5,-1
        CALL SMOOTH(vsig,nx,ny,ksigh,smth)
      ENDDO
    ENDIF
    
    ! Try to use V10 for surface V if available
    vsfc(:,:) = 0.
    IF (mtype.EQ.'mm5')THEN
      CALL get_mm5_2d(current_lun, 'V10      ', time_to_proc, vsfc, &
                      'D   ', status)
    ELSEIF(mtype(1:3) .EQ. 'wrf') THEN
      ! Get WRF V10
      CALL get_wrfnc_2d(current_lun,'V10','A',nx,ny,1,vsfc,status)
    ENDIF
    IF ((status .ne. 0.).OR. &
        ((maxval(vsfc).eq.0.).AND.(minval(vsfc).eq.0))) THEN
      PRINT *, 'V10 not available, using lowest sigma layer for vsfc'
      vsfc(:,:) = vsig(:,:,1)
    ELSE 
      IF (do_smoothing) THEN
        DO smth = 0.5,-0.5,-1
          CALL SMOOTH(vsfc,nx,ny,1,smth)
        ENDDO
      ENDIF
    ENDIF
    below_ground = vsfc
    CALL vinterp_3d(vsig, trap_bot_ind, trap_top_ind, &
                    weight_top_lin, below_ground,vprs ,&
                    nx, ny, ksigh, kprs)

    ! Get w-component of wind
    ALLOCATE (wsig ( nx , ny , ksigh ) )
    ALLOCATE (wsigf( nx , ny , ksigf ) )

    IF (mtype .EQ. 'mm5') THEN
      CALL get_mm5_3d(current_lun, 'W        ', time_to_proc, wsigf, &
                      'D    ', status)
    ELSEIF (mtype(1:3) .EQ. 'wrf') THEN
      ! Get WRF W on full layers
      CALL get_wrfnc_3d(current_lun,'W','A',nx,ny,ksigf,1,wsigf,status)
    ENDIF
    IF (status.NE.0) THEN
      PRINT '(A)', 'Problem getting w-wind!  Aborting...'
      CALL abort
    ENDIF
    DO k = 1, ksigh
      wsig(:,:,k) = 0.5*(wsigf(:,:,k)+wsigf(:,:,k+1))
    ENDDO
    DEALLOCATE (wsigf)
    wsfc(:,:) = wsig(:,:,1)
    IF (do_smoothing) THEN
      DO smth = 0.5,-0.5,-1
        CALL smooth(wsig,nx,ny,ksigh,smth)
        CALL smooth(wsfc,nx,ny,1,smth)
      ENDDO
    ENDIF
    below_ground(:,:) = 0.0
    CALL vinterp_3d(wsig,trap_bot_ind, trap_top_ind, &
                    weight_top_lin, below_ground, wprs, &
                    nx, ny, ksigh, kprs)
    DEALLOCATE (wsig)

    ! Compute omega for LAPS
    DO k = 1,kprs
      omprs(:,:,k)=-(prslvl(k)*wprs(:,:,k)*grav)/(r*tvprs(:,:,k))
    ENDDO
   
    ! Compute the storm relative helicity

    CALL helicity(usig, vsig, zsig, terdot, nx, ny, ksigh, srhel)
    PRINT *, 'Min/Max Helicity: ', MINVAL(srhel), MAXVAL(srhel)

    ! Compute ventilation index
    CALL ventilation(usig,vsig,zsig,pblhgt,terdot,nx,ny,ksigh,upbl,vpbl,vnt_index)
    PRINT *, 'Min/Max Ventilation : ', MINVAL(vnt_index), MAXVAL(vnt_index)    
 
    ! Compute TKE
    ALLOCATE(tkesigcomp(nx,ny,ksigh))  
    CALL compute_tke(psig,tsig,usig,vsig,zsig,terdot,nx,ny,ksigh,tkesigcomp)
    PRINT '(A,2F10.3)','Min/Max TKE computed from DTF3: ', &
       MINVAL(tkesigcomp),MAXVAL(tkesigcomp)
    DEALLOCATE (usig)
    DEALLOCATE (vsig)
    ! Get TKE
    PRINT *, 'Getting TKE....'
    ALLOCATE(tkesig(nx,ny,ksigh))
    below_ground = 0.0
    IF (mtype .EQ. 'mm5') THEN
      CALL get_mm5_3d(current_lun, 'TKE      ', time_to_proc, tkesig, &
                      'D    ', status)
    ELSEIF (mtype(1:3) .EQ. 'wrf') THEN
      CALL get_wrfnc_3d(current_lun,'TKE_MYJ','A',nx,ny,ksigf,1,tkesig,status)
    ENDIF
    IF (status .NE. 0) THEN
       tkesig = tkesigcomp
    ELSE
      PRINT '(A,2F10.3)', 'Min/Max TKE from model output: ', &
         MINVAL(tkesig),MAXVAL(tkesig)
    ENDIF
    DEALLOCATE(tkesigcomp)
    PRINT *, 'Vertically interpolating TKE'
    CALL vinterp_3d(tkesig,trap_bot_ind,trap_top_ind,&
                    weight_top_lin, below_ground, tkeprs, &
                    nx,ny,ksigh,kprs)
    DEALLOCATE(tkesig)
    DEALLOCATE(below_ground)
    
    PRINT *, '      Min/Max U wind (m/s):   ', minval(uprs),maxval(uprs)
    PRINT *, '      Min/Max V wind (m/s):   ', minval(vprs),maxval(vprs)
    PRINT *, '      Min/Max W wind (m/s):   ', minval(wprs),maxval(wprs)
    PRINT *, '      Min/Max Omega (Pa/s):   ', minval(omprs),maxval(omprs)
    PRINT *, '      Min/Max SR Helicity:    ', minval(srhel),maxval(srhel)
    PRINT *, '      Min/Max TKE:            ', minval(tkeprs),maxval(tkeprs)
    ! If making Vis5D output, compute vorticity
    IF (make_v5d(domain_num)) THEN
       recipdx = 1.0 / (2.0 * grid_spacing)
       DO k = 1,kprs
         DO j = 2,ny-1
           DO i = 2,nx-1
             dvdx = (vprs(i+1,j,k) - vprs(i-1,j,k)) * recipdx
             dudy = (uprs(i,j+1,k) - uprs(i,j-1,k)) * recipdx
             abs_vort(i,j,k) = (mapfac_d(i,j) * (dvdx-dudy) ) + coriolis(i,j)
           ENDDO
         ENDDO
       ENDDO
    ENDIF
    RETURN
 
  END SUBROUTINE interp_winds 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_clouds_reflectivity
    ! Gets the cloud parameters on sigma surfaces, computes various
    ! parameters, then interpolates to sigma
 
    IMPLICIT NONE
    REAL, ALLOCATABLE             :: below_ground ( : , : )
    REAL, ALLOCATABLE             :: cldliqmr_sig ( : , : , : )
    REAL, ALLOCATABLE             :: cldicemr_sig ( : , : , : )
    REAL, ALLOCATABLE             :: condmr_sig   ( : , : , : )
    REAL, ALLOCATABLE             :: rainmr_sig   ( : , : , : )
    REAL, ALLOCATABLE             :: snowmr_sig   ( : , : , : )
    REAL, ALLOCATABLE             :: graupelmr_sig( : , : , : )
    REAL, ALLOCATABLE             :: refl_sig     ( : , : , : )
    REAL                          :: cldliqthresh
    REAL                          :: cldicethresh
    REAL                          :: snowthresh
    INTEGER                       :: status
    REAL                          :: zero_thresh

    ! Set up the zero_thresh (threshold below which arrays
    ! are zeroed out for microphysics species) for MM5
    ! and WRF separately

    IF (mtype .EQ. 'mm5') THEN
       zero_thresh = 1.e-10
    ELSE
       zero_thresh = 1.e-6
    ENDIF
    ! Get each of the cloud species mixing ratio arrays
    ALLOCATE (cldliqmr_sig (nx,ny,ksigh))
    ALLOCATE (cldicemr_sig (nx,ny,ksigh))
    ALLOCATE (rainmr_sig (nx,ny,ksigh))
    ALLOCATE (snowmr_sig (nx,ny,ksigh))
    ALLOCATE (graupelmr_sig(nx,ny,ksigh))

    ! Cloud liquid and rain water
    IF (clwflag) THEN
      IF (mtype(1:3) .EQ. 'mm5') THEN
        CALL get_mm5_3d(current_lun, 'CLW      ', time_to_proc, cldliqmr_sig,&
                      'D    ', status)
      ELSEIF(mtype(1:3).EQ.'wrf') THEN
        CALL get_wrfnc_3d(current_lun,'QCLOUD','A',nx,ny,ksigh,1,cldliqmr_sig,&
                          status)
      ENDIF
      IF (status.NE.0) THEN
        PRINT '(A)', 'Problem getting CLW!  Setting to 0.'
        cldliqmr_sig = 0.0
      ENDIF
      WHERE(cldliqmr_sig .LT. zero_thresh) cldliqmr_sig = 0.0
      clwmrsfc(:,:) = cldliqmr_sig(:,:,1)

      IF (mtype .EQ. 'mm5') THEN 
        CALL get_mm5_3d(current_lun, 'RNW      ', time_to_proc, rainmr_sig,&
                      'D    ', status)
      ELSEIF (mtype(1:3) .EQ. 'wrf') THEN
         CALL get_wrfnc_3d(current_lun,'QRAIN','A',nx,ny,ksigh,1,rainmr_sig,&
                          status)
      ENDIF
      IF (status.NE.0) THEN
        PRINT '(A)', 'Problem getting RNW!  Setting to 0.'
        rainmr_sig = 0.0
      ENDIF
      WHERE(rainmr_sig .LT. zero_thresh) rainmr_sig = 0.0
      rainmrsfc(:,:) = rainmr_sig(:,:,1)

    ELSE
      cldliqmr_sig = 0.0
      rainmr_sig = 0.0
      clwmrsfc = 0.0
      rainmrsfc = 0.0
    ENDIF
     
    ! Snow and Ice
    IF (iceflag) THEN
      IF (mtype .EQ. 'mm5')THEN
        CALL get_mm5_3d(current_lun, 'ICE      ', time_to_proc, cldicemr_sig,&
                      'D    ', status)
      ELSEIF(mtype(1:3).eq.'wrf') THEN
        CALL get_wrfnc_3d(current_lun,'QICE','A',nx,ny,ksigh,1,cldicemr_sig,&
                          status)
      ENDIF
      IF (status.NE.0) THEN
        PRINT '(A)', 'Problem getting ICE!  Setting to 0.'
        cldicemr_sig = 0.0
      ENDIF
      WHERE(cldicemr_sig .LT. zero_thresh) cldicemr_sig = 0.0
      icemrsfc(:,:) = cldicemr_sig(:,:,1)

      IF (mtype .EQ. 'mm5') THEN  
        CALL get_mm5_3d(current_lun, 'SNOW     ', time_to_proc, snowmr_sig,&
                      'D    ', status)
      ELSEIF (mtype(1:3) .EQ. 'wrf') THEN
         CALL get_wrfnc_3d(current_lun,'QSNOW','A',nx,ny,ksigh,1,snowmr_sig,&
                          status)
      ENDIF

      IF (status.NE.0) THEN
        PRINT '(A)', 'Problem getting SNOW!  Setting to 0.'
        snowmr_sig = 0.0
      ENDIF
      WHERE(snowmr_sig .LT. zero_thresh) snowmr_sig = 0.0
      snowmrsfc(:,:) = snowmr_sig(:,:,1)
    ELSE
      cldicemr_sig = 0.0
      snowmr_sig = 0.0
      icemrsfc = 0.0
      snowmrsfc = 0.0
    ENDIF  

    ! Graupel
    IF (graupelflag) THEN
      IF (mtype.EQ.'mm5') THEN
        CALL get_mm5_3d(current_lun, 'GRAUPEL  ', time_to_proc, graupelmr_sig,&
                      'D    ', status)
      ELSEIF(mtype(1:3) .EQ. 'wrf') THEN
        CALL  get_wrfnc_3d(current_lun,'QGRAUP','A',nx,ny,ksigh,1, &
                          graupelmr_sig,status)
      ENDIF
      IF (status.NE.0) THEN
        PRINT '(A)', 'Problem getting GRAUPEL!  Setting to 0.'
        graupelmr_sig = 0.0
      ENDIF
      WHERE(graupelmr_sig .LT. zero_thresh) graupelmr_sig = 0.0
      graupmrsfc(:,:) = graupelmr_sig(:,:,1)
    ELSE
      graupelmr_sig = 0.0
      graupmrsfc = 0.0
    ENDIF 

    ! Now compute the bases, tops, coverage, and ceiling
     
    CALL clouds(nx,ny,ksigh,grid_spacing,cldliqmr_sig,cldicemr_sig, &
                snowmr_sig, zsig, terdot, cldbase, cldtop, ceiling, cldamt)


    ! Compute total condensate so we can compute integrated liquid water
    ALLOCATE (condmr_sig(nx,ny,ksigh))
    condmr_sig = cldliqmr_sig + cldicemr_sig + snowmr_sig + rainmr_sig &
                 + graupelmr_sig
    print *, 'Calling integrated_liquid'
    print *, 'min/max condmr_sig = ',minval(condmr_sig),maxval(condmr_sig)
    print *, 'min/max mrsig = ',minval(mrsig),maxval(mrsig)
    print *, 'min/max rhodrysig = ', minval(rhodrysig),maxval(rhodrysig)
    print *, 'min/max zsig = ', minval(zsig),maxval(zsig)
    CALL integrated_liquid(nx,ny,ksigh,condmr_sig,mrsig, rhodrysig,zsig,terdot, &
                           intliqwater,totpcpwater) 
    DEALLOCATE(condmr_sig)   
 
    ! Compute reflectivities
    ALLOCATE (refl_sig(nx,ny,ksigh))
    CALL reflectivity(nx,ny,ksigh,rhomoistsig, zsig, &
                      rainmr_sig,cldicemr_sig,snowmr_sig,graupelmr_sig,&
                      refl_sig, max_refl, echo_tops)

    ! Interpolate our 3D cloud mixing ratios and reflectivities to
    ! the pressure surfaces
    ALLOCATE (below_ground (nx,ny ))
    
    below_ground(:,:)=0.0
    CALL vinterp_3d(cldliqmr_sig, trap_bot_ind, trap_top_ind, &
                    weight_top_lin, below_ground, cldliqmr_prs, &
                    nx, ny, ksigh, kprs)
 
    CALL vinterp_3d(cldicemr_sig, trap_bot_ind, trap_top_ind, &
                    weight_top_lin, below_ground, cldicemr_prs, &
                    nx, ny, ksigh, kprs)

    CALL vinterp_3d(rainmr_sig, trap_bot_ind, trap_top_ind, &
                    weight_top_lin, below_ground, rainmr_prs, &
                    nx, ny, ksigh, kprs)

    CALL vinterp_3d(snowmr_sig,trap_bot_ind, trap_top_ind, &
                    weight_top_lin, below_ground,snowmr_prs, &
                    nx, ny, ksigh, kprs)
 
    CALL vinterp_3d(graupelmr_sig, trap_bot_ind, trap_top_ind, &
                    weight_top_lin, below_ground, graupelmr_prs, &
                    nx, ny, ksigh, kprs)
 
    CALL vinterp_3d(refl_sig, trap_bot_ind, trap_top_ind, &
                    weight_top_lin, below_ground,refl_prs, &
                    nx, ny, ksigh, kprs)
 
    DEALLOCATE (below_ground) 
    
    ! Diagnostic Print of top 4 layers

    PRINT *, 'Mean Condensate Values:'
    PRINT *, 'Level   Species   Min     Max    Mean'
    PRINT *, '(mb)             (g/kg)  (g/kg) (g/kg)'
    PRINT *, '------ -------- ------- ------- -------'
    DO k = kprs - 3, kprs
      print '(f6.1,1x,a8,1x,f7.4,1x,f7.4,1x,f7.4)', &
        prslvl(k)*0.01, 'ICE     ', MINVAL(cldicemr_prs(:,:,k))*1000., &
        MAXVAL(cldicemr_prs(:,:,k))*1000., &
        SUM( cldicemr_prs(:,:,k) )*1000./(nx*ny)
      print '(f6.1,1x,a8,1x,f7.4,1x,f7.4,1x,f7.4)', &
        prslvl(k)*0.01, 'SNOW    ', MINVAL(snowmr_prs(:,:,k))*1000., &
        MAXVAL(snowmr_prs(:,:,k))*1000., &
        SUM( snowmr_prs(:,:,k) )*1000./(nx*ny)
      print '(f6.1,1x,a8,1x,f7.4,1x,f7.4,1x,f7.4)', &
        prslvl(k)*0.01, 'CLW     ', MINVAL(cldliqmr_prs(:,:,k))*1000., &
        MAXVAL(cldliqmr_prs(:,:,k))*1000., &
        SUM( cldliqmr_prs(:,:,k) )*1000./(nx*ny)
      print '(f6.1,1x,a8,1x,f7.4,1x,f7.4,1x,f7.4)', &
        prslvl(k)*0.01, 'RAIN    ', MINVAL(rainmr_prs(:,:,k))*1000., &
        MAXVAL(rainmr_prs(:,:,k))*1000., &
        SUM( rainmr_prs(:,:,k) )*1000./(nx*ny)
      print '(f6.1,1x,a8,1x,f7.4,1x,f7.4,1x,f7.4)', &
        prslvl(k)*0.01, 'GRAUPEL ', MINVAL(graupelmr_prs(:,:,k))*1000., &
        MAXVAL(graupelmr_prs(:,:,k))*1000., &
        SUM( graupelmr_prs(:,:,k) )*1000./(nx*ny)
    ENDDO

    ! Make  precip codes
    pcptype_sfc = nonecode
    pcptype_prs = nonecode
    DO j = 1 , ny
      DO i = 1 , nx
        DO k = 1 , kprs
           IF (rainmr_prs(i,j,k).GT.0) THEN
             IF (snowmr_prs(i,j,k) .GT. 0.) THEN
               IF (graupelmr_prs(i,j,k).GT.snowmr_prs(i,j,k)) THEN
                 pcptype_prs(i,j,k)=rainicecode
               ELSE
                 pcptype_prs(i,j,k) = rainsnowcode
               ENDIF
             ELSE
               IF (graupelmr_prs(i,j,k) .GT. 0) THEN
                 pcptype_prs(i,j,k) = rainicecode
               ELSE
                 pcptype_prs(i,j,k) = raincode
               ENDIF
             ENDIF
           ELSE
             IF (snowmr_prs(i,j,k).GT.0) THEN
               IF (graupelmr_prs(i,j,k) .GT. snowmr_prs(i,j,k)) THEN
                 pcptype_prs(i,j,k) = sleetcode
               ELSE
                 pcptype_prs(i,j,k) = snowcode
               ENDIF
             ELSE
               IF (graupelmr_prs(i,j,k).GT. 0) THEN
                 pcptype_prs(i,j,k)= sleetcode
               ELSE
                 pcptype_prs(i,j,k) = nonecode
               ENDIF
             ENDIF
           ENDIF
         ENDDO
         IF (rainmr_sig(i,j,1).GT.0) THEN
           IF (snowmr_sig(i,j,1) .GT. 0.) THEN
             IF (graupelmr_sig(i,j,1).GT.snowmr_sig(i,j,1)) THEN
                 pcptype_sfc(i,j)=rainicecode
             ELSE
                 pcptype_sfc(i,j) = rainsnowcode
             ENDIF
           ELSE
             IF (graupelmr_sig(i,j,1) .GT. 0) THEN
               pcptype_sfc(i,j) = rainicecode
             ELSE
               pcptype_sfc(i,j) = raincode
             ENDIF
           ENDIF
         ELSE
           IF (snowmr_sig(i,j,1).GT.0) THEN
             IF (graupelmr_sig(i,j,1) .GT. snowmr_sig(i,j,1)) THEN
               pcptype_sfc(i,j) = sleetcode
             ELSE
               pcptype_sfc(i,j) = snowcode
             ENDIF
           ELSE
             IF (graupelmr_sig(i,j,1).GT. 0) THEN
               pcptype_sfc(i,j)= sleetcode
             ELSE
               pcptype_sfc(i,j) = nonecode
             ENDIF
           ENDIF
         ENDIF
       ENDDO
     ENDDO
     ! Set surface reflectivity to be the maximum value in the lowest
     ! third of the atmosphere.
     ! 
     ! refl_sfc(:,:) = MAXVAL(refl_sig(:,:,1:ksigh/3),DIM=3)
   
     ! No...go back to using only the lowest level for consistency with
     ! the precip type icons.

     refl_sfc(:,:) = refl_sig(:,:,1)

     DEALLOCATE(refl_sig)
     DEALLOCATE(cldliqmr_sig)
     DEALLOCATE(cldicemr_sig)
     DEALLOCATE(rainmr_sig)
     DEALLOCATE(snowmr_sig)
     DEALLOCATE(graupelmr_sig)
     PRINT *, '      Min/Max Sfc Reflectivity (dBZ): ',minval(refl_sfc), &
                                                       maxval(refl_sfc)
     PRINT *, '      Min/Max 3d Reflectivity  (dBZ): ',minval(refl_prs), &
                                                       maxval(refl_prs)
     PRINT *, '      Cloud Base (m) at domain center: ', cldbase(nx/2,ny/2)
     PRINT *, '      Cloud Top (m) at domain center:  ', cldtop(nx/2,ny/2)
     PRINT *, '      Cloud cover (fraction) at center:', cldamt(nx/2,ny/2)
     PRINT *, '      Ceiling (ALG m) at center:       ', ceiling(nx/2,ny/2)
     PRINT *, '      Precip Type codes at center:'
     PRINT *, '      LEVEL(mb)       CODE'
     PRINT'(6x,"    SFC  ",9x,F2.0)',pcptype_sfc(nx/2,ny/2)
     DO k = 1 , kprs
       PRINT '(6x,F9.0,9x,F2.0)',prslvl(k)*0.01,pcptype_prs(nx/2,ny/2,k)
     ENDDO
     PRINT *, '      Min/Max Int Liquid Water (kg/m2):',minval(intliqwater), &
                                                        maxval(intliqwater)
     PRINT *, '      Min/Max Tot Precip Water (kg/m2):',minval(totpcpwater), &
                                                        maxval(totpcpwater)
     
    RETURN

  END SUBROUTINE get_clouds_reflectivity   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_precip
    ! Gets the total precipitation fields from the run and subtracts off
    ! intial and current total values to produce incremental amounts as well
    ! as total amounts.  Also derives snowfall.

    IMPLICIT NONE

    REAL, ALLOCATABLE                 :: raincon ( : , : )
    REAL, ALLOCATABLE                 :: rainnon ( : , : )
    INTEGER                           :: status
    INTEGER, ALLOCATABLE              :: fallen_precip_type ( : , : )

    ALLOCATE ( raincon ( nx , ny ) )
    ALLOCATE ( rainnon ( nx , ny ) ) 
    ALLOCATE ( fallen_precip_type ( nx , ny ) )

    IF (mtype .EQ. 'mm5') THEN
      CALL get_mm5_2d(current_lun, 'RAIN CON ', time_to_proc, raincon, &
                      'D   ', status)
    ELSEIF(mtype(1:3) .EQ. 'wrf') THEN
      CALL get_wrfnc_2d(current_lun,'RAINC','A',nx,ny,1,raincon,status)
      ! Convert to cm
      raincon = raincon * 0.1
    ENDIF
    IF (status.NE.0) THEN
      PRINT '(A)', 'Problem RAIN CON'
      raincon = 0.0
    ENDIF
    WHERE(raincon .LT. .00001) raincon = 0.

    IF (mtype .EQ. 'mm5') THEN
      CALL get_mm5_2d(current_lun, 'RAIN NON ', time_to_proc, rainnon, &
                      'D   ', status)
    ELSEIF (mtype(1:3) .EQ. 'wrf') THEN
      CALL get_wrfnc_2d(current_lun,'RAINNC','A',nx,ny,1,rainnon,status)
      ! Convert to cm
      rainnon = rainnon * 0.1
    ENDIF

    IF (status.NE.0) THEN
      PRINT '(A)', 'Problem RAIN NON'
      rainnon = 0.0
    ENDIF  
    WHERE(rainnon .LT. .00001) rainnon = 0.

    IF (make_v5d(domain_num)) THEN
      IF (mtype .EQ. 'mm5') THEN
         CALL get_mm5_2d(current_lun, 'SNOWCOVR ', time_to_proc, snowcover, &
                      'D   ', status)
      ELSEIF(mtype(1:3).EQ.'wrf') THEN
         CALL get_wrfnc_2d(current_lun,'ACSNOW','A',nx,ny,1,snowcover,status)
      ENDIF
      IF (status.NE.0) THEN
        PRINT '(A)', 'Problem SNOWCOVR'
        snowcover = 0.0
      ENDIF
    ENDIF

    
    ! Convert from cm to m
    rainnon = rainnon * 0.01
    raincon = raincon * 0.01
  
    IF (initialize) THEN
      pcp_init = raincon + rainnon
      pcp_inc = 0.0
      pcp_tot = 0.0
      con_pcp_init = raincon
      con_pcp_inc = 0.0
      con_pcp_tot = 0.0
      snow_inc = 0.0
      snow_tot  = 0.0
    ELSE
      pcp_inc = 0.0
      con_pcp_inc = 0.0
      snow_inc = 0.0
      con_pcp_inc = raincon - con_pcp_init - con_pcp_tot
      pcp_inc = raincon + rainnon - pcp_init - pcp_tot
      WHERE (pcp_inc .LT. 0) pcp_inc = 0.
      WHERE (con_pcp_inc .LT. 0) con_pcp_inc = 0.
      pcp_tot = pcp_tot + pcp_inc
      con_pcp_tot = con_pcp_tot + con_pcp_inc
      CALL wintprec (tsig, zsig, zprs, psfc, tsfc, terdot, pcp_inc, nx, ny, &
                     ksigh, kprs, k700, k850, k1000, raincon, &
                     fallen_precip_type)
      CALL snowfall (tsfc, pcp_inc, fallen_precip_type, nx, ny, &
                     snow_inc, snow_tot)
    ENDIF
    DEALLOCATE (raincon)
    DEALLOCATE (rainnon)
    DEALLOCATE (fallen_precip_type)
    PRINT *, '      Min/Max Inc. Precip (m): ',minval(pcp_inc),maxval(pcp_inc)
    PRINT *, '      Min/Max Conv. Precip   : ',minval(con_pcp_inc), &
                                               maxval(con_pcp_inc)
    PRINT *, '      Min/Max Acc. Precip (m): ',minval(pcp_tot),maxval(pcp_tot)
    PRINT *, '      Min/Max Conv. Precip   : ',minval(con_pcp_tot), &
                                               maxval(con_pcp_tot)
    PRINT *, '      Min/Max Inc. Snow   (m): ',minval(snow_inc),maxval(snow_inc)
    PRINT *, '      Min/Max Acc. Snow   (m): ',minval(snow_tot),maxval(snow_tot)
    
    RETURN
  END SUBROUTINE get_precip
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE make_stability

    ! Make severe weather stability indices.  Currently just creates CAPE, CIN,
    ! and lifted index

    IMPLICIT NONE
    REAL, ALLOCATABLE                     :: thetaesig ( : , : , : )
    REAL, ALLOCATABLE                     :: zsigf      ( : , : , : )
    INTEGER                               :: itest, jtest
    ALLOCATE(thetaesig ( nx , ny , ksigh ) )
    ALLOCATE(zsigf      ( nx , ny , ksigf ) )

    itest = nx/2
    jtest = ny/2
   
    ! Compute thetae on sigma
    print *, ' --- computing thetae on sigma --- ' 
    DO k = 1 , ksigh
      DO j = 1 , ny
        DO i = 1 , nx
           thetaesig(i,j,k) = eq_potential_temp(tsig(i,j,k),psig(i,j,k),&
                                                mrsig(i,j,k),rhsig(i,j,k)*0.01)
        ENDDO
      ENDDO
    ENDDO
    ! Compute heights on full sigma surfaces
    print *, ' --- computing heights on full sigmas ---'
    DO k = 2 , ksigh
      zsigf(:,:,k) = ( zsig(:,:,k-1) + zsig(:,:,k) ) * 0.5
    ENDDO
    zsigf(:,:,1) = terdot
    zsigf(:,:,ksigf) = 2 * zsig(:,:,ksigh) - zsigf(:,:,ksigh)
   
    ! Call routine to compute CAPE, CIN, and LI

    print *, ' --- computing CAPE/CIN/LI ---' 
    CALL capecin(psig*0.01,tsig,thetaesig,thetasig,rhsig*0.01, &
                 zsigf,tprs,liftedind,cape,cin,k500,nx,ny,ksigh,kprs)
    cin = -cin
    print *, 'Stability Stuff'
    print *, 'K   PRESSURE  T      THETA  THETAE RH   Z'
    DO k = 1, ksigh
      print '(I4,F10.1,3F7.1,F5.0,F8.0)', &
        k, psig(itest,jtest,k), tsig(itest,jtest,k),thetasig(itest,jtest,k), &
        thetaesig(itest,jtest,k),rhsig(itest,jtest,k),zsigf(itest,jtest,k)
    ENDDO
    print *, 'CAPE =',cape(itest,jtest)
    print *, 'CIN = ',cin(itest,jtest)
    print *, 'LI = ', liftedind(itest,jtest)
    DEALLOCATE (thetaesig)
    DEALLOCATE (zsigf)

    PRINT *, '      Min/Max CAPE (J/kg):  ', MINVAL(cape),MAXVAL(cape)
    PRINT *, '      Min/Max CIN  (J/kg):  ', MINVAL(cin),MAXVAL(cin)
    PRINT *, '      Min/Max LI   (K):     ', MINVAL(liftedind),MAXVAL(liftedind)
    RETURN
  END SUBROUTINE make_stability 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE make_misc
    ! Makes miscellaneous derived parameters
   
    IMPLICIT NONE
    REAL, EXTERNAL       :: heatindex
    REAL, ALLOCATABLE    :: tdsig(:,:,:)
    ! Compute surface visibility from relative humidity and temp/dewpoint

    visibility = 6000.0 * (tsfc - tdsfc) / ( rhsfc**1.75)
 
    ! Convert to meters from km
    visibility = visibility * 1000.
    WHERE(visibility .GT. 99990.) visibility = 99990.
    DO j = 1 , ny
      DO i = 1 , nx
        ! Compute heat index if temp is above 80F (300K)
        IF (tsfc(i,j) .GE.300.) THEN
          heatind(i,j) = HEATINDEX(tsfc(i,j),rhsfc(i,j))
        ELSE
          heatind(i,j) = tsfc(i,j)
        ENDIF
    
      ENDDO
    ENDDO
   
    ! Compute fire indices

    ! We need dewpoint on sigma for the Haines indices
    ALLOCATE (tdsig(nx,ny,ksigh))
    tdsig = tsig/((-rvolv*ALOG(rhsig*0.01)*tsig)+1.0)

    ! Mid-level Haines Index
    CALL haines_layer(psig*0.01,tsig,tdsig,ham_index,nx,ny,ksigh, &
                      850., 700.)

    ! High-level Haines Index
    CALL haines_layer(psig*0.01,tsig,tdsig,hah_index,nx,ny,ksigh, &
                      700., 500.)

    DEALLOCATE(tdsig)

    ! Fosberg FWI

    CALL fosberg_fwi(tsfc,rhsfc,psfc*0.01,usfc,vsfc,nx,ny,fwi_index)
    PRINT *, '      Min/Max Surface Visibility (m): ',MINVAL(visibility), &
                                                      MAXVAL(visibility)
    PRINT *, '      Min/Max Surface Heat Index (K): ',MINVAL(heatind), &
                                                      MAXVAL(heatind)
    RETURN
  END SUBROUTINE make_misc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE deallocate_internal

    ! This routine deallocates all of the arrays internal to this module that
    ! may be shared among more than one routine in this module.  It does not
    ! deallocate the zsig array, however, which gets reused with each 
    ! time step!!

    IMPLICIT NONE
    
    DEALLOCATE ( psig )
    DEALLOCATE ( tsig )
    DEALLOCATE ( tvsig )
    DEALLOCATE ( thetasig )
    DEALLOCATE ( mrsig )
    DEALLOCATE ( rhsig )
    DEALLOCATE ( rhodrysig ) 
    DEALLOCATE ( rhomoistsig )
    DEALLOCATE ( trap_top_ind )
    DEALLOCATE ( trap_bot_ind )
    DEALLOCATE ( weight_top_lin )
    DEALLOCATE ( weight_top_log )
    DEALLOCATE ( tvprs ) 
    DEALLOCATE ( mrprs )
    RETURN
  
  END SUBROUTINE deallocate_internal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE postproc_lfm
