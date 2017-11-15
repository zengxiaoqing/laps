MODULE wrf_lga

! A module containing necessary items to convert a WRF forecast
! file to LGA or to SWIM for the current LAPS domain

  USE map_utils
  USE wrf_netcdf
  USE time_utils
  USE constants
  USE horiz_interp
  USE mem_namelist, ONLY: r_missing_data
  USE mem_allsky
  IMPLICIT NONE

  PRIVATE
  integer                :: icentw, jcentw
  integer                :: icentl, jcentl
  REAL                   :: rmissingflag
  ! LAPS Pressure Levels 
  REAL, ALLOCATABLE      :: pr_laps(:)

  ! LGA variables
  REAL, ALLOCATABLE      :: ht(:,:,:)
  REAL, ALLOCATABLE      :: t3(:,:,:)
  REAL, ALLOCATABLE      :: sh(:,:,:)
  REAL, ALLOCATABLE      :: u3(:,:,:)
  REAL, ALLOCATABLE      :: v3(:,:,:)
  REAL, ALLOCATABLE      :: om(:,:,:)
  ! LGB Variables
  REAL, ALLOCATABLE      :: usf(:,:)
  REAL, ALLOCATABLE      :: vsf(:,:)
  REAL, ALLOCATABLE      :: tsf(:,:)
  REAL, ALLOCATABLE      :: tsk(:,:) ! surface skin temp. (added by Wei-Ting 130312)
  REAL, ALLOCATABLE      :: dsf(:,:)
  REAL, ALLOCATABLE      :: slp(:,:)
  REAL, ALLOCATABLE      :: psf(:,:)
  REAL, ALLOCATABLE      :: rsf(:,:)
  REAL, ALLOCATABLE      :: p(:,:)
  REAL, ALLOCATABLE      :: pcp(:,:) ! RAINNC+RAINC (added by Wei-Ting 130312)
  ! LAPS static variables
  REAL, ALLOCATABLE      :: topo_laps(:,:)
  REAL, ALLOCATABLE      :: lat(:,:)
  REAL, ALLOCATABLE      :: lon(:,:)
! INTEGER                :: nxl, nyl, nzl
  CHARACTER(LEN=200)     :: laps_data_root
  CHARACTER(LEN=10)      :: laps_domain_name  
  REAL                   :: redp_lvl
  ! WRF on pressure levels
  REAL, ALLOCATABLE      :: ht_wrfp(:,:,:)
  REAL, ALLOCATABLE      :: t3_wrfp(:,:,:)
  REAL, ALLOCATABLE      :: sh_wrfp(:,:,:)
  REAL, ALLOCATABLE      :: u3_wrfp(:,:,:)
  REAL, ALLOCATABLE      :: v3_wrfp(:,:,:)
  REAL, ALLOCATABLE      :: om_wrfp(:,:,:)
  REAL, ALLOCATABLE      :: qc_wrfp(:,:,:)
  REAL, ALLOCATABLE      :: qi_wrfp(:,:,:)
  REAL, ALLOCATABLE      :: qr_wrfp(:,:,:)
  REAL, ALLOCATABLE      :: qs_wrfp(:,:,:)
  REAL, ALLOCATABLE      :: aod_wrfp(:,:,:) ! extinction coefficient
  REAL, ALLOCATABLE      :: usf_wrf(:,:)
  REAL, ALLOCATABLE      :: vsf_wrf(:,:)
  REAL, ALLOCATABLE      :: tsf_wrf(:,:)
  REAL, ALLOCATABLE      :: tsk_wrf(:,:) ! surface skin temp. (added by Wei-Ting 130312)
  REAL, ALLOCATABLE      :: dsf_wrf(:,:)
  REAL, ALLOCATABLE      :: slp_wrf(:,:)
  REAL, ALLOCATABLE      :: psf_wrf(:,:)
  REAL, ALLOCATABLE      :: rsf_wrf(:,:)
  REAL, ALLOCATABLE      :: lmk_wrf(:,:)
  REAL, ALLOCATABLE      :: snc_wrf(:,:)
  REAL, ALLOCATABLE      :: sna_wrf(:,:)
  REAL, ALLOCATABLE      :: p_wrf(:,:)
  REAL, ALLOCATABLE      :: pcp_wrf(:,:) ! RAINNC+RAINC (added by Wei-Ting 130312)
  REAL, ALLOCATABLE      :: tvb_wrf(:,:)  ! Mean virtual temperature in lowest 60mb
  ! WRF on native variables
  REAL, ALLOCATABLE      :: pr_wrfs(:,:,:)
  REAL, ALLOCATABLE      :: ht_wrfs(:,:,:)
  REAL, ALLOCATABLE      :: dz_wrfs(:,:,:) ! If we want layer thicknesses
  REAL, ALLOCATABLE      :: aod_wrfs(:,:,:)
  REAL, ALLOCATABLE      :: t3_wrfs(:,:,:)
  REAL, ALLOCATABLE      :: sh_wrfs(:,:,:)
  REAL, ALLOCATABLE      :: u3_wrfs(:,:,:)
  REAL, ALLOCATABLE      :: v3_wrfs(:,:,:)
  REAL, ALLOCATABLE      :: om_wrfs(:,:,:)
  REAL, ALLOCATABLE      :: rho_wrfs(:,:,:) ! Density
  REAL, ALLOCATABLE      :: mr_wrfs(:,:,:) ! Mixing Ratio
  REAL, ALLOCATABLE      :: qc_wrfs(:,:,:) ! Cloud Liquid Mixing Ratio
  REAL, ALLOCATABLE      :: qi_wrfs(:,:,:) ! Cloud Ice Mixing Ratio
  REAL, ALLOCATABLE      :: qr_wrfs(:,:,:) ! Rain  Mixing Ratio
  REAL, ALLOCATABLE      :: qs_wrfs(:,:,:) ! Snow Mixing Ratio
  ! WRF static variables 
  INTEGER                :: cdf,cdp ! added cdp by Wei-Ting (130312)
  TYPE(proj_info)        :: wrfgrid
  REAL, ALLOCATABLE      :: topo_wrf(:,:)
  CHARACTER(LEN=19)      :: reftime
  INTEGER                :: tau_hr, tau_min,tau_sec
  INTEGER                :: itimestep,projcode,istat_aod
  INTEGER                :: nxw,nyw,nzw
  REAL                   :: dx_m, dy_m,dt
  REAL                   :: lat1_wrf, lon1_wrf
  REAL                   :: truelat1_wrf, truelat2_wrf, stdlon_wrf

  PUBLIC wrf2lga, wrf2swim
CONTAINS

  SUBROUTINE wrf2swim(wrffile_in,i4time,nxl,nyl,nzl,latl,lonl,pres_1d,land_frac,snow_cover,istatus)

     IMPLICIT NONE

     CHARACTER(LEN=150)           :: wrffile_in 
     CHARACTER(LEN=256)           :: wrffile
     INTEGER, INTENT(IN)          :: i4time
     INTEGER                      :: i4reftime
     CHARACTER(LEN=13)            :: reftime13
     INTEGER, INTENT(OUT)         :: istatus
     INTEGER                      :: k,k1000,bg_valid,icaller
     INTEGER,EXTERNAL             :: cvt_wfo_fname13_i4time
     INTEGER                      :: nxl,nyl,nzl
     REAL                         :: i_ll, j_ll, i_ul, j_ul, i_ur, j_ur, i_lr, j_lr 
     REAL                         :: latl(nxl,nyl),lonl(nxl,nyl),land_frac(nxl,nyl),pres_1d(nzl)
     REAL                         :: snow_cover(nxl,nyl)
     LOGICAL                      :: need_hinterp
      istatus = 1

      wrffile = wrffile_in

     
     ! Get some LAPS setup stuff
     rmissingflag = r_missing_data
     CALL find_domain_name(laps_data_root,laps_domain_name,istatus)
!    print *, "LAPS_DATA_ROOT = ", TRIM(laps_data_root)
!    print *, "DOMAIN NAME = ", TRIM(laps_domain_name)
     print *, "Dims:   ", nxl,nyl,nzl
     print *, "rmissingflag:   ", rmissingflag

    ! Allocate static fields

     ALLOCATE(pr_laps(nzl))
     ALLOCATE(lat(nxl,nyl))
     ALLOCATE(lon(nxl,nyl))
     ALLOCATE(topo_laps(nxl,nyl))
     CALL get_laps_domain(nxl,nyl,laps_domain_name,lat,lon,topo_laps,istatus)
     IF (istatus .NE. 1) THEN
       PRINT *, "Error reading LAPS static info."
       RETURN
     ENDIF
     pr_laps = pres_1d
     find_k1000: DO k = 1,nzl
      IF (NINT(pres_1d(k)) .EQ. 100000) THEN
        k1000 = k
        EXIT find_k1000
      ENDIF
     ENDDO find_k1000

     ! Print some config stuff
!    print *, "LAPS_DATA_ROOT = ", TRIM(laps_data_root)
!    print *, "DOMAIN NAME = ", TRIM(laps_domain_name)
     print *, "Dims:   ", nxl,nyl,nzl

     ! Get the WRF config
     CALL open_wrfnc(wrffile,cdf,istatus) 
     CALL get_wrf2_timeinfo(cdf,reftime,dt,itimestep,tau_hr,tau_min,tau_sec,istatus)
     reftime13 = reftime(1:4) // reftime(6:7) // reftime(9:13) // reftime(15:16)
     i4reftime = cvt_wfo_fname13_i4time(reftime13)
     bg_valid = i4reftime + tau_hr * 3600 + tau_min * 60 + tau_sec
     CALL get_wrf2_map(cdf,'T',projcode,lat1_wrf,lon1_wrf,stdlon_wrf, &
             truelat1_wrf,truelat2_wrf,dx_m,dy_m,nxw,nyw,nzw,istatus)
     CALL map_set(projcode,lat1_wrf,lon1_wrf,dx_m,stdlon_wrf,truelat1_wrf,truelat2_wrf, &
                  nxw,nyw,wrfgrid)
 
     ! Make sure LAPS domain covers the WRF domain, and see if it is an exact match
     CALL latlon_to_ij(wrfgrid,lat(1,1),lon(1,1),i_ll,j_ll)
     CALL latlon_to_ij(wrfgrid,lat(1,nyl),lon(1,nyl),i_ul,j_ul)
     CALL latlon_to_ij(wrfgrid,lat(nxl,nyl),lon(nxl,nyl),i_ur,j_ur)
     CALL latlon_to_ij(wrfgrid,lat(nxl,1),lon(nxl,1),i_lr,j_lr)
     print *, "Location of LAPS corners in WRF domain:"
     print *, i_ll,j_ll
     print *, i_ul,j_ul
     print *, i_ur,j_ur
     print *, i_lr,j_lr
     IF (NINT(i_ll) .LT. 1 .OR. NINT(j_ll) .LT. 1  .OR.  &
         NINT(i_ul) .LT. 1 .OR. NINT(j_ul) .GT. nyw .OR. &
         NINT(i_ur) .GT. nxw .OR. NINT(j_ur) .GT. nyw .OR. &
         NINT(i_lr) .GT. nxw .OR. NINT(j_lr) .LT. 1) THEN
       PRINT *, "LAPS Domain exceeds bounds of WRF background!"
       istatus = 0 
       RETURN
     ELSE
       need_hinterp = .true.
       IF (NINT(i_ll) .EQ. 1 .AND. NINT(j_ll) .EQ. 1 .AND. &
           NINT(i_ul) .EQ. 1 .AND. NINT(j_ul) .EQ. nyw .AND. &
           NINT(i_ur) .EQ. nxw .AND. NINT(j_ur) .EQ. nyw .AND. &
           NINT(i_lr) .EQ. nxw .AND. NINT(j_lr) .EQ. 1 ) THEN
         PRINT *, "Exact match between LAPS and background.  No hinterp needed!"
         need_hinterp = .false.
       ENDIF
     ENDIF

     icentl = nxl/2
     jcentl = nyl/2
     icentw = nxw/2
     jcentw = nyw/2 
     ! Get WRF on sigma
     ALLOCATE (pr_wrfs(nxw,nyw,nzw)) 
     ALLOCATE (ht_wrfs(nxw,nyw,nzw))
     ALLOCATE (dz_wrfs(nxw,nyw,nzw))
     ALLOCATE (aod_wrfs(nxw,nyw,nzw))
     ALLOCATE (t3_wrfs(nxw,nyw,nzw))
     ALLOCATE (sh_wrfs(nxw,nyw,nzw)) 
     ALLOCATE (mr_wrfs(nxw,nyw,nzw))
     ALLOCATE (rho_wrfs(nxw,nyw,nzw))
     ALLOCATE (qc_wrfs(nxw,nyw,nzw))
     ALLOCATE (qi_wrfs(nxw,nyw,nzw))
     ALLOCATE (qr_wrfs(nxw,nyw,nzw))
     ALLOCATE (qs_wrfs(nxw,nyw,nzw))
     ALLOCATE (usf_wrf(nxw,nyw))
     ALLOCATE (vsf_wrf(nxw,nyw))
     ALLOCATE (tsf_wrf(nxw,nyw))
     ALLOCATE (tsk_wrf(nxw,nyw)) ! surface skin temp. (added by Wei-Ting 130312)
     ALLOCATE (rsf_wrf(nxw,nyw))
     ALLOCATE (lmk_wrf(nxw,nyw))
     ALLOCATE (snc_wrf(nxw,nyw))
     ALLOCATE (sna_wrf(nxw,nyw))
     ALLOCATE (dsf_wrf(nxw,nyw))
     ALLOCATE (slp_wrf(nxw,nyw))
     ALLOCATE (psf_wrf(nxw,nyw))
     ALLOCATE (p_wrf(nxw,nyw))
     ALLOCATE (pcp_wrf(nxw,nyw)) ! RAINNC+RAINC (added by Wei-Ting 130312)
     ALLOCATE (topo_wrf(nxw,nyw))
     ALLOCATE (tvb_wrf(nxw,nyw))

     icaller = 1
     CALL fill_wrfs(icaller,istatus)
     IF (istatus .NE. 1) THEN
       PRINT *, "Problem getting WRF data"
       RETURN
     ENDIF

     ! Vertically interpolate to pressure levels
     PRINT *, "Allocating arrays for WRF on Pressure"
       ! Allocate wrfp
     ALLOCATE (ht_wrfp(nxw,nyw,nzl))
     ALLOCATE (t3_wrfp(nxw,nyw,nzl))
     ALLOCATE (sh_wrfp(nxw,nyw,nzl))
     ALLOCATE (qc_wrfp(nxw,nyw,nzl))
     ALLOCATE (qi_wrfp(nxw,nyw,nzl))
     ALLOCATE (qr_wrfp(nxw,nyw,nzl))
     ALLOCATE (qs_wrfp(nxw,nyw,nzl))
     ALLOCATE (aod_wrfp(nxw,nyw,nzl))
     ht_wrfp = rmissingflag
     t3_wrfp = rmissingflag
     sh_wrfp = rmissingflag
     qc_wrfp = rmissingflag
     qi_wrfp = rmissingflag
     qr_wrfp = rmissingflag
     qs_wrfp = rmissingflag
     aod_wrfp = rmissingflag
  
     ! Vertically interpolate
     PRINT *, "Calling vinterp_wrfarw2p ",nzl,pr_laps(1)
     CALL vinterp_wrfarw2p(icaller,nzl,istatus)
     IF (istatus .NE. 1) THEN
       PRINT *, "Problem vertically interpolating WRF data"
       RETURN
     ENDIF

     ! dealloc wrfs
     print *, "Deallocating WRF sigma var"
     DEALLOCATE(pr_wrfs,ht_wrfs,dz_wrfs,aod_wrfs,sh_wrfs,t3_wrfs,mr_wrfs, &  
        rho_wrfs,qc_wrfs,qi_wrfs,qr_wrfs,qs_wrfs)

     ! Horizontally interpolate to LAPS grid
       ! allocate lga/lgb

     print *, "Allocating LGA variables"
     ALLOCATE(ht(nxl,nyl,nzl))
     ALLOCATE(sh(nxl,nyl,nzl))
     ALLOCATE(t3(nxl,nyl,nzl))
     ! Allocate LGB (2d) variables
     print *, "Allocating LGB variables"
     ALLOCATE (usf(nxl,nyl))
     ALLOCATE (vsf(nxl,nyl))
     ALLOCATE (tsf(nxl,nyl))
     ALLOCATE (tsk(nxl,nyl)) ! surface skin temp. (added by Wei-Ting 130312)
     ALLOCATE (rsf(nxl,nyl))
     ALLOCATE (dsf(nxl,nyl))
     ALLOCATE (slp(nxl,nyl))
     ALLOCATE (psf(nxl,nyl))
     ALLOCATE (p  (nxl,nyl))
     ALLOCATE (pcp(nxl,nyl)) ! precitation (added by Wei-Ting 130312)
     usf = rmissingflag
     vsf = rmissingflag
     tsf = rmissingflag
     tsk = rmissingflag ! surface skin temp. (added by Wei-Ting 130312)
     dsf = rmissingflag
     slp = rmissingflag
     psf = rmissingflag
     p   = rmissingflag
     pcp = rmissingflag ! RAINNC+RAINC (added by Wei-Ting 130312)

     IF (need_hinterp) THEN
       PRINT *, "Problem in wrf2swim: hinterp indicated as needed" 
       istatus = 0
       RETURN
     ELSE
       ht = ht_wrfp
       t3 = t3_wrfp
       sh = sh_wrfp
       clwc_3d = qc_wrfp
       cice_3d = qi_wrfp
       rain_3d = qr_wrfp
       snow_3d = qs_wrfp
       if(istat_aod .eq. 1 .and. mode_aero_cld .eq. 3)then
         write(6,*)' Transferring AOD-3D to all-sky array'
         aod_3d(:,:,:) = aod_wrfp(:,:,:)
       endif
       psf = psf_wrf
       tsf = tsf_wrf
       tsk = tsk_wrf ! surface skin temp. (added by Wei-Ting 130312)
       dsf = dsf_wrf
       rsf = rsf_wrf
       usf = usf_wrf
       vsf = vsf_wrf
       pcp = pcp_wrf ! RAINNC+RAINC (added by Wei-Ting 130312)
       land_frac = lmk_wrf ! land mask
       snow_cover = snc_wrf * sna_wrf
     ENDIF 
!     pcp = 0 ! since pcp hasn't be used for now, assume that the value is 0

     print *, "Deallocating WRF press vars"
     DEALLOCATE (ht_wrfp,t3_wrfp,sh_wrfp,qc_wrfp,qi_wrfp,qr_wrfp,qs_wrfp,aod_wrfp)

     print *, "Deallocating WRF sfc vars"
     DEALLOCATE (usf_wrf,vsf_wrf,tsf_wrf,tsk_wrf,rsf_wrf,lmk_wrf,snc_wrf,sna_wrf,dsf_wrf,slp_wrf,&
         psf_wrf,p_wrf,pcp_wrf,tvb_wrf) ! added tsk_wrf & pcp_wrf by Wei-Ting (130312)
 
     ! Create MSLP and reduced pressure
     print *, "topo/slp/psf/p/tsf/tsk/dsf/rsf/usf/vsf",topo_laps(icentl,jcentl), &
       slp(icentl,jcentl),psf(icentl,jcentl),p(icentl,jcentl),tsf(icentl,jcentl), &
       tsk(icentl,jcentl),dsf(icentl,jcentl),rsf(icentl,jcentl),usf(icentl,jcentl), &
       vsf(icentl,jcentl) ! added tsk by Wei-Ting (130312)



     print *, "Deallocating LAPS static vars"
     DEALLOCATE (lat,lon,topo_laps,topo_wrf)

     ! Deallocate memory
     print *, "Deallocating lga vars"
     DEALLOCATE(ht,t3,sh,pr_laps)
     print *, "Deallocating lgb vars"
     DEALLOCATE(usf,vsf,tsf,tsk,rsf,dsf,slp,psf,p,pcp) ! added tsk & pcp by Wei-Ting (130312)
     PRINT *, "Successful processing of ", trim(wrffile) ! modified by Wei-Ting
     PRINT *, "Return from wrf2swim"
     RETURN 
  END SUBROUTINE wrf2swim

  SUBROUTINE wrf2lga(wrffile,i4time,cmodel,istatus )

     IMPLICIT NONE

     CHARACTER(LEN=256), INTENT(IN) :: wrffile(2) ! add 1 dim. and (LEN=255 -> LEN=256) by Wei-Ting (130312) to contain previous time
     CHARACTER(LEN=12), INTENT(IN) :: cmodel
     INTEGER, INTENT(IN)          :: i4time
     INTEGER                      :: i4reftime
     INTEGER                      :: nxl,nyl,nzl
     CHARACTER(LEN=13)            :: reftime13
     INTEGER, INTENT(OUT)         :: istatus
     INTEGER                      :: k,k1000,bg_valid,icaller
     INTEGER,EXTERNAL             :: cvt_wfo_fname13_i4time
     REAL                         :: i_ll, j_ll, i_ul, j_ul, i_ur, j_ur, i_lr, j_lr 
     LOGICAL                      :: need_hinterp
      istatus = 1

     
     ! Get some LAPS setup stuff
     CALL find_domain_name(laps_data_root,laps_domain_name,istatus)
     print *, "LAPS_DATA_ROOT = ", TRIM(laps_data_root)
     print *, "DOMAIN NAME = ", TRIM(laps_domain_name)
    
     CALL get_grid_dim_xy(nxl,nyl,istatus)
     IF (istatus .NE. 1) THEN
       PRINT *, "Could not get LAPS xy dims"
       RETURN
     ENDIF
     CALL get_laps_dimensions(nzl,istatus)
     IF (istatus .NE. 1) THEN
       PRINT *, "Could not get LAPS z dim"
       RETURN
     ENDIF
     print *, "Dims:   ", nxl,nyl,nzl
     CALL get_r_missing_data(rmissingflag,istatus) 
     CALL get_laps_redp(redp_lvl,istatus)
    ! Allocate static fields

     ALLOCATE(pr_laps(nzl))
     ALLOCATE(lat(nxl,nyl))
     ALLOCATE(lon(nxl,nyl))
     ALLOCATE(topo_laps(nxl,nyl))
     CALL get_laps_domain(nxl,nyl,laps_domain_name,lat,lon,topo_laps,istatus)
     IF (istatus .NE. 1) THEN
       PRINT *, "Error reading LAPS static info."
       RETURN
     ENDIF
     CALL get_pres_1d(i4time,nzl,pr_laps,istatus)
     IF (istatus .NE. 1) THEN
       PRINT *, "Error reading LAPS pressure levels"
       RETURN
     ENDIF
     find_k1000: DO k = 1,nzl
      IF (NINT(pr_laps(k)) .EQ. 100000) THEN
        k1000 = k
        EXIT find_k1000
      ENDIF
     ENDDO find_k1000
     ! Print some config stuff
     print *, "LAPS_DATA_ROOT = ", TRIM(laps_data_root)
     print *, "DOMAIN NAME = ", TRIM(laps_domain_name)
     print *, "Dims:   ", nxl,nyl,nzl


     ! Get the WRF config
     CALL open_wrfnc(wrffile(1),cdf,istatus) ! modified by Wei-Ting to use right time ( wrffile -> wrffile(1) )
     CALL open_wrfnc(wrffile(2),cdp,istatus) ! added by Wei-Ting to use previous time (just for calculate PCP)
     CALL get_wrf2_timeinfo(cdf,reftime,dt,itimestep,tau_hr,tau_min,tau_sec,istatus)
     reftime13 = reftime(1:4) // reftime(6:7) // reftime(9:13) // reftime(15:16)
     i4reftime = cvt_wfo_fname13_i4time(reftime13)
     bg_valid = i4reftime + tau_hr * 3600 + tau_min * 60 + tau_sec
     CALL get_wrf2_map(cdf,'T',projcode,lat1_wrf,lon1_wrf,stdlon_wrf, &
             truelat1_wrf,truelat2_wrf,dx_m,dy_m,nxw,nyw,nzw,istatus)
     CALL map_set(projcode,lat1_wrf,lon1_wrf,dx_m,stdlon_wrf,truelat1_wrf,truelat2_wrf, &
                  nxw,nyw,wrfgrid)
 
     ! Make sure LAPS domain covers the WRF domain, and see if it is an exact match
     CALL latlon_to_ij(wrfgrid,lat(1,1),lon(1,1),i_ll,j_ll)
     CALL latlon_to_ij(wrfgrid,lat(1,nyl),lon(1,nyl),i_ul,j_ul)
     CALL latlon_to_ij(wrfgrid,lat(nxl,nyl),lon(nxl,nyl),i_ur,j_ur)
     CALL latlon_to_ij(wrfgrid,lat(nxl,1),lon(nxl,1),i_lr,j_lr)
     print *, "Location of LAPS corners in WRF domain:"
     print *, i_ll,j_ll
     print *, i_ul,j_ul
     print *, i_ur,j_ur
     print *, i_lr,j_lr
     IF (NINT(i_ll) .LT. 1 .OR. NINT(j_ll) .LT. 1  .OR.  &
         NINT(i_ul) .LT. 1 .OR. NINT(j_ul) .GT. nyw .OR. &
         NINT(i_ur) .GT. nxw .OR. NINT(j_ur) .GT. nyw .OR. &
         NINT(i_lr) .GT. nxw .OR. NINT(j_lr) .LT. 1) THEN
       PRINT *, "LAPS Domain exceeds bounds of WRF background!"
       istatus = 0 
       RETURN
     ELSE
       need_hinterp = .true.
       IF (NINT(i_ll) .EQ. 1 .AND. NINT(j_ll) .EQ. 1 .AND. &
           NINT(i_ul) .EQ. 1 .AND. NINT(j_ul) .EQ. nyw .AND. &
           NINT(i_ur) .EQ. nxw .AND. NINT(j_ur) .EQ. nyw .AND. &
           NINT(i_lr) .EQ. nxw .AND. NINT(j_lr) .EQ. 1 ) THEN
         PRINT *, "Exact match between LAPS and background.  No hinterp needed!"
         need_hinterp = .false.
       ENDIF
     ENDIF

     icentl = nxl/2
     jcentl = nyl/2
     icentw = nxw/2
     jcentw = nyw/2 
     ! Get WRF on sigma
     ALLOCATE (pr_wrfs(nxw,nyw,nzw)) 
     ALLOCATE (ht_wrfs(nxw,nyw,nzw))
     ALLOCATE (t3_wrfs(nxw,nyw,nzw))
     ALLOCATE (sh_wrfs(nxw,nyw,nzw)) 
     ALLOCATE (u3_wrfs(nxw,nyw,nzw))
     ALLOCATE (v3_wrfs(nxw,nyw,nzw))
     ALLOCATE (om_wrfs(nxw,nyw,nzw)) 
     ALLOCATE (mr_wrfs(nxw,nyw,nzw))
     ALLOCATE (rho_wrfs(nxw,nyw,nzw))
     ALLOCATE (usf_wrf(nxw,nyw))
     ALLOCATE (vsf_wrf(nxw,nyw))
     ALLOCATE (tsf_wrf(nxw,nyw))
     ALLOCATE (tsk_wrf(nxw,nyw)) ! surface skin temp. (added by Wei-Ting 130312)
     ALLOCATE (rsf_wrf(nxw,nyw))
     ALLOCATE (dsf_wrf(nxw,nyw))
     ALLOCATE (slp_wrf(nxw,nyw))
     ALLOCATE (psf_wrf(nxw,nyw))
     ALLOCATE (p_wrf(nxw,nyw))
     ALLOCATE (pcp_wrf(nxw,nyw)) ! RAINNC+RAINC (added by Wei-Ting 130312)
     ALLOCATE (topo_wrf(nxw,nyw))
     ALLOCATE (tvb_wrf(nxw,nyw))

     icaller = 2
     CALL fill_wrfs(icaller,istatus)
     IF (istatus .NE. 1) THEN
       PRINT *, "Problem getting WRF data"
       RETURN
     ENDIF

     ! Vertically interpolate to pressure levels
     PRINT *, "Allocating arrays for WRF on Pressure"
       ! Allocate wrfp
     ALLOCATE (ht_wrfp(nxw,nyw,nzl))
     ALLOCATE (t3_wrfp(nxw,nyw,nzl))
     ALLOCATE (sh_wrfp(nxw,nyw,nzl))
     ALLOCATE (u3_wrfp(nxw,nyw,nzl))
     ALLOCATE (v3_wrfp(nxw,nyw,nzl))
     ALLOCATE (om_wrfp(nxw,nyw,nzl))
     ht_wrfp = rmissingflag
     t3_wrfp = rmissingflag
     u3_wrfp = rmissingflag
     v3_wrfp = rmissingflag
     om_wrfp = rmissingflag
  
     ! Vertically interpolate
     PRINT *, "Calling vinterp_wrfarw2p ",nzl,pr_laps(1)
     CALL vinterp_wrfarw2p(icaller,nzl,istatus)
     IF (istatus .NE. 1) THEN
       PRINT *, "Problem vertically interpolating WRF data"
       RETURN
     ENDIF

     ! dealloc wrfs
     print *, "Deallocating WRF sigma var"
     DEALLOCATE(pr_wrfs,ht_wrfs,sh_wrfs,t3_wrfs,u3_wrfs,v3_wrfs,om_wrfs,mr_wrfs, &
        rho_wrfs)

     ! Horizontally interpolate to LAPS grid
       ! allocate lga/lgb

     print *, "Allocating LGA variables"
     ALLOCATE(ht(nxl,nyl,nzl))
     ALLOCATE(sh(nxl,nyl,nzl))
     ALLOCATE(t3(nxl,nyl,nzl))
     ALLOCATE(u3(nxl,nyl,nzl))
     ALLOCATE(v3(nxl,nyl,nzl))
     ALLOCATE(om(nxl,nyl,nzl))
     ! Allocate LGB (2d) variables
     print *, "Allocating LGB variables"
     ALLOCATE (usf(nxl,nyl))
     ALLOCATE (vsf(nxl,nyl))
     ALLOCATE (tsf(nxl,nyl))
     ALLOCATE (tsk(nxl,nyl)) ! surface skin temp. (added by Wei-Ting 130312)
     ALLOCATE (rsf(nxl,nyl))
     ALLOCATE (dsf(nxl,nyl))
     ALLOCATE (slp(nxl,nyl))
     ALLOCATE (psf(nxl,nyl))
     ALLOCATE (p  (nxl,nyl))
     ALLOCATE (pcp(nxl,nyl)) ! precitation (added by Wei-Ting 130312)
     usf = rmissingflag
     vsf = rmissingflag
     tsf = rmissingflag
     tsk = rmissingflag ! surface skin temp. (added by Wei-Ting 130312)
     dsf = rmissingflag
     slp = rmissingflag
     psf = rmissingflag
     p   = rmissingflag
     pcp = rmissingflag ! RAINNC+RAINC (added by Wei-Ting 130312)

     IF (need_hinterp) THEN
       ht = rmissingflag
       sh = rmissingflag
       t3 = rmissingflag
       u3 = rmissingflag
       v3 = rmissingflag
       om = rmissingflag
       usf = rmissingflag
       vsf = rmissingflag
       tsf = rmissingflag
       tsk = rmissingflag ! surface skin temp. (added by Wei-Ting 130312)
       dsf = rmissingflag
       slp = rmissingflag
       psf = rmissingflag
       p   = rmissingflag
       pcp = rmissingflag ! RAINNC+RAINC (added by Wei-Ting 130312)
                                                                                                                            

       ! call hinterplga
       PRINT *, "Calling hinterp_wrf2lga"
       CALL hinterp_wrf2lga(nxl,nyl,nzl,istatus)
       IF (istatus .NE. 1) THEN 
         PRINT *, "Problem in hinterp_wrf2lga" 
         istatus = 0
         RETURN
       ENDIF 
     ELSE
       ht = ht_wrfp
       t3 = t3_wrfp
       sh = sh_wrfp
       u3 = u3_wrfp
       v3 = v3_wrfp
       om = om_wrfp
       psf = psf_wrf
       tsf = tsf_wrf
       tsk = tsk_wrf ! surface skin temp. (added by Wei-Ting 130312)
       dsf = dsf_wrf
       rsf = rsf_wrf
       usf = usf_wrf
       vsf = vsf_wrf
       pcp = pcp_wrf ! RAINNC+RAINC (added by Wei-Ting 130312)
     ENDIF 
!     pcp = 0 ! since pcp hasn't be used for now, assume that the value is 0

     print *, "Deallocating WRF press vars"
     DEALLOCATE (ht_wrfp,t3_wrfp,sh_wrfp,u3_wrfp,v3_wrfp,om_wrfp)

     print *, "Deallocating WRF sfc vars"
     DEALLOCATE (usf_wrf,vsf_wrf,tsf_wrf,tsk_wrf,rsf_wrf,dsf_wrf,slp_wrf,&
         psf_wrf,p_wrf,pcp_wrf,tvb_wrf) ! added tsk_wrf & pcp_wrf by Wei-Ting (130312)
 
     ! Create MSLP and reduced pressure
     print *, "Creating reduced pressure arrays"
     CALL make_derived_pressures(nxl,nyl,nzl)
     print *, "topo/slp/psf/p/tsf/tsk/dsf/rsf/usf/vsf",topo_laps(icentl,jcentl), &
       slp(icentl,jcentl),psf(icentl,jcentl),p(icentl,jcentl),tsf(icentl,jcentl), &
       tsk(icentl,jcentl),dsf(icentl,jcentl),rsf(icentl,jcentl),usf(icentl,jcentl), &
       vsf(icentl,jcentl) ! added tsk by Wei-Ting (130312)



     print *, "Deallocating LAPS static vars"
     DEALLOCATE (lat,lon,topo_laps,topo_wrf)

     ! Write the data out
     print *, "Writing LGA"
     CALL write_lga(nxl,nyl,nzl,i4reftime,bg_valid,cmodel,rmissingflag, &
           pr_laps*0.01,ht,t3,sh,u3,v3,om,istatus)
     IF(istatus .NE. 1) THEN
       print *, "Error writing LGA"
       RETURN
     ENDIF

     print *, "Writing LGB" ! modified by Wei-Ting (LGA -> LGB)
     CALL write_lgb(nxl,nyl,i4reftime,bg_valid,cmodel,rmissingflag, &
           usf,vsf,tsf,tsk,rsf,psf,slp,dsf,p,pcp,istatus) ! added tsk & pcp by Wei-Ting (130312)
     IF(istatus .NE. 1) THEN
       print *, "Error writing LGB"
       RETURN
     ENDIF


     ! Deallocate memory
     print *, "Deallocating lga vars"
     DEALLOCATE(ht,t3,sh,u3,v3,om,pr_laps)
     print *, "Deallocating lgb vars"
     DEALLOCATE(usf,vsf,tsf,tsk,rsf,dsf,slp,psf,p,pcp) ! added tsk & pcp by Wei-Ting (130312)
     PRINT *, "Successful processing of ", trim(wrffile(1)) ! modified by Wei-Ting
     RETURN 
  END SUBROUTINE wrf2lga
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE fill_wrfs(icaller,istatus)
 
    IMPLICIT NONE
    INTEGER, INTENT(OUT)  :: istatus
    INTEGER               :: status,i,j,k,icaller
    REAL, ALLOCATABLE     :: dum3d(:,:,:)
    REAL, ALLOCATABLE     :: dum3df(:,:,:)
    REAL, ALLOCATABLE     :: dum3df2(:,:,:)
    REAL, ALLOCATABLE     :: dum2dt1(:,:) ! added by Wei-Ting (130312) to get RAINNC
    REAL, ALLOCATABLE     :: dum2dt2(:,:) ! added by Wei-Ting (130312) to get RAINNC
    REAL, EXTERNAL        :: mixsat, relhum, dewpt2
    REAL                  :: rh 
    REAL                  :: tvbar, tvbar_nlevs
    REAL, PARAMETER       :: tvbar_thick = 6000.
    ! Varialbles have already been allocated by our driver routine
    ! so just start getting them
    istatus = 1
    PRINT *, " Allocating arrays ",icaller
    ALLOCATE(dum3d(nxw,nyw,nzw))
    ALLOCATE(dum3df(nxw,nyw,nzw+1))
    ALLOCATE(dum3df2(nxw,nyw,nzw+1))
    ALLOCATE(dum2dt1(nxw,nyw))
    ALLOCATE(dum2dt2(nxw,nyw))

    ! Get 3D pressure array
    ! Get pressures
    PRINT *, "Getting PB"
    CALL get_wrfnc_3d(cdf,"PB","T",nxw,nyw,nzw,1,dum3d,status)
    IF (status .GT. 0) THEN
      PRINT *, "+++CRITICAL:  Could not get base pressure!"
      istatus = 0
      RETURN
    ENDIF
    pr_wrfs = dum3d
                        
    PRINT *, "Getting P"                                                                    
    CALL get_wrfnc_3d(cdf,"P","T",nxw,nyw,nzw,1,dum3d,status)
    IF (status .GT. 0) THEN
      PRINT *, "+++CRITICAL:  Could not get pert pressure!"
      istatus = 0
      RETURN
    ENDIF
    pr_wrfs = pr_wrfs + dum3d
    print *, "Min/Max WRF 3D Pressure: ",minval(pr_wrfs),maxval(pr_wrfs)
  
    ! Get heights
    print *, "Getting PHB" 
    CALL get_wrfnc_3d(cdf,"PHB","T",nxw,nyw,nzw+1,1,dum3df,status)
    IF (status.NE.0) THEN
      PRINT *, 'Could not properly obtain WRF base-state geopotential.'
      istatus = 0
      RETURN
    ENDIF
    dum3df2 = dum3df
  
    PRINT *, "Getting PH"
    CALL get_wrfnc_3d(cdf,"PH","T",nxw,nyw,nzw+1,1,dum3df,status)
    IF (status.NE.0) THEN
      PRINT *, 'Could not properly obtain WRF geopotential.'
      istatus = 0
      RETURN
    ENDIF
    PRINT *,"Destaggering (vertically) heights"
    dum3df2 = (dum3df2 + dum3df) / grav
    DO k = 1,nzw
      ht_wrfs(:,:,k) = 0.5 * (dum3df2(:,:,k) + dum3df2(:,:,k+1))
      if(icaller .eq. 1)then
        dz_wrfs(:,:,k) = dum3df2(:,:,k+1) - dum3df2(:,:,k)
      endif
    ENDDO
  
    PRINT *, "Getting AOD"
    CALL get_wrfnc_3d(cdf,"TAOD5503D","T",nxw,nyw,nzw,1,dum3df,status)
    IF (status.NE.0) THEN
      PRINT *, 'Could not properly obtain WRF AOD - set to missing'
      aod_wrfs(:,:,:) = r_missing_data
      istat_aod = 0
    ELSE
      PRINT *, 'Success Reading WRF AOD - divide by dz'
      aod_wrfs(:,:,:) = dum3df(:,:,:) / dz_wrfs(:,:,:)
      print *, "Min/Max WRF 3D AOD: ",minval(aod_wrfs),maxval(aod_wrfs)
      istat_aod = 1
    ENDIF

    ! Get theta and convert to temperature
    PRINT *, "Getting Theta"
    CALL get_wrfnc_3d(cdf, "T","T",nxw,nyw,nzw,1,dum3d,status)
    IF (status.NE.0) THEN
      PRINT *, 'Could not properly obtain WRF perturbation theta.'
      istatus = 0
      RETURN
    ENDIF

    PRINT *, "Computing temp"
    dum3d = dum3d + 300.
    DO k = 1, nzw
      DO j = 1, nyw
        DO i = 1,nxw
          t3_wrfs(i,j,k) = dum3d(i,j,k)/ ((100000./pr_wrfs(i,j,k))**kappa)
        ENDDO
      ENDDO
    ENDDO

    ! Get Q on sigma
    PRINT *, "Getting Q"
    CALL get_wrfnc_3d(cdf, "QVAPOR","T",nxw,nyw,nzw,1,mr_wrfs,status)
    IF (status.NE.0) THEN
      PRINT *, 'Could not properly obtain WRF mixing ratio.'
      istatus = 0
      RETURN
    ENDIF

    ! Derive specific humidity
    PRINT *, "Computing SH"
    sh_wrfs = mr_wrfs / (1. + mr_wrfs)

    if(icaller .eq. 1)then
      PRINT *, "Getting QC"
      CALL get_wrfnc_3d(cdf, "QCLOUD","T",nxw,nyw,nzw,1,qc_wrfs,status)
      IF (status.NE.0) THEN
        PRINT *, 'Could not properly obtain WRF qc mixing ratio.'
        istatus = 0
        RETURN
      ENDIF

      PRINT *, "Getting QI"
      CALL get_wrfnc_3d(cdf, "QICE","T",nxw,nyw,nzw,1,qi_wrfs,status)
      IF (status.NE.0) THEN
        PRINT *, 'Could not properly obtain WRF qi mixing ratio.'
        istatus = 0
        RETURN
      ENDIF

      PRINT *, "Getting QR"
      CALL get_wrfnc_3d(cdf, "QRAIN","T",nxw,nyw,nzw,1,qr_wrfs,status)
      IF (status.NE.0) THEN
        PRINT *, 'Could not properly obtain WRF qr mixing ratio.'
        istatus = 0
        RETURN
      ENDIF

      PRINT *, "Getting QS"
      CALL get_wrfnc_3d(cdf, "QSNOW","T",nxw,nyw,nzw,1,qs_wrfs,status)
      IF (status.NE.0) THEN
        PRINT *, 'Could not properly obtain WRF qs mixing ratio.'
        istatus = 0
        RETURN
      ENDIF

    endif

    if(icaller .eq. 2)then
    ! Get U on sigma   
      PRINT *, "Getting U"
      CALL get_wrfnc_3d(cdf, "U","T",nxw,nyw,nzw,1,u3_wrfs,status)
      IF (status.NE.0) THEN
        PRINT *, 'Could not properly obtain WRF U-comp.'
        istatus = 0
        RETURN
      ENDIF

    ! Get V on sigma
      PRINT *, "Getting V"
      CALL get_wrfnc_3d(cdf, "V","T",nxw,nyw,nzw,1,v3_wrfs,status)
      IF (status.NE.0) THEN
        PRINT *, 'Could not properly obtain WRF V-comp.'
        istatus = 0
        RETURN
      ENDIF

    ! Get W on sigma
      PRINT *, "Getting W"
      CALL get_wrfnc_3d(cdf, "W","T",nxw,nyw,nzw+1,1,dum3df,status)
      IF (status.NE.0) THEN
        PRINT *, 'Could not properly obtain WRF W-comp'
        istatus = 0
        RETURN 
      ENDIF
      PRINT*, "Destaggering (vertically) w"
      DO k = 1,nzw
        om_wrfs(:,:,k) = 0.5*(dum3df(:,:,k)+dum3df(:,:,k+1))
      ENDDO
    endif 

    ! Now, derive density, virtual potential temp, and omega
    PRINT *, "Computing Rho, Theta-V, and Omega"
    DO k = 1,nzw
      DO j=1,nyw
        DO i = 1,nxw
          rho_wrfs(i,j,k) = pr_wrfs(i,j,k) / ( r * t3_wrfs(i,j,k)*(1.+0.61*sh_wrfs(i,j,k)))
          if(icaller .eq. 2)then
            om_wrfs(i,j,k) = -1. * rho_wrfs(i,j,k) * grav * om_wrfs(i,j,k)
          endif
        ENDDO
      ENDDO
    ENDDO

    if(icaller .eq. 1)then ! convert mixing ratio to density
      PRINT *, "Convert QC,QI,QR,QS to density"
      qc_wrfs(:,:,:) = qc_wrfs(:,:,:) * rho_wrfs(:,:,:)
      qi_wrfs(:,:,:) = qi_wrfs(:,:,:) * rho_wrfs(:,:,:)
      qr_wrfs(:,:,:) = qr_wrfs(:,:,:) * rho_wrfs(:,:,:)
      qs_wrfs(:,:,:) = qs_wrfs(:,:,:) * rho_wrfs(:,:,:)
   endif
   
    ! Get surface fields
    
    PRINT *, "Getting WRF TOPO"
    CALL get_wrfnc_2d(cdf, "HGT","T",nxw,nyw,1,topo_wrf,status)
    IF (status .NE. 0) THEN
      PRINT *, "Could not get topo from WRF"
      istatus = 0 
      RETURN
    ENDIF
    
    PRINT *, "Getting WRF LANDMASK"
    CALL get_wrfnc_2d(cdf, "LANDMASK","T",nxw,nyw,1,lmk_wrf,status)
    IF (status .NE. 0) THEN
      PRINT *, "Could not get land mask from WRF"
      istatus = 0 
      RETURN
    ENDIF
    
    PRINT *, "Getting WRF SNOWC"
    CALL get_wrfnc_2d(cdf, "SNOWC","T",nxw,nyw,1,snc_wrf,status)
    IF (status .NE. 0) THEN
      PRINT *, "Could not get snow cover from WRF"
      istatus = 0 
      RETURN
    ENDIF
    
    PRINT *, "Getting WRF SNOALB"
    CALL get_wrfnc_2d(cdf, "SNOALB","T",nxw,nyw,1,sna_wrf,status)
    IF (status .NE. 0) THEN
      PRINT *, "Could not get snow albedo from WRF"
      istatus = 0 
      RETURN
    ENDIF

    if(icaller .eq. 2)then
      PRINT *, "Getting USF"
      usf_wrf = u3_wrfs(:,:,1)
      PRINT *, "Getting VSF"
      vsf_wrf = v3_wrfs(:,:,1)
    endif

    PRINT *, "Getting PSF"
    CALL get_wrfnc_2d(cdf, "PSFC","T",nxw,nyw,1,psf_wrf,status)
    IF ((status .NE. 0).OR.(MAXVAL(psf_wrf) .LT. 10000.))THEN
      PRINT *, "Could not get PSFC, using lowest sigma level"
      psf_wrf = pr_wrfs(:,:,1)
    ENDIF
 
    PRINT *, "Getting T2" 
    CALL get_wrfnc_2d(cdf, "T2","T",nxw,nyw,1,tsf_wrf,status)
    IF ((status .NE. 0).OR.(MAXVAL(tsf_wrf) .LT. 100.))THEN
      PRINT *, "Could not get T2, using lowest sigma level"
      tsf_wrf = t3_wrfs(:,:,1)
    ENDIF
    
    ! added tsk(skin temp.) by Wei-Ting (130312)
    PRINT *, "Getting TSK" 
    CALL get_wrfnc_2d(cdf, "TSK","T",nxw,nyw,1,tsk_wrf,status)
    IF ((status .NE. 0).OR.(MAXVAL(tsk_wrf) .LT. 100.))THEN
      PRINT *, "Could not get TSK, using lowest sigma level"
      tsk_wrf = t3_wrfs(:,:,1)
    ENDIF
    
    ! added PCP (RAINNC+RAINC) by Wei-Ting (130312) & Modified (130326)
    PRINT *, "Getting Precipitation"
    PRINT *, "!!!!! This Precipitaion is an accumulation per N hours. !!!!!"
    PRINT *, "!!!!! N depends on the time difference of each WRFOUT.  !!!!!"
    PRINT *, "   Getting RAINNC(t)"
    CALL get_wrfnc_2d(cdf,"RAINNC","T",nxw,nyw,1,dum2dt2,status)
    IF (status .NE. 0) THEN
      PRINT *, "   Could not get RAINNC(t), setting the value = 0"
      dum2dt2 = 0
    ENDIF
    pcp_wrf = dum2dt2
    PRINT *, "   Getting RAINC(t)"
    CALL get_wrfnc_2d(cdf,"RAINC","T",nxw,nyw,1,dum2dt2,status)
    IF (status .NE. 0) THEN
      PRINT *, "   Could not get RAINC(t), setting the value = 0"
      dum2dt2 = 0
    ENDIF
    pcp_wrf = pcp_wrf+dum2dt2

    PRINT *, "   Getting RAINNC(t-1)"
    IF (cdp .LT. 0 .AND. cdf .GT. 0) THEN
      PRINT *, "   Could not get RAINNC(t-1), maybe result from t = initial time!"
      PRINT *, "   Set RAINNC(t-1) = 0"
      dum2dt1 = 0
    ELSE
      CALL get_wrfnc_2d(cdp,"RAINNC","T",nxw,nyw,1,dum2dt1,status)
      IF (status .NE. 0) THEN
         PRINT *, "   Could not get RAINNC(t-1), setting the value = 0"
         dum2dt1 = 0
      ENDIF
    ENDIF
    pcp_wrf = pcp_wrf-dum2dt1
    PRINT *, "   Getting RAINC(t-1)"
    IF (cdp .LT. 0 .AND. cdf .GT. 0) THEN
      PRINT *, "   Could not get RAINC(t-1), maybe result from t = initial time!"
      PRINT *, "   Set RAINC(t-1) = 0"
      dum2dt1 = 0
    ELSE
      CALL get_wrfnc_2d(cdp,"RAINC","T",nxw,nyw,1,dum2dt1,status)
      IF (status .NE. 0) THEN
         PRINT *, "   Could not get RAINC(t-1), setting the value = 0"
         dum2dt1 = 0
      ENDIF
    ENDIF
    pcp_wrf = pcp_wrf-dum2dt1
    where ( pcp_wrf < 0 ) ; pcp_wrf = 0 ; endwhere ! keep pcp >= 0
    print *, "Min/Max WRF Precipitation : ",minval(pcp_wrf),maxval(pcp_wrf)
    ! End of reading RAINNC+RAINC

    ! qvapor at 2m
    PRINT *, "Getting Q2"
    CALL get_wrfnc_2d(cdf, "Q2","T",nxw,nyw,1,rsf_wrf,status)
    IF ((status .NE. 0).OR.(MAXVAL(rsf_wrf) .LT. 0.0001))THEN
      PRINT *, "Could not get Q2, using lowest sigma level"
      rsf_wrf = sh_wrfs(:,:,1)
    ELSE 

      ! Because 2m qv and T are derived from the PBL scheme and
      ! the WRF is apparently not checking for saturation, clean this
      ! up now
      PRINT *, "Checking Q2 for supersaturation"
      DO j = 1, nyw
        DO i= 1, nxw
          rsf_wrf(i,j) = MIN(rsf_wrf(i,j),mixsat(rsf_wrf(i,j),psf_wrf(i,j)))
          ! Compute dewpoint
          rh = relhum(tsf_wrf(i,j),rsf_wrf(i,j),psf_wrf(i,j))
          dsf_wrf(i,j) = dewpt2(tsf_wrf(i,j),rh)
          ! Compute tvbar
          tvbar = 0.
          tvbar_nlevs = 0.
          comptvb:  DO k = 1,nzw
            IF ((psf_wrf(i,j)-pr_wrfs(i,j,k)) .LE. tvbar_thick) THEN
               tvbar = tvbar + (t3_wrfs(i,j,k)*(1.+0.61*sh_wrfs(i,j,k)) )
               tvbar_nlevs = tvbar_nlevs + 1.
            ELSE
              exit comptvb
            ENDIF
          ENDDO comptvb
          tvb_wrf(i,j) = tvbar / tvbar_nlevs
        ENDDO
      ENDDO
    ENDIF 
    ! Smooth the tvb_wrf field
    CALL smooth2(nxw,nyw,4,tvb_wrf)
    ! Convert mr to sh
    rsf_wrf(:,:) = rsf_wrf(:,:)/(1. + rsf_wrf(:,:))

    ! Diagnostics
    PRINT *, "WRF Sigma data from center of WRF domain"
    if(icaller .eq. 1)then ! called from 'wrf2swim'
      print *, "K   PRESS     HEIGHT   DZ     TEMP   SH         QC        QI       QR       QS      AOD"
      print *, "--- --------  -------  -----  -----  -------    -------   ------   -----    -----   ------"
    else                   ! called from 'wrf2lga'
      print *, "K   PRESS     HEIGHT   TEMP   SH       U       V      OM"
      print *, "--- --------  -------  -----  -------  ------  ------ -----------"
    endif
    DO k = 1,nzw
      if(icaller .eq. 1)then ! called from 'wrf2swim'
        PRINT ('(I3,1x,F8.1,2x,F7.0,2x,F6.1,2x,f5.1,2x,F7.5,2x,4F9.5,2x,F8.6)'), &
          k,pr_wrfs(icentw,jcentw,k),ht_wrfs(icentw,jcentw,k),dz_wrfs(icentw,jcentw,k), t3_wrfs(icentw,jcentw,k), &
          sh_wrfs(icentw,jcentw,k),qc_wrfs(icentw,jcentw,k),qi_wrfs(icentw,jcentw,k), &
          qr_wrfs(icentw,jcentw,k),qs_wrfs(icentw,jcentw,k),aod_wrfs(icentw,jcentw,k)
      else                   ! called from 'wrf2lga'
        PRINT ('(I3,1x,F8.1,2x,F7.0,2x,F5.1,2x,F7.5,2x,F6.1,2x,F6.1,2x,F11.8)'), &
          k,pr_wrfs(icentw,jcentw,k),ht_wrfs(icentw,jcentw,k), t3_wrfs(icentw,jcentw,k), &
          sh_wrfs(icentw,jcentw,k),u3_wrfs(icentw,jcentw,k),v3_wrfs(icentw,jcentw,k),&
          om_wrfs(icentw,jcentw,k)
      endif
    ENDDO
    PRINT *, "SFC:"
    PRINT ('(F8.1,2x,F7.0,2x,F5.1,2x,F5.1,2x,F5.1,2x,F5.1,2x,F7.5,2x,F6.1,2x,F6.1)'), &
        psf_wrf(icentw,jcentw),topo_wrf(icentw,jcentw),tsf_wrf(icentw,jcentw), &
        tsk_wrf(icentw,jcentw),dsf_wrf(icentw,jcentw),tvb_wrf(icentw,jcentw), &
        rsf_wrf(icentw,jcentw), usf_wrf(icentw,jcentw), vsf_wrf(icentw,jcentw)
        ! added tsk_wrf by Wei-Ting (130312)
    
    PRINT *, "Deallocating arrays"
    DEALLOCATE(dum3d,dum3df,dum3df2,dum2dt1,dum2dt2)
    RETURN
  END SUBROUTINE fill_wrfs 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE vinterp_wrfarw2p(icaller,nzl,istatus)

    IMPLICIT NONE

    INTEGER, INTENT(OUT)  :: istatus
    INTEGER  :: i,j,k,ks,kp, ksb,kst,icaller,nzl
    REAL     :: lpb,lpt,lp, wgtb, wgtt
    REAL,PARAMETER     :: dTdlnPBase = 50.0
    REAL               :: tvbot, tvbar, deltalnp
    REAL               :: dz
    istatus = 1

    write(6,*)' Subroutine vinterp_wrfarw2p: ',icaller,nzl,pr_laps(1)

    ! Loop over horizontal domain, interpolating vertically for each
    ! column
    DO j = 1, nyw
      DO i = 1, nxw   

        pressloop: DO kp = 1,nzl
  
          ! Initialize kst and ksb, which will hold the vertical sigma
          ! index values of the top and bottom bounding layers
          kst = 0
          ksb = 0 

          ! Find bounding levels in raw data
          sigmaloop: DO ks = 1,nzw
            IF (pr_wrfs(i,j,ks) .LE. pr_laps(kp)) THEN   

              kst = ks
              ksb = kst - 1
              EXIT sigmaloop
            ENDIF
          ENDDO sigmaloop

          IF (kst .GT. 1) THEN ! Interpolate between two bounding points
            lp = ALOG(pr_laps(kp))
            lpt = ALOG(pr_wrfs(i,j,kst))
            lpb = ALOG(pr_wrfs(i,j,ksb))
            wgtb = (lpt - lp) / (lpt - lpb)
            wgtt = 1.0 - wgtb

            ! Height
            ht_wrfp(i,j,kp) = wgtb * ht_wrfs(i,j,ksb) + &
                              wgtt * ht_wrfs(i,j,kst)

            ! Temp
            t3_wrfp(i,j,kp) = wgtb * t3_wrfs(i,j,ksb) + &
                              wgtt * t3_wrfs(i,j,kst)

            ! SH
            sh_wrfp(i,j,kp) = wgtb * sh_wrfs(i,j,ksb) + &
                              wgtt * sh_wrfs(i,j,kst)

            if(icaller .eq. 1)then ! called from 'wrf2swim'
              qc_wrfp(i,j,kp) = wgtb * qc_wrfs(i,j,ksb) + &
                                wgtt * qc_wrfs(i,j,kst)

              qi_wrfp(i,j,kp) = wgtb * qi_wrfs(i,j,ksb) + &
                                wgtt * qi_wrfs(i,j,kst)

              qr_wrfp(i,j,kp) = wgtb * qr_wrfs(i,j,ksb) + &
                                wgtt * qr_wrfs(i,j,kst)

              qs_wrfp(i,j,kp) = wgtb * qs_wrfs(i,j,ksb) + &
                                wgtt * qs_wrfs(i,j,kst)

              if(istat_aod .eq. 1)then
                aod_wrfp(i,j,kp) = wgtb * aod_wrfs(i,j,ksb) + &
                                   wgtt * aod_wrfs(i,j,kst)
              endif

            else

            ! U3
              u3_wrfp(i,j,kp) = wgtb * u3_wrfs(i,j,ksb) + &
                                wgtt * u3_wrfs(i,j,kst)

            ! V3
              v3_wrfp(i,j,kp) = wgtb * v3_wrfs(i,j,ksb) + &
                                wgtt * v3_wrfs(i,j,kst)
 
            ! OM
              om_wrfp(i,j,kp) = wgtb * om_wrfs(i,j,ksb) + &
                                wgtt * om_wrfs(i,j,kst)
            endif


          ELSEIF (kst .EQ. 1) THEN ! Extrapolate downward
            lpt = ALOG(pr_wrfs(i,j,kst))
            lpb = ALOG(pr_laps(kp))
            deltalnp = lpb - lpt
            tvbot = tvb_wrf(i,j) + deltalnp*dTdlnPBase
            tvbar = 0.5*(tvb_wrf(i,j) + tvbot)
            ! Height
 
            IF ((pr_laps(kp) - pr_wrfs(i,j,1)).LT. 500.) THEN
              ! Very small difference in pressures, so
              ! assume 10 m per 100 Pa, because
              ! hypsometric eq breaks down in these cases
              dz =  0.1 * (pr_laps(kp) - pr_wrfs(i,j,1))
            ELSE
              dz = tvbar * rog * ALOG(pr_laps(kp)/pr_wrfs(i,j,1))
            ENDIF
            ht_wrfp(i,j,kp) = ht_wrfs(i,j,1) - dz
            ! Temp
            t3_wrfp(i,j,kp) = tvbot/(1.+0.61*mr_wrfs(i,j,1))
            ! SH
            sh_wrfp(i,j,kp) = sh_wrfs(i,j,1)

            if(icaller .eq. 1)then ! called from 'wrf2swim'
              qc_wrfp(i,j,kp) = qc_wrfs(i,j,1)
              qi_wrfp(i,j,kp) = qi_wrfs(i,j,1)
              qr_wrfp(i,j,kp) = qr_wrfs(i,j,1)
              qs_wrfp(i,j,kp) = qs_wrfs(i,j,1)
              if(istat_aod .eq. 1)then
                aod_wrfp(i,j,kp) = aod_wrfs(i,j,1)
              endif
              
            else
            
            ! U3
              u3_wrfp(i,j,kp) = u3_wrfs(i,j,1)

            ! V3
              v3_wrfp(i,j,kp) = v3_wrfs(i,j,1)

            ! OM
              om_wrfp(i,j,kp) = 0.

            endif

          ELSE ! kst never got set .. extrapolate upward
            ! Assume isothermal (above tropopause)
            t3_wrfp(i,j,kp) = t3_wrfs(i,j,nzw)

            if(icaller .eq. 1)then ! called from 'wrf2swim'
              qc_wrfp(i,j,kp) = qc_wrfs(i,j,nzw)
              qi_wrfp(i,j,kp) = qi_wrfs(i,j,nzw)
              qr_wrfp(i,j,kp) = qr_wrfs(i,j,nzw)
              qs_wrfp(i,j,kp) = qs_wrfs(i,j,nzw)
              if(istat_aod .eq. 1)then
                aod_wrfp(i,j,kp) = aod_wrfs(i,j,nzw)
              endif
            else
              u3_wrfp(i,j,kp) = u3_wrfs(i,j,nzw)
              v3_wrfp(i,j,kp) = v3_wrfs(i,j,nzw)
              om_wrfp(i,j,kp) = 0.
            endif

            dz = t3_wrfs(i,j,nzw) * rog * ALOG(pr_wrfs(i,j,nzw)/pr_laps(kp))
            ht_wrfp(i,j,kp) = ht_wrfs(i,j,nzw) + dz 
            ! Reduce moisture toward zero
            IF ( pr_laps(kp-1) .GT. pr_wrfs(i,j,nzw) ) THEN
              sh_wrfp(i,j,kp) = 0.5*sh_wrfs(i,j,nzw)
            ELSE
              sh_wrfp(i,j,kp) = 0.5*sh_wrfp(i,j,kp-1)
            ENDIF
          ENDIF

        ENDDO pressloop
      ENDDO
    ENDDO

    ! If we do smoothing, do it here
    
    ! Print some diagnostics
    PRINT *, "WRF Press data from center of WRF domain ",nzl
    if(icaller .eq. 1)then ! called from 'wrf2swim'
      print *, "KP  PRESS     HEIGHT   TEMP   SH        QC       QI      QR       QS"
      print *, "--- --------  -------  -----  -------   ------   ------  -----    -----"
    else
      print *, "KP  PRESS     HEIGHT   TEMP   SH       U       V      OM"
      print *, "--- --------  -------  -----  -------  ------  ------ -----------"
    endif
    DO k = 1,nzl
      if(icaller .eq. 1)then ! called from 'wrf2swim'
        PRINT ('(I3,1x,F8.1,2x,F7.0,2x,F5.1,2x,F7.5,2x,4F8.5,2x,F6.1,2x,F11.8)'), &
          k,pr_laps(k),ht_wrfp(icentw,jcentw,k), t3_wrfp(icentw,jcentw,k), &
          sh_wrfp(icentw,jcentw,k),qc_wrfp(icentw,jcentw,k),qi_wrfp(icentw,jcentw,k), &
          qr_wrfp(icentw,jcentw,k),qs_wrfp(icentw,jcentw,k)
      else
        PRINT ('(I3,1x,F8.1,2x,F7.0,2x,F5.1,2x,F7.5,2x,F6.1,2x,F6.1,2x,F11.8)'), &
          k,pr_laps(k),ht_wrfp(icentw,jcentw,k), t3_wrfp(icentw,jcentw,k), &
          sh_wrfp(icentw,jcentw,k),u3_wrfp(icentw,jcentw,k),v3_wrfp(icentw,jcentw,k),&
          om_wrfp(icentw,jcentw,k)
      endif
    ENDDO
    
    RETURN
  END SUBROUTINE vinterp_wrfarw2p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE hinterp_wrf2lga(nxl,nyl,nzl,istatus)

  IMPLICIT NONE
  INTEGER, INTENT(OUT)  :: nxl,nyl,nzl,istatus
  INTEGER               :: i,j,k
  REAL, ALLOCATABLE     :: xloc(:,:), yloc(:,:),dum2d(:,:)
  REAL, ALLOCATABLE     :: topo_wrfl(:,:)  ! WRF Topo interpolated to LAPS grid
  REAL                  :: ri,rj,dtopo,dz, wgt1, wgt2,lp
  istatus = 1 

  ALLOCATE(dum2d(nxl,nyl))
  ! First, generate xloc/yloc locations
  ALLOCATE(xloc(nxl,nyl))
  ALLOCATE(yloc(nxl,nyl))
  DO j = 1, nyl
    DO i = 1, nxl
      CALL latlon_to_ij(wrfgrid,lat(i,j),lon(i,j),ri,rj)
      xloc(i,j) = ri
      yloc(i,j) = rj
    ENDDO
  ENDDO

  ! Now, horizontally interpolate level by level
  DO k=1,nzl
    ! Height
    CALL interpolate_standard(nxw,nyw,ht_wrfp(:,:,k), &
                              nxl,nyl,xloc,yloc,METHOD_LINEAR, &
                              dum2d)
    ht(:,:,k) = dum2d


    ! Temp
    CALL interpolate_standard(nxw,nyw,t3_wrfp(:,:,k), &
                              nxl,nyl,xloc,yloc,METHOD_LINEAR, &
                              dum2d)
    t3(:,:,k) = dum2d

    ! SH
    CALL interpolate_standard(nxw,nyw,sh_wrfp(:,:,k), &
                              nxl,nyl,xloc,yloc,METHOD_LINEAR, &
                              dum2d)
    sh(:,:,k) = dum2d

    ! U3
    CALL interpolate_standard(nxw,nyw,u3_wrfp(:,:,k), &
                              nxl,nyl,xloc,yloc,METHOD_LINEAR, &
                              dum2d)
    u3(:,:,k) = dum2d

    ! v3
    CALL interpolate_standard(nxw,nyw,v3_wrfp(:,:,k), &
                              nxl,nyl,xloc,yloc,METHOD_LINEAR, &
                              dum2d)
    v3(:,:,k) = dum2d

    ! om
    CALL interpolate_standard(nxw,nyw,om_wrfp(:,:,k), &
                              nxl,nyl,xloc,yloc,METHOD_LINEAR, &
                              dum2d)
    om(:,:,k) = dum2d

  ENDDO

  PRINT *, "LGA Press data from center of LAPS domain"
  print *, "KP  PRESS     HEIGHT   TEMP   SH       U       V      OM"
  print *, "--- --------  -------  -----  -------  ------  ------ -----------"
  DO k = 1,nzl
    PRINT ('(I3,1x,F8.1,2x,F7.0,2x,F5.1,2x,F7.5,2x,F6.1,2x,F6.1,2x,F11.8)'), &
          k,pr_laps(k),ht(icentl,jcentl,k), t3(icentl,jcentl,k), &
          sh(icentl,jcentl,k),u3(icentl,jcentl,k),v3(icentl,jcentl,k),&
          om(icentl,jcentl,k)
  ENDDO

  DEALLOCATE(dum2d)

  ALLOCATE(topo_wrfl(nxl,nyl))
  ! Do the surface variables
  CALL interpolate_standard(nxw,nyw,usf_wrf, nxl,nyl,xloc,yloc, &
      METHOD_LINEAR,usf)
  CALL interpolate_standard(nxw,nyw,vsf_wrf, nxl,nyl,xloc,yloc, &
      METHOD_LINEAR,vsf)
  CALL interpolate_standard(nxw,nyw,tsf_wrf, nxl,nyl,xloc,yloc, &
      METHOD_LINEAR,tsf)
  CALL interpolate_standard(nxw,nyw,tsk_wrf, nxl,nyl,xloc,yloc, &
      METHOD_LINEAR,tsk) ! added tsk by Wei-Ting (130312)
  CALL interpolate_standard(nxw,nyw,dsf_wrf, nxl,nyl,xloc,yloc, &
      METHOD_LINEAR,dsf)
  CALL interpolate_standard(nxw,nyw,rsf_wrf, nxl,nyl,xloc,yloc, &
      METHOD_LINEAR,rsf)
  CALL interpolate_standard(nxw,nyw,psf_wrf, nxl,nyl,xloc,yloc, &
      METHOD_LINEAR,psf)
  CALL interpolate_standard(nxw,nyw,topo_wrf, nxl,nyl,xloc,yloc, &
      METHOD_LINEAR,topo_wrfl)
  CALL interpolate_standard(nxw,nyw,pcp_wrf, nxl,nyl,xloc,yloc, &
      METHOD_LINEAR,pcp) ! added pcp by Wei-Ting (130312)
  where ( pcp < 0 ) ; pcp = 0 ; endwhere ! keep pcp >= 0 added by Wei-Ting (130312)

  ! Adjust for terrain differences between WRF and LAPS
  PRINT *, "Adjusting WRF surface to LAPS surface"
  DO j = 1 , nyl
    DO i = 1 , nxl
      dtopo = topo_wrfl(i,j) - topo_laps(i,j)
      IF (ABS(dtopo) .GT. 10.) THEN
        IF (dtopo .GT. 0) THEN ! Move downward to LAPS level
          downloop:  DO k = nzl , 1, -1
            IF (ht(i,j,k) .LE. topo_laps(i,j)) THEN
              dz = topo_wrfl(i,j) - ht(i,j,k)
              wgt1 = (topo_wrfl(i,j) - topo_laps(i,j)) / dz
              wgt2 = 1.0 - wgt1
              lp = wgt1*ALOG(pr_laps(k)) + wgt2*ALOG(psf(i,j))
              psf(i,j) = EXP(lp)
              tsf(i,j) = wgt1*t3(i,j,k)  + wgt2*tsf(i,j)
              tsk(i,j) = wgt1*t3(i,j,k)  + wgt2*tsk(i,j) ! added tsk by Wei-Ting (130312)
              rsf(i,j) = wgt1*sh(i,j,k)  + wgt2*rsf(i,j)
              usf(i,j) = wgt1*u3(i,j,k)  + wgt2*usf(i,j)
              vsf(i,j) = wgt1*v3(i,j,k)  + wgt2*vsf(i,j)

              EXIT downloop
            ENDIF
          ENDDO downloop
        ELSE  ! Move upward to LAPS level
          uploop:  DO k = 1, nzl
            IF (ht(i,j,k) .GE. topo_laps(i,j)) THEN
              dz = ht(i,j,k) - topo_wrfl(i,j)
              wgt2 = (topo_laps(i,j) - topo_wrfl(i,j)) / dz
              wgt1 = 1.0 - wgt2
              psf(i,j) = wgt1 * psf(i,j) + wgt2*pr_laps(k)
              tsf(i,j) = wgt1 * tsf(i,j) + wgt2*t3(i,j,k)
              tsk(i,j) = wgt1 * tsk(i,j) + wgt2*t3(i,j,k) ! added tsk by Wei-Ting (130312)
              rsf(i,j) = wgt1 * rsf(i,j) + wgt2*sh(i,j,k)
              usf(i,j) = wgt1 * usf(i,j) + wgt2*u3(i,j,k)
              vsf(i,j) = wgt1 * vsf(i,j) + wgt2*v3(i,j,k) 
 
              EXIT uploop
            ENDIF
          ENDDO uploop
        ENDIF
      ENDIF
   
    ENDDO
  ENDDO
  DEALLOCATE(xloc,yloc,topo_wrfl)

  END SUBROUTINE hinterp_wrf2lga

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE make_derived_pressures(nxl,nyl,nzl)

    IMPLICIT NONE
    INTEGER :: i,j,k,nxl,nyl,nzl
    REAL    :: wgt1, wgt2,lp1,lp2,dz,lp
    DO j = 1, nyl 
      DO i =  1, nxl

        ! Sea-level pressure
        slploop: DO k = 1, nzl
          IF (ht(i,j,k) .GE. 0) THEN
            lp1 = ALOG(pr_laps(k-1))
            lp2 = ALOG(pr_laps(k))
            dz = ht(i,j,k) - ht(i,j,k-1)
            wgt1 = ht(i,j,k)/dz
            wgt2 = 1. - wgt1
            lp = wgt1 * lp1 + wgt2 * lp2
            slp(i,j) = EXP(lp) 
            EXIT slploop
          ENDIF
        ENDDO slploop

        ! LAPS reduced pressure
        redploop: DO k = 1, nzl
          IF (ht(i,j,k) .GE. redp_lvl) THEN
            lp1 = ALOG(pr_laps(k-1))
            lp2 = ALOG(pr_laps(k))
            dz = ht(i,j,k) - ht(i,j,k-1)
            wgt1 = (ht(i,j,k)-redp_lvl)/dz
            wgt2 = 1. - wgt1
            lp = wgt1 * lp1 + wgt2 * lp2
            p(i,j) = EXP(lp)
            EXIT redploop
          ENDIF
        ENDDO redploop
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE make_derived_pressures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE wrf_lga
