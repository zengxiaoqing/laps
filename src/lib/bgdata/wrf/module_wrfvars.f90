MODULE basicvars

!  Contains the arrays and routines needed for allocating
!  and populating the basic variables from WRF.

  USE constants
  USE wrf_netcdf
  USE wrfpost_config
  IMPLICIT NONE

  PUBLIC
  ! 3d
  REAL, ALLOCATABLE        :: p3 (:,:,:)
  REAL, ALLOCATABLE        :: t3 (:,:,:)
  REAL, ALLOCATABLE        :: z3 (:,:,:)
  REAL, ALLOCATABLE        :: z3f(:,:,:)
  REAL, ALLOCATABLE        :: qv3(:,:,:)
  REAL, ALLOCATABLE        :: rh3(:,:,:)
  REAL, ALLOCATABLE        :: u3 (:,:,:)
  REAL, ALLOCATABLE        :: v3 (:,:,:)
  REAL, ALLOCATABLE        :: w3 (:,:,:)
 
  ! 2d
  REAL, ALLOCATABLE        :: p_sfc(:,:)
  REAL, ALLOCATABLE        :: mslp(:,:)
  REAL, ALLOCATABLE        :: z_1000(:,:)
  REAL, ALLOCATABLE        :: tv_sl(:,:)
  REAL, ALLOCATABLE        :: t_2m (:,:)
  REAL, ALLOCATABLE        :: q_2m (:,:)
  REAL, ALLOCATABLE        :: td_2m(:,:)
  REAL, ALLOCATABLE        :: rh_2m(:,:)
  REAL, ALLOCATABLE        :: u_10m (:,:)
  REAL, ALLOCATABLE        :: v_10m (:,:)
  REAL, ALLOCATABLE        :: wdir_10m (:,:)
  REAL, ALLOCATABLE        :: wspd_10m (:,:)
  REAL, ALLOCATABLE        :: qpf_sum(:,:)
  REAL, ALLOCATABLE        :: qpf_inc(:,:)
  REAL, ALLOCATABLE        :: pcprate(:,:)
  INTEGER, ALLOCATABLE     :: conv_flag(:,:) 
  REAL, ALLOCATABLE        :: xlat(:,:),xlon(:,:),topo(:,:)
  REAL, ALLOCATABLE        :: landmask(:,:)
  REAL, ALLOCATABLE        :: t_skin(:,:)
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE wrfpost_get_basic()

  IMPLICIT NONE
  include 'netcdf.inc'
  REAL, ALLOCATABLE               :: d1(:,:,:)
  REAL, ALLOCATABLE               :: d2(:,:,:)
  REAL, ALLOCATABLE               :: d1f(:,:,:)
  REAL, ALLOCATABLE               :: d2f(:,:,:)
  REAL, ALLOCATABLE               :: dum2d (:,:)
  REAL, EXTERNAL                  :: relhum, dewpt2, wspd, wdir,mixsat
  INTEGER :: status, i, j, k, rcode
 
  ALLOCATE (d1 (nx,ny,nz) )
  ALLOCATE (d2 (nx,ny,nz) )
  ALLOCATE (d1f (nx,ny,nz+1) )
  ALLOCATE (d2f (nx,ny,nz+1) )

  IF (.NOT. ALLOCATED(p3)) ALLOCATE ( p3 (nx,ny,nz) )
  IF (.NOT. ALLOCATED(z3)) ALLOCATE ( z3 (nx,ny,nz) )
  IF (.NOT. ALLOCATED(z3f)) ALLOCATE ( z3f (nx,ny,nz+1) )
  IF (.NOT. ALLOCATED(t3)) ALLOCATE ( t3 (nx,ny,nz) )
  IF (.NOT. ALLOCATED(qv3)) ALLOCATE ( qv3(nx,ny,nz) )
  IF (.NOT. ALLOCATED(u3))ALLOCATE ( u3 (nx,ny,nz) )
  IF (.NOT. ALLOCATED(v3)) ALLOCATE ( v3 (nx,ny,nz) )
  IF (.NOT. ALLOCATED(w3)) ALLOCATE ( w3 (nx,ny,nz) )
  IF (.NOT. ALLOCATED(p_sfc)) ALLOCATE ( p_sfc (nx,ny) )
  IF (.NOT. ALLOCATED(t_2m)) ALLOCATE ( t_2m(nx,ny) )
  IF (.NOT. ALLOCATED(q_2m))ALLOCATE ( q_2m (nx,ny) )
  IF (.NOT. ALLOCATED(u_10m)) ALLOCATE ( u_10m (nx,ny) )
  IF (.NOT. ALLOCATED(v_10m)) ALLOCATE ( v_10m (nx,ny) )
  IF (.NOT. ALLOCATED(qpf_sum)) ALLOCATE ( qpf_sum (nx,ny) )
  IF (.NOT. ALLOCATED(pcprate)) ALLOCATE ( pcprate (nx,ny) )
  IF (.NOT. ALLOCATED(xlat)) ALLOCATE ( xlat (nx,ny) )
  IF (.NOT. ALLOCATED(xlon)) ALLOCATE ( xlon (nx,ny) )
  IF (.NOT. ALLOCATED(topo)) ALLOCATE ( topo (nx,ny) )
  IF (.NOT. ALLOCATED(t_skin)) ALLOCATE ( t_skin (nx,ny) )
  IF (.NOT. ALLOCATED(landmask)) ALLOCATE ( landmask (nx,ny) )
  IF (.NOT. ALLOCATED(conv_flag)) ALLOCATE(conv_flag (nx,ny) )
  conv_flag(:,:) = 0
  ! Get pressures
  CALL get_wrfnc_3d(ncfile,"PB","T",nx,ny,nz,1,d1,status)
  IF (status .GT. 0) THEN
    PRINT *, "+++CRITICAL:  Could not get base pressure!"
    CALL ABORT
  ENDIF
  p3 = d1

  CALL get_wrfnc_3d(ncfile,"P","T",nx,ny,nz,1,d1,status)
  IF (status .GT. 0) THEN
    PRINT *, "+++CRITICAL:  Could not get pert pressure!"
    CALL ABORT
  ENDIF
  p3 = p3 + d1

  
  ! Get heights
  CALL get_wrfnc_3d(ncfile,"PHB","T",nx,ny,nz+1,1,d1f,status)
  IF (status.NE.0) THEN
    PRINT *, 'Could not properly obtain WRF base-state geopotential.'
    CALL ABORT
  ENDIF
  d2f = d1f
  CALL get_wrfnc_3d(ncfile,"PH","T",nx,ny,nz+1,1,d1f,status)
  IF (status.NE.0) THEN
    PRINT *, 'Could not properly obtain WRF geopotential.'
    CALL ABORT
  ENDIF
  d2f = (d2f + d1f) / grav
  z3f = d2f
  DO k = 1,nz
    z3(:,:,k) = 0.5 * (z3f(:,:,k) + z3f(:,:,k+1))
  ENDDO

  ! Get theta and convert to temperature
  CALL get_wrfnc_3d(ncfile, "T","T",nx,ny,nz,1,d1,status)
  IF (status.NE.0) THEN
    PRINT *, 'Could not properly obtain WRF perturbation theta.'
    CALL ABORT
  ENDIF
  d1 = d1 + 300.
  DO k = 1, nz
    DO j = 1, ny
      DO i = 1,nx
        t3(i,j,k) = d1(i,j,k)/ ((100000./p3(i,j,k))**kappa)
      ENDDO
    ENDDO
  ENDDO
 
   ! Get Q on sigma
   CALL get_wrfnc_3d(ncfile, "QVAPOR","T",nx,ny,nz,1,qv3,status)
   IF (status.NE.0) THEN
     PRINT *, 'Could not properly obtain WRF mixing ratio.'
     CALL ABORT
   ENDIF

   ! Get U on sigma
   CALL get_wrfnc_3d(ncfile, "U","T",nx,ny,nz,1,u3,status)
   IF (status.NE.0) THEN
     PRINT *, 'Could not properly obtain WRF U-comp.'
     CALL ABORT
   ENDIF

   ! Get U on sigma
   CALL get_wrfnc_3d(ncfile, "V","T",nx,ny,nz,1,v3,status)
   IF (status.NE.0) THEN
     PRINT *, 'Could not properly obtain WRF V-comp.'
     CALL ABORT
   ENDIF

   ! Get W on sigma
   CALL get_wrfnc_3d(ncfile, "W","T",nx,ny,nz+1,1,d1f,status)
   IF (status.NE.0) THEN
     PRINT *, 'Could not properly obtain WRF W-comp.'
     CALL ABORT
   ENDIF
   DO k = 1,nz
     w3(:,:,k) = 0.5*(d1f(:,:,k)+d1f(:,:,k+1))
   ENDDO

  
  DEALLOCATE(d1)
  DEALLOCATE(d2)
  DEALLOCATE(d1f)
  DEALLOCATE(d2f)

  ! Get the 2D stuff
  
  ! Pressure at surface
  p_sfc(:,:) = 0.
  CALL get_wrfnc_2d(ncfile, "PSFC","T",nx,ny,1,p_sfc,status)
  IF ((status .NE. 0).OR.(MAXVAL(p_sfc) .LT. 10000.))THEN
    PRINT *, "Could not get PSFC, using lowest sigma level"
    p_sfc = p3(:,:,1)
  ENDIF
  
  ! 2m Temp
  t_2m(:,:) = 0.
  CALL get_wrfnc_2d(ncfile, "T2","T",nx,ny,1,t_2m,status)
  IF ((status .NE. 0).OR.(MAXVAL(t_2m) .LT. 100.))THEN
    PRINT *, "Could not get T2, using lowest sigma level"
    t_2m = t3(:,:,1)
  ENDIF
  
  ! qvapor at 2m
  q_2m(:,:) = 0.
  CALL get_wrfnc_2d(ncfile, "Q2","T",nx,ny,1,q_2m,status)
  IF ((status .NE. 0).OR.(MAXVAL(q_2m) .LT. 0.0001))THEN
    PRINT *, "Could not get Q2, using lowest sigma level"
    q_2m = qv3(:,:,1)
  ENDIF
  ! Because 2m qv and T are derived from the PBL scheme and
  ! the WRF is apparently not checking for saturation, clean this
  ! up now
  DO j = 1, ny
    DO i= 1, nx
       q_2m(i,j) = MIN(q_2m(i,j),mixsat(t_2m(i,j),p_sfc(i,j)))
    ENDDO
  ENDDO
  
  ! 10m wind
  !u_10m(:,:) = 0. 
  !CALL get_wrfnc_2d(ncfile, "U10","T",nx,ny,1,u_10m,status) 
  !IF ((status .NE. 0).OR.(MAXVAL(u_10m) .EQ. 0.))THEN
   ! PRINT *, "Could not get U10, using lowest sigma level"
    u_10m = u3(:,:,1)
  !ENDIF
  
  !v_10m(:,:) = 0.
  !CALL get_wrfnc_2d(ncfile, "V10","T",nx,ny,1,v_10m,status)
  !IF ((status .NE. 0).OR.(MAXVAL(v_10m) .EQ. 0.))THEN
  !  PRINT *, "Could not get V10, using lowest sigma level"
    v_10m = v3(:,:,1)
  !ENDIF  
  
  ! Precip
  ALLOCATE(dum2d (nx,ny))
  qpf_sum(:,:) = 0.
  CALL get_wrfnc_2d(ncfile, "RAINNC","T",nx,ny,1,dum2d,status)
  IF (status .NE. 0)THEN
    PRINT *, "Could not get RAINNC, assuming 0."
  ELSE
    qpf_sum = dum2d
    CALL get_wrfnc_2d(ncfile, "RAINC","T",nx,ny,1,dum2d,status)
    IF (status .NE. 0)THEN
       PRINT *, "Could not get RAINC, assuming 0."
    ELSE
      ! Set convective flag for any point that received parameterized precip
       DO j = 1, ny
         DO i = 1, nx
           IF(dum2d(i,j) .GT. 0) conv_flag(i,j) = 1
            ! Zero out small values
           !IF (dum2d(i,j) .LT. 1e-4) dum2d(i,j) = 0
          ENDDO
        ENDDO
       qpf_sum = qpf_sum + dum2d
    ENDIF
  ENDIF
  WHERE(qpf_sum .LT. 1.e-4) qpf_sum = 0.
  ! Lets get incremental precipitation
  dum2d(:,:) = 0.
  IF (.NOT. ALLOCATED (qpf_inc)) ALLOCATE(qpf_inc(nx,ny))
  IF (have_prev) THEN
    IF (pretime_str .NE. reftime_str) THEN
      ! There is a previous output time available
      ! and it is not the 0-h tau, so let's compute
      ! get the qpf from the previous time and subtract it out
      CALL open_wrfnc(filename_prev, ncfile_prev, status)
      IF (status .EQ. 0) THEN
        CALL get_wrfnc_2d(ncfile_prev, "RAINNC","T",nx,ny,1,dum2d,status)
	IF (status .NE. 0) dum2d(:,:) = 0.
	qpf_inc = dum2d
	CALL get_wrfnc_2d(ncfile_prev, "RAINC","T",nx,ny,1,dum2d,status)
	IF (status .EQ. 0) THEN 
          ! If we had convective precip at any grid point in the previous time
          ! go ahead and set the flag for convection
          DO j = 1, ny
            DO i = 1, nx
              IF(dum2d(i,j) .GT. 0) conv_flag(i,j) = 1
              ! Zero out small values
              !IF (dum2d(i,j) .LT. 1e-4) dum2d(i,j)= 0.
            ENDDO
          ENDDO
          qpf_inc = qpf_inc + dum2d
        ENDIF
        WHERE(qpf_inc .LT. 1.e-4) qpf_inc = 0.
        qpf_inc = qpf_sum - qpf_inc
        CALL close_wrfnc(ncfile_prev)
      ELSE 
        PRINT *, "ERROR opening previous file"
        qpf_inc = 0.
      ENDIF
    ELSE
      qpf_inc = qpf_sum 
    ENDIF
  ELSE
    qpf_inc(:,:) = 0.
  ENDIF
  
  PRINT *, "   Min/Max QPF total: ", MINVAL(qpf_sum),MAXVAL(qpf_sum)
  PRINT *, "   Min/Max QPF incre: ", MINVAL(qpf_inc),MAXVAL(qpf_inc)

  ! Preciptation rate (get the timestep resolved and convective
  ! precip for this period, divide by the timestep (Seconds) and
  ! multiply by 3600 to get mm per hour.  If we cannot obtain
  ! it, then deallocate it.
  PRINT *, '  -- Precip Rate'
  ALLOCATE(dum2d(nx,ny))
  CALL get_wrfnc_2d(ncfile, "RAINNCV", "T",nx,ny,1,dum2d,status)
  IF (status.GT.0) THEN
    PRINT *, "Did not get RAINNCV"
    pcprate = 0.
  ELSE
    WHERE(dum2d .LT. 1.e-07) dum2d = 0.
    pcprate = dum2d
    CALL get_wrfnc_2d(ncfile, "RAINCV", "T",nx,ny,1,dum2d,status)
    IF (status .EQ. 0) THEN
      WHERE(dum2d .LT. 1.e-06) dum2d =0.
      pcprate = pcprate + dum2d
      DO j = 1, ny
        DO i = 1, nx
          IF(dum2d(i,j) .GT. 0) conv_flag(i,j) = 1
        ENDDO
      ENDDO
    ENDIF
  ENDIF
  rcode = NF_GET_ATT_REAL(ncfile,0,"DT",dt)
  IF (rcode .EQ. NF_NOERR) THEN
    pcprate = pcprate / dt   ! Converts to mm per sec
    WHERE (pcprate .LT. 1.e-8) pcprate = 0.
  ELSE
    print *, "Could not get timestep to compute pcprate!"
    DEALLOCATE(pcprate)
  ENDIF
  DEALLOCATE(dum2d)

  CALL get_wrfnc_2d(ncfile, "XLAT","T",nx,ny,1,xlat,status)
  CALL get_wrfnc_2d(ncfile, "XLONG","T",nx,ny,1,xlon,status)
  CALL get_wrfnc_2d(ncfile, "HGT","T",nx,ny,1,topo,status)
  CALL get_wrfnc_2d(ncfile, "LANDMASK","T",nx,ny,1,landmask,status)
  CALL get_Wrfnc_2d(ncfile, "TSK", "T", nx,ny,1,t_skin,status)


  ! Compute a few more things
  ! RH and dewpoint at surface are nice to have
  PRINT *, "  Computing Td_2m and RH_2m"
  ALLOCATE(rh_2m(nx,ny))
  ALLOCATE(td_2m(nx,ny))
  ALLOCATE(rh3(nx,ny,nz))
 
  DO k= 1, nz
    DO j = 1, ny
      DO i = 1, nx 
        rh3(i,j,k) = relhum(t3(i,j,k),qv3(i,j,k),p3(i,j,k)) * 100.
      ENDDO
    ENDDO
  ENDDO
  DO j = 1, ny
    DO i = 1, nx
      rh_2m(i,j) = relhum(t_2m(i,j),q_2m(i,j),p_sfc(i,j))
      IF (rh_2m(i,j) .GT. 1.0) THEN
        print *,"i/j/t/p/q/rh = ",i,j,t_2m(i,j),p_sfc(i,j),q_2m(i,j),rh_2m(i,j)
      ENDIF
      td_2m(i,j) = dewpt2(t_2m(i,j),rh_2m(i,j))
      rh_2m(i,j) = rh_2m(i,j) * 100.
    ENDDO
  ENDDO
  PRINT *, '    Min/Max Td_2m = ', MINVAL(td_2m),MAXVAL(td_2m)
  PRINT *, '    Min/Max RH_2m = ', MINVAL(rh_2m),MAXVAL(rh_2m)

  ! Lets also get a surface true wind direction and speed
  ALLOCATE(wdir_10m(nx,ny))
  ALLOCATE(wspd_10m(nx,ny))
  DO j = 1, ny
    DO i = 1, nx
      wdir_10m(i,j) = wdir(u_10m(i,j),v_10m(i,j),xlon(i,j),grid%stdlon, &
                     grid%cone)
      wspd_10m(i,j) = wspd(u_10m(i,j),v_10m(i,j))
    ENDDO
  ENDDO

  ! Mean sea-level pressure
  ALLOCATE(mslp(nx,ny))
  ALLOCATE(tv_sl(nx,ny))
  ALLOCATE(z_1000(nx,ny))
  PRINT *, "WRFPOST_BASICVARS:  Calling derive_mslp"
  CALL derive_mslp(nx,ny,nz,t3,qv3,p3,p_sfc,topo,mslp,tv_sl,z_1000)

  RETURN

END SUBROUTINE wrfpost_get_basic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE basicvars
