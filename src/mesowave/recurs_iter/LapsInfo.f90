SUBROUTINE LapsInfo

!==========================================================
!  This routine gets the necesary LAPS information.
!
!  HISTORY: MAY. 2004 by YUANFU XIE.
!==========================================================

  CALL get_laps_info(nx,ny,stanlat,stanlat2,stanlon, &
                     laps_cycle_time,badflag ,maxstations,maproj)

  ! Debug: PRINT*,'LAPS: ',nx,ny,stanlat,stanlat2,stanlon, &
  !                   laps_cycle_time,badflag ,maxstations,maproj(1:6)

  l(1:4) = (/mx,my,mt,mv/)

  ! Set number of fictitious grid points surrounding the grid:
  nfic = 20

  n(1) = nx+2*nfic
  n(2) = ny+2*nfic
  n(3) = ncycles
  n(4) = nvlaps

  PRINT*,'N = ',n
  PRINT*,'Nxy without fictitous: ',nx,ny

  ! QC: threshold value check:
  qc_cons(1) = 10.0		! T
  qc_cons(5) = 10.0		! DT
  qc_cons(2) = 5.0		! U
  qc_cons(3) = 5.0		! V
  qc_cons(4) = 400.0		! MSL P
  qc_cons(6) = 400.0		! RED P
  
  ! Check:
  IF ((n(1) .GT. mx) .OR. (n(2) .GT. my) .OR. &
      (n(3) .GT. mt) .OR. (n(4) .GT. mv)) THEN
     PRINT*,'LapsInfo: Analysis array is too small!'
     PRINT*,'X: ',n(1),mx
     PRINT*,'Y: ',n(2),my
     PRINT*,'T: ',n(3),mt
     PRINT*,'V: ',n(4),mv
     STOP
  ENDIF

  ! Read in recursive filter parameters:
  call get_directory('static', dir_s, len)
  name(1:len) = dir_s(1:len)
  name(len+1:len+12) = 'Namelist.txt'

  PRINT*,'Path to NAMELIST: ',name(1:len+12)

  OPEN(unit=13,file=name(1:len+12),form='formatted')
  DO i1=1,4
     READ(13,*)
  ENDDO
  DO i1=1,n(4)
     READ(13,*) al(1:3,i1),np(1:3,i1)
  ENDDO
  ! Number of minimization iterations:
  READ(13,*) maxitr

  ! Number of recursive filter iterations:
  READ(13,*) nrf(1:n(4))
  CLOSE(13)

END SUBROUTINE
