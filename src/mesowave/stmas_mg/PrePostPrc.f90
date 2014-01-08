!dis Forecast Systems Laboratory
!dis NOAA/OAR/ERL/FSL
!dis 325 Broadway
!dis Boulder, CO 80303
!dis
!dis Forecast Research Division
!dis Local Analysis and Prediction Branch
!dis LAPS
!dis
!dis This software and its documentation are in the public domain and
!dis are furnished "as is." The United States government, its
!dis instrumentalities, officers, employees, and agents make no
!dis warranty, express or implied, as to the usefulness of the software
!dis and documentation for any purpose. They assume no responsibility
!dis (1) for the use of the software and documentation; or (2) to provide
!dis technical support to users.
!dis
!dis Permission to use, copy, modify, and distribute this software is
!dis hereby granted, provided that the entire disclaimer notice appears
!dis in all copies. All modifications to this software must be clearly
!dis documented, and are solely the responsibility of the agent making
!dis the modifications. If significant modifications or enhancements
!dis are made to this software, the FSL Software Policy Manager
!dis (softwaremgr@fsl.noaa.gov) should be notified.
!dis

MODULE PrePostPrc

!==========================================================
! This module sets up STMAS pre/post processes.
!
! HISTORY:
! Creation: YUANFU XIE 8-2005
!==========================================================

  USE Definition

CONTAINS

SUBROUTINE PrPstNLs

!==========================================================
! This routine reads in namelists for STMAS multigrid.
!
! HISTORY:
! Creation: YUANFU XIE 8-2005
!==========================================================

  IMPLICIT NONE

  CHARACTER*200 :: dir
  INTEGER :: nam,ios

  ! Get namelist for STMAS:
  CALL get_directory('static',dirstc,dirlen)
  dir = dirstc(1:dirlen)//'stmas_mg.nl'
  OPEN(unit=11,file=dir(1:dirlen+12),form='formatted')
  READ(11,NML=STMAS,IOSTAT=ios)
  CLOSE(11)
  ! Check the numtmf consistent to numgrd
  IF (numtmf .NE. numgrd(3)) THEN
    PRINT*,'STMAS>PrPstNLs: Error in stmas_mg.nl, numtmf does not match numgrd'
    STOP
  ENDIF

  ! Check the limit for number of variables to be analyzed:
  IF (numvar .GT. MAXVAR) THEN
    WRITE(*,*) 'STMAS>PrPstNLs: Error: too many variables!'
    WRITE(*,*) 'STMAS>PrPstNLs: Increase MAXVAR and rerun!'
    STOP
  ENDIF

  ! Get names of analyzing variables:
  dir = dirstc(1:dirlen)//'stmas_mg.vr'
  OPEN(unit=11,file=dir(1:dirlen+12),form='formatted')
  DO nam=1,numvar
    READ(11,*) varnam(nam),thresh(nam),needbk(nam),bounds(nam),radius(nam),pnlt_v(nam),lndsea(nam),slevel(nam)
  ENDDO
  CLOSE(11)

END SUBROUTINE PrPstNLs

SUBROUTINE PrPstLSX

!==========================================================
! This routine writes out the analyses into LSX gridded
! data in NetCDF format.
!
! HISTORY:
! Creation: YUANFU XIE 6-2005
! Modified: YUANFU XIE 6-2006: Write out more time frame.
! Modified: YUANFU XIE 6-2006: Add sfc pressure and change
!			       the equations for theta and
!			       thetae from REDP to SFCP.
!==========================================================

  USE Definition

  IMPLICIT NONE

  ! Local variables:
  CHARACTER :: ext*3,dir*200
  CHARACTER*125 :: cmt(LSXVAR)
  CHARACTER*3 :: vnm(LSXVAR)		! Variable names
  CHARACTER*3 :: vun(LSXVAR)		! Units
  CHARACTER*3 :: crd(LSXVAR)		! Coordinates

  INTEGER :: i,j,ncm
  INTEGER :: ngd(2)	! Actual grid points
  INTEGER :: lvl(LSXVAR) ! Number of levels of each field
  INTEGER :: itm	! Time frame to write out
  INTEGER :: i4t	! i4 time corresponding to itm
  INTEGER :: idx(MAXVAR)! Indices for derived variables;
  INTEGER :: iwv	! Index of V; 
  INTEGER :: nvr,mvr	! Number of variables;
  INTEGER :: len
  INTEGER :: sts	! Return status

  REAL, EXTERNAL :: EPT
  REAL :: tmp,dew,td_c,pr_m
  REAL :: gdx(numgrd(1),numgrd(2))	! Grid spacing (X)
  REAL :: gdy(numgrd(1),numgrd(2))	! Grid spacing (Y)
  REAL :: gdp(numgrd(1),numgrd(2))	! Pressure in mb
  REAL :: gdt(numgrd(1),numgrd(2))	! Theta
  REAL :: dum(numgrd(1),numgrd(2))	! Unused value
  REAL :: dmm(numgrd(1),numgrd(2))	! Unused value
  REAL :: mrc(numgrd(1),numgrd(2))
  REAL :: dat(numgrd(1)-2*numfic(1),numgrd(2)-2*numfic(2),LSXVAR)

  ! added by min-ken hsieh
  ! declare L1S array to store 1 hr and 24 hr precip.
  ! then output to L1S file via write_laps_data
  REAL :: L1S(numgrd(1)-2*numfic(1),numgrd(2)-2*numfic(2),2)
  CHARACTER*125 :: L1S_cmt(2)
  CHARACTER*3 :: L1S_vnm(2)		! Variable names
  CHARACTER*3 :: L1S_vun(2)		! Units
  CHARACTER*3 :: L1S_crd(2)		! Coordinates
  CHARACTER*3 :: LWM_vnm(2)		! Variable names

  INTEGER :: L1S_lvl(2)			! Number of levels of each field

  ! Mixing ratio function from TD and P:
  REAL :: WMR

  ! Time frame to write out:
  !DO itm = numgrd(3)-numfic(3),max0(1,numgrd(3)-numfic(3)-2),-1	! Time frame
  ! Yuanfu: change output to current time frame only for forecast:
  DO itm = numgrd(3)-numfic(3),max0(1,numgrd(3)-numfic(3)),-1	! Time frame

  i4t = i4time-(numgrd(3)-numfic(3)-itm)*lapsdt	! i4time corresponding to itm
  mvr = LSXVAR
  lvl = 0
  crd = 'AGL'

  ! Number of gridpoints without fictitous points:
  ngd = numgrd(1:2)-2*numfic(1:2)
  gdx = phydxy(1)
  gdy = phydxy(2)

  ! Parsing the variable names:
  nvr = 0
  DO i=1,numvar

    SELECT CASE (varnam(i))
    CASE ("TEMP")	! Temperature
      nvr = nvr+1
      IF (nvr .GT. LSXVAR) THEN
	WRITE(*,2)
        STOP
      ENDIF
      vnm(nvr) = 'T  '
      vun(nvr) = 'K  '
      cmt(nvr) = 'Surface temperature'
      dat(1:ngd(1),1:ngd(2),nvr) = &
	analys(numfic(1)+1:numgrd(1)-numfic(1), &
	       numfic(2)+1:numgrd(2)-numfic(2),itm,i)
    CASE ("VISB")	! Visibility
      nvr = nvr+1
      IF (nvr .GT. LSXVAR) THEN
	WRITE(*,2)
        STOP
      ENDIF
      vnm(nvr) = 'VIS'
      vun(nvr) = 'M  '
      cmt(nvr) = 'Visibilty'
      dat(1:ngd(1),1:ngd(2),nvr) = &
	analys(numfic(1)+1:numgrd(1)-numfic(1), &
	       numfic(2)+1:numgrd(2)-numfic(2),itm,i)
    CASE ("CEIL")	! Cloud ceiling
      nvr = nvr+1
      IF (nvr .GT. LSXVAR) THEN
	WRITE(*,2)
        STOP
      ENDIF
      vnm(nvr) = 'CC'
      vun(nvr) = 'M  '
      cmt(nvr) = 'Cloud ceiling'
      dat(1:ngd(1),1:ngd(2),nvr) = &
	analys(numfic(1)+1:numgrd(1)-numfic(1), &
	       numfic(2)+1:numgrd(2)-numfic(2),itm,i)
    CASE ("DEWP")	! Dew point temperature
      nvr = nvr+1
      IF (nvr .GT. LSXVAR) THEN
	WRITE(*,2)
        STOP
      ENDIF
      vnm(nvr) = 'TD '
      vun(nvr) = 'K  '
      cmt(nvr) = 'Dewpoint'
      dat(1:ngd(1),1:ngd(2),nvr) = &
	analys(numfic(1)+1:numgrd(1)-numfic(1), &
	       numfic(2)+1:numgrd(2)-numfic(2),itm,i)
    CASE ("TGD ")	! Ground temperature
      nvr = nvr+1
      IF (nvr .GT. LSXVAR) THEN
	WRITE(*,2)
        STOP
      ENDIF
      vnm(nvr) = 'TGD'
      vun(nvr) = 'K  '
      cmt(nvr) = 'Ground Temperature'
      dat(1:ngd(1),1:ngd(2),nvr) = &
	analys(numfic(1)+1:numgrd(1)-numfic(1), &
	       numfic(2)+1:numgrd(2)-numfic(2),itm,i)
    CASE ("MSLP")	! Mean Sea Level Pressure
      nvr = nvr+1
      IF (nvr .GT. LSXVAR) THEN
	WRITE(*,2)
        STOP
      ENDIF
      vnm(nvr) = 'MSL'
      vun(nvr) = 'PA '
      cmt(nvr) = 'MSL Pressure'
      dat(1:ngd(1),1:ngd(2),nvr) = &
	analys(numfic(1)+1:numgrd(1)-numfic(1), &
	       numfic(2)+1:numgrd(2)-numfic(2),itm,i)
    CASE ("REDP")	! Reduced pressure
      IF (nvr+2 .GT. LSXVAR) THEN
	WRITE(*,2)
        STOP
      ENDIF
      nvr = nvr+1
      vnm(nvr) = 'P  '
      vun(nvr) = 'PA '
      write(cmt(nvr),11) int(rdplvl),' M REDUCED PRESSURE'
 11   format(i5,a19)
      ! cmt(nvr) = '0  M REDUCED PRESSURE'
      dat(1:ngd(1),1:ngd(2),nvr) = &
	analys(numfic(1)+1:numgrd(1)-numfic(1), &
	       numfic(2)+1:numgrd(2)-numfic(2),itm,i)
      IF (press_pert .EQ. 1) THEN
        nvr = nvr+1
        vnm(nvr) = 'PP '
        vun(nvr) = 'PA '
        cmt(nvr) = '0  M REDUCED PRESSURE CHANGE'
        CALL PresChng(dat(1,1,nvr-1),ngd,dat(1,1,nvr))
      ENDIF
    CASE ("SFCP")	! Surface pressure
      IF (nvr+1 .GT. LSXVAR) THEN
	WRITE(*,2)
        STOP
      ENDIF
      nvr = nvr+1
      vnm(nvr) = 'PS '
      vun(nvr) = 'PA '
      cmt(nvr) = 'SURFACE PRESSURE'
      dat(1:ngd(1),1:ngd(2),nvr) = &
	analys(numfic(1)+1:numgrd(1)-numfic(1), &
	       numfic(2)+1:numgrd(2)-numfic(2),itm,i)
    CASE ("WNDU") 	! U wind binded with WNDV
      iwv = 0
      DO j=1,numvar
        IF (varnam(j) .EQ. 'WNDV') iwv = j
      ENDDO
      IF (iwv .EQ. 0) THEN
	WRITE(*,3)
        STOP
      ENDIF
      IF (nvr+4 .GT. LSXVAR) THEN
	WRITE(*,2)
        STOP
      ENDIF
      nvr = nvr+1
      vnm(nvr) = 'U  '
      vun(nvr) = 'M/S'
      cmt(nvr) = 'U component'
      dat(1:ngd(1),1:ngd(2),nvr) = &
	analys(numfic(1)+1:numgrd(1)-numfic(1), &
	       numfic(2)+1:numgrd(2)-numfic(2),itm,i)
      nvr = nvr+1
      vnm(nvr) = 'V  '
      vun(nvr) = 'M/S'
      cmt(nvr) = 'V component'
      dat(1:ngd(1),1:ngd(2),nvr) = &
	analys(numfic(1)+1:numgrd(1)-numfic(1), &
	       numfic(2)+1:numgrd(2)-numfic(2),itm,iwv)

      ! Write surface wind to lwm for plot purpose:
      ! LWM_vnm(1) = 'SU '
      ! LWM_vnm(2) = 'SV '
      ! CALL PUT_LAPS_MULTI_2D(i4t,'lwm',LWM_vnm,vun(nvr-1),cmt(nvr-1), &
      !                        dat(1,1,nvr-1),ngd(1),ngd(2),2,sts)

      nvr = nvr+1
      vnm(nvr) = 'VOR'
      vun(nvr) = '/S '
      cmt(nvr) = 'Vorticity'
      ! Interior:
      dat(2:ngd(1)-1,2:ngd(2)-1,nvr) = 0.5/phydxy(2)*( &
	analys(numfic(1)+2:numgrd(1)-numfic(1)-1, &
	       numfic(2)+1:numgrd(2)-numfic(2)-2,itm,i)- &
	analys(numfic(1)+2:numgrd(1)-numfic(1)-1, &
	       numfic(2)+3:numgrd(2)-numfic(2)  ,itm,i))+ &
		 		       0.5/phydxy(1)*( &
	analys(numfic(1)+3:numgrd(1)-numfic(1)  , &
	       numfic(2)+2:numgrd(2)-numfic(2)-1,itm,iwv)- &
	analys(numfic(1)+1:numgrd(1)-numfic(1)-2, &
	       numfic(2)+2:numgrd(2)-numfic(2)-1,itm,iwv))
      ! Extrapolate boundary values:
      CALL Extraplt(dat(1,1,nvr),ngd)
      nvr = nvr+1
      vnm(nvr) = 'DIV'
      vun(nvr) = '/S '
      cmt(nvr) = 'Divergence'
      ! Interior:
      dat(2:ngd(1)-1,2:ngd(2)-1,nvr) = -0.5/phydxy(2)*( &
	analys(numfic(1)+2:numgrd(1)-numfic(1)-1, &
	       numfic(2)+1:numgrd(2)-numfic(2)-2,itm,iwv)- &
	analys(numfic(1)+2:numgrd(1)-numfic(1)-1, &
	       numfic(2)+3:numgrd(2)-numfic(2)  ,itm,iwv))+ &
		 		       0.5/phydxy(1)*( &
	analys(numfic(1)+3:numgrd(1)-numfic(1)  , &
	       numfic(2)+2:numgrd(2)-numfic(2)-1,itm,i)- &
	analys(numfic(1)+1:numgrd(1)-numfic(1)-2, &
	       numfic(2)+2:numgrd(2)-numfic(2)-1,itm,i))
      ! Extrapolate boundary values:
      CALL Extraplt(dat(1,1,nvr),ngd)
    CASE ("WNDV")	! V wind dealt with U
    CASE ("PCP1")	! 1 hour precip accumulation		added by min-ken hsieh
      ! make a copy in L1S array
      L1S(1:ngd(1),1:ngd(2),1) = &
        analys(numfic(1)+1:numgrd(1)-numfic(1), &
               numfic(2)+1:numgrd(2)-numfic(2),itm,i)

    CASE ("PC24")	! 24 hour precip accumulation		added by min-ken hsieh
      ! make a copy in L1S array
      L1S(1:ngd(1),1:ngd(2),2) = &
        analys(numfic(1)+1:numgrd(1)-numfic(1), &
               numfic(2)+1:numgrd(2)-numfic(2),itm,i)
    CASE ("PCP3")
    CASE ("PCP6")
    CASE DEFAULT
      WRITE(*,1) varnam(i)
    END SELECT
  ENDDO

1 FORMAT('PrPstPrc>PrPstLSX: no such variable to write: ',a4)
2 FORMAT('PrPstPrc>PrPstLSX: too much variables to write; aborts')
3 FORMAT('PrPstPrc>PrPstLSX: V wind is missing')

  !======================
  ! Other derived fields:
  !======================
    
  ! Theta and Theta_e:	potential/equivalent potential temperature:
  ! Search necessary basic variables:
  ncm = 0
  idx = 0
  DO j=1,numvar
    IF (varnam(j) .EQ. "TEMP") THEN
      ncm = ncm+1
      idx(1) = j
    ENDIF
    IF (varnam(j) .EQ. "DEWP") THEN
      ncm = ncm+1
      idx(2) = j
    ENDIF
    IF (varnam(j) .EQ. "SFCP") THEN
      ncm = ncm+1
      idx(3) = j
    ENDIF
  ENDDO
  ! Found necessary variables for potential temperature:
  IF ((idx(1) .NE. 0) .AND. (idx(3) .NE. 0)) THEN
    nvr = nvr+1
    IF (nvr .GT. LSXVAR) THEN
      WRITE(*,2)
      STOP
    ENDIF
    vnm(nvr) = 'TH '
    vun(nvr) = 'K  '
    cmt(nvr) = 'Potential temperature'
    ! Theta:
    gdt(1:numgrd(1),1:numgrd(2)) = &
      analys(1:numgrd(1),1:numgrd(2),itm,idx(1))* &
	(1.0e5/analys(1:numgrd(1),1:numgrd(2),itm,idx(3))) &
	**(gascnt/spheat)
    ! Remove fictitous point and convert to r4:
    dat(1:ngd(1),1:ngd(2),nvr) = &
      gdt(numfic(1)+1:numgrd(1)-numfic(1), &
	  numfic(2)+1:numgrd(2)-numfic(2))
  ENDIF
  ! Found necessary variables for mixing ratio:
  IF ((idx(2) .NE. 0) .AND. (idx(3) .NE. 0)) THEN
    nvr = nvr+1
    IF (nvr .GT. LSXVAR) THEN
      WRITE(*,2)
      STOP
    ENDIF
    vnm(nvr) = 'MR '
    vun(nvr) = 'G/KG'
    cmt(nvr) = 'mixing ratio'
    ! Mixing ratio from td and p:
    DO j=1,ngd(2)
      DO i=1,ngd(1)
        td_c = analys(i,j,itm,idx(2))-temp_0	! Dewpoint in celsius
        pr_m = analys(i,j,itm,idx(3))/100.0	! Pressure in mb
        dat(i,j,nvr) = WMR(pr_m,td_c)
      ENDDO
    ENDDO
  ENDIF
  ! Found necessary variables for relative humidity:
  IF ((idx(1) .NE. 0) .AND. (idx(2) .NE. 0)) THEN
    nvr = nvr+1
    IF (nvr .GT. LSXVAR) THEN
      WRITE(*,2)
      STOP
    ENDIF
    vnm(nvr) = 'RH '
    vun(nvr) = 'M  '
    cmt(nvr) = 'relative humidity'
    ! RH from T and td: HUM returns RH in [0,1]
    CALL HUM(analys(1,1,itm,idx(1)),analys(1,1,itm,idx(2)), &
              dat(1,1,nvr),numgrd(1),numgrd(2),dum,dmm)
    dat(1:ngd(1),1:ngd(2),nvr) = dat(1:ngd(1),1:ngd(2),nvr)*100		! RH in [0,100]
  ENDIF
  ! Found necessary variables for equivalent potential temperature:
  IF (ncm .EQ. 3) THEN
    nvr = nvr+1
    IF (nvr .GT. LSXVAR) THEN
      WRITE(*,2)
      STOP
    ENDIF
    vnm(nvr) = 'THE'
    vun(nvr) = 'K  '
    cmt(nvr) = 'Equivalent potential temperature'
    ! Use LAPS function EPT to compute theta_e:
    gdp = analys(1:numgrd(1),1:numgrd(2),itm,idx(3))/mb2pas! Pressure(mb)
    DO j=1,ngd(2)
      DO i=1,ngd(1)
	tmp = analys(numfic(1)+i,numfic(2)+j,itm,idx(1))-temp_0
	dew = analys(numfic(1)+i,numfic(2)+j,itm,idx(2))-temp_0
        dat(i,j,nvr) = EPT(tmp,dew,gdp(i+numfic(1),j+numfic(2)))+temp_0
      ENDDO
    ENDDO
  ENDIF

  ! Moisture convergence:
  ! Search necessary basic variables:
  ncm = 0
  idx = 0
  DO j=1,numvar
    IF (varnam(j) .EQ. "TEMP") THEN
      ncm = ncm+1
      idx(1) = j
    ENDIF
    IF (varnam(j) .EQ. "DEWP") THEN
      ncm = ncm+1
      idx(2) = j
    ENDIF
    IF (varnam(j) .EQ. "REDP") THEN
      ncm = ncm+1
      idx(3) = j
    ENDIF
    IF (varnam(j) .EQ. "WNDU") THEN
      ncm = ncm+1
      idx(4) = j
    ENDIF
    IF (varnam(j) .EQ. "WNDV") THEN
      ncm = ncm+1
      idx(5) = j
    ENDIF
  ENDDO
  ! Found necessary variables:
  IF (ncm .EQ. 5) THEN
    nvr = nvr+1
    IF (nvr .GT. LSXVAR) THEN
      WRITE(*,2)
      STOP
    ENDIF
    vnm(nvr) = 'MRC'
    vun(nvr) = '/S '
    cmt(nvr) = 'Moisture convergence'
    CALL MESO_ANL(analys(1,1,itm,idx(4)),analys(1,1,itm,idx(5)), &
      gdp,analys(1,1,itm,idx(1)),analys(1,1,itm,idx(2)), &
      gdt,gdx,gdy,dum,mrc,dum,dum,dum,numgrd(1),numgrd(2))
    dat(1:ngd(1),1:ngd(2),nvr) = &
      mrc(numfic(1)+1:numgrd(1)-numfic(1), &
	  numfic(2)+1:numgrd(2)-numfic(2))
  ENDIF

  !====================
  !  Write to LSX file:
  !====================

  ! Get the directory for LSX file:
  ext = 'lsx'
  CALL GET_DIRECTORY(ext,dir,len)

  ! Write data to a lsx file:
  CALL WRITE_LAPS_DATA(i4t,dir,ext,ngd(1),ngd(2),nvr,nvr, &
     		       vnm,lvl,crd,vun,cmt,dat,sts)

  ! Write data to a lsx file under balance:
  dir(len-3:len+8) = 'balance/lsx/'
  CALL WRITE_LAPS_DATA(i4t,dir(1:len+8),ext,ngd(1),ngd(2),nvr,nvr, &
     		       vnm,lvl,crd,vun,cmt,dat,sts)

  !====================
  !  Write to L1S file:
  !====================

  ! Get the directory for L1S file:
  ext = 'l1s'
  CALL GET_DIRECTORY(ext,dir,len)
  L1S_cmt(1) = 'LAPS 60 Minute Snow Accumulation'
  L1S_cmt(2) = 'storm total precip. accum.'
  L1S_vnm(1) = 'R01'
  L1S_vnm(2) = 'RTO'
  L1S_vun = 'M'
  L1S_crd = 'MSL'
  L1S_lvl = 0

  ! Write data to a lsx file:
  ! Temporarily turn off L1S output as soil analysis needs snow as well
  ! however, Mile's work is only on rain.
  ! CALL WRITE_LAPS_DATA(i4t,dir,ext,ngd(1),ngd(2),2,2, &
  !   		       L1S_vnm,L1S_lvl,L1S_crd,L1S_vun,L1S_cmt,L1S,sts)
  ! End of writing a series of time frames:
  ENDDO

END SUBROUTINE PrPstLSX

SUBROUTINE PresChng(pres,ngrd,chng)

!==========================================================
!  This routine computes pressure change from pressure.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIe.
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ngrd(2)
  REAL, INTENT(IN) :: pres(ngrd(1),ngrd(2))
  REAL, INTENT(OUT) :: chng(ngrd(1),ngrd(2))

  ! Local variables:
  INTEGER :: i,j,nbi,nbj,nbs,ips
  REAL :: gam,kap,val,wgt,wsm
  REAL :: prs(ngrd(1),ngrd(2))

  ! Barnes parameters:
  gam = 0.2
  kap = 2.0e3/gam	! Kapa_0
  nbs = 60		! Neighbors to analyzed
  
  ! Initial:
  chng = 0.0
  prs = 0.0

  ! Barnes analysis:
  DO ips=1,1

    kap = kap*gam	! Adjust kappa

    ! Every gridpoint:
    DO j=1,ngrd(2)
      DO i=1,ngrd(1)
	
	! For all considered neighbors: every other one to save time
	val = 0.0
	wsm = 0.0
	DO nbj=-nbs,nbs,4
	  DO nbi=-nbs,nbs,4

	    wgt = EXP(-(FLOAT(nbi)**2+FLOAT(nbj)**2)/kap)
	    IF ((i+nbi .GE. 1) .AND. (i+nbi .LE. ngrd(1)) .AND. &
	        (j+nbj .GE. 1) .AND. (j+nbj .LE. ngrd(2))) THEN
	      val = val+(pres(i+nbi,j+nbj)-prs(i+nbi,j+nbj))*wgt
	      wsm = wsm + wgt
	    ENDIF
	  
	  ENDDO
	ENDDO

        ! Update grid value:
        chng(i,j) = chng(i,j)+val/wsm

      ENDDO
    ENDDO

    ! Save iterated result:
    prs = chng

  ENDDO

  ! Compute pressure changes:
  chng = pres-prs

END SUBROUTINE PresChng

SUBROUTINE Extraplt(grid,ngrd)

!==========================================================
!  This routine extrapolates interior points to its boundary
!  values assuming ngrd > 3.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ngrd(2)
  REAL, INTENT(INOUT) :: grid(ngrd(1),ngrd(2))

  ! X:
  grid(1,2:ngrd(2)-1) = &
    2.0*grid(2,2:ngrd(2)-1)-grid(3,2:ngrd(2)-1)
  grid(ngrd(1),2:ngrd(2)-1) = &
    2.0*grid(ngrd(1)-1,2:ngrd(2)-1)-grid(ngrd(1)-2,2:ngrd(2)-1)

  ! Y:
  grid(1:ngrd(1),1) = &
    2.0*grid(1:ngrd(1),2)-grid(1:ngrd(1),3)
  grid(1:ngrd(1),ngrd(2)) = &
    2.0*grid(1:ngrd(1),ngrd(2)-1)-grid(1:ngrd(1),ngrd(2)-2)

END SUBROUTINE Extraplt

END MODULE PrePostPrc
