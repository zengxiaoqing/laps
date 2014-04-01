!dis    Forecast Systems Laboratory
!dis    NOAA/OAR/ERL/FSL
!dis    325 Broadway
!dis    Boulder, CO     80303
!dis
!dis    Forecast Research Division
!dis    Local Analysis and Prediction Branch
!dis    LAPS
!dis
!dis    This software and its documentation are in the public domain and
!dis    are furnished "as is."  The United States government, its
!dis    instrumentalities, officers, employees, and agents make no
!dis    warranty, express or implied, as to the usefulness of the software
!dis    and documentation for any purpose.  They assume no responsibility
!dis    (1) for the use of the software and documentation; or (2) to provide
!dis    technical support to users.
!dis
!dis    Permission to use, copy, modify, and distribute this software is
!dis    hereby granted, provided that the entire disclaimer notice appears
!dis    in all copies.  All modifications to this software must be clearly
!dis    documented, and are solely the responsibility of the agent making
!dis    the modifications.  If significant modifications or enhancements
!dis    are made to this software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis

MODULE LAPSDatSrc

!==========================================================
!  This module defines LAPS data retrieval functionality.
!
!  HISTORY:
!	Creation: YUANFU XIE	8-2005
!==========================================================

  USE Definition
  USE MEM_NAMELIST

CONTAINS

SUBROUTINE LAPSInfo

!==========================================================
!  This routine configures the STMAS for its analyses.
!  Tasks:
!	1. 
!
!  NOTE: three letter variables are local; six global;
!
!  HISTORY: 
! 	Creation: YUANFU XIE	6-2005
!==========================================================

  IMPLICIT NONE

  ! Local variables:
  CHARACTER*9 :: fnm
  INTEGER :: err	! Error indicator
  INTEGER :: i

  ! VARIABLES FOR ACCESSING WIND PARAMETERS:
  CHARACTER*150 :: STATIC,FILENM
  INTEGER       :: LENGTH

  !*********************
  ! LAPS configuration:
  !*********************

  ! Get number of gridpoints:
  CALL GET_GRID_DIM_XY(numgrd(1),numgrd(2),err)
  IF (err .ne. 1) PRINT*, 'STMAS>LAPSInfo: Error getting numgrd'

  ! Get LAPS cycle time:
  CALL GET_LAPS_CYCLE_TIME(lapsdt,err)
  IF (err .ne. 1) PRINT*, 'STMAS>LAPSInfo: Error getting cycle time'

  ! Get current system time:
  CALL GET_SYSTIME(i4time,fnm,err)
  WRITE(6,1) fnm(1:9), 'LSO'	! Assume STMAS using LSO instead of LSO_QC
1 format('STMAS>LAPSInfo: Getting surface data at: ',a9,' from ',a6)

  ! Get a flag for missing data:
  CALL GET_R_MISSING_DATA(mising,err)
  IF (err .ne. 1) PRINT*, 'STMAS>LAPSInfo: Error getting mising flag'

  ! Get a flag for bad surface data:
  CALL GET_SFC_BADFLAG(badsfc,err)
  IF (err .ne. 1) PRINT*, 'STMAS>LAPSInfo: Error getting badsfc flag'

  CALL GET_MAXSTNS(mxstts,err)
  IF (err .ne. 1) print*, 'STMAS>LAPSInfo: Error getting maxstations'

  ! USE A NEW LAPS WIND PARAMETER SCHEME:
  CALL GET_DIRECTORY('static',STATIC,LENGTH)
  FILENM = STATIC(1:LENGTH)//'/wind.nl'
  CALL READ_NAMELIST_LAPS ('wind',FILENM)

  ! Maximum number of stations includes both surface and SND surface:
  mxstts = mxstts+max_pr

  ! Check:
  IF (verbal .EQ. 1) THEN
    WRITE(*,2) numgrd(1:3),lapsdt,mxstts,mxstts-max_pr,max_pr,mising,badsfc
  ENDIF
2 FORMAT('STMAS>LAPSInfo: Num  gridpoints: ',3I6,/, &
	 'STMAS>LAPSInfo: LAPS cycle time: ',I6,/, &
	 'STMAS>LAPSInfo: Maxnumber sites: ',I6,/, &
	 'STMAS>LAPSInfo: Maxnumber surfs: ',I6,/, &
	 'STMAS>LAPSInfo: Maxnumber sonde: ',I6,/, &
	 'STMAS>LAPSInfo: Missing/bad ids: ',e16.8,f16.5)

END SUBROUTINE LAPSInfo

SUBROUTINE LAPSConf

!==========================================================
!  This routine reads in necessary configuration grid data.
!
!  HISTORY: 
!	Creation: 8-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  ! Local variables:
  CHARACTER*9 :: fnm,ext,var
  CHARACTER :: unt*60,com*60
  INTEGER :: i,j,err	! Error indicator

  ! Physical grid points and spacing info:
  ext = 'nest7grid'
  var = 'LAT'
  CALL RD_LAPS_STATIC(dirstc,ext,numgrd(1),numgrd(2),1, &
		      var,unt,com,latgrd,phydxy(1),err)
  IF (err .NE. 1) PRINT*,'STMAS>LAPSInfo: Error getting LAT'
  var = 'LON'
  CALL RD_LAPS_STATIC(dirstc,ext,numgrd(1),numgrd(2),1, &
		      var,unt,com,longrd,phydxy(2),err)
  IF (err .NE. 1) PRINT*,'STMAS>LAPSInfo: Error getting LON'

  ! Topography:
  CALL READ_STATIC_GRID(numgrd(1),numgrd(2),'AVG',topogr,err)
  IF (err .NE. 1) THEN
    WRITE(6,*) 'LAPSConf: error get LAPS AVG'
    STOP
  ENDIF

  var = 'LDF'
  CALL RD_LAPS_STATIC(dirstc,ext,numgrd(1),numgrd(2),1, &
		      var,unt,com,lndfac,grdspc(2),err)
  IF (err .NE. 1) PRINT*,'STMAS>LAPSInfo: Error getting LON'
  ! Removing meaningless land factors:
  DO j=1,numgrd(2)
    DO i=1,numgrd(1)
      IF (lndfac(i,j) .GT. 1.0) lndfac(i,j) = 1.0
      IF (lndfac(i,j) .LT. 0.0) lndfac(i,j) = 0.0
    ENDDO
  ENDDO

  ! Get reduced pressure level:
  call GET_LAPS_REDP(rdplvl,err)

  ! time step:
  grdspc(3) = lapsdt

  !*********************
  ! STMAS configuration:
  !*********************

  ! Analysis grid numbers:
  numgrd = numgrd+2*numfic

  ! Analysis domain:
  domain(1,1:2) = 1.0-numfic(1:2)	! X/Y: use grid numbers
  domain(2,1:2) = FLOAT(numgrd(1:2)-numfic(1:2))
  ! domain(1,3) = MOD(i4time-(numtmf-1)*lapsdt,86400) 
  ! domain(2,3) = domain(1,3)+(numtmf-1)*lapsdt
  ! New time window using i4time instead of HHMM:
  i4wndw(1) = i4time-(numtmf-1)*lapsdt
  i4wndw(2) = i4time
  domain(1,3) = 0.0
  domain(2,3) = i4wndw(2)-i4wndw(1)

  grdspc(1:2) = 1.0			! Based the domain setting

END SUBROUTINE LAPSConf

SUBROUTINE LAPSBKGD

!==========================================================
!  This routine reads into LGA background fields from LAPS.
!
!  HISTORY: 
! 	Creation: YUANFU XIE	8-2005
!==========================================================

  IMPLICIT NONE

  CHARACTER*31 :: ext	! extension of bkg file used (LAPS)
  INTEGER :: i,j,k,err	! err: 1 normal; 0 file not found
  INTEGER :: tim	! time of bkg file used (LAPS)
  INTEGER :: iwv,id,it	! index for v component of wind

  ! Set LAPS circle time frames:
  DO i=0,numtmf-1
    i4prev(numtmf-i) = i4time-i*lapsdt
  ENDDO

  ! Check LAPS time frames:
  IF (verbal .EQ. 1) THEN
    DO i=1,numtmf
      WRITE(*,11) i,i4prev(i),MOD(i4prev(i),86400)
    ENDDO
  ENDIF
11 FORMAT('STMAS>LAPSBKGD: LAPS Time stamp',I3,':',I11,I7)

  ! Read background fields:
  DO j=1,numvar

    ! Wind:
    IF (varnam(j) .EQ. 'WNDU') THEN
      ! Search index for v component of wind:
      iwv = 0
      DO i=1,numvar
        IF (varnam(i) .EQ. 'WNDV') iwv = i
      ENDDO
      IF (iwv .EQ. 0) THEN
	WRITE(*,12)
	STOP
      ENDIF
12 FORMAT('STMAS>LAPSBKGD: V component of wind is missing!')

      ! Get wind fields:
      IF (needbk(j) .EQ. 1) THEN
        DO i=1,numtmf
          CALL GET_BKGWIND_SFC(i4prev(i),ext,tim, &
            bkgrnd(1,1,i,j),bkgrnd(1,1,i,iwv),lapsdt, &
            numgrd(1),numgrd(2),err)
          IF (err .EQ. 0) WRITE(*,13) varnam(j),i4prev(i),i,j
        ENDDO
      ELSE
        bkgrnd(1:numgrd(1),1:numgrd(2),i,j) = 0.0
        bkgrnd(1:numgrd(1),1:numgrd(2),i,iwv) = 0.0
      ENDIF
    !PCP1:
    ELSE IF (varnam(j) .EQ. 'PCP1') THEN        ! precip 1hr bkg        added by min-ken hsieh
      ! Get pcp1hr fields:
      IF (needbk(j) .EQ. 1) THEN
        DO i=1,numtmf
          CALL GET_MODELFG_2D(i4prev(i),'PCP',numgrd(1),numgrd(2),bkgrnd(1,1,i,j),err)
          IF(err .NE. 1)THEN
            WRITE(6,*)' No model first guess preicp, using zero field'
            bkgrnd(1:numgrd(1),1:numgrd(2),i,j) = 0.0
          ENDIF
        ENDDO
      ELSE
        bkgrnd(1:numgrd(1),1:numgrd(2),1:numtmf,j) = 0.0
      ENDIF
    ELSE IF (varnam(j) .EQ. 'PCP3') THEN        ! precip 3hr bkg        added by min-ken hsieh
      ! assign zero fields:
      DO i=1,numtmf
        bkgrnd(1:numgrd(1),1:numgrd(2),i,j) = 0.0
      ENDDO
    ELSE IF (varnam(j) .EQ. 'PCP6') THEN        ! precip 6hr bkg        added by min-ken hsieh
      ! assign zero fields:
      DO i=1,numtmf
        bkgrnd(1:numgrd(1),1:numgrd(2),i,j) = 0.0
      ENDDO
    ELSE IF (varnam(j) .EQ. 'PC24') THEN        ! precip 24hr bkg       added by min-ken hsieh
      ! assign zero fields:
      DO i=1,numtmf
        bkgrnd(1:numgrd(1),1:numgrd(2),i,j) = 0.0
      ENDDO
    ELSE IF (varnam(j) .EQ. 'TGD ') THEN	! Ground/Skin temperature added by Yuanfu
      IF (needbk(j) .EQ. 1) THEN
        DO i=1,numtmf
          CALL GET_MODELFG_2D(i4prev(i),varnam(j),numgrd(1),numgrd(2),bkgrnd(1,1,i,j),err)
          ! Temporarily add a test of TGD. If it is small (<1.0), assume
          ! laps does not provide good TGD data. When Paula or Steve fix LGA/LGB
          ! we can change back to test err only!!!
          ! IF (err .NE. 1) THEN
          IF (err .NE. 1 .or. bkgrnd(1,1,i,j) .lt. 1.0) THEN
            ! Use t and td and landfactor as used in LAPS sfc:
            it = 0	! search temp read in already
            id = 0	! search dewp read in already
            DO k=1,j-1
              IF (varnam(k) .EQ. 'TEMP') it = k
              IF (varnam(k) .EQ. 'DEWP') id = k
            ENDDO
            IF (it .EQ. 0 .AND. id .EQ. 0) THEN
              PRINT*,'LAPSBKGD: Place skin/ground temperature analysis after TEMP/DEWP in stmas_mg.vr'
              STOP
            ENDIF
            PRINT*,'LAPSBKGD: No background for skin/ground temp; use sfc temp and dewp'
            bkgrnd(1:numgrd(1),1:numgrd(2),i,j) = &
              bkgrnd(1:numgrd(1),1:numgrd(2),i,it)*lndfac(1:numgrd(1),1:numgrd(2))+ &
              bkgrnd(1:numgrd(1),1:numgrd(2),i,id)*(1.0-lndfac(1:numgrd(1),1:numgrd(2)))
          ENDIF
        ENDDO
      ELSE
        bkgrnd(1:numgrd(1),1:numgrd(2),1:numtmf,j) = 0.0
      ENDIF
    !Other fields:
    ELSE IF ((varnam(j) .NE. 'WNDV') .AND. &    ! V in with U
             (varnam(j) .NE. 'CEIL')) THEN      ! No ceiling bkg
      IF (needbk(j) .EQ. 1) THEN
        DO i=1,numtmf
          CALL GET_BACKGROUND_SFC(i4prev(i),varnam(j),ext,tim, &
            bkgrnd(1,1,i,j),lapsdt,numgrd(1),numgrd(2),err)
          IF (err .EQ. 0) THEN
            WRITE(*,13) varnam(j),i4prev(i),i,j
            STOP
          ENDIF
        ENDDO
      ELSE
        bkgrnd(1:numgrd(1),1:numgrd(2),1:numtmf,j) = 0.0
      ENDIF
    ENDIF
  ENDDO
13 FORMAT('STMAS>LAPSBKGD: Background is not found for: ',A4,i16,2i3)

  ! Save background if requested for debugging:
  IF (savdat .EQ. 1) THEN
    OPEN(unit=10,file='STMAS_bkg.dat',form='formatted')
    WRITE(10,*) numgrd(1:2),numtmf,domain
    WRITE(10,*) bkgrnd(1:numgrd(1),1:numgrd(2),1:numtmf,saveid)
    CLOSE(10)
  ENDIF

END SUBROUTINE LAPSBKGD

SUBROUTINE LAPSOBSV(m)

!==========================================================
!  This routine reads in LSO observation data from LAPS.
!
!  HISTORY:
! 	Creation: YUANFU XIE	8-2005
!	Modified: 6-2006 by YUANFU XIE: add sfc pressure.
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: m	! Maximum number of sites

  ! Local variables:
  CHARACTER*24 :: tim		! Obs file time
  CHARACTER :: stn(m)*20	! Station names
  CHARACTER :: prd(m)*11	! Provider names
  CHARACTER :: pwx(m)*25	! Present weather
  CHARACTER :: rtp(m)*6		! Report type
  CHARACTER :: stp(m)*6		! Station type (manual/auto)
  CHARACTER :: amt(m,5)*4	! cloud amount

  INTEGER :: nog,nob		! Number obs over grid/box
  INTEGER :: nss		! Sfc sonde data start number = nob+1
  INTEGER :: nsg,nsb		! Number sfc sonde obs over grid/box
  INTEGER :: wid(m)		! WMO id
  INTEGER :: otm(m)		! Observation time
  INTEGER :: cld(m)		! Number of cloud layers

  REAL :: lat(m),lon(m), &	! Lat/Lon
	    elv(m)		! Elevation
  REAL :: tmp(m),tmpea(m), &	! Temperature/expected accuracy
	    dew(m),dewea(m), &	! Dewpoint/EA
	    rhd(m),rhdea(m), &	! Relative Humidity/EA
	    wdi(m),wdiea(m), &	! Wind direction/EA
	    spd(m),spdea(m), &	! Wind speed/EA
	    gdi(m),	     &	! Gust wind direction
	    gsp(m),          &	! Gust wind speed
	    alt(m),altea(m), &	! Altimeter/EA
	    spr(m),prsea(m), &	! Station pressure/EA
	    msp(m), &		! Mean sea level pressure/EA
	    pcc(m),pccea(m), &	! 3-hour pressure change character
	    pch(m),pchea(m), &	! 3-hour pressure change
	    vis(m),visea(m), &	! Visibility/EA
	    sol(m),solea(m), &	! Solar/EA
	    slt(m),sltea(m), &	! Soil/water temperature/EA
	    slm(m),slmea(m), &	! Soil moist/EA
	    pc1(m),pcpea(m), &	! 1-hour precipitation/EA
	    pc3(m),pc6(m), &
	    p24(m), &		! 3,6,24-hour precipitation
	    snw(m),snwea(m), &	! Snow depth/EA
	    mxt(m),mnt(m)	! 24-hour maximum/minimum temperature
  REAL :: cht(m,5)		! cloud layer heights
  
  INTEGER :: i,j,k,err,iwv,kmax ! kmax need for reading sonde sfc data
  ! Old obs time treatment:
  ! INTEGER :: hrs,mns,nit	! Time: hours, minutes and mid-night

  ! No matter what sts value, get_sfc_obtime converts the obstime:
  INTEGER :: ot,sts		! Obstime conversion status 

  INTEGER :: jmin,jmax          ! For checking the min/max obs
  REAL :: vmin,vmax

  REAL :: xyt(3)		! X, Y and T
  REAL :: prs,ALT_2_SFC_PRESS

  REAL :: t_c,d_c,f_to_c,c_to_f,dwpt	! for converting rh to td

  ! Terrain interpolation variables:
  INTEGER :: ix,ix1,jy,jy1
  REAL :: alpha,beta,topo

  ! Read observation data by LAPS time frames:
  numobs = 0
  DO i=1,numtmf

    ! Frame by frame: READ_SURFACE_DATAQC or READ_SURFACE_DATA
    CALL READ_SURFACE_DATA(i4prev(i),tim,nog,nob, &
	otm,wid,stn,prd,pwx,rtp,stp,lat,lon,elv, &
	tmp,dew,rhd,wdi,spd,gdi,gsp,alt,spr,msp,pcc,pch, &
	vis,sol,slt,slm,pc1,pc3,pc6,p24,snw,cld,mxt,mnt, &
	tmpea,dewea,rhdea,wdiea,spdea,altea,prsea,visea, &
	solea,sltea,slmea,pcpea,snwea,amt,cht,mxstts,err)

    IF (err .NE. 1) THEN
      PRINT*,'LAPSOBSV: error in reading LSO data, check!',i
    ENDIF

    ! Read in surface data from sondes: 
    ! (LAPS surface data is lso and snd both)
    nss = nob+1
    nsg = 0	! READ_SFC_SND does not initialize nsg
    nsb = 0	! READ_SFC_SND does not initialize nsb
    CALL GET_LAPS_DIMENSIONS(kmax,err)
    CALL READ_SFC_SND(i4prev(i),tim,nsg,nsb, &
        otm(nss),wid(nss),stn(nss),prd(nss),pwx(nss), &
        rtp(nss),stp(nss),lat(nss),lon(nss),elv(nss), &
        tmp(nss),dew(nss),rhd(nss),wdi(nss),spd(nss), &
        gdi(nss),gsp(nss),alt(nss),spr(nss),msp(nss),pcc(nss),pch(nss), &
        vis(nss),sol(nss),slt(nss),slm(nss),pc1(nss), &
        pc3(nss),pc6(nss),p24(nss),snw(nss),cld(nss),mxt(nss),mnt(nss), &
        tmpea(nss),dewea(nss),rhdea(nss),wdiea(nss),spdea(nss),altea(nss), &
        prsea(nss),visea(nss),solea(nss),sltea(nss),slmea(nss),pcpea(nss), &
        snwea(nss),amt(nss,1),cht(nss,1),mxstts,latgrd,longrd, &
        numgrd(1),numgrd(2),kmax,max_pr,max_pr_levels,topogr,sts)
    ! Combine surface data together:
    nob = nob+nsb	! Add snd data to the total

    IF (nob .EQ. 0) THEN
      PRINT*,'LAPSOBSV: No sfc obs found!'
      STOP
    ELSE
        ! Convert LAPS surface obs time to i4time:
      DO j=1,nob
        ot = otm(j)
        ! OTM from LSO is in HHMM:
        IF (ot .GE. 0 .AND. ot .LE. 2400) THEN
          CALL GET_SFC_OBTIME(ot,i4prev(i),otm(j),sts)
        ELSE
          ! OBS time error:
          otm(j) = -100
          IF (verbal .EQ. 1) &
            PRINT*,'LAPSOBS: Invalid Obs time: set to bad value: ',ot,j,i
        ENDIF
      ENDDO

    ENDIF

    IF (err .NE. 1 .AND. sts .NE. 1) THEN
      ! LSO data cannot be read in:
      WRITE(*,21) i
    ELSE
      ! Assign LSO data to the corresponding arrays:

      ! Retrieve location and time sequence:
      DO j=1,nob	! Through all obs sites
	! X and Y:
	CALL LATLON_TO_RLAPSGRID(lat(j),lon(j), &
		latgrd,longrd,numgrd(1),numgrd(2), &
		xyt(1),xyt(2),err)

        ! Good station locations:
        IF (err .EQ. 1) THEN

	  ! T: from LAPS time form: HHMM to seconds
	  ! hrs = otm(j)/100
	  ! mns = otm(j)-hrs*100
	  ! xyt(3) = hrs*3600+mns*60
	  ! IF (otm(j) .LT. 0) xyt(3) = 2*86400	! Void: Bad data

          ! Use i4time to handle obs times:
          xyt(3) = otm(j)-i4wndw(1)

	  ! Adjust the time when crossing the midnight:
	  IF ((xyt(3)+86400 .GE. domain(1,3)) .AND. &
	      (xyt(3)+86400 .LE. domain(2,3)) ) &
	    xyt(3) = xyt(3)+86400

	  ! Pass the location/time to obs arrays:
          DO k=1,numvar
	    rawobs(2:4,j+numobs(k),k) = xyt(1:3)
	  ENDDO

        ELSE  ! Bad station location
          DO k=1,numvar
            rawobs(1:4,j+numobs(k),k) = badsfc
          ENDDO
        ENDIF

      ENDDO

      ! Place the observations into right variables:
      DO j=1,numvar
	stanam(1+numobs(j):nob+numobs(j),j) = stn(1:nob)	!Added by min-ken,hsieh:stanam for STMASVer
	
	SELECT CASE (varnam(j))
	CASE ("TEMP")
	  rawobs(1,1+numobs(j):nob+numobs(j),j) = tmp(1:nob)
	  weight(1+numobs(j):nob+numobs(j),j) = tmpea(1:nob)
	CASE ("DEWP")
          ! Use td obs or rh obs if td missing:
          DO k=1,nob
            IF (dew(k) .NE. mising .AND. &
                dew(k) .NE. badsfc) THEN
              rawobs(1,numobs(j)+k,j) = dew(k)
	      weight(numobs(j)+k,j) = dewea(k)
            ELSEIF (rhd(k) .NE. mising .AND. &
                    rhd(k) .NE. badsfc .AND. &
                    tmp(k) .NE. mising .AND. &
                    tmp(k) .NE. badsfc) THEN
                ! obs value:
                t_c = f_to_c(tmp(k))
                d_c = dwpt(t_c,rhd(k))
                rawobs(1,numobs(j)+k,j) = c_to_f(d_c)
                dew(k) = rawobs(1,numobs(j)+k,j) ! Replace dew with one from RH
                ! obs expect accuracy:
                t_c = f_to_c(tmpea(k))
                d_c = dwpt(t_c,rhdea(k))
                weight(numobs(j)+k,j) = c_to_f(d_c)
            ELSE
              rawobs(1,numobs(j)+k,j) = badsfc
	      weight(numobs(j)+k,j) = badsfc
            ENDIF
          ENDDO
	CASE ("VISB")
	  rawobs(1,1+numobs(j):nob+numobs(j),j) = vis(1:nob)
	  weight(1+numobs(j):nob+numobs(j),j) = visea(1:nob)
        CASE ("CEIL")
	  rawobs(1,1+numobs(j):nob+numobs(j),j) = cht(1:nob,1)
	  weight(1+numobs(j):nob+numobs(j),j) = 1.0
        CASE ("MSLP")
	  rawobs(1,1+numobs(j):nob+numobs(j),j) = msp(1:nob)
	  weight(1+numobs(j):nob+numobs(j),j) = prsea(1:nob)
        CASE ("TGD ")
	  rawobs(1,1+numobs(j):nob+numobs(j),j) = slt(1:nob)	! Soil temp
	  weight(1+numobs(j):nob+numobs(j),j) = sltea(1:nob)
	CASE ("REDP")
	  DO k=1,nob
	    ! Collect either station pressure or altimeter:
	    IF ((spr(k) .NE. mising) .AND. &
		(spr(k) .NE. badsfc)) THEN
	      prs = spr(k)
	    ELSEIF ((alt(k) .NE. mising) .AND. &
		    (alt(k) .NE. badsfc) .AND. &
                    (rawobs(2,numobs(j)+k,j) .NE. mising) .AND. &
                    (rawobs(3,numobs(j)+k,j) .NE. badsfc)) THEN

              ! Grid indices and interpolation coefficiences
              ix = INT(rawobs(2,numobs(j)+k,j))
              ix1 = MIN(ix+1,numgrd(1))
              jy = INT(rawobs(3,numobs(j)+k,j))
              jy1 = MIN(jy+1,numgrd(2))

              alpha = rawobs(2,numobs(j)+k,j)-ix
              beta  = rawobs(3,numobs(j)+k,j)-jy

              topo = (1.0-alpha)*(1.0-beta)*topogr(ix ,jy )+ &
                          alpha *(1.0-beta)*topogr(ix1,jy )+ &
                     (1.0-alpha)*     beta *topogr(ix ,jy1)+ &
                          alpha *     beta *topogr(ix1,jy1)

	      prs = ALT_2_SFC_PRESS(alt(k),topo)
	    ELSE
	      prs = badsfc
	    ENDIF

	    ! Convert to reduced pressure at the reduced pressure level, rdplvl:
            IF (prs .NE. badsfc) THEN
              IF (ABS(elv(k)-rdplvl) .LE. 10) THEN 
                ! Elevation is close (10m hardcoded) to the reduced level, 
                ! use pressure as reduced pressure:
                rawobs(1,numobs(j)+k,j) = prs
              ELSE
                ! Use reduce_p to reduce prs to the reduced pressure level:
	        IF (tmp(k) .NE. badsfc .AND. dew(k) .NE. badsfc) THEN
	          CALL REDUCE_P(tmp(k),dew(k),prs,topo, &
		                lapses(1),lapses(2), &
		                rawobs(1,numobs(j)+k,j),rdplvl,badsfc)
                ELSE
	          rawobs(1,numobs(j)+k,j) = badsfc
                ENDIF
              ENDIF
	    ELSE
	      rawobs(1,numobs(j)+k,j) = badsfc
            ENDIF
	  ENDDO
	  weight(1+numobs(j):nob+numobs(j),j) = 1.0 ! altea(1:nob)
	CASE ("SFCP")
	  DO k=1,nob
	    ! Collect either station pressure or altimeter:
	    IF ((spr(k) .NE. mising) .AND. &
		(spr(k) .NE. badsfc)) THEN
	      rawobs(1,numobs(j)+k,j) = spr(k)
	    ELSEIF ((alt(k) .NE. mising) .AND. &
		    (alt(k) .NE. badsfc) .AND. &
                    (rawobs(2,numobs(j)+k,j) .NE. mising) .AND. &
                    (rawobs(3,numobs(j)+k,j) .NE. badsfc) ) THEN

              ! Grid indices and interpolation coefficiences
              ix = INT(rawobs(2,numobs(j)+k,j))
              ix1 = MIN(ix+1,numgrd(1))
              jy = INT(rawobs(3,numobs(j)+k,j))
              jy1 = MIN(jy+1,numgrd(2))

              alpha = rawobs(2,numobs(j)+k,j)-ix
              beta  = rawobs(3,numobs(j)+k,j)-jy

              topo = (1.0-alpha)*(1.0-beta)*topogr(ix ,jy )+ &
                          alpha *(1.0-beta)*topogr(ix1,jy )+ &
                     (1.0-alpha)*     beta *topogr(ix ,jy1)+ &
                          alpha *     beta *topogr(ix1,jy1)

	      rawobs(1,numobs(j)+k,j) = ALT_2_SFC_PRESS(alt(k),topo)
if (stn(k) .EQ. 'KSPR') then
print*,'E1452: ',topo, elv(k),alt(k),rawobs(1,numobs(j)+k,j)
endif

	    ELSE
	      rawobs(1,numobs(j)+k,j) = badsfc
	    ENDIF
	  ENDDO
	  weight(1+numobs(j):nob+numobs(j),j) = 1.0 ! altea(1:nob)
	CASE ("WNDU")
	  ! Find the index for v component:
	  iwv = 0
	  DO k=1,numvar
	    IF (varnam(k) .EQ. "WNDV") iwv = k
	  ENDDO
	  IF (iwv .EQ. 0) THEN
	    WRITE(*,24)
	  ENDIF
24 FORMAT('STMAS>LAPSOBS: Warning: no v component wind analysis!')

	  ! Convert wind from direction/speed to U/V:
	  DO k=1,nob
	    IF ((wdi(k) .EQ. mising) .OR. &
	        (wdi(k) .EQ. badsfc) .OR. &
	        (spd(k) .EQ. mising) .OR. &
		(spd(k) .EQ. badsfc)) THEN
	      rawobs(1,numobs(j)+k,j) = badsfc
	      rawobs(1,numobs(j)+k,iwv) = badsfc
	    ELSE
	      ! Conversion:
	      CALL DISP_TO_UV(wdi(k),spd(k),xyt(1),xyt(2))
	      CALL UVTRUE_TO_UVGRID(xyt(1),xyt(2), &
		rawobs(1,numobs(j)+k,j), &
		rawobs(1,numobs(j)+k,iwv),lon(k))
	    ENDIF
	  ENDDO
	  weight(1+numobs(j):nob+numobs(j),j) = 1.0
	  weight(1+numobs(j):nob+numobs(j),iwv) = 1.0
	CASE ("WNDV")
	  ! V should be in already. See CASE ("WNDU").
	CASE ("PCP1")						!added by min-ken hsieh
	  rawobs(1,1+numobs(j):nob+numobs(j),j) = pc1(1:nob)
	  weight(1+numobs(j):nob+numobs(j),j) = pcpea(1:nob)
	CASE ("PCP3")						!added by min-ken hsieh
	  rawobs(1,1+numobs(j):nob+numobs(j),j) = pc3(1:nob)
	  weight(1+numobs(j):nob+numobs(j),j) = pcpea(1:nob)
	CASE ("PCP6")						!added by min-ken hsieh
	  rawobs(1,1+numobs(j):nob+numobs(j),j) = pc6(1:nob)
	  weight(1+numobs(j):nob+numobs(j),j) = pcpea(1:nob)
	CASE ("PC24")						!added by min-ken hsieh
	  rawobs(1,1+numobs(j):nob+numobs(j),j) = p24(1:nob)
	  weight(1+numobs(j):nob+numobs(j),j) = pcpea(1:nob)
	CASE DEFAULT
	  WRITE(*,22) varnam(j)
	END SELECT
      ENDDO
    ENDIF

    ! Update frm:
    numobs(1:numvar) = numobs(1:numvar)+nob
  ENDDO
21 FORMAT('STMAS>LAPSOBS: Cannot read in LSO data and SND data: ',i8)
22 FORMAT('STMAS>LAPSOBS: No such var in LSO data: ',A4)

  ! Check number of obs:
  DO i=1,numvar
    WRITE(*,23) varnam(i),numobs(i)
    ! IF (verbal .EQ. 1) WRITE(*,23) varnam(i),numobs(i)
  ENDDO
23 FORMAT('STMAS>LAPSOBSV: NumObs of (raw) ',A4,': ',I8)

  ! Remove invalid data:
  CALL RmvInvld

  ! Remove redundant observations:
  ! CALL RmvDupls

  ! Check the data ranges:
  IF (verbal .EQ. 1) THEN
    DO i=1,numvar
      vmin = 1000.0
      vmax = -1000.0
      DO j=1,numobs(i)
        IF (rawobs(1,j,i) .LT. vmin) THEN
          vmin = rawobs(1,j,i)
          jmin = j
        ENDIF
        IF (rawobs(1,j,i) .GT. vmax) THEN
          vmax = rawobs(1,j,i)
          jmax = j
        ENDIF
      ENDDO
       
      WRITE(*,26) varnam(i),vmin,jmin,stanam(jmin,i),vmax,jmax,stanam(jmax,i)
!		  MINVAL(rawobs(1,1:numobs(i),i)), &
!		  MAXVAL(rawobs(1,1:numobs(i),i))
    ENDDO
  ENDIF
26 FORMAT('STMAS>LAPSOBSV: ',A4,' min/max values: ', F11.2,i4,1x,a6,F11.2,i4,1x,a6)

  ! Write out requested obs for testing:
  IF (savdat .EQ. 1) THEN
    OPEN(unit=10,file='STMAS_ob1.dat',form='formatted')
    WRITE(10,*) numobs(saveid),numtmf,domain,grdspc(3)
    WRITE(10,*) rawobs(1:4,1:numobs(saveid),saveid)
    CLOSE(10)
  ENDIF

END SUBROUTINE LAPSOBSV

SUBROUTINE RmvInvld

!==========================================================
!  This routine removes the invalid data (missing/badsfc)
!  data from observations.
!
!  HISTORY:
!	Creation: 8-2005 by YUANFU XIE
!==========================================================

  IMPLICIT NONE

  INTEGER :: i,j,num

  ! For every variables:
  DO i=1,numvar

    num = numobs(i)
    numobs(i) = 0
    ! For every obs:
    DO j=1,num
      IF ((rawobs(1,j,i) .NE. mising) .AND. &
	  (rawobs(1,j,i) .NE. badsfc)) THEN
	! Valid data:
	numobs(i) = numobs(i)+1
	rawobs(1:4,numobs(i),i) = rawobs(1:4,j,i)
	weight(numobs(i),i) = weight(j,i)
	stanam(numobs(i),i) = stanam(j,i)		!Added by min-ken,hsieh:stanam for STMASVer
      ENDIF
    ENDDO
  ENDDO
 
  ! Check numbers of obs left:
  DO i=1,numvar
    IF (verbal .EQ. 1) WRITE(*,31) varnam(i),numobs(i)
  ENDDO
31 FORMAT('STMAS>LAPS_QCs: NumObs of (Vld) ',A4,': ',I8)

END SUBROUTINE RmvInvld

SUBROUTINE RmvDupls

!==========================================================
!  This routine removes the duplicated observations.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER :: i,j,k,l,nob
  REAL :: dis

  obsspc = 1.0e10
  DO i=1,numvar

    ! Use weight as flag:
    DO j=1,numobs(i)
      DO k=j+1,numobs(i)
	dis = 0.0
	DO l=2,4
	  dis = dis+(rawobs(l,j,i)-rawobs(l,k,i))* &
		    (rawobs(l,j,i)-rawobs(l,k,i))
	ENDDO

	! Mark those redundants:
	IF (dis .LT. epsiln) THEN
	  weight(k,i) = 0.0
	  IF (verbal .EQ. 1) THEN
	    WRITE(*,11) rawobs(1:4,j,i),rawobs(1:4,k,i), &
		varnam(i),j,k
	  ENDIF
	ELSE
	  ! Minimal observation spacing:
          IF ((ABS(rawobs(2,j,i)-rawobs(2,k,i)) .GT. 0.0) .AND. &
              (ABS(rawobs(2,j,i)-rawobs(2,k,i)) .LT. obsspc(1,i))) &
            obsspc(1,i) = ABS(rawobs(2,j,i)-rawobs(2,k,i))
          IF ((ABS(rawobs(3,j,i)-rawobs(3,k,i)) .GT. 0.0) .AND. &
              (ABS(rawobs(3,j,i)-rawobs(3,k,i)) .LT. obsspc(2,i))) &
            obsspc(2,i) = ABS(rawobs(3,j,i)-rawobs(3,k,i))
          IF ((ABS(rawobs(4,j,i)-rawobs(4,k,i)) .GT. 0.0) .AND. &
              (ABS(rawobs(4,j,i)-rawobs(4,k,i)) .LT. obsspc(3,i))) &
            obsspc(3,i) = ABS(rawobs(4,j,i)-rawobs(4,k,i))
	ENDIF
      ENDDO
11 FORMAT('STMAS>RmvDupls: Redundant data: ',/,4F14.4,/, &
	4F14.4,/,A6,2I8)

    ENDDO
    WRITE(*,12) varnam(i),obsspc(1:3,i)
12 FORMAT('STMAS>RmvDupls: Minimal obs (',A4,') spacing: ',3E12.4)

    ! Remove redundants:
    nob = 0
    DO j=1,numobs(i)
      IF (weight(j,i) .GT. 0.0) THEN
        nob = nob+1
        rawobs(1:4,nob,i) = rawobs(1:4,j,i)
	weight(nob,i) = weight(j,i)
	stanam(nob,i) = stanam(j,i)		!Added by Yuanfu: stanam for STMASVer
      ENDIF
    ENDDO
    numobs(i) = nob

  ENDDO

END SUBROUTINE RmvDupls

SUBROUTINE LAPSUnit

!==========================================================
!  This routine converts LSO observation units into a unit
!  consistent with the background.
!
!  HISTORY:
!	Creation: 8-2005 by YUANFU XIE.
!	Modified: 6-2006 by YUANFU XIE: add sfc pressure.
!==========================================================

  IMPLICIT NONE

  INTEGER :: i

  ! Check all variables:
  DO i=1,numvar

    ! Find necessary conversion:
    SELECT CASE (varnam(i))
    CASE ("TEMP")
      ! Convert to Kelvin from Fahrenheit:
      rawobs(1,1:numobs(i),i) = &
	(rawobs(1,1:numobs(i),i)-32.0)*5.0/9.0+temp_0
    CASE ("DEWP")
      ! Convert to Kelvin from Fahrenheit:
      rawobs(1,1:numobs(i),i) = &
	(rawobs(1,1:numobs(i),i)-32.0)*5.0/9.0+temp_0
    CASE ("TGD ")
      ! Convert to Kelvin from Fahrenheit:
      rawobs(1,1:numobs(i),i) = &
	(rawobs(1,1:numobs(i),i)-32.0)*5.0/9.0+temp_0
    CASE ("WNDU")
      ! Convert to m/s from knots:
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*knt2ms
    CASE ("WNDV")
      ! Convert to m/s from knots:
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*knt2ms
    CASE ("VISB")
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*mile2m
    CASE ("REDP")
      ! Convert to pascal from mb:
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*mb2pas
    CASE ("SFCP")
      ! Convert to pascal from mb:
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*mb2pas
    CASE ("MSLP")
      ! Convert to pascal from mb:
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*mb2pas
    CASE ("PCP1")
      ! Convert to meter from inch:
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*inch2m
    CASE ("PCP3")
      ! Convert to meter from inch:
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*inch2m
    CASE ("PCP6")
      ! Convert to meter from inch:
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*inch2m
    CASE ("PC24")
      ! Convert to meter from inch:
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*inch2m
    END SELECT
  ENDDO

END SUBROUTINE LAPSUnit

SUBROUTINE LAPS_QCs

!==========================================================
!  This routine runs quality control over data by threshold
!  values and standard deviation.
!
!  HISTORY:
! 	Creation: YUANFU XIE	6-2005
!==========================================================

  IMPLICIT NONE

  ! Interpolation indices and coefficients:
  CALL LAPSIntp

  ! Optional QCs:
  IF (qc_val .EQ. 1) CALL Thrshold

  ! Save QCed obs:
  CALL CpyQCObs

END SUBROUTINE LAPS_QCs

SUBROUTINE CpyQCObs

!==========================================================
!  This routine copies QCed observation data from rawobs to
!  qc_obs after all QC is done. This routine can be avoid
!  if the qc_obs array keeps the same structure as rawobs.
!
!  HISTORY:
!	Creation: 8-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER :: i,j

  ! Copy:
  DO i=1,numvar
    DO j=1,numobs(i)
      qc_obs(1:4,j,i) = rawobs(1:4,j,i)

      ! Save innovations:
      IF (needbk(i) .EQ. 1) &
        qc_obs(1,j,i) = qc_obs(1,j,i)-bkgobs(j,i)

    ENDDO
  ENDDO

END SUBROUTINE CpyQCObs

SUBROUTINE Thrshold

!==========================================================
!  This routine does the threshold value QC checks.
!
!  HISTORY:
!	Creation: 8-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER :: i,j,num

  ! Check:
  DO i=1,numvar
    IF (needbk(i) .EQ. 1) THEN
      num = numobs(i)
      numobs(i) = 0
      DO j=1,num

        ! QC check: avoid bkg = mising with roundoff error:
        IF (ABS(rawobs(1,j,i)-bkgobs(j,i)) .LE. thresh(i)) THEN
	  numobs(i) = numobs(i)+1
	  rawobs(1:4,numobs(i),i) = rawobs(1:4,j,i)
	  weight(numobs(i),i) = weight(j,i)
	  stanam(numobs(i),i) = stanam(j,i)		!Added by min-ken,hsieh:stanam for STMASVer
	  indice(1:6,numobs(i),i) = indice(1:6,j,i)
	  coeffs(1:6,numobs(i),i) = coeffs(1:6,j,i)
	  bkgobs(numobs(i),i) = bkgobs(j,i)
        ELSE
          print*,'Thresh out: ',rawobs(1,j,i),bkgobs(j,i),thresh(i),j,i,stanam(j,i)
        ENDIF
      ENDDO
    ENDIF
  ENDDO
 
  ! Check numbers of obs left:
  DO i=1,numvar
    IF (verbal .EQ. 1) WRITE(*,31) varnam(i),numobs(i)
  ENDDO
31 FORMAT('STMAS>LAPS_QCs: NumObs of (Vlu) ',A4,': ',I8)

END SUBROUTINE Thrshold

SUBROUTINE LAPSIntp

!==========================================================
!  This routine interpolates gridpoints to observation site
!  and saves the indices and coefficients.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!       Modification:
!                 25-08-2008 by min-ken hsieh
!                 pass stanam into Grid2Obs to map each obs its stn name for STMASVer
!==========================================================

  IMPLICIT NONE

  INTEGER :: i,j,ix,iy,it

  DO i=1,numvar
    CALL Grid2Obs(indice(1,1,i),coeffs(1,1,i), &
	rawobs(1,1,i),numobs(i),weight(1,i),stanam(1,i),numgrd, &
	grdspc,domain)

    ! Compute background values at observation sites:
    DO j=1,numobs(i)

      ! Interpolate background to the obs site:
      bkgobs(j,i) = 0.0
      DO it=3,6,3
        DO iy=2,5,3
          DO ix=1,4,3
            bkgobs(j,i) = bkgobs(j,i) + &
	      bkgrnd(indice(ix,j,i), &
		     indice(iy,j,i), &
		     indice(it,j,i),i)* &
	      coeffs(ix,j,i)*coeffs(iy,j,i)*coeffs(it,j,i)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  ENDDO

END SUBROUTINE LAPSIntp

SUBROUTINE STMASVer
!==========================================================
!  This routine prepare all parameters and pass them to
!  verify subroutine in src/lib/laps_routine.f.
!  it will output verify file in log/qc directory, just
!  like what laps_sfc.x dose.
!  HISTORY:
!       Creation: 22-8-2008 by min-ken hsieh.
!==========================================================

  IMPLICIT NONE

  !Local vaiables
  INTEGER :: i, j, len, istatus, err
  INTEGER :: iunit     		 !log file handle
  INTEGER :: ii(mxstts),jj(mxstts)

  ! since we use bilinear interpolation, these arrays are dummy..
  REAL :: x1a(numgrd(1)), x2a(numgrd(2)), y2a(numgrd(1),numgrd(2)) 
  REAL :: ea

  CHARACTER*60  :: title
  CHARACTER*256 :: ver_file
  CHARACTER*9   :: a9time

  !for time loop
  INTEGER :: kt,nn
  INTEGER :: obstime
  REAL :: obsOfThisTime(mxstts)
  CHARACTER*20 :: staOfThisTime(mxstts)
  CHARACTER*1 :: tmtag

  ! open log file
  iunit = 11
  CALL make_fnam_lp(i4time,a9time,istatus)
  IF(istatus .eq. 0) GOTO 999
  CALL get_directory('log', ver_file, len)
  ver_file = ver_file(1:len)//'qc/stmas.ver.'//a9time(6:9)
  CALL s_len(ver_file, len)
  !PRINT*, "min-ken",ver_file
  OPEN(iunit,file=ver_file(1:len),status='unknown',err=999)

  ! variable loop
  DO i=1,numvar
    !because we only care about real obs
    !we do not apply verify to those obs made by bkgrnd

    ! qc_obs array actually store obs - obsbkg field (done by CpyQCObs)
    ! we need to add obsbkg back to qc_obs
    IF (needbk(i) .EQ. 1) &
      qc_obs(1,1:numobs(i),i) = qc_obs(1,1:numobs(i),i)+bkgobs(1:numobs(i),i)



    SELECT CASE (varnam(i))
      CASE ("TEMP")

	!time loop
	DO kt = 1,numtmf
	  write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          DO j=1,numobs(i)
            IF(INT(qc_obs(4,j,i)).EQ.obstime) THEN
              nn = nn + 1
	      obsOfThisTime(nn) = qc_obs(1,j,i)
	      staOfThisTime(nn) = stanam(j,i)
	      ii(nn) = qc_obs(2,j,i)
	      jj(nn) = qc_obs(3,j,i)
            ENDIF
          ENDDO
	
          title = 'Temperature background verification of tmf = '//tmtag//' (deg C)'
          ea = 1.50*5.0/9.0
          CALL verify(bkgrnd(1:numgrd(1),1:numgrd(2),kt,i),obsOfThisTime(1:nn),		&
      	    	    staOfThisTime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
 
          title = 'Temperature verification of tmf = '//tmtag//' (deg C)'
          CALL verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsOfThisTime(1:nn),		&
      	    	    staOfThisTime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	ENDDO

      CASE ("WNDU")

        !time loop
        DO kt = 1,numtmf
          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          DO j=1,numobs(i)
            IF(INT(qc_obs(4,j,i)).EQ.obstime) THEN
              nn = nn + 1
              obsOfThisTime(nn) = qc_obs(1,j,i)
              staOfThisTime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            ENDIF
          ENDDO

          title = 'U Wind Component background verification of tmf = '//tmtag//' (m/s)'
          ea = 2.00*knt2ms
          CALL verify(bkgrnd(1:numgrd(1),1:numgrd(2),kt,i),obsOfThisTime(1:nn),		&
      	    	    staOfThisTime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)

          title = 'U Wind Component verification of tmf = '//tmtag//' (m/s)'
          CALL verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsOfThisTime(1:nn),		&
      	    	    staOfThisTime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	ENDDO

      CASE ("WNDV")

        !time loop
        DO kt = 1,numtmf
          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          DO j=1,numobs(i)
            IF(INT(qc_obs(4,j,i)).EQ.obstime) THEN
              nn = nn + 1
              obsOfThisTime(nn) = qc_obs(1,j,i)
              staOfThisTime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            ENDIF
          ENDDO

          title = 'V Wind Component background verification of tmf = '//tmtag//' (m/s)'
          ea = 2.00*knt2ms
          CALL verify(bkgrnd(1:numgrd(1),1:numgrd(2),kt,i),obsOfThisTime(1:nn),		&
      	    	    staOfThisTime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)

          title = 'V Wind Component verification of tmf = '//tmtag//' (m/s)'
          CALL verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsOfThisTime(1:nn),		&
      	    	    staOfThisTime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	ENDDO

      CASE ("DEWP")

        !time loop
        DO kt = 1,numtmf
          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          DO j=1,numobs(i)
            IF(INT(qc_obs(4,j,i)).EQ.obstime) THEN
              nn = nn + 1
              obsOfThisTime(nn) = qc_obs(1,j,i)
              staOfThisTime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            ENDIF
          ENDDO

          title = 'Dew Point background verification of tmf = '//tmtag//' (deg C)'
          ea = 2.00*5.0/9.0
          CALL verify(bkgrnd(1:numgrd(1),1:numgrd(2),kt,i),obsOfThisTime(1:nn),		&
      	    	    staOfThisTime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
 
          title = 'Dew Point verification of tmf = '//tmtag//' (deg C)'
          CALL verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsOfThisTime(1:nn),		&
      	    	    staOfThisTime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	ENDDO

      CASE ("TGD ")

        !time loop
        DO kt = 1,numtmf
          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          DO j=1,numobs(i)
            IF(INT(qc_obs(4,j,i)).EQ.obstime) THEN
              nn = nn + 1
              obsOfThisTime(nn) = qc_obs(1,j,i)
              staOfThisTime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            ENDIF
          ENDDO

          title = 'TGD background verification of tmf = '//tmtag//' (deg C)'
          ea = 2.00*5.0/9.0
          CALL verify(bkgrnd(1:numgrd(1),1:numgrd(2),kt,i),obsOfThisTime(1:nn),		&
      	    	    staOfThisTime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
 
          title = 'TDG verification of tmf = '//tmtag//' (deg C)'
          CALL verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsOfThisTime(1:nn),		&
      	    	    staOfThisTime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	ENDDO

      CASE ("REDP")
	!Convert pascal to mb
	qc_obs(1,1:numobs(i),i) = qc_obs(1,1:numobs(i),i)/mb2pas

        !time loop
        DO kt = 1,numtmf
 	  !Convert pascal to mb
	  analys(1:numgrd(1),1:numgrd(2),kt,i) = analys(1:numgrd(1),1:numgrd(2),kt,i)/mb2pas
	  bkgrnd(1:numgrd(1),1:numgrd(2),kt,i) = bkgrnd(1:numgrd(1),1:numgrd(2),kt,i)/mb2pas

          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          DO j=1,numobs(i)
            IF(INT(qc_obs(4,j,i)).EQ.obstime) THEN
              nn = nn + 1
              obsOfThisTime(nn) = qc_obs(1,j,i)
              staOfThisTime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            ENDIF
          ENDDO

          title = 'Reduced pressure background verification of tmf = '//tmtag//' (mb)'
          ea = 0.68
          CALL verify(bkgrnd(1:numgrd(1),1:numgrd(2),kt,i),obsOfThisTime(1:nn),		&
      	    	    staOfThisTime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)

          title = 'Reduced pressure verification of tmf = '//tmtag//' (mb)'
          CALL verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsOfThisTime(1:nn),		&
      	    	    staOfThisTime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	ENDDO


      CASE ("SFCP")
	!Convert pascal to mb
	qc_obs(1,1:numobs(i),i) = qc_obs(1,1:numobs(i),i)/mb2pas

        !time loop
        DO kt = 1,numtmf
	  !Convert pascal to mb
	  analys(1:numgrd(1),1:numgrd(2),kt,i) = analys(1:numgrd(1),1:numgrd(2),kt,i)/mb2pas
	  bkgrnd(1:numgrd(1),1:numgrd(2),kt,i) = bkgrnd(1:numgrd(1),1:numgrd(2),kt,i)/mb2pas

          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          DO j=1,numobs(i)
            IF(INT(qc_obs(4,j,i)).EQ.obstime) THEN
              nn = nn + 1
              obsOfThisTime(nn) = qc_obs(1,j,i)
              staOfThisTime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            ENDIF
          ENDDO

          title = 'SFC pressure background verification of tmf = '//tmtag//' (mb)'
          ea = 0.68
          CALL verify(bkgrnd(1:numgrd(1),1:numgrd(2),kt,i),obsOfThisTime(1:nn),		&
      	    	    staOfThisTime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)

          title = 'SFC pressure verification of tmf = '//tmtag//' (mb)'
          CALL verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsOfThisTime(1:nn),		&
      	    	    staOfThisTime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	ENDDO


      CASE ("MSLP")
	!Convert pascal to mb
	qc_obs(1,1:numobs(i),i) = qc_obs(1,1:numobs(i),i)/mb2pas

        !time loop
        DO kt = 1,numtmf
	  !Convert pascal to mb
	  analys(1:numgrd(1),1:numgrd(2),kt,i) = analys(1:numgrd(1),1:numgrd(2),kt,i)/mb2pas
	  bkgrnd(1:numgrd(1),1:numgrd(2),kt,i) = bkgrnd(1:numgrd(1),1:numgrd(2),kt,i)/mb2pas

          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          DO j=1,numobs(i)
            IF(INT(qc_obs(4,j,i)).EQ.obstime) THEN
              nn = nn + 1
              obsOfThisTime(nn) = qc_obs(1,j,i)
              staOfThisTime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            ENDIF
          ENDDO

          title = 'MSL pressure background verification of tmf = '//tmtag//' (mb)'
          ea = 0.68
          CALL verify(bkgrnd(1:numgrd(1),1:numgrd(2),kt,i),obsOfThisTime(1:nn),		&
      	    	    staOfThisTime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)

          title = 'MSL pressure verification of tmf = '//tmtag//' (mb)'
          CALL verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsOfThisTime(1:nn),		&
      	    	    staOfThisTime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	ENDDO


      CASE ("PCP1")
	!Convert meter to mm 
	qc_obs(1,1:numobs(i),i) = qc_obs(1,1:numobs(i),i)*1000.0

        !time loop
        DO kt = 1,numtmf
	  !Convert meter to mm
	  analys(1:numgrd(1),1:numgrd(2),kt,i) = analys(1:numgrd(1),1:numgrd(2),kt,i)*1000.0
	  bkgrnd(1:numgrd(1),1:numgrd(2),kt,i) = bkgrnd(1:numgrd(1),1:numgrd(2),kt,i)*1000.0

          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          DO j=1,numobs(i)
            IF(INT(qc_obs(4,j,i)).EQ.obstime) THEN
              nn = nn + 1
              obsOfThisTime(nn) = qc_obs(1,j,i)
              staOfThisTime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            ENDIF
          ENDDO

	  !on precipitation we only verify obs and analysis
          title = '1hr Precipitation verification of tmf = '//tmtag//' (meter)'
          CALL verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsOfThisTime(1:nn),		&
      	    	    staOfThisTime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	ENDDO


      CASE ("PCP3")
	!Convert meter to mm 
	qc_obs(1,1:numobs(i),i) = qc_obs(1,1:numobs(i),i)*1000.0

        !time loop
        DO kt = 1,numtmf
	  !Convert meter to mm
	  analys(1:numgrd(1),1:numgrd(2),kt,i) = analys(1:numgrd(1),1:numgrd(2),kt,i)*1000.0
	  bkgrnd(1:numgrd(1),1:numgrd(2),kt,i) = bkgrnd(1:numgrd(1),1:numgrd(2),kt,i)*1000.0

          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          DO j=1,numobs(i)
            IF(INT(qc_obs(4,j,i)).EQ.obstime) THEN
              nn = nn + 1
              obsOfThisTime(nn) = qc_obs(1,j,i)
              staOfThisTime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            ENDIF
          ENDDO

	  !on precipitation we only verify obs and analysis
          title = '3hr Precipitation verification of tmf = '//tmtag//' (meter)'
          CALL verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsOfThisTime(1:nn),		&
      	    	    staOfThisTime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	ENDDO


      CASE ("PCP6")
	!Convert meter to mm 
	qc_obs(1,1:numobs(i),i) = qc_obs(1,1:numobs(i),i)*1000.0

        !time loop
        DO kt = 1,numtmf
	  !Convert meter to mm
	  analys(1:numgrd(1),1:numgrd(2),kt,i) = analys(1:numgrd(1),1:numgrd(2),kt,i)*1000.0
	  bkgrnd(1:numgrd(1),1:numgrd(2),kt,i) = bkgrnd(1:numgrd(1),1:numgrd(2),kt,i)*1000.0

          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          DO j=1,numobs(i)
            IF(INT(qc_obs(4,j,i)).EQ.obstime) THEN
              nn = nn + 1
              obsOfThisTime(nn) = qc_obs(1,j,i)
              staOfThisTime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            ENDIF
          ENDDO

	  !on precipitation we only verify obs and analysis
          title = '6hr Precipitation verification of tmf = '//tmtag//' (meter)'
          CALL verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsOfThisTime(1:nn),		&
      	    	    staOfThisTime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	ENDDO


      CASE ("PC24")
	!Convert meter to mm 
	qc_obs(1,1:numobs(i),i) = qc_obs(1,1:numobs(i),i)*1000.0

        !time loop
        DO kt = 1,numtmf
	  !Convert meter to mm
	  analys(1:numgrd(1),1:numgrd(2),kt,i) = analys(1:numgrd(1),1:numgrd(2),kt,i)*1000.0
	  bkgrnd(1:numgrd(1),1:numgrd(2),kt,i) = bkgrnd(1:numgrd(1),1:numgrd(2),kt,i)*1000.0

          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          DO j=1,numobs(i)
            IF(INT(qc_obs(4,j,i)).EQ.obstime) THEN
              nn = nn + 1
              obsOfThisTime(nn) = qc_obs(1,j,i)
              staOfThisTime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            ENDIF
          ENDDO

	  !on precipitation we only verify obs and analysis
          title = '24hr Precipitation verification of tmf = '//tmtag//' (meter)'
          CALL verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsOfThisTime(1:nn),		&
      	    	    staOfThisTime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	ENDDO


      CASE DEFAULT
        WRITE(*,22) varnam(i)

    END SELECT
  ENDDO	! end of variable loop

  CLOSE(iunit)

  !Remember to free memory of bkgobs
  !
  DEALLOCATE(bkgobs,STAT=err)
  IF (err .NE. 0) THEN
    PRINT*,'STMAS>STMASVer: cannot deallocate bkgobs memory!'
    STOP
  ENDIF
  DEALLOCATE(stanam,STAT=err)
  IF (err .NE. 0) THEN
    PRINT*,'STMAS>STMASVer: cannot deallocate stanam memory!'
    STOP
  ENDIF


  !everything is done.
  PRINT *,' Normal completion of STMASVer'
  RETURN 

999  PRINT *,'ERROR opening ',ver_file(1:len)
  RETURN

22 FORMAT('STMAS>STMASVer: Do not know to verify such var: ',A4)

END SUBROUTINE STMASVer

SUBROUTINE AddBkgrd
!==========================================================
!  This routine mask out obs covered area defined by radius 
!  and then add background grid data into obs vector as if 
!  we have obs in those areas. 
!
!  HISTORY:
!       Creation: 26-08-2008 by min-ken hsieh
!       Modification:
!                 11-2008 by min-ken hsieh
!                            Using Jb term in STMASAna instead of adding bkg to obs here.
!                            We only mark uncovered areas here.
!==========================================================

  IMPLICIT NONE

  !Local variable
  INTEGER :: i,j,kx,ky,kt,nn
  INTEGER :: FirstCoveredGrid(2),LastCoveredGrid(2)
  INTEGER :: obstime
  LOGICAL :: land(numgrd(1),numgrd(2))		!land mask
  LOGICAL :: sameAsStn(numgrd(1),numgrd(2))	!grid land/sea is the same as stn
  LOGICAL :: stnOverLand

  CHARACTER*1 :: tmtag

  !initialize land array
  DO j = 1,numgrd(2)
    DO i = 1,numgrd(1)
      land(i,j) = (lndfac(i,j).GT.0.0)
    ENDDO
  ENDDO

  !variable loop
  DO i=1,numvar
    ! find out areas have been covered by obs for each time frame
    DO kt=1,numtmf
      uncovr(:,:,kt,i) = .TRUE.
      nn= 0
      obstime = domain(1,3)+(kt-1)*lapsdt
      DO j=1,numobs(i)
        IF(INT(qc_obs(4,j,i)).EQ.obstime) THEN
          nn = nn + 1
          FirstCoveredGrid(1:2) = MAX0(1,FLOOR(qc_obs(2:3,j,i))-radius(i))
          LastCoveredGrid(1:2) = MIN0(numgrd(1:2),FLOOR(qc_obs(2:3,j,i))+radius(i)+1)

          !land/water process
          IF(lndsea(i).EQ.1) THEN
            stnOverLand = (lndfac(INT(qc_obs(2,j,i)),INT(qc_obs(3,j,i))).GT.0.0)
            DO ky=FirstCoveredGrid(2),LastCoveredGrid(2)
	      DO kx=FirstCoveredGrid(1),LastCoveredGrid(1)
                sameAsStn(kx,ky)= (stnOverLand .EQV. land(kx,ky))
              ENDDO
            ENDDO
          ELSE
            sameAsStn = .TRUE.
          ENDIF
    
          ! mask out
          DO ky=FirstCoveredGrid(2),LastCoveredGrid(2)
            DO kx=FirstCoveredGrid(1),LastCoveredGrid(1)
	      uncovr(kx,ky,kt,i) = .NOT.sameAsStn(kx,ky)
            ENDDO
          ENDDO
 
        ENDIF
      ENDDO

    ENDDO

  ENDDO ! end of variable loop
  RETURN

END SUBROUTINE AddBkgrd

SUBROUTINE JbGridpt

!==========================================================
!  This routine identifies those gridpoints over which the
!  background are needed for J_b term in the cost function.
!  Note. 
!  (1) This intends to improve the efficient of Min-Ken's
!  AddBkgrd routine which takes 5 minutes over CONUS domain.
!  (2) Consider spatial only, i.e., no temporal variation
!  as the namelist, stmas_mg.vr, now provides spatial radius
!  of data coverage.
!
!  HISTORY:
!       Creation: Nov, 2008 by Yuanfu Xie.
!==========================================================

  ! Local variables:
  INTEGER :: i,j,ic,jc,kc,iv,jo,is,ie,js,je,km,kp
  LOGICAL :: o,g
  REAL    :: r2,ol

  uncovr = .TRUE.

  !variable loop:
  DO iv=1,numvar
    ! Time frames:
    r2 = radius(iv)*radius(iv)
    DO jo=1,numobs(iv)

      ! Center gridpoint:
      ic = INT(qc_obs(2,jo,iv))
      jc = INT(qc_obs(3,jo,iv))
      kc = (INT(qc_obs(4,jo,iv))-domain(1,3))/lapsdt+1
      km = max0(kc-1,1)
      kp = min0(kc+1,numgrd(3))
      IF ((verbal .EQ. 1) .AND. ((kc .LT. 1) .OR. (kc .GT. numgrd(3)))) THEN
        PRINT*,'Obs out of range: ',iv,jo,qc_obs(4,jo,iv),domain(1,3)
        CYCLE
      ENDIF
      IF (kc .EQ. 1) kc = 2
      IF (kc .EQ. numgrd(3)) kc = kc-1

      ! Landfactor:
      ol = lndfac(ic,jc)+lndfac(ic+1,jc)+lndfac(ic,jc+1)+lndfac(ic+1,jc+1)
      o = .FALSE.
      IF ((ol .GT. 0.0) .OR. (lndsea(iv) .EQ. 0)) o = .TRUE.

      ! Covered circle:
      is = MAX0(ic-radius(iv),1)
      ie = MIN0(ic+radius(iv),numgrd(1))
      js = MAX0(jc-radius(iv),1)
      je = MIN0(jc+radius(iv),numgrd(2))

      DO j=js,je
        DO i=is,ie
          g = .FALSE.
          IF ((lndfac(i,j) .GT. 0.0) .OR. (lndsea(iv) .EQ. 0)) g = .TRUE.
          ! Within the influence radius and landfactors the same:
          IF (((i-ic)*(i-ic)+(j-jc)*(j-jc) .LT. r2) .AND. (o .EQV. g)) &
            uncovr(i,j,km:kp,iv) = .FALSE.
        ENDDO
      ENDDO

    ENDDO
  ENDDO
END SUBROUTINE JbGridpt



SUBROUTINE JbInGaus

!==========================================================
!  This routine identifies those gridpoints over which the
!  background are needed for J_b term in the cost function.
!  Note. 
!  This is a modified version of JbGridpt.
!
!  HISTORY:
!       Creation: Nov, 2008 by Yuanfu Xie.
!==========================================================

  ! Local variables:
  INTEGER :: i,j,ic,jc,kc,iv,jo,is,ie,js,je,iex
  LOGICAL :: o,g
  REAL    :: r2,ol

  diagnl = 1.0

  iex = 10		! Uncover gridpoint with Gaussian decay

  !variable loop:
  DO iv=1,numvar
    ! Time frames:
    r2 = radius(iv)*radius(iv)
    DO jo=1,numobs(iv)

      ! Center gridpoint:
      ic = INT(qc_obs(2,jo,iv))
      jc = INT(qc_obs(3,jo,iv))
      kc = (INT(qc_obs(4,jo,iv))-domain(1,3))/lapsdt+1
      IF ((kc .LT. 1) .OR. (kc .GT. numgrd(3))) THEN
        PRINT*,'Obs out of range: ',iv,jo,qc_obs(4,jo,iv),domain(1,3)
        CYCLE
      ENDIF
      IF (kc .EQ. 1) kc = 2
      IF (kc .EQ. numgrd(3)) kc = kc-1

      ! Landfactor:
      ol = lndfac(ic,jc)*lndfac(ic+1,jc)*lndfac(ic,jc+1)*lndfac(ic+1,jc+1)
      o = .FALSE.
      IF ((ol .GT. 0.0) .OR. (lndsea(iv) .EQ. 0)) o = .TRUE.

      ! Covered circle:
      is = MAX0(ic-radius(iv)-iex,1)
      ie = MIN0(ic+radius(iv)+iex,numgrd(1))
      js = MAX0(jc-radius(iv)-iex,1)
      je = MIN0(jc+radius(iv)+iex,numgrd(2))

      DO j=js,je
        DO i=is,ie
          g = .FALSE.
          IF ((lndfac(i,j) .GT. 0.0) .OR. (lndsea(iv) .EQ. 0)) g = .TRUE.
          ! Within the influence radius and landfactors the same:
          IF (((i-ic)*(i-ic)+(j-jc)*(j-jc) .LT. r2) .AND. (o .EQV. g)) THEN
            diagnl(i,j,iv) = 0.0
          ELSE
            IF (o .EQV. g) diagnl(i,j,iv) = &
              1.0-exp(r2-(i-ic)*(i-ic)-(j-jc)*(j-jc))
          ENDIF
        ENDDO
      ENDDO

    ENDDO
  ENDDO
END SUBROUTINE JbInGaus

END MODULE LAPSDatSrc
