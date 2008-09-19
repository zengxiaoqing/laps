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

PROGRAM STMAS_MG

!==========================================================
!  This program is to analyze surface data using time and
!  space observation through a multigrid technique.
!
!  HISTORY: 
! 	Creation: YUANFU XIE	6-2005
!==========================================================

  USE Definition
  USE LAPSDatSrc
  USE MemoryMngr
  USE PrePostPrc
  USE STMASAnalz

  IMPLICIT NONE

  INTEGER :: i

  ! Namelist:
  CALL PrPstNLs

  ! LAPS parameters:
  CALL LAPSInfo

  ! Allocate memory for LAPS and interactive vars:
  CALL LAPSMemo		! For LAPS usage
  CALL IntrMemo		! For interactive variables

  ! LAPS grid configuration:
  CALL LAPSConf

  ! Background fields:
  CALL LAPSBKGD

  ! Surface LSO obs:
  CALL LAPSOBSV(mxstts)

  ! Convert to units consistent with background:
  CALL LAPSUnit

  ! STMAS QC:
  CALL LAPS_QCs

  ! Add background to obs where obs is spare !by min-ken hseih
  CALL AddBkgrd

  ! Release memory for LAPS:
  CALL LAPSRels

  ! Allocate memory for STMAS:
  CALL STMASMem

  ! Test analytic function:
  ! CALL STMASTst

  ! STMAS Analyses:
  DO i=1,numvar
    WRITE(*,*) 'STMAS_MAIN: Start analyzing ',varnam(i)
    ! Check if there is any obs for analysis:
    IF (numobs(i) .GT. 0) THEN
      ! modified by min-ken hsieh
      ! pass stanam into STMASAna
      !
      CALL STMASAna(analys(1,1,1,i),numgrd,grdspc, &
	domain,bkgrnd(1,1,1,i),numtmf, &
	qc_obs(1,1,i),numobs(i),weight(1,i), stanam(1,i),&
	obsspc(1,i),indice(1,1,i),coeffs(1,1,i), &
	bounds(i),stmasi,stmasr)
    ELSE
      ! No analysis:
      analys(1:numgrd(1),1:numgrd(2),1:numgrd(3),i) = 0.0
    ENDIF
  ENDDO

  ! Add increment to the background:
  CALL STMASInc


  ! Write analyses to LSX
  IF (savdat .EQ. 1) THEN
    WRITE(*,*) numgrd,numvar
    WRITE(11,*) analys
  ENDIF
  CALL PrPstLSX

  ! modified by min-ken hsieh
  ! Verify STMAS Analysis
  CALL STMASVer

  ! Release dynamic memory:
  CALL IntrRels
  CALL STMARels

  ! End of analysis:
  WRITE(*,*) 'STMAS analysis succeeds'

END PROGRAM STMAS_MG
