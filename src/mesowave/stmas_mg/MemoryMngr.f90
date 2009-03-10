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

MODULE MemoryMngr

!==========================================================
!  This module manages memory allocation and deallocation.
!
!  HISTORY:
!	Creation: YUANFU XIE	8-2005
!==========================================================

  USE Definition

CONTAINS

SUBROUTINE LAPSMemo

!==========================================================
!  This routine allocates necessary memory for STMAS usage.
!
!  HISTORY: MAY. 2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER :: err

  ALLOCATE(i4prev(numtmf), &
	   latgrd(numgrd(1),numgrd(2)), &
	   longrd(numgrd(1),numgrd(2)), &
	   topogr(numgrd(1),numgrd(2)), &
	   rawobs(4,numtmf*mxstts,numvar), &
	   bkgobs(numtmf*mxstts,numvar), &
	   stanam(numtmf*mxstts,numvar), &		!Added by min-ken.hsieh: stanam for STMASVer
	   STAT=err)
  IF (err .NE. 0) THEN
    PRINT*,'STMAS>GetMemos: cannot allocate enough memory!'
    STOP
  ENDIF

END SUBROUTINE LAPSMemo

SUBROUTINE LAPSRels

!==========================================================
!  This routine releases necessary memory for STMAS usage.
!
!  HISTORY: MAY. 2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER :: err

!
! modified by min-ken hsieh
! because we still need bkgobs for verify
! we deallocate it in STMASVer subroutine
!
  DEALLOCATE(i4prev, &
	     latgrd, &
	     longrd, &
	     topogr, &
	     rawobs, &
	     STAT=err)

  IF (err .NE. 0) THEN
    PRINT*,'STMAS>RelsMemo: cannot delete allocated memory!'
    STOP
  ENDIF

END SUBROUTINE LAPSRels

SUBROUTINE IntrMemo

!==========================================================
!  This routine allocates necessary memory for interactive
!  variables.
!
!  HISTORY: AUG. 2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER :: err

  ALLOCATE(bkgrnd(numgrd(1),numgrd(2),numtmf,numvar), &
           uncovr(numgrd(1),numgrd(2),numtmf,numvar), &	! Add by min-ken, uncovered used in 
	   diagnl(numgrd(1),numgrd(2),numvar), &	! Add B's diagonal array for J_b term
	   lndfac(numgrd(1),numgrd(2)), &		! AddBkgrd and STMASAna
	   qc_obs(4,numtmf*mxstts,numvar), &
	   weight(numtmf*mxstts,numvar), &
	   indice(6,numtmf*mxstts,numvar), &
	   coeffs(6,numtmf*mxstts,numvar), &
	   STAT=err)
  IF (err .NE. 0) THEN
    PRINT*,'STMAS>IntrMemo: cannot allocate memory!'
    STOP
  ENDIF

END SUBROUTINE IntrMemo

SUBROUTINE IntrRels

!==========================================================
!  This routine releases necessary memory for interactive
!  variables.
!
!  HISTORY: AUG. 2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER :: err

  DEALLOCATE(bkgrnd, &
	     uncovr, &
             diagnl, &
	     lndfac, &
	     qc_obs, &
	     weight, &
	     indice, &
	     coeffs, &
	     STAT=err)
  IF (err .NE. 0) THEN
    PRINT*,'STMAS>IntrMemo: cannot release memory!'
    STOP
  ENDIF

END SUBROUTINE IntrRels

SUBROUTINE STMASMem

!==========================================================
!  This routine allocates necessary memory for STMAS 
!  analysis variables.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER :: err

  ALLOCATE(analys(numgrd(1),numgrd(2),numgrd(3),numvar), &
	   STAT=err)
  IF (err .NE. 0) THEN
    PRINT*,'STMAS>STMASMem: cannot allocate memory!'
    STOP
  ENDIF

END SUBROUTINE STMASMem

SUBROUTINE STMARels

!==========================================================
!  This routine allocates necessary memory for STMAS 
!  analysis variables.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER :: err

  DEALLOCATE(analys, &
	     STAT=err)
  IF (err .NE. 0) THEN
    PRINT*,'STMAS>STMASMem: cannot release memory!'
    STOP
  ENDIF

END SUBROUTINE STMARels

END MODULE MemoryMngr
