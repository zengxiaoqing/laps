MODULE INPUT_BG_OBS

  USE PRMTRS_STMAS
  USE GENERALTOOLS

  IMPLICIT NONE
  INTEGER(KIND=4) ,  PARAMETER :: LURAO=11,LUNDX=21,LUNEW=50,NOBSMAX=1000000
  REAL(KIND=8)                 :: DX0,DY0
  INTEGER(KIND=4) ,ALLOCATABLE :: STT(:)
  INTEGER(KIND=4) ,ALLOCATABLE :: GRIDMASK(:,:,:,:,:)
  REAL(KIND=4)    ,ALLOCATABLE :: X00(:,:),Y00(:,:),P00(:),Z00(:,:,:,:)
  REAL(KIND=8)    ,ALLOCATABLE :: TPG(:,:),BK0(:,:,:,:,:),C00(:,:),D00(:,:),S00(:,:,:,:)
  REAL(KIND=8)    ,ALLOCATABLE :: OBRADIUS(:,:)
  REAL(KIND=4)    ,ALLOCATABLE :: OBP(:,:),OBS(:),OBE(:)

CONTAINS

SUBROUTINE BKGRNDOBS
!*************************************************
! MAIN ROUTINE OF THIS PREPARATION CODE
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
  PRINT*,'READNMLST'
  CALL READNMLST
  PRINT*,'BKGMEMALC'
  CALL BKGMEMALC
  PRINT*,'MEMORYALC'
  CALL MEMORYALC
  PRINT*,'RDBCKGRND'
!  CALL RDBCKGRND
  CALL RDLAPSBKG
  PRINT*,'GETBKGRND'
  CALL GETBKGRND
  PRINT*,'RDBUFROBS'
  CALL RDBUFROBS
  IF(IFBKGND.EQ.1)THEN
    PRINT*,'ADDBKGRND'
    CALL ADDBKGRND
  ENDIF
  PRINT*,'OBSMEMALC'
  CALL OBSMEMALC
  PRINT*,'MEMORYRLS'
  CALL MEMORYRLS
  RETURN
END SUBROUTINE BKGRNDOBS

SUBROUTINE READNMLST
!*************************************************
! READ NAME LIST
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER(KIND=4) :: S,N,ER,NU
! --------------------
  NU=2
  OPEN(NU,FILE='namelist.txt',ACTION='READ',STATUS='OLD')
  READ(NU,*)IFBKGND
  READ(NU,*)IFBOUND
  READ(NU,*)IFPCDNT
  READ(NU,*)NUMSTAT
  ALLOCATE(SL0(NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'SL0 ALLOCATE WRONG'
  DO S=1,NUMSTAT
    READ(NU,*)SL0(S)
  ENDDO
  ALLOCATE(PENAL0X(NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'PENAL0X ALLOCATE WRONG'
  ALLOCATE(PENAL0Y(NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'PENAL0Y ALLOCATE WRONG'
  ALLOCATE(PENAL0Z(NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'PENAL0Z ALLOCATE WRONG'
  ALLOCATE(PENAL0T(NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'PENAL0T ALLOCATE WRONG'
  ALLOCATE(PENAL_X(NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'PENAL_X ALLOCATE WRONG'
  ALLOCATE(PENAL_Y(NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'PENAL_Y ALLOCATE WRONG'
  ALLOCATE(PENAL_Z(NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'PENAL_Z ALLOCATE WRONG'
  ALLOCATE(PENAL_T(NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'PENAL_T ALLOCATE WRONG'
  DO S=1,NUMSTAT
    READ(NU,*)PENAL0X(S)
    READ(NU,*)PENAL0Y(S)
    READ(NU,*)PENAL0Z(S)
    READ(NU,*)PENAL0T(S)
  ENDDO
  READ(NU,*)PNLT0PU
  READ(NU,*)PNLT0PV
  READ(NU,*)NUMDIMS
  READ(NU,*)NUMGRID(1)
  READ(NU,*)NUMGRID(2)
  READ(NU,*)NUMGRID(3)
  READ(NU,*)NUMGRID(4)
  READ(NU,*)GRDSPAC(3)
  READ(NU,*)GRDSPAC(4)
  READ(NU,*)ORIPSTN(1)
  READ(NU,*)ORIPSTN(2)
  READ(NU,*)ORIPSTN(3)
  READ(NU,*)ORIPSTN(4)
  READ(NU,*)MAXGRID(1)
  READ(NU,*)MAXGRID(2)
  READ(NU,*)MAXGRID(3)
  READ(NU,*)MAXGRID(4)
  READ(NU,*)FNSTGRD
  ALLOCATE(OBRADIUS(MAXDIMS,NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'OBRADIUS ALLOCATE WRONG'
  DO S=1,NUMSTAT
    DO N=1,MAXDIMS
      READ(NU,*)OBRADIUS(N,S)
    ENDDO
  ENDDO
  READ(NU,*)U_CMPNNT
  READ(NU,*)V_CMPNNT
  READ(NU,*)W_CMPNNT
  READ(NU,*)PRESSURE
  READ(NU,*)TEMPRTUR
  READ(NU,*)XSL
  READ(NU,*)YSL
  READ(NU,*)PSL
  READ(NU,*)CSL
  READ(NU,*)DSL
  READ(NU,*)COSSTEP
  READ(NU,*)MIDGRID
  READ(NU,*)FINSTEP
  CLOSE(NU)
  NGPTOBS=2**NUMDIMS
  RETURN
END SUBROUTINE READNMLST

SUBROUTINE BKGMEMALC
!*************************************************
! ALLOCATE MEMORY FOR BACKGROUND ARRAY
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER(KIND=4) :: ER
! --------------------
  ALLOCATE(GRDBKGD0(MAXGRID(1),MAXGRID(2),MAXGRID(3),MAXGRID(4),NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'GRDBKGD0 ALLOCATE WRONG'
  ALLOCATE(XX0(MAXGRID(1),MAXGRID(2)),STAT=ER)
  IF(ER.NE.0)STOP 'XX0 ALLOCATE WRONG'
  ALLOCATE(YY0(MAXGRID(1),MAXGRID(2)),STAT=ER)
  IF(ER.NE.0)STOP 'YY0 ALLOCATE WRONG'
  ALLOCATE(CR0(MAXGRID(1),MAXGRID(2)),STAT=ER)
  IF(ER.NE.0)STOP 'CR0 ALLOCATE WRONG'
  ALLOCATE(DG0(MAXGRID(1),MAXGRID(2)),STAT=ER)
  IF(ER.NE.0)STOP 'DG0 ALLOCATE WRONG'
  ALLOCATE(DN0(MAXGRID(1),MAXGRID(2),MAXGRID(3),MAXGRID(4)),STAT=ER)
  IF(ER.NE.0)STOP 'DN0 ALLOCATE WRONG'
  IF(IFPCDNT.EQ.0)THEN
    ALLOCATE(ZZ0(MAXGRID(1),MAXGRID(2),MAXGRID(3),MAXGRID(4)),STAT=ER)
    IF(ER.NE.0)STOP 'ZZ0 ALLOCATE WRONG'
  ELSEIF(IFPCDNT.EQ.1)THEN
    ALLOCATE(PP0(MAXGRID(3)),STAT=ER)
    IF(ER.NE.0)STOP 'PP0 ALLOCATE WRONG'
  ENDIF
  RETURN
END SUBROUTINE BKGMEMALC

SUBROUTINE OBSMEMALC
!*************************************************
! ALLOCATE OBSERVATION MEMORY AND READ IN DATA AND SCALE
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER(KIND=4) :: N,O,ER
! --------------------
  IF(NALLOBS.EQ.0)RETURN
  PRINT*,'NALLOBS=',NALLOBS
  ALLOCATE(OBSPOSTN(NUMDIMS,NALLOBS),STAT=ER)
  IF(ER.NE.0)STOP 'OBSPOSTN ALLOCATE WRONG'
  ALLOCATE(OBSCOEFF(NGPTOBS,NALLOBS),STAT=ER)
  IF(ER.NE.0)STOP 'OBSCOEFF ALLOCATE WRONG'
  ALLOCATE(OBSIDXPC(MAXDIMS,NALLOBS),STAT=ER)
  IF(ER.NE.0)STOP 'OBSIDXPC ALLOCATE WRONG'
  ALLOCATE(OBSVALUE(NALLOBS),STAT=ER)
  IF(ER.NE.0)STOP 'OBSVALUE ALLOCATE WRONG'
  ALLOCATE(OBSERROR(NALLOBS),STAT=ER)
  IF(ER.NE.0)STOP 'OBSERROR ALLOCATE WRONG'
  ALLOCATE(OBSSTATE(NALLOBS),STAT=ER)
  IF(ER.NE.0)STOP 'OBSSTATE ALLOCATE WRONG'
!  ALLOCATE(OBSEINF1(NALLOBS),STAT=ER)
!  IF(ER.NE.0)STOP 'OBSEINF1 ALLOCATE WRONG'
!  ALLOCATE(OBSEINF2(NALLOBS),STAT=ER)
!  IF(ER.NE.0)STOP 'OBSEINF2 ALLOCATE WRONG'
!  ALLOCATE(OBSEINF3(NALLOBS),STAT=ER)
!  IF(ER.NE.0)STOP 'OBSEINF3 ALLOCATE WRONG'
  OPEN(111,FILE='TMPOBS_U.DAT')
  OPEN(222,FILE='TMPOBS_V.DAT')
  OPEN(333,FILE='TMPOBS_Z.DAT')
  OPEN(444,FILE='TMPOBS_T.DAT')
  DO O=1,NALLOBS
    DO N=1,NUMDIMS
      OBSPOSTN(N,O)=OBP(N,O)
    ENDDO
    OBSVALUE(O)=OBS(O)
    OBSERROR(O)=OBE(O)
    OBSSTATE(O)=STT(O)
    IF(STT(O).EQ.1)WRITE(111,'(5F15.4,I3)')(OBP(N,O),N=1,NUMDIMS),OBS(O),OBE(O),STT(O)
    IF(STT(O).EQ.2)WRITE(222,'(5F15.4,I3)')(OBP(N,O),N=1,NUMDIMS),OBS(O),OBE(O),STT(O)
    IF(STT(O).EQ.3)WRITE(333,'(5F15.4,I3)')(OBP(N,O),N=1,NUMDIMS),OBS(O),OBE(O),STT(O)
    IF(STT(O).EQ.4)WRITE(444,'(5F15.4,I3)')(OBP(N,O),N=1,NUMDIMS),OBS(O),OBE(O),STT(O)
  ENDDO
  CLOSE(111)
  CLOSE(222)
  CLOSE(333)
  CLOSE(444)
  RETURN
END SUBROUTINE OBSMEMALC

SUBROUTINE MEMORYALC
!*************************************************
! MEMORY ALLOCATE
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER(KIND=4) :: N,ER
  INTEGER :: ST		! STATUS
! --------------------
  DO N=1,MAXDIMS
    FCSTGRD(N)=1
  ENDDO
  ! OPEN(11,FILE='fort.11',STATUS='OLD',ACTION='READ')
  ! READ(11,*)FCSTGRD(1:NUMDIMS)
  ! PRINT*,FCSTGRD(1),FCSTGRD(2),FCSTGRD(3)
  ! SPATIAL X Y DIMENSIONS:
  CALL GET_GRID_DIM_XY(FCSTGRD(1),FCSTGRD(2),ST)
  CALL GET_LAPS_DIMENSIONS(FCSTGRD(3),ST)
  GRDSPAC(1)=((FCSTGRD(1)-1)*1.0)/((NUMGRID(1)-1)*1.0)
  GRDSPAC(2)=((FCSTGRD(2)-1)*1.0)/((NUMGRID(2)-1)*1.0)
  ALLOCATE(BK0(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),FCSTGRD(4),NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'BK0 ALLOCATE WRONG'
  ALLOCATE(X00(FCSTGRD(1),FCSTGRD(2)),STAT=ER)
  IF(ER.NE.0)STOP 'X00 ALLOCATE WRONG'
  ALLOCATE(Y00(FCSTGRD(1),FCSTGRD(2)),STAT=ER)
  IF(ER.NE.0)STOP 'Y00 ALLOCATE WRONG'
  ALLOCATE(TPG(FCSTGRD(1),FCSTGRD(2)),STAT=ER)
  IF(ER.NE.0)STOP 'TPG ALLOCATE WRONG'
  ALLOCATE(C00(FCSTGRD(1),FCSTGRD(2)),STAT=ER)
  IF(ER.NE.0)STOP 'C00 ALLOCATE WRONG'
  ALLOCATE(D00(FCSTGRD(1),FCSTGRD(2)),STAT=ER)
  IF(ER.NE.0)STOP 'D00 ALLOCATE WRONG'
  ALLOCATE(P00(FCSTGRD(3)),STAT=ER)
  IF(ER.NE.0)STOP 'P00 ALLOCATE WRONG'
  ALLOCATE(S00(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),FCSTGRD(4)),STAT=ER)
  IF(ER.NE.0)STOP 'S00 ALLOCATE WRONG'
  ALLOCATE(GRIDMASK(FCSTGRD(1),FCSTGRD(2),MAXGRID(3),FCSTGRD(4),NUMSTAT),STAT=ER) !!!!!ATTENTION
  IF(ER.NE.0)STOP 'GRIDMASK ALLOCATE WRONG'
  ALLOCATE(OBP(NUMDIMS,NOBSMAX),STAT=ER)
  IF(ER.NE.0)STOP 'OBP ALLOCATE WRONG'
  ALLOCATE(OBS(NOBSMAX),STAT=ER)
  IF(ER.NE.0)STOP 'OBS ALLOCATE WRONG'
  ALLOCATE(OBE(NOBSMAX),STAT=ER)
  IF(ER.NE.0)STOP 'OBE ALLOCATE WRONG'
  ALLOCATE(STT(NOBSMAX),STAT=ER)
  IF(ER.NE.0)STOP 'STT ALLOCATE WRONG'
  RETURN
END SUBROUTINE MEMORYALC

SUBROUTINE MEMORYRLS
!*************************************************
! MEMORY RELEASE
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
  DEALLOCATE(BK0)
  DEALLOCATE(X00)
  DEALLOCATE(Y00)
  DEALLOCATE(TPG)
  DEALLOCATE(C00)
  DEALLOCATE(D00)
  DEALLOCATE(P00)
  DEALLOCATE(S00)
  DEALLOCATE(OBP)
  DEALLOCATE(OBS)
  DEALLOCATE(OBE)
  DEALLOCATE(STT)
  DEALLOCATE(GRIDMASK)
  DEALLOCATE(OBRADIUS)
  RETURN
END SUBROUTINE MEMORYRLS

SUBROUTINE RDLAPSBKG
!*************************************************
! READ IN BACKGROUND
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!	   OCTOBER   2007, by YUANFU XIE (USE LAPS)
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER(KIND=4) :: I,J,K
  REAL(KIND=8) :: UU(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3))
  REAL(KIND=8) :: VV(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3))
  REAL(KIND=8) :: ZZ(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3))
  REAL(KIND=8) :: TT(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3))
  REAL(KIND=8) :: OX,OY,EX,EY

  CHARACTER*4  :: VN		! VARNAME ADDED BY YUANFU
  CHARACTER*9  :: TM		! ASCII TIME ADDED BY YUANFU
  INTEGER      :: ST		! STATUS OF LAPS CALLS ADDED BY YUANFU
  INTEGER*4    :: I4		! SYSTEM I4 TIME ADDED BY YUANFU
  REAL         :: HT(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3)), &
                  T3(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3)), &
                  U3(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3)), &
                  V3(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3)), &
                  LA(FCSTGRD(1),FCSTGRD(2)),LO(FCSTGRD(1),FCSTGRD(2)), &
                  AG(FCSTGRD(1),FCSTGRD(2)),P1(FCSTGRD(3)),DS(2)
! --------------------
!  OPEN(11,FILE='fort.11',STATUS='OLD',ACTION='READ')
!  READ(11,*)
!  READ(11,*)ZZ,TT,UU,VV
!  READ(11,*)Y00,X00,TPG,P00
!  CLOSE(11)
  ! USE LAPS INGEST TO READ BACKGROUND FIELDS: ADDED BY YUANFU

  ! SYSTEM TIME:
  CALL GET_SYSTIME(I4,TM,ST)

  ! HEIGHT FIELD:
  VN = 'HT'
  CALL GET_MODELFG_3D(I4,VN,FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),HT,ST)
  ! TEMPERATURE FIELD:
  VN = 'T3'
  CALL GET_MODELFG_3D(I4,VN,FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),T3,ST)
  ! WIND U:
  VN = 'U3'
  CALL GET_MODELFG_3D(I4,VN,FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),U3,ST)
  ! WIND V:
  VN = 'V3'
  CALL GET_MODELFG_3D(I4,VN,FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),V3,ST)
  ! DOMAIN CONFIGURATION:
  CALL READ_STATIC_GRID(FCSTGRD(1),FCSTGRD(2),'LAT',&
                        LA,ST)
  CALL READ_STATIC_GRID(FCSTGRD(1),FCSTGRD(2),'LON',&
                        LO,ST)
  CALL READ_STATIC_GRID(FCSTGRD(1),FCSTGRD(2),'AVG',&
                        AG,ST)
  ! PRESSURE LEVELS:
  CALL GET_PRES_1D(I4,FCSTGRD(3),P1,ST)

  ! GRID SPACING:
  CALL GET_GRID_SPACING_ACTUAL(LA((FCSTGRD(1)-1)/2+1, &
                                  (FCSTGRD(2)-1)/2+1), &
			       LO((FCSTGRD(1)-1)/2+1, &
                                  (FCSTGRD(2)-1)/2+1),DS,ST)
  DX0=DS(1)
  DY0=DX0 !DS(2)

  ! CONVERT TO REAL 8:
  ZZ = HT
  TT = T3
  UU = U3
  VV = V3
  Y00 = LA
  X00 = LO
  TPG = AG
  P00 = P1

  ! END OF YUANFU'S MODIFICATION

  DO K=1,FCSTGRD(3)
  DO J=1,FCSTGRD(2)
  DO I=1,FCSTGRD(1)
    S00(I,J,K,1)=1.0D0
  ENDDO
  ENDDO
  ENDDO
  DO J=1,FCSTGRD(2)
  DO I=1,FCSTGRD(1)
    D00(I,J)=0.0D0
    C00(I,J)=2.0*7.29E-5*SIND(Y00(I,J))
  ENDDO
  ENDDO
  DO K=1,FCSTGRD(3)
  DO J=1,FCSTGRD(2)
  DO I=1,FCSTGRD(1)
    BK0(I,J,K,1,1)=UU(I,J,K)
    BK0(I,J,K,1,2)=VV(I,J,K)
    BK0(I,J,K,1,3)=ZZ(I,J,K)
    BK0(I,J,K,1,4)=TT(I,J,K)
  ENDDO
  ENDDO
  ENDDO
  OPEN(20,FILE='VERTICAL.DAT',ACTION='WRITE')
  DO K=1,FCSTGRD(3)
    WRITE(20,*)P00(K)
  ENDDO
  CLOSE(20)
  OX=X00(1,1)
  EX=X00(FCSTGRD(1),FCSTGRD(2))
  OY=Y00(1,1)
  EY=Y00(FCSTGRD(1),FCSTGRD(2))
  CALL DRCONTOUR(OX,EX,OY,EY,NUMDIMS,FCSTGRD,UU,'UB.GRD')
  CALL DRCONTOUR(OX,EX,OY,EY,NUMDIMS,FCSTGRD,VV,'VB.GRD')
  CALL DRCONTOUR(OX,EX,OY,EY,NUMDIMS,FCSTGRD,TT,'TB.GRD')
  CALL DRCONTOUR(OX,EX,OY,EY,NUMDIMS,FCSTGRD,ZZ,'ZB.GRD')
  CALL DRCONTOUR_2D(OX,EX,OY,EY,FCSTGRD(1),FCSTGRD(2),TPG,'GG.GRD')
  RETURN
END SUBROUTINE RDLAPSBKG

SUBROUTINE RDBCKGRND
!*************************************************
! READ IN BACKGROUND
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER(KIND=4) :: I,J,K
  REAL(KIND=8) :: UU(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3))
  REAL(KIND=8) :: VV(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3))
  REAL(KIND=8) :: ZZ(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3))
  REAL(KIND=8) :: TT(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3))
  REAL(KIND=8) :: OX,OY,EX,EY
! --------------------
  OPEN(11,FILE='fort.11',STATUS='OLD',ACTION='READ')
  READ(11,*)
  READ(11,*)ZZ,TT,UU,VV
  READ(11,*)Y00,X00,TPG,P00
  CLOSE(11)
  DO K=1,FCSTGRD(3)
  DO J=1,FCSTGRD(2)
  DO I=1,FCSTGRD(1)
    S00(I,J,K,1)=1.0D0
  ENDDO
  ENDDO
  ENDDO
  DO J=1,FCSTGRD(2)
  DO I=1,FCSTGRD(1)
    D00(I,J)=0.0D0
    C00(I,J)=2.0*7.29E-5*SIND(Y00(I,J))
  ENDDO
  ENDDO
  DO K=1,FCSTGRD(3)
  DO J=1,FCSTGRD(2)
  DO I=1,FCSTGRD(1)
    BK0(I,J,K,1,1)=UU(I,J,K)
    BK0(I,J,K,1,2)=VV(I,J,K)
    BK0(I,J,K,1,3)=ZZ(I,J,K)
    BK0(I,J,K,1,4)=TT(I,J,K)
  ENDDO
  ENDDO
  ENDDO
  OPEN(20,FILE='VERTICAL.DAT',ACTION='WRITE')
  DO K=1,FCSTGRD(3)
    WRITE(20,*)P00(K)
  ENDDO
  CLOSE(20)
  OX=X00(1,1)
  EX=X00(FCSTGRD(1),FCSTGRD(2))
  OY=Y00(1,1)
  EY=Y00(FCSTGRD(1),FCSTGRD(2))
  CALL DRCONTOUR(OX,EX,OY,EY,NUMDIMS,FCSTGRD,UU,'UB.GRD')
  CALL DRCONTOUR(OX,EX,OY,EY,NUMDIMS,FCSTGRD,VV,'VB.GRD')
  CALL DRCONTOUR(OX,EX,OY,EY,NUMDIMS,FCSTGRD,TT,'TB.GRD')
  CALL DRCONTOUR(OX,EX,OY,EY,NUMDIMS,FCSTGRD,ZZ,'ZB.GRD')
  CALL DRCONTOUR_2D(OX,EX,OY,EY,FCSTGRD(1),FCSTGRD(2),TPG,'GG.GRD')
  RETURN
END SUBROUTINE RDBCKGRND

SUBROUTINE GETBKGRND
!*************************************************
! GET BACKGROUND FOR ANALYSIS
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER(KIND=4) :: I,J,K,T,S
  REAL(KIND=8)    :: XB(MAXGRID(1)),YB(MAXGRID(2)),ZB(MAXGRID(3)),TB(MAXGRID(4))
  REAL(KIND=8)    :: XF(FCSTGRD(1)),YF(FCSTGRD(2)),ZF(FCSTGRD(3)),TF(FCSTGRD(4))
  INTEGER(KIND=4) :: FG(MAXDIMS),MG(MAXDIMS)
  REAL(KIND=8)    :: Z1(1),T1(1),Z2(1),T2(1),DX,DY
! --------------------
  DO I=1,FCSTGRD(1)
    XF(I)=(I-1)*1.0D0
  ENDDO
  DO J=1,FCSTGRD(2)
    YF(J)=(J-1)*1.0D0
  ENDDO
  DX=((FCSTGRD(1)-1)*1.0D0)/((MAXGRID(1)-1)*1.0D0)
  DO I=1,MAXGRID(1)
    XB(I)=XF(1)+(I-1)*DX
  ENDDO
  DY=((FCSTGRD(2)-1)*1.0D0)/((MAXGRID(2)-1)*1.0D0)
  DO J=1,MAXGRID(2)
    YB(J)=YF(1)+(J-1)*DY
  ENDDO
  DO K=1,FCSTGRD(3)
    ZF(K)=P00(K)
  ENDDO
  OPEN(2,FILE='P_LEVEL.DAT',STATUS='OLD',ACTION='READ')
  DO K=1,MAXGRID(3)
    READ(2,*)ZB(K)
    PP0(K)=ZB(K)
  ENDDO
  CLOSE(2)
  DO T=1,FCSTGRD(4)
    TF(T)=0.0D0
  ENDDO
  DO T=1,MAXGRID(4)
    TB(T)=0.0D0
  ENDDO
  CALL FCST2BKGD(NUMDIMS,NGPTOBS,NUMSTAT,FCSTGRD,XF,YF,ZF,TF,MAXGRID,XB,YB,ZB,TB,BK0,GRDBKGD0)
  CALL FCST2BKGD(NUMDIMS,NGPTOBS,1,FCSTGRD,XF,YF,ZF,TF,MAXGRID,XB,YB,ZB,TB,S00,DN0)
  FG(1)=FCSTGRD(1)
  FG(2)=FCSTGRD(2)
  FG(3)=1
  FG(4)=1
  Z1(1)=0.0
  T1(1)=0.0
  MG(1)=MAXGRID(1)
  MG(2)=MAXGRID(2)
  MG(3)=1
  MG(4)=1
  Z2(1)=0.0
  T2(1)=0.0
  CALL FCST2BKGD(2,4,1,FG,XF,YF,Z1,T1,MG,XB,YB,Z2,T2,C00,CR0)
  CALL FCST2BKGD(2,4,1,FG,XF,YF,Z1,T1,MG,XB,YB,Z2,T2,D00,DG0)
  DX=((FCSTGRD(1)-1)*DX0)/((MAXGRID(1)-1)*1.0D0)
  DY=((FCSTGRD(2)-1)*DY0)/((MAXGRID(2)-1)*1.0D0)
  DO I=1,MAXGRID(1)
  DO J=1,MAXGRID(2)
    XX0(I,J)=0.0D0+(I-1)*DX
    YY0(I,J)=0.0D0+(J-1)*DY
  ENDDO
  ENDDO
  RETURN
END SUBROUTINE GETBKGRND

SUBROUTINE RDBUFROBS
!*************************************************
! READ IN OBSERVATION, MODIFIED FROM 'raob2dwl.f'
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  REAL(KIND=8),PARAMETER :: MS=10e10
  INTEGER(KIND=4),PARAMETER :: NH=4,NR=15,LN=255
  CHARACTER(LEN=8)       :: SS
  INTEGER(KIND=4)        :: L,NL,IR,DT,NC,NW,O,OS,IS,IP
  REAL(KIND=8)           :: AD,BD,HD(NH),P1
  REAL(KIND=8)           :: RA(NR,LN)
  REAL(KIND=8)           :: X,Y,P,T,UU,VV,ZZ,TT,QQ,UE,VE,ZE,TE,QE
  REAL(KIND=4)           :: OP(NUMDIMS),OB,OE
  INTEGER     , EXTERNAL :: IREADSB,IREADMG,I4DY

  ! VARIABLES FOR LAPS OBS: BY YUANFU
  CHARACTER*150          :: OD
  CHARACTER*9            :: A9
  INTEGER                :: LD,I4,N4,ST

  CALL GET_SYSTIME(I4,A9,ST)
  CALL GET_FILESPEC('bufr',2,OD,ST)
  CALL GET_FILE_TIME(OD,I4,N4)
  CALL S_LEN(OD,LD)
  CALL MAKE_FNAM_LP(N4,A9,ST)
  OD(LD+4:LD+9) = OD(LD-4:LD)
  OD(LD-5:LD+3) = A9
  LD = LD+9

  !OPEN(UNIT=LURAO,FILE='072711900.bufr',FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
  OPEN(UNIT=LURAO,FILE=OD(1:LD),FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')

  ! END LAPS OBS INGEST BY YUANFU

!  OPEN(UNIT=LUNDX,FILE='prepobs_prep.bufrtable',STATUS='OLD',ACTION='READ')
!  OPEN(UNIT=LUNEW,FILE='readresult.dat')
  CALL OPENBF(LURAO,'IN',LURAO)
  NC=0
  NW=0
  CALL GET_CONFIG(IS)
  IF(IS.NE.1)STOP 'LAPS PARAMETERS ARE WRONG!!!'
  O=0
  DO WHILE(IREADMG(LURAO,SS,DT).EQ.0)
    DO WHILE(IREADSB(LURAO).EQ.0)
      CALL UFBINT(LURAO,HD,NH,1,IR,'XOB YOB ELV DHR')
      CALL UFBINT(LURAO,RA,NR,LN,NL, &
                 'XDR YDR PRSS POB POE HRDR UOB VOB WOE ZOB ZOE TOB TOE QOB QOE')
      DO L=1,NL
        NW=NW+1
! FOR X LOCATION
        X=HD(1)
        IF(RA(1,L).LT.MS)X=RA(1,L)
        IF(X.GE.MS)CYCLE
        IF(X.LT.0.0)X=X+360.0
! FOR Y LOCATION
        Y=HD(2)
        IF(RA(2,L).LT.MS)Y=RA(2,L)
        IF(Y.GE.MS)CYCLE
! FOR T LOCATION
        T=HD(4)
        IF(RA(6,L).LT.MS)T=RA(6,L)
!        IF(T.GE.MS)CYCLE
        AD=I4DY(DT)
        CALL GETOBDATE(AD,T,BD)
        T=BD
! FOR P LOCATION
        IF(RA(3,L).GE.MS.AND.RA(4,L).GE.MS)P=MS !(RA(4,L).GE.MS.OR.RA(5,L).GE.MS)
        IF(RA(4,L).LT.MS)P=RA(4,L) !.AND.RA(5,L).LT.MS
        IF(RA(3,L).LT.MS)P=RA(3,L)
        IF(RA(3,L).LT.MS.AND.RA(4,L).LT.MS.AND. &
        DABS(RA(3,L)-RA(4,L)).GT.0.001)STOP 'WRONG' !(RA(4,L).LT.MS.AND.RA(5,L).LT.MS)
        IF(P.LT.MS)P=P*100.0D0
! FOR ZZ OBSERVATION
        ZZ=HD(3)
        IF(RA(10,L).LT.MS)ZZ=RA(10,L)
        ZE=RA(11,L)
!        IF(ZE.GE.MS)ZZ=MS
! CHECK P AND ZZ
        IF(P.GE.MS.AND.ZZ.GE.MS)CYCLE
! TRANSFORM
        OP(1)=X
        OP(2)=Y
        IF(P.LT.MS.AND.ZZ.LT.MS)THEN
          CALL HT_TO_PRS(OP(1),OP(2),ZZ,P1,IS)
          IF(IS.EQ.1)WRITE(555,*)OP(1),OP(2),P,P1,P-P1,ZZ
        ENDIF
        IP=0
        IF(P.GE.MS)THEN
          IP=1
          CALL HT_TO_PRS(OP(1),OP(2),ZZ,P,IS)
          IF(IS.NE.1)CYCLE
        ENDIF
        OP(3)=P
! FOR U COMPONENT OBSERVATION
        UU=RA(7,L)
        UE=RA(9,L)
! FOR V COMPONENT OBSERVATION
        VV=RA(8,L)
        VE=RA(9,L)
! FOR TEMPERATURE OBSERVATION
        TT=RA(12,L)
        TE=RA(13,L)
        IF(TT.LT.MS)TT=TT+273.15D0
! FOR SPECIFIC HUMIDITY OBSERVATION
        QQ=RA(14,L)
        QE=RA(15,L)
! OUTPUT
        UE=1.0
        VE=1.0
        ZE=2.0
        TE=0.5
        IF(UU.LT.MS.AND.UE.LT.MS)THEN
          OB=UU
          OE=UE
          OS=1
          CALL HANDLEOBS(OP,OB,OE,OS,O,IP)
        ENDIF
        IF(VV.LT.MS.AND.VE.LT.MS)THEN
          OB=VV
          OE=VE
          OS=2
          CALL HANDLEOBS(OP,OB,OE,OS,O,IP)
        ENDIF
        IF(ZZ.LT.MS.AND.ZE.LT.MS)THEN
          OB=ZZ
          OE=ZE
          OS=3
          CALL HANDLEOBS(OP,OB,OE,OS,O,IP)
        ENDIF
        IF(TT.LT.MS.AND.TE.LT.MS)THEN
          OB=TT
          OE=TE
          OS=4
          CALL HANDLEOBS(OP,OB,OE,OS,O,IP)
        ENDIF
        NC=NC+1
      ENDDO
    ENDDO
  ENDDO
  NALLOBS=O
  PRINT*,'THE NUMBER OF LOCATION IS',NW,'AND',NC,'AVAILABLE'
  PRINT*,'NALLOBS=',NALLOBS
  CALL CLOSBF(LURAO)
  CLOSE(LURAO)
!  CLOSE(LUNDX)
!  CLOSE(LUNEW)
  RETURN
END SUBROUTINE RDBUFROBS

SUBROUTINE HANDLEOBS(OP,OB,OE,OS,O,IP)
!*************************************************
! CALCULATE THE DIFFERENCE BETWEEN OBSERVATION AND BACKGROUND
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER(KIND=4) :: I,J,K,T,M,N,NP(MAXDIMS),NN(MAXDIMS),IS,O,OS,IP
  REAL(KIND=4) :: X,Y,P
  REAL(KIND=8) :: AC(NUMDIMS,NGPTOBS),OC(NUMDIMS),CO(NGPTOBS),HT
  REAL(KIND=4) :: OP(NUMDIMS),OB,OE
! --------------------
  CALL LATLON_TO_RLAPSGRID(OP(2),OP(1),Y00,X00,FCSTGRD(1),FCSTGRD(2),X,Y,IS)
  IF(X.LT.1.0.OR.Y.LT.1.0.OR.X.GT.FCSTGRD(1).OR.Y.GT.FCSTGRD(2).OR.IS.NE.1)RETURN
  CALL VRTCLPSTN(FCSTGRD(3),P00,OP(3),P,IS)
  IF(IS.NE.1)RETURN
  DO N=1,MAXDIMS
    NP(N)=1
  ENDDO
  NP(1)=INT(X)
  NP(2)=INT(Y)
  NP(3)=INT(P)
  DO N=1,MAXDIMS
    IF(NP(N).EQ.FCSTGRD(N).AND.FCSTGRD(N).NE.1)NP(N)=FCSTGRD(N)-1
  ENDDO
  OC(1)=X
  OC(2)=Y
  OC(3)=P
  M=0
  DO T=NP(4),MIN0(NP(4)+1,FCSTGRD(4))
  DO K=NP(3),MIN0(NP(3)+1,FCSTGRD(3))
  DO J=NP(2),MIN0(NP(2)+1,FCSTGRD(2))
  DO I=NP(1),MIN0(NP(1)+1,FCSTGRD(1))
    NN(1)=I
    NN(2)=J
    NN(3)=K
    NN(4)=T
    M=M+1
    DO N=1,NUMDIMS
      AC(N,M)=NN(N)*1.0D0
    ENDDO
  ENDDO
  ENDDO
  ENDDO
  ENDDO
  CALL INTERPLTN(NUMDIMS,NGPTOBS,CO,AC,OC)
  HT=0.0
  M=0
  DO T=NP(4),MIN0(NP(4)+1,FCSTGRD(4))
  DO K=NP(3),MIN0(NP(3)+1,FCSTGRD(3))
  DO J=NP(2),MIN0(NP(2)+1,FCSTGRD(2))
  DO I=NP(1),MIN0(NP(1)+1,FCSTGRD(1))
    M=M+1
    HT=HT+CO(M)*BK0(I,J,K,T,OS)
  ENDDO
  ENDDO
  ENDDO
  ENDDO
  IF(OS.EQ.4.AND.ABS(OB-HT).GE.8.0)RETURN
  IF(OS.EQ.3.AND.ABS(OB-HT).GE.50.0)RETURN
  IF(IP.EQ.1.AND.OS.EQ.3)RETURN
  OB=OB-HT
  CALL VRTCLPSTN8(MAXGRID(3),PP0,OP(3),P,IS)
  IF(IS.NE.1)RETURN
  OC(3)=P
  O=O+1
  IF(O.GT.NOBSMAX)STOP 'NUMBER OF OBSERVATIONS EXCEEDED'
  DO N=1,NUMDIMS
    OBP(N,O)=OC(N)-1.0D0
  ENDDO
  OBS(O)=OB
  OBE(O)=OE
  STT(O)=OS
  RETURN
END SUBROUTINE HANDLEOBS

SUBROUTINE HT_TO_PRS(X0,Y0,Z,P,IS)
!*************************************************
! CONVERT HEIGHT TO PRESSURE TO USE AS MANY DATA AS POSSIBLE
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER(KIND=4) :: IS,I,J,K,T,N,M,NP(2)
  REAL(KIND=4) :: X,Y,X0,Y0
  REAL(KIND=8) :: Z,P,AC(2,4),OC(2),CO(4),HT(FCSTGRD(3))
! --------------------
  CALL LATLON_TO_RLAPSGRID(Y0,X0,Y00,X00,FCSTGRD(1),FCSTGRD(2),X,Y,IS)
  IF(X.LT.1.0.OR.Y.LT.1.0.OR.X.GT.FCSTGRD(1).OR.Y.GT.FCSTGRD(2).OR.IS.NE.1)THEN
    IS=0
    RETURN
  ENDIF
  OC(1)=X
  OC(2)=Y
  NP(1)=INT(X)
  NP(2)=INT(Y)
  DO N=1,2
    IF(NP(N).EQ.FCSTGRD(N).AND.FCSTGRD(N).NE.1)NP(N)=FCSTGRD(N)-1
  ENDDO
  M=0
  DO J=NP(2),MIN0(NP(2)+1,FCSTGRD(2))
  DO I=NP(1),MIN0(NP(1)+1,FCSTGRD(1))
    M=M+1
	AC(1,M)=I*1.0D0
    AC(2,M)=J*1.0D0
  ENDDO
  ENDDO
  CALL INTERPLTN(2,4,CO,AC,OC)
  T=1
  DO K=1,FCSTGRD(3)
    HT(K)=0.0
    M=0
    DO J=NP(2),MIN0(NP(2)+1,FCSTGRD(2))
    DO I=NP(1),MIN0(NP(1)+1,FCSTGRD(1))
      M=M+1
      HT(K)=HT(K)+CO(M)*BK0(I,J,K,T,3)
    ENDDO
    ENDDO
  ENDDO
  CALL ZPCONVERT(FCSTGRD(3),HT,P00,Z,P,IS)
  RETURN
END SUBROUTINE HT_TO_PRS

SUBROUTINE ADDBKGRND
!*************************************************
! SET BACKGROUND AS OBSERVATIONS AT THE GRID POINTS WHERE OBSERVATIONS CAN NOT AFFECT
! HISTORY: OCTOBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER(KIND=4) :: I,J,K,T,S,O,N,RX,RY,NP(MAXDIMS),NS(NUMSTAT)
  REAL(KIND=8)    :: RH,RZ,OE(NUMSTAT),OC(NUMDIMS),P,PP(MAXGRID(3))
! --------------------
  IF(NALLOBS.EQ.0)RETURN
  DO S=1,NUMSTAT
    DO T=1,FCSTGRD(4)
    DO K=1,MAXGRID(3)      !!!!!!!ATTENTION, DUE TO SPECIAL CHARACTER OF VERTICAL COORDINATE
    DO J=1,FCSTGRD(2)
    DO I=1,FCSTGRD(1)
      GRIDMASK(I,J,K,T,S)=1
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDDO
  DO S=1,NUMSTAT
    OE(S)=0.0D0
  ENDDO
  DO S=1,NUMSTAT
    NS(S)=0
  ENDDO
  DO O=1,NALLOBS
    S=STT(O)
    NS(S)=NS(S)+1
    DO N=1,NUMDIMS
      OC(N)=OBP(N,O)+1.0D0
    ENDDO
    IF(IFPCDNT.EQ.1)THEN
      DO K=1,MAXGRID(3)
        PP(K)=PP0(K)
      ENDDO
    ENDIF
    IF(MAXGRID(3).GE.2)THEN
      K=IDINT(OC(3))
      P=(OC(3)-K)*(PP(K+1)-PP(K))+PP(K)
    ENDIF
    RX=0
    RY=0
    IF(FCSTGRD(1).GE.2)RX=IDINT(OBRADIUS(1,S)/DX0)+1
    IF(FCSTGRD(2).GE.2)RY=IDINT(OBRADIUS(2,S)/DY0)+1
    DO T=1,FCSTGRD(4)
    DO K=1,MAXGRID(3)      !!!!!!!ATTENTION, DUE TO SPECIAL CHARACTER OF VERTICAL COORDINATE
    DO J=MAX0(IDINT(OC(2))-RY,1),MIN0(IDINT(OC(2))+RY+1,FCSTGRD(2))
    DO I=MAX0(IDINT(OC(1))-RX,1),MIN0(IDINT(OC(1))+RX+1,FCSTGRD(1))
      RH=DSQRT(((OC(1)-I)*DX0)**2+((OC(2)-J)*DY0)**2)
      RZ=DABS((PP(K)-P))
      IF(RH.LE.OBRADIUS(1,S).AND.RZ.LE.OBRADIUS(3,S))GRIDMASK(I,J,K,T,S)=0
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    OE(S)=DMAX1(OE(S),OBE(O)*1.0D0)
  ENDDO
  O=NALLOBS
  DO S=1,NUMSTAT
    DO T=1,FCSTGRD(4)
    DO K=1,MAXGRID(3)      !!!!!!!ATTENTION, DUE TO SPECIAL CHARACTER OF VERTICAL COORDINATE
    DO J=1,FCSTGRD(2),3
    DO I=1,FCSTGRD(1),3
      IF(GRIDMASK(I,J,K,T,S).EQ.1.AND.NS(S).GE.1)THEN
        NP(1)=I-1
        NP(2)=J-1
        NP(3)=K-1
        NP(4)=T-1
        O=O+1
        IF(O.GT.NOBSMAX)STOP 'NUMBER OF OBSERVATIONS EXCEEDED'
        DO N=1,NUMDIMS
          OBP(N,O)=NP(N)
        ENDDO
        OBS(O)=0.0D0
        OBE(O)=OE(S)
        STT(O)=S
      ENDIF
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDDO
  NALLOBS=O
  RETURN
END SUBROUTINE ADDBKGRND

SUBROUTINE DRCONTOUR(OX,EX,OY,EY,ND,NG,UV,FN)
!*************************************************
! DRAW RADIAL WIND (AFFILIATE)
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER(KIND=4) :: ND,NG(NUMDIMS),IU,I,J,K
  REAL(KIND=8) :: MN,MX,OX,EX,OY,EY
  REAL(KIND=8) :: UV(NG(1),NG(2),NG(3))
  CHARACTER(LEN=6) :: FN
! --------------------
  K=2 !(NG(3)+1)/2
  IU=2
  OPEN(IU,FILE=FN,STATUS='UNKNOWN')
  MX=-100000000.0
  MN=100000000.0
  DO I=1,NG(1)
  DO J=1,NG(2)
    IF(UV(I,J,K).LT.9000.0)THEN
      IF(UV(I,J,K).GT.MX)MX=UV(I,J,K)
      IF(UV(I,J,K).LT.MN)MN=UV(I,J,K)
    ENDIF
  ENDDO
  ENDDO
  WRITE(IU,'(A4)')'DSAA'
  WRITE(IU,'(2I4)')NG(1),NG(2)
  WRITE(IU,*)OX,EX
  WRITE(IU,*)OY,EY
  WRITE(IU,*)MN,MX
  DO I=1,NG(1)
  DO J=1,NG(2)
    IF(UV(I,J,K).GT.9000.0)UV(I,J,K)=2.E38
  ENDDO
  ENDDO
  DO J=1,NG(2)
    WRITE(IU,*)(UV(I,J,K),I=1,NG(1))
  ENDDO
  CLOSE(IU)
  RETURN
END SUBROUTINE DRCONTOUR

SUBROUTINE DRCONTOUR_2D(OX,EX,OY,EY,IM,JM,UV,FN)
!*************************************************
! DRAW RADIAL WIND (AFFILIATE)
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER(KIND=4) :: IU,I,J,K,IM,JM
  REAL(KIND=8) :: MN,MX,OX,EX,OY,EY
  REAL(KIND=8) :: UV(IM,JM)
  CHARACTER(LEN=6) :: FN
! --------------------
  IU=2
  OPEN(IU,FILE=FN,STATUS='UNKNOWN')
  MX=-100000000.0
  MN=100000000.0
  DO I=1,IM
  DO J=1,JM
    IF(UV(I,J).LT.9000.0)THEN
      IF(UV(I,J).GT.MX)MX=UV(I,J)
      IF(UV(I,J).LT.MN)MN=UV(I,J)
    ENDIF
  ENDDO
  ENDDO
  WRITE(IU,'(A4)')'DSAA'
  WRITE(IU,'(2I4)')IM,JM
  WRITE(IU,*)OX,EX
  WRITE(IU,*)OY,EY
  WRITE(IU,*)MN,MX
  DO I=1,IM
  DO J=1,JM
    IF(UV(I,J).GT.9000.0)UV(I,J)=2.E38
  ENDDO
  ENDDO
  DO J=1,JM
    WRITE(IU,*)(UV(I,J),I=1,IM)
  ENDDO
  CLOSE(IU)
  RETURN
END SUBROUTINE DRCONTOUR_2D

END MODULE INPUT_BG_OBS
