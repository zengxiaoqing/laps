
!====================================================
!  This header file declares necessary variables for
!  the optimization package: LBFGSB.
!
!  Author: YUANFU XIE   NOAA/FSL
!  Date:   AUGUST 2000
!
!====================================================
      
!.....Variables required by LBFGS_B:

INTEGER, PARAMETER :: msave=7
INTEGER, PARAMETER :: mvar = mx*my*mt
CHARACTER*60       :: ctask,csave
DOUBLE PRECISION   :: wk(mvar*(2*msave+4)+ &
                         12*msave*msave+12*msave)
DOUBLE PRECISION   :: bdlow(mvar),bdupp(mvar)
DOUBLE PRECISION   :: factr,dsave(29)
INTEGER            :: iprnt,nbund(mvar),iwrka(3*mvar),isbmn
LOGICAL            :: isave(44),lsave(4)

! Put them in a common block:
COMMON /lbfgs/ctask,csave,wk,bdlow,bdupp, &
              factr,dsave,iprnt,nbund,iwrka, &
              isave,lsave,isbmn

!.....End of LBFGS_B declarations.
