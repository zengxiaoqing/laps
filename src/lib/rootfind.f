cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS
cdis
cdis    This software and its documentation are in the public domain and
cdis    are furnished "as is."  The United States government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  They assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  If significant modifications or enhancements
cdis    are made to this software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis
      FUNCTION ROOTFIND(FUNC,X1,X2,TOL)

C     THIS ROUTINE RETURNS ROOTS BASED ON A METHOD
C     BY VAN WIJNGAARDEN-DEKKER-BRENT AND IS A CUT
C     ABOVE MANY CONTEMPORARY METHODS.

C     INPUTS ARE:

C     FUNC -- A DEFINED FUNCTION WITH A ROOT KNOWN
C     TO LIE BETWEEN X1 AND X2.

C     X1 -- BRACKET VALUE

C     X2 -- BRACKET VALUE

C     TOL -- THE TOLERANCE IS THE DESIRED ACCURACY

C     CODED UP SEPTEMBER 29 1993 -- DAN BIRKENHEUER
C     CODED TO COMPLY WITH F77 STANDARDS.  THIS ROUTINE
C     SHOULD BE FULLY PORTABLE.

      PARAMETER (ITMAX=100,EPS=3.E-8)
      EXTERNAL FUNC
      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)
      IF(FB*FA.GT.0.) WRITE (6,*) 'Root must be bracketed for ROOTFIND.'
      FC=FB
      DO 2 ITER=1,ITMAX
        IF(FB*FC.GT.0.) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.*EPS*ABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          ROOTFIND=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.*XM*S
            Q=1.-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B)
2     CONTINUE
      WRITE (6,*) 'module ROOTFIND exceeding maximum iterations.'
      ROOTFIND=B
      RETURN
      END
