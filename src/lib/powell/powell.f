      SUBROUTINE POWELL(P,XI,N,NP,FTOL,ITER,FRET,func)
      PARAMETER (NMAX=40,ITMAX=50)
      EXTERNAL FUNC
      real func ! funciton type
      DIMENSION P(NP),XI(NP,NP),PT(NMAX),PTT(NMAX),XIT(NMAX)
      FRET=FUNC(P)
      DO 11 J=1,N
        PT(J)=P(J)
11    CONTINUE
      ITER=0
1     ITER=ITER+1
      FP=FRET
      IBIG=0
      DEL=0.
      DO 13 I=1,N
        DO 12 J=1,N
          XIT(J)=XI(J,I)
12      CONTINUE
        FPTT=FRET
        CALL LINMIN(P,XIT,N,FRET)
        IF(ABS(FPTT-FRET).GT.DEL)THEN
          DEL=ABS(FPTT-FRET)
          IBIG=I
        ENDIF
13    CONTINUE
      IF(2.*ABS(FP-FRET).LE.FTOL*(ABS(FP)+ABS(FRET)))RETURN
c      IF(ITER.EQ.ITMAX) PAUSE 'Powell exceeding maximum iterations.'
      IF(ITER.EQ.ITMAX) then
c		type*, 'Powell exceeding maximum iterations.'
		return
	endif
      DO 14 J=1,N
        PTT(J)=2.*P(J)-PT(J)
        XIT(J)=P(J)-PT(J)
        PT(J)=P(J)
14    CONTINUE
      FPTT=FUNC(PTT)
      IF(FPTT.GE.FP)GO TO 1
      T=2.*(FP-2.*FRET+FPTT)*(FP-FRET-DEL)**2-DEL*(FP-FPTT)**2
      IF(T.GE.0.)GO TO 1
      CALL LINMIN(P,XIT,N,FRET)
      DO 15 J=1,N
        XI(J,IBIG)=XIT(J)
15    CONTINUE
      GO TO 1
      END
