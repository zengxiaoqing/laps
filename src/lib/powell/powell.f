cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis 

      SUBROUTINE POWELL(P,XI,N,NP,FTOL,ITER,FRET,func,print_switch)
      PARAMETER (NMAX=40,ITMAX=5)
      EXTERNAL FUNC
      integer print_switch
      real func                 ! funciton type
      DIMENSION P(NP),XI(NP,NP),PT(NMAX),PTT(NMAX),XIT(NMAX)
      FRET=FUNC(P)
      if(fret.eq.0.0) then      !notify
         if (print_switch .eq. 1) then
            write(6,*)'POWELL:fret = 0.0'
         endif
      endif
      DO 11 J=1,N
         PT(J)=P(J)
 11   CONTINUE
      ITER=0
 1    ITER=ITER+1
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
 13   CONTINUE
      IF(2.*ABS(FP-FRET).LE.FTOL*(ABS(FP)+ABS(FRET)))then
c         write(6,*) 'POWELL difference less than FTOL'
c         write(6,*) fp, fret, ftol,'fp, fret,ftol'
         RETURN
      endif
c     IF(ITER.EQ.ITMAX) PAUSE 'Powell exceeding maximum iterations.'
      IF(ITER.EQ.ITMAX) then
         if (print_switch .eq. 1) then
            write(6,*) 'Powell exceeding maximum iterations.'
         endif
         return
      endif
      DO 14 J=1,N
         PTT(J)=2.*P(J)-PT(J)
         XIT(J)=P(J)-PT(J)
         PT(J)=P(J)
 14   CONTINUE
      FPTT=FUNC(PTT)
      IF(FPTT.GE.FP)GO TO 1
      T=2.*(FP-2.*FRET+FPTT)*(FP-FRET-DEL)**2-DEL*(FP-FPTT)**2
      IF(T.GE.0.)GO TO 1
      CALL LINMIN(P,XIT,N,FRET)
      DO 15 J=1,N
         XI(J,IBIG)=XIT(J)
 15   CONTINUE
      GO TO 1
      return
      END
      
