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

      SUBROUTINE LINMIN(P,XI,N,FRET)
      PARAMETER (NMAX=50,TOL=1.E-4)
      EXTERNAL F1DIM
      DIMENSION P(N),XI(N)
      COMMON /F1COM/ NCOM,PCOM(NMAX),XICOM(NMAX)
      NCOM=N
      DO 11 J=1,N
        PCOM(J)=P(J)
        XICOM(J)=XI(J)
11    CONTINUE
      AX=0.
      XX=1.
      CALL MNBRAK(AX,XX,BX,FA,FX,FB,F1DIM)
      FRET=BRENT(AX,XX,BX,F1DIM,TOL,XMIN)
      DO 12 J=1,N
        XI(J)=XMIN*XI(J)
        P(J)=P(J)+XI(J)
12    CONTINUE
      RETURN
      END
