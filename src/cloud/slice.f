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
cdis
cdis   
cdis

        SUBROUTINE SLICE (A,IMAX,JMAX,KMAX,B,ib,jb,HZ,EW,NS,KK,JJ,II)
C THIS SUBROUTNE ALLOWS A SLICE TO BE TAKEN THROUGH A 3-D ARRAY
C CONTROL IS BY SETTING INTEGERS HZ,EW,NS, TO 1 OR 0.  KK,JJ,II
C THEN DETERMINE WHICH SLICE IS TAKEN. RESULT IS PUT INTO B.
        INTEGER HZ,EW,NS
        DIMENSION A(IMAX,JMAX,KMAX),B(ib,jb)

        IF (HZ.EQ.1) THEN
             DO J=1,JMAX
             DO I=1,IMAX
             B(I,J)=A(I,J,KK)
             ENDDO
             ENDDO
        endif

        IF(EW.EQ.1) THEN
                DO K=1,KMAX
                DO I=1,IMAX
                B(I,K)=A(I,JJ,K)
                ENDDO
                ENDDO
        endif

        IF(NS.EQ.1) THEN
                DO K=1,KMAX
                DO J=1,JMAX
                B(J,K)=A(II,J,K)
                ENDDO
                ENDDO
        endif

        RETURN
        END
