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
