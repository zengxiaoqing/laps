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
C
        SUBROUTINE make_fnam_lp (I4TIME, FILE_NAME, ISTATUS)
C
cdoc    make_fnam_lp CONSTRUCTS THE FILE NAME STRING 'yyjjjhhmm' FOR
cdoc    THE TIME CORRESPONDING TO I4TIME.
C
C================================================================
C
        INTEGER*4 I4TIME, ISTATUS
        CHARACTER*9 FILE_NAME
C
        INTEGER*4 NYEAR, NMONTH, NDAY, NHOUR, NMIN, NSEC, NJULIAN, I
C
        INTEGER*4 NJUL_DAYS(12)
        DATA NJUL_DAYS/0,31,59,90,120,151,181,212,243,273,304,334/
C
C================================================================
C
C Convert I4 time to string yyjjjhhmm.
C
        CALL CV_I4TIM_INT_LP (I4TIME, NYEAR, NMONTH, NDAY, NHOUR, NMIN,
     1                      NSEC, ISTATUS)
        IF (istatus .ne. 1) GO TO 100

        NJULIAN = NJUL_DAYS(NMONTH) + NDAY

        IF (NMONTH .GT. 2 .AND. MOD (NYEAR,4) .EQ. 0) NJULIAN=NJULIAN +
     11

        NYEAR = mod(NYEAR,100) ! Steve Albers 1997

        WRITE(FILE_NAME,1001,ERR=90) NYEAR, NJULIAN, NHOUR, NMIN
1001    FORMAT (I2.2,I3.3,I2.2,I2.2)

        DO I = 1, 9
           IF (FILE_NAME(I:I) .EQ. ' ') FILE_NAME(I:I) = '0'
        END DO

        ISTATUS = 1
        RETURN
C
C Error in ENCODE.
C
90      ISTATUS = 0
        WRITE( 6,* ) 'Error in make_fnam_lp: error in encode.'
        RETURN
C
C Error in subroutine...
C
100     ISTATUS = 0
        WRITE( 6,* ) 'Error in make_fnam_lp: error in subroutine.'
        RETURN
        END
C
        SUBROUTINE CV_I4TIM_INT_LP (I4TIME,NYEAR,NMONTH,NDAY,NHOUR,
     1                     NMIN,NSEC,ISTATUS)
C
cdoc    CV_I4TIM_INT_LP CONVERTS I4 TIME TO SIX INTEGERS
C
C================================================================
C
C       INTEGER*4 I4TIME, NYEAR, NMONTH, NDAY, NHOUR, NMIN, NSEC, ISTATUS
C
        INTEGER*4 NSECMO(12), NSECYR, NSECDA, NSECHR, NSECMN
        INTEGER*4 IBASE, I, NN, LFTOVR, NDAYS

        PARAMETER (IBASE=60)

        DATA NSECYR/31536000/NSECDA/86400/NSECHR/3600/NSECMN/60/
        DATA NSECMO/2678400,0,2678400,2592000,2678400,
     1  2592000,2678400,2678400,2592000,2678400,2592000,2678400/
C
C================================================================
C
        ISTATUS = 1
C
C Verify input.
C
        IF (I4TIME .LT. 0) THEN
           ISTATUS = 0
           WRITE(6,*) 'Error in CV_I4TIM_INT_LP: negative time ', I4TIME
           RETURN
        END IF

C
C Subtract out number of years.
C
        LFTOVR = I4TIME
        NN = LFTOVR
        DO I = 0, 70
           IF (MOD(I,4) .EQ. 0) THEN
              NN = NN - (NSECYR + NSECDA)
           ELSE
              NN = NN - NSECYR
           END IF
           IF (NN .LT. 0) GO TO 8
           LFTOVR = NN
        END DO
8       NYEAR = IBASE + I

C
C Subtract out number of months.
C
        NSECMO(2) = 2419200
        IF (MOD(NYEAR,4) .NE. 0) GO TO 10
           NSECMO(2) = 2505600
10      NN=LFTOVR
        DO I=1,12
           NN = NN - NSECMO(I)
           IF (NN .LT. 0) GO TO 30
           LFTOVR = NN
        END DO
30      NMONTH = I

C
C Subtract out number of days.
C
        NDAYS = LFTOVR / NSECDA
        LFTOVR = LFTOVR - (NDAYS * NSECDA)
        NDAY = NDAYS + 1

C
C Subtract out number of hours.
C
        NHOUR = LFTOVR / NSECHR
        LFTOVR = LFTOVR - (NHOUR * NSECHR)

C
C Subtract out number of minutes.
C
        NMIN = LFTOVR / NSECMN
        LFTOVR = LFTOVR - (NMIN * NSECMN)

C
C What's left over is number of seconds.
C
        NSEC = LFTOVR
        RETURN
        END

      subroutine make_fnam13_lp(initial_i4time,forecast_time,filename,
     +     status)

cdoc  Converts initial time and forecast time to a 13 character filename

      integer initial_i4time, forecast_time, status
      character*13 filename

      call make_fnam_lp(initial_i4time,filename,status)
      write(filename(10:13),'(i2.2,i2.2)') forecast_time/3600,
     +     mod(forecast_time,60)
      return
      end

