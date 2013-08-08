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
cdoc    THE TIME CORRESPONDING TO I4TIME (seconds since 1-1-1960).
C
C================================================================
C
        INTEGER I4TIME, ISTATUS
        CHARACTER*9 FILE_NAME
C
        INTEGER NYEAR, NMONTH, NDAY, NHOUR, NMIN, NSEC, NJULIAN, I
C
        INTEGER NJUL_DAYS(12)
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
cdoc    CV_I4TIM_INT_LP CONVERTS I4TIME (seconds since 1-1-1960) TO SIX INTEGERS
cdoc    Note that NYEAR is the number of years since 1900
C
C================================================================
C
C       INTEGER I4TIME, NYEAR, NMONTH, NDAY, NHOUR, NMIN, NSEC, ISTATUS
C
        INTEGER NSECMO(12), NSECYR, NSECDA, NSECHR, NSECMN
        INTEGER IBASE, I, NN, LFTOVR, NDAYS

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

      subroutine c_time2fname(utime,a9time)

cdoc  Convert utime to a9time. Jacket routine that calls 'make_fnam_lp'.

      integer utime, i4time, istatus
      character*(*) a9time

      i4time = utime + 315619200
      call make_fnam_lp (i4time,a9time,istatus)

      return
      end


      subroutine afwa_julhr_i4time(I_A1JUL,I_A1MIN,i4time)

cdoc  I_A1JUL is number of hours since Dec 31, 1967 at 00z
cdoc  This is converted to i4time, number of sec since Jan 1, 1960 at 00z
      i4time_hr  = I_A1JUL * 3600 + (8*365 - 1 + 2) * 86400
      i4time_min = I_A1MIN*60
      i4time     = i4time_hr + i4time_min

      return
      end


        subroutine cv_asc_i4time(ascii_time,I4time)

cdoc    Converts 9 character ascii time into i4time (seconds since 1-1-1960)

        character*9 ascii_time
        integer i4time ! seconds since 1-1-1960

        read(ascii_time(1:2),2,err=900)iyear
        read(ascii_time(3:5),3,err=900)idoy
        read(ascii_time(6:7),2,err=900)ihour
        read(ascii_time(8:9),2,err=900)imin

2       format(i2)
3       format(i3)

!       Valid for years 1960-2060
        if(iyear .lt. 60)iyear = iyear + 100

        lp = (iyear + 3 - 60) / 4

        i4time =  (iyear-60) * 31536000
     1  + (idoy-1)   * 86400
     1  + ihour      * 3600
     1  + imin       * 60

        i4time = i4time + 86400 * lp

        return

!       Error return
900     write(6,*)' Error in cv_asc_i4time: ascii_time = ',ascii_time

        return

        end

C
        subroutine      cv_i4tim_asc_lp(i4time,atime,istatus)
C
cdoc  Takes in an i4time and returns the time as an ASCII string
cdoc  (e.g. 27-MAR-1990 12:30:00.00 ).  The i4time is assumed to
cdoc  be a 1960-relative time, although the starting year is easily
cdoc  changed in the code.
C
C     IMPORTS - i4time ! seconds since 1-1-1960
C
C     EXPORTS - atime, istatus
C
C================================================================
C

        implicit        none

        integer       i4time,
     1          istatus,
     1          rmndr,
     1          nsec,
     1          monthsec(12),
     1          year,
     1          month,
     1          day,
     1          hour,
     1          min,
     1          sec

        character*24    atime
        character*4     ayear
        character*3     amonth(12)
        character*2     aday
        character*2     ahour
        character*2     amin
        character*2     asec

        data            monthsec/2678400,2419200,2678400,2592000,
     1                   2678400,2592000,2678400,2678400,
     1                   2592000,2678400,2592000,2678400/

        data            amonth/'JAN','FEB','MAR','APR','MAY','JUN',
     1                 'JUL','AUG','SEP','OCT','NOV','DEC'/

C
C================================================================
C

        if (i4time .lt. 0) then
           istatus=0
           write (6,*) 'Error in input to cv_i4tim_asc_lp: negative time
     1'
           return
        endif

        rmndr=i4time
        do year=1960,2100
                if (mod(year,4) .eq. 0) then
                        nsec=31622400
                else
                        nsec=31536000
                endif
                if (rmndr .lt. nsec) goto 10
                rmndr=rmndr-nsec
        enddo

10      do month=1,12
                nsec=monthsec(month)
                if (mod(year,4) .eq. 0 .and. month .eq. 2) nsec=nsec+864
     100
                if (rmndr .lt. nsec) goto 20
                rmndr=rmndr-nsec
        enddo

20      do day=1,31
                if (rmndr .lt. 86400) goto 30
                rmndr=rmndr-86400
        enddo

30      do hour=0,23
                if (rmndr .lt. 3600) goto 40
                rmndr=rmndr-3600
        enddo

40      do min=0,59
                if (rmndr .lt. 60) goto 50
                rmndr=rmndr-60
        enddo

50      sec=rmndr

!       encode(4,900,ayear) year
!       encode(2,901,aday) day
!       encode(2,901,ahour) hour
!       encode(2,901,amin) min
!       encode(2,901,asec) sec

        write(ayear,900) year
        write(aday,901) day
        write(ahour,901) hour
        write(amin,901) min
        write(asec,901) sec
900     format(i4)
901     format(i2)

!       if (day .lt. 10) aday(1:1)='0'
        if (hour .lt. 10) ahour(1:1)='0'
        if (min .lt. 10) amin(1:1)='0'
        if (sec .lt. 10) asec(1:1)='0'

        atime=aday//'-'//amonth(month)//'-'//ayear//' '//ahour//':'//
     1      amin//':'//asec//'.00 '

        istatus=1
        return
        end

      subroutine jd_to_i4time(jd,i4time,istatus)

      double precision jd, jd_1960
      integer i4time ! seconds since 1-1-1960

cdoc  Converts Julian Day (Number of days since Jan 1 4713BCE to i4time)

!     Author: Steve Albers 2005

      jd_1960 = 2436934.5

      days_since_1960 = jd - jd_1960

      i4time = int(days_since_1960 * 86400.)

      istatus = 1

      return
      end

      subroutine i4time_to_jd(i4time,jd,istatus)

      double precision jd, jd_1960, days_since_1960
      integer i4time,istatus ! seconds since 1-1-1960

cdoc  Converts i4time to Julian Day (Number of days since Jan 1 4713BCE)

!     Author: Steve Albers 2013

      jd_1960 = 2436934.5

      days_since_1960 = dble(i4time) / 86400.D0

!     write(6,*)' days since 1960 is ',days_since_1960

      jd = jd_1960 + days_since_1960              

      istatus = 1

      return
      end
