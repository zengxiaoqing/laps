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
        SUBROUTINE i4time_fname_lp (FNAME_IN, I4TIME, ISTATUS)
C
C       This routine converts a standard PROFS file name, i.e. yydddhhmm,
C       into the corresponding I4 time.
C
C       ON INPUT
C          FILE_NAME - The PROFS style file name.
C
C       ON OUTPUT
C          I4TIME - The corresponding I4 time of the file name.
C          ISTATUS -  The return status.
C
C================================================================
C
        CHARACTER*9 wfo_fname13_to_fname9
        CHARACTER*9 FILE_NAME
        CHARACTER*(*) FNAME_IN
        INTEGER*4 I4TIME, ISTATUS
C
        INTEGER*4 INT_FILE(9), I, NYEAR, JDAY, NHOUR, MIN, MONTH, NDAY
        INTEGER*4 I4TIME_INT_LP
C
C================================================================
C
C       Check the length of FILE_NAME to be sure that it is valid.
C
        call get_filetime_length(fname_in,lenf)

        if(lenf.ne.9.and.lenf.ne.13)then
           write(6,*)' error in i4time_fname_lp: ',fname_in,lenf
           go to 1000
        endif

c
c added 9-3-96 JSmart for use in conversions of wfo filenames of
c the type yyyymmdd_hhmm (13 characters).
c

        if(fname_in(9:9) .eq. '_')then ! We have a WFO file
           file_name = wfo_fname13_to_fname9(fname_in)

        else   !assume 9 character filename type
           file_name = fname_in

        endif
C
C       Split the file name into individual integers, checking for invalid
C       characters in the file name.
C
        DO I = 1,9
           INT_FILE(I) = ICHAR(FILE_NAME(I:I))-48
           IF (INT_FILE(I) .LT. 0 .OR. INT_FILE(I) .GT.9) THEN
              GO TO 1500
           END IF
        END DO
C
C       Convert these numbers in year, day of year (julian day), hours, 
C       and minutes.
        NYEAR = 10*INT_FILE(1) + INT_FILE(2)
        JDAY = 100*INT_FILE(3) + 10*INT_FILE(4) + INT_FILE(5)
        NHOUR = 10*INT_FILE(6) + INT_FILE(7)
        MIN = 10*INT_FILE(8) + INT_FILE(9)

C       Get the actual year, here is where we assume what century it is...
C       Steve Albers 1997
        if(NYEAR .lt. 60)then
            NYEAR = NYEAR + 2000
        else
            NYEAR = NYEAR + 1900
        endif
C
C       Convert the day of year (julian day) into month and day.
C
        CALL CV_JUL_MMDD_LP (JDAY, NYEAR, MONTH, NDAY, ISTATUS)
        IF (0 .EQ. ISTATUS)
     1    GO TO 2000
C
C       Convert these integers into I4 time.
C
        I4TIME = I4TIME_INT_LP (NYEAR, MONTH, NDAY, NHOUR, MIN, 0, ISTAT
     1US)
        IF (0 .EQ. ISTATUS)
     1    GO TO 2000
C
C       Successful return
C
        ISTATUS = 1
        RETURN
C
C       Error detected in this routine.
C
 1500   CONTINUE
        ISTATUS = 0
        WRITE( 6,* ) 'Error in i4time_fname_lp: bad digit in file name.'
        WRITE( 6,* ) FILE_NAME
        RETURN
C
C       Error detected in called routine
C
 2000   CONTINUE
        ISTATUS = 0
        WRITE( 6,* ) 'Error in i4time_fname_lp: error in subroutine.'
        WRITE( 6,* ) FILE_NAME
        RETURN

 1000   CONTINUE
        ISTATUS = 0
        WRITE( 6,* ) 'Error in i4time_fname_lp: wrong length for file na
     1me.'
        WRITE( 6,* ) FILE_NAME

        RETURN
C
        END
C
        SUBROUTINE CV_JUL_MMDD_LP (JULIAN_DAY, YEAR, MONTH, DAY, ISTATUS
     1)
C
C       This routine converts from day of year (Julian days) to month and 
C       day in integer format.
C
C       ON INPUT
C          JULIAN_DAY - The Day of Year (Julian date) to be converted.
C          YEAR - The year of the Day of Year (Julian date), if .lt. 100, 
C                 it is assumed to be the last two digits of 19xx.
C
C       ON OUTPUT
C          MONTH - The integer representation for the month (1-12)
C          DAY - The integer value of the day of the month.
C          ISTATUS - The return status.
C
C================================================================
C
        INTEGER*4 JULIAN_DAY, YEAR, DAY, ISTATUS
C
        INTEGER*4 MNTH(12), TEMP_YEAR, MONTH, MAX_DAY
        LOGICAL*1 LEAP
        DATA MNTH/31,0,31,30,31,30,31,31,30,31,30,31/
C
C================================================================
C
        IF (YEAR.LT.100) THEN
           TEMP_YEAR = YEAR + 1900
        ELSE
           TEMP_YEAR = YEAR
        ENDIF
C
C       Check to see if the year is a leap year
C
        LEAP = ( MOD(TEMP_YEAR,4)  .EQ. 0 )
C
        IF ( LEAP ) THEN
           MAX_DAY = 366
           MNTH(2) = 29
        ELSE
           MAX_DAY = 365
           MNTH(2) = 28
        END IF
C
C       Check day of year (JULIAN_DAY) for being too large or too small
C
        IF (JULIAN_DAY .LE. 0 .OR. JULIAN_DAY .GT. MAX_DAY) THEN
           ISTATUS = 0
           WRITE( 6,* ) 'Error in cv_jul_mmdd_LP.f: bad Julian day', JUL
     1IAN_DAY
           RETURN
        END IF
C
C       Convert to month and day
C
        MONTH = 1
        DAY = JULIAN_DAY
        DO WHILE (DAY .GT. MNTH(MONTH))
           DAY = DAY - MNTH(MONTH)
           MONTH = MONTH + 1
        END DO
        ISTATUS = 1
C
        RETURN
        END

C
        FUNCTION I4TIME_INT_LP (NYEAR,NMONTH,NDAY,NHOUR,NMIN,NSEC,ISTATU
     1S)
C
C       I4TIME_INT_LP RETURNS I4 TIME (# OF SECONDS SINCE 00:00 01-JAN-60)
C         GIVEN A 6 INTEGERS CONTAINING YEAR, MONTH, DAY, HOUR,
C         MINUTE, SECOND
C
C================================================================
C
        INTEGER*4 I4TIME_INT_LP
        INTEGER*4 NYEAR, NMONTH, NDAY, NHOUR, NMIN, NSEC, ISTATUS
C
        INTEGER*4 NSECYR, NSECDA, NSECHR, NSECMN
        INTEGER*4 NYR, NYRS, NLEAP
        INTEGER*4 NDAYS(12), NSECMO(12)
        INTEGER*4 IBASE, ISUM, I
        PARAMETER (IBASE=60)

        DATA NSECYR/31536000/NSECDA/86400/NSECHR/3600/NSECMN/60/
        DATA NSECMO/2678400,2419200,2678400,2592000,2678400,
     1  2592000,2678400,2678400,2592000,2678400,2592000,2678400/

        DATA NDAYS/31,29,31,30,31,30,31,31,30,31,30,31/
C
C================================================================
C
        ISTATUS = 1
C
C Sum the number of years.
C
        NYR=NYEAR
        IF (NYR.GT.1900) NYR=NYR-1900
        NYRS=NYR-IBASE
        IF (NYRS.LT.0.OR.NYRS.GT.67) GO TO 1000
        ISUM=NYRS*NSECYR
        NLEAP=(NYRS/4)+1                        ! Account for leap years.
        ISUM=ISUM+(NLEAP*NSECDA)
C
C Sum in the number of months.
C
        IF (NMONTH.LT.1.OR.NMONTH.GT.12) GO TO 1000
C
        IF (NMONTH.NE.1) THEN
                DO I=1,NMONTH-1
                        ISUM=ISUM + NSECMO(I)
                END DO
        END IF
C
C Correct for Jan or Feb of a leap yer.
C
!       IF (MOD(NYRS,4).EQ.0) THEN
        IF (MOD(NYR,4).EQ.0) THEN       ! Steve Albers, Linda Wharton 1993
                IF (NMONTH.LE.2) ISUM=ISUM-NSECDA
        END IF
C
C Sum in the number of days.
C
        IF (NDAY.LT.1.OR.NDAY.GT.NDAYS(NMONTH))
     1  GO TO 1000
        ISUM=ISUM+((NDAY-1)*NSECDA)
C
C Sum in the number of hours.
C
        IF (NHOUR.LT.0.OR.NHOUR.GT.23) GO TO 1000
        ISUM=ISUM+(NHOUR*NSECHR)
C
C Sum in the number of minutes.
C
        IF (NMIN.LT.0.OR.NMIN.GT.59) GO TO 1000
        ISUM=ISUM+(NMIN*NSECMN)
C
C Sum in the number of seconds.
C
        IF (NSEC.LT.0.OR.NSEC.GT.59) GO TO 1000
        ISUM=ISUM+NSEC
C
        I4TIME_INT_LP=ISUM
        RETURN
C
C Return error code.
C
1000    ISTATUS = 0
        I4TIME_INT_LP=0
        WRITE(6,*) 'Error in I4TIME_INT_LP'
        RETURN
C
        END

