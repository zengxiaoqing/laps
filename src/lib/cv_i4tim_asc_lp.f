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
        subroutine      cv_i4tim_asc_lp(i4time,atime,istatus)
C
C     Takes in an i4time and returns the time as an ASCII string
C     (e.g. 27-MAR-1990 12:30:00.00 ).  The i4time is assumed to
C     be a 1960-relative time, although the starting year is easily
C     changed in the code.
C
C     IMPORTS - i4time
C
C     EXPORTS - atime, istatus
C
C================================================================
C

        implicit        none

        integer*4       i4time,
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
