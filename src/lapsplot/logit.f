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

c
c..............................  logit  ...............................
c
c  logit makes a log of user selections by appending to getprd_log.dat.
c
        subroutine logit(str)
        implicit none
        character*(*) str
        integer lenin
        character*80 reqin
c
        open(23,status='old',err=5004,
     .    file='User_Dev:lapsplot.log')

!5001   read(23,5002,end=5003,err=5004) lenin,reqin
!5002   format(q,a)
!       go to 5001

5003    write(23,*) str
        close(23)
5004    continue
        return
        end
c
