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
        integer*4 j_status(20),iprod_number(20),i4time_array(20)
        character*9 a9_time

        call get_systime(i4time,a9_time,istatus)
        if(istatus .ne. 1)go to 999

!      (-1) DUMMY PROCESS
!       (0) Normal full Cloud Analysis
!       (1) Calculate only main fields,
!           derived fields were moved elsewhere
!       (2) Reread data, then calc derived fields
!           (for testing)
!       (3) means derived prods only

        isplit = 1

        call laps_cloud(i4time,
     1                  i_diag,
     1                  n_prods,
     1                  iprod_number,
     1                  isplit,
     1                  j_status)

999     continue

        end

