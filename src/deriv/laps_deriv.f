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

!                   S. Albers - Original Version

!       1997 Jul 31 K. Dritz  - Added call to get_grid_dim_xy to get the
!                               values of NX_L, NY_L.
!       1997 Jul 31 K. Dritz  - Now pass NX_L, NY_L as arguments to laps_cloud.
!       1997 Jul 31 K. Dritz  - Added call to get_meso_sao_pirep to get the
!                               value of N_PIREP, which is passed to
!                               laps_cloud.
!       1997 Jul 31 K. Dritz  - Added call to get_maxstns, and pass the value
!                               of maxstns to laps_cloud.
!       1997 Jul 31 K. Dritz  - Compute max_cld_snd as maxstns + N_PIREP and
!                               pass to laps_cloud.

        integer j_status(20),iprod_number(20),i4time_array(20)
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

        isplit = 3

        call get_grid_dim_xy(NX_L,NY_L,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting horizontal domain dimensions'
           go to 999
        endif

        call get_laps_dimensions(NZ_L,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting vertical domain dimension'
           go to 999
        endif

        call get_meso_sao_pirep(N_MESO,N_SAO,N_PIREP,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting N_PIREP'
           go to 999
        endif

        call get_maxstns(maxstns,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting maxstns'
           go to 999
        endif

        max_cld_snd = maxstns + N_PIREP
          
        call laps_cloud_deriv(i4time,
     1                  NX_L,NY_L,
     1                  NZ_L,
     1                  N_PIREP,
     1                  maxstns,
     1                  max_cld_snd,
     1                  i_diag,
     1                  n_prods,
     1                  iprod_number,
     1                  isplit,
     1                  j_status)


        call put_derived_wind_prods(NX_L,NY_L,NZ_L           ! Input
     1          ,NX_L,NY_L,NZ_L                              ! Input (sic)
     1          ,max_radars,r_missing_data                   ! Input
     1          ,i4time)                                     ! Input


999     continue

        end

