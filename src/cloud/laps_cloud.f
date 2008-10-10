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

        use mem_namelist, ONLY: read_namelist_laps

        use mem_namelist, ONLY: NX_L, NY_L, nk_laps, maxstns, N_PIREP       

        integer j_status(20),iprod_number(20)
        character*150 static_dir,filename
        character*9 a9time

!       Read global parameters into module memory structure
        call get_directory('static',static_dir,len_dir)
        filename = static_dir(1:len_dir)//'/nest7grid.parms'
        call read_namelist_laps('lapsparms',filename)

!       Read cloud parameters into module memory structure
        filename = static_dir(1:len_dir)//'/cloud.nl'
        call read_namelist_laps('cloud_anal',filename)

        call get_systime(i4time,a9time,istatus)
        if(istatus .ne. 1)go to 999

        write(6,*)' systime = ',a9time

        isplit = 1

        max_cld_snd = maxstns + N_PIREP
          
        call laps_cloud_sub(i4time,
     1                  NX_L,NY_L,
     1                  nk_laps,
     1                  N_PIREP,
     1                  maxstns,
     1                  max_cld_snd,
     1                  i_diag,
     1                  n_prods,
     1                  iprod_number,
     1                  isplit,
     1                  j_status)

999     continue

        end

 



