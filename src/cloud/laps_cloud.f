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

        integer*4 j_status(20),iprod_number(20)
        character*9 a9_time

        call get_systime(i4time,a9_time,istatus)
        if(istatus .ne. 1)go to 999

        write(6,*)' systime = ',a9_time

        isplit = 1

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
          
        call laps_cloud(i4time,
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

999     continue

        end

 
       subroutine get_cloud_parms(l_use_vis,l_use_vis_add                ! O
     1                           ,l_use_vis_partial                      ! O
     1                           ,l_use_39,latency_co2                   ! O
     1                           ,pct_req_lvd_s8a                        ! O
     1                           ,i4_sat_window,i4_sat_window_offset     ! O
     1                           ,istatus)                               ! O

       logical l_use_vis,l_use_vis_add,l_use_vis_partial,l_use_39 
       namelist /cloud_nl/ l_use_vis, l_use_vis_add, l_use_vis_partial       
     1                    ,l_use_39, latency_co2
     1                    ,pct_req_lvd_s8a
     1                    ,i4_sat_window,i4_sat_window_offset
 
       character*150 static_dir,filename

!      Default value that can be overridden in namelist
       l_use_vis_add = .false.
       l_use_vis_partial = .true.
 
       call get_directory('static',static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/cloud.nl'
 
       open(1,file=filename,status='old',err=900)
       read(1,cloud_nl,err=901)
       close(1)

       print*,'success reading cloud_nl in ',filename
       write(*,cloud_nl)

       istatus = 1
       return

  900  print*,'error opening file ',filename
       istatus = 0
       return

  901  print*,'error reading cloud_nl in ',filename
       write(*,cloud_nl)
       istatus = 0
       return

       end




