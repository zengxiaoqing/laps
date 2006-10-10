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

!       1997 Jul      Ken Dritz        Added call to get_grid_dim_xy.
!       1997 Jul      Ken Dritz        Pass NX_L, NY_L to ingest_pro.

        character*9 a9_time
        character*8 c8_project,c8_blp_format,c8_blp_format_in
        character*3 ext

        call GETENV('LAPS_A9TIME',a9_time)
        call s_len(a9_time,ilen)

        lun_out = 1

        call get_c8_project(c8_project,istatus)
        if(istatus .ne. 1)goto999

        call get_c8_blpfmt(c8_blp_format_in,istatus)
        if(istatus .ne. 1)goto999

        if(ilen .eq. 9)then
            write(6,*)' systime (from env) = ',a9_time
            call i4time_fname_lp(a9_time,i4time,istatus)
        else
            call get_systime(i4time,a9_time,istatus)
            if(istatus .ne. 1)go to 999
            write(6,*)' systime = ',a9_time
        endif

        call get_grid_dim_xy(NX_L,NY_L,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting horizontal domain dimensions'
           go to 999
        endif

!       if(i4time .eq. (i4time / 3600) * 3600)then
            write(6,*)
            write(6,*)' Running WPDN (NIMBUS/WFO) profiler ingest'       
            call ingest_pro(i4time,NX_L,NY_L,lun_out,j_status)
            write(6,*)' Return from WPDN (NIMBUS/WFO) profiler ingest'              

!       else
!           write(6,*)' Not on the hour, no WPDN profiler ingest run'       

!       endif

        if(c8_blp_format_in .eq. 'default')then
            c8_blp_format = c8_project
        else
            c8_blp_format = c8_blp_format_in
        endif            

        if(c8_blp_format .eq. 'RSA')then
            write(6,*)
            write(6,*)' Running RSA/LDAD local wind profile ingest '
            call ingest_rsapro(i4time,NX_L,NY_L,lun_out,j_status)
            write(6,*)' Return from RSA/LDAD local wind profile ingest'
        elseif(c8_blp_format .eq. 'WFO' .or. 
     1         c8_blp_format .eq. 'MADIS')then
            write(6,*)
            write(6,*)' Running MADIS (WFO) Multi-agency profile ingest'       
            ext = 'pro'
            call ingest_madis_map(i4time,NX_L,NY_L,ext,lun_out
     1                           ,j_status)
            write(6,*)' Return from MADIS (WFO) MAP ingest'
        else
            write(6,*)
            write(6,*)' Running BLP (NIMBUS) local profiler ingest'
            call ingest_blppro(i4time,NX_L,NY_L,lun_out,j_status)
            write(6,*)' Return from BLP (NIMBUS) local profiler ingest'
        endif

        write(6,*)
        write(6,*)' Running VAD (NIMBUS) ingest'
        call ingest_vad(istatus)
        write(6,*)' Return from VAD (NIMBUS) ingest'

 999    continue

        close(lun_out)

        end


       subroutine get_pro_parms(c8_blp_format,istatus)

!      This subroutine and namelist isn't being used at the present time.

       character*8 c8_blp_format

       namelist /pro_nl/ c8_blp_format
 
       character*150 static_dir,filename
 
       call get_directory('static',static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/ingest_pro.nl'
 
       open(51,file=filename,status='old',err=900)
       read(51,pro_nl,err=901)
       close(51)

       print*,'success reading pro_nl in ',filename
       write(*,pro_nl)

       istatus = 1
       return

  900  print*,'error opening file ',filename
       istatus = 0
       return

  901  print*,'error reading pro_nl in ',filename
       write(*,pro_nl)
       istatus = 0
       return

       end
