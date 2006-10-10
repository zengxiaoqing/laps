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
!       1997 Jul      Ken Dritz        Pass NX_L, NY_L to ingest_lrs.

        character*9 a9_time
        character*8 c8_project, c8_blpfmt, c8_blpfmt_in

        call get_systime(i4time,a9_time,istatus)
        if(istatus .ne. 1)go to 999

        call get_grid_dim_xy(NX_L,NY_L,istatus)
	if (istatus .ne. 1) then
           write (6,*) 'Error getting horizontal domain dimensions'
           go to 999
        endif

        write(6,*)
        write(6,*)' Running WPDN (NIMBUS) Rass ingest'
        call ingest_lrs(i4time,NX_L,NY_L,j_status)
        write(6,*)' Return from WPDN (NIMBUS) Rass ingest'

        call get_c8_project(c8_project,istatus)
        if(istatus .ne. 1)goto999

        call get_c8_blpfmt(c8_blpfmt_in,istatus)
        if(istatus .ne. 1)goto999

        if(c8_project .eq. 'RSA')then
            write(6,*)
            write(6,*)' Running RSA/LDAD local RASS ingest '
            call ingest_rsalrs(i4time,NX_L,NY_L,j_status)
            write(6,*)' Return from RSA/LDAD local RASS ingest'
        elseif(c8_project   .eq. 'WFO' .or. 
     1         c8_blpfmt_in .eq. 'MADIS'    )then
            write(6,*)
            write(6,*)' Running MADIS (WFO) Multi-agency profile ingest'       
            call ingest_madis_map(i4time,NX_L,NY_L,'lrs',j_status)
            write(6,*)' Return from MADIS (WFO) MAP ingest'
        else
            write(6,*)
            write(6,*)' Running BLP (NIMBUS) Rass ingest'
            call ingest_blplrs(i4time,NX_L,NY_L,j_status)
            write(6,*)' Return from BLP (NIMBUS) Rass ingest'
        endif

999     continue
        end

