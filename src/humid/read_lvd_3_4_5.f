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
cdis cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis


        subroutine read_lvd_3_4_5(dir,i4time,ch3,ch4,ch5,
     1   ii,jj,kk,ngoes, istatus)


c       $log: read_lvd_3_4_5.for,v $
c revision 1.1  1996/08/30  20:57:17  birk
c initial revision
c

        implicit none

c parameter variables

        integer ii,jj,kk
        real ch3(ii,jj),ch4(ii,jj),ch5(ii,jj)
        integer i4time
        integer ngoes
        integer istatus
        character*256 dir



        character*31 ext
        integer lapsp(1)
        character*3     var(1)
        character*4     lvl_coord(1)
        character*10    units(1)
        character*125   comment(1)
c        integer len

        integer i,j


        do i = 1,ii
        do j = 1,jj
        ch3(i,j) = 0.0
        ch4(i,j) = 0.0
        ch5(i,j) = 0.0
        enddo
        enddo


c real laps data for chan 3

        

        ext = 'lvd'
c        call get_directory('lvd',dir,len)
        var(1) = 's4a'
        lapsp(1) = 0


        call read_laps(i4time,i4time, dir,
     1  ext,ii,jj,
     1  1,1,var,lapsp,
     1  lvl_coord,units,comment,ch3,istatus)

        if (istatus.ne.1) then
           write(6,*) 'Error reading channel 3'
           return
        endif


c fill ngoes parameter from common

        if(comment(1)(5:5) .eq. '8') ngoes = 8
        if(comment(1)(5:5) .eq. '9') ngoes = 9
        if(comment(1)(5:5) .eq. 'a') ngoes = 10


c real laps data for chan 4

c        dir = '../lapsprd/lvd/'
        ext = 'lvd'
c        call get_directory('lvd',dir,len)
        var(1) = 's8a'
        lapsp(1) = 0


        call read_laps(i4time,i4time, dir,
     1  ext,ii,jj,
     1  1,1,var,lapsp,
     1  lvl_coord,units,comment,ch4,istatus)

        if (istatus.ne.1) then
           write(6,*) 'Error reading channel 4'
           return
        endif

c real laps data for chan 5

c        dir = '../lapsprd/lvd/'
        ext = 'lvd'
c        call get_directory('lvd',dir,len)

        var(1) = 'sca'
        lapsp(1) = 0


        call read_laps(i4time,i4time, dir,
     1  ext,ii,jj,
     1  1,1,var,lapsp,
     1  lvl_coord,units,comment,ch5,istatus)


        if (istatus.ne.1) then
           write(6,*) 'Error reading channel 5'
           return
        endif





        return
        end



