cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps
cdis
cdis    this software and its documentation are in the public domain and
cdis    are furnished "as is."  the united states government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  they assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  if significant modifications or enhancements
cdis    are made to this software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis


        subroutine read_lvd_3_4_5(i4time,ch3,ch4,ch5,
     1   ii,jj,kk,istatus)


c       $log: read_lvd_3_4_5.for,v $
c revision 1.1  1996/08/30  20:57:17  birk
c initial revision
c

        implicit none

c parameter variables

        integer ii,jj,kk
        real ch3(ii,jj),ch4(ii,jj),ch5(ii,jj)
        integer i4time
        integer istatus



        character*50 dir
        character*31 ext
        integer lapsp(1)
        character*3     var(1)
        character*4     lvl_coord(1)
        character*10    units(1)
        character*125   comment(1)
        integer len

        integer i,j


        do i = 1,ii
        do j = 1,jj
        ch3(i,j) = 0.0
        ch4(i,j) = 0.0
        ch5(i,j) = 0.0
        enddo
        enddo


c real laps data for chan 3

        
c        dir = '../lapsprd/lvd/'
        ext = 'lvd'
        call get_directory('lvd',dir,len)
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

c real laps data for chan 4

c        dir = '../lapsprd/lvd/'
        ext = 'lvd'
        call get_directory('lvd',dir,len)
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
        call get_directory('lvd',dir,len)

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



