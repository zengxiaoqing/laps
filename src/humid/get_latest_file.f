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
        subroutine get_latest_file (path,i4time,filename,istatus)


c       $log: get_latest_file.for,v $
c revision 1.1  1996/08/30  20:45:49  birk
c initial revision
c

c this routine returns the latest filename in a given path sans the .ext.
c it is useful to determine i4times from existing files on the system.

c author dan birkenheuer

c date: 8/14/96

        implicit none

        character*256 path
        character*9 filename,f_save
        integer i4time
        integer istatus

        integer numoffiles,maxfiles
        parameter (maxfiles = 3000)
        character*256 c_filenames(maxfiles),dum
        integer i4time_test,i,mintime




        call get_file_names (path, numoffiles,
     1  c_filenames, maxfiles, istatus)

        if (istatus.eq.1 .and. numoffiles.gt.0) then

        mintime = 1000000000

        do i = 1,numoffiles

        dum = c_filenames(i)

        filename = dum (index(dum,' ')-13:index(dum,' ')-4)

        call i4time_fname_lp (filename, i4time_test, istatus)

        mintime = min(abs(i4time_test-i4time),mintime)

        if(mintime.eq.abs(i4time_test-i4time) ) f_save = filename

        enddo

        filename = f_save

        if (mintime.lt.3600) return

        endif

        istatus = 0

        return
        end

