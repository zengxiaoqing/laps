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
     1       c_filenames, maxfiles, istatus)
        
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

           write (6,*) 'File found is ',mintime,'(sec) old'
           
           if (mintime.lt.3600) return
           
        endif

        write (6,*) 'File found is ',mintime,'(sec) old'
        
        istatus = 0
        
        return
        end

