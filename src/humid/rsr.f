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

       subroutine rsr (i4time_in, rad, ii,jj,kk,ngoes, 
     1                 istatus)


       implicit none


       integer ii,jj,kk, i4time_in, istatus
       real rad (ii,jj,kk)
       real local_rad_data (ii,jj,19)
       integer ngoes
       character*150 ext
       character*256 dir
       integer len
       character*3 var_lsr(19)
       integer lvl_lsr(19)
       character*4 lvl_coord(19)
       character*10 units(19)
       character*125 comment(19)
       character*9 filename1
       integer  i,j,k
       integer i4time



       istatus = 0 ! bad




       ext = 'lsr'
       call get_directory('lsr',dir,len)

c     install new changes for revise satellite path

 
       if (ngoes.eq.8) then
          dir = dir(1:len)//'goes08/'
          len = len + 7
       elseif (ngoes.eq.9) then
          dir = dir(1:len)//'goes09/'
          len = len + 7
       elseif (ngoes.eq.10) then
          dir = dir(1:len)//'goes10/'
          len = len + 7
       elseif (ngoes.eq.11) then
          dir = dir(1:len)//'goes11/'
          len = len + 7
       elseif (ngoes.eq.12) then
          dir = dir(1:len)//'goes12/'
          len = len + 7
       endif

c     insertion of wait for data   10/8/02 db
       call wait_for_data (
     1      dir,                !path to data
     1      i4time_in-3600/2,   !data lookback (use data after this time)
     1      5,                  !frequency to check for data arrival (sec)
     1      600,                !amount of time to wait before giving up 
     1      3600*4,             !if data are not there at this lookback
                                !dont bother to even wait, assume no chance
     1      istatus)

       if (istatus.ne.1) then
          write (6,*) 'rsr.f:: Failure in waitfordata'
       endif


       call get_latest_file (dir,i4time_in,filename1,istatus)

       if (istatus.ne.1)then
          write(6,*) 'rsr.f:: Failure in call to get_latest_file'
          return
       endif

       write (6,*) ' ' 
       write (6,*) ' ' 
       write (6,*) ' ' 
       write (6,*) 'Directory path is: ', dir
       write (6,*) 'Attempting: ', filename1
c     convert filename to i4time
       call i4time_fname_lp (filename1,i4time,istatus)

       write (6,*) 'Sounder time difference (sec)', i4time-i4time_in

c     put in test for old data, reject data if older than 1 hour.

       if( i4time-i4time_in .gt. 3600) then
         write (6,*) 'Sounder data is too old to use, GT one hour old'
         return
       endif

       do k = 1,19
          var_lsr(k) = ' '
          write (var_lsr(k), 23) k
 23       format ('s',i2.2)
          units(k) = 'Radiance'
          lvl_lsr(k) = 0
          lvl_coord(k) = 'AGL'
       enddo


       call read_laps(i4time,i4time, dir,
     1      ext,ii,jj,
     1      19,19,var_lsr,lvl_lsr,
     1      lvl_coord,units,comment,local_rad_data,istatus)

       if (istatus.ne.1) then
          write(6,*) 'Error reading radiance data ... returning'
          return
       endif

c fill ngoes from comment line

       if(comment(1)(5:5) .eq. '8') ngoes = 8
       if(comment(1)(5:5) .eq. '9') ngoes = 9
       if(comment(1)(5:5) .eq. 'a') ngoes = 10
       if(comment(1)(5:5) .eq. 'b') ngoes = 11
       if(comment(1)(5:5) .eq. 'c') ngoes = 12


       do k = 1,kk
          do j = 1,jj
             do i = 1,ii
                rad(i,j,k) = local_rad_data(i,j,k)
                if (rad(i,j,k) .le. 0.0) then
                   write(6,*) 'Zero error ',rad(i,j,k), i,j,k,filename1
                endif
             enddo
          enddo
       enddo

       istatus = 1

       return
       end

