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
       endif


       call get_latest_file (dir,i4time_in,filename1,istatus)

       if (istatus.ne.1) return

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

       do k = 1,kk
          do j = 1,jj
             do i = 1,ii
                rad(i,j,k) = local_rad_data(i,j,k)
             enddo
          enddo
       enddo

       istatus = 1

       return
       end

