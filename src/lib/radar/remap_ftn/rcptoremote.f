cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
      Subroutine rcp_to_remote(i4time,ext
     1                   ,destination_node,destination_dir)

!     You guessed it, copy forbidden data outside the LAB

      character*150 dir
      character*(*) destination_dir
      character*(*) destination_node
      character*200 command
!     character*200 remote_command

      character*9 gtime
      character*31 ext,ext_in

      character*91 full_fname
      integer*4 len_fname

      write(6,*)' Subroutine rcp_to_remote'

      call s_len (destination_node,len_node)

      call s_len (destination_dir,len_dir)

!     Get ascii 9 time
      call make_fnam_lp(i4time,gtime,istatus)

!     Get directory to copy file from
      call downcase(ext,ext_in)
      call s_len(ext_in,len_ext)
      dir='../lapsprd/'//ext_in(1:len_ext)//'/'

!     Get full file name of file being copied
      call cvt_fname_data(dir,gtime,ext_in,full_fname,istatus)
      call s_len(full_fname,len_fname)

      write(6,*)full_fname

!     Generate full command

!     Option 1 - Use simple rcp command
      command = 'gzip '//full_fname(1:len_fname)//';'
     1             //' rcp '//full_fname(1:len_fname)//'* '
     1             //destination_node(1:len_node)
     1             //destination_dir(1:len_dir)

!     Option 2 - Use Dan's tar command
!     remote_command = '(cd '//destination_dir(1:len_dir)
!    1                       //'; tar -xf -)'
!     command = 'gzip '//full_fname(1:len_fname)
!    1                 //' | 'tar -cf | rsh '
!    1                 //destination_node(1:len_node)
!    1                 //remote_command 

      write(6,*)' Calling system with this command:'
      write(6,*)command

      call system(command)

      return
      end
