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
      integer len_fname

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
