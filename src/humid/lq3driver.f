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
      program lq3driver

c     routine to automatically execute routines lpw_driver1a & b for the laps
c     scheduler.
c
c     note, it is essential that this module be executed in the laps exe
c     path position ... its path location is vital to its proper execution
c
c
c     it is also strongly recommended that this and the associated moisture
c     modules be compiled with the highest permissible optimization.
c     especially the satellite structure inclusion can take up to 1/2 hour
c     of cpu time when run without optimization.

      implicit none
      
      include 'lapsparms.cmn'
      include 'grid_fname.cmn'
      
      integer*4
     1     ii,jj,kk,
     1     i4time,
     1     istatus,
     1     jstatus(3)
      
      real mdf
      integer lct
      
      character*9 filename
      character*200 fname
      integer len

c      grid_fnam_common = 'nest7grid'
      call get_laps_config(grid_fnam_common,istatus)
c     call get_laps_config('nest7grid',istatus)
      
      ii = nx_l_cmn
      jj = ny_l_cmn
      kk = nk_laps
      mdf = r_missing_data_cmn
      lct = laps_cycle_time_cmn
      
      call get_directory('etc',fname,len)
      print *,fname(1:len)
      open(11,file=fname(1:len)//'systime.dat',status='unknown')
      read (11,*)i4time
      read (11,22) filename
 22   format (1x,a9)
      close (11)
      
c     convert filename to i4time
 
      call i4time_fname_lp (filename,i4time,istatus)
      
      call lq3_driver1a (i4time,ii,jj,kk,mdf,lct,jstatus)
      
      write(6,*) 'lq3, lh3, and lh4 (1=success)'
      write(6,*) jstatus, ' output matrix'
      
      stop
      end
      
