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
      
