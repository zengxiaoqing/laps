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
      subroutine read_gvap (filename, nstations, path_to_gvap8,
     1     path_to_gvap10,
     1     lat,lon, wt,w1,w2,w3,
     1     nn, istatus)

c     routine to read goes watervapor for system testing
c     author Dan Birkenheuer
c     2/8/99
c     
c     paramter list
c     filename, name of the wisconsin file
c     nstations, dimension of variables
c     lat, latitude of each station
c     lon, longitude of each station
c     wt, total water
c     w1, layer 1 water
c     w2, layer 2 water
c     w3, layer 3 water
c     nn, number of stations with real data for consideration
c
      implicit none

c     input variables
      character*9 filename
      integer nstations,nn,istatus,idummy,idummy2
      character*256 path_to_gvap8,path_to_gvap10
      real lat(nstations)
      real lon(nstations)
      real wt(nstations)
      real w1(nstations)
      real w2(nstations)
      real w3(nstations)

c     internal variables

      integer i
      real dummy
      integer ptg_index
      character*512 const_file
      integer cf
      character*9 filefound
      character*120 extension
      integer extension_index
      integer time_diff         !time allowed for latency (sec)

      time_diff = 3600*4        ! 4-hour latency allowed

      call s_len(path_to_gvap8, ptg_index, istatus)

c     reading goes 8

c     get most recent file in directory

      call get_newest_file (filename, time_diff,
     1     path_to_gvap8,ptg_index,filefound,
     1     extension, extension_index, istatus)

      if(istatus.ne.1) then     !failure in getting file
         return
      endif

      const_file = path_to_gvap8(1:ptg_index)//'20'//filefound
     1     //'.'//extension(1:extension_index)

      call s_len(const_file, cf, istatus)

      write(6,*) 'opening file ',const_file(1:cf)
      
      open(22, file=const_file(1:cf),form='formatted',
     1     status='old',err = 668)
      read(22,*,end=668,err=668) ! first header line is ignored
      do i = 1,nstations
         read(22,*,end=665,err=665) idummy,idummy,lat(i),lon(i),
     1        idummy,idummy,wt(i), w1(i),w2(i),w3(i)

      enddo

 665  close (22)
      nn = i-1
      if (nn .eq. 0) go to 666
      write(6,*) nn, ' number of records read GOES 8'
      istatus = 1
c     write(6,*) (wt(i),i=1,nn)

      go to 669

 668  close (22)
      write(6,*) 'failed reading GOES 8'
      nn = 1

 669  continue

c     reading goes 10

c     get most recent file in directory

      call get_newest_file (filename, time_diff,
     1     path_to_gvap10,ptg_index,filefound,
     1     extension, extension_index, istatus)

      if(istatus.ne.1) then     !failure in getting file
         return
      endif

      const_file = path_to_gvap10(1:ptg_index)//'20'//filefound
     1     //'.'//extension(1:extension_index)

      call s_len(const_file, cf, istatus)

      write(6,*) 'opening file ',const_file(1:cf)
      
      open(23, file=const_file(1:cf),form='formatted',
     1     status='old',err = 666)
      read(23,*,end=666,err=666) ! first header line is ignored
      do i = 1,nstations
         read(23,*,end=667,err=667) idummy,idummy,lat(i),lon(i),
     1        idummy,idummy,wt(i), w1(i),w2(i),w3(i)

      enddo


 667  close(23)

      nn = i-1
      if (nn .eq. 0) go to 666
      write(6,*) nn, ' number of records read GOES 10'
      istatus = 1
c      write(6,*) (wt(i),i=1,nn)

      return

 666  if(nn.eq.0) then
         istatus = 0
         write(6,*) ' no available gvap data'
      else
         write(6,*) 'GOES 10 not avail but GOES 8 data detetected'
         istatus = 1
      endif

      return

      end

