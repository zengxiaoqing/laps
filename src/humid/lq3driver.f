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
      program lq3driver

c     routine to automatically execute routines lpw_driver1a for the laps
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
c     this software has recently been modified to place all of the 
c     state variable I/O routines in the driver or calling routine, thus
c     passing all state variables to the subroutine lq3_driver1a.f via
c     a parameter list.  
c
c     the purpose of this modification is to allow a future large scale driver
c     routine to call this algroritm along with others so that the entire
c     laps system can be in one executable, a prerequisite for parallel 
c     processing.
c
c     there are some current caveats that one should be aware of.  some 
c     variables in the calling routine such as this driver, contain information
c     that is changed in the called routine.  for this reason it is recommended
c     that duplicate variables be used for the call so that original state
c     variables can be kept in a pristine manner.  At this time the following
c     variables are potentially modified by the call to lq3_driver1a.f the
c     main moisture subroutine:
c     t (surface temperature)
c     td (surface dew point)
c     p (surface pressure)
c     at this time all other variables are returned in their original state.

      USE module_sfc_structure

      implicit none
      
c      include 'lapsparms.cmn'
      include 'grid_fname.cmn'

      
      integer
     1     ii,jj,kk,
     1     i4time,
     1     istatus,
     1     jstatus(3)
      
      real mdf
      integer lct

c     lat lon variables

      real, dimension (:,:), allocatable :: lat, lon

c     pressure
      real, dimension (:,:,:), allocatable :: p_3d

c     cloud variables
      
      integer c_istatus,i4timep
      real, dimension (:,:,:), allocatable :: cg
      

c     lt1 variables

      character, dimension(:), allocatable :: lvl_coordlt1*4
      character, dimension(:), allocatable :: commentlt1*125
      character, dimension(:), allocatable :: unitslt1*10
      character, dimension(:), allocatable :: varlt1*3
      real, dimension(:,:,:),  allocatable :: lt1dat
      integer t_istatus, save_i4time
      integer,   dimension(:), allocatable :: lvllm

c     specific humidity background variable "data"

      real, dimension(:,:,:), allocatable :: data

      character
     1     dirlt1*250,dir*250,rhdir*250,dirpw*250,dir3*250,
     1     extlt1*31,ext*50,rhext*50,extpw*50,ext3*50

      data extlt1/'lt1'/

c     surface variables preserved and passed
      real, dimension(:,:), allocatable :: gt
      real, dimension(:,:), allocatable :: gp
      real, dimension(:,:), allocatable :: gtd

   
      character*256  directory
      real rspacing_dum
      character*125 comment_2d
      character*10  units_2d
      character*3 var_2d
      integer len_dir
      character*200 fname
      real factor
      
      character*9 filename

      integer len

c  code



      call get_directory(extpw,dirpw,len)
      call get_directory(ext3,dir3,len)
      call get_directory(extlt1,dirlt1,len)
      call get_directory(ext,dir,len)
      call get_directory(rhext,rhdir,len)

      call get_laps_config(grid_fnam_common,istatus)

      if(istatus .ne. 1)then
         write(6,*)' error in get_laps_config'
         stop
      endif







c     get horizontal dimensions ii,jj

      call get_grid_dim_xy(ii,jj,istatus)
      if (istatus.ne.1) then
         write (6,*) 'error in grid_dim_xy'
         stop
      endif






c     get vertical dimension kk

      call get_laps_dimensions(kk,istatus)
      if (istatus.ne.1) then
         write (6,*) 'error in laps_dimensions (k)'
         stop
      endif



c     allocate arrays now dimensions known
      allocate (lat(ii,jj),lon(ii,jj))
      allocate (lt1dat(ii,jj,kk))
      allocate (lvl_coordlt1(kk))
      allocate (commentlt1(kk))
      allocate (unitslt1(kk))
      allocate (varlt1(kk))
      allocate (lvllm(kk))
      allocate (p_3d(ii,jj,kk)) ! pressure allocated here (needed for lt1 get)
      allocate (data(ii,jj,kk))
      allocate (cg(ii,jj,kk))
      allocate (gt(ii,jj))
      allocate (gp(ii,jj))
      allocate (gtd(ii,jj))
c     everything after this point can "go to 666" for deallocate stop









c     get the missing data flag

      call get_r_missing_data(mdf,istatus)
      if (istatus.ne.1) then
         write (6,*) 'error in get_r_missing_data value'
         go to 666
      endif










c     get laps cycle time variable

      call get_laps_cycle_time(lct,istatus)
      if (istatus.ne.1) then
         write (6,*) 'error in get_laps_cycle_time'
         go to 666
      endif








c     **** obtain lat lons for domain


            
c      grid_fnam_common = 'nest7grid' ! used in get_directory to modify
                                ! extension based on the grid domain
      ext = 'nest7grid'
      
c     get the location of the static grid directory
      call get_directory(ext,directory,len_dir)
      
      var_2d='lat'
      call rd_laps_static (directory,ext,ii,jj,1,var_2d,
     1     units_2d,comment_2d,
     1     lat,rspacing_dum,istatus)
      if(istatus .ne. 1)then
         write(6,*)' error reading laps static-lat'
         go to 666
      endif

      call check_nan2 (lat,ii,jj,istatus)
      if (istatus.ne.1) then
         write(6,*) 'NaNs in lat file  abort'
         stop
      endif
      
      var_2d='lon'
      call rd_laps_static (directory,ext,ii,jj,1,var_2d,
     1     units_2d,comment_2d,
     1     lon,rspacing_dum,istatus)
      if(istatus .ne. 1)then
         write(6,*)' error reading laps static-lon'
         go to 666
      endif
      
      call check_nan2 (lon,ii,jj,istatus)
      if (istatus.ne.1) then
         write(6,*) 'NaNs in lon file  abort'
         go to 666
      endif   











c     get the filename and build the i4time

      call get_directory('etc',fname,len)
      print *,fname(1:len)
      open(11,file=fname(1:len)//'systime.dat',status='unknown')
      read (11,*)i4time
      read (11,22) filename
 22   format (1x,a9)
      close (11)
      
c     convert filename to i4time
 
      call i4time_fname_lp (filename,i4time,istatus)












c     read in laps 3-d pressure and upper air data

      save_i4time = i4time  !   save the i4time for the background call
      varlt1 = 't3 '







c     read in the 3-d pressure field
      call get_pres_3d (i4time,ii,jj,kk,p_3d,istatus)

      if(istatus.ne.1) then
         write (6,*) 'error in getting 3d pressures'
         write (6,*) 'aborting'
         istatus = 0
         go to 666
      endif

c     dependence here now remains to one dimension

      p_3d = p_3d * 0.01

      lvllm (1:kk)  = int( p_3d (1,1,1:kk))  ! convert to hpa

      p_3d = p_3d / 0.01  ! restore p_3d to pristine form








c     read in the laps 3-d temperature field

      call read_laps (i4time,i4time,
     1     dirlt1,
     1     extlt1,
     1     ii,jj,kk,kk,
     1     varlt1,
     1     lvllm,
     1     lvl_coordlt1,
     1     unitslt1,
     1     commentlt1,
     1     lt1dat,
     1     t_istatus)
      if (t_istatus.ne.1) then
         print*, 'no lt1 quality control performed...'
         print*, 'missing 3-d temp data'
         write(6,*) 'ABORTING MOISTURE RUN...!!!'
         istatus = 0            ! failure
         go to 666
      endif
      















c     read in background humidity data (sh)

c     initialize data field with special missing data flag for sh

      data = -1.e+30
      
      call get_modelfg_3d(i4time,'sh ',ii,jj,kk
     1     ,data,istatus)
      
      if (istatus.ne.1) then 
         write (6,*) 'getting background field failed... abort'
         go to 666
      endif

      call check_nan3 (data,ii,jj,kk,istatus)
      if (istatus.ne.1) then
         write(6,*) 'NaN detected from RUC/MAPS...abort'
         go to 666
      endif














c     ****   laps cloud data. used for cloud, bl, goes

      call mak_cld_grid (i4time,i4timep,cg,ii,jj,kk,
     1        lct,c_istatus)

      call check_nan3 (cg,ii,jj,kk,istatus)
      if (istatus.ne.1) then
         write(6,*) 'NaN detected from Cloud Grid...ABORT'
         go to 666
      endif










c     get surface temperature and save pristine array (gt)

      call glst(i4time,gt,ii,jj,istatus)
      if(istatus.ne.1) go to 666

      call check_nan2 (gt,ii,jj,istatus)

      if(istatus.ne.1) then
         write(6,*) 'NaN detected in var:t  routine:lsin.f'
         go to 666
      endif
     














c     get surface pressure  (gp)

      call glsp(i4time,gp,ii,jj,istatus)
      if(istatus.ne.1) go to 666

      call check_nan2 (gp,ii,jj,istatus)

      if(istatus.ne.1) then
         write(6,*) 'NaN detected in var:p  routine:lsin.f'
         go to 666
      endif









c     obtain surface td value

      call glstd(i4time,gtd,ii,jj,istatus)
      if(istatus.ne.1) go to 666

      call check_nan2 (gtd,ii,jj,istatus)

      if(istatus.ne.1) then
         write(6,*) 'NaN detected in var:td  routine:lsin.f'
         go to 666
      endif

  







      

c     restore i4time from saved time

      i4time = save_i4time












c     call the main humidity algorighm with filled state variables

      call lq3_driver1a (i4time,ii,jj,kk,mdf,lat,lon,p_3d,
     1     lt1dat,lvllm,data,cg,c_istatus,i4timep,gt,gp,gtd,
     1     lct,t_istatus,jstatus)
      
      write(6,*) 'lq3, lh3, and lh4 (1=success)'
      write(6,*) jstatus, ' output matrix'



c     FINISHED

 666  continue  !  stop statements go here to deallocate

c     deallocate all allocated arrays
      deallocate (lvl_coordlt1,commentlt1,unitslt1,varlt1)
      deallocate (lat,lon)
      deallocate (lt1dat)
      deallocate (p_3d,lvllm)
      deallocate (data)
      deallocate (cg)
      deallocate (gt)
      deallocate (gp)
      deallocate (gtd)
      
      end
      
