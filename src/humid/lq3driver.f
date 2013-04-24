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
     1     ii,jj,kk,iout,
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

      data extpw/'lh1'/
      data ext3/'lh2'/
      data ext/'lq3'/
      data rhext/'lh3'/

c     create common block for cloud_sat data (direct insert)
c     via Steve Albers new function for humidification
c     the variable with the _nl suffix is from the namelist and will be
c     assigned here since a namelist variable cannot be in a 
c     common block also
      common /cloud_sat_insert/ max_cdelrh,cf_set
      real max_cdelrh
      real cf_set   
c     common block for cloud weighting .. direct insert into func_o.f
      common /func_o_insert/ cloud_weight
      real cloud_weight

c     direct insert for radiometer weight in func_o.f
      common /radiometer/ radio_wt
      real radio_wt

c     namelist data

      integer covar_switch
      integer print_switch
      integer  raob_switch
      integer radiometer_switch
      integer  raob_lookback
      real raob_radius
      integer endian  ! 1 = big, 0 = little , big default
      integer goes_switch
      integer cloud_switch
      integer cloud_d
      real    max_cdelrh_nl
      real    cf_set_nl
      real    cloud_weight_nl
      real    radio_wt_nl
      integer tiros_switch
      integer sounder_switch
      integer sat_skip
      integer gvap_switch
      integer IHOP_flag
      integer time_diff         !time allowed for latency (sec)
      integer sfc_mix
      integer mod_4dda_1
      real    mod_4dda_factor
      real    t_ref
      integer gps_switch
      character*256 path_to_gvap12,path_to_gvap10,path_to_gps,path2covar
      namelist /moisture_switch_nl/ covar_switch,print_switch,
     1     raob_switch,radiometer_switch,
     1     raob_lookback, endian,
     1     raob_radius, goes_switch, cloud_switch, cloud_d,
     1     max_cdelrh_nl,cf_set_nl,cloud_weight_nl,radio_wt_nl
     1     ,tiros_switch, sounder_switch, sat_skip
     1     ,gvap_switch, IHOP_flag, time_diff, gps_switch
     1     ,sfc_mix, mod_4dda_1,mod_4dda_factor,
     1     t_ref,path_to_gvap12,path_to_gvap10,path_to_gps,
     1     path2covar

      
      integer i,j,k

c  code

      write (6,*) 'LQ3 Starting now'
      write (6,*)


      iout = 1 ! have subroutine generate output files



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


      write (6,*) ' ALL directory information has been read'



c     fill namelist variables---------------------------------------------
c     
c     set namelist parameters to defaults 
      covar_switch = 0
      print_switch = 0
      cloud_switch = 0
      cloud_d = 0
      max_cdelrh_nl = 0.11 ! test default; the 0.11 value should not appear in code
      cf_set_nl = 0.6 ! test default: 0.3 is better value
      cloud_weight = 0.5 ! normal default up to 3/8/12 DB
      radio_wt_nl = 1.0 ! current default
      raob_switch = 0
      radiometer_switch = 2000
      raob_lookback = 0
      endian = 1 ! big endian is default
      raob_radius = 45000.0  ! meters (45km)
      goes_switch = 0
      sounder_switch = 0
      tiros_switch = 0
      sat_skip = 0
      gvap_switch = 1
      IHOP_flag = 0 
      time_diff = 0
      gps_switch = 1
      sfc_mix = 0
      mod_4dda_1 = 0
      mod_4dda_factor = 0.02
      t_ref = -132.0
      path_to_gvap12 = ' '
      path_to_gvap10 = ' '
      path_to_gps = ' '
      path2covar = ' '
      
      call get_directory('static',fname,len)
      open (23, file=fname(1:len)//'moisture_switch.nl',
     1     status = 'old', err = 666)
      
      read(23,moisture_switch_nl,end=666)
      
      
      close (23)

      write (6,*)
      write (6,*) 'NAMELIST has been read'


c-----------------------end filling namelist ---





c     assign name list cloud_sat parameter to the variable in the common block
c     namelist variable suggested by Steve Albers for routine Cloud_sat.f only
c     this is a direct insert from the namelist into that routine via a special
c     common block that was deviced in Feb 2012
      max_cdelrh = max_cdelrh_nl
      cf_set = cf_set_nl
      cloud_weight = cloud_weight_nl
      radio_wt = radio_wt_nl


c     get horizontal dimensions ii,jj

      call get_grid_dim_xy(ii,jj,istatus)
      if (istatus.ne.1) then
         write (6,*) 'error in grid_dim_xy'
         stop
      endif

      write (6,*) 
      write (6,*) ' grid dimensions are, ',ii,jj




c     get vertical dimension kk

      call get_laps_dimensions(kk,istatus)
      if (istatus.ne.1) then
         write (6,*) 'error in laps_dimensions (k)'
         stop
      endif

      write (6,*)
      write (6,*) 'vertical dimension is: ',kk
      write (6,*)


      write(6,*) 'ALLOCATING ARRAYS NOW top level driver'
c     allocate arrays now dimensions known
      allocate (lat(ii,jj),lon(ii,jj))
      write(6,*) 'lat lon done'
      allocate (lt1dat(ii,jj,kk))
      write (6,*) 'lt1dat done'
      allocate (lvl_coordlt1(kk))
      write(6,*) 'lvl_coordlt1 done'
      allocate (commentlt1(kk))
      write(6,*) 'commentlt1 done'
      allocate (unitslt1(kk))
      write (6,*) 'unitslt1 done'
      allocate (varlt1(kk))
      write (6,*) 'varlt1 done'
      allocate (lvllm(kk))
      write (6,*) 'lvllm done'
      allocate (p_3d(ii,jj,kk)) ! pressure allocated here (needed for lt1 get)
      write (6,*) 'p_3d done'
      allocate (data(ii,jj,kk))
      write (6,*) 'data done'
      allocate (cg(ii,jj,kk))
      write (6,*) 'cg done'
      allocate (gt(ii,jj))
      write (6,*) 'gt done'
      allocate (gp(ii,jj))
      write (6,*) 'gp done'
      allocate (gtd(ii,jj))
      write (6,*) 'gtd done'
      write (6,*)
      write (6,*)
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
!     open(11,file=fname(1:len)//'systime.dat',status='unknown')
!     read (11,*)i4time
!     read (11,22) filename
!22   format (1x,a9)
!     close (11)
      call get_systime(i4time,filename,istatus)
      print *,' systime returned from get_systime = ',filename
      
c     convert filename to i4time
 
!     call i4time_fname_lp (filename,i4time,istatus)












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


c     check for "constant" input field and abort"
      
      do k = 1,kk

c     test levels independently

        do i = 1,ii
          do j = 1,jj

            if (data(1,1,k).ne.data(i,j,k)) then
c     test is passed
               go to 555 ! jump out of kill loop
            endif

          enddo ! jj
        enddo ! ii
         write(6,*) 'Data in SH input array are all EQUAL in value'
         write (6,*) 'Value detected = ', data(1,1,k), ' level ', k
         write(6,*) 'Assuming such a constant field is an ERROR'
         write(6,*) 'Aborting LQ3 run'
         go to 666
555      continue  ! made it past kill loop for all same data
         write (6,*) 'Success in variable SH field, level ',k
      enddo ! kk














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
     1     covar_switch,
     1     print_switch,
     1     cloud_switch,
     1     cloud_d,
     1     raob_switch,
     1     radiometer_switch,
     1     raob_lookback,
     1     endian,
     1     raob_radius,
     1     goes_switch,
     1     sounder_switch,
     1     tiros_switch,
     1     sat_skip,
     1     gvap_switch,
     1     IHOP_flag, 
     1     time_diff,
     1     gps_switch,
     1     sfc_mix,
     1     mod_4dda_1,
     1     mod_4dda_factor,
     1     t_ref,
     1     path_to_gvap12,
     1     path_to_gvap10,
     1     path_to_gps,
     1     path2covar,
     1     lct,iout,t_istatus,jstatus)
      
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
      
