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
      subroutine lq3_driver1a (i4time,ii,jj,kk,mdf,lct,jstatus)



      implicit none

      include 'grid_fname.cmn'

c     parameter variables

      integer ii,jj,kk
      real mdf
      integer lct
      
      
      integer*4
     1     jstatus(3)
      
      
      integer
     1     istatus,
     1     t_istatus,
     1     i4time,
     1     i4timep,
     1     c_istatus,
     1     save_i4time,
     1     ramsi4time
      
      real
     1     ssh2                 !function
c     1        make_ssh !function type
      
      
      real*4
     1     data(ii,jj,kk),
     1     tpw(ii,jj),
     1     tempsh
      
      integer kstart (ii,jj)
      real qs(ii,jj)
      real ps(ii,jj)
      
      real bias_one

      real*4 mask (ii,jj),cg(ii,jj,kk)
      
      real*4 lt1dat(ii,jj,kk)
      
      character
     1     dirlt1*250,dir*250,rhdir*250,dirpw*250,dir3*250,
     1     extlt1*31,ext*50,rhext*50,extpw*50,ext3*50,
     1     varlt1(kk)*3,
     1     lvl_coordlt1(kk)*4,
     1     unitslt1(kk)*10,
     1     commentlt1(kk)*125
      
c     lat lon variables
   
      character*256  directory
c      character*256 grid_fnam_common
      real lat(ii, jj), lon(ii, jj)
      real rspacing_dum
      character*125 comment_2d
      character*10  units_2d
      character*3 var_2d
      integer len_dir
      character*200 fname
      real factor
      
c     rams stuff--------
      character*9
     1     filename,savefilename,ramsfile
      character ramsvar(kk)*3, ramslvlcoord(kk)*4,
     1     ramsunits(kk)*10, ramscomments(kk)*125
      character rams_dir*250, rams_ext*31
      
c     ------------------
      
      real data_in(ii,jj,kk), delta_moisture(kk), avg_moisture(kk)
      real data_start(ii,jj,kk)
      real data_pre_bound (ii,jj,kk)
      real diff_data(ii*jj)
      real ave,adev,sdev,var,skew,curt
      
      
      character*125 commentline
      
      
      integer*4
     1     i,j,k
      
      integer*4 counter
      
      integer*4 lvllm(kk)
      
      real*4 maps_rh(ii,jj,kk)

      character*3 desired_field
   
      real*4 plevel(kk)
      integer*4 mlevel(kk)

c    
c     gps variables

      real gps_data (ii,jj)
      real gps_w (ii,jj)
      integer istatus_gps
      
c     
c     gvap variables
c     
      real gvap_data(ii,jj)
      real gvap_w (ii,jj)
      integer istatus_gvap
      
      real pressure_of_level    !function call
      
      integer  raob_switch
      integer  raob_lookback
      integer goes_switch
      integer cloud_switch
      integer tiros_switch
      integer sounder_switch
      integer sat_skip
      integer gvap_switch
      integer time_diff         !time allowed for latency (sec)
      integer sfc_mix
      integer mod_4dda_1
      real    mod_4dda_factor
      real    t_ref
      integer gps_switch
      character*256 path_to_gvap8,path_to_gvap10,path_to_gps
      namelist /moisture_switch_nl/ raob_switch,
     1     raob_lookback, goes_switch, cloud_switch
     1     ,tiros_switch, sounder_switch, sat_skip
     1     ,gvap_switch, time_diff, gps_switch
     1     ,sfc_mix, mod_4dda_1,mod_4dda_factor,
     1     t_ref,path_to_gvap8,path_to_gvap10,path_to_gps
      
      integer len
      character*200 cdomain
      
      data extpw/'lh1'/
      data ext3/'lh2'/
      data extlt1/'lt1'/
      data ext /'lq3'/
      data rhext /'lh3'/
      
c     real eslo,esice
      
      
      
c----------------------code   ------------------
c     initialize laps field
      
      write (6,*) 'version 1.32; 7/8/99; increase cloud thresh (.6)'
      
c     call get_laps congif to fill common block used in pressure assignment
c     routine
      
      call get_directory(extpw,dirpw,len)
      call get_directory(ext3,dir3,len)
      call get_directory(extlt1,dirlt1,len)
      call get_directory(ext,dir,len)
      call get_directory(rhext,rhdir,len)

      call get_laps_config(grid_fnam_common,istatus)

c      call get_laps_config('nest7grid',istatus)
      if(istatus .ne. 1)then
         write(6,*)' error in get_laps_config'
         return
      endif


c     
c     set namelist parameters to defaults 
      cloud_switch = 1
      raob_switch = 0
      raob_lookback = 0
      goes_switch = 0
      sounder_switch = 0
      tiros_switch = 0
      sat_skip = 0
      gvap_switch = 0
      time_diff = 0
      gps_switch = 0
      sfc_mix = 0
      mod_4dda_1 = 0
      mod_4dda_factor = 0.02
      t_ref = -132.0
      path_to_gvap8 = ' '
      path_to_gvap10 = ' '
      path_to_gps = ' '
      
      call get_directory('static',fname,len)
      open (23, file=fname(1:len)//'moisture_switch.nl',
     1     status = 'old', err = 24)
      
      read(23,moisture_switch_nl,end=24)
      
      
      close (23)
      
      
      if (cloud_switch.eq.0) then
         write(6,*) 'Cloud switch off, ignore clouds'
         write(6,*) 'If available, clouds will be used in GOES adjust'
      else
         write (6,*) 'Clouds will be used in the analysis'
      endif
      
      if (raob_switch.eq.0) then
         write(6,*) 'raob switch off, ignoring raobs (.snd files)'
      else
         write (6,*) 'Raob switch on... will use raobs if present'
      endif
      
      write(6,*) 'RAOB look back set to ', raob_lookback, 'seconds'
      
      if (goes_switch.eq.0) then
         write(6,*) 'GOES switch off, ignoring goes data'
      else
         write(6,*) 'GOES switch on, attempting GOES ', goes_switch
      endif
      
      if (sounder_switch.eq.0) then
         write(6,*) 'Sounder switch off'
         write(6,*) 'Using IMAGER data only'
      else
         write(6,*) 'Sounder ON using Sounder data'
      endif
      
      if (tiros_switch.eq.0) then
         write(6,*) 'tiros switch off'
      else 
         write (6,*) 'Attempting to run tiros data from NOAA-',  
     1        tiros_switch
      endif
      
      if (tiros_switch.ne.0 .and. goes_switch.ne.0) then
         write(6,*) 'USING BOTH TIROS AND GOES DATA IN THIS RUN'
      endif
      
      if (sat_skip .eq. 1) then
         write(6,*) 'Use full resolution satellite'
      else
         write(6,*) 'Using partial satellite resolution ',sat_skip
      endif
      
      if (gvap_switch .eq. 1) then
         write(6,*) 'Using goes derived pw, assume data connection'
      else
         write(6,*) 'GVAP not used... nominal state'
      endif

      if (gps_switch .eq. 1) then
         write(6,*) 'Using GPS IPW data'
      else
         write(6,*) 'GPS data not used'
      endif
      
      if (time_diff .ne. 0) then
         write(6,*) 'GVAP latency assigned to ',time_diff, 'seconds'
      else
         write(6,*) 'NO latency assinged to GVAP data'
      endif
      
      if (sfc_mix .eq. 1) then
         write(6,*) 'Mixing moisture from sfc'
      else
         write(6,*) 'Sfc moisture field ignored'
      endif
      
      if (mod_4dda_1 .eq.1) then
         write(6,*) 'Mod 4dda active, modifying moisture on output'
         write(6,*) 'Mod 4dda factor is set to, ',mod_4dda_factor
      else
         write(6,*) 'Mod 4dda turned off ... nominal state'
      endif
      
      write(6,*) 'T_ref is set to: ',t_ref
      
      if (path_to_gvap8 .eq. ' '.and. path_to_gvap10 .eq. ' ')then
         write(6,*) 'Path to gvap not assigned, assigning gvap switch 0'
         gvap_switch = 0
      else
         write(6,*) 'Gvap switch assigned, using assigned switch'
         write(6,*) 'Path is ', path_to_gvap8, ' ',path_to_gvap10
         write(6,*) 'GVAP switch is set to ',gvap_switch
      endif
     
      
c     initialize field to lq3 internal missing data flag.
c     initialize total pw to laps missing data flag
      
      do i = 1,ii
         do j = 1,jj
            do k = 1,kk
               data (i,j,k) = -1e+30
            enddo
            tpw(i,j) = mdf
         enddo
      enddo

      write(6,*) 'running NEW T_ref change' 



      jstatus(1) = 0
      jstatus(2) = 0            !%loc(rtsys_abort_prod)
      jstatus(3) = 0



      do k = 1,kk
         lvllm(k) = nint( pressure_of_level(k)  * .01 )
      enddo


      do k = 1,kk
         plevel(k) = float ( lvllm(k)  )
         mlevel(k) = plevel(k)
      enddo

c     mark the maps gridpoints

      do j = 1,jj,(jj-1)/(jj-1)
         do  i = 1,ii,(ii-1)/(ii-1)

            mask (i,j) = 1

         enddo
      enddo

c     translate the i4time into filename to be used for this run
c     and store

      call make_fnam_lp (i4time,filename,istatus)

      savefilename = filename

      write(6,*) 'FILENAME = ',filename

c     preserve the i4time

      save_i4time = i4time

c     Get background field
      
      call get_modelfg_3d(i4time,'sh ',ii,jj,kk
     1     ,data,istatus)
      
      if (istatus.ne.1) then 
         write (6,*) 'getting background field failed... abort'
         return
      endif
      
      
      
      
      call check_nan3 (data,ii,jj,kk,istatus)
      if (istatus.ne.1) then
         write(6,*) 'NaN detected from RUC/MAPS...abort'
         return
      endif
      
      
      i4time = save_i4time
      filename = savefilename
      
      
      
c     check for negative input and warn
      
      do k = 1,kk
         do j = 1,jj
            do i = 1,ii
               
               if (data (i,j,k).lt.0.0) then
                  write(6,*) 'neg. input found, data,i,j,k'
                  write(6,*) data(i,j,k),i,j,k
                  data(i,j,k) = 0.0
               endif
               
            enddo
         enddo
         
         
      enddo
      
      
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
         return
      endif

      call check_nan2 (lat,ii,jj,istatus)
      if (istatus.ne.1) then
         write(6,*) 'NaNs in lat file  abort'
         return
      endif
      
      var_2d='lon'
      call rd_laps_static (directory,ext,ii,jj,1,var_2d,
     1     units_2d,comment_2d,
     1     lon,rspacing_dum,istatus)
      if(istatus .ne. 1)then
         write(6,*)' error reading laps static-lon'
         return
      endif
      
      call check_nan2 (lon,ii,jj,istatus)
      if (istatus.ne.1) then
         write(6,*) 'NaNs in lon file  abort'
         return
      endif      
      
      
c     open file for laps temp data
      do k = 1,kk
         varlt1(k) = 't3 '
      enddo
      
      call read_laps (save_i4time,save_i4time,
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
         return
      endif
      
      call check_nan3 (lt1dat,ii,jj,kk,istatus)
      if (istatus.ne.1) then
         write(6,*) 'NaN detected from lt1...ABORT'
         return
      endif
      
      
      
      
c     perform initialquality control check for supersaturation after ingest
      write(6,*)  'perform qc for supersaturation'
      counter = 0
      do k = 1,kk
         write (6,*)
         write (6,*) 'Level ', k, '   ', plevel(k)
         write (6,*)
         
         do j = jj,1,-1
            do i = 1,ii
               tempsh = ssh2( float(lvllm(k)) ,lt1dat(i,j,k)-273.15,
     1              lt1dat(i,j,k)-273.15, t_ref )/1000.
               
               if ( data(i,j,k)/tempsh .ge. 1.0) then
                  cdomain(i:i) = 'x'
                  if(data(i,j,k)/tempsh .gt. 1.01) then
                     cdomain(i:i) = 's'
                  endif
                  counter = counter + 1
                  data(i,j,k) = tempsh
               else
                  write (cdomain(i:i),34) int(data(i,j,k)/tempsh*10.)
 34               format (i1)
               endif
            enddo
            write(6,*) cdomain(1:ii)
         enddo
      enddo
      
      if(counter.gt.0) then
         
         write(6,*) ' '
         write(6,*) 'Questionable INPUT DATA DETECTED'
         write(6,*)  counter,' times.'
         write(6,*) ' '
      endif

c     initial check for computed nans
      call check_nan3(data,ii,jj,kk,istatus)
      if(istatus.ne.1) then
         write(6,*) 'initial data corrupt after checking for'
         write(6,*) 'supsaturation... hosed before beginning'
         write(6,*) 'var:data  routine:lq3driver1a.f'
         return
      endif

      
      
c     record total moisture
      
      do k = 1,kk
         do i = 1,ii
            do j = 1,jj
               data_in(i,j,k) = data(i,j,k)
               data_start(i,j,k) = data(i,j,k)
            enddo
         enddo
      enddo
      
      
c     ****  execute raob step if switch is on
      
      
      if(raob_switch.eq.1) then
         write (6,*) 'begin raob insertion'
         call raob_step (i4time,data,plevel, raob_lookback,
     1        lat,lon, lt1dat, ii,jj,kk)
         
         
         write(6,*) 'Reporting effects of RAOB insertion'
         
         call report_change (data_in, data, plevel,mdf,ii,jj,kk)

         do i = 1,ii
            do j = 1,jj
               do k  = 1,kk
                  data_in(i,j,k) = data(i,j,k)
               enddo
            enddo
         enddo
         
c     end report moisture change block
         
      else
         write(6,*) 'the raob switch is off... raobs skipped'
      endif
      
c     ****  get laps cloud data. used for cloud, bl, goes
      
      call mak_cld_grid (i4time,i4timep,cg,ii,jj,kk,
     1     lct,c_istatus)
      
      c_istatus = 0
      if (i4time.eq.i4timep) c_istatus = 1
      
      call check_nan3 (cg,ii,jj,kk,istatus)
      if (istatus.ne.1) then
         write(6,*) 'NaN detected from Cloud Grid...ABORT'
         return
      endif

c     ***   insert bl moisture

      do k = 1,kk
         do i = 1,ii
            do j = 1,jj
               data_pre_bound(i,j,k) = data(i,j,k)
            enddo
         enddo
      enddo

      print*, 'calling lsin'
c     insert boundary layer data
      call lsin (i4time,plevel,lt1dat,data,cg,tpw,bias_one,
     1     kstart,qs,ps,lat,lon,ii,jj,kk,istatus)

c     check fields after lsin call
      call check_nan3(data,ii,jj,kk,istatus)
      if(istatus.ne.1) then
         write (6,*) 'Nan generated in lsin'
         write (6,*) 'var:data  routine:lq3driver1a.f'
         return
      endif
      call check_nan2(tpw,ii,jj,istatus)
      if(istatus.ne.1) then 
         write(6,*) 'Nan generated in lsin'
         write(6,*) 'var:tpw   routine:lq3driver1a.f'
         return
      endif
      
      write(6,*) 'finished with routine lsin'
      if( sfc_mix.eq.1)then
         write(6,*) 'Lsin allowed to modify data field'
         
         write(6,*) 'Reporting incremental boundary layer effects'
         
         call report_change (data_in, data, plevel,mdf,ii,jj,kk)
         
         write(6,*) 'Reporting net change'
         call report_change (data_start, data, plevel, mdf, ii,jj,kk)

         do i = 1,ii
            do j = 1,jj
               do k  = 1,kk
                  data_in(i,j,k) = data(i,j,k)
               enddo
            enddo
         enddo
         
c     end report moisture change block
         
      else
         write(6,*) 'Lsin and sfc mixing step skipped'
         
         do k = 1,kk
            do i = 1,ii
               do j = 1,jj
                  data(i,j,k) = data_pre_bound(i,j,k)
                  data_in(i,j,k) = data(i,j,k)
               enddo
            enddo
         enddo
    
      endif

c     gps data inserstion step

      istatus_gps = 0

      if (gps_switch .eq. 1) then

         call process_gps (ii,jj,gps_data,gps_w,
     1        tpw,lat,lon,time_diff,
     1        path_to_gps,filename,istatus_gps)
         
c     gvap data insertion step
         
      endif

      istatus_gvap = 0
      
      if (gvap_switch.eq.1) then

         write(6,*) 
         write(6,*) 'Begin GVAP insertion setep'
         write(6,*) 
         
         call process_gvap(ii,jj,gvap_data,gvap_w,tpw,
     1        lat,lon,time_diff,
     1        path_to_gvap8,path_to_gvap10,filename,istatus_gvap)

      endif

      call check_nan3(data,ii,jj,kk,istatus)
      if(istatus.ne.1) then
         write(6,*) 'Failed in nan prior to adjust'
         write(6,*) 'var:data  routine lq3driver'
         return
      endif
         
      if(istatus_gps.eq.1 .or. istatus_gvap.eq.1) then ! apply gvap and/or gps 

      do k = 1,kk
         do j = 1,jj
            do i = 1,ii
               if(data(i,j,k).ge.0.0 .and. 
     1              (gvap_w(i,j)+gps_w(i,j)) .ne. 0.0 ) then
                  data(i,j,k) = 
     1                 gvap_w(i,j)/(gvap_w(i,j)+gps_w(i,j))
     1                 *(data(i,j,k) * gvap_data(i,j))
     1                 +
     1                 gps_w(i,j)/(gvap_w(i,j)+gps_w(i,j))
     1                 *(data(i,j,k) *gps_data(i,j))
     1                 + data(i,j,k)
               endif
               call check_nan(data(i,j,k), istatus)
               if (istatus.ne.1) then ! nan generated in computation
c                  write(6,*) 'Correcting Nan value, var:data',i,j,k
                  data(i,j,k) = data_in(i,j,k)
               endif
            enddo
         enddo
      enddo


c     CHECKING PROCESS OUTPUT

      call check_nan2(gvap_w,ii,jj,istatus)
      if(istatus.ne.1) then
         write(6,*) 'Failed in nan after adjust'
         write(6,*) 'var:gvap_w  routine lq3driver'
         return
      endif
         

      call check_nan2(gps_w,ii,jj,istatus)
      if(istatus.ne.1) then
         write(6,*) 'Failed in nan after adjust'
         write(6,*) 'var:gps_w  routine lq3driver'
         return
      endif
         

      call check_nan2(gps_data,ii,jj,istatus)
      if(istatus.ne.1) then
         write(6,*) 'Failed in nan after adjust'
         write(6,*) 'var:gps_data  routine lq3driver'
         return
      endif
         

      call check_nan2(gvap_data,ii,jj,istatus)
      if(istatus.ne.1) then
         write(6,*) 'Failed in nan after adjust'
         write(6,*) 'var:gsp_data  routine lq3driver'
         return
      endif
         

      call check_nan3(data,ii,jj,kk,istatus)
      if(istatus.ne.1)then
         write(6,*) 'Nan report from TPW processing'
         write(6,*) 'var:data  routine lq3driver1a.f'
         return
      endif
      
      write(6,*) 'Reporting changes from TPW data types'
      
      call report_change (data_in, data, plevel,mdf,ii,jj,kk)
      
      write(6,*) 'Reporting net change'
      call report_change (data_start, data, plevel, mdf, ii,jj,kk)
      
      do i = 1,ii
         do j = 1,jj
            do k  = 1,kk
               data_in(i,j,k) = data(i,j,k)
            enddo
         enddo
      enddo
      

      else
         write(6,*) 'gvap weights not applied, istatus = 0'
      endif
            

      write(6,*) 

c     make call to TIROS moisture insertion


      if (tiros_switch .ne. 0) then

         if(c_istatus.eq.1 .and. t_istatus.eq.1) then

        
            write (6,*) 'begin TIROS insertion step'
      
            call tiros (
     1           data,          ! specific humidity g/g
     1           lat,lon,       ! lat and longitude (deg)
     1           i4time,        ! i4time of run (seconds)
     1           plevel,        ! pressure hpa (laps vert grid)
     1           cg,            ! cloud array
     1           lt1dat,        ! lt1 (laps 3d temps)
     1           14,            ! satellite number
     1           ii,jj,kk       ! grid dimensions
     1           )

         else
            
            write(6,*)
            write(6,*)
            write(6,*)
            write(6,*) 'tiros moisture insertion step skipped'
            write(6,*) 'cloud or lt1 data not current'
            write(6,*) 'cannot assume clear conditions or'
            write(6,*) 'use alternate lt1.. this will create'
            write(6,*) 'forward model problems....'
            write(6,*) 'tiros moisture insertion step skipped'
            write(6,*)
            write(6,*)
            write(6,*)

         endif

      else

         write(6,*) 'tiros switch is off... tiros step skipped...'
         
      endif


c     make call to goes moisture insertion


      if (goes_switch .ne. 0 ) then

         if(c_istatus.eq.1 .and. t_istatus.eq.1) then

            write (6,*) 'begin goes insertion step'
            call goes_sbn (
     1           data,          ! 3-d specific humidity g/g
     1           lat,lon,       ! 2-d lat and longitude
     1           i4time,        ! i4time of run
     1           plevel,        ! pressure mb
     1           cg,            ! 3-e cloud field 0-1 (1=cloudy)
     1           lt1dat,        ! laps lt1 (3-d temps)
     1           goes_switch,   ! goes switch and satellite number
     1           sounder_switch, ! sounder switch, 0=imager,1=sndr
     1           sat_skip,      ! normally 1 for full resolution
     1           ii,jj,kk
     1           )
            
            write (6,*) 'GOES step complete, effects logged.'
            
            call report_change (data_in, data, plevel,mdf,ii,jj,kk)
            
            write(6,*) 'Reporting net change'
            call report_change (data_start, data, plevel, mdf, ii,jj,kk)
            
            do i = 1,ii
               do j = 1,jj
                  do k  = 1,kk
                     data_in(i,j,k) = data(i,j,k)
                  enddo
               enddo
            enddo
         

            
            
c     end report moisture change block
        

         else

            write(6,*)
            write(6,*)
            write(6,*)
            write(6,*) 'goes moisture insertion step skipped'
            write(6,*) 'cloud or lt1 data not current'
            write(6,*) 'cannot assume clear conditions or'
            write(6,*) 'use alternate lt1.. this will create'
            write(6,*) 'forward model problems....'
            write(6,*) 'goes moisture insertion step skipped'
            write(6,*)
            write(6,*)
            write(6,*)
            
         endif
         
      else
         
         write(6,*) 'goes switch is off... goes step skipped...'
         
      endif


c     *** insert cloud moisture, this section now controled by a switch

      if(cloud_switch.eq.0) then
         write(6,*) ' '
         write(6,*) 'Cloud switch  '
         write(6,*) 'Skipping cloud moistening here'
         write(6,*) ' '
         write(6,*) ' '
         c_istatus = 0          !force this to skip here
      endif
      
      if(c_istatus.ne.1 .or. i4time.ne.i4timep)then
         c_istatus = 0
         write(6,*) 'Cloud field not available for exact time'
         write(6,*) 'assume ALL CLEAR values used'
      else                      ! increase moisture based on cloud amount
         
         write(6,*) 'Saturate in cloudy areas'

         call check_nan3(cg,ii,jj,kk,istatus)
         if(istatus.ne.1)then
            write(6,*) 'Nan in cg data prior to saturation'
            return
         endif
         call check_nan3(lt1dat,ii,jj,kk,istatus)
         if(istatus.ne.1) then
            write(6,*) 'NaN in lt1dat prior to saturation'
            return
         endif

         call check_nan3(data,ii,jj,kk,istatus)
         if(istatus.ne.1)then
            write(6,*) 'Nan in data prior to saturation'
            return
         endif


         do k = 1,kk
            write(6,*) lvllm(k),'checking lvllm prior to sat'
            do j = 1,jj
               do i = 1,ii
                  
                  if(cg(i,j,k) .gt. 0.6 .and. cg(i,j,k) .lt. 1.0) then !cloudy
                     
                     tempsh = ssh2( float(lvllm(k)),
     1                    lt1dat(i,j,k)-273.15,
     1                    lt1dat(i,j,k)-273.15, t_ref )/1000.
                     data(i,j,k) = cg(i,j,k)* tempsh
     1                    +(1.-cg(i,j,k))*data(i,j,k)
                     
                  elseif (cg(i,j,k).ge.1.0) then 
c     ! still cloudy...put in for albers
                     
                     tempsh = ssh2( float(lvllm(k)) 
     1                    ,lt1dat(i,j,k)-273.15,
     1                    lt1dat(i,j,k)-273.15, t_ref )/1000.
                     data(i,j,k) = tempsh
                     
                  endif
                  
                  
               enddo
            enddo
         enddo

         call check_nan3(data,ii,jj,kk,istatus)
         if(istatus.ne.1)then
            write(6,*) 'Nan in data AFTER call to saturation'
            return
         endif

         
         write (6,*) 'Reporting cloud effects on analysis'
         call report_change (data_in, data, plevel,mdf,ii,jj,kk)
         
         write(6,*) 'Reporting net change'
         call report_change (data_start, data, plevel, mdf, ii,jj,kk)
         
         do i = 1,ii
            do j = 1,jj
               do k  = 1,kk
                  data_in(i,j,k) = data(i,j,k)
               enddo
            enddo
         enddo
         
      endif
      
  
c     mod_4dda_1 to decrease overall water in 4dda mode running at AFWA
      
      if(mod_4dda_1 .eq. 1) then ! act to decrease overall water
         
         do k=1,kk
            factor=1.-(float(k)*mod_4dda_factor)
            do j=1,jj
               do i=1,ii
                  data(i,j,k)=data(i,j,k)*factor
               enddo
            enddo
         enddo
         
         write(6,*) ' mod_4dda loop complete'
         
         call report_change (data_in, data, plevel,mdf,ii,jj,kk)
         
      endif
      
      
c     repeat quality control check for supersaturation after pre-analysis
      write (6,*)  'perform qc for supersaturation'
      counter = 0
      do k = 1,kk
         write(6,*) 'Level ',k, '    ', plevel(k)
         do j = jj,1,-1
            do i = 1,ii
               
               tempsh = ssh2( float(lvllm(k)) ,lt1dat(i,j,k)-273.15,
     1              lt1dat(i,j,k)-273.15,t_ref )/1000.
               
               if ( data(i,j,k)/tempsh .ge. 1.0) then
                  cdomain(i:i) = 'x'
                  if(data(i,j,k)/tempsh .gt. 1.01) cdomain(i:i) = 's'
                  counter = counter + 1
                  data(i,j,k) = tempsh
               elseif (data(i,j,k) .lt. 0.0) then
                  cdomain(i:i) = 'M'
                  
               else
                  write (cdomain(i:i),35) int(data(i,j,k)/tempsh*10.)
 35               format (i1)
                  
               endif
               
               
            enddo
            write(6,*) cdomain(1:ii)
         enddo
      enddo
      
      if (counter.ne.0) then
         write (6,*) 'supersaturation has been corrected, 
     1        ',counter,' times.'
      endif
      
      
c     recompute tpw including clouds and supersat corrections
      
      call int_tpw(data,kstart,qs,ps,plevel,tpw,ii,jj,kk)
      
c     place the accepted missing data flag in output field
c     sum over the entire grid for a total water sum value for 
c     QC study.
      
      tempsh = 0.0
      
      do i = 1,ii
         do j = 1,jj
            do k = 1,kk
               
               if(data(i,j,k) .lt.0.0) then
                  data(i,j,k) = mdf !  put in missing data flag if missing
               else
                  tempsh = tempsh + data(i,j,k) ! sum if good data
               endif
               
            enddo
         enddo
      enddo
      
      
      
      
c     log the amount of water vapor
      
      write (6,*) ' '
      write (6,*) ' '
      write (6,*) '***************************** '
      write (6,*) 'Average water in volume (g/g)*10000'
      write (6,*) tempsh/float(ii)/float(jj)/float(kk)*10000.
      write (6,*) '***************************** '
      write (6,*) ' '
      write (6,*) ' '
      
c     check for NaN values and Abort if found
      
      call check_nan3(data,ii,jj,kk,istatus)
      if(istatus.ne.1) then
         write(6,*) 'NaN values detected (sh array)... aborting'
         return
      endif
      
      call check_nan2(tpw,ii,jj,istatus)
      if(istatus.ne.1) then
         write(6,*) 'NaN values detected (tpw array)... aborting'
         return
      endif
      
      
c     write final 3-d sh field to disk
      commentline = 'maps with clouds and surface effects only'
      call writefile (save_i4time,commentline,mlevel,data,
     1     ii,jj,kk,istatus)
      if(istatus.eq.1)        jstatus(1) = 1
      
c     write total precipitable water field
      call write_lh4 (save_i4time,tpw,bias_one,ii,jj,istatus)
      if(istatus.eq.1) jstatus(3) = 1
      
      
c     generate lh3 file (RH true, RH liquid)
      if (t_istatus.eq.1) then
         call lh3_compress(data,lt1dat,save_i4time,lvllm,t_ref,
     1        ii,jj,kk,istatus)
         if(istatus.eq.1)        jstatus(2) = 1
      else
         print*, 'no lh3 or rh data produced...'
         print*, 'no laps 3-d temp data avail'
         jstatus(2) = 0
      endif
      
      write (6,*) 'Reporting overall changes to moisture'
      
      call report_change (data_start, data, plevel,mdf,ii,jj,kk)
      
      return
      
 24   write(6,*) 'error finding moisture switch file'
      write(6,*) 'check to see it is under'
      write(6,*) fname(1:len)//'moisture_switch.nl'
      write(6,*) 'aborting'
      
      return
      
      end
      
