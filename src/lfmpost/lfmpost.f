!dis   
!dis    Open Source License/Disclaimer, Forecast Systems Laboratory
!dis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
!dis    
!dis    This software is distributed under the Open Source Definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis    
!dis    In particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis    
!dis    - Redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis    
!dis    - Redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis    
!dis    - All modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis    
!dis    - If significant modifications or enhancements are made to this
!dis    software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis    
!dis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
!dis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
!dis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
!dis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
!dis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
!dis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
!dis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
!dis   
!dis 

PROGRAM lfmpost

! F90 program to post-process local forecast model output files into
!  LAPS-format (netCDF fsf/fua), Vis5D, GRIB, etc.  Currently supports
!  MM5v3 and WRFv1.2.1 for input model data.  
! 
!  This program assumes you are processing model using a grid configuration
!  identical to the grid in the specified LAPS_DATA_ROOT!
!
! Brent Shaw, Sep 2002
!   Base on original MM5POST program by same author
!
!

  USE setup
  USE mm5v3_io
  USE wrfv1_netcdf
  USE postproc_lfm
  USE constants
  USE time_utils
  USE vis5d
  USE map_utils
  USE grib
 
  IMPLICIT NONE
  INTEGER                     :: lun_data
  INTEGER                     :: t,k
  INTEGER                     :: time_index
  INTEGER                     :: status
  LOGICAL                     :: file_ready
  INTEGER                     :: funit,nbytes
  CHARACTER(LEN=255)          :: gribdir
  INTEGER                     :: gdir_len,ddir_len,gfile_len
  CHARACTER(LEN=255)          :: gribfile,gribdone
  CHARACTER (LEN=24)          :: current_time
  INTEGER                     :: laps_reftime
  INTEGER                     :: laps_valtime,laps_valtime_prev
  INTEGER                     :: period_sec
  INTEGER                     :: pass
  REAL                        :: smth
  REAL, ALLOCATABLE           :: rhodry ( : , : , : )
  REAL, ALLOCATABLE           :: cldliqcon_prs ( : , : , : )
  REAL, ALLOCATABLE           :: cldicecon_prs ( : , : , : )
  REAL, ALLOCATABLE           :: raincon_prs   ( : , : , : )
  REAL, ALLOCATABLE           :: snowcon_prs   ( : , : , : )
  REAL, ALLOCATABLE           :: graupelcon_prs( : , : , : )
  ! Declarations added for make_points
  INTEGER                     :: year
  INTEGER                     :: month
  INTEGER                     :: day
  INTEGER                     :: hour
  INTEGER                     :: minute
  INTEGER                     :: second
  INTEGER                     :: igds(18)
  INTEGER                     :: ip
  INTEGER                     :: jp
  INTEGER                     :: dir_pt
  INTEGER                     :: spd_pt
  INTEGER                     :: ceiling_pt
  INTEGER                     :: pt
  REAL                        :: t_pt
  REAL                        :: td_pt
  INTEGER                     :: rh_pt
  REAL                        :: u_pt
  REAL                        :: v_pt
  REAL                        :: dir_pt_r
  REAL                        :: spd_pt_r
  REAL                        :: cld_pt
  REAL                        :: vis_pt
  REAL                        :: pcp_pt
  REAL                        :: snow_pt
  CHARACTER(LEN=16)           :: date_str
  CHARACTER(LEN=8)            :: wx_pt
  REAL, EXTERNAL              :: bint
  REAL                        :: ceiling_pt_r
  REAL                        :: u_pt_g
  REAL                        :: v_pt_g
  INTEGER                     :: t_v5d
  CHARACTER(LEN=255)          :: flagfile
  CHARACTER(LEN=9)            :: a9time
  CHARACTER(LEN=4)            :: fcst_hhmm
  INTEGER                     :: istatus
  INTEGER                     :: startb
 
  CALL setup_lfmpost
 
  CALL mm5_to_lapstime(cycle_date,laps_reftime)

  ! Start the loop over all times

  CALL s_len(lfmprd_dir, ddir_len)

  WRITE(gribdir,'(A,"/d",I2.2,"/grib")') &
    lfmprd_dir(1:ddir_len),domain_num
  CALL s_len(gribdir,gdir_len)
    
  laps_valtime_prev = 0
  time_loop:  DO t = 1, num_times_to_proc, time_index_inc
    
    time_index = t + start_time_index - 1
    current_time = times_to_proc(t)
    ! Compute LAPS reference time
    CALL mm5_to_lapstime(current_time,laps_valtime)
    ! Open this times file for reading
    IF (mtype .EQ. 'mm5') THEN
      CALL make_data_file_name(mm5_data_root,domain_num_str,split_output,&
                             time_index-1, data_file,file_num3)
      IF (split_output) THEN 
        IF (realtime) THEN
          flagfile = TRIM(data_file) // '.done'
          INQUIRE(FILE=flagfile, EXIST=file_ready)
          IF (.NOT.file_ready) THEN
             CALL io_wait(data_file,max_wait_sec) 
          ENDIF
          OPEN(89,FILE=flagfile,STATUS='OLD')
        ENDIF
        CALL open_mm5v3(data_file,lun_data,status)     
      ELSE IF (t.EQ.1) THEN
        CALL open_mm5v3(data_file,lun_data,status)
      ENDIF

    ELSEIF(mtype.EQ.'wrf') THEN
      CALL make_wrf_file_name(lfmprd_dir,domain_num,sim_minutes(t),data_file)
      INQUIRE(FILE=data_file,EXIST=file_ready)
      IF (.NOT.file_ready) THEN
        CALL io_wait(data_file,max_wait_sec)
      ENDIF
      IF (realtime) THEN
        CALL sleep(30)  ! To make sure it is really done!
      ENDIF
      CALL open_wrfnc(data_file,lun_data,status)
    ENDIF
    IF (status .NE. 0) THEN
      PRINT *, 'Stopping due to error opening ', TRIM(data_file)
      CALL ABORT
    ENDIF
    
    ! Allocate the variables.  Allocation/deallocation of most of the arrays
    ! is done inside the time loop in case a user is running this program
    ! concurrently with the model (using split output) on a shared system.  
   
    ALLOCATE ( psfc         ( nx , ny ) )
    ALLOCATE ( tsfc         ( nx , ny ) )
    ALLOCATE ( thetasfc     ( nx , ny ) )
    ALLOCATE ( thetaesfc    ( nx , ny ) )
    ALLOCATE ( rhsfc        ( nx , ny ) )
    ALLOCATE ( tdsfc        ( nx , ny ) )
    ALLOCATE ( thetaprs     ( nx , ny , kprs ) )
    ALLOCATE ( zprs         ( nx , ny , kprs ) )
    ALLOCATE ( rhprs        ( nx , ny , kprs ) )
    ALLOCATE ( tprs         ( nx , ny , kprs ) )
    ALLOCATE ( tdprs        ( nx , ny , kprs ) )
    ALLOCATE ( shprs        ( nx , ny , kprs ) )
    ALLOCATE ( redp         ( nx , ny ) )
    ALLOCATE ( pmsl         ( nx , ny ) )
    ALLOCATE ( usfc         ( nx , ny ) )
    ALLOCATE ( uprs         ( nx , ny , kprs ) )
    ALLOCATE ( vsfc         ( nx , ny ) )
    ALLOCATE ( vprs         ( nx , ny , kprs ) )
    ALLOCATE ( wsfc         ( nx , ny ) )
    ALLOCATE ( clwmrsfc       ( nx , ny ) )
    ALLOCATE ( icemrsfc       ( nx , ny ) )
    ALLOCATE ( snowmrsfc       ( nx , ny ) )
    ALLOCATE ( rainmrsfc       ( nx , ny ) )
    ALLOCATE ( graupmrsfc       ( nx , ny ) )
    ALLOCATE ( wprs         ( nx , ny , kprs ) )
    ALLOCATE ( omprs        ( nx , ny , kprs ) )
    ALLOCATE ( cldbase      ( nx , ny ) )
    ALLOCATE ( cldtop       ( nx , ny ) )   
    ALLOCATE ( cldamt       ( nx , ny ) )
    ALLOCATE ( ceiling      ( nx , ny ) )
    ALLOCATE ( intliqwater  ( nx , ny ) )
    ALLOCATE ( totpcpwater  ( nx , ny ) )
    ALLOCATE ( max_refl     ( nx , ny ) )
    ALLOCATE ( echo_tops    ( nx , ny ) )
    ALLOCATE ( cldliqmr_prs ( nx , ny , kprs ) )
    ALLOCATE ( cldicemr_prs ( nx , ny , kprs ) )
    ALLOCATE ( rainmr_prs   ( nx , ny , kprs ) )
    ALLOCATE ( snowmr_prs   ( nx , ny , kprs ) )
    ALLOCATE ( graupelmr_prs( nx , ny , kprs ) )
    ALLOCATE ( refl_prs     ( nx , ny , kprs ) )
    ALLOCATE ( refl_sfc     ( nx , ny ) )
    ALLOCATE ( pcptype_sfc  ( nx , ny ) )
    ALLOCATE ( pcptype_prs  ( nx , ny , kprs) )
    IF (.NOT.ALLOCATED(pcp_init)) ALLOCATE ( pcp_init     ( nx , ny ) )
    ALLOCATE ( pcp_inc      ( nx , ny ) )
    ALLOCATE ( con_pcp_inc      ( nx , ny ) )
    IF (.NOT.ALLOCATED(con_pcp_init)) ALLOCATE ( con_pcp_init     ( nx , ny ) )
    IF (.NOT.ALLOCATED(pcp_tot))       ALLOCATE ( pcp_tot      ( nx , ny ) )
    IF (.NOT.ALLOCATED(con_pcp_tot))  ALLOCATE ( con_pcp_tot  ( nx , ny ) )
    IF (.NOT.ALLOCATED(snow_init))     ALLOCATE ( snow_init    ( nx , ny ) )
    ALLOCATE ( snow_inc     ( nx , ny ) )
    IF (.NOT.ALLOCATED(snow_tot))      ALLOCATE ( snow_tot     ( nx , ny ) )
    IF (.NOT.ALLOCATED(srhel))         ALLOCATE ( srhel        ( nx , ny ) )
    ALLOCATE ( cape         ( nx , ny ) )
    ALLOCATE ( cin          ( nx , ny ) )
    ALLOCATE ( liftedind    ( nx , ny ) )
    ALLOCATE ( visibility   ( nx , ny ) )
    ALLOCATE ( heatind      ( nx , ny ) )
    ALLOCATE ( lwout        ( nx , ny ) )
    ALLOCATE ( swout        ( nx , ny ) )
    ALLOCATE ( pblhgt       ( nx , ny ) )
    ALLOCATE ( shflux       ( nx , ny ) )
    ALLOCATE ( lhflux       ( nx , ny ) )
    ALLOCATE ( ground_t     ( nx , ny ) )
    IF (make_v5d) ALLOCATE ( abs_vort (nx,ny,kprs))
    IF (make_v5d) ALLOCATE ( thick_10_5 (nx,ny)) 
    IF (make_v5d) ALLOCATE ( snowcover (nx,ny))

    PRINT '("PROCESSING TIME ",A24," FILE = ",A," UNIT # = ",I4)', &
             current_time, TRIM(data_file), lun_data
    CALL process_one_lfm(lun_data,current_time)
    IF (mtype .EQ. 'mm5') THEN
      IF (split_output) CLOSE (lun_data)
    ELSEIF (mtype .EQ. 'wrf') THEN
      IF (split_output) CALL close_wrfnc(lun_data)
    ENDIF 

    ! If smoothing is desired, call the smooth routine.  Two pass
    ! for 2dx filter on surface fields, 20 pass for 4dx on UA
    
    IF (do_smoothing) THEN
      print *, 'Performing extra smoothing on selected variables...'
      ! Extra smoothing for upper-air heights, etc.
      DO pass = 1,10
        CALL SMOOTH(pmsl,nx,ny,1,1.)
        CALL SMOOTH(redp,nx,ny,1,1.)
        CALL SMOOTH(zprs,nx,ny,kprs,1.)
        CALL SMOOTH(srhel,nx,ny,1,1.) 
        IF (make_v5d) THEN
           CALL SMOOTH(abs_vort,nx,ny,kprs,1.)
           CALL SMOOTH(thick_10_5,nx,ny,1,1.)
        ENDIF
      ENDDO
    ENDIF
    ! Call output routines.  This is where you can insert a call
    ! to a routine to write whatever format you want

  
    IF (make_laps) THEN 
      PRINT *, 'Making LAPS output files (fua/fsf)...'
      ! For LAPS, some of the variables need to be manipulated a bit

      WHERE (refl_prs .LE. 0.) refl_prs = -10.
      WHERE (refl_sfc .LE. 0.) refl_sfc = -10.
      WHERE (max_refl .LE. 0.) max_refl = -10. 

      ! Convert intliqwater and totpcpwater to meters (MKS)
      intliqwater = intliqwater*0.001
      totpcpwater = totpcpwater*0.001

      ! Convert cloud species from mixing ratios to mass per volume. Since
      ! they are mixing ratios, we need to simply multiply by the density
      ! of *dry* air
      ALLOCATE (rhodry (nx,ny,kprs))
      DO k = 1,kprs
        rhodry(:,:,k) = prslvl(k) / (r * tprs(:,:,k))
      ENDDO
      ALLOCATE(cldliqcon_prs (nx,ny,kprs))
      ALLOCATE(cldicecon_prs (nx,ny,kprs))
      ALLOCATE(raincon_prs (nx,ny,kprs))
      ALLOCATE(snowcon_prs (nx,ny,kprs))
      ALLOCATE(graupelcon_prs (nx,ny,kprs))
      cldliqcon_prs = cldliqmr_prs * rhodry
      cldicecon_prs = cldicemr_prs * rhodry
      raincon_prs = rainmr_prs * rhodry
      snowcon_prs = snowmr_prs * rhodry
      graupelcon_prs = graupelmr_prs * rhodry
      DEALLOCATE(rhodry)

      CALL output_laps_format(zprs,uprs,vprs,wprs,omprs,tprs, shprs,rhprs, &
                            cldliqcon_prs,cldicecon_prs,raincon_prs, &
                            snowcon_prs,graupelcon_prs, &
                            refl_prs,pcptype_prs,usfc,vsfc,wsfc, &
                            tsfc,tdsfc,rhsfc,cldbase,cldtop,pmsl,redp ,psfc, &
                            intliqwater,totpcpwater,pcp_inc,pcp_tot, &
                            snow_inc,snow_tot,thetasfc,thetaesfc,cape,cin,&
                            cldamt,ceiling,echo_tops,max_refl,refl_sfc, &
                            pcptype_sfc,srhel,liftedind,heatind,visibility,&
                            terdot,lwout,swout,shflux,lhflux,pblhgt,ground_t,&
                            prslvl*0.01,lfmprd_dir,laps_data_root,domain_num,&
                            laps_reftime,laps_valtime,nx,ny,kprs,realtime)

      DEALLOCATE(cldliqcon_prs)
      DEALLOCATE(cldicecon_prs)
      DEALLOCATE(raincon_prs)
      DEALLOCATE(snowcon_prs)
      DEALLOCATE(graupelcon_prs)
    ENDIF

    ! Make GRIB output as required

    IF ( (gribsfc).OR.(gribua)) THEN
      IF (t .EQ. 1) THEN
        CALL make_igds(proj,igds)
      ENDIF
      IF (laps_valtime_prev .GT. 0) THEN
        period_sec = laps_valtime - laps_valtime_prev
      ELSE
        period_sec = 0
      ENDIF
      CALL make_fnam_lp(laps_reftime,a9time,istatus)
      CALL make_fcst_time(laps_valtime,laps_reftime,fcst_hhmm,istatus) 
      gribfile = gribdir(1:gdir_len) // '/' //a9time//'00'//fcst_hhmm(1:2) &
                           //'.grib'
      CALL s_len(gribfile, gfile_len)
      print *, 'Opening ', gribfile(1:gfile_len)
      CALL open_grib_c(gribfile,funit)
      PRINT *, 'OPENED FUNIT ', funit, ' for ',TRIM(gribfile)
      startb = 1
      nbytes = 0
      IF (gribsfc) THEN
        CALL grib_sfc_vars(table_version,center_id,subcenter_id, &
              process_id, laps_reftime,laps_valtime,period_sec,igds,nx,ny, &
              tsfc,tdsfc,rhsfc,usfc,vsfc, &
              wsfc,pmsl,psfc,totpcpwater,pcp_inc,pcp_tot,snow_inc, &
              snow_tot,thetasfc,thetaesfc,cape,cin,srhel, &
              liftedind,terdot,lwout,swout,shflux,lhflux,pblhgt,ground_t, &
              clwmrsfc,icemrsfc,rainmrsfc,snowmrsfc,graupmrsfc, &
              cldamt,cldbase,cldtop,visibility, &
              ceiling,echo_tops,max_refl,refl_sfc, &
              pcptype_sfc, funit,startb,nbytes)
        startb = startb + nbytes
      ENDIF
      IF (gribua) THEN
  
        CALL grib_ua_vars(table_version,center_id,subcenter_id, &
             process_id,laps_reftime,laps_valtime,period_sec,igds,nx,ny,kprs, &
             prslvl*0.01,zprs,uprs,vprs,wprs,omprs,tprs,shprs,rhprs,&
             cldliqmr_prs,cldicemr_prs,rainmr_prs,snowmr_prs,graupelmr_prs, &
             pcptype_prs,refl_prs,funit,startb,nbytes)

        startb = startb+nbytes
      ENDIF
      print *, 'Number of bytes written in GRIB: ',nbytes
 
      CALL close_grib_c(funit)
      gribdone = gribfile(1:gfile_len) // '.done'
      OPEN(99,FILE=gribdone,STATUS='UNKNOWN')
      CLOSE(99)
    ENDIF
    ! Make Vis5D output if requested in namelist

    IF (make_v5d) THEN 

      PRINT *, 'Making Vis5D output...'
      ! If this is the first output time, then the Vis5D file
      ! needs to be initialized
      IF (t .EQ. 1) THEN
                        
         PRINT '(A)', 'Initializing Vis5D Data File...'
         CALL v5dinit(v5d_compress) 

      ENDIF
      WHERE(refl_prs .LT. 0) refl_prs = 0.
      WHERE(refl_sfc .LT. 0) refl_sfc = 0.
      WHERE(max_refl .LT. 0) max_refl = 0.
      ! Call v5dout to fill data for this time period
      IF (time_index_inc .NE. 1) THEN
        t_v5d = t / time_index_inc + 1
      ELSE
        t_v5d = t
      ENDIF
      PRINT '(A,I3)', 'Populating Vis5D file for time period ', t_v5d
      CALL v5dout(t_v5d,zprs,tprs,tdprs,rhprs,uprs,vprs,wprs,omprs, &
                  abs_vort,shprs,cldliqmr_prs,cldicemr_prs,rainmr_prs, &
                  snowmr_prs,graupelmr_prs,refl_prs,pcptype_prs, &
                  thick_10_5,tsfc,tdsfc,rhsfc,usfc,vsfc,wsfc, &
                  pmsl,psfc,cldbase,cldtop,cldamt,ceiling, &
                  intliqwater,totpcpwater,pcp_inc,pcp_tot, &
                  con_pcp_inc, con_pcp_tot, &
                  snow_inc,snow_tot,pcptype_sfc,thetasfc,thetaesfc, &
                  cape,cin,liftedind,srhel,refl_sfc,max_refl,echo_tops, &
                  heatind,visibility,snowcover,lwout,swout,shflux,lhflux, &
                  pblhgt, ground_t)
      IF (t .EQ. num_times_to_proc) CALL v5dend
   
    ENDIF

    IF (make_points) THEN
      CALL split_date_char(times_to_proc(t),year,month,day,hour,minute,second)
      WRITE(date_str, '(I2.2,"/",I2.2,"/",I4.4,1x,I2.2,":",I2.2)') month, &
        day, year,hour,minute
      point_loop: DO pt = 1, num_points
        ip = NINT(point_rec(pt)%i)
        jp = NINT(point_rec(pt)%j)
        t_pt = bint(point_rec(pt)%i, point_rec(pt)%j,tsfc,nx,ny)
        t_pt = (t_pt-273.15)*9./5.+32. ! Convert to F
        td_pt = bint(point_rec(pt)%i,point_rec(pt)%j,tdsfc,nx,ny)
        td_pt =(td_pt-273.15)*9./5.+32.
        rh_pt = NINT(bint(point_rec(pt)%i,point_rec(pt)%j,rhsfc,nx,ny))
        u_pt_g = bint(point_rec(pt)%i,point_rec(pt)%j,usfc,nx,ny)
        v_pt_g = bint(point_rec(pt)%i,point_rec(pt)%j,vsfc,nx,ny)
        CALL gridwind_to_truewind(londot(ip,jp),proj,u_pt_g,v_pt_g, &
                                  u_pt, v_pt)
        CALL uv_to_disp(u_pt, v_pt, dir_pt_r, spd_pt_r)
        dir_pt = NINT(dir_pt_r/10.)*10  ! Integer to nearest 10 degrees
        spd_pt = NINT(spd_pt_r * 1.9425) ! Convert to knots
        vis_pt = bint(point_rec(pt)%i,point_rec(pt)%j,visibility,nx,ny) &
                 *0.00062317
        cld_pt = bint(point_rec(pt)%i,point_rec(pt)%j,cldamt,nx,ny)
        ceiling_pt_r = bint(point_rec(pt)%i,point_rec(pt)%j,ceiling,nx,ny)
        IF (ceiling_pt_r.LT.50000.) THEN
          ceiling_pt = NINT(ceiling_pt_r*3.2808/100.)  ! Hundreds of feet
        ELSE
          ceiling_pt = 999
        ENDIF
        SELECT CASE (INT(pcptype_sfc(ip,jp)))
          CASE (0)
            IF (vis_pt .GE. 7.) THEN
              IF (cld_pt .GE. 0.75) THEN
                wx_pt = 'CLOUDY  '
              ELSE IF ( (cld_pt .LT. 0.75).AND.(cld_pt .GE. 0.25) ) THEN
                wx_pt = 'PT CLDY '
              ELSE
                wx_pt = 'CLEAR   '
              ENDIF
            ELSE IF ((vis_pt .LT. 7.).AND.(vis_pt .GE. 3.))THEN
              wx_pt = 'HAZE    '
            ELSE
              wx_pt = 'FOG     '
            ENDIF
              
          CASE (1)
            wx_pt = 'RAIN    '
          CASE (2)
            wx_pt = 'SNOW    '
          CASE (3)
            wx_pt = 'FRZRAIN '
          CASE (4)
            wx_pt = 'SLEET   '
          CASE (5)
            wx_pt = 'HAIL    '
          CASE (6)
            wx_pt = 'DRIZZLE '
          CASE (7)
            wx_pt = 'MIXED   '
          CASE (8)
            wx_pt = 'MIXED   '
          CASE (9)
            wx_pt = 'RAIN/ICE'
          CASE DEFAULT
            wx_pt = 'UNKNOWN '
        END SELECT
        pcp_pt = bint(point_rec(pt)%i,point_rec(pt)%j,pcp_inc,nx,ny)*39.37
        snow_pt =bint(point_rec(pt)%i,point_rec(pt)%j,snow_inc,nx,ny)*39.37

        WRITE(point_rec(pt)%output_unit, &
        '(A,2F7.1,1x,I3,1x,I3.3,"/",I2.2,1x,I3.3,1x,F5.2,1x,A,2F6.2)') &
        date_str,t_pt,td_pt,rh_pt,dir_pt,spd_pt,ceiling_pt,vis_pt,wx_pt, &
        pcp_pt,snow_pt
   
        IF (t_pt .GT. point_rec(pt)%hi_temp)  THEN
           point_rec(pt)%hi_temp = t_pt
           point_rec(pt)%hi_temp_time = date_str
        ENDIF
        IF (t_pt .LT. point_rec(pt)%lo_temp) THEN
           point_rec(pt)%lo_temp = t_pt
           point_rec(pt)%lo_temp_time = date_str 
        ENDIF
        point_rec(pt)%total_pcp = point_rec(pt)%total_pcp + pcp_pt
        point_rec(pt)%total_snow = point_rec(pt)%total_snow + snow_pt
        point_rec(pt)%avg_temp = point_rec(pt)%avg_temp + t_pt
        point_rec(pt)%avg_dewpt = point_rec(pt)%avg_dewpt + td_pt   
        IF (t .EQ. num_times_to_proc) THEN
          ! Finish summary and close the file
         WRITE(point_rec(pt)%output_unit,'(80x)')
         WRITE(point_rec(pt)%output_unit,'(80x)') 
         WRITE(point_rec(pt)%output_unit, '("SUMMARY INFORMATION FOR PERIOD")')
         WRITE(point_rec(pt)%output_unit,&
           '("--------------------------------------------------")')
         WRITE(point_rec(pt)%output_unit, &
          '("HIGH TEMPERATURE:  ",F6.1,2x,"AT ",A)') point_rec(pt)%hi_temp, &
          point_rec(pt)%hi_temp_time
         WRITE(point_rec(pt)%output_unit, &
          '("LOW TEMPERATURE:   ",F6.1,2x,"AT ",A)') point_rec(pt)%lo_temp, &
          point_rec(pt)%lo_temp_time 
         point_rec(pt)%avg_temp = point_rec(pt)%avg_temp / &
                      FLOAT(num_times_to_proc/time_index_inc)
         WRITE(point_rec(pt)%output_unit, &
          '("AVG TEMPERATURE:   ",F6.1)') point_rec(pt)%avg_temp
         point_rec(pt)%avg_dewpt = point_rec(pt)%avg_dewpt / &
                 FLOAT(num_times_to_proc/time_index_inc)
         WRITE(point_rec(pt)%output_unit, &
          '("AVG DEWPOINT:      ",F6.1)') point_rec(pt)%avg_dewpt
         WRITE(point_rec(pt)%output_unit, &
          '("TOTAL PRECIP:      ",F6.2)') point_rec(pt)%total_pcp
         WRITE(point_rec(pt)%output_unit, &
          '("TOTAL SNOW:        ",F6.2)') point_rec(pt)%total_snow
         CLOSE(point_rec(pt)%output_unit)
       ENDIF
      ENDDO point_loop
    ENDIF
    ! Deallocate all variables except pcp/snow init/total
    DEALLOCATE ( psfc )
    DEALLOCATE ( tsfc )
    DEALLOCATE ( thetasfc )
    DEALLOCATE ( thetaesfc )
    DEALLOCATE ( rhsfc )
    DEALLOCATE ( tdsfc )
    DEALLOCATE ( thetaprs )
    DEALLOCATE ( zprs )
    DEALLOCATE ( rhprs )
    DEALLOCATE ( tprs )
    DEALLOCATE ( tdprs )
    DEALLOCATE ( shprs )
    DEALLOCATE ( redp )
    DEALLOCATE ( pmsl )
    DEALLOCATE ( usfc )
    DEALLOCATE ( uprs )
    DEALLOCATE ( vsfc )
    DEALLOCATE ( vprs )
    DEALLOCATE ( wsfc )
    DEALLOCATE ( wprs )
    DEALLOCATE ( omprs )
    DEALLOCATE ( cldbase )
    DEALLOCATE ( cldtop )   
    DEALLOCATE ( cldamt )
    DEALLOCATE ( ceiling )
    DEALLOCATE ( intliqwater )
    DEALLOCATE ( totpcpwater )
    DEALLOCATE ( max_refl )
    DEALLOCATE ( echo_tops )
    DEALLOCATE ( cldliqmr_prs )
    DEALLOCATE ( cldicemr_prs )
    DEALLOCATE ( rainmr_prs )
    DEALLOCATE ( snowmr_prs )
    DEALLOCATE ( graupelmr_prs )
    DEALLOCATE ( refl_prs )
    DEALLOCATE ( refl_sfc )
    DEALLOCATE ( pcptype_sfc )
    DEALLOCATE ( pcptype_prs )
    DEALLOCATE ( pcp_inc )
    DEALLOCATE ( con_pcp_inc)
    DEALLOCATE ( snow_inc )
    DEALLOCATE ( cape )
    DEALLOCATE ( cin )
    DEALLOCATE ( liftedind  )
    DEALLOCATE ( visibility )
    DEALLOCATE ( heatind   )
    DEALLOCATE ( lwout )
    DEALLOCATE ( swout )
    DEALLOCATE ( shflux )
    DEALLOCATE ( lhflux )
    DEALLOCATE ( pblhgt )
    DEALLOCATE ( ground_t )
    DEALLOCATE ( clwmrsfc )
    DEALLOCATE ( icemrsfc )
    DEALLOCATE ( rainmrsfc )
    DEALLOCATE ( snowmrsfc )
    DEALLOCATE ( graupmrsfc )
    IF (make_v5d) THEN
      DEALLOCATE (abs_vort)
      DEALLOCATE (thick_10_5)
      DEALLOCATE (snowcover)
    ENDIF
    IF (realtime) THEN
      IF (split_output) THEN
       ! CLOSE(89,STATUS='DELETE')
      ENDIF
    ENDIF
    laps_valtime_prev = laps_valtime
  ENDDO time_loop
 
  ! Deallocate the precip arrays
  DEALLOCATE ( pcp_init )
  DEALLOCATE ( con_pcp_tot )
  DEALLOCATE ( con_pcp_init )
  DEALLOCATE ( pcp_tot )
  DEALLOCATE ( snow_init )
  DEALLOCATE ( snow_tot )
  DEALLOCATE ( latdot )
  DEALLOCATE ( londot )
  DEALLOCATE ( terdot )
END PROGRAM lfmpost

