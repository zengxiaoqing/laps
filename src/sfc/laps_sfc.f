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
c
c
	program laps_sfc
c
        include 'lapsparms.cmn'
        include 'grid_fname.cmn'

	character laps_domain*9
c
	call get_config(istatus)
	if(istatus .ne. 1) then
	   write(6,*) 'LAPS_SFC: ERROR getting domain dimensions'
	   stop
	endif
c
	laps_domain = grid_fnam_common
	call laps_sfc_sub(nx_l_cmn,ny_l_cmn,nk_laps,maxstns_cmn,
     &                    laps_cycle_time_cmn,grid_spacing_m_cmn,
     &                    laps_domain)
c
	end
c
c
	subroutine laps_sfc_sub(ni,nj,nk,mxstn,laps_cycle_time,
     &                          grid_spacing,laps_domain)
c
c
c*****************************************************************************
c
c	Driver program for the LAPS variational surface analysis.  This 
c	program gets the correct time and passes it to the routines (former
c	programs) that get the mesonet, SAO, and VAS data; grid and quality
c	control it; and perform the variational analysis.
c
c
c	History:
c		P. Stamus	03-09-90  Original version.
c				03-22-90  Check for BATCH or INTERACTIVE.
c				03-29-90  Added check for old VAS data.
c				04-06-90  Pass del,gam,ak to mdat for hdr.
c				04-11-90  Pass lat/lon & topo to routines
c				12-10-90  Change to RT_DEV for lat/lon, topo.
c				11-08-91  Changes for new grids.
c				10-15-92  Add Steve's get_laps_config call.
c				01-06-93  New version...new lso and lvd stuff.
c				04-19-93  Common for fnorm (Barnes routine)
c				08-05-93  Replace fnorm with fnorm2. 
c                               12-08-93  Change to format in 'systime.dat'.
c                               03-04-94  Change static file read (read to rd)
c                               07-20-94  Add include file.
c                               02-03-95  Background reads here...then pass.
c                               07-20-95  Bag fnorm2 for new Barnes wt method.
c                               08-08-95  Add code for verify routine.
c                               03-22-96  Fixes for 30 min cycle.
c                               04-10-96  More 30 min...bkgwts file.
c                               11-07-96  Use 3d sfc wind for bkg, adj wts.
c                               12-13-96  More porting changes...common for
c                                         sfc data, LGS grids. Bag stations.in
c                               03-26-97  Add ability to do interactive runs.
c                                         (Removes need for 'laps_sfci')
c                                         Remove equivs.
c                               09-11-97  Changes for dynamic LAPS.
c                               01-20-98  Move wt calcs to correct place.
c	                        07-11-98  New background calls.  Still temp fix
c                                           until can remove weight code (new 
c                                           spline doesn't use).
c                               09-24-98  Carry background flags for each var.
c                                           Rm ceil QC check.
c                               09-30-98  Housekeeping.
c                               12-02-98  Remove status check for LT1.
c                               07-07-99  General upgrades and cleanup.  Change
c                                           read_surface_obs to read_surface_data.
c                                           Rm *4 from all declarations.
c                               09-19-99  Check T/Td bkgs until LGB can do it.
c                               12-01-99  Rotate bkg winds to grid north.
c                               12-17-99  Add option to use either LSO or the
c                                           Kalman estimate LSO_QC.
c                               01-23-00  Make that option readable via namelist.
c                               01-26-00  Skip qcdata call if using LSO_QC.
c
c
c       Notes:
c
c       1. Running 'laps_sfc.x' makes the analysis use the time from 
c          '../sched/systime.dat'.  Running 'laps_sfc.x -i' allows the
c          user to enter the desired start time.
c
c       2. The background weights for the station pressure (wt_stnp) are
c          not currently used.  They are included because they might be in
c          the future.
c
c*****************************************************************************
c
	real lat(ni,nj), lon(ni,nj), topo(ni,nj), ldf(ni,nj)
	real x1a(ni), x2a(nj), y2a(ni,nj)
	real grid_spacing
c
	integer*4 i4time
	integer jstatus(20)		! 20 is standard for prodgen drivers
	integer narg, iargc
c
	character atime*24, filename*9, filename_last*9
	character infile_last*256
	character dir_s*256,ext_s*31,units*10,comment*125,var_s*3
	character laps_domain*9, use*6
c
c.....  Stuff for backgrounds.
c
	real u_bk(ni,nj), v_bk(ni,nj), t_bk(ni,nj), td_bk(ni,nj) ! Kt & Deg F
	real wt_u(ni,nj), wt_v(ni,nj), wt_t(ni,nj), wt_td(ni,nj)
	real rp_bk(ni,nj), mslp_bk(ni,nj), stnp_bk(ni,nj)
	real wt_rp(ni,nj), wt_mslp(ni,nj), wt_stnp(ni,nj)
	real vis_bk(ni,nj), wt_vis(ni,nj)
	real tgd_bk_f(ni,nj), wt_tgd(ni,nj)
	real wt(ni,nj)
        real dum_2d(ni,nj)
c
	integer back_t, back_td, back_rp, back_uv, back_vis, back_sp
	integer back_mp, back_tgd
	character var_req*4, ext_bk*31, back*9
c
c..... Stuff for the sfc data and other station info (LSO +)
c
        include 'sfcob.inc'
        type (sfcob) obs(mxstn)

	real lat_s(mxstn), lon_s(mxstn), elev_s(mxstn)
	real t_s(mxstn), t_ea(mxstn), max24t(mxstn), min24t(mxstn)
        real td_s(mxstn), td_ea(mxstn), rh(mxstn), rh_ea(mxstn)

        real dd_s(mxstn), ddg_s(mxstn), dd_ea(mxstn)
        real ff_s(mxstn), ffg_s(mxstn), ff_ea(mxstn)

        real alt_s(mxstn), alt_ea(mxstn), delp(mxstn)
	real pstn_s(mxstn), pmsl_s(mxstn), p_ea(mxstn)

	real store_hgt(mxstn,5) 

        real vis_s(mxstn), vis_ea(mxstn)
        real solar_s(mxstn), solar_ea(mxstn)

        real sfct(mxstn), sfct_ea(mxstn)
        real sfcm(mxstn), sfcm_ea(mxstn)
        real pcp1(mxstn), pcp3(mxstn), pcp6(mxstn), pcp24(mxstn)
        real snow(mxstn), snow_ea(mxstn), pcp_ea(mxstn)

	real rii(mxstn), rjj(mxstn)
c
	integer kloud_s(mxstn), obstime(mxstn)
        integer wmoid(mxstn), delpch(mxstn)
	integer ii(mxstn), jj(mxstn)
c
	character atime_s*24
	character store_amt(mxstn,5)*4
        character stations(mxstn)*20, provider(mxstn)*11,stn20*20      
        character reptype(mxstn)*6, autostntype(mxstn)*6
        character wx_s(mxstn)*25 
c
c.....  Work arrays for the QC routine.
c
        integer rely(26,mxstn)
	character stn3(mxstn)*3
c
c.....  Stuff for intermediate grids (old LGS file)
c
	real u1(ni,nj), v1(ni,nj)
	real t1(ni,nj), td1(ni,nj), tb81(ni,nj)
	real rp1(ni,nj), sp1(ni,nj), mslp1(ni,nj)
	real vis1(ni,nj), elev1(ni,nj)
c
c.....  Namelist stuff
c
	integer use_lso_qc, skip_internal_qc, itheta
	character nl_file*256
c
	namelist /surface_analysis/ use_lso_qc,skip_internal_qc
     1                             ,itheta, redp_lvl, del, gam, ak       
c
c*************************************************************
c.....	Start here.  First see if this is an interactive run.
c*************************************************************
c
	call tagit('laps_sfc',20000126)
        call get_sfc_badflag(badflag,istatus)
	narg = iargc()
cc	print *,' narg = ', narg
c
c.....  Now get the analysis time from the scheduler or the user.
c
	if(narg .eq. 0) then
c
	   ihours = 1	! default # of hrs back for time-tendencies
c
	   call get_systime(i4time,filename,istatus)
c
	else
c
 970	   write(6,973)
 973	   format(' Enter input filename (yydddhhmm): ',$)
	   read(5,972) filename
 972	   format(a)
c
 974	   write(6,975)
 975	   format(
     &    ' Hours to go back for time-tendencies (1..6)[1]: ',$)
	   read(5,976) ihours
 976	   format(i1)
	   if(ihours .eq. 0) then
	      print *, ' Will skip backgrounds.'
	   elseif(ihours.lt.1 .or. ihours.gt.6) then
	      print *, ' ERROR: Hrs out of bounds.  Try again.'
	      go to 974
	   endif
	   call i4time_fname_lp(filename,i4time,status)
c
	endif
c
c
500	dt = ihours * laps_cycle_time    ! multiples of time to get bkgs.
	i4time_last = i4time - ifix( dt )
	call cv_i4tim_asc_lp(i4time, atime, status)	   ! get the atime
	call make_fnam_lp(i4time_last,filename_last,istatus)  ! make earlier filename	
c
        del = 1.3e5
	gam = .000004
	ak = 1.e-6
c
c.....  Set the namelist variables to their defaults.  If there's a problem reading
c.....  the namelist, at least we can continue.
c
	use_lso_qc = 0        !use normal LSO
	skip_internal_qc = 0  !use internal QC routine
c
c.....  Read the namelist and get that info, then get the LAPS lat/lon and topo 
c.....  data so we can pass them to the routines that need them. 
c
cc	dir_s = '../static/' 
	call get_directory('static', dir_s, len)
	ext_s = laps_domain
	nl_file = dir_s(1:len) // 'surface_analysis.nl'
	call s_len(nl_file, len)
c
	open(20,file=nl_file(1:len),status='old',err=930)
	read(20,surface_analysis,end=930)
	close(20)
	go to 501
c
c.....  Skip here if there are problems with the namelist.
c
 930	continue
	print *,' '
	print *,' WARNING.  Problem reading the surface analysis ',
     &          'namelist file.'
	print *,'    Check to see if ', nl_file(1:len)
        print *,'         is there and correct.'
	print *,'    Continuing with default values for namelist ',
     &          'variables.'
	print *,' '
c
c.....  Continue...get static stuff.
c
 501	continue
	var_s = 'LAT'
        call rd_laps_static(dir_s,ext_s,ni,nj,1,var_s,units,comment,
     &                      lat ,grid_spacing,istatus)
	var_s = 'LON'
        call rd_laps_static(dir_s,ext_s,ni,nj,1,var_s,units,comment,
     &                      lon ,grid_spacing,istatus)
	var_s = 'AVG'
        call rd_laps_static(dir_s,ext_s,ni,nj,1,var_s,units,comment,
     &                      topo ,grid_spacing,istatus)
	var_s = 'LDF'
        call rd_laps_static(dir_s,ext_s,ni,nj,1,var_s,units,comment,
     &                      ldf  ,grid_spacing,istatus)
c
c.....  Read in the obs and calculate a weight based on distance to each
c.....  station.
c
c.....	READ IN THE SURFACE OBS:  dd/ff in deg/kt, t and td in F, elev in m,
c.....	                          and the pressure variable. cld hts are msl.
c
c
	if(use_lso_qc .eq. 1) then
	   use = 'LSO_QC'
	else
	   use = 'LSO   '
	endif
c
	write(6,305) filename(1:9), use
 305	format(' Getting surface data at: ',a9,' from ',a6)
c
	if(use_lso_qc .ne. 1) then
	   call read_surface_data(i4time,atime_s,n_obs_g,n_obs_b,      !regular LSO
     &       obstime,wmoid,stations,provider,wx_s,reptype,autostntype,
     &       lat_s,lon_s,elev_s,t_s,td_s,rh,dd_s,ff_s,ddg_s,ffg_s,
     &       alt_s,pstn_s,pmsl_s,delpch,delp,vis_s,solar_s,sfct,sfcm,
     &       pcp1,pcp3,pcp6,pcp24,snow,kloud_s,max24t,min24t,t_ea,
     &       td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,
     &       sfct_ea,sfcm_ea,pcp_ea,snow_ea,store_amt,store_hgt,mxstn,
     &       istatus)
	else
	   call read_surface_dataqc(i4time,atime_s,n_obs_g,n_obs_b,    !QC'd LSO
     &       obstime,wmoid,stations,provider,wx_s,reptype,autostntype,
     &       lat_s,lon_s,elev_s,t_s,td_s,rh,dd_s,ff_s,ddg_s,ffg_s,
     &       alt_s,pstn_s,pmsl_s,delpch,delp,vis_s,solar_s,sfct,sfcm,
     &       pcp1,pcp3,pcp6,pcp24,snow,kloud_s,max24t,min24t,t_ea,
     &       td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,
     &       sfct_ea,sfcm_ea,pcp_ea,snow_ea,store_amt,store_hgt,mxstn,
     &       istatus)
	endif
c
        if(istatus .ne. 1 .and. istatus .ne. -1)then
            write(6,*)' bad istatus from read_surface_data..',istatus       
            stop
        endif

	do i=1,n_obs_b ! Put obs into data structure 
            obs(i)%sfct_f = sfct(i)
            obs(i)%t_ea_f = t_ea(i)
        enddo ! i

	print *,' '
	write(6,320) use,atime_s,n_obs_g,n_obs_b
320	format(1x,a6,' data valid time: ',a24,' Num obs: ',2i6)
c
	print *,' '
c
c.....  Copy 3 characters of the station name to another array for use
c.....  in the qc_data routine.  This can be removed once the QC is changed
c.....  to the stand-alone Kalman.  Blank out the array first.
c
	do i=1,mxstn
	   stn3(i)(1:3) = '   '
	enddo !i
c
	do i=1,n_obs_b ! Rearrange station strings
           stn20 = stations(i)
           call right_justify(stn20)
	   stn3(i) = stn20(18:20)
           call left_justify(stations(i))
	enddo !i
c
c.....	Find the i,j location of each station, then calculate the
c.....  background weights (based on station density).
c
	call find_ij(lat_s,lon_s,lat,lon,n_obs_b,mxstn,
     &               ni,nj,ii,jj,rii,rjj)
c
        call get_r_missing_data(r_missing_data,istatus)

        do ista=1,n_obs_b ! Write station info including elevation
                          ! Insert info into obs data structure
           obs(ista)%ri = rii(ista)
           obs(ista)%rj = rjj(ista)
           obs(ista)%i = ii(ista)
           obs(ista)%j = jj(ista)

           if(       ii(ista) .ge. 1 .and. ii(ista) .le. ni
     1         .and. jj(ista) .ge. 1 .and. jj(ista) .le. nj )then ! in domain

               elev_diff = elev_s(ista) - topo(ii(ista),jj(ista))
               write(6,999)ista,stations(ista)(1:5), reptype(ista)(1:6)
     &                    ,autostntype(ista)(1:6), rii(ista), rjj(ista)
     &                    ,elev_s(ista),topo(ii(ista),jj(ista))
     &                    ,elev_diff

               obs(ista)%elev_diff = elev_diff

           else
               write(6,999)ista,stations(ista)(1:5), reptype(ista)(1:6)
     &                    ,autostntype(ista)(1:6), rii(ista), rjj(ista)

               obs(ista)%elev_diff = r_missing_data

           endif

        enddo !ista
 999    format(i5,': ',a5,2x,a6,2x,a6,' is at i,j: ',f5.1,',',f5.1
     1        ,2x,3f8.0)
c
        call zero(wt, ni,nj)
        if(n_obs_b .gt. 0.)then
	    call bkgwts(lat,lon,topo,n_obs_b,lat_s,lon_s,elev_s,
     &                  rii,rjj,wt,ni,nj,mxstn)
        endif
c
c.....  Zero the weights then get the background data.  Try for an FSF
c.....  forecast first, except for surface winds from the 3d analysis, then
c.....  a previous LAPS analysis. If both LAPS and FSF are missing, print a
c.....  warning. Ground temp is read from 'get_modelfg_2d', that looks for
c.....  an FSF forecast according to the fdda parameter, otherwise LGB.
c
        call zero(wt_u, ni,nj)
        call zero(wt_v, ni,nj)
        call zero(wt_t, ni,nj)
        call zero(wt_td, ni,nj)
        call zero(wt_rp, ni,nj)
        call zero(wt_mslp, ni,nj)
        call zero(wt_stnp, ni,nj)
        call zero(wt_vis, ni,nj)
        call zero(wt_tgd, ni,nj)
c
	back_t = 0
	back_td = 0
	back_rp = 0
	back_sp = 0
	back_mp = 0
	back_uv = 0
	back_vis = 0
	back_tgd = 0
c
	if(ihours .eq. 0) then     ! skip the backgrounds
	   ilaps_bk = 0
	   print *,' Skipping backgrounds...'
	   go to 600
	endif
c
c.....  Get the backgrounds.  Convert units while we're here.
c
	call get_bkgwind_sfc(i4time,ext_bk,ibkg_time,u_bk,v_bk,
     &            laps_cycle_time,ni,nj,istatus)
	if(istatus .eq. 1) then
	   call make_fnam_lp(ibkg_time,back,istatus)
	   write(6,951) ext_bk(1:6), back
	   call move(wt, wt_u, ni,nj)
	   call move(wt, wt_v, ni,nj)
c
c.....  Rotate background winds from true north to grid north.
c
	   do j=1,nj
	   do i=1,ni
	      utrue = u_bk(i,j)
	      vtrue = v_bk(i,j)
	      call uvtrue_to_uvgrid(utrue, vtrue, ugrid, vgrid,lon(i,j))
	      u_bk(i,j) = ugrid
	      v_bk(i,j) = vgrid
	   enddo !i
	   enddo !j
c
	   call conv_ms2kt(u_bk,u_bk,ni,nj)
	   call conv_ms2kt(v_bk,v_bk,ni,nj)
	   back_uv = 1
	else
	   print *,'     No background available'
	   call zero(u_bk,ni,nj)
	   call zero(v_bk,ni,nj)
	endif
c
	print *,' '
	print *,' Getting temperature background....'
	var_req = 'TEMP'
	call get_background_sfc(i4time,var_req,ext_bk,ibkg_time,dum_2d,
     &          laps_cycle_time,ni,nj,istatus)
	if(istatus .eq. 1) then
	   ilaps_bk = 1
	   call make_fnam_lp(ibkg_time,back,istatus)
!          write(6,951) ext_bk(1:6), back
	   call conv_k2f(dum_2d,t_bk,ni,nj) ! conv K to deg F
	   call move(wt, wt_t, ni,nj)
	   back_t = 1
	else
	   print *,'     No background available for ',var_req
	   call zero(t_bk,ni,nj)
	endif
c
	print *,' '
	print *,' Getting MSL pressure background....'
	var_req = 'MSLP'
	call get_background_sfc(i4time,var_req,ext_bk,ibkg_time,mslp_bk,
     &          laps_cycle_time,ni,nj,istatus)
	if(istatus .eq. 1) then
	   call make_fnam_lp(ibkg_time,back,istatus)
!	   write(6,951) ext_bk(1:6), back
	   call multcon(mslp_bk,0.01,ni,nj) ! conv Pa to mb
	   call move(wt, wt_mslp, ni,nj)
	   back_mp = 1
	else
	   print *,'     No background available for ',var_req
	   call zero(mslp_bk,ni,nj)
	endif
c
	print *,' '
	print *,' Getting station pressure background....'
	var_req = 'SFCP'
	call get_background_sfc(i4time,var_req,ext_bk,ibkg_time,stnp_bk,
     &          laps_cycle_time,ni,nj,istatus)
	if(istatus .eq. 1) then
	   call make_fnam_lp(ibkg_time,back,istatus)
!          write(6,951) ext_bk(1:6), back
	   call move(wt, wt_stnp, ni,nj)
	   call multcon(stnp_bk,0.01,ni,nj) ! conv Pa to mb
	   back_sp = 1
	else
	   print *,'     No background available for ',var_req
	   call zero(stnp_bk,ni,nj)
	endif
c
	print *,' '
	print *,' Getting visibility background....'
	var_req = 'VISB'
	call get_background_sfc(i4time,var_req,ext_bk,ibkg_time,vis_bk,
     &          laps_cycle_time,ni,nj,istatus)
	if(istatus .eq. 1) then
	   call make_fnam_lp(ibkg_time,back,istatus)
!	   write(6,951) ext_bk(1:6), back
	   call move(wt, wt_vis, ni,nj)
	   call conv_m2miles(vis_bk,vis_bk,ni,nj)
	   call visg2log(vis_bk,ni,nj,badflag) ! conv miles to log(miles)
	   back_vis = 1
	else
	   print *,'     No background available for ',var_req
	   call zero(vis_bk,ni,nj)
	endif
c
	print *,' '
	print *,' Getting dew point temperature background....'
	var_req = 'DEWP'
	call get_background_sfc(i4time,var_req,ext_bk,ibkg_time,dum_2d,
     &          laps_cycle_time,ni,nj,istatus)
	if(istatus .eq. 1) then
	   call make_fnam_lp(ibkg_time,back,istatus)
!	   write(6,951) ext_bk(1:6), back
	   call move(wt, wt_td, ni,nj)
	   call conv_k2f(dum_2d,td_bk,ni,nj) ! conv K to deg F
	   back_td = 1
	else
	   print *,'     No background available for ',var_req
	   call zero(td_bk,ni,nj)
	endif
	print *,' '
c
	print *,' '
	print *,' Getting reduced pressure background....'
	var_req = 'REDP'
	call get_background_sfc(i4time,var_req,ext_bk,ibkg_time,rp_bk,
     &          laps_cycle_time,ni,nj,istatus)
	if(istatus .eq. 1) then
	   call make_fnam_lp(ibkg_time,back,istatus)
!	   write(6,951) ext_bk(1:6), back
	   call move(wt, wt_rp, ni,nj)
	   call multcon(rp_bk,0.01,ni,nj) ! conv Pa to mb
	   back_rp = 1
	else
	   print *,'     No background available for ',var_req
	   call zero(rp_bk, ni,nj)
	endif
c
	print *,' '
	print *,' Getting ground temperature background....'
	var_req = 'tgd'
        call get_modelfg_2d(i4time,var_req,ni,nj,dum_2d,istatus)       
!       istatus = 0 ! Temporary statement until 'get_modelfg_2d' is fixed
	if(istatus .eq. 1) then
	   call make_fnam_lp(i4time,back,istatus)
	   write(6,*) 'Read successful for ',var_req
	   call move(wt, wt_tgd, ni,nj)
           call conv_k2f(dum_2d,tgd_bk_f,ni,nj) ! conv K to deg F
	   back_tgd = 1
	elseif(back_t .eq. 1 .and. back_td .eq. 1)then
           write(6,*)' Using the temp & dewpoint bkg as a substitute'    
           do j = 1,nj
           do i = 1,ni
               tgd_bk_f(i,j) = t_bk(i,j)  *       ldf(i,j) 
     1                       + td_bk(i,j) * (1. - ldf(i,j))       
           enddo ! i
           enddo ! j
	else
	   print *,'     No background available for ',var_req
	   call zero(tgd_bk_f, ni,nj)
	endif
c
 951	format(3x,'Using background from ',a6,' at ',a9)
c
c.....  Adjust wts for winds; if above 2500 m, increase wt by order of mag
c.....  so background winds will have more influence.
c
	do j=1,nj
	do i=1,ni
	   if(topo(i,j) .gt. 2500.) then
	      wt_u(i,j) = wt_u(i,j) * 10.
	      wt_v(i,j) = wt_v(i,j) * 10.
	   endif
	enddo !i
	enddo !j
c
 600	continue
c
c.....	QC the surface data.
c
	if(use_lso_qc .eq. 1) then        !data already qc'd
	   print *,' ** Using pre-QCd LSO_QC...skipping LSX QC step.'       
	   go to 521
	endif
	if(skip_internal_qc .eq. 1) then  !skip QC routine
	   print *, ' **  Skipping LSX internal QC routine **  '
	   go to 521
	endif
c
	call get_directory('lso', infile_last, len)
	infile_last = infile_last(1:len) // filename_last // '.lso'
c
	call qcdata(filename,infile_last,rely,mxstn,
     &     t_s, td_s, dd_s, ff_s, ddg_s, ffg_s, pstn_s, pmsl_s, alt_s, 
     &     vis_s, stn3, rii, rjj, ii, jj, n_obs_b, n_sao_b, n_sao_g,
     &     ni,nj,mslp_bk,back_mp,
     &     istatus)
c
	if(istatus .eq. 1) then
	   jstatus(2) = 1
	elseif(istatus .eq. 0) then
	   jstatus(2) = 0
	   print *, ' +++ No data for QC routine. +++'
	   go to 521
	else
	   print *,
     &        ' +++ ERROR.  Problem in QC routine. +++'
	   jstatus(2) = -2
	   go to 521
	endif
c
c.....	Check each of the primary analysis variables.
c
	do mm=1,n_obs_b
	  if(rely(7,mm) .lt. 0) then	! temperature
	     print *, 'QC: Bad T at ',stations(mm),' with rely/value '       
     1              ,rely(7,mm),t_s(mm)
	     t_s(mm) = badflag
	  endif
 121	enddo  !mm

	do mm=1,n_obs_b
	  if(rely(8,mm) .lt. 0) then	! dewpt
	      print *, 'QC: Bad TD at ',stations(mm)
     1               ,' with rely/value ',rely(8,mm),td_s(mm)
	      td_s(mm) = badflag
	  endif
 122	enddo  !mm

	do mm=1,n_obs_b
	  if(rely(9,mm) .lt. 0) then	! wind direction
	   print *, 'QC: Bad DIR at ',stations(mm),' with rely/value '       
     1            ,rely(9,mm),dd_s(mm)
	    dd_s(mm) = badflag
	  endif
 123	enddo  !mm

	do mm=1,n_obs_b
	  if(rely(10,mm) .lt. 0) then	! wind speed
	      print *, 'QC: Bad SPD at ',stations(mm)
     1               ,' with rely/value ',rely(10,mm),ff_s(mm)
	      ff_s(mm) = badflag
	  endif
 124	enddo  !mm

	do mm=1,n_obs_b
          if(rely(14,mm) .lt. 0) then	! mslp
	      print *, 'QC: Bad MSLP at ',stations(mm)
     1               ,' with rely/value ',rely(14,mm),pmsl_s(mm)
	      pmsl_s(mm) = badflag
	  endif
 125	enddo  !mm

	do mm=1,n_obs_b
	  if(rely(15,mm) .lt. 0) then	! altimeter 
	      print *, 'QC: Bad ALT at ',stations(mm)
     1               ,' with rely/value ',rely(15,mm),alt_s(mm)
	      alt_s(mm) = badflag
	  endif
 126	enddo  !mm

	do mm=1,n_obs_b
          if(rely(25,mm) .lt. 0) then	! visibility 
	      print *, 'QC: Bad VIS at ',stations(mm)
     1               ,' with rely/value ',rely(25,mm),vis_s(mm)
	      vis_s(mm) = badflag
	  endif
 128	enddo  !mm

 521	continue                          
c
c.....  QC the backgrounds.  If Td > T, set Td = T...temporary fix until LGB
c.....  can do this....
c
	do j=1,nj
	do i=1,ni
	   if(td_bk(i,j) .gt. t_bk(i,j)) td_bk(i,j) = t_bk(i,j)
	enddo !i
	enddo !j
c
c.....  Set up arrays for the verify routine.
c
	do i=1,ni
	   x1a(i) = float(i)
	enddo !i
	do j=1,nj
	   x2a(j) = float(j)
	enddo !j
c
c
cx....	Now call MDATLAPS to put the data on the grid.
c
	call mdat_laps(i4time,atime,ni,nj,mxstn,laps_cycle_time,lat,
     &     lon,topo,x1a,x2a,y2a,redp_lvl,
     &     lon_s, elev_s, t_s, td_s, dd_s, ff_s, pstn_s, pmsl_s, 
     &     alt_s, vis_s, stations, rii, rjj, ii, jj, n_obs_b, n_sao_g,
     &     u_bk, v_bk, t_bk, td_bk, rp_bk, mslp_bk, stnp_bk, vis_bk,
     &     wt_u, wt_v, wt_rp, wt_mslp, ilaps_bk, 
     &     u1, v1, rp1, t1, td1, sp1, tb81, mslp1, vis1, elev1,
     &     back_t,back_td,back_uv,back_sp,back_rp,back_mp,back_vis,
     &     jstatus) 
c
	if(jstatus(1) .ne. 1) then
	   print *,' From MDAT_LAPS:  Error Return.  Stop.'
	   stop 
	endif
c
c
c.....	Call LAPSVANL to do the actual variational analysis, and calculate
c.....	derived variables, etc.  The output file goes to the lapsprd 
c.....	directory (machine dependent) and has the extension '.lsx'.
c
c       Units at this stage are: p  in mb
c                                t  in F
c                                td in F
c                                u  in kt
c                                v  in kt
	call laps_vanl(i4time,filename,ni,nj,nk,mxstn,
     &     itheta,redp_lvl,laps_cycle_time,
     &     dt,del,gam,ak,lat,lon,topo,ldf,grid_spacing, laps_domain,
     &     lat_s, lon_s, elev_s, t_s, td_s, ff_s, pstn_s, pmsl_s,
     &     vis_s, stations, n_obs_b, n_sao_b, n_sao_g, obs,
     &     u_bk,v_bk,t_bk,td_bk,rp_bk,mslp_bk,stnp_bk,vis_bk,tgd_bk_f,   
     &     wt_u, wt_v, wt_t, wt_td, wt_rp, wt_mslp, wt_vis, ilaps_bk, 
     &     back_t,back_td,back_uv,back_sp,back_rp,back_mp,back_vis,
     &     u1, v1, rp1, t1, td1, sp1, tb81, mslp1, vis1, elev1,
     &     x1a,x2a,y2a,ii,jj,jstatus)
c
	if(jstatus(3) .ne. 1) then
	  print *,' From LAPS_VANL: Error Return.' 
	endif
c
c.....	That's about it...let's go home.
c
	print *,' End of LAPS Surface Analysis'
	return
c
	end

