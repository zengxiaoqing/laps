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
	character laps_domain*9
c
	laps_domain = 'nest7grid'
	call get_laps_config(laps_domain,istatus)
	if(istatus .ne. 1) then
	   write(6,*) 'LAPS_SFC: ERROR getting domain dimensions'
	   stop
	endif
c
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
	include 'laps_sfc.inc'
c
	real*4 lat(ni,nj), lon(ni,nj), topo(ni,nj)
	real*4 x1a(ni), x2a(nj), y2a(ni,nj)
	real*4 grid_spacing
c
	integer*4 jstatus(20)	! 20 is standard for prodgen drivers
	integer narg, iargc
c
	character atime*24, filename*9, filename_last*9
	character infile1*256, infile_last*256
	character dir_s*256,ext_s*31,units*10,comment*125,var_s*3
	character laps_domain*9
c
c.....  Stuff for backgrounds.
c
	real*4 u_bk(ni,nj), v_bk(ni,nj), t_bk(ni,nj), td_bk(ni,nj)
	real*4 wt_u(ni,nj), wt_v(ni,nj), wt_t(ni,nj), wt_td(ni,nj)
	real*4 rp_bk(ni,nj), mslp_bk(ni,nj), stnp_bk(ni,nj)
	real*4 wt_rp(ni,nj), wt_mslp(ni,nj), wt_stnp(ni,nj)
	real*4 vis_bk(ni,nj), wt_vis(ni,nj)
	real*4 wt(ni,nj)  
c
	real*4 background(ni,nj,8)
	integer*4 lvl_bk(8)
	character var_bk(8)*3, back*9, unitsl(8)*10, lvlcl(8)*4
	character coml(8)*125, dir_bk*256, ext_bk*31
c
	real*4 bk_rams(ni,nj,5)
	integer*4 lvl_bkr(5)
	character var_bkr(5)*3, backrams*9, dir_rams*256, ext_rams*31
	character unitsr(5)*10, lvlcr(5)*4, comr(5)*125
c
	real*4 bk_sw3d(ni,nj,2)
	integer*4 lvl_sw3d(2)
	character var_sw3d(2)*3, backsw3d*9, dir_sw3d*256, ext_sw3d*31
	character units3(2)*10, lvlc3(2)*4, com3(2)*125
c
c..... Stuff for the sfc data and other station info (LSO +)
c
	real*4 lat_s(mxstn), lon_s(mxstn), elev_s(mxstn)
	real*4 t_s(mxstn), td_s(mxstn), dd_s(mxstn), ff_s(mxstn)
	real*4 ddg_s(mxstn), ffg_s(mxstn)
	real*4 pstn_s(mxstn), pmsl_s(mxstn), alt_s(mxstn)
	real*4 cover_s(mxstn), hgt_ceil(mxstn), hgt_low(mxstn)
	real*4 solar_s(mxstn), store_hgt(mxstn,5), vis_s(mxstn)
	real*4 rii(mxstn), rjj(mxstn)
c
	integer*4 kloud_s(mxstn),idp3_s(mxstn),obstime(mxstn)
	integer*4 ii(mxstn), jj(mxstn)
c
	character atime_s*24
	character stn(mxstn)*3,obstype(mxstn)*8,wx_s(mxstn)*8
	character store_emv(mxstn,5)*1, store_amt(mxstn,5)*4
c
c.....  Work arrays for the QC routine.
c
        integer*4 rely(26,mxstn), ivals1(mxstn)
c
c.....  Stuff for intermediate grids (old LGS file)
c
	real*4 u1(ni,nj), v1(ni,nj)
	real*4 t1(ni,nj), td1(ni,nj), tb81(ni,nj)
	real*4 rp1(ni,nj), sp1(ni,nj), mslp1(ni,nj)
	real*4 vis1(ni,nj), elev1(ni,nj)
c
c
c.....	Start here.  First see if this is an interactive run.
c
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
	del = 1.e6
	gam = .0008
	ak = 1.e-6
c
c.....	Get the LAPS lat/lon and topo data here so we can pass them to the 
c.....	routines that need them.
c
cc	dir_s = '../static/' 
	call get_directory('static', dir_s, len)
	ext_s = laps_domain
	var_s = 'LAT'
        call rd_laps_static(dir_s,ext_s,ni,nj,1,var_s,units,comment,
     &                      lat ,grid_spacing,istatus)
	var_s = 'LON'
        call rd_laps_static(dir_s,ext_s,ni,nj,1,var_s,units,comment,
     &                      lon ,grid_spacing,istatus)
	var_s = 'AVG'
        call rd_laps_static(dir_s,ext_s,ni,nj,1,var_s,units,comment,
     &                      topo ,grid_spacing,istatus)
c
c.....	Find east/west and north/south sides of grid (max extension of grid)
c
	grid_east = -999.
	grid_west = 0.
	grid_north = 0.
	grid_south = 90.
	do i=1,ni
	  if(lat(i,nj) .gt. grid_north) grid_north = lat(i,nj)
	  if(lat(i,1) .lt. grid_south) grid_south = lat(i,1)
	enddo !i
	do j=1,nj	
	  if(lon(ni,j) .gt. grid_east) grid_east = lon(ni,j)
	  if(lon(1,j) .lt. grid_west) grid_west = lon(1,j)
	enddo !j
c
c.....  Read in the obs and calculate a weight based on distance to each
c.....  station.
c
c.....	READ IN THE SURFACE OBS:  dd/ff in deg/kt, t and td in F, elev in m,
c.....	                          and the pressure variable. cld hts are msl.
c
c 	infile1 = '../lapsprd/lso/'//filename//'.lso' 
	call get_directory('lso',infile1,len)
	infile1 = infile1(1:len) // filename(1:9) // '.lso'
c
	write(6,305) infile1
 305	format(' Getting surface data from: ',a70)
	call read_surface_obs(infile1,mxstn,atime_s,n_meso_g,
     &   n_meso_pos,n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,
     &   n_obs_pos_g,n_obs_b,n_obs_pos_b,stn,obstype,lat_s,lon_s,elev_s,
     &   wx_s,t_s,td_s,dd_s,ff_s,ddg_s,ffg_s,pstn_s,pmsl_s,alt_s,
     &   kloud_s,hgt_ceil,hgt_low,cover_s,solar_s,idp3_s,store_emv,
     &   store_amt,store_hgt,vis_s,obstime,istatus)
c
	if(istatus.ne.1 .or. n_obs_b.eq.0) then	  !surface obs not available
	  jstatus(1) = 0
	  stop 'No sfc obs from LSO'
	endif
c
	print *,' '
	write(6,320) atime_s,n_sao_g,n_sao_b,n_obs_b
320	format(' LSO data vaild time: ',a24,' Num obs: ',3i6)
c
	if(n_sao_b .lt. 10) then
	   jstatus(1) = -2
	   print *,' Insufficient number of surface observations'
	   print *,' for a clean analysis.  Stopping.'
	   stop 
	endif
	print *,' '
c
c.....	Find the i,j location of each station, then calculate the
c.....  background weights (based on station density).
c
	call find_ij(lat_s,lon_s,lat,lon,n_obs_b,mxstn,
     &               ni,nj,ii,jj,rii,rjj)
c
        do ista=1,n_obs_b
           write(6,999) ista,stn(ista),rii(ista),rjj(ista)
        enddo !ista
 999    format(i4,': ',a3,' is at i,j: ',f5.1,',',f5.1)
c
        call zero(wt, ni,nj)
	call bkgwts(lat,lon,topo,n_obs_b,lat_s,lon_s,elev_s,
     &              rii,rjj,wt,ni,nj,mxstn)
c
c.....  Zero the weights then get the background data.  Try for a RAMS 
c.....  forecast first, then the surface winds from the 3d analysis, then
c.....  a previous LAPS analysis.  If RAMS or the 3d wind are missing just 
c.....  use LAPS.  If both LAPS and RAMS are missing, print a warning.
c
        call zero(wt_u, ni,nj)
        call zero(wt_v, ni,nj)
        call zero(wt_t, ni,nj)
        call zero(wt_td, ni,nj)
        call zero(wt_rp, ni,nj)
        call zero(wt_mslp, ni,nj)
        call zero(wt_stnp, ni,nj)
        call zero(wt_vis, ni,nj)
c
	if(ihours .eq. 0) then     ! skip the backgrounds
	   irams_bk = 0
	   ilaps_bk = 0
	   isw3d_bk = 0
	   print *,' Skipping backgrounds...'
	   go to 600
	endif
c
	irams_bk = 1
	irams_loop = 1
	do i=1,5
	   lvl_bkr(i) = 0
	enddo !i
c
	var_bkr(1) = 'U'       ! u-wind (m/s)
	var_bkr(2) = 'V'       ! v-wind (m/s)
	var_bkr(3) = 'RP'      ! reduced pressure (Pa) (1500 m in CO)
	var_bkr(4) = 'T'       ! temp (K)
	var_bkr(5) = 'TD'      ! dewpt (K)
	i4time_rams = i4time_last + 60     ! get the 1-h fcst
c	dir_rams = '../lapsprd/rsf/'
	call get_directory('rsf',dir_rams,len)
	ext_rams = 'rsf'
 100	call read_laps_data(i4time_rams,dir_rams,ext_rams,ni,nj,5,5,
     &   var_bkr,lvl_bkr,lvlcr,unitsr,comr,bk_rams,istatus)
	if(istatus .ne. 1) then
	   if(irams_loop .gt. 5) then     ! give up
	      print *,' ++ No RAMS backgound available. ++'
	      irams_bk = 0
	      go to 300
	   else
	      irams_loop = irams_loop + 1
	      print *,' RAMS forecast not available...trying earlier.'
	      i4time_rams = (i4time_rams - 3600) + 60
	      go to 100
	   endif
	endif
c
c.....	Get the previous LSX.
c
 300	ilaps_bk = 1
	ilaps_loop = 1
	do i=1,8
	   lvl_bk(i) = 0
	enddo !i
	var_bk(1) = 'U'     ! u-wind (m/s)
	var_bk(2) = 'V'     ! v-wind (m/s)
	var_bk(3) = 'P'     ! reduced pressure (Pa) (1500 m in CO)
	var_bk(4) = 'T'     ! temp (K)
	var_bk(5) = 'TD'    ! dew point (K)
	var_bk(6) = 'MSL'   ! msl pressure (Pa)
	var_bk(7) = 'VIS'   ! visibility (m)
	var_bk(8) = 'PS'    ! station pressure (Pa)
c	dir_bk = '../lapsprd/lsx/'
	call get_directory('lsx',dir_bk,len)
	ext_bk = 'lsx'
	i4time_bk = i4time_last 
 200	call read_laps_data(i4time_bk,dir_bk,ext_bk,ni,nj,8,8,var_bk,
     &     lvl_bk,lvlcl,unitsl,coml,background,istatus)
	if(istatus .ne. 1) then
	   if(ilaps_loop .gt. 3) then     ! give up
	      print *,' ++ No  LSX background available. ++'
	      ilaps_bk = 0
	      go to 350
	   else
	      ilaps_loop = ilaps_loop + 1
	      ihours = ihours + 1
	      dt = 3600. * ihours
	      i4time_bk = i4time - ifix( dt )
	      go to 200
	   endif
	endif
c
 350	continue
 	isw3d_bk = 1
	isw3d_loop = 1
	do i=1,2
	   lvl_sw3d(i) = 0
	enddo !i
	var_sw3d(1) = 'SU'     ! u-wind (m/s)
	var_sw3d(2) = 'SV'     ! v-wind (m/s)
c	dir_sw3d = '../lapsprd/lwm/'
	call get_directory('lwm',dir_sw3d,len)
	ext_sw3d = 'lwm'
	i4time_sw3d = i4time_last
 202	call read_laps_data(i4time_sw3d,dir_sw3d,ext_sw3d,ni,nj,2,2,
     &     var_sw3d,lvl_sw3d,lvlc3,units3,com3,bk_sw3d,istatus)
	if(istatus .ne. 1) then
	   if(isw3d_loop .gt. 3) then     ! give up
	      print *,' ++ No  LWM background available. ++'
	      isw3d_bk = 0
	      go to 370
	   else
	      isw3d_loop = isw3d_loop + 1
	      ihours = ihours + 1
	      dt = 3600. * ihours
	      i4time_sw3d = i4time - ifix( dt )
	      go to 202
	   endif
	endif
c
 370	continue
c
c.....  Now convert units and move stuff to the proper arrays.
c      
	call constant(u_bk,badflag,ni,nj)
	call constant(v_bk,badflag,ni,nj)
	call constant(t_bk,badflag,ni,nj)
	call constant(td_bk,badflag,ni,nj)
	call constant(rp_bk,badflag,ni,nj)
	call constant(mslp_bk,badflag,ni,nj)
	call constant(stnp_bk,badflag,ni,nj)
	call constant(vis_bk,badflag,ni,nj)
c
	if(ilaps_bk.eq.0 .and. irams_bk.eq.0) then
	   print *,
     &     ' ++WARNING. No backgrounds available for the analysis.++'
	   go to 600
	endif
	if(irams_bk .eq. 1) then     ! have RAMS for background
	   call move_3dto2d(bk_rams,1, u_bk, ni,nj,nk)    !rams u
	   call move_3dto2d(bk_rams,2, v_bk, ni,nj,nk)    !rams v
	   call move_3dto2d(bk_rams,3,rp_bk, ni,nj,nk)    !rams 1500 p
	   call move_3dto2d(bk_rams,4, t_bk, ni,nj,nk)    !rams t
	   call move_3dto2d(bk_rams,5,td_bk, ni,nj,nk)    !rams td
c
           call move(wt, wt_u, ni,nj)
           call move(wt, wt_v, ni,nj)
           call move(wt, wt_rp, ni,nj)
           call move(wt, wt_t, ni,nj)
           call move(wt, wt_td, ni,nj)
	endif
c
	if(ilaps_bk .eq. 1) then
	   if(irams_bk .eq. 1) then  ! have RAMS, just fill in the others
	      call move_3dto2d(background,6, mslp_bk, ni,nj,nk)  !laps msl p
	      call move_3dto2d(background,7,  vis_bk, ni,nj,nk)  !laps vis
	      call move_3dto2d(background,8, stnp_bk, ni,nj,nk)  !laps stn p
c
              call move(wt, wt_mslp, ni,nj)
              call move(wt, wt_vis, ni,nj)
              call move(wt, wt_stnp, ni,nj)
c
	   else                      ! only previous laps for background
c
	      if(isw3d_bk .eq. 1) then      !but...have 3d sfc wind
		 call move_3dto2d(bk_sw3d,1, u_bk, ni,nj,nk) !lwm 3d sfc u
		 call move_3dto2d(bk_sw3d,2, v_bk, ni,nj,nk) !lwm 3d sfc v
	      else
		 call move_3dto2d(background,1, u_bk, ni,nj,nk) !laps u
		 call move_3dto2d(background,2, v_bk, ni,nj,nk) !laps v
	      endif
c
	      call move_3dto2d(background,3,   rp_bk, ni,nj,nk)  !laps reduced p (co: 1500m)
	      call move_3dto2d(background,4,    t_bk, ni,nj,nk)  !laps t
	      call move_3dto2d(background,5,   td_bk, ni,nj,nk)  !laps td
	      call move_3dto2d(background,6, mslp_bk, ni,nj,nk)  !laps msl p
	      call move_3dto2d(background,7,  vis_bk, ni,nj,nk)  !laps vis
	      call move_3dto2d(background,8, stnp_bk, ni,nj,nk)  !laps stn p
c
              call move(wt, wt_u, ni,nj)
              call move(wt, wt_v, ni,nj)
              call move(wt, wt_rp, ni,nj)
              call move(wt, wt_t, ni,nj)
              call move(wt, wt_td, ni,nj)
              call move(wt, wt_mslp, ni,nj)
              call move(wt, wt_vis, ni,nj)
              call move(wt, wt_stnp, ni,nj)
	   endif
	endif
c
	call conv_ms2kt(u_bk,u_bk,ni,nj)
	call conv_ms2kt(v_bk,v_bk,ni,nj)
	call conv_k2f(t_bk,t_bk,ni,nj)
	call conv_k2f(td_bk,td_bk,ni,nj)
	call multcon(rp_bk,0.01,ni,nj)          ! conv Pa to mb
	call multcon(mslp_bk,0.01,ni,nj)        ! conv Pa to mb
	call multcon(stnp_bk,0.01,ni,nj)        ! conv Pa to mb
	call conv_m2miles(vis_bk,vis_bk,ni,nj)
	call visg2log(vis_bk,ni,nj,badflag)     ! conv miles to log(miles)
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
c.....  Print the background times.
c
	print *,' '
	if(ilaps_bk .eq. 1) then
	   call make_fnam_lp(i4time_bk,back,istatus)
	   write(6,400) back
 400	   format(' Using  LSX background from: ',a9)
	else
	   print *,' ++ No  LSX background available. ++'
	endif
	if(irams_bk .eq. 1) then
	   call make_fnam_lp(i4time_rams,backrams,istatus)
	   write(6,405) backrams
 405	   format(' Using RAMS background from: ',a9)
	else
	   print *,' ++ No RAMS background available. ++'
	endif
	if(isw3d_bk .eq. 1) then
	   call make_fnam_lp(i4time_sw3d,backsw3d,istatus)
	   write(6,407) backsw3d
 407	   format(' Using LWM sfc wind background from: ',a9)
	else
	   print *,' ++ No LWM sfc wind background available. ++'
	endif
	print *,' '
c
 600	continue
c
c.....	QC the surface data.
c
	if(iskip .gt. 0) then  !check QC flag
	  print *, ' **  omit qc of data  **  '
	  goto 521
	endif
	call get_directory('lso', infile_last, len)
	infile_last = infile_last(1:len) // filename_last // '.lso'
c
	call qcdata(filename,infile_last,rely,ivals1,mxstn,
     &     lat_s, lon_s, elev_s, t_s, td_s, dd_s, ff_s, ddg_s, 
     &     ffg_s, pstn_s, pmsl_s, alt_s, cover_s, hgt_ceil, 
     &     hgt_low, solar_s, store_hgt, vis_s, kloud_s, idp3_s, 
     &     obstime, stn, obstype, wx_s, store_emv, store_amt,
     &     rii, rjj, ii, jj, n_obs_b, n_sao_b, n_sao_g,
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
     &    ' +++ ERROR.  Problem in QC routine. +++'
	  jstatus(2) = -2
	  go to 521
	endif
c
c.....	Check each of the primary analysis variables.
c
	do mm=1,n_obs_b
	  nn = ivals1(mm)
	  if(nn .lt. 1) go to 121
	  if(rely(7,nn) .lt. 0) then	! temperature
	    print *, 'QC: Bad T at ',stn(mm),' with value ',t_s(mm)
	    t_s(mm) = badflag
	  endif
 121	enddo  !mm
	do mm=1,n_obs_b
	  nn = ivals1(mm)
	  if(nn .lt. 1) go to 122
	  if(rely(8,nn) .lt. 0) then	! dewpt
	    print *, 'QC: Bad TD at ',stn(mm),' with value ',td_s(mm)
	    td_s(mm) = badflag
	  endif
 122	enddo  !mm
	do mm=1,n_obs_b
	  nn = ivals1(mm)
	  if(nn .lt. 1) go to 123
	  if(rely(9,nn) .lt. 0) then	! wind direction
	    print *, 'QC: Bad DIR at ',stn(mm),' with value ',dd_s(mm)
	    dd_s(mm) = badflag
	  endif
 123	enddo  !mm
	do mm=1,n_obs_b
	  nn = ivals1(mm)
	  if(nn .lt. 1) go to 124
	  if(rely(10,nn) .lt. 0) then	! wind speed
	    print *, 'QC: Bad SPD at ',stn(mm),' with value ',ff_s(mm)
	    ff_s(mm) = badflag
	  endif
 124	enddo  !mm
	do mm=1,n_obs_b
	  nn = ivals1(mm)
	  if(nn .lt. 1) go to 126
	  if(rely(15,nn) .lt. 0) then	! altimeter 
	    print *, 'QC: Bad ALT at ',stn(mm),' with value ',alt_s(mm)
	    alt_s(mm) = badflag
	  endif
 126	enddo  !mm
	do mm=1,n_obs_b
	  nn = ivals1(mm)
	  if(nn .lt. 1) go to 127
	  if(rely(17,nn) .lt. 0) then	! ceiling 
	   print *, 'QC: Bad CEIL at ',stn(mm),' with value ',hgt_ceil(mm)
	   hgt_ceil(mm) = badflag
	  endif
 127	enddo  !mm
	do mm=1,n_obs_b
	  nn = ivals1(mm)
	  if(nn .lt. 1) go to 128
	  if(rely(25,nn) .lt. 0) then	! visibility 
	    print *, 'QC: Bad VIS at ',stn(mm),' with value ',vis_s(mm)
	    vis_s(mm) = badflag
	  endif
 128	enddo  !mm
 521	continue                          
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
     &     lon,grid_east,grid_west,grid_north,grid_south,topo,x1a,x2a,
     &     y2a, lon_s, elev_s, t_s, td_s, dd_s, ff_s, pstn_s, pmsl_s, 
     &     alt_s, vis_s, stn, rii, rjj, ii, jj, n_obs_b, n_sao_g,
     &     u_bk, v_bk, t_bk, td_bk, rp_bk, mslp_bk, stnp_bk, vis_bk, 
     &     wt_u, wt_v, wt_rp, wt_mslp, ilaps_bk, irams_bk,
     &     u1, v1, rp1, t1, td1, sp1, tb81, mslp1, vis1, elev1,
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
	call laps_vanl(i4time,filename,ni,nj,nk,mxstn,laps_cycle_time,
     &        dt,del,gam,ak,lat,lon,topo,grid_spacing, laps_domain,
     &        lat_s, lon_s, elev_s, t_s, td_s, ff_s, pstn_s, vis_s, 
     &        stn, n_obs_b, n_sao_b, n_sao_g,
     &        u_bk, v_bk, t_bk, td_bk, rp_bk, mslp_bk, vis_bk, 
     &        wt_u, wt_v, wt_t, wt_td, wt_rp, wt_mslp, wt_vis, 
     &        ilaps_bk, irams_bk,
     &        u1, v1, rp1, t1, td1, sp1, tb81, mslp1, vis1, elev1,
     &        x1a,x2a,y2a,ii,jj,jstatus)
c
	if(jstatus(3) .ne. 1) then
	  print *,' From LAPS_VANL: Error Return.' 
	endif
	if(jstatus(4) .ne. 1) then
	  print *,' From LAPS_VANL LT1:  Error Return.' 
	endif
c
c.....	That's about it...let's go home.
c
	return
	end

