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
	subroutine laps_vanl(i4time,infile,dir_in,ext_in,dir,ext,
     &        ihrs,dt,del,gam,ak,lat,lon,topo,grid_spacing,jstatus)
c
c******************************************************************************
c
c	Mesoanalysis Model developed by J.McGinley
c
c	Uses variational constraints to impose approximations of the
c	equations of motion in the analysis of u,v,and p.  Input
c	required is a first guess uo,vo,and po and the same fields 
c	for 1 to 6 hours previous to analysis time, um2, and vm2. 
c	The analysis is controlled by three parameters: 1. gam which 
c	determines the extent to which pressure is modified 2. del, 
c	which determines the magnitude of the residual of the momentum 
c	equations. 3.  tfact, the delta t factor.  This is the time 
c	scale of the input data, the delta t of the wind change.
c	Artificiallly varying this parameter can modulate the effects of
c	local velocity change and provide a less or more geostrophic
c	analysis.  Suggested values: 
c	gam .01 to .00001 m**2/sec**2/pa**2
c	should be the squared error of the winds divided by the 
c	squared error of the pressure measurements.  For the surface 
c	gam=1m**2/sec**2/2500pa**2=.0004
c	del 1.e6 to 1.e8 sec**-2
c	controls the rms error of the equation of motion over the 
c	domain.  best value 1.e7.  For gam=const increasing del 
c	will result in increased adjustment to winds and pressure
c	tfact= 7200 sec (2 hours)
c	nonlinear and friction terms added.  Non linear terms are
c	computed from the first guess fields.  Frictional term assumes
c	ak=1.e-5.  Fu=U*absU*ak.  The algorithm can be iterated on
c	the non linear terms if desired.  set npass>1.  Omega is computed
c	at the top of a 50 mb deep boundary layer.(pa/sec)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c	Changes made:
c		P.A. Stamus	11-19-86  Added graphics for fields.
c				07-07-87  Added read for del,gam,ci's.
c				09-28-87  Added wind barb plotting stuff.
c				10-02-87  Included wind change vector plot.
c				10-06-87  Added code to plot del,gam,dt,ak.
c				10-26-87  Added Coriolis term (oops!)
c				11-17-87  Added dp code.
c				01-19-88  LAPS version.
c				03-22-88  Changed lat/lon to LAPS standard.
c				06-02-88  Changes for new LAPS grid.
c				06-27-88  More changes for new grid.
c				08-31-88  Take out graphics (now in LAPSplot)
c					   and store analysis fields; clean
c					   things up a bit.
c				10-11-88  make filenames time dependent
c-----------------------------------------------------------------------------
c		P. A. Stamus 	01-10-89  Rewritten and combined with the
c					  LAPSint program for batch mode.
c				04-12-89  Added flags to header, changed time
c					  file method to call i4time.
c				04-28-89  Added writes to PROD_DEV for primary
c					  variables and cld top & base.
c				05-03-89  Made sure everything was < 72 cols
c				05-08-89  Added meso_anl derived fields.
c				05-11-89  Add Moisture advect., LI, and write
c					  outputs in LAPS standards.
c				06-01-89  New grid -- add nest6grid.
c				11-07-89  Add nummeso/sao to headers.
c				12-06-89  Change cld hts to msl in, both out.
c				12-08-89  Change beta in cld splines to .3
c				02-26-90  Write all fields to PROD_DEV.
c				----------------------------------------------
c				03-14-90  Subroutine version.
c				04-06-90  Pass in del,gam,ak from driver.
c				04-10-90  Loosen spline (make beta=1.)
c				04-11-90  Pass in LAPS lat/lon and topo.
c				04-17-90  Bag clouds, strmfct, vel pot; add 
c					  RH, Tempadv, MSL P.  keep ceil.
c				04-25-90  Add data quality array, etc.
c				05-01-90  Changed data quality array...
c				06-15-90  New version of temperature calcs
c				06-29-90  Test of new dewpoints
c				07-13-90  Ck for earlier hrs in 1 hr missing.
c				10-09-90  Bring in MAPS 700 T for anl.
c				10-22-90  Install background fields.
c				10-30-90  Bag quality array..assume good input.
c				04-10-91  Correct units probs with meso_anl.
c				05-03-91  Add CSSI,P/NBE code.
c				09-30-91  Add dummy arrays for routines.
c				11-13-91  Changes for new grids.
c				12-16-91  Changes for put_tmp_anl.
c				01-15-92  Add visibility analysis.
c				03-23-92  Fix pressure in get_hts_hydro.
c				07-28-93  New Barnes routines. 
c                               12-08-93  Changes for new put_tmp_anl.
c                               02-24-94  Upgrade enhance_vis routine,bag ceil
c                               03-02-94  Install RAMS backgrnd stuff.
c                               04-14-94  Changes for CRAY port.
c                               07-20-94  Add include file.
c                               08-31-94  Change to LGA from LMA for bkgs.
c                               02-02-95  Background calls to driver routine.
c                               07-20-95  Heat Indx, chng for new Barnes mthd.
c                               08-08-95  Add verify code, LSO read.
c                               05-23-96  Convert RH to % for some reason.
c                               08-08-96  Change LGA read for Upper variables.
c                               08-12-96  Add 500 theta check for sfc temps.
c                               09-10-96  Subroutine Channel corrected
c                               10-09-96  Gridded stn elevs for temp anl.
c                               11-07-96  Num obs ck.
c                               12-17-96  More porting changes...common for
c                                           sfc data, LGS grids. 
c
c*****************************************************************************
cx
	include 'laps_sfc.inc'
c
	parameter(roi = 20)	  !radius of influence for Barnes 
	parameter( 	          !QC parameters: # of standard deviations 
     &            bad_p  = 5.0,		! for reduced pressure
     &            bad_t  = 3.5,		! for temperature
     &            bad_td = 2.5,		! for dewpoint
     &            bad_u  = 4.0,		! for u-wind
     &            bad_v  = 4.0,		! for v-wind
     &            bad_th = 3.5,		! for theta
     &            bad_the = 2.5,	! for theta-e
     &            bad_vis = 500.,	! for visibility
     &            bad_tb8 = 5.0)	! for tb8 Brightness temps.
c
	parameter(bad = 1.e6 - 2.)	! larger than 'bad' are.
c
c.....	LAPS lat/lon and terrain grids, and Coriolis.
c
	integer*4 istatus
	real*4 lat(ni,nj), lon(ni,nj), topo(ni,nj)
	real*4 fo(ni,nj), fo2(ni,nj), akk(ni,nj), pbl_top(ni,nj)
c
c.....  Stuff for intermediate grids (old LGS file)
c
	real*4 u1(ni,nj), v1(ni,nj)
	real*4 t1(ni,nj), td1(ni,nj), tb81(ni,nj)
	real*4 rp1(ni,nj), sp1(ni,nj), mslp1(ni,nj)
	real*4 vis1(ni,nj), elev1(ni,nj)
c
	common/LGS_grids/
     &     u1, v1, rp1, t1, td1, sp1, tb81, mslp1, vis1, elev1
c
c.....	Grids for the first data's analyses.
c
	real*4 u(ni,nj), v(ni,nj)
	real*4 rp(ni,nj), psfc(ni,nj), vis(ni,nj)
	real*4 t(ni,nj), theta(ni,nj), thetae(ni,nj), tb8(ni,nj)
	real*4 td(ni,nj), ceil(ni,nj), mslp(ni,nj)
c
c.....	Grids for the variational analyses of rp, u, v
c
	real*4 p_a(ni,nj), u_a(ni,nj), v_a(ni,nj)
c
c.....	Grids for the derived quantities.
c
	real*4 du(ni,nj), dv(ni,nj), spd(ni,nj)
	real*4 drp(ni,nj), w(ni,nj)
	real*4 tt(ni,nj), ttd(ni,nj)
	real*4 li(ni,nj), qadv(ni,nj), rh(ni,nj)
	real*4 cssi(ni,nj), fire(ni,nj), hi(ni,nj)
	real*4 p_1d_pa(nk)
	real*4 pbe_2d(ni,nj), nbe_2d(ni,nj)
c
c.....  Grids for PUT_TEMP_ANAL and related stuff.
c
        real*4 LT1_out(ni,nj,nk,2), grid_spacing
	real*4 t_3d_k(ni,nj,nk), ht_3d_m(ni,nj,nk), rh_3d(ni,nj,nk)
cc        equivalence(t_3d_k,   LT1_out(1,1,1,1))
cc        equivalence(ht_3d_m,  LT1_out(1,1,1,2))
        integer*4 lni_m, lnj_m
        parameter(lni_m = ((ni-1)/iden_ratio)+1)
        parameter(lnj_m = ((nj-1)/iden_ratio)+1)
        integer*4 adbox(lni_m,lnj_m), ndbox(lni_m,lnj_m)
c
c.....	Grids for variables derived by the MESO_ANL subroutine.
c
	real*4 q(ni,nj), qcon(ni,nj), thadv(ni,nj), tadv(ni,nj)
c
c.....	Grids for the background fields.
c
        real*4 u_bk(ni,nj), v_bk(ni,nj), t_bk(ni,nj), td_bk(ni,nj)
        real*4 wt_u(ni,nj), wt_v(ni,nj), wt_t(ni,nj), wt_td(ni,nj)
        real*4 rp_bk(ni,nj), mslp_bk(ni,nj), stnp_bk(ni,nj)
        real*4 wt_rp(ni,nj), wt_mslp(ni,nj), wt_stnp(ni,nj)
        real*4 vis_bk(ni,nj), wt_vis(ni,nj)
c
        common/backgrnd/
     &     u_bk, v_bk, t_bk, td_bk, rp_bk, mslp_bk, stnp_bk, vis_bk, 
     &     wt_u, wt_v, wt_t, wt_td, wt_rp, wt_mslp, wt_stnp, wt_vis, 
     &     ilaps_bk, irams_bk
c
c.....	Grids for other stuff.
c
	real*4 ddiv(ni,nj), vort(ni,nj), st(ni,nj), sm(ni,nj)
	real*4 f(ni,nj), fu(ni,nj), fv(ni,nj), div(ni,nj)
	real*4 a(ni,nj), z(ni,nj), dx(ni,nj), dy(ni,nj)
	real*4 nu(ni,nj),nv(ni,nj), t7(ni,nj), h7(ni,nj), td7(ni,nj)
	real*4 t5(ni,nj)
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
	character infile_o*70,atime_s*24
	character stn(mxstn)*3,obstype(mxstn)*8,wx_s(mxstn)*8
	character store_emv(mxstn,5)*1, store_amt(mxstn,5)*4
c
	common/LSO_sfc_obs/
     &     lat_s, lon_s, elev_s, t_s, td_s, dd_s, ff_s, ddg_s, 
     &     ffg_s, pstn_s, pmsl_s, alt_s, cover_s, hgt_ceil, 
     &     hgt_low, solar_s, store_hgt, vis_s, kloud_s, idp3_s, 
     &     obstime, stn, obstype, wx_s, store_emv, store_amt,
     &     rii, rjj, ii, jj, n_obs_b, n_sao_b, n_sao_g
c
        character atime_s*24, title*40, ver_file*200
c
c.....	dummy work arrays
c
	real*4 d1(ni,nj),d2(ni,nj),d3(ni,nj),d4(ni,nj),d5(ni,nj)
        real*4 d6(ni,nj), d7(ni,nj)
	real*4 dm1(ni_maps,nj_maps,nk), dm2(ni_maps,nj_maps,nk)
	real*4 dm3(ni_maps,nj_maps,nk), dm4(ni,nj,nk)
c
        include 'laps_cloud.inc'
        real*4 d_c(ni,nj,kcloud)
        integer*4 kdum(ni,nj)
c
	real*4 lapse_t, lapse_td
	real make_td
	character infile*200, name*10, filename*9, atime*24
	character var_fire*3, com_fire*125, units_fire*10, ext_f*31
	character var_lga*3, ext_lga*31
	integer*4 jstatus(20)
	logical l_fill
c
	character*80 grid_fnam_common
	common/ grid_fnam_cmn / grid_fnam_common
c
c.....	Stuff for LAPS outputs (i.e., standard forms).
c
	real*4 data(ni,nj,27)
	integer*4 imax,jmax,lvl(27)
	character dir*50,ext*31,var(27)*3,lvl_coord(27)*4,units(27)*10
	character comment(27)*125,dir_in*50,ext_in*31
	equivalence(u_a,    data(1,1,1))
	equivalence(v_a,    data(1,1,2))
	equivalence(p_a,    data(1,1,3))
	equivalence(t,      data(1,1,4))
	equivalence(td,     data(1,1,5))
	equivalence(w,      data(1,1,6))
	equivalence(rh,     data(1,1,7))
	equivalence(ceil,   data(1,1,8))
	equivalence(mslp,   data(1,1,9))
	equivalence(tadv,   data(1,1,10))
	equivalence(theta,  data(1,1,11))
	equivalence(thetae, data(1,1,12))
	equivalence(psfc,   data(1,1,13))
	equivalence(vort,   data(1,1,14))
	equivalence(q,      data(1,1,15))
	equivalence(qcon,   data(1,1,16))
	equivalence(div,    data(1,1,17))
	equivalence(thadv,  data(1,1,18))
	equivalence(qadv,   data(1,1,19))
	equivalence(li,     data(1,1,20))
	equivalence(spd,    data(1,1,21))
	equivalence(cssi,   data(1,1,22))
	equivalence(pbe_2d, data(1,1,23))
	equivalence(nbe_2d, data(1,1,24))
	equivalence(vis,    data(1,1,25))
	equivalence(fire,   data(1,1,26))
	equivalence(hi,     data(1,1,27))
c
c
c.....	Start...set up constants, initialize arrays, etc.
c
	ibt = 1      !assume have sat data...code cks later.
	grid_fnam_common = laps_domain
	jstatus(3) = -1		 ! start w/this until changed
	pi = 4. * atan(1.)
	smsng = badflag
	imax = ni
	jmax = nj
	kmax = nk
	imax_m = ni_maps
	jmax_m = nj_maps
	beta = 3.
	betac = .3
	itmax = 100
	alf = 100.
	alfg = 1000.
	err = .1
	ovr = 1.4
	scale = 0.
	ci = 0.
	npass = 1
	dpbl = 5000.
	rho = 1.25
	rho2 = rho * rho
	omg2 = 2. * 7.292e-5
	rdpdg = pi / 180.
	re = 6371222.
	gor = 9.808 / 287.04
	fon = 5. / 9.
	anof = 9. / 5.
	call zero(z,imax,jmax)
        call zero(ceil,imax,jmax)   ! dummy ceil ht array
	pblht = 500.		! pbl height in meters
	call make_fnam_lp(i4time,filename,istatus)
c
c.....  calculate Coriolis term for each grid point
c
  	do j=1,jmax
	do i=1,imax
	   fo(i,j) = omg2 * sin( rdpdg * lat(i,j) )
	   fo2(i,j) = fo(i,j) * fo(i,j)
	enddo !i
	enddo !j
c
c.....	set up some other stuff (including grid intervals, etc.)
c
	do j=1,jmax-1
	do i=1,imax
	  dlat = lat(i,j+1) - lat(i,j)
	  alat = (lat(i,j+1) + lat(i,j)) * .5
	  dy(i,j) = dlat * rdpdg * re
	  dx(i,j) = dy(i,j) * cos( rdpdg * alat )
	  a(i,j) = -(1. + fo2(i,j) * del) * gam * rho2 / del
	enddo !i
	enddo !j
c
	do i=1,imax
	  dy(i,jmax) = dy(i,jmax-1)
	  dx(i,jmax) = dx(i,jmax-1)
	  a(i,jmax) = -(1. + fo2(i,jmax) * del) * gam * rho2 / del
	enddo !i
c
cz..... Compute T on the surface using the LGA (or equiv) 700 T and HT.
c
	i4time_tol = 21600
	ext_lga = 'lga'
	l_fill = .true.
c
c.....  Get the latest 3d fields, pull out the var/lvls needed.
c
	print *,' Get LGA 700 T'
	itheta7 = 1
	var_lga = 'T3 '
	call get_maps_laps_4d(i4time,var_lga,ni,nj,nk,ni_maps,
     &                nj_maps,dm1,dm2,dm3,dm4,l_fill,istatus)
c
	if(istatus .ne. 1)  then
	   print *,' LGA 700 T not available. Using constant 5C.'
	   call constant(t7,278.15,imax,jmax)
	   itheta7 = 0
	else
	   do j=1,nj
	   do i=1,ni
	      t7(i,j) = dm4(i,j,9) ! lvl 9 = 700 hPa
	   enddo !i
	   enddo !j
	endif
c
	print *,' Get LGA 700 HT'
	var_lga = 'HT '
	call get_maps_laps_4d(i4time,var_lga,ni,nj,nk,ni_maps,
     &                nj_maps,dm1,dm2,dm3,dm4,l_fill,istatus)
c
	if(istatus .ne. 1)  then
	 print *,' LGA 700 HT not available. Using constant 3000 m.'
	 call constant(h7,3000.,imax,jmax)
	else
	   do j=1,nj
	   do i=1,ni
	      h7(i,j) = dm4(i,j,9) ! lvl 9 = 700 hPa
	   enddo !i
	   enddo !j
	endif
c
	print *,' Get LGA 700 TD'
	var_lga = 'SH '  ! specific humidity 
	call get_maps_laps_4d(i4time,var_lga,ni,nj,nk,ni_maps,
     &                nj_maps,dm1,dm2,dm3,dm4,l_fill,istatus)
c
	if(istatus .ne. 1)  then
	   print *,' LGA 700 Td not available. Using constant -5C.'
	   call constant(td7,268.15,imax,jmax) 
	else
	   do j=1,nj
	   do i=1,ni
	    t7c = t7(i,j) - 273.15  !K to C
	    qgkg = dm4(i,j,9) * 1000.   !lvl 9 = 700 hPa
	    td7(i,j) = make_td(700., t7c, qgkg, 0.) + 273.15 !in K
	   enddo !i
	   enddo !j
	endif
c
c.....  Get the 500 Temps while we're here
c
	print *,' Get LGA 500 T'
	itheta5 = 1
	var_lga = 'T3 '
	call get_maps_laps_4d(i4time,var_lga,ni,nj,nk,ni_maps,
     &                nj_maps,dm1,dm2,dm3,dm4,l_fill,istatus)
c
	if(istatus .ne. 1)  then
	 print *,' LGA 500 T not available.'
	 itheta5 = 0
	else
	   do j=1,nj
	   do i=1,ni
	      t5(i,j) = dm4(i,j,13) ! lvl 13 = 500 hPa
	   enddo !i
	   enddo !j
	endif
c
c.....  Convert units.
c
	call conv_k2f(td7,td7,imax,jmax)  ! conv K to F
	call conv_k2f(t7,t7,imax,jmax)	  ! conv K to F
c
c.....  Get lapse rate (usually std), and mean pressure.
c
	call mean_lapse(n_obs_b,elev_s,t_s,td_s,a_t,lapse_t,a_td,
     &                  lapse_td,hbar)
	call mean_pres(n_obs_b,pstn_s,pbar)
c
c.....  Calculate the terrain est temps and pressure.
c
	do j=1,jmax
	do i=1,imax
	   ter = elev1(i,j)     ! departures at stn elev
	   dz = ter - h7(i,j)
c
	   if(t1(i,j) .ne. 0.) then
	      tts = t7(i,j) + (lapse_t * dz)
	      t1(i,j) = t1(i,j) - tts
	   endif
	   if(td1(i,j) .ne. 0.) then
	      ttds = td7(i,j) + (lapse_td * dz)
	      td1(i,j) = td1(i,j) - ttds
	   endif
c
	   ter = topo(i,j)     ! departures on laps topo
	   dz = ter - h7(i,j)
	   tt(i,j) = t7(i,j) + (lapse_t * dz)   
	   ttd(i,j) = td7(i,j) + (lapse_td * dz)	
	   if(tb81(i,j) .ne. 0.) then
	      tb81(i,j) = tb81(i,j) - tt(i,j)
	   endif
	   if(t_bk(i,j) .ne. 0.) then
	      t_bk(i,j) = t_bk(i,j) - tt(i,j)
	   endif
	   if(td_bk(i,j) .ne. 0.) then
	      td_bk(i,j) = td_bk(i,j) - ttd(i,j)
	   endif
c
	   dz = ter - hbar     ! calc a sfc pressure
	   tbar = a_t + (lapse_t * (hbar + (dz * .5)))
	   tbar = (tbar - 32.) * fon + 273.15 ! conv F to K
	   psfc(i,j) = pbar * exp(-dz * gor / tbar)
	enddo !i
	enddo !j
c
c.....	Now call the solution algorithm for the tb8 data.
c
	call zero(w,imax,jmax)
	print *,'  At spline call for tb8'
	name = 'TB8   '	
        call spline(tb8,tb81,w,alf,z,beta,0.,st,cormax,.3,imax,jmax,
     &        roi,bad_tb8,imiss,d1,d2,d3,name)	! z is a zero array here...
	if(imiss .ne. 0) ibt = 0 ! all zeros in tb8 array
c
c.....	Now force the t analysis with the tb8 and background data.
c
	name = 'NOPLOT'	
	gamma = 5.
	if(ibt .eq. 0) gamma = 0.		! data not there
	print *,'  At spline call for t'
        call spline(t,t1,t_bk,alf,wt_t,beta,gamma,tb8,cormax,.3,imax,
     &        jmax,roi,bad_t,imiss,d1,d2,d3,name)
c
c.....	Now call the solution algorithm for the dew point.
c
c	beta_td = 3.0
!	write(6,19999)
!19999	format(' Enter beta for dewpoints: ',$)
!	read(5,*) beta_td
!	write(6,19998) beta_td
19998	format(' Using beta_td of: ',f10.2)
	print *,'  At spline call for td'
        call spline(td,td1,td_bk,alf,wt_td,beta_td,0.,sm,cormax,.3,
     &        imax,jmax,roi,bad_td,imiss,d1,d2,d3,name)
c
c.....	Convert the analysed perturbations back to t, td, and tb8 (don't
c.....	bother with the t and td backgrounds since we're done with them).
c
	call add(t,tt,t,imax,jmax)
	call add(td,ttd,td,imax,jmax)
	call add(tb8,tt,tb8,imax,jmax)
c
c.....	Check to make sure that td is not greater than t...1st time.
c
	do j=1,jmax
	do i=1,imax
	  if(td(i,j) .gt. t(i,j)) td(i,j) = t(i,j)
	enddo !i
	enddo !j
c
c.....	Find an estimated terrain theta and theta-e while checking
c.....  the sfc theta vs an upper level theta.  If the sfc theta
c.....  is greater than the upper theta, set the sfc theta to the upper
c.....  theta, then calculate and store a corrected temperature.  For
c.....  the upper theta, pick a level that is above your highest data,
c.....  that takes into account mixing, etc.  For Colorado, use 500 mb;
c.....  lower elevations try 700 mb.
c
	itheta_all = 1
	if(itheta .eq. 0) then
	   print *,' Skipping theta check for surface temperatures.'
	   itheta_all = 0
	elseif(itheta .eq. 7) then
	   print *,' Checking sfc temperatures with 700 mb theta.'
	   if(itheta7 .eq. 0) then
	      print *,' 700 mb Temps missing. No theta check.'
	      itheta_all = 0
	   else
	      call zero(d1,imax,jmax)
	      call conv_f2k(t7,t7,imax,jmax)
	      ratio = (1000. / 700.) ** (.286)
	      do j=1,jmax
	      do i=1,imax
		 d1(i,j) = t7(i,j) * ratio
	      enddo !i
	      enddo !j
	   endif
	elseif(itheta .eq. 5) then
	   print *,' Checking sfc temperatures with 500 mb theta.'
	   if(itheta5 .eq. 0) then
	      print *,' 500 mb Temps missing. No theta check.'
	      itheta_all = 0
	   else
	      call zero(d1,imax,jmax)
	      ratio = (1000. / 500.) ** (.286)
	      do j=1,jmax
	      do i=1,imax
		 d1(i,j) = t5(i,j) * ratio
	      enddo !i
	      enddo !j
	   endif
	else
	   print *,' Bad itheta value. Skipping sfc temp check.'
	   itheta_all = 0
	endif
c
c
	icnt_th = 0
	do j=1,jmax
	do i=1,imax
	   torg_f = t(i,j)
	   tc = (t(i,j) - 32.) * fon            ! sfc T in F to C
	   theta_c = o(tc,psfc(i,j))            ! sfc Th in C
	   theta_old = theta_c
	   theta_k = theta_c + 273.15           ! sfc Th in K
	   if(itheta_all .ne. 0) then
	     if(theta_k .gt. d1(i,j)) then      ! if sfc Th > 500 Th...
		theta_k = d1(i,j)               ! set sfc Th = 500 Th
		theta_c = theta_k - 273.15      ! adj sfc Th in C
		tc = tda(theta_c,psfc(i,j))     ! adj sfc T in C
		tnew_f = (tc * anof) + 32.      ! replace sfc T in F
		t(i,j) = tnew_f
		icnt_th = icnt_th + 1
		dff_t = tnew_f - torg_f
	 write(6,2244) i,j,theta_old,theta_c,torg_f,tnew_f,dff_t
	     endif
	  endif
	  theta(i,j) = theta_c                  ! sfc Th in C
	  tdc = (td(i,j) - 32.) * fon           ! sfc Td in F to C
	  thetae(i,j) = oe(tc,tdc,psfc(i,j))	! in C
	enddo !i
	enddo !j
 2244	format(' Adjusting sfc TH at ',2i4,/,'  Old/New TH(C): ',
     &         2f10.2,/,'  Old/New Sfc Temp(F): ',2f10.2,
     &         '   Difference: ',f10.2,/)
	print *,' Changed sfc TH/temp at ',icnt_th,' points.'
c
c.....	Check again (since we changed t) to make sure that td is not 
c.....  greater than t, so thermo stuff won't blow up later.
c
	do j=1,jmax
	do i=1,imax
	  if(td(i,j) .gt. t(i,j)) td(i,j) = t(i,j)
	enddo !i
	enddo !j
c
c..... Call the solution algorithm for the rest of the fields.
c
	print *,'  At spline call for u'
	call spline(u,u1,u_bk,alf,wt_u,beta,0.,st,cormax,.1,imax,jmax,
     &        roi,bad_u,imiss,d1,d2,d3,name)
	print *,'  At spline call for v'
	call spline(v,v1,v_bk,alf,wt_v,beta,0.,st,cormax,.1,imax,jmax,
     &        roi,bad_v,imiss,d1,d2,d3,name)
	print *,'  At spline call for red_p'
	call spline(rp,rp1,rp_bk,alf,wt_rp,beta,0.,st,cormax,.1,imax,
     &        jmax,roi,bad_p,imiss,d1,d2,d3,name)
	print *,'  At spline call for msl p'
	call spline(mslp,mslp1,mslp_bk,alf,wt_mslp,beta,0.,st,cormax,
     &      .1,imax,jmax,roi,bad_p,imiss,d1,d2,d3,name)
	print *,'  At spline call for visibility'
	call spline(vis,vis1,vis_bk,alf,wt_vis,beta,0.,st,cormax,10.,
     &        imax,jmax,roi,bad_vis,imiss,d1,d2,d3,name)
c
c.....	If no background fields are available, skip over the variational
c.....	section.  Fields will be Barnes/splines, and derived values will be
c.....	calculated.  The fields may not be very good....
c
	if(ilaps_bk.eq.0 .and. irams_bk.eq.0) then
	  call move(rp,p_a,imax,jmax)
	  call multcon(p_a,100.,imax,jmax)	! conv mb to Pa
	  call conv_kt2ms(u,u_a,imax,jmax)	! conv kt to m/s and move array
	  call conv_kt2ms(v,v_a,imax,jmax)	! conv kt to m/s and move array
	  go to 500
	endif
c
cv....	This is the where the variational analysis stuff starts.
c.....	First compute the pressure change
c
!	call diff(rp,rp_bk,drp,imax,jmax)	! change that drove analysis
c
c.....	compute wind changes   
c
	do 100 n=1,npass
c
	  if(n .gt. 1) then
	    call move(u_a,u,imax,jmax)
	    call move(v_a,v,imax,jmax)
	  else
	    call multcon(rp,100.,imax,jmax) 	!convert pressure to Pa.
	    call move(rp,p_a,imax,jmax)
	    call conv_kt2ms(u,u,imax,jmax)	!convert winds from kt to m/sec
	    call conv_kt2ms(v,v,imax,jmax)	!   "      "     "   "     "
	    call conv_kt2ms(u_bk,u_bk,imax,jmax)!   "      "     "   "     "
	    call conv_kt2ms(v_bk,v_bk,imax,jmax)!   "      "     "   "     "
	  endif
c
c.....  calculate and save the wind changes. first ones are off the bkg.
c
cc	  if(n .eq. 1) then
cc	    call zero(d1,imax,jmax)
cc	    call zero(d2,imax,jmax)
cc	    do j=1,jmax
cc	    do i=1,imax
cc	      if(u1(i,j) .ne. 0.) d1(i,j) = u1(i,j) - u_bk(i,j)  ! diff = 
cc	      if(v1(i,j) .ne. 0.) d2(i,j) = v1(i,j) - v_bk(i,j)  !   ob - bk
cc	    enddo !i
cc	    enddo !j
cc	    nbpass = 1
cc	    kdim = 5 	! rad of infl of 0.01 
cc	    call barnes2(du,imax,jmax,d1,bad,kdim,nbpass,d3,d4)
cc	    call barnes2(dv,imax,jmax,d2,bad,kdim,nbpass,d3,d4) 
cc	  else
	    call diff(u,u_bk,du,imax,jmax)
	    call diff(v,v_bk,dv,imax,jmax)
cc	  endif
c
c.....	compute divergence change and vorticity
c
	  do j=2,jmax
	  do i=2,imax
	    ddiv(i,j) =((du(i,j-1) - du(i-1,j-1)) / dx(i,j) +
     &                     (dv(i-1,j) - dv(i-1,j-1)) / dy(i,j)) / dt
	    vort(i,j) = (v(i,j) - v(i-1,j)) / dx(i,j) -
     &                     (u(i,j) - u(i,j-1)) / dy(i,j)
	  enddo !i
	  enddo !j
c
c.....	compute the nonlinear and friction terms
c
          call get_directory('static',infile,len)
          infile=infile(1:len)//'drag_coef.dat'
          call s_len(infile,len)
	  open(51,file=infile(1:len),
     &         form='unformatted',status='old')
c	  open(51,file='../static/surface/drag_coef.dat',
c     &         form='unformatted',status='old')
	  read(51) akk
	  close(51)
	  call nonlin(nu,nv,u,v,u_bk,v_bk,imax,jmax,dx,dy)
	  call frict(fu,fv,u,v,u_bk,v_bk,imax,jmax,ak,akk)
c
c.....	compute forcing function
c
	  do j=2,jmax-1
	  do i=2,imax-1
	    ddiva = (ddiv(i,j) + ddiv(i+1,j+1) + 
     &                      ddiv(i,j+1) + ddiv(i+1,j)) * 0.25
	    anux = ( nu(i+1,j) - nu(i-1,j) +
     &         nu(i+1,j-1) - nu(i-1,j-1) ) / dx(i,j) 
	    anvy = ( nv(i,j+1) - nv(i,j-1) + 
     &         nv(i-1,j+1) - nv(i-1,j-1) ) / dy(i,j)
	    anonlinterm = (anux + anvy) * 0.25
	    frictterm = ((fu(i+1,j) - fu(i-1,j) + fu(i+1,j-1) - 
     &         fu(i-1,j-1))/ dx(i,j) + (fv(i,j+1) - fv(i,j-1) + 
     &         fv(i-1,j+1) - fv(i-1,j-1)) / dy(i,j)) * .25
	    f(i,j) = rho * (-ddiva + fo(i,j) * vort(i,j)) +
     &         a(i,j) * rp(i,j) - rho * anonlinterm -
     &         rho * frictterm
	    f(imax,j) = f(imax-1,j)
	    f(i,jmax) = f(i,jmax-1)
	    u_a(i,j) = u(i,j)
	    v_a(i,j) = v(i,j)
	  enddo !i
	  enddo !j
	  scale=0.
	  ci=0.
c
	  call leib(p_a,f,60,.1,imax,jmax,z,z,z,z,a,dx,dy,1.)
c
	  write(6,1200) filename,gam,del,dt
1200	  format(1x,'wind and pressure analysis for ',a9/1x,' with gam,
     &           del, and dt = ',3e12.4)
c
	  do j=2,jmax-1
	  do i=2,imax-1
	    dpdy = (p_a(i,j+1) - p_a(i,j)) / dy(i,j)
	    dpdx = (p_a(i+1,j) - p_a(i,j)) / dx(i,j)
	    dvdtnf= (dv(i,j)+dv(i,j+1)+dv(i-1,j+1)+dv(i-1,j))/4./dt
     &             +(nv(i,j)+nv(i,j+1)+nv(i-1,j+1)+nv(i-1,j))*.25
     &             +(fv(i,j)+fv(i,j+1)+fv(i-1,j+1)+fv(i-1,j))*.25
	    u_a(i,j) = (u(i,j) - del * fo(i,j)*(dvdtnf + dpdy / rho)) /
     &                 (1. + del * fo2(i,j))
	    dudtnf= (du(i,j)+du(i+1,j)+du(i+1,j-1)+du(i,j-1))/4./dt
     &             +(nu(i,j)+nu(i+1,j)+nu(i+1,j-1)+nu(i,j-1))*.25
     &             +(fu(i,j)+fu(i+1,j)+fu(i+1,j-1)+fu(i,j-1))*.25
            v_a(i,j) = (v(i,j) + del * fo(i,j)*(dudtnf + dpdx / rho)) /
     &                 (1. + del * fo2(i,j))
	  enddo !i
	  enddo !j
c
c.....	fill in boundaries of u_a and v_a for finite diff calcs.
c
	call bounds(u_a,imax,jmax)
	call bounds(v_a,imax,jmax)
c
100	continue	! end loop on npass
500	continue	! skip to here if cold starting or no backgrnd
c
c.....	Channel the winds around the terrain
c
       
	call get_directory('static',infile,len)
	infile=infile(1:len)//'pbl_top.dat'
	call s_len(infile,len)
	open(52,file=infile(1:len),
     &         form='unformatted',status='old')
c	open(52,file='../static/surface/pbl_top.dat',
c     &       form='unformatted',status='old')
	read(52) pbl_top
	close(52)
cc	call vortdiv(u_a,v_a,vort,div,imax,jmax,dx,dy)
cc	call channel(u_a,v_a,topo,imax,jmax,pbl_top,pblht,dx,dy,z,
cc   &               d1,d2,d3,d4,d5,d6,d7,div)
c
c.....	Calculate the vorticity and divergence.
c
	call vortdiv(u_a,v_a,vort,div,imax,jmax,dx,dy)
c
c.....	Compute w by integrating the surface winds over some pbl and
c.....	allow for terrain lift.
c
	do j=2,jmax-1
	do i=2,imax-1
	   dterdx = (topo(i,j)+topo(i,j-1)-topo(i-1,j)-topo(i-1,j-1)
     1                ) * .5 / dx(i,j)
	   dterdy = (topo(i,j)+topo(i-1,j)-topo(i-1,j-1)-topo(i,j-1)
     1                ) * .5 / dy(i,j)
	   ubar = (u_a(i-1,j-1) + u_a(i,j-1)) * .5
	   vbar = (v_a(i-1,j) + v_a(i-1,j-1)) * .5
	   dvh = ubar * dterdx + vbar * dterdy
	   w(i,j) = (dvh - div(i,j) * pblht) * 100.	! cm/sec
	enddo !i
	enddo !j
	call bounds(w,imax,jmax)
c
c.....	Now convert some stuff, then call the thermo routines.
c
	call multcon(p_a,0.01,imax,jmax)
	call zero(d1,ni,nj)
	sflag = 0.
!	sflag = 2.e6	! this means li_laps will read 500 T from a file
	call li_laps(t,td,psfc,d1,i4time,imax,jmax,li,sflag,istatus)
c
	call make_cssi(t,td,mslp,u_a,v_a,cssi,imax,jmax,d1,d2)
c
	call conv_f2k(t,t,imax,jmax)		! conv F to K
	call conv_f2k(td,td,imax,jmax)		! conv F to K
	call addcon(theta,273.15,theta,imax,jmax)	! C to K
	call meso_anl(u_a,v_a,psfc,t,td,theta,dx,dy,q,qcon,qadv,
     &                thadv,tadv,ni,nj)
c
c.....	Convert stuff not already converted to MKS units.
c
        call zero(d1,imax,jmax)
        call zero(d2,imax,jmax)
	call multcon(p_a,100.,imax,jmax)	! conv mb to Pa
	call multcon(psfc,100.,imax,jmax)	! conv mb to Pa
	call multcon(mslp,100.,imax,jmax)	! conv mb to Pa
	call multcon(w,.01,imax,jmax)		! conv cm/s to m/s
	call addcon(thetae,273.15,thetae,imax,jmax)	! C to  K
	call windspeed(u_a,v_a,spd,imax,jmax)	! calc windspeed
	call hum(t,td,rh,imax,jmax,d1,d2)	! calc rel hum.
	call vlog2vis(vis,vis,imax,jmax)	! conv log(vis) to vis-miles
c
c.....  Adjust visibility analysis.
c
        call enhance_vis(i4time,vis,rh,topo,imax,jmax,d1,d_c,
     &                                                kcloud,kdum)
	call conv_miles2m(vis,vis,imax,jmax)	! conv miles to meters
c
c.....  Call the Fire index routine.
c
	print *,' '
	print *,' Fire Wx index...'
	isnow = 1
	ismoist = 1
	call zero(fire,imax,jmax)
	call zero(d6,imax,jmax)
	call zero(d7,imax,jmax)
c
	i4time_tol = 7200
	i4time_f = i4time - 3600
	ilev = -1
	var_fire = 'LSM'                   ! soil moisture
	ext_f = 'lm1'
	call get_2d_field(i4time,i4time_tol,i4time_near,ext_f,var_fire,
     &                    ilev,d6,imax,jmax,istatus)
	if(istatus .ne. 1) then
	   print *,' Error getting soil moisture.'
	   ismoist = 0
	endif
c
	var_fire = 'SC '                   ! snow cover
	ext_f = 'lm2'
	call get_laps_2d(i4time_f,ext_f,var_fire,units_fire,com_fire,
     &                   imax,jmax,d7,istatus)
	if(istatus .ne. 1) then
	   print *,' Error getting snow cover.'
	   isnow = 0
	endif
c
	call lp_fire_danger(rh,t,spd,d6,d7,topo,ismoist,isnow,
     &                                            fire,istatus)
	print *,' Fire: ',ismoist, isnow, istatus
	print *,' '
c
c.....  Calculate Heat Index
c
	print *,' Heat Index...'
	call heat_index(t,rh,hi,imax,jmax)
	call conv_f2k(hi,hi,imax,jmax)		! conv F to K
c
c.....	Now calculate PBE and NBE.  First fix up the .LT1 file....
c
	print *,' PBE/NBE calcs...'
	no_3d_t = 1	! flag for if 3-d temps missing: 1-there,0-missing
	jstatus(4) = 1
        LT1_write_flag = 0     ! flag for writing the LT1 file: 1-yes,0-no
c
	call move(t,d1,imax,jmax)     ! move temps into dummy array
	call zero(d2,imax,jmax)
	call zero(d3,imax,jmax)
	call zero(d4,imax,jmax)
	call zero(d5,imax,jmax)
	call zero(d6,imax,jmax)
	call zero(d7,imax,jmax)
c
        call put_temp_anal(i4time,imax,jmax,kmax,imax_m,jmax_m,
     &                 dm1,dm2,dm3,
     &                 lni_m,lnj_m,iden_ratio,dm4,
     &                 adbox,ndbox,d2,d3,d4,d5,d6,d7,
     &                 ht_3d_m,LT1_out,
     &                 lat,lon,topo,d1,psfc,
     &                 LT1_write_flag,laps_cycle_time,grid_spacing,
     &                 rh_3d,t_3d_k,istatus)
c
	print *,' From put_temp_anal: istatus = ', istatus
	if(istatus .ne. 1) then
	   no_3d_t = 0	   
	   jstatus(4) = 2	! written but incompl. data
	endif
c
	if(no_3d_t .ne. 1) then		! skip p/nbe calc...no upper air temps
	  print *,' ** Upper temps missing.  Skipping BE calc. **'
	  call constant(pbe_2d,1.e6,imax,jmax)
	  call constant(nbe_2d,1.e6,imax,jmax)
	  go to 888
	endif
c
	do k = 1,nk ! Calculate pressure at each level
	    p_1d_pa(k) = pressure_of_level(k) ! Pressure at each level
	enddo !k
c
c	Calculate a 3-D Height Field
c	call get_heights_hydrostatic(t_3d_k,psfc,topo,d1,d2,d3,d4,
c     &                               ni,nj,nk,ht_3d_m)
c
c	Get PBE and NBE - Make sure t_sfc_k(i,j) >= td_sfc_k(i,j)
	call laps_be(ni,nj,nk,t,td,psfc,
     &  t_3d_k,ht_3d_m,p_1d_pa,topo,pbe_2d,nbe_2d)
c
 888	continue
c
c.....  Change RH to %
c
	call multcon(rh,100.,imax,jmax)
c
c.....  Write out some stats.
c
	print *,' ======================================================='
	print *,' u-wind (m/s):'
	call stats(u_a,imax,jmax)
	print *,' -------------------------------'
	print *,' v-wind (m/s):'
	call stats(v_a,imax,jmax)
	print *,' -------------------------------'
	print *,' wind spd (m/s):'
	call stats(spd,imax,jmax)
	print *,' -------------------------------'
	print *,' vert vel (m/s):'
	call stats(w,imax,jmax)
	print *,' -------------------------------'
	print *, redp_lvl,' m pressure (pa):'
	call stats(p_a,imax,jmax)
	print *,' -------------------------------'
	print *,' temp (k):'
	call stats(t,imax,jmax)
	print *,' -------------------------------'
	print *,' dewpt (k):'
	call stats(td,imax,jmax)
	print *,' -------------------------------'
	print *,' rel hum:'
	call stats(rh,imax,jmax)
	print *,' -------------------------------'
	print *,' theta (k):'
	call stats(theta,imax,jmax)
	print *,' -------------------------------'
	print *,' LI (K):'
	call stats(li,imax,jmax)
	print *,' -------------------------------'
	print *,' vis (m):'
	call stats(vis,imax,jmax)
	print *,' -------------------------------'
	print *,' fire :'
	call stats(fire,imax,jmax)
	print *,' -------------------------------'
	print *,' heat index :'
	call stats(hi,imax,jmax)
	print *,' ======================================================='

c.....	Now write out the grids to PROD_DEV.
c
	print *,' Saving primary fields.'
	do i=1,27
	  lvl(i) = 0
	  lvl_coord(i) = 'AGL'
	enddo !i
	lvl_coord(10) = 'MSL'
	lvl_coord(26) = 'MSL'
	var(1) = 'U'	! u-wind (m/s)
	var(2) = 'V'	! v-wind (m/s)
	var(3) = 'P'	! reduced press (Pa)
	var(4) = 'T'	! temp (K)
	var(5) = 'TD'	! dew point (K)
	var(6) = 'VV'	! vert. vel (m/s)
	var(7) = 'RH'	! relative humidity (%)
	var(8) = 'CCE'	! ceiling ht. msl (m)
	var(9) = 'MSL'	! MSL pressure (Pa)
	var(10) = 'TAD'	! temperature advection (K/sec)
	var(11) = 'TH'	! potential temp (K)
	var(12) = 'THE' ! equiv pot temp (K)
	var(13) = 'PS'	! surface press (Pa)
	var(14) = 'VOR'	! sfc vorticity (/s)
	var(15) = 'MR'	! mixing ratio (g/kg)
	var(16) = 'MRC'	! moisture convergence (g/kg/s)
	var(17) = 'DIV'	! sfc divergence (/s)
	var(18) = 'THA' ! pot temp adv (K/s)
	var(19) = 'MRA' ! moisture adv (g/kg/s)
	var(20) = 'LI'	! lifted index (K)
	var(21) = 'SPD'	! wind speed (m/s)
	var(22) = 'CSS' ! CSSI 
	var(23) = 'PBE' ! Pos Bouyant Energy (J/kg)
	var(24) = 'NBE' ! Neg Bouyant Energy (J/kg)
	var(25) = 'VIS' ! Visibility (m)
	var(26) = 'FWX' ! Fire threat index (integer)
	var(27) = 'HI'  ! Heat Index (K)
	units(1) = 'M/S'
	units(2) = 'M/S'
	units(3) = 'PA'
	units(4) = 'K'
	units(5) = 'K'
	units(6) = 'M/S'
	units(7) = '%'
	units(8) = 'M'
	units(9) = 'PA'
	units(10) = 'K/S'
	units(11) = 'K'
	units(12) = 'K'
	units(13) = 'PA'
	units(14) = '/S'
	units(15) = 'G/KG'
	units(16) = 'G/KG/S'
	units(17) = '/S'
	units(18) = 'K/S'
	units(19) = 'G/KG/S'
	units(20) = 'K'
	units(21) = 'M/S'
	units(23) = 'J/KG'
	units(24) = 'J/KG'
	units(25) = 'M'
	units(26) = ' '
	units(27) = 'K'
c
        do i=1,27
           write(comment(i)(50:52),180) n_sao_g
           write(comment(i)(54:56),180) n_sao_b
           write(comment(i)(58:60),180) n_obs_g
180        format(i3)
        enddo !i
	write(comment(3)(1:4),181) ifix(redp_lvl)
 181	format(i4)
	comment(3)(5:23) = ' M REDUCED PRESSURE'
	comment(9) = 'MSL PRESSURE'
	comment(22)= 'CSSI - COLORADO SEVERE STORM INDEX'
	comment(26)= 'LAPS FIRE THREAT INDEX'
	comment(26)(62:120) = 
     &      'INDEX: 0-NONE, 5-SLGT, 10-MDT, 15-HI, 20-EXTREME'
	comment(27)= 'HEAT INDEX'
c

	call write_laps_data(i4time,dir,ext,imax,jmax,27,27,var,
     &                  lvl,lvl_coord,units,comment,data,istatus)
c
	jstatus(3) = 1		! everything ok...
c
c.....  Now finish up with some verification.
c
	iunit = 11
        call get_directory('log',ver_file,len)
	ver_file = ver_file(1:len)//'qc/laps_sfc.ver.'//filename(6:9)
        call s_len(ver_file,len)
	call zero(d1,imax,jmax)

	open(iunit,file=ver_file(1:len),status='unknown')
c
	title = 'Temperature'
	call conv_k2f(t,d1,imax,jmax)
	call verify(d1,t_s,stn,n_obs_b,title,iunit)
c	
	title = 'Dew Point'
	call conv_k2f(td,d1,imax,jmax)
	call verify(d1,td_s,stn,n_obs_b,title,iunit)
c
	title = 'Wind Speed'
	call conv_ms2kt(spd,d1,imax,jmax)
	call verify(d1,ff_s,stn,n_obs_b,title,iunit)
c
	close(iunit)
c
c.....  That's it.  Let go home.
c
	print *,' Normal completion of LAPSVANL'
	return
c
	end

