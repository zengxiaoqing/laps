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
	subroutine laps_vanl(i4time,filename,ni,nj,nk,mxstn,
     &     itheta_in,redp_lvl,sfc_nl_parms,
     &     laps_cycle_time,dt,del,gam,ak,lat,lon,topo,ldf,grid_spacing, 
     &     lat_s,lon_s,elev_s,t_s,td_s,dd_s,ff_s,pstn_s,
     &     mslp_s,pred_s,vis_s,stn,n_obs_b,n_sao_b,n_sao_g,obs,
     &     u_bk,v_bk,t_bk,td_bk,rp_bk,mslp_bk,sp_bk,vis_bk,tgd_bk_f,
     &     wt_u, wt_v, wt_t, wt_td, wt_rp, wt_mslp, wt_vis, ilaps_bk, 
     &     back_t,back_td,back_uv,back_sp,back_rp,back_mp,back_vis,
     &     bad_t,bad_td,bad_u,bad_v,bad_p,bad_mp,bad_th,bad_the,
     &     bad_vis,bad_tb8,
     &     u1, v1, rp1, t1_f, td1_f, sp1, tb81, mslp1, vis1, elev1,
     &     x1a,x2a,y2a,ii,jj,jstatus)
c
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
c	Artificially varying this parameter can modulate the effects of
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
c                               08-07-97  Dynamic changes..rm equivs.
c                               08-27-97  Changes for dynamic LAPS.
c                               11-06-97  New puttmp call;rm some work arrays.
c                               05-12-98  Added MSL press verification.
c                               07-11-98  Added NaN checks on output.
c                               09-25-98  Add bkg flags for each bkg variable.
c	                        10-07-98  Range ck for pressures...until can
c	                                    fix better.
c                               12-02-98  Remove slot for ceiling in output;
c                                           reduce num_var by 1 (to 26). Also
c                                           remove vis verification, added 
c                                           verif of T, Td, P backgrounds.
c                               01-14-99  Remove spline for tb8 (sat data).
c                                           Add vis field ck. 
c                               01-19-99  Temp. turn off all splines, replace
c                                           with Barnes.
c                               01-28-99  Rm Barnes calls. Add check for t,td
c                                           bkgs, if none use tt, ttd.
c                               06-18-99  Remove LI, CAPE, CIN calcs.
c                               07-08-99  Remove LI/CAPE/CIN arrays.  Change 
c                                           stn char array.  Rm *4 from all
c                                           declarations. Don't calc subsitute
c                                           T/Td backgrounds here. Don't use
c                                           sat data in T spline if no bkg.
c                               08-13-99  Change spline call to allow diff wts
c                                           for diff variables.  Fix for SGI..
c                                           constants passed to subroutines.
c                               09-19-99  Turn off range check for p_a on output.
c                               11-23-99  Change pbar calc...check for a bad one.
c                                           Fix error with grid spacing.
c
c*****************************************************************************
c 
        use mem_namelist, ONLY: iwrite_output, rms_pres

        use mem_sfcanl, ONLY: alloc_sfcanl_arrays, point_sfcanl_arrays       

        use mem_sfcanl, ONLY: u_a,v_a,p_a,t,td,vv,rh,hi,mslp,tadv
     +                       ,theta,thetae,psfc,vort,q,qcon,div,thadv
     +                       ,qadv,spd,cssi,vis,fire,tgd_k    

	include 'laps_sfc.inc'
        include 'laps_cloud.inc'

	parameter(                !Expected observation error, ea. var.
     &            obs_error_redp  = 0.1,  ! for reduced pressure
     &            obs_error_t     = 0.1,  ! for temperature
     &            obs_error_tb8   = 0.1,  ! for Brighness temperatures
     &            obs_error_td    = 0.1,  ! for dew point
     &            obs_error_mslp  = 0.1,  ! for MSL pressure
     &            obs_error_wind  = 0.1,  ! for wind
     &            obs_error_vis   = 0.1)  ! for visibility
c
	parameter(bad = 1.e6 - 2.)	! larger than 'bad' are.

c
c
c.....	LAPS lat/lon and terrain grids, and Coriolis.
c
	integer istatus
        real grid_spacing
	real lat(ni,nj), lon(ni,nj), topo(ni,nj), ldf(ni,nj)
	real fo(ni,nj), fo2(ni,nj), akk(ni,nj), pbl_top(ni,nj)
c
c.....  Stuff for intermediate grids (old LGS file)
c
	real u1(ni,nj), v1(ni,nj)
	real t1_f(ni,nj), td1_f(ni,nj), tb81(ni,nj)                 ! Deg F
	real rp1(ni,nj), sp1(ni,nj), mslp1(ni,nj)
	real vis1(ni,nj), elev1(ni,nj)
c
c.....	Grids for the first data's analyses.
c
 	real u(ni,nj), v(ni,nj)
!       real t(ni,nj), td(ni,nj)
!       real theta(ni,nj), thetae(ni,nj), psfc(ni,nj), vis(ni,nj)
	real rp(ni,nj)
        real tb8(ni,nj)
!       real mslp(ni,nj), tgd_k(ni,nj)
c
c.....	Grids for the variational analyses of rp, u, v (grid north)
c
!       real p_a(ni,nj), u_a(ni,nj), v_a(ni,nj) ! Post Variational Pa & M/S
	real p_a_orig(ni,nj), u_a_orig(ni,nj), v_a_orig(ni,nj)   ! Reference
c                                                                ! Pa & M/S
c.....	Grids for the derived quantities.
c
	real du(ni,nj), dv(ni,nj)
	real drp(ni,nj)
!       real vv(ni,nj), spd(ni,nj)
	real tt(ni,nj), ttd(ni,nj)
!       real qadv(ni,nj)
!       real rh(ni,nj), hi(ni,nj)
!	real cssi(ni,nj), fire(ni,nj)
	real p_1d_pa(nk)
c
c.....	Grids for variables derived by the MESO_ANL subroutine.
c
!	real q(ni,nj), qcon(ni,nj), thadv(ni,nj), tadv(ni,nj)
c
c.....	Grids for the background fields.
c
        real u_bk(ni,nj), v_bk(ni,nj), t_bk(ni,nj), td_bk(ni,nj)  ! kt
        real wt_u(ni,nj), wt_v(ni,nj), wt_t(ni,nj), wt_td(ni,nj)
        real rp_bk(ni,nj), mslp_bk(ni,nj), sp_bk(ni,nj)
        real wt_rp(ni,nj), wt_mslp(ni,nj)
        real vis_bk(ni,nj), wt_vis(ni,nj)
	real tgd_bk_f(ni,nj)
	real tb8_bk(ni,nj)
	real dum1(ni,nj), dum2(ni,nj)
        real u_bk_ms(ni,nj), v_bk_ms(ni,nj)                       ! m/s
        integer back_t, back_td, back_rp, back_uv, back_vis, back_sp
        integer back_mp
        real wt_bkg_a(ni,nj)                         
c
c.....	Grids for other stuff.
c
        real fnorm(0:ni-1,0:nj-1)
	real ddiv(ni,nj)
!       real vort(ni,nj), div(ni,nj)
	real f(ni,nj), fu(ni,nj), fv(ni,nj)
	real a(ni,nj), z(ni,nj), dx(ni,nj), dy(ni,nj)
	real nu(ni,nj),nv(ni,nj), h7(ni,nj)
	real t5(ni,nj), t7(ni,nj), td7(ni,nj)                     ! Deg K
	real pres_3d(ni,nj,nk)
c
c..... Stuff for the sfc data and other station info (LSO +)
c
        include 'sfcob.inc'
        type (sfcob) obs(mxstn)

	real lat_s(mxstn), lon_s(mxstn), elev_s(mxstn)
	real t_s(mxstn), td_s(mxstn)
        real dd_s(mxstn), ff_s(mxstn) ! true north (knots)
	real pstn_s(mxstn), mslp_s(mxstn), pred_s(mxstn), vis_s(mxstn)
        real ob_full(mxstn)
        real u_s(mxstn), v_s(mxstn) ! grid north (knots)
c
	character stn(mxstn)*20
        character title*60, ver_file*256
c
c.....	dummy work arrays
c
	real d1(ni,nj), d2(ni,nj)
	real dm1(ni,nj,nk)
	real dums(mxstn)
c
	real lapse_t, lapse_td
	real make_td
	real x1a(ni), x2a(nj), y2a(ni,nj)
	integer ii(mxstn), jj(mxstn)
	integer jstatus(20)
	character name*10, filename*9, infile*256
	character var_fire*3, com_fire*125, units_fire*10, ext_f*31
	character var_lga*3, ext_lga*31
c
c
c.....	Stuff for LAPS outputs (i.e., standard forms).
c
	parameter(num_var = 24)
	real data(ni,nj,num_var)
	integer imax,jmax,lvl(num_var)
	character dir*256,ext*31,var(num_var)*3,lvl_coord(num_var)*4
	character units(num_var)*10, comment(num_var)*125
c
c
c.....	Start...set up constants, initialize arrays, etc.
c
        call alloc_sfcanl_arrays(ni,nj)
        call point_sfcanl_arrays()

        I4_elapsed = ishow_timer()

	call tagit('laps_vanl', 19991123)
        rms_thresh_norm = 1.0    ! used for barnes_multivariate
	zcon = 0.
	ibt = 1      !assume have sat data...code cks later.
	jstatus(3) = -1		 ! start w/this until changed
	pi = 4. * atan(1.)
	imax = ni
	jmax = nj
	kmax = nk
	itmax = 60               !max num interations for leib
	erf = 0.1                !error in leib
	err = .0003
	ovr = 1.4
	scale = 0.
	npass = 1
	rho = 1.25
	rho2 = rho * rho
	omg2 = 2. * 7.292e-5
	rdpdg = pi / 180.
	re = 6371222.
	gor = 9.808 / 287.04
	fon = 5. / 9.
	anof = 9. / 5.
	pblht = 500.		! pbl height in meters
	fill_val = 1.e37

        call get_r_missing_data(r_missing_data,istatus)
c
	call zero(z,imax,jmax)
	call zero(p_a,imax,jmax)
	call zero(u_a,imax,jmax)
	call zero(v_a,imax,jmax)
	call zero(nu,imax,jmax)
	call zero(nv,imax,jmax)
	call zero(fu,imax,jmax)
	call zero(fv,imax,jmax)
	call zero(ddiv,imax,jmax)
	call zero(vort,imax,jmax)
	call zero(f,imax,jmax)
	call zero(du,imax,jmax)
	call zero(dv,imax,jmax)

        weight_bkg_const = 5e28 ! a la wind.nl

        wt_bkg_a = weight_bkg_const
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
	do j=1,jmax
	do i=1,imax
	   call get_grid_spacing_actual_xy(lat(i,j),lon(i,j)
     &             ,grid_spacing_actual_mx
     &             ,grid_spacing_actual_my
     &             ,istatus)
	   dy(i,j) = grid_spacing_actual_my
	   dx(i,j) = grid_spacing_actual_mx
           if(del .gt. 0.)then ! we will not divide by zero
	       a(i,j) = - gam * rho2 / del
           endif
	enddo !i
	enddo !j
c
cz..... Compute T on the surface using the LGA (or equiv) 700 T and HT.
c
	i4time_tol = 21600
	ext_lga = 'lga'

c       Determine 3D grid levels closest to 700mb and 500mb. This makes the
c       approximation that the 3D grid is on a constant pressure grid.

        icen = ni/2
        jcen = nj/2

        call get_pres_3d(i4time,ni,nj,nk,pres_3d,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error: Bad status returned from get_pres_3d'       
            return
        endif

        arg = rlevel_of_field(70000.
     1                        ,pres_3d,ni,nj,nk,icen,jcen,istatus)
        if(istatus .ne. 1)return
        k_700 = nint(arg)

        arg = rlevel_of_field(50000.
     1                        ,pres_3d,ni,nj,nk,icen,jcen,istatus)
        if(istatus .ne. 1)return
        k_500 = nint(arg)

c
c.....  Get the latest 3d fields, pull out the var/lvls needed.
c
	print *,' Get LGA 700 T, level ',k_700
	itheta7 = 1
	var_lga = 'T3 '
	call get_modelfg_3d(i4time,var_lga,ni,nj,nk,dm1,istatus)
c
	if(istatus .ne. 1)  then
	   print *,' LGA 700 T not available. Using constant 5C.'
	   call constant(t7,278.15,imax,jmax)
	   itheta7 = 0
	else
	   call move_3dto2d(dm1,k_700,t7,ni,nj,nk)  
	endif
c
	print *,' Get LGA 700 HT, level ',k_700
	var_lga = 'HT '
	call get_modelfg_3d(i4time,var_lga,ni,nj,nk,dm1,istatus)
c
	if(istatus .ne. 1)  then
	   print *,' LGA 700 HT not available. Using constant 3000 m.'
	   call constant(h7,3000.,imax,jmax)
	else
	   call move_3dto2d(dm1,k_700,h7,ni,nj,nk)  ! lvl 9 = 700 hPa
	endif
c
	print *,' Get LGA 700 TD, level ',k_700
	var_lga = 'SH '  ! specific humidity 
	call get_modelfg_3d(i4time,var_lga,ni,nj,nk,dm1,istatus)
c
	if(istatus .ne. 1)  then
	   print *,' LGA 700 Td not available. Using constant -5C.'
	   call constant(td7,268.15,imax,jmax) 
	else
	   do j=1,nj
	   do i=1,ni
	    t7_c = t7(i,j) - 273.15  !K to C
	    qgkg = dm1(i,j,k_700) * 1000.   
	    td7(i,j) = make_td(700., t7_c, qgkg, 0.) + 273.15   ! in K
	   enddo !i
	   enddo !j
	endif
c
c.....  Get the 500 Temps while we're here
c
	print *,' Get LGA 500 T, level ',k_500
	itheta5 = 1
	var_lga = 'T3 '
	call get_modelfg_3d(i4time,var_lga,ni,nj,nk,dm1,istatus)
c
	if(istatus .ne. 1)  then
 	   print *,' LGA 500 T not available.'
	   call constant(t5,badflag,imax,jmax) 
	   itheta5 = 0
	else
	   call move_3dto2d(dm1,k_500,t5,ni,nj,nk)  
	endif
c
c.....  Get lapse rate (usually std), and mean pressure.
c
	print *,' '
	call mean_lapse(n_obs_b,elev_s,t_s,td_s,a_t,lapse_t,a_td,
     &                  lapse_td,hbar,badflag)
	print *,' '
cc	call mean_pres(n_obs_b,pstn_s,pbar)
	call mean_pressure(pstn_s,n_obs_b,sp_bk,imax,jmax,badflag,pbar)
	if(pbar .le. 0.) then
	   print *,'  ERROR. Mean pressure is: ', pbar
	   print *,'    Setting pbar to 925 mb so we can continue...',
     &             'the analyses may not be any good.'
	   print *,'    Check your pressure obs and backgrounds.'
	   pbar = 925.
	else
	   print *,'  Mean pressure from the obs or bkg is: ', pbar
	endif

	do j=1,jmax
	do i=1,imax
c
c.....     Correct gridded obs (t1_f,td1_f) using lapse rate and difference
c          between gridded stn elevations (elev1) and laps terrain (topo)

!          Here it seems to be important that the stations are mapped onto the
!          grid with rounding up allowed to get the best possible departures

	   ter_s = elev1(i,j)     ! departures at stn elev
c
	   dz = ter_s - topo(i,j)
	   if(t1_f(i,j) .ne. 0.) then
	      tts = (lapse_t * dz)
	      t1_f(i,j) = t1_f(i,j) - tts
	   endif
	   if(td1_f(i,j) .ne. 0.) then
	      ttds = (lapse_td * dz)
	      td1_f(i,j) = td1_f(i,j) - ttds
	   endif
c
c.....    Departures of other stuff on laps topo
c
	   ter_g = topo(i,j)     ! departures on laps topo
c
	   dz = ter_g - hbar     ! calc a sfc pressure
	   tbar = a_t + (lapse_t * (hbar + (dz * .5)))
	   tbar = (tbar - 32.) * fon + 273.15        ! conv F to K
	   psfc(i,j) = pbar * exp(-dz * gor / tbar)  ! Has no meteorological 
                                                     ! variation       

!          Use LGA 700mb t/td and lapse rates to estimate sfc t/td
	   dz = 0. ! ter_g - ter_g
	   tt(i,j) = 0. ! (lapse_t * dz)   
	   ttd(i,j) = 0. ! (lapse_td * dz)	
c
!          if(tb81(i,j) .ne. 0.) then
!              tb81(i,j) = tb81(i,j) - tt(i,j)
!          endif

           if(sp_bk(i,j) .ne. 0. .and. back_sp .eq. 1)then ! psfc from bkgnd  
                                                           ! can be used as it
                                                           ! is on the LAPS trn
               psfc(i,j) = sp_bk(i,j)
           endif

	enddo !i
	enddo !j

        I4_elapsed = ishow_timer()
c
cc	do j=1,jmax
c	   do i=1,imax
c	      if(elev1(i,j) .ne. 0.) then
c		 write(6,7119) i,j,elev1(i,j)
c	      endif
c	   enddo
c	enddo
c
c.....	Now call the solution algorithm for the tb8 data.
c
c	print *,'  Fill in tb8 field using smooth Barnes'
c
	n_obs_var = 0

c       fill_val = 1.e37
c       smsng = 1.e37
c	npass = 1
c	rom2 = 0.005
c
	call zero(tb8, imax,jmax)
	if(back_t .ne. 1) call zero(tb81, imax,jmax) !if no bkg, 
                                                     ! don't use sat data
c
c	call dynamic_wts(imax,jmax,n_obs_var,rom2,d,fnorm)
c	call barnes2(tb8,imax,jmax,tb81,smsng,mxstn,npass,fnorm)
c       print *,' Got the weights'
c	call check_field_2d(tb8,imax,jmax,fill_val,istatus)
c	if(istatus .eq. 0) ibt = 0       !empty field
c
c	name = 'TB8   '	
c       call spline(tb8,tb81,tb8_bk,alf,z,beta,0.,z,cormax,.3,imax,jmax,
c     &        rms_thresh_norm,bad_tb8,imiss,mxstn,obs_error_tb8,name)       
c	if(imiss .ne. 0) ibt = 0 ! all zeros in tb8 array
c
c.....	Now force the t analysis with the tb8 and background data.
c
	name = 'NOPLOT'	
	bad_tm = bad_t
	if(back_t .ne. 1) bad_tm = bad_t * 2.
	print *,' '
	print *,'  At spline call for t (F)'

 7119	format(2i5,f10.2)

        if(.true.)then
	    alf = 10000.
	    alf2a = 0.
	    beta = 100.
	    gamma = 5.
	    if(ibt .eq. 0) gamma = 0.
            call spline(t,t1_f,t_bk,alf,alf2a,beta,gamma,tb81,cormax,
     &        err,imax,jmax,sfc_nl_parms%rms_temp,
     &        bad_tm,imiss,mxstn,obs_error_t,name,topo,ldf,wt_bkg_a)
        else ! use data structures for handling obs
            bad_tm_land  = bad_tm
            bad_tm_water = bad_tm
            call barnes_multivariate_sfc_jacket('t',obs,mxstn,t_bk
     1                                     ,badflag,imax,jmax
     1                                     ,rms_thresh_norm
     1                                     ,bad_tm_land
     1                                     ,bad_tm_water
     1                                     ,topo,ldf,wt_bkg_a,t,istatus)       
        endif

	print *,' '
	print *,'  At spline call for td (F)'
	bad_tmd = bad_td
	if(back_t .ne. 1) bad_tmd = bad_td * 2.

        if(.true.)then
	    alf = 10000.
	    alf2a = 0.
	    beta = 100.  
            call spline(td,td1_f,td_bk,alf,alf2a,beta,zcon,z,cormax,
     &        err,imax,jmax,sfc_nl_parms%rms_dewpoint,
     &        bad_tmd,imiss,mxstn,obs_error_td,name,topo,ldf,wt_bkg_a)
        else ! use data structures for handling obs
            bad_tmd_land  = bad_tmd
            bad_tmd_water = bad_tmd
            call barnes_multivariate_sfc_jacket('td',obs,mxstn,td_bk
     1                                    ,badflag,imax,jmax
     1                                    ,rms_thresh_norm
     1                                    ,bad_tmd_land
     1                                    ,bad_tmd_water
     1                                    ,topo,ldf,wt_bkg_a,td,istatus)       
        endif
c
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
	print *,' '

        if(itheta_in .eq. -1)then
           if(topo(imax/2,jmax/2) .ge. 1000.)then
              itheta = 5
           else
              itheta = 7
           endif
	   print *,' itheta auto-set based on terrain = ',itheta
        else
           itheta = itheta_in
        endif

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
!             call conv_f2k(t7,t7,imax,jmax)
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
	dff_tmx = -9.e30
	dff_tmn =  9.e30
	do j=1,jmax
	do i=1,imax
	   torg_f = t(i,j)
	   t_c = (t(i,j) - 32.) * fon           ! sfc T in F to C
	   theta_c = o(t_c,psfc(i,j))           ! sfc Th in C
	   theta_old = theta_c
	   theta_k = theta_c + 273.15           ! sfc Th in K
	   if(itheta_all .ne. 0) then
	     if(theta_k .gt. d1(i,j)) then      ! if sfc Th > Upper Th...
		theta_k = d1(i,j)               ! set sfc Th = Upper Th
		theta_c = theta_k - 273.15      ! adj sfc Th in C
		t_c = tda(theta_c,psfc(i,j))    ! adj sfc T in C
		tnew_f = (t_c * anof) + 32.     ! replace sfc T in F
		t(i,j) = tnew_f
		icnt_th = icnt_th + 1
		dff_t = tnew_f - torg_f
		if(dff_t .gt. dff_tmx) then
		   dff_tmx = dff_t
		   i_mx = i
		   j_mx = j
		endif
		if(dff_t .lt. dff_tmn) then
		   dff_tmn = dff_t
		   i_mn = i
		   j_mn = j
		endif
!	 write(6,2244) i,j,theta_old,theta_c,torg_f,tnew_f,dff_t
	     endif
	  endif
	  theta(i,j) = theta_c                  ! sfc Th in C
	  td_c = (td(i,j) - 32.) * fon          ! sfc Td in F to C
	  thetae(i,j) = oe(t_c,td_c,psfc(i,j))	! in C
	enddo !i
	enddo !j
 2244	format(' Adjusting sfc TH at ',2i4,/,'  Old/New TH(C): ',
     &         2f10.2,/,'  Old/New Sfc Temp(F): ',2f10.2,
     &         '   Difference: ',f10.2,/)
	print *,' '
	print *,' Changed sfc TH/temp at ',icnt_th,' points.'
	if(icnt_th .gt. 0) then
	   print *,'   Max change at ',i_mx,',',j_mx,': ',dff_tmx
	   print *,'   Min change at ',i_mn,',',j_mn,': ',dff_tmn
	endif
	print *,' '
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
c..... Note: 'gamma' (satellite weight) is zero for these.
c
	print *,' '
	print *,'  At spline call for u (kt)'
	bad_uw = bad_u
	if(back_uv .ne. 1) bad_uw = bad_u * 2.
	alf = 10000.
	alf2a = 0.
	beta = 100.
	call spline(u,u1,u_bk,alf,alf2a,beta,zcon,z,cormax,err,imax,jmax,
     &              sfc_nl_parms%rms_wind,bad_uw,imiss,
     &              mxstn,obs_error_wind,name,topo,ldf,wt_bkg_a)       
c
	print *,' '
	print *,'  At spline call for v (kt)'
	bad_vw = bad_v
	if(back_uv .ne. 1) bad_vw = bad_v * 2.
	alf = 10000.
	alf2a = 0.
	beta = 100.
	call spline(v,v1,v_bk,alf,alf2a,beta,zcon,z,cormax,err,imax,jmax,
     &              sfc_nl_parms%rms_wind,bad_vw,imiss,
     &              mxstn,obs_error_wind,name,topo,ldf,wt_bkg_a)        
c
	print *,' '
	print *,'  At spline call for red_p (mb)'
	bad_rp = bad_p
	if(back_rp .ne. 1) bad_rp = bad_p * 2.
	alf = 10000.
	alf2a = 0.
	beta = 100.
C
C TH: 29 November 2002 Begin hack.
C We set 'name' to be 'PRESSURE' so the spline routine knows how to set
C the mask_sea flag.
C
        name = 'PRESSURE'
	call spline(rp,rp1,rp_bk,alf,alf2a,beta,zcon,z,cormax,err,imax,
     &        jmax,rms_thresh_norm,bad_rp,imiss,mxstn,obs_error_redp,
     &        name,topo,ldf,wt_bkg_a)     
	name = 'NOPLOT'	
c
	print *,' '
	print *,'  At spline call for msl p (mb)'
cc	if(back_mp .ne. 1) bad_mp = bad_p * 2.
	alf = 10000.
	alf2a = 0.
	beta = 100.
        name = 'PRESSURE'
	call spline(mslp,mslp1,mslp_bk,alf,alf2a,beta,zcon,z,cormax,
     &      err,imax,jmax,rms_pres,bad_mp,imiss,mxstn,
     &      obs_error_mslp,name,topo,ldf,wt_bkg_a)
C
C TH: End hack.
C
c
!       Call routine to check pres arrays and adjust psfc 
        if(.false.)then ! adjust based on mslp/mslp_bk
            call pstn_anal(back_mp,back_sp,mslp_bk,mslp,imax,jmax
     1                    ,sp_bk,psfc)       
        else           ! adjust based on rp/rp_bk
            call pstn_anal(back_rp,back_sp,rp_bk,rp,imax,jmax
     1                    ,sp_bk,psfc)       
        endif
	name = 'NOPLOT'	

	print *,' '
	print *,'  At spline call for visibility (log)'
	bad_vs = bad_vis
	if(back_vis .ne. 1) bad_vs = bad_vis * 2.
	alf = 10000.
	alf2a = 0.
	beta = 100.
	call spline(vis,vis1,vis_bk,alf,alf2a,beta,zcon,z,cormax,err,
     &        imax,jmax,rms_thresh_norm,bad_vs,imiss,mxstn,
     &        obs_error_vis,name,topo,ldf,wt_bkg_a)

        write(6,*)' Analyze TGD observations'
        call barnes_multivariate_sfc_jacket('tgd',obs,mxstn
     1                                 ,tgd_bk_f
     1                                 ,badflag,imax,jmax
     1                                 ,rms_thresh_norm
     1                                 ,sfc_nl_parms%bad_tgd_land
     1                                 ,sfc_nl_parms%bad_tgd_water
     1                                 ,topo,ldf
     1                                 ,wt_bkg_a
     1                                 ,d2,istatus)
	call conv_f2k(d2,tgd_k,imax,jmax)                  ! conv F to K
c
c.....	If no background fields are available, skip over the variational
c.....	section.  Fields will be Barnes/splines, and derived values will be
c.....	calculated.  The fields may not be very good....
c
	if(ilaps_bk .eq. 0 .or. del .eq. 0.) then
	  call move(rp,p_a,imax,jmax)
	  call multcon(p_a,100.,imax,jmax)	! conv mb to Pa
	  call conv_kt2ms(u,u_a,imax,jmax)	! conv kt to m/s and move array
	  call conv_kt2ms(v,v_a,imax,jmax)	! conv kt to m/s and move array
          write(6,*)' Skipping variational adjustment of u,v,p '
     1             ,ilaps_bk,del
	  go to 500
	endif

cv....	This is the where the variational analysis stuff starts.
c
c.....	Compute and save the wind changes.
c
	do 100 n=1,npass
c
	  if(n .gt. 1) then
	    call move(u_a,u,imax,jmax)
	    call move(v_a,v,imax,jmax)
	  else
	    call multcon(rp,100.,imax,jmax) 	!convert pressure to Pa.
	    call move(rp,p_a,imax,jmax)
	    call conv_kt2ms(u,u,imax,jmax)	   ! convert winds (kt -> m/s)
	    call conv_kt2ms(v,v,imax,jmax)	   !    "      "    "      "  
	    call conv_kt2ms(u_bk,u_bk_ms,imax,jmax)!    "      "    "      "  
	    call conv_kt2ms(v_bk,v_bk_ms,imax,jmax)!    "      "    "      " 

!           Save the fields prior to the variational step for reference
            u_a_orig = u      ! m/s
            v_a_orig = v      ! m/s
            p_a_orig = p_a    ! Pascals

	  endif
c
	  call diff(u,u_bk_ms,du,imax,jmax)
	  call diff(v,v_bk_ms,dv,imax,jmax)
c
c.....	Compute divergence change and vorticity
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
	  call get_directory('static', infile, len)
	  infile = infile(1:len) // '/drag_coef.dat'
	  call s_len(infile, len)
c	  open(51,file='../static/drag_coef.dat',
	  open(51,file=infile(1:len),
     &         form='unformatted',status='old')
	  read(51) akk
	  close(51)
          ro=.667  ! ro is V/fL where L is ave data distance *4
	  call nonlin_2d(nu,nv,u,v,u_bk_ms,v_bk_ms,imax,jmax,dx,dy)
	  call frict_2d(fu,fv,u,v,u_bk_ms,v_bk_ms,imax,jmax,ak,akk)
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
     &         a(i,j) * rp(i,j) - rho * ro* anonlinterm + 
     &         rho * frictterm
	    f(imax,j) = f(imax-1,j)
	    f(i,jmax) = f(i,jmax-1)
	    u_a(i,j) = u(i,j)
	    v_a(i,j) = v(i,j)
	  enddo !i
	  enddo !j
	  scale=0.
c
	  call zero(z, ni,nj)
	  call leib_2d(p_a,f,itmax,erf,imax,jmax,z,z,z,z,a,dx,dy)
c
	  write(6,1200) filename,ro,gam,del,dt
1200	  format(1x,'wind and pressure analysis for ',a9/1x,' with ro, 
     &      gam,  del, and dt = ',4e12.4)
c
	  do j=2,jmax-1
	  do i=2,imax-1
	    dpdy = (p_a(i,j+1) - p_a(i,j)) / dy(i,j)
	    dpdx = (p_a(i+1,j) - p_a(i,j)) / dx(i,j)
	    dvdtnf= ro*(dv(i,j)+dv(i,j+1)+dv(i-1,j+1)+dv(i-1,j))/4./dt
     &             +ro*(nv(i,j)+nv(i,j+1)+nv(i-1,j+1)+nv(i-1,j))*.25
     &             -(fv(i,j)+fv(i,j+1)+fv(i-1,j+1)+fv(i-1,j))*.25
	    u_a(i,j) = (u(i,j) - del * fo(i,j)*(dvdtnf + dpdy / rho)) /
     &                 (1. + del * fo2(i,j))
	    dudtnf=ro* (du(i,j)+du(i+1,j)+du(i+1,j-1)+du(i,j-1))/4./dt
     &             +ro*(nu(i,j)+nu(i+1,j)+nu(i+1,j-1)+nu(i,j-1))*.25
     &             -(fu(i,j)+fu(i+1,j)+fu(i+1,j-1)+fu(i,j-1))*.25
            v_a(i,j) = (v(i,j) + del * fo(i,j)*(dudtnf + dpdx / rho)) /
     &                 (1. + del * fo2(i,j))
	  enddo !i
	  enddo !j
c
c.....	Fill in boundaries of u_a and v_a for finite diff calcs.
c
	call bounds(u_a,imax,jmax)
	call bounds(v_a,imax,jmax)
c
100	continue	! end loop on npass

!       Compare before and after adjustments
        call diff(u_a,u_a_orig,d1,imax,jmax)                   ! m/s
        write(6,*)
        write(6,*)' Stats for variational adjustment on U (m/s)'
        call stats(d1,imax,jmax)

        call diff(v_a,v_a_orig,d1,imax,jmax)                   ! m/s
        write(6,*)
        write(6,*)' Stats for variational adjustment on V (m/s)'
        call stats(d1,imax,jmax)

        call diff(p_a,p_a_orig,d1,imax,jmax)                   ! Pascals
        write(6,*)
        write(6,*)' Stats for variational adjustment on P (Pa)'
        call stats(d1,imax,jmax)

        write(6,*)

500	continue	! skip to here if cold starting or no backgrnd
c
c.....	Channel the winds around the terrain
c
	call get_directory('static', infile, len)
!	infile = infile(1:len) // '/pbl_top.dat'
!	call s_len(infile, len)
c	open(52,file='../static/surface/pbl_top.dat',
!	open(52,file=infile(1:len),
!    &       form='unformatted',status='old')
!	read(52) pbl_top
!	close(52)
cc	call vortdiv(u_a,v_a,vort,div,imax,jmax,dx,dy)
cc	call channel(u_a,v_a,topo,imax,jmax,pbl_top,pblht,dx,dy,z,div)
c
c.....	Calculate the final vorticity and divergence.
c
	call vortdiv(u_a,v_a,vort,div,imax,jmax,dx,dy)
c
c.....	Compute a vertical velocity by integrating the surface winds over 
c.....  some pbl and allow for terrain lift.
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
	   vv(i,j) = (dvh - div(i,j) * pblht) * 100.	! cm/sec
	enddo !i
	enddo !j
	call bounds(vv,imax,jmax)
c
c.....	Now convert some stuff and call the derived routines.
c
	call multcon(p_a,0.01,imax,jmax)
	call make_cssi(t,td,mslp,u_a,v_a,cssi,imax,jmax,badflag)
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
c
!!	print *,' '
!!	print *,' Plotting MSL pressure (mb): '
!!	call aplot(mslp, ni, nj)
c
	call multcon(p_a,100.,imax,jmax)	! conv mb to Pa
	call multcon(psfc,100.,imax,jmax)	! conv mb to Pa
	call multcon(mslp,100.,imax,jmax)	! conv mb to Pa
	call multcon(vv,.01,imax,jmax)		! conv cm/s to m/s
	call addcon(thetae,273.15,thetae,imax,jmax)	! C to  K
	call windspeed(u_a,v_a,spd,imax,jmax)	! calc windspeed
	call vlog2vis(vis,vis,imax,jmax)	! conv log(vis) to vis-miles
c
c.....  Calculate RH and change to %
c
	call hum(t,td,rh,imax,jmax,d1,d2)	! calc rel hum.
	call multcon(rh,100.,imax,jmax)
c
c.....  Adjust visibility analysis.
c
	call enhance_vis(i4time,vis,rh,topo,imax,jmax,kcloud)
	call conv_miles2m(vis,vis,imax,jmax)	! conv miles to meters
	vis_mx = -1.e25
	vis_mn =  1.e25
	do j=1,jmax
	do i=1,imax
	   if(vis(i,j) .gt. vis_mx) vis_mx = vis(i,j)
	   if(vis(i,j) .lt. vis_mn) vis_mn = vis(i,j)
	enddo !i
	enddo !j
	print *,' Vis check: Max, Min (meters)= ', vis_mx, vis_mn
	if(vis_mx.ge.400000. .or. vis_mn.lt.0.) then
	   print *,' Vis max or min out of range. Bag field.'
	   call constant(vis,badflag,imax,jmax)
	endif
c
c.....  Call the Fire index routine.
c
	print *,' '
	print *,' Fire Wx index...'
	isnow = 1
	ismoist = 1
	call zero(fire,imax,jmax)
	call zero(d1,imax,jmax)
	call zero(d2,imax,jmax)
c
	i4time_tol = 7200
	i4time_f = i4time - 3600
	ilev = -1
	var_fire = 'LSM'                   ! soil moisture
	ext_f = 'lm1'
	call get_laps_2dvar(i4time,i4time_tol,i4time_near,lat,lon,
     &                      dum1,dum2,
     &                      ext_f,var_fire,units_fire,com_fire,
     &                      imax,jmax,d1,ilev,istatus)
	if(istatus .ne. 1) then
	   print *,' Error getting soil moisture.'
	   ismoist = 0
	endif
c
        ilev = 0
	var_fire = 'SC '                   ! snow cover
	ext_f = 'lm2'
	call get_laps_2dvar(i4time,i4time_tol,i4time_near,lat,lon,
     &                      dum1,dum2,
     &                      ext_f,var_fire,units_fire,com_fire,
     &                      imax,jmax,d2,ilev,istatus)
	if(istatus .ne. 1) then
	   print *,' Error getting snow cover.'
	   isnow = 0
	endif
c
	call lp_fire_danger(imax,jmax,rh,t,spd,d1,d2,topo,ldf,ismoist,
     &                                            isnow,fire,istatus)
	print *,' Fire: ',ismoist, isnow, istatus
	print *,' '
c
c.....  Calculate Heat Index
c
	print *,' Heat Index...'
	call heat_index(t,rh,hi,imax,jmax,r_missing_data)
c
c.....	Now write out the grids to PROD_DEV.
c
 888	continue
	print *,' Saving primary fields.'
	do i=1,num_var
	  lvl(i) = 0
	  lvl_coord(i) = 'AGL'
	enddo !i
	lvl_coord(9) = 'MSL'
c
        do i=1,num_var
           write(comment(i),180) n_sao_g,n_sao_b
        enddo !i
 180	   format(49x,3i4)
        print*,comment(1)
c
	do k=1,num_var  ! fill the output array
	do j=1,jmax
	do i=1,imax
	   data(i,j,k) = fill_val
	enddo !i
	enddo !j
	enddo !k
c
c.....  Move the 2-d analyses to the 3-d storage array for writing.
c.....  Check the fields for NaN's and other bad stuff first.
c
	print *,
     1  ' ======================================================='
	print *,' u-wind (m/s):'
	var(1) = 'U'		! u-wind (m/s)
	units(1) = 'M/S'
	comment(1)= 'U (10m AGL)'
	call check_field_2d(u_a, imax,jmax,fill_val,istatus)
	if(istatus .eq. 1)
     &     call move_2dto3d(   u_a, data,  1, imax, jmax, num_var)
c
	print *,' -------------------------------'
	print *,' v-wind (m/s):'
	var(2) = 'V'		! v-wind (m/s)
	units(2) = 'M/S'
	comment(2)= 'V (10m AGL)'
	call check_field_2d(v_a, imax,jmax,fill_val,istatus)
	if(istatus .eq. 1)
     &     call move_2dto3d(   v_a, data,  2, imax, jmax, num_var)
c
	print *,' -------------------------------'
	print *, redp_lvl,' m pressure (pa):'
	var(3) = 'P'		! reduced press (Pa)
	units(3) = 'PA'
	write(comment(3)(1:4),181) ifix(redp_lvl)
 181	format(i4)
	comment(3)(5:23) = ' M REDUCED PRESSURE'
	call check_field_2d(p_a, imax,jmax,fill_val,istatus)
cc	do j=1,jmax
cc	do i=1,imax
cc	  if(p_a(i,j).lt.80000. .or. p_a(i,j).gt.105000.) then
cc	     istatus = 0
cc	     print *,' Value out of range at ',i,j
cc	     go to 1181
cc	  endif
cc	enddo !i
cc	enddo !j
1181    continue
	print *, 'rp istatus = ', istatus
	if(istatus .eq. 1)
     &     call move_2dto3d(   p_a, data,  3, imax, jmax, num_var)
c
	print *,' -------------------------------'
	print *,' temp (k):'
	var(4) = 'T'		! temp (K)
	units(4) = 'K'
	comment(4)= 'T (2m AGL)'
	call check_field_2d(t, imax,jmax,fill_val,istatus)
	if(istatus .eq. 1)
     &     call move_2dto3d(     t, data,  4, imax, jmax, num_var)
c
	print *,' -------------------------------'
	print *,' dewpt (k):'
	var(5) = 'TD'		! dew point (K)
	units(5) = 'K'
	comment(5)= 'TD (2m AGL)'
	call check_field_2d(td, imax,jmax,fill_val,istatus)
	if(istatus .eq. 1)
     &     call move_2dto3d(    td, data,  5, imax, jmax, num_var)
c
	print *,' -------------------------------'
	print *,' vert vel (m/s):'
	var(6) = 'VV'		! vert. vel (m/s)
	units(6) = 'M/S'
	call check_field_2d(vv, imax,jmax,fill_val,istatus)
	if(istatus .eq. 1)
     &     call move_2dto3d(    vv, data,  6, imax, jmax, num_var)
c
	print *,' -------------------------------'
	print *,' rel hum:'
	var(7) = 'RH'		! relative humidity (%)
	units(7) = '%'
	comment(7)= 'RH (2m AGL)'
	call check_field_2d(rh, imax,jmax,fill_val,istatus)
	if(istatus .eq. 1)
     &     call move_2dto3d(    rh, data,  7, imax, jmax, num_var)
c
	var(8) = 'HI'		! Heat Index (K)
	units(8) = 'K'
	comment(8)= 'HEAT INDEX'
	call move_2dto3d(    hi, data, 8, imax, jmax, num_var)
c
	print *,' -------------------------------'
	print *,' MSL pressure (pa) :'
	var(9) = 'MSL'		! MSL pressure (Pa)
	units(9) = 'PA'
	comment(9) = 'MSL PRESSURE'
	call check_field_2d(mslp, imax,jmax,fill_val,istatus)
	do j=1,jmax
	do i=1,imax
	  if(mslp(i,j).lt.85000. .or. mslp(i,j).gt.110000.) then
	     istatus = 0
             write(6,*)' ERROR: mslp out of range at ',i,j,mslp(i,j)       
	     go to 1182
	  endif
	enddo !i
	enddo !j
1182    continue
	if(istatus .eq. 1)
     &     call move_2dto3d(  mslp, data,  9, imax, jmax, num_var)
c
	print *,' -------------------------------'
	print *,' temp adv (K/s):'
	var(10) = 'TAD'		! temperature advection (K/sec)
	units(10) = 'K/S'
	call check_field_2d(tadv, imax,jmax,fill_val,istatus)
	if(istatus .eq. 1)
     &     call move_2dto3d(  tadv, data, 10, imax, jmax, num_var)
c
	print *,' -------------------------------'
	print *,' theta (k):'
	var(11) = 'TH'		! potential temp (K)
	units(11) = 'K'
	call check_field_2d(theta, imax,jmax,fill_val,istatus)
	if(istatus .eq. 1)
     &     call move_2dto3d( theta, data, 11, imax, jmax, num_var)
c
	print *,' -------------------------------'
	print *,' thetae'
	var(12) = 'THE'		! equiv pot temp (K)
	units(12) = 'K'
	call check_field_2d(thetae, imax,jmax,fill_val,istatus)
	if(istatus .eq. 1)
     &     call move_2dto3d(thetae, data, 12, imax, jmax, num_var)
c
	print *,' -------------------------------'
	print *,' sfc p (pa)'
	var(13) = 'PS'		! surface press (Pa)
	units(13) = 'PA'
	comment(13) = 'UNREDUCED SURFACE PRESSURE (0m AGL)'
	call check_field_2d(psfc, imax,jmax,fill_val,istatus)
	if(istatus .eq. 1)
     &     call move_2dto3d(  psfc, data, 13, imax, jmax, num_var)
c
	print *,' -------------------------------'
	print *,' vort (/s)'
	var(14) = 'VOR'		! sfc vorticity (/s)
	units(14) = '/S'
	comment(14) = 'VORTICITY (10m AGL)'
	call check_field_2d(vort, imax,jmax,fill_val,istatus)
	if(istatus .eq. 1)
     &     call move_2dto3d(  vort, data, 14, imax, jmax, num_var)
c
	print *,' -------------------------------'
	print *,' mix ratio (g/kg)'
	var(15) = 'MR'		! mixing ratio (g/kg)
	units(15) = 'G/KG'
	comment(15) = 'MR (2m AGL)'
	call check_field_2d(q, imax,jmax,fill_val,istatus)
	if(istatus .eq. 1)
     &     call move_2dto3d(     q, data, 15, imax, jmax, num_var)
c
	print *,' -------------------------------'
	print *,' moist conv (g/kg/s)'
	var(16) = 'MRC'		! moisture convergence (g/kg/s)
	units(16) = 'G/KG/S'
	call check_field_2d(qcon, imax,jmax,fill_val,istatus)
	if(istatus .eq. 1)
     &     call move_2dto3d(  qcon, data, 16, imax, jmax, num_var)
c
	print *,' -------------------------------'
	print *,' div (/s)'
	var(17) = 'DIV'		! sfc divergence (/s)
	units(17) = '/S'
	comment(17) = 'DIV (10m AGL)'
	call check_field_2d(div, imax,jmax,fill_val,istatus)
	if(istatus .eq. 1)
     &     call move_2dto3d(   div, data, 17, imax, jmax, num_var)
c
	print *,' -------------------------------'
	print *,' theta adv (K/s)'
	var(18) = 'THA'		! pot temp adv (K/s)
	units(18) = 'K/S'
	call check_field_2d(thadv, imax,jmax,fill_val,istatus)
	if(istatus .eq. 1)
     &     call move_2dto3d( thadv, data, 18, imax, jmax, num_var)
c
	print *,' -------------------------------'
	print *,' moist adv (g/kg/s)'
	var(19) = 'MRA'		! moisture adv (g/kg/s)
	units(19) = 'G/KG/S'
	call check_field_2d(qadv, imax,jmax,fill_val,istatus)
	if(istatus .eq. 1)
     &     call move_2dto3d(  qadv, data, 19, imax, jmax, num_var)
c
	print *,' -------------------------------'
	print *,' wind spd (m/s):'
	var(20) = 'SPD'		! wind speed (m/s)
	units(20) = 'M/S'
	comment(20) = 'SPD (10m AGL)'
	call check_field_2d(spd, imax,jmax,fill_val,istatus)
	if(istatus .eq. 1)
     &     call move_2dto3d(   spd, data, 20, imax, jmax, num_var)
c
	var(21) = 'CSS'		! CSSI 
	units(21) = ' '
	comment(21)= 'CSSI - COLORADO SEVERE STORM INDEX'
        call move_2dto3d(  cssi, data, 21, imax, jmax, num_var)
c
	print *,' -------------------------------'
	print *,' vis (m):'
	var(22) = 'VIS'		! Visibility (m)
	units(22) = 'M'
	call check_field_2d(vis, imax,jmax,fill_val,istatus)
	if(istatus .eq. 1)
     &     call move_2dto3d(   vis, data, 22, imax, jmax, num_var)
c
	var(23) = 'FWX'		! Fire threat index (integer)
	units(23) = ' '
	comment(23)(1:22)= 'LAPS FIRE THREAT INDEX'
	comment(23)(63:121) = 
     &      'INDEX: 0-NONE, 5-SLGT, 10-MDT, 15-HI, 20-EXTREME'
	call move_2dto3d(  fire, data, 23, imax, jmax, num_var)
c
	var(24) = 'TGD'		! Ground Temperature
	units(24) = 'K'
	comment(24) = 'TGD (0m AGL)'
	call move_2dto3d(  tgd_k, data, 24, imax, jmax, num_var)       
c
	print *,
     1  ' ======================================================='
c
c.....  Now actually write the LSX file.
c
        if(iwrite_output .ge. 0)then
	    call get_directory('lsx', dir, len)
	    ext = 'lsx'
	    call write_laps_data(i4time,dir,ext,imax,jmax,num_var,
     &             num_var,var,lvl,lvl_coord,units,comment,data,istatus)      
            write(6,*)' LSX file write completed, istatus = ',istatus
        endif
c
	jstatus(3) = 1		! everything ok...

        I4_elapsed = ishow_timer()
c
c.....  Now finish up with some verification.  Expected accuracys
c.....  based on FMH-1 Appendix C, but fixed estimates for normal 
c.....  conditions.
c
	iunit = 11
	call get_directory('log', ver_file, len)
c	ver_file = '../log/qc/laps_sfc.ver.'//filename(6:9)
	ver_file = ver_file(1:len)//'qc/laps_sfc.ver.'//filename(6:9)
	call s_len(ver_file, len)
	open(iunit,file=ver_file(1:len),status='unknown',err=999)
c	
	title = 'Temperature background verification (deg F)'
	ea = 1.50
	call zero(d1,imax,jmax)
	call move(t_bk,d1,imax,jmax)
	call verify(d1,t_s,stn,n_obs_b,title,iunit,
     &              ni,nj,mxstn,x1a,x2a,y2a,ii,jj,ea,badflag)
c
	title = 'Temperature verification (deg F)'
	ea = 1.50
	call zero(d1,imax,jmax)
	call conv_k2f(t,d1,imax,jmax)
	call verify(d1,t_s,stn,n_obs_b,title,iunit,
     &              ni,nj,mxstn,x1a,x2a,y2a,ii,jj,ea,badflag)
c
	title = 'Dew Point background verification  (deg F)'
	ea = 2.00
	call zero(d1,imax,jmax)
	call move(td_bk,d1,imax,jmax)
	call verify(d1,td_s,stn,n_obs_b,title,iunit,
     &              ni,nj,mxstn,x1a,x2a,y2a,ii,jj,ea,badflag)
c	
	title = 'Dew Point verification (deg F)'
	ea = 2.00
	call zero(d1,imax,jmax)
	call conv_k2f(td,d1,imax,jmax)
	call verify(d1,td_s,stn,n_obs_b,title,iunit,
     &              ni,nj,mxstn,x1a,x2a,y2a,ii,jj,ea,badflag)
c
        do ista = 1,n_obs_b
            if(dd_s(ista) .ne. badflag .and. 
     1         ff_s(ista) .ne. badflag)then      
                call disptrue_to_uvgrid(dd_s(ista),ff_s(ista)
     1                                 ,u_s(ista),v_s(ista),lon_s(ista))       
            else
                u_s(ista) = badflag
                v_s(ista) = badflag
            endif
        enddo ! ista

        title = 'U Wind Component background verification (kt)'
 	ea = 2.00
!	call zero(d1,imax,jmax)
!	call conv_ms2kt(u_bk,d1,imax,jmax)
        call verify(u_bk,u_s,stn,n_obs_b,title,iunit,
     &              ni,nj,mxstn,x1a,x2a,y2a,ii,jj,ea,badflag)
c
 	title = 'U Wind Component verification (kt)'
 	ea = 2.00
 	call zero(d1,imax,jmax)
 	call conv_ms2kt(u_a,d1,imax,jmax)
 	call verify(d1,u_s,stn,n_obs_b,title,iunit,
     &              ni,nj,mxstn,x1a,x2a,y2a,ii,jj,ea,badflag)

        title = 'V Wind Component background verification (kt)'
 	ea = 2.00
!	call zero(d1,imax,jmax)
!	call conv_ms2kt(v_bk,d1,imax,jmax)
        call verify(v_bk,v_s,stn,n_obs_b,title,iunit,
     &              ni,nj,mxstn,x1a,x2a,y2a,ii,jj,ea,badflag)
c
 	title = 'V Wind Component verification (kt)'
 	ea = 2.00
 	call zero(d1,imax,jmax)
 	call conv_ms2kt(v_a,d1,imax,jmax)
 	call verify(d1,v_s,stn,n_obs_b,title,iunit,
     &              ni,nj,mxstn,x1a,x2a,y2a,ii,jj,ea,badflag)
c
	title = 'Wind Speed background verification (kt)'
	ea = 2.00
	call zero(d1,imax,jmax)
	call windspeed(u_bk,v_bk,d1,imax,jmax)	! calc windspeed (kt)
	call verify(d1,ff_s,stn,n_obs_b,title,iunit,
     &              ni,nj,mxstn,x1a,x2a,y2a,ii,jj,ea,badflag)
c
	title = 'Wind Speed verification (kt)'
	ea = 2.00
	call zero(d1,imax,jmax)
	call conv_ms2kt(spd,d1,imax,jmax)
	call verify(d1,ff_s,stn,n_obs_b,title,iunit,
     &              ni,nj,mxstn,x1a,x2a,y2a,ii,jj,ea,badflag)
c
	title = 'MSL pressure background verification (mb)'   
	ea = 0.68
	call zero(d1,imax,jmax)
	call move(mslp_bk,d1,imax,jmax)
!	call multcon(d1,.01,imax,jmax)
	call verify(d1,mslp_s,stn,n_obs_b,title,iunit,
     &              ni,nj,mxstn,x1a,x2a,y2a,ii,jj,ea,badflag)
c
	title = 'MSL pressure verification (mb)'   
	ea = 0.68
	call zero(d1,imax,jmax)
	call move(mslp,d1,imax,jmax)
	call multcon(d1,.01,imax,jmax)
	call verify(d1,mslp_s,stn,n_obs_b,title,iunit,
     &              ni,nj,mxstn,x1a,x2a,y2a,ii,jj,ea,badflag)
c
	title = 'Reduced pressure background verification (mb)'   
	ea = 0.68
	call zero(d1,imax,jmax)
	call move(rp_bk,d1,imax,jmax)
!	call multcon(d1,.01,imax,jmax)
	call verify(d1,pred_s,stn,n_obs_b,title,iunit,
     &              ni,nj,mxstn,x1a,x2a,y2a,ii,jj,ea,badflag)
c
	title = 'Reduced pressure verification (mb)'   
	ea = 0.68
	call zero(d1,imax,jmax)
	call move(p_a,d1,imax,jmax)
	call multcon(d1,.01,imax,jmax)
	call verify(d1,pred_s,stn,n_obs_b,title,iunit,
     &              ni,nj,mxstn,x1a,x2a,y2a,ii,jj,ea,badflag)
c	
	title = 'Ground Temperature background verification (deg F)'
	ea = 1.50
        call get_sfcob_field(obs,mxstn,'tgd',ob_full,istatus)
	call verify(tgd_bk_f,ob_full,stn,n_obs_b,title,iunit,
     &              ni,nj,mxstn,x1a,x2a,y2a,ii,jj,ea,badflag)
c
	title = 'Ground Temperature verification (deg F)'
	ea = 1.50
	call conv_k2f(tgd_k,d1,imax,jmax)
	call verify(d1,ob_full,stn,n_obs_b,title,iunit,
     &              ni,nj,mxstn,x1a,x2a,y2a,ii,jj,ea,badflag)
c
	close(iunit)
c
c.....  That's it.  Let go home.
c
	print *,' Normal completion of LAPSVANL'
        I4_elapsed = ishow_timer()
	return

 999	print *,'ERROR opening ',ver_file(1:len)
	return
c
	end

