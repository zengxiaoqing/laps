cdis
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
        program locpost
        integer istatus

        real, allocatable, dimension(:,:) :: lat
        real, allocatable, dimension(:,:) :: lon
        real, allocatable, dimension(:,:) :: topo
        real, allocatable, dimension(:,:) :: ldf

        call get_grid_dim_xy(NX_L, NY_L, istatus)
        if (istatus .ne. 1) then
            write(6,*) 'return get_grid_dim_xy, status: ', istatus
            stop
        endif

!       Allocate static arrays (lat, lon, topo, ldf)
        allocate( lat(NX_L,NY_L), STAT=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' ERROR: Could not allocate lat'
            stop
        endif

        allocate( lon(NX_L,NY_L), STAT=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' ERROR: Could not allocate lon'
            stop
        endif

        allocate( topo(NX_L,NY_L), STAT=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' ERROR: Could not allocate topo'
            stop
        endif

        allocate( ldf(NX_L,NY_L), STAT=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' ERROR: Could not allocate ldf'
            stop
        endif

!       Read static arrays (lat, lon, topo, ldf)
        call read_static_grid(NX_L,NY_L,'LAT',lat,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error getting LAPS LAT'
            stop
        endif

        call read_static_grid(NX_L,NY_L,'LON',lon,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error getting LAPS LON'
            stop
        endif

        call read_static_grid(NX_L,NY_L,'AVG',topo,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error getting LAPS AVG'
            stop
        endif

        call read_static_grid(NX_L,NY_L,'LDF',ldf,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error getting LAPS LDF'
            stop
        endif

        call build_sfc_static_sub(NX_L,NY_L,lat,lon,topo)

        lun = 31
        lun_out = 19

        call locpost_radar(NX_L,NY_L,lat,lon,topo,ldf
     1                    ,lun,lun_out,istatus)

        end

	subroutine build_sfc_static_sub(ni,nj,lat,lon,topo)
c
c*****************************************************************************
c
c	Program to build the necessary static files for the LAPS surface
c       analysis.  Files created include:
c
c           REMOVED JUL 95 --> 'fnorm2.lut' (the barnes look-up table), 
c           REMOVED JAN 97 --> 'bkgwts.wts' (background weights), 
c                 'drag_coef.dat' (surface roughness), 
c             and 'pbl_top.dat' (boundary layer depth).
c
c	Changes:
c 	P.A. Stamus	04-12-94  Original from the individual programs.
c                       09-01-94  Porting changes.
c                       07-26-95  Remove fnorm2 calcs.
c                       01-08-97  Porting changes. Remove 'stations.in'
c                                  requirement, general clean up.
c                       04-09-97  Remove equivalences.
c
c****************************************************************************
c
        implicit none
ccccc	include 'laps_sfc.inc'
c
c.....	Grids for the outputs, weights, and stuff 
c
        integer ni,nj
	real akk(ni,nj)
	real top(ni,nj), d1(ni,nj), d2(ni,nj)
        real fnorm(0:ni-1,0:nj-1)
c
        character*80 grid_fnam_common
        common/grid_fnam_cmn/ grid_fnam_common
c
c..... LAPS Lat/lon grids.
c
	real lat(ni,nj), lon(ni,nj), topo(ni,nj)
	integer grid_spacing, len
	character dir_s*150,ext_s*31,units*10,comment*125,var_s*3
        integer imx,jmx, imn, jmn, icnt, npass, i_tmn, j_tmn
        integer i_tmx, j_tmx, n_obs_var, i, j, imaxm1, jmaxm1
        integer imax, jmax, istatus
        real slope1, slope2, slp, std, ave, smin, smax, con, scale
        real co_max_scl,co_topo_range , termx, x, pi, fno, smsng, rom2
        real topo_range, d, pbldz, top_adj, topmx, topo_mx, topo_mn
c
c
c.....  START:
c.....	Get the LAPS lat/lon and topo data here so we can pass them to the 
c.....	routines that need them.
c

!        call get_laps_config(laps_domain,istatus)
!	dir_s = '../static/' 
!	dir_s = './' 
        call get_directory('static',dir_s,len)
        ext_s = 'nest7grid' 
c
	imax = ni
	jmax = nj
	imaxm1 = imax - 1
	jmaxm1 = jmax - 1
c
	topo_mx = -1.e30
	topo_mn =  1.e30
	do j=1,jmax
	do i=1,imax
	   if(topo(i,j) .gt. topo_mx) then
	      topo_mx = topo(i,j)
	      i_tmx = i
	      j_tmx = j
	   endif
	   if(topo(i,j) .lt. topo_mn) then
	      topo_mn = topo(i,j)
	      i_tmn = i
	      j_tmn = j
	   endif
	enddo !i
	enddo !j
	topo_range = topo_mx - topo_mn
c
	print *,' '
	print *,' Topography: '
	write(6,211) topo_mx, i_tmx, j_tmx
 211	format(1x,'Max of ',f10.1,' at ',i3,',',i3)
	write(6,212) topo_mn, i_tmn, j_tmn
 212	format(1x,'Min of ',f10.1,' at ',i3,',',i3)
	write(6,213) topo_range
 213	format(1x,'Range of ',f10.1)
	print *,' '
c
	pi = 4. * atan(1.0) 
	x = 1.e-3
	fno = 1. / (sqrt(2. * pi)) 
c
c..........................................................................
c
c	This section calculates the drag coefficient for the LAPS domain.
c
c	Original by John McGinley.  Date unknown.
c	Changes:
c		P. Stamus	11-12-19  Changes for new LAPS grids.
c				11-20-91  Changed scaling to reduce akk.
c                               03-17-95  Auto scaling code for other domains.
c
c..........................................................................
c
	co_max_scl = 51.
	co_topo_range = 3099.8
	call constant(akk,1.,imax,jmax)
c
	smax = -1.e30
	smin = 1.e30
	icnt = 0
	std = 0.
	ave = 0.
	do j=2,jmaxm1
	do i=2,imaxm1
           slope1 = topo(i,j+1) - topo(i,j-1)
           slope2 = topo(i+1,j) - topo(i-1,j)
           if(slope1 .eq. 0.) slope1 = 0.001
           if(slope2 .eq. 0.) slope2 = 0.001
           slp = sqrt((slope1)**2 + (slope2)**2)
           icnt = icnt + 1
           akk(i,j) = slp
           ave = ave + slp
           std = std + slp ** 2
           if(slp .gt. smax) then
	      smax = slp
	      imx = i
	      jmx = j
	   endif
           if(slp .lt. smin) then
	      smin = slp
	      imn = i
	      jmn = j
	   endif
	enddo !i
	enddo !j
c
	ave = ave / float(icnt)
	std = sqrt(std / float(icnt))
c
c.....	We want this factor to serve as a  multiplier on the drag coef 
c.....  constant in lapsanal.  The scale factor was arbitrarilly set so
c.....  'akk' ranged between 51 and 1.1.  Here, the range of the topography
c.....  in the domain is scaled against the Colorado version, so the scaling
c.....  will adjust to smoother (or rougher) domains.
c
	con = co_max_scl * (topo_range / co_topo_range)
	scale = (con / smax) * (smin + 1.)
c
	do j=2,jmaxm1
	do i=2,imaxm1
	   akk(i,j) = scale * akk(i,j) / (smin + 1.)  !AKK GOES FROM 51 TO 1.1
	   if(akk(i,j) .lt. 1.) akk(i,j) = 1.
	enddo !i
	enddo !j
c
         
	open(19,file=dir_s(1:len)//'drag_coef.dat',
     &                    form='unformatted',status='unknown')
	print *,' Drag coef.: '
	write(6,1000)ave,std,smax,imx,jmx,smin,imn,jmn,scale
1000	format(1x,'Ave: ',f12.4,' Std Dev: ',f12.4,/,' Max: ',f12.4,
     &         ' at ',2i3,/,' Min: ',f12.4,' at ',2i3,/,
     &         ' Calc scale factor: ',f12.4)
	write(19) akk
        close(19)
	print *,' '
	print *,' Normal completion of DRAG_COEF.'



        goto 9999 ! Bypass PBL top code

c
c
c...........................................................................
c
c	This section calculates the top of the pbl across the LAPS domain.
c
c	Original by John McGinley (somewhere in the depths of time).
c	Changes:
c		P. Stamus	11-12-91  Changes for new LAPS grids.
c				11-15-91  Incr plbtop to prevent hi wnds
c
c...........................................................................
c
	termx = 0.
	do j=1,nj
	do i=1,ni
	  if(topo(i,j) .gt. termx) termx = topo(i,j)
        enddo !i
        enddo !j
c
        npass = 1
	smsng = 1.e37
c
	rom2 = 0.01
	call dynamic_wts(imax,jmax,n_obs_var,rom2,d,fnorm)
	call barnes2(top,imax,jmax,topo,smsng,npass,d1,d2,fnorm)
c
	topmx = 0.
	do j=1,nj
	do i=1,ni
	  if(top(i,j) .gt. topmx) topmx = top(i,j)
	enddo !i
	enddo !j
c
        write(6,302) termx, topmx
 302   format(1x,' Max. terrain ',f8.0,'; Max. smooth terrain ',f8.0)
c
!	pbldz=475.	! pbl depth...old value 1000 m
!	pbldz=775.	! pbl depth...old value 1000 m
!	pbldz=925.	! pbl depth...old value 1000 m
	pbldz = 1225.
c
	top_adj = termx - topmx + pbldz
	do j=1,nj
	do i=1,ni
	  top(i,j) = top(i,j) + top_adj
	  if((top(i,j) - topo(i,j)) .lt. 100.) then
	    write (6,300) top(i,j),topo(i,j),i,j
	    top(i,j) = topo(i,j) + 100.
	  endif
	enddo !i
	enddo !j
 300    format(1x,'top pbl ',f8.0,' top terrain ',f8.0,' at ',2i6)
c
	open(10,file=dir_s(1:len)//'/pbl_top.dat',
     &          form='unformatted',status='unknown')
	write(10) top
	close(10)
	print *,' '
	print *,' Normal completion of PBL_TOP'
c
 9999   print *,' '
	print *,' Normal completion of build_sfc_static'
c
	return
	end
c
c
	subroutine barnes2(t,imax,jmax,to,smsng,npass,h1,h2, fnorm)
c
ccc	include '../../source/sfc/laps_sfc.inc'
        implicit none
        integer imax, jmax
	real to(imax,jmax), t(imax,jmax), val(imax*jmax)
	real h1(imax,jmax), h2(imax,jmax)
	integer iob(imax*jmax), job(imax*jmax), dx, dy
        real badd
	parameter(badd = 1.e6 - 2.)	! bad data value
        real sum, sumwt, sum2, sumwt2, smsng
        integer i,j,n, ipass, ncnt, npass 
	real fnorm(0:imax-1,0:jmax-1)

	call zero(h1,imax,jmax)
	call zero(h2,imax,jmax)
c
c
c.....	loop over field npass times 
c
!	print *,' *** In BARNES2 ***'
        ncnt = 0
        do j=1,jmax
        do i=1,imax
          if (to(i,j).ne.0. .and. to(i,j).lt.badd) then
            ncnt = ncnt + 1
            iob(ncnt) = i
            job(ncnt) = j 
            val(ncnt) = to(i,j)
          endif
        enddo !i
        enddo !j
	if(ncnt .eq. 0) then
	  print *,' *** NCNT = 0 in BARNES2. ***'
	  return
	endif
c
	do ipass=1,npass
c

	  do j=1,jmax
	  do i=1,imax
	    sum = 0.
	    sumwt = 0.
	    sum2 = 0.
	    sumwt2 = 0.
	    do n=1,ncnt
	      dy = abs(j - job(n))
	      dx = abs(i - iob(n))
cc	      sum = fnorm2(dx,dy,kdim) * val(n) + sum
	      sum2 = fnorm(dx,dy) * val(n) + sum2
cc	      sumwt = sumwt + fnorm2(dx,dy,kdim)
	      sumwt2 = sumwt2 + fnorm(dx,dy)
	    enddo !n
c
	    if(sumwt2 .eq. 0.) then
	      if(ipass .eq. 1) then
cc	        if(kdim .eq. 1) then
cc	          print *,
cc     &              ' *** WARNING. sumwt = 0, kdim = 1 in BARNES2 ***'
cc	          go to 500
cc	        else
		   print *,' got into wierd loop.............'
cc	          do kk=kdim,1,-2
	            sum2 = 0.
	            sumwt2 = 0.
	            do n=1,ncnt
	              dx = abs(i - iob(n))
	              dy = abs(j - job(n))
cc	              sum = fnorm2(dx,dy,kk) * val(n) + sum
	              sum2 = (fnorm(dx,dy) + .01) * val(n) + sum2
	              sumwt2 = sumwt2 + (fnorm(dx,dy) + .01)
	            enddo !n
	            if(sumwt2 .ne. 0.) go to 490
cc	          enddo !kk
cc	        endif
	      else
	        go to 500
	      endif
	    endif
c
490	    continue
cc	    t(i,j) = sum / sumwt
	    t(i,j) = sum2 / sumwt2
c
500	  continue
          enddo !i
	  enddo !j
c
	  if(ipass .eq. 2) then
	    call move(h2,to,imax,jmax)
	    call diff(h1,t,t,imax,jmax)
!	    write(9,915)
!915	    format(' after 2nd barnes pass')
!	    write(9,909) rom2
!909	    format(' radm2 = ',f8.4)
	    go to 550
	  endif
!	  write(9,912)
912	  format(' after 1st barnes pass')
	  if(npass .eq. 1) go to 550
	  call move(t,h1,imax,jmax)
 	  call move(to,h2,imax,jmax)
	  do n=1,ncnt
	    val(n) = t(iob(n),job(n)) - val(n)  
	  enddo !n
550	continue
        enddo !ipass
c
!	print *,' *** BARNES2 Done. ***'
	return
	end
c
c
	subroutine dynamic_wts(imax,jmax,n_obs_var,rom2,d, fnorm)
c
c=====================================================================
c
c     Routine to calculate the weights to be used in the Barnes
c     analysis.  The data density is used to set the cutoff for
c     the response function.  Then that cutoff is used to calculate
c     the exp, based on differences so that no additional distance
c     calculations are required in the Barnes routine.  All of this
c     is done in gridpoint space.
c
c     Original:  07-14-95  P. Stamus, NOAA/FSL
c     Changes:
c
c     Notes:
c
c       1.  If variable 'rom2' is passed in as zero, it is calculated
c           from the data density.  Otherwise, the value passed in is
c           used in the weight calculation.
c
c       2.  The response for 2d waves is hard-wired in this routine.
c           This is the 'con' variable, and comes from the eqn:
c                     D = exp -(pi**2 R**2)/lamba**2
c           If we set D (the response) to our desired cutoff, set 
c           lamba to the desired wavelength in gridpt space (2d),
c           then solve for R in terms of d, we get the 'con' value
c           (i.e.,  R = (con)*d).  Here are some values for different
c           cutoffs:
c                     D = 0.01     R = 1.36616d
c                         0.10     R = 0.96602d
c                         0.25     R = 0.74956d
c                         0.50     R = 0.53002d
c
c=====================================================================
c
	implicit none
        integer imax,jmax 
        real fnorm(0:imax-1,0:jmax-1)
        real pi, con, area, fno, rr, d, rom2
        integer n_obs_var,  dx, dy
          
ccc	include '../../source/sfc/laps_sfc.inc'
c
c.... First, find the area that each ob covers in gridpt space (this
c.... of course assumes a uniform coverage).
c
cc	con = 0.96602     ! resp of 0.10
cc	con = 0.74956     ! resp of 0.25
	con = 0.53002     ! resp of 0.50
	if(rom2 .eq. 0.) then
	   area = float(imax * jmax) / n_obs_var
	   d = sqrt( area )
	   rom2 = 1. / ((con * con) * (d * d))
	   write(6,900) n_obs_var, area, d, rom2
 900	   format(1x,'Num obs: ',i5,'  Area: ',f8.2,'  d: ',f8.2,
     &       '  rom2: ',f8.5)
	else
	   d = sqrt(1./(con * con * rom2))
	   write(6,902) rom2, d
 902	   format(' Using preset rom2 of: ',f8.5,'  Calc d: ',f8.2)
	endif
c
c.... Now calculate the weights for all the possible distances.
c
	pi = 4. * atan(1.)
	fno = 1. / (sqrt(2. * pi))
c
	do dy=0,jmax-1
	do dx=0,imax-1
	   rr = dx*dx + dy*dy
	   fnorm(dx,dy) = fno * (exp( -(rr * rom2)))
	enddo !dx
	enddo !dy
c
c.... That's it.
c
	return
	end
