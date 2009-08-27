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
	subroutine qcdata(filename,infile_l,rely,mxstn,
     &     t_s, td_s, dd_s, ff_s, ddg_s, ffg_s, pstn_s, pmsl_s, alt_s, 
     &     vis_s, stn, rii, rjj, ii, jj, n_obs_b, n_sao_b, n_sao_g,
     &     ni, nj, t_bk_f, td_bk_f, mslp_bk_mb, 
     &     thresh_t,thresh_td,thresh_mslp,
     &     iback_t, iback_td, iback_mp,
     &     istatus)
c
c=========================================================================
c
c       LAPS Quality Control routine.
c
c       Original by C. Hartsough, FSL  c. 1992
c       Changes:  
c           P.Stamus  20 Dec 1996  Porting changes for go anywhere LAPS
c                     25 Aug 1997  Changes for dynamic LAPS.
c                     17 Nov 1997  Fill empty arrays (HP compiler prob).
c                     13 Jul 1999  Remove some variables not used in sfc.
c                                  Rm *4 from all declarations.
c
c=========================================================================
c
c
c..... Stuff for the sfc data and other station info (LSO +)
c
        use mem_namelist, ONLY: l_dev_ck

	real t_s(mxstn), td_s(mxstn), dd_s(mxstn), ff_s(mxstn)
	real ddg_s(mxstn), ffg_s(mxstn), vis_s(mxstn)
	real pstn_s(mxstn), pmsl_s(mxstn), alt_s(mxstn)
	real rii(mxstn), rjj(mxstn)
c
	integer ii(mxstn), jj(mxstn)
c
	character stn(mxstn)*3, c_field*4
c
c..... Model background 
        real t_bk_f(ni,nj)
        real td_bk_f(ni,nj)
        real mslp_bk_mb(ni,nj)
c
c..... Arrays for the prev hour's OBS file input data
c                            
	real lat_l(mxstn),lon_l(mxstn),elev_l(mxstn)
	real td_l(mxstn),t_l(mxstn)
	real dd_l(mxstn),ff_l(mxstn),ddg_l(mxstn),ffg_l(mxstn)
	real pstn_l(mxstn),pmsl_l(mxstn),alt_l(mxstn)
	real store_hgt_l(mxstn,5),ceil_l(mxstn),lowcld_l(mxstn)
	real cover_l(mxstn),vis_l(mxstn),rad_l(mxstn)
	Integer obstime_l(mxstn),kloud_l(mxstn),idp3_l(mxstn)
	Character  infile_l*256, atime_l*24, stn_l(mxstn)*3
	character  obstype_l(mxstn)*8,wx_l(mxstn)*8
	character  store_emv_l(mxstn,5)*1, store_amt_l(mxstn,5)*4
c
c	Other files for internal use...
c
	integer rely(26,mxstn)            ! QC flags for current cycle
	integer rely_l(26,mxstn)          ! QC flags for previous cycle
        integer ivals(mxstn)              
        integer ivals_l(mxstn)            ! Ob indices for previous cycle
	integer istatus, jstatus
	character filename*9, outfile*256
c
c
c.....  Start here.
c
        write(6,*)
        write(6,*)' Subroutine qc_data...'

	n_obs_curr = n_obs_b
!       missing = -99.9
        badflag = -99.9
	imissing = -99
	jstatus = 0
c
	n_meso_g = 0          !old mesonet gone...
	n_meso_pos = 0
c
	do n=1,mxstn
	   ivals_l(n) = imissing
	do i=1,26
	   rely(i,n)   = imissing
	   rely_l(i,n) = imissing
	enddo !i
	enddo !n
c
c.....	OPEN QC FILE
c
	call get_directory('log', outfile, len)
	outfile = outfile(1:len) // 'sfcqc.log.' // filename(6:9)
	print *, ' opening log file for sfc qc:', outfile
	open(60,file=outfile,status='unknown')
c
	write(60,*)' ** LAPS Surface Quality Control **'
	write(60,*)' ** begin qc on surface data ',filename ,' ** '
c
c.....  Get previous sfc data (current passed in via argument list)
c
	call read_surface_old(infile_l,mxstn,atime_l,n_meso_g_l,
     &  n_meso_pos_l,n_sao_g_l,n_sao_pos_g_l,n_sao_b_l,n_sao_pos_b_l,
     &	n_obs_g_l,n_obs_pos_g_l,n_obs_b_l,n_obs_pos_b_l,stn_l,obstype_l,
     &	lat_l,lon_l,elev_l,wx_l,t_l,td_l,dd_l,ff_l,ddg_l,ffg_l,pstn_l,
     &	pmsl_l,alt_l,kloud_l,ceil_l,lowcld_l,cover_l,rad_l,idp3_l,
     &	store_emv_l,store_amt_l,store_hgt_l,vis_l,obstime_l,istatus)
	if(istatus .ne. 1) then
	 write(60,*) ' ERROR: could not read in data. Abort.'
	 goto 999
	end if
c
c.....  climatological extreme test            
c                 
	write(60,*) ' # of stns available: ',
     &                         n_obs_b, n_obs_b_l, n_obs_pos_b_l
	if(n_obs_b.le.0 .or. n_obs_b_l.le.0) then
	 write(60,*) ' No data for current and/or previous hr..Abort QC'      
	 return
	endif
c
	call time_ck2(stn,n_obs_b,stn_l,n_obs_b_l,ivals_l,stn,       
     &	              n_obs_curr)
	do 10 n = 1, n_obs_b
	 m = n 
	 if(m .lt. 1) go to 10
c
c	 ** no climatological check on STN **       rely(01,m)=imissing 
c	 ** no climatological check on OBSTYPE **   rely(02,m)=imissing
c	 ** no climatological check on LAT_S **     rely(03,m)=imissing
c	 ** no climatological check on LON_S **     rely(04,m)=imissing
c	 ** no climatological check on ELEV_S **    rely(05,m)=imissing
c	 ** no climatological check on WX **        rely(06,m)=imissing
c
	 if(t_s(n)    .ge.  -50.          .and.
     &      t_s(n)    .le.  130.               ) rely(07,m)=10 
c
	 if(td_s(n)   .ge.  -50.          .and.
     &      td_s(n)   .le.   90.               ) rely(08,m)=10
c
	 if(dd_s(n)   .ge.    0.          .and.
     &      dd_s(n)   .le.  360.               ) rely(09,m)=10
c
	 if(ff_s(n)   .ge.    0.          .and.
     &      ff_s(n)   .le.  120.               ) rely(10,m)=10
c
	 if(ddg_s(n)  .ge.    0.          .and.
     &      ddg_s(n)  .le.  360.               ) rely(11,m)=10
c
	 if(ffg_s(n)  .ge.    0.          .and.
     &      ffg_s(n)  .le.  200.               ) rely(12,m)=10
c
	 if(pstn_s(n) .ge.  500.          .and.
     &      pstn_s(n) .le. 1100.               ) rely(13,m)=10
c
	 if(pmsl_s(n) .ge.  900.          .and.
     &      pmsl_s(n) .le. 1100.               ) rely(14,m)=10
c
	 if(alt_s(n)  .ge.  900.          .and.
     &      alt_s(n)  .le. 1100.               ) rely(15,m)=10 
c
c	 ** no climatological check on KLOUD_S **   rely(16,m)=imissing 
c	 ** no climatological check on HGT_CEIL **  rely(17,m)=imissing 
c	 ** no climatological check on HGT_LOW **   rely(18,m)=imissing 
c	 ** no climatological check on COVER_S **   rely(19,m)=imissing 
c	 ** no climatological check on SOLAR_S **   rely(20,m)=imissing 
c	 ** no climatological check on IDP3_S **    rely(21,m)=imissing 
c	 ** no climatological check on STORE_EMV ** rely(22,m)=imissing 
c	 ** no climatological check on STORE_AMT ** rely(23,m)=imissing
c	 ** no climatological check on STORE_HGT ** rely(24,m)=imissing
c
	 if(vis_s(n)  .ge.    0.          .and.
     &      vis_s(n)  .le.  200.               ) rely(25,m) = 10
c
c	 ** no climatological check on OBSTIME **   rely(26,m)=imissing
c
 10	continue
c
	do 12 n=1,n_obs_b_l
	 m = ivals_l(n)
	 if (m .lt. 1) go to 12
c
c	 ** no climatological check on STN_L **      rely_l(01,m)=imissing 
c	 ** no climatological check on OBSTYPE_L **  rely_l(02,m)=imissing
c	 ** no climatological check on LAT_L **      rely_l(03,m)=imissing
c	 ** no climatological check on LON_L **      rely_l(04,m)=imissing
c	 ** no climatological check on ELEV_L **     rely_l(05,m)=imissing
c	 ** no climatological check on WX_L **       rely_l(06,m)=imissing
c
	 if(t_l(n)        .ge.   -50.      .and.
     &      t_l(n)        .le.   130.           ) rely_l(07,m)=10 
	 if(td_l(n)       .ge.   -50.      .and.
     &      td_l(n)       .le.    90.           ) rely_l(08,m)=10
	 if(dd_l(n)       .ge.     0.      .and.
     &      dd_l(n)       .le.   360.           ) rely_l(09,m)=10
	 if(ff_l(n)       .ge.     0.      .and.
     &      ff_l(n)       .le.   120.           ) rely_l(10,m)=10
	 if(ddg_l(n)      .ge.     0.      .and.
     &      ddg_l(n)      .le.   360.           ) rely_l(11,m)=10
	 if(ffg_l(n)      .ge.     0.      .and.
     &      ffg_l(n)      .le.   200.           ) rely_l(12,m)=10
	 if(pstn_l(n)     .ge.   500.      .and.
     &      pstn_l(n)     .le.  1100.           ) rely_l(13,m)=10
	 if(pmsl_l(n)     .ge.   900.      .and.
     &      pmsl_l(n)     .le.  1100.           ) rely_l(14,m)=10
	 if(alt_l(n)      .ge.   900.      .and.
     &      alt_l(n)      .le.  1100.           ) rely_l(15,m)=10 
c
c	 ** no climatological check on KLOUD_L **     rely_l(16,m)=imissing 
c	 ** no climatological check on CEIL_L **      rely_l(17,m)=imissing 
c	 ** no climatological check on LOWCLD_L **    rely_l(18,m)=imissing 
c	 ** no climatological check on COVER_L **     rely_l(19,m)=imissing 
c	 ** no climatological check on RAD_L **       rely_l(20,m)=imissing 
c	 ** no climatological check on IDP3_L **      rely_l(21,m)=imissing 
c	 ** no climatological check on STORE_EMV_L ** rely_l(22,m)=imissing 
c	 ** no climatological check on STORE_AMT_L ** rely_l(23,m)=imissing
c	 ** no climatological check on STORE_HGT_L ** rely_l(24,m)=imissing
c
	 if(vis_l(n)      .ge.     0.      .and.
     &      vis_l(n)      .le.   200.           ) rely_l(25,m) = 10 
c 
c        ** no climatological check on OBSTIME_L ** rely_l(26,m)=imissing c
 12	continue
c
	write(60,*) '  '
	write(60,*) ' --- Reliability after climatalogical test --- '
	write(60,*) '  '
	do 100 j=1,n_obs_curr
	 write(60,900) j, stn(j), (rely(i,j),i=1,26)
 100	continue
 900	format(i5,a5,15i4,/,10x,11i4)
c
c.....  standard deviation check
c
        if(l_dev_ck)then
	  call time_ck(stn,n_obs_b,stn_l,n_obs_b_l,ivals)
cc	  call dev_ck( 5, n_meso_pos, n_obs_b, elev_s, rely, 
cc     +	       n_obs_b_l, elev_l, rely_l, ivals)
	  call dev_ck( 7, n_meso_pos, n_obs_b, t_s, rely, 
     +	       n_obs_b_l, t_l, rely_l, ivals)
	  call dev_ck( 8, n_meso_pos, n_obs_b, td_s, rely, 
     +	       n_obs_b_l, td_l, rely_l, ivals)
	  call dev_ck(13, n_meso_pos, n_obs_b, pstn_s, rely, 
     +	       n_obs_b_l, pstn_l, rely_l, ivals)
	  call dev_ck(14, n_meso_pos, n_obs_b, pmsl_s, rely, 
     +	       n_obs_b_l, pmsl_l, rely_l, ivals)
	  call dev_ck(15, n_meso_pos, n_obs_b, alt_s, rely, 
     +	       n_obs_b_l, alt_l, rely_l, ivals)
        endif
c
c.....  model background checks
c
        c_field = 'T'
        ifld = 7
        thresh = thresh_t

        call bkg_chk(iback_t,c_field,ifld,thresh,t_s
     1              ,ii,jj,t_bk_f,rely,badflag,n_obs_curr
     1              ,stn,mxstn,ni,nj,istatus)

        c_field = 'TD'
        ifld = 8
        thresh = thresh_td

        call bkg_chk(iback_t,c_field,ifld,thresh,td_s
     1              ,ii,jj,td_bk_f,rely,badflag,n_obs_curr
     1              ,stn,mxstn,ni,nj,istatus)

        c_field = 'MSLP'
        ifld = 14
        thresh = thresh_mslp

        call bkg_chk(iback_mp,c_field,ifld,thresh,pmsl_s
     1              ,ii,jj,mslp_bk_mb,rely,badflag,n_obs_curr
     1              ,stn,mxstn,ni,nj,istatus)

c
	write(60,*) '  '
	write(60,*) ' --- Reliability after std dev & bkg tests --- '
	write(60,*) '  '
	do 101 j = 1,n_obs_curr
	 write(60,900) j, stn(j), (rely(i,j),i=1,26)
 101	continue

c
c.....  qc finished
c
	jstatus = 1
 999	write(60,*) ' ** qc of surface data complete. ** '
        write(6,*)  ' ** qc of surface data complete. ** '
c
	return
	end

c ---------------------------------------------------------------------------
c
c

        subroutine bkg_chk(iback,c_field,ifld,thresh,ob_s
     1              ,ii,jj,background_a,rely,badflag,n_obs_curr
     1              ,stn,mxstn,ni,nj,istatus)

 	real ob_s(mxstn), background_a(ni,nj)
        real badflag
	integer rely(26,mxstn), ii(mxstn), jj(mxstn)
        character*(*) c_field
	character stn(mxstn)*3

        if(iback .eq. 1)then
            write(6,*)' BKG '//c_field//' check'
            do j = 1,n_obs_curr
                observation = ob_s(j)

                if(observation .ne. badflag)then
                    ista = ii(j)
                    jsta = jj(j)

                    if(ista .ge. 1 .and. ista .le. ni .and. 
     1                 jsta .ge. 1 .and. jsta .le. nj)then ! Inside the domain
                        background = background_a(ista,jsta)
                        diff = observation - background

                        if(abs(diff) .le. thresh)then
!                           write(6,911)j, stn(j), observation, 
!    1                                  background, diff       
 911                        format('QC: ',i5,1x,a3,1x,3f9.1)     

                        else
                            write(6,912)j, stn(j), observation, 
     1                                  background, diff       
 912                        format('QC: ',i5,1x,a3,1x,3f9.1
     1                            ,' Exceeds threshold')     

                            rely(ifld,j) = -25 ! Set to reject

                        endif

                    endif ! Inside the domain

                endif ! Station has missing data

            enddo ! j

        else
            write(6,*)' Skipping BKG '//c_field//' check'

        endif

        return
        end

c ---------------------------------------------------------------------------
c
c
	subroutine time_ck(stn1,num1,stn2,num2,ivals)

c       Returns array of indices in stn2 array corresponding to each index
c       in the stn1 array. This matches the stations in the two arrays.

	dimension ivals(num1)
	character stn1(num1)*3,stn2(num2)*3
	do 10 n=1,num1
	 ivals(n) = -99
	 do 20 m=1,num2
	  if(stn2(m) .ne. stn1(n)) goto 20
	  ivals(n) = m
20	 continue
10	continue	 
	return
	end
c ---------------------------------------------------------------------------
c
c
	subroutine time_ck2(stn1,num1,stn2,num2,ivals_l,
     &	                    stn_all,n_obs_curr)
	dimension ivals_l(num2)
	character stn1(num1)*3,stn2(num2)*3,stn_all(n_obs_curr)*3

!       Calculate array of indices in stn_all array corresponding to each 
!       index in the stn2 array. This matches the stations in the two arrays.
	do 30 n=1,num2
	 ivals_l(n) = -99
	 do 40 m=1,n_obs_curr
	  if(stn_all(m) .ne. stn2(n)) goto 40
	  ivals_l(n) = m
40	 continue
	 if(ivals_l(n) .lt. 0) 
     &          write(60,*) n,' cannot find previous ',stn2(n)
30	continue	 

	return
	end
c -----------------------------------------------------------------------------
c
c
	subroutine dev_ck(ifld,n_meso_g,n_obs_b,  aa_s,rely,  
     +	            n_obs_b_l,aa_l,rely_l,ivals)
	real aa_s(n_obs_b),aa_l(n_obs_b_l)
	real diff(n_obs_b),stdev(n_obs_b)  !work arrays
	integer rely(26,n_obs_b)  
	integer rely_l(26,n_obs_b),ivals(n_obs_b), qc
	missing = -99.
	imissing = -99
c
c	compute avg change from prev hr (omit mesonet pressures)
c
	sumdiff = 0.
	numdiff = 0
	do 20 n=1,n_obs_b
	 diff(n) = missing
	 if(ifld .eq. 13 .and. n .le. n_meso_g) goto 20
	 m  = ivals(n)
	 if(m .eq. imissing) goto 20
	 if(aa_s(n) .le. missing .or. aa_l(m) .le. missing) goto 20
	 diff(n) = aa_s(n) - aa_l(m)
	 numdiff = numdiff + 1
	 sumdiff = sumdiff + diff(n)
 20	continue
	if (numdiff .eq. 0) then
	 write(60,*) ' numdiff=0...no dev_ck done.'
	 return
	else
	 avgdiff = sumdiff/float(numdiff)
	end if
c
c	compute std dev of each stn value from avg value for field
c
	sumdev = 0.
	do 210 n=1,n_obs_b
	 if(diff(n) .ne. missing) sumdev = sumdev + (diff(n)-avgdiff)**2
 210	continue
	std = sqrt(sumdev/float(numdiff))
        if(std .eq. 0.) then
	 write(60,*) ' std dev = 0., dev_ck not completed...'
	 return
	endif
c
	do 220 n=1,n_obs_b
	 stdev(n) = missing
	 if(diff(n) .ne. missing) stdev(n) = diff(n)/std
 220	continue
c
c	add or subtract reliability points based on std dev
c
	do 230 n=1,n_obs_b
	 m1 = n 
	 if(m1 .lt. 1) go to 230
	 qc = +25
	 if(abs(stdev(n)) .ge. 5.0 .and. abs(diff(n)) .gt. 5.0) qc = -25
	 if(diff(n) .ne. missing) rely(ifld,m1) = rely(ifld,m1) + qc
 230	continue
	write(60,*) ' subr dev_ck complete ', ifld
	return
	end
