c
c
       Program QC_main
c
       include 'lapsparms.cmn'
       include 'surface_obs.inc'
       character laps_domain*9
c
       laps_domain = 'nest7grid'
       call get_laps_config(laps_domain,istatus)
       if(istatus .ne. 1) then
          write(6,*) 'QC_main: Error getting domain dimensions'
          stop
       endif
c
       call QC_main_sub(NX_L_CMN,NY_L_CMN,maxstns_cmn,
     &                  laps_cycle_time_cmn,badflag)
c     
       end
c
c
       subroutine QC_main_sub(ni, nj, m, laps_cycle_time, badflag)
c
c*********************************************************************
c
c     LAPS QC driver...does quality control on surface observations
c     read in from the raw LS2 files.
c     
c     Original: John McGinley, NOAA/FSL  Spring 1998 (QC_kalmanVII)
c     Changes:
c       25 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       07 Oct 1998  Peter Stamus, NOAA/FSL
c          Finished 1st set of changes.  Change 'ran' to 'ran1' for 
c          portability.
c
c     Notes:
c
c*********************************************************************
c
c  set up tabled exponential and freq used constants
c
        common tab(10000),pi,re,rdpdg,reorpd
        parameter( nvar = 10 )
c
c  arrays for reading and converting observations
c
        real*4   lata (m),lona(m),eleva(m)
        real*4   latb (m),lonb(m),elevb(m)
        real*4   pstn(m),pmsl(m),alt(m)
        real*4   pstna(m),pmsla(m),alta(m)
        real*4   pstnb(m),pmslb(m),altb(m)
        real*4   pmslc(m),altc(m)
        real*4   pmsle(m),alte(m)
        real*4   pmslf(m),altf(m)
        real*4   u(m),v(m)
c
c ....lso reader arrays
c
        real*4 lat(m), lon(m), elev(m)
        real*4 t(m), t_ea(m), max24t(m), min24t(m)
        real*4 td(m), td_ea(m), rh(m), rh_ea(m)
        real*4 dd(m), ddg(m), dd_ea(m)
        real*4 ff(m), ffg(m), ff_ea(m)
        real*4 alt_ea(m), delp(m)
        real*4 p_ea(m)
        real*4 vis(m), vis_ea(m)
        real*4 solar(m), solar_ea(m)
        real*4 sfct(m), sfct_ea(m)
        real*4 sfcm(m), sfcm_ea(m)
        real*4 pcp1(m), pcp3(m), pcp6(m), pcp24(m)
        real*4 snow(m), snow_ea(m), pcp_ea(m)
        real*4 store_cldht(m,5)
c
        integer*4 i4time, wmoid(m), jstatus
        integer*4 time(m), delpch(m), kkk_s(m)
c
        character infile*256, store_cldamt(m,5)*4
        character stations(m)*20, provider(m)*11
        character reptype(m)*6, autostntype(m)*6
c
c.....  LSO write arrays
c
        real*4  store_1(m,4), 
     &          store_2(m,3), store_2ea(m,3),
     &          store_3(m,4), store_3ea(m,2),
     &          store_4(m,5), store_4ea(m,2),
     &          store_5(m,4), store_5ea(m,4),
     &          store_6(m,5), store_6ea(m,2),
     &          store_7(m,3)
c
c  two sets of past obs a and b
c
        real*4   ta(m),tda(m),dda(m),ffa(m),ua(m),va(m)
        real*4   tb(m),tdb(m),ddb(m),ffb(m),ub(m),vb(m)
        character providera(m)*11, providerb(m)*11
        character reptypea(m)*6, reptypeb(m)*6
c
c  utility arrays
c
        real*4   tc(m),tdc(m),uc(m),vc(m)
        real*4   te(m),tde(m),ue(m),ve(m)
        real*4   tf(m),tdf(m),uf(m),vf(m)
c
c  arrays for Kalman
c     error covariances for each variable
c
        real*4   w(m,m),vv(m,m),pt(m,m),ptd(m,m),pu(m,m),pv(m,m)
        real*4   ppm(m,m),pal(m,m)
c
c  input arrays-model and analysis weighting arrays
c
        real dt(m),dta(m),dtd(m),dtda(m)
        real du(m),dua(m),dv(m),dva(m),dpm(m),dpma(m),dalt(m),dalta(m)
        real F(m,m),mwt(m,m),dwt(m,m)
c
c  output arrays
c
        real*4 xt(m),xtt(m),xtd(m),xttd(m),xu(m),xtu(m),xv(m),
     &     xtv(m),xal(m),xtal(m),xpm(m),xtpm(m)
c
c  fourier and observation estimates 
c
        real   monster(m,m,nvar),yta(m),ytda(m),yua(m)
        real*4   yva(m),ypmsla(m),yalta(m),fcf(m,m,nvar)
c
c   arrays for error processing for each variable by iteration
c
        real ar(m),wr(m),vr(m),wit(m,m),witd(m,m),witu(m,m),witv(m,m)
        real witpm(m,m),wital(m,m),vit(m,m),vitd(m,m),vitu(m,m)
        real vitv(m,m),vitpm(m,m),vital(m,m)
c
c theta arrays for times a and b
c
        real theta(m),thetaa(m),thetab(m),thetaf(m)
c
c utility arrays
c
        real b(m),c(m,m),zot(m)
c
c error variables
c
        real oberr,moderr
c
c integer arrays for data input
c       Integer*4   obstime(m),kloud(m),idp3(m),mm
c quality control status for each variable
c
        integer qcstat(m),qcstatd(m),qcstapm(m),qcstal(m),qcstauv(m)
        integer qcd(m)
c
c switches 
c
        integer on,off
        logical exists, flagstart
        data exists/.false./
        data flagstart/.false./
c        
c character arrys for file names, station names, time, wx symbols
c

	integer     i4prev1, i4prev2, cycle
	character*9 filename, fname1, fname2
        Character   atime*24, atime_cur*24, stn(m)*5, wx(m)*25
        character   stna(m)*5, stnb(m)*5
        character   monfile*256, atime_mon*24
        character   dir_s*180, dir_out*256, ext_out*31
c
c
c..... Start here.  First get the time from the user or scheduler.
c
        narg = iargc()
c
        if(narg .eq. 0) then
           call get_systime(i4time,filename,istatus)
        else
           write(6,973)
 973       format(' Enter input filename (yyjjjhhmm): ',$)
           read(5,972) filename
 972       format(a)
           call i4time_fname_lp(filename,i4time,status)
        endif
c
c.....  Set the data cycle in seconds and figure out the 
c.....  previous time variables and filenames.
c
        cycle = laps_cycle_time
	i4prev1 = i4time - cycle
	i4prev2 = i4prev1 - cycle
        call make_fnam_lp(i4prev1, fname1, istatus)
        call make_fnam_lp(i4prev2, fname2, istatus)
c
c  observation and model error assumptions
c *** expected obs accuracy for each variable at each station
c      is calculated in the LSO file, and is passed in by
c      the read_surface_data routine.  Can use these someday.***
c
        oberr=1.0*9./5.!C error converted to F
        oberrt=1.0*9./5.
        oberrtd=2.0*9./5.
        oberru=1.0/.515!m/sec error conv to knts
        oberrpm=1.0
        oberral=1.0
c
        moderr=4.0
        moderrt=4.0
        moderrtd=4.0
        moderru=4.0
        moderrpm=2.0
        moderral=2.0
        grosst=25.
        grosstd=30.
        grossuv=25.
        grosspm=20.
        grossal=20.0
c
c  offset for model operator (>> than ob-fg)
c
        offset=1000.
        on=1
        off=0
c
c  constants
c
        re=6371122.
        pi = 4.0 * atan(1.0)
        rdpdg=pi/180.
        reorpd=re*rdpdg  
c
c  averaging period for error for w and vv
c
        ia=24
c
c  set up exponential lookup table
c
        do i=1,10000
           tab(i)= exp (-float((i-1)*(i-1))/1000.)
        enddo !i
c
c.....  Check for a 'monster' file in the 'LSQ' output
c.....  directory.  If there's one there, read it.
c.....  If not, fill the error arrays here.
c
        call get_directory('lsq', infile, len)
        infile = infile(1:len) // 'monster.dat'
        call s_len(infile, len)
        INQUIRE(FILE=infile(1:len), EXIST=exists)
        if(.not. exists) then
           flagstart = .true.  !no monster file
        else 
           flagstart = .false. !yes monster file
        endif
c
 500    continue
        call get_directory('lsq', dir_s, len)
        monfile = dir_s(1:len) // 'monster.dat'

        IF(flagstart) THEN
c
c.....  No monster file, initilize error arrays
c
           call zero(w,m,m)
cc           call zero(p,m,m)  !**** no p array...???
           call zero(pt,m,m)
           call zero(ptd,m,m)
           call zero(pu,m,m)
           call zero(pv,m,m)
           call zero(ppm,m,m)
           call zero(pal,m,m)
           call zero(vv,m,m)

           do i=1,m
              zot(i)=0. 
              do j=1,m
              do k=1,nvar
                 monster(i,j,k)=badflag
              enddo !k
              enddo !j
              do j=1,2
                 iiiii=i*j*k
c
c random number intial seed
c              
                 iiiii=ran1(iiiii)*2000000.-1000000.
                 vit(i,j)=oberrt*ffz(iiiii,20)
                 wit(i,j)=moderr*ffz(iiiii,20)
                 witd(i,j)=moderr*ffz(iiiii,20)
                 vitd(i,j)=oberrtd*ffz(iiiii,20)
                 witu(i,j)=moderr*ffz(iiiii,20)
                 witv(i,j)=moderr*ffz(iiiii,20)
                 vitv(i,j)=oberru*ffz(iiiii,20)
                 vitu(i,j)=oberru*ffz(iiiii,20)
                 witpm(i,j)=moderr*ffz(iiiii,20)
                 vitpm(i,j)=oberrpm*ffz(iiiii,20)
                 wital(i,j)=moderr*ffz(iiiii,20)
                 vital(i,j)=oberral*ffz(iiiii,20)
              enddo !on j
              Pt(i,i)=moderrt*moderr
              Ptd(i,i)=moderrtd*moderr
              Pu(i,i)=moderru*moderr
              Pv(i,i)=moderru*moderr
              Ppm(i,i)=moderrpm*moderr
              Pal(i,i)=moderral*moderr
           enddo !i
           on=1
           off=0
           it=1 
c
c.....  Get the LSO data file for 2 cycles ago.
c
           call read_surface_data(i4prev2,atime,n_obs_g,n_obs_b,time,
     &  wmoid,stations,providerb,wx,reptypeb,autostntype,latb,lonb,
     &  elevb,tb,tdb,rh,ddb,ffb,ddg,ffg,altb,pstnb,pmslb,delpch,delp,
     &  vis,solar,sfct,sfcm,pcp1,pcp3,pcp6,pcp24,snow,kkk_s,max24t,
     &  min24t,t_ea,td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,
     &  sfct_ea,sfcm_ea,pcp_ea,snow_ea,store_cldamt,store_cldht,
     &  m,jstatus)
c
           if(jstatus .ne. 1) then
              print *,
     &           ' No LSO data for (Hr-2). No Kalman QC for ',filename
              stop
           else
              print *,' Found LSO data (Hr-2) at ', atime
           endif
c
           maxstab=n_obs_b
           do k=1,maxstab
              stnb(k)(1:5)=stations(k)(1:5)
              print*,stnb(k) 
           enddo !k
c
           print*, maxstab, 'obs read for cycle '    
c
           it=2
c
c.....  Now get LSO data for 1 cycle ago.
c
           call read_surface_data(i4prev1,atime,n_obs_g,n_obs_b,time,
     &  wmoid,stations,providera,wx,reptypea,autostntype,lata,lona,
     &  eleva,ta,tda,rh,dda,ffa,ddg,ffg,alta,pstna,pmsla,delpch,delp,
     &  vis,solar,sfct,sfcm,pcp1,pcp3,pcp6,pcp24,snow,kkk_s,max24t,
     &  min24t,t_ea,td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,
     &  sfct_ea,sfcm_ea,pcp_ea,snow_ea,store_cldamt,store_cldht,
     &  m,jstatus)
c
           if(jstatus .ne. 1) then
              print *,
     &          ' No LSO data for (Hr-1). No Kalman QC for ',filename
              stop
           else
              print *,' Found LSO data for (Hr-1) at ', atime
           endif
c
           maxstaa=n_obs_b
           do k=1,maxstaa
              stna(k)(1:5)=stations(k)(1:5)
           enddo !k
           print*, maxstaa, 'obs read for cycle '
           call reorder(ta,tda,dda,ffa,lata,lona,eleva,pstna,pmsla,alta,
     &       stna,providera,reptypea,
     &       tb,tdb,ddb,ffb,latb,lonb,elevb,pstnb,pmslb,altb,stnb,
     &       providerb,reptypeb,maxstaa,maxstab,m,badflag)
           call convuv(ddb,ffb,ub,vb,maxstab,m,badflag)
           call convuv(dda,ffa,ua,va,maxstaa,m,badflag)
c
c  convert to theta from raw temps
c
           call thet(tb,thetab,maxstaa,eleva,m)
c
c  replace missing values of theta from neighbors
c
           call fill(stnb,thetab,latb,lonb,elevb,thetab,sl,
     &       maxstaa,m,8./1000.,-1./1000. )         
           call thet2T(tb,thetab,maxstaa,eleva,m)
           call writemon(tb,tdb,ub,vb,pmslb,altb,nvar,maxstaa,m,
     &          monster,(it  )) 
           call thet(ta,thetaa,maxstaa,eleva,m)
           call fill(stna,thetaa,lata,lona,eleva,thetaa,sl,
     &             maxstaa,m,8./1000.,-1./1000.)          
c
c  now for first guess ta, use thetas and interp to stn locations
c
           call thet2T(ta,thetaa,maxstaa,eleva,m)
c
c same thing for dewpoint and winds
c
           call fill(stna,tda,lata,lona,eleva,thetaa,vl,
     &          maxstaa,m,15./1000.,-35./1000.)        
           call fill(stna,ua,lata,lona,eleva,thetaa,vl,
     &          maxstaa,m,10./1000.,-10./1000.)        
           call fill(stna,va,lata,lona,eleva,thetaa,vl,
     &       maxstaa,m,10./1000.,-10./1000.)         
           call fill(stna,pmsla,lata,lona,eleva,thetaa,vl,
     &       maxstaa,m,2./1000.,-2./1000.)        
           call fill(stna,alta,lata,lona,eleva,thetaa,vl,
     &       maxstaa,m,2./1000.,-2./1000.)         
           icnt=1
c
        ELSE


           open(15,file=monfile(1:len),
     &         form='unformatted',status='old')
           read(15) atime_mon, i4time_mon
c
c.....  Check the time in the monster file.  If its too old,
c.....  bag it and start over.
c
           print *,' Found monster file valid for: ', atime_mon
           if( abs(i4time-i4time_mon) .gt. laps_cycle_time) then
              print *,'    But monster file too old. Starting over...'
              flagstart = .true.
              close(15)
              go to 500
           endif
c
           read(15) vit,vitd,vitu,vitv,vitpm,vital
           read(15) wit,witd,witu,witv,witpm,wital
           read(15) pt,ptd,pu,pv,ppm,pal
           read(15) monster
           read(15) ta,tda,ua,va,pmsla,alta
           close(15)
        ENDIF
c
c start of routine
c
        maxstab=maxstaa
        badthr=3.5
c
c put the obs in the monster stack
c
        call writemon(ta,tda,ua,va,pmsla,alta,nvar,maxstaa,m,
     &          monster,(it  )) 



	go to 1234


        if(it.ge.26) then
           icyc=24
           nn=32
        else
           go to 334
        endif
        if(it.ge.66) then
           icyc=48
           nn=64
        endif
c
c     perform fourier analysis on the last 24 or 48 obs
c
 333    call fouranal(monster,fcf,maxstaa,m,6,nvar,nn,icyc) !1 is maxvar
c
c     use seriesn to estimate next ob/or use taylor series
c
 334    continue
        call project(fcf,yta,monster,1,nvar,
     &          maxstaa,m,nn,atime,it,icyc,oberrt,badthr)
        call project(fcf,ytda,monster,2,nvar,
     &          maxstaa,m,nn,atime,it,icyc,oberrtd,badthr)
        call project(fcf,yua,monster,3,nvar,
     &          maxstaa,m,nn,atime,it,icyc,oberru,badthr)
        call project(fcf,yva,monster,4,nvar,
     &          maxstaa,m,nn,atime,it,icyc,oberru,badthr)
        call project(fcf,ypmsla,monster,5,nvar,
     &          maxstaa,m,nn,atime,it,icyc,oberrpm,badthr)
        call project(fcf,yalta,monster,6,nvar,
     &          maxstaa,m,nn,atime,it,icyc,oberral,badthr)
c     
c  input ob estimates as a departure from background ta
c
        call perturb(yta,ta,yta,maxstaa,m,offset,on)      
        call perturb(ta,ta,dta,maxstaa,m,offset,on)      
        call perturb(ytda,tda,ytda,maxstaa,m,offset,on)      
        call perturb(tda,tda,dtda,maxstaa,m,offset,on)      
        call perturb(yua,ua,yua,maxstaa,m,offset,on)      
        call perturb(ua,ua,dua,maxstaa,m,offset,on)      
        call perturb(yva,va,yva,maxstaa,m,offset,on)      
        call perturb(va,va,dva,maxstaa,m,offset,on)      
        call perturb(ypmsla,pmsla,ypmsla,maxstaa,m,offset,on)      
        call perturb(pmsla,pmsla,dpma,maxstaa,m,offset,on)      
        call perturb(yalta,alta,yalta,maxstaa,m,offset,on)      
        call perturb(alta,alta,dalta,maxstaa,m,offset,on)      
c     
c compute distance/wind/theta defined weighting functions for error cov
c
        call weights(mwt,dwt,ua,va,thetaa,lata,lona,
     &       cycle,maxstaa,m)
c
        write(*,*) 'Kalman for TEMP'
*****MODEL IS UNDER REVISION - ARGUEMENTS WILL CHANGE**********
        call model(F,dta,thetaa,ua,va,thetaa,mwt,lata,lona,eleva
     &          ,maxstaa,m,ni,nj,timer,cycle,14.,4.,24.,6.)
c
c  create the v and w matrices using the error history
c
        Call avgdiagerr(wit,B,c,W,maxstaa,ia,m,it-1)
        Call avgdiagerr(Vit,B,c,vV,maxstaa,ia,m,it-1)
c  kalman....
        call kalman(dta,F,YTA,pt,it, w,vv,xt,xtt,maxstaa,m,atime)
        write(*,*) 'Kalman for DEWPOINT'
        call model(F,dtda,tda,ua,va,thetaa,mwt,lata,lona,eleva
     1          ,maxstaa,m,ni,nj,timer,cycle,4.,2.,12.,18.)
c
c  create the v and w matrices using the error history
c
        Call avgdiagerr(witd,B,c,W,maxstaa,ia,m,it-1)
        Call avgdiagerr(Vitd,B,c,vV,maxstaa,ia,m,it-1)
        call kalman(dtda,F,YTDA,ptd,it,w,vv,xtd,xttd,maxstaa,m,atime)
        write(*,*) 'Kalman for U-WIND'
        call model(F,dua,thetaa,ua,va,thetaa,mwt,lata,lona,eleva
     1          ,maxstaa,m,ni,nj,timer,cycle,2.,0.,20.,6.)
c
c  create the v and w matrices using the error history
c
        Call avgdiagerr(witu,B,c,W,maxstaa,ia,m,it-1)
        Call avgdiagerr(Vitu,B,c,vV,maxstaa,ia,m,it-1)
c  kalman....
        call kalman(dua,F,yua,pu,it, w,vv,xu,xtu,maxstaa,m,atime)       
        write(*,*) 'Kalman for V-WIND'
        call model(F,dva,thetaa,ua,va,thetaa,mwt,lata,lona,eleva
     1          ,maxstaa,m,ni,nj,timer,cycle,0.,0.,24.,6.)
c
c  create the v and w matrices using the error history
c
        Call avgdiagerr(witv,B,c,W,maxstaa,ia,m,it-1)
        Call avgdiagerr(Vitv,B,c,vV,maxstaa,ia,m,it-1)
c
c  kalman....
c
        call kalman(dva,F,yva,pv,it, w,vv,xv,xtv,maxstaa,m,atime)
        write(*,*) 'Kalman for PMSL'
        call model(F,dpma,thetaa,ua,va,thetaa,mwt,lata,lona,eleva
     1          ,maxstaa,m,ni,nj,timer,cycle,0.,3.,16.,16.)
c
c  create the v and w matrices using the error history
c
        Call avgdiagerr(witpm,B,c,W,maxstaa,ia,m,it-1)
        Call avgdiagerr(Vitpm,B,c,vV,maxstaa,ia,m,it-1)
c  kalman....
        call kalman(dpma,F,YpmslA,ppm,it, w,vv,xpm,xtpm,maxstaa,m,atime)
        write(*,*) 'Kalman for ALSTG'
        call model(F,dalta,thetaa,ua,va,thetaa,mwt,lata,lona,eleva
     1          ,maxstaa,m,ni,nj,timer,cycle,0.,.1,24.,16.)
c
c  create the v and w matrices using the error history
c
        Call avgdiagerr(wital,B,c,W,maxstaa,ia,m,it-1)
        Call avgdiagerr(Vital,B,c,vV,maxstaa,ia,m,it-1)
c  kalman....
        call kalman(dalta,F,yalta,pal,it, w,vv,xal,xtal,maxstaa,m,atime)
c
c read in new obs.....maxstaa will go ue
c



1234    continue



        call read_surface_data(i4time,atime_cur,n_obs_g,n_obs_b,time,
     &     wmoid,stations,provider,wx,reptype,autostntype,lat,lon,elev,
     &     t,td,rh,dd,ff,ddg,ffg,alt,pstn,pmsl,delpch,delp,vis,solar,
     &     sfct,sfcm,pcp1,pcp3,pcp6,pcp24,snow,kkk_s,max24t,min24t,t_ea,
     &     td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,sfct_ea,
     &     sfcm_ea,pcp_ea,snow_ea,store_cldamt,store_cldht,
     &     m,jstatus)
c
        if(jstatus .ne. 1) then
           print *,' No LSO data for ', filename,'.  Stopping'
           stop
        else
           print *,' Found LSO data (current) at ', atime_cur
        endif
        atime = atime_cur
c
        maxsta=n_obs_b
        do k=1,maxsta
           stn(k)(1:5)=stations(k)(1:5)
        enddo !k
        print*, maxsta,'obs read for cycle ',it
        write(*,1066) atime
 1066   format(1x,a24)
        call reorder(t,td,dd,ff,lat,lon,elev,pstn,pmsl,alt,
     &          stn,provider,reptype,
     &          ta,tda,dda,ffa,lata,lona,eleva,pstna,pmsla,alta,
     &          stna,providera,reptypea,maxsta,maxstaa,m,badflag)
        print*, maxsta,maxstaa,maxstab
c
c redo the theta conversion, interpolation, and back to ta for new stns
c
        call thet(ta,thetaa,maxstaa,eleva,m)
        write(*,*) 'FILL after new  read'
        call fill(stna,thetaa,lata,lona,eleva,thetaa,sl,
     &       maxstaa,m,8./1000.,-1./1000.)          
        call thet2T(ta,thetaa,maxstaa,eleva,m)
        call fill(stna,tda,lata,lona,eleva,thetaa,sl,
     &          maxstaa,m,8./1000.,-1./1000.)          
        call fill(stna,ua,lata,lona,eleva,thetaa,sl,
     &          maxstaa,m,10./1000.,-10./1000.)          
        call fill(stna,va,lata,lona,eleva,thetaa,sl,
     &          maxstaa,m,10./1000.,-10./1000.)          
        call fill(stna,pmsla,lata,lona,eleva,thetaa,sl,
     &          maxstaa,m,2./1000.,-2./1000.)          
        call fill(stna,alta,lata,lona,eleva,thetaa,sl,
     &          maxstaa,m,2./1000.,-2./1000.)          
c
c convert new winds to u,v
c
        call convuv(dd,ff,u,v,maxstaa,m,badflag)
c
c set qc status and gross error check
c

	go to 12345



        write(*,*) 'gross err T'
        call qcset(t,ta,qcstat,m,maxstaa,badflag,grosst)
        write(*,*) 'gross err TD'
        call qcset(td,tda,qcstatd,m,maxstaa,badflag,grosstd)
        write(*,*) 'gross err U'
        call qcset(u,ua,qcstauv,m,maxstaa,badflag,grossuv)
        write(*,*) 'gross err V'
        call qcset(v,va,qcstauv,m,maxstaa,badflag,grossuv)
        write(*,*) 'gross err PMSL'
        call qcset(pmsl,pmsla,qcstapm,m,maxstaa,badflag,grosspm)
        write(*,*) 'gross err ALT'
        call qcset(alt,alta,qcstal,m,maxstaa,badflag,grossal)
c     
c find observed departure from background 
c
        call perturb(t,ta,dt,maxstaa,m,0.,on)      
        call perturb(td,tda,dtd,maxstaa,m,0.,on)      
        call perturb(u,ua,du,maxstaa,m,0.,on)      
        call perturb(v,va,dv,maxstaa,m,0.,on)      
        call perturb(pmsl,pmsla,dpm,maxstaa,m,0.,on)      
        call perturb(alt,alta,dalt,maxstaa,m,0.,on)      
c     
c perform error analysis based on truth estimate write into history array
c remove offset from obestimates
c
        call perturb(yta,zot,yta,maxstab,m,offset,off)      
        call perturb(ytda,zot,ytda,maxstab,m,offset,off)      
        call perturb(xt,zot,xt,maxstab,m,offset,off)      
        call perturb(xtt ,zot,xtt ,maxstab,m,offset,off)      
        call perturb(xtd,zot,xtd,maxstab,m,offset,off)      
        call perturb(xttd ,zot,xttd,maxstab,m,offset,off)      
        call perturb(yua,zot,yua,maxstab,m,offset,off)      
        call perturb(yva,zot,yva,maxstab,m,offset,off)      
        call perturb(xu,zot,xu,maxstab,m,offset,off)      
        call perturb(xtu ,zot,xtu ,maxstab,m,offset,off)      
        call perturb(xv,zot,xv,maxstab,m,offset,off)      
        call perturb(xtv ,zot,xtv,maxstab,m,offset,off)      
        call perturb(ypmsla,zot,ypmsla,maxstab,m,offset,off)      
        call perturb(yalta,zot,yalta,maxstab,m,offset,off)      
        call perturb(xpm,zot,xpm,maxstab,m,offset,off)      
        call perturb(xtpm ,zot,xtpm ,maxstab,m,offset,off)      
        call perturb(xal,zot,xal,maxstab,m,offset,off)      
        call perturb(xtal ,zot,xtal,maxstab,m,offset,off)      
        write(*,*) 'ERROR PROCESSING FOR T'
        call errorproc(dt,yta,xtt,xt,wr,vr,ar,maxstab,m,qcstat
     &          ,oberrt,badthr,atime,icnt,1)  
        call writemon(wr,wr,wr,wr,wr,wr,1,maxstab,m,wit,it)
        call writemon(vr,vr,vr,vr,vr,vr,1,maxstab,m,vit,it)
        write(*,*) 'ERROR PROCESSING FOR TD'
        call errorproc(dtd,ytda,xttd,xtd,wr,vr,ar,maxstab,m,qcstatd  
     &          ,oberrtd,badthr,atime,icnt,2)  
        call writemon(wr,wr,wr,wr,wr,wr,1,maxstab,m,witd,it)
        call writemon(vr,vr,vr,vr,vr,vr,1,maxstab,m,vitd,it)
        write(*,*) 'ERROR PROCESSING FOR U'
        call errorproc(du,yua,xtu,xu,wr,vr,ar,maxstab,m,qcstauv
     &          ,oberru,badthr,atime,icnt,3)  
        call writemon(wr,wr,wr,wr,wr,wr,1,maxstab,m,witu,it)
        call writemon(vr,vr,vr,vr,vr,vr,1,maxstab,m,vitu,it)
        write(*,*) 'ERROR PROCESSING FOR v'
        call errorproc(dv,yva,xtv,xv,wr,vr,ar,maxstab,m,qcstauv  
     &          ,oberru,badthr,atime,icnt,4)  
        call writemon(wr,wr,wr,wr,wr,wr,1,maxstab,m,witv,it)
        call writemon(vr,vr,vr,vr,vr,vr,1,maxstab,m,vitv,it)
        write(*,*) 'ERROR PROCESSING FOR PMSL'
        call errorproc(dpm,ypmsla,xtpm,xpm,wr,vr,ar,maxstab,m,qcstapm
     &          ,oberrpm,badthr,atime,icnt,5)  
        call writemon(wr,wr,wr,wr,wr,wr,1,maxstab,m,witpm,it)
        call writemon(vr,vr,vr,vr,vr,vr,1,maxstab,m,vitpm,it)
        write(*,*) 'ERROR PROCESSING FOR ALSTG' 
        call errorproc(dalt,yalta,xtal,xal,wr,vr,ar,maxstab,m,qcstal  
     &          ,oberral,badthr,atime,icnt,6)  
        call writemon(wr,wr,wr,wr,wr,wr,1,maxstab,m,wital,it)
        call writemon(vr,vr,vr,vr,vr,vr,1,maxstab,m,vital,it)
c     
c  recover optimum fields   offset has been removed already 
c here a: model, b: ob, e: com 
c
        call perturb(te,ta,xt,maxstab,m,0.,off)
        call perturb(tde,tda,xtd ,maxstab,m,0.,off)
        call perturb(tb,ta,yta,maxstab,m,0.,off)
        call perturb(tdb,tda,ytda ,maxstab,m,0.,off)
        call perturb(ta,ta,xtt,maxstab,m,0.,off)
        call perturb(tda,tda,xttd ,maxstab,m,0.,off)
        call perturb(ue,ua,xu,maxstab,m,0.,off)
        call perturb(ve,va,xv ,maxstab,m,0.,off)
        call perturb(ub,ua,yua,maxstab,m,0.,off)
        call perturb(vb,va,yva ,maxstab,m,0.,off)
        call perturb(ua,ua,xtu,maxstab,m,0.,off)
        call perturb(va,va,xtv ,maxstab,m,0.,off)
        call perturb(pmsle,pmsla,xpm,maxstab,m,0.,off)
        call perturb(alte,alta,xal ,maxstab,m,0.,off)
        call perturb(pmslb,pmsla,ypmsla,maxstab,m,0.,off)
        call perturb(altb,alta,yalta ,maxstab,m,0.,off)
        call perturb(pmsla,pmsla,xtpm,maxstab,m,0.,off)
        call perturb(alta,alta,xtal ,maxstab,m,0.,off)
c     
c store optimum estimates in the "c" arrays
c c is now comb
c
        call replace(te,tc,maxstab,1,m,1)
        call replace(tde,tdc,maxstab,1,m,1)
        call replace(ue,uc,maxstab,1,m,1)
        call replace(ve,vc,maxstab,1,m,1)
        call replace(pmsle,pmslc,maxstab,1,m,1)
        call replace(alte,altc,maxstab,1,m,1)
c
c  QC check and ob replacement...we have ta(maxstaa),te(maxstab)
c     and t(maxstaa) 
c stage 1 ...buddy check. get value without using ob...fill one
c
        do i=1,maxstaa
           tf(i)=t(i)
           tdf(i)=td(i)
           uf(i)=u(i)
           vf(i)=v(i)
           pmslf(i)=pmsl(i)
           altf(i)=alt(i)
        enddo !i
c     
c stage 2....latest buddy values 
c
        write(*,*) 'FILL for F series'
        write(*,*) 't fill'
        call thet(tf,thetaf,maxstaa,eleva,m)
        call fillone(stna,thetaf,lata,lona,eleva,thetaf,sl,
     &       maxstaa,m,8./1000.,-1./1000.)          
        call thet2T(tf,thetaf,maxstaa,eleva,m)
c
c same thing for dewpoint and winds
c
        write(*,*) 'td fill'
        call fillone(stna,tdf,lata,lona,eleva,thetaa,vl,
     &       maxstaa,m,15./1000.,-35./1000.)        
        write(*,*) 'uf fill'
        call fillone(stna,uf,lata,lona,eleva,thetaa,vl,
     &          maxstaa,m,10./1000.,-10./1000.)        
        write(*,*) 'vf fill'
        call fillone(stna,vf,lata,lona,eleva,thetaa,vl,
     &          maxstaa,m,10./1000.,-10./1000.)         
        write(*,*) 'pmslf fill'
        call fillone(stna,pmslf,lata,lona,eleva,thetaa,vl,
     &          maxstaa,m,2./1000.,-2./1000.)        
        write(*,*) 'altf fill'
        call fillone(stna,altf,lata,lona,eleva,thetaa,vl,
     &          maxstaa,m,2./1000.,-2./1000.)         
c
c do buddy check flat
c
        do i=1,maxstaa
           if(t(i).eq.badflag) go to 70
           if(abs(tf(i)-t(i)).gt.badthr*oberrt) then
              qcstat(i)=qcstat(i)+100
              print*,'stn ',i,t(i),' fails buddy check temp ',tf(i)
              write(*,*) 
     &             'stn ',i,t(i),' fails buddy check temp ',tf(i)
           endif
 70        if(td(i).eq.badflag) go to 71
           if(abs(tdf(i)-td(i)).gt.badthr*oberrtd) then 
              qcstatd(i)=qcstatd(i)+100
              print*,'stn ',i,td(i),' fails buddy check dewp ',tdf(i)
              write(*,*) 
     &             'stn ',i,td(i),' fails buddy check dewp ',tdf(i)
           endif
 71        if(u(i).eq.badflag) go to 72
           if(abs(uf(i)-u(i)).gt.badthr*oberru) then
              qcstauv(i)=qcstauv(i)+100
              write(*,*) 
     &             'stn ',i,u(i),' fails buddy check u-wind ',uf(i)
           endif
           uf(i)=u(i)
 72        if(v(i).eq.badflag) go to 73
           if(abs(vf(i)-v(i)).gt.badthr*oberru) then 
              qcstauv(i)=qcstauv(i)+100
              print*,'stn ',i,v(i),' fails buddy check v-wind ',vf(i)
              write(*,*) 
     &             'stn ',i,v(i),' fails buddy check v-wind ',vf(i)
           endif
 73        if(pmsl(i).eq.badflag) go to 74
           if(abs(pmslf(i)-pmsl(i)).gt.badthr*oberrpm)then 
              qcstapm(i)=qcstapm(i)+100
              print*,
     &             'stn ',i,pmsl(i),' fails buddy check pmsl ',pmslf(i)
              write(*,*) 'stn ',i,pmsl(i),' fails buddy check pmsl ',
     &             pmslf(i)
           endif
           pmslf(i)=pmsl(i)
 74        if( alt(i).eq.badflag) go to 75
           if(abs(altf(i)-alt(i)).gt.badthr*oberral) then 
              qcstal(i)=qcstal(i)+100
              print*,
     &             'stn ',i, alt(i),' fails buddy check pmsl ', altf(i)
              write(*,*) 'stn ',i,alt (i),' fails buddy check pmsl ',
     &             altf(i)
           endif
 75     enddo !i
c     
c stage 3 determine what variable gets passed to ta for next cycle 
c
        do i=1,maxstab
           call flag(qcstat,qcd,m,6)
           if (t(i).eq.badflag) then
              te(i)=tc(i)
           else
              if(qcd(2).eq.1) te(i)=tc(i)
              if(qcstat(i).ne.111100) te(i)=t(i)
           endif
           call flag(qcstatd,qcd,m,6)
           if (td(i).eq.badflag) then
              tde(i)=tdc(i)
           else
              if(qcd(2).eq.1) tde(i)=tdc(i)
              if(qcstat(i).ne.111100) tde(i)=td(i)
           endif
           call flag(qcstauv,qcd,m,6)
           if (u(i).eq.badflag) then
              ue(i)=uc(i)
           else
              if(qcd(2).eq.1) ue(i)=uc(i)
              if(qcd(6).lt.1.or.qcd(5).lt.1.or.qcd(4).lt.1)ue(i)=u(i)
           endif
           call flag(qcstauv,qcd,m,6)
           if (v(i).eq.badflag) then
              ve(i)=vc(i)
           else
              if(qcd(2).eq.1) ve(i)=vc(i)
              if(qcd(6).lt.1.or.qcd(5).lt.1.or.qcd(4).lt.1)ve(i)=v(i)
           endif
           call flag(qcstapm,qcd,m,6)
           if (pmsl(i).eq.badflag) then
              pmsle(i)=pmslc(i)
           else
              if(qcd(2).eq.1) pmsle(i)=pmslc(i)
              if(qcstapm(i).ne.111100) pmsle(i)=pmsl(i)
           endif
           call flag(qcstal,qcd,m,6)
           if (alt(i).eq.badflag) then
              alte(i)=altc(i)
           else
              if(qcd(2).eq.1) alte(i)=altc(i)
              if(qcstal(i).ne.111100) alte(i)=alt(i)
           endif
        enddo !i
c
c  special treatment of new obs
c    if new ob is missing, ta contains analyzed value for i>maxstab
c

12345   continue


        do i=maxstab+1,maxstaa
           tb(i)=badflag
           te(i)=ta(i)
           ta(i)=badflag
           tc(i)=badflag
           tdb(i)=badflag
           tde(i)=tda(i)
           tda(i)=badflag
           tdc(i)=badflag
           ub(i)=badflag
           ue(i)=ua(i)
           ua(i)=badflag
           uc(i)=badflag
           vb(i)=badflag
           ve(i)=va(i)
           va(i)=badflag
           vc(i)=badflag
           pmslb(i)=badflag
           pmsle(i)=pmsla(i)
           pmsla(i)=badflag
           pmslc(i)=badflag
           altb(i)=badflag
           alte(i)=alta(i)
           alta(i)=badflag
           altc(i)=badflag
        enddo !i
c
c write out observations
c
        call writec(stn,maxstaa,m,' Stations   ',atime)
        call writev(t,maxstaa,1,m, ' RAW T OBS  ',atime,on,0.)
        call writev(tb,maxstab,1,m, ' OBS GUESS T',atime,on,0.)
        call writev(ta,maxstab,1,m, ' MOD GUESS T',atime,on,0.)
        call writei(qcstat,maxstaa,m,
     &   ' QC Status : 100000 failed comb, 10000 MD , 1000 Ob chk,',
     &   ' 100 bd ck,10 gs err, 1 Ob msng    ', atime)
        call writev(tc,maxstab,1,m,' OPTIMUM T ',atime,on,0.)
        call writev(te,maxstaa,1,m,' QCd T OBS    ',atime,on,0.)
c
c  write out QC status and obs
c
        call writec(stn,maxstaa,m,' Stations   ',atime)
        call writev(td,maxstaa,1,m, ' RAW TD OBS  ',atime,on,0.)
        call writev(tdb,maxstab,1,m, 'OBS TD GUESS',atime,on,0.)
        call writev(tda,maxstab,1,m, 'MOD TD GUESS',atime,on,0.)
        call writei(qcstatd,maxstaa,m,
     &   ' QC Status : 100000 failed comb, 10000 MD , 1000 Ob chk,',
     &   ' 100 bd ck,10 gs err, 1 Ob msng    ',atime)
        call writev(tdc,maxstab,1,m,' OPTIMUM TD ',atime,on,0.)
        call writec(stn,maxstaa,m,' Stations   ',atime)
        call writev(u,maxstaa,1,m, ' RAW U OBS  ',atime,on,0.)
        call writev(v,maxstaa,1,m, ' RAW V OBS  ',atime,on,0.)
        call writev(ub,maxstab,1,m, ' OBS GUESS U',atime,on,0.)
        call writev(vb,maxstab,1,m, 'OBS V GUESS',atime,on,0.)
        call writev(ua,maxstab,1,m, 'MOD U GUESS',atime,on,0.)
        call writev(va,maxstab,1,m, 'MOD V GUESS',atime,on,0.)
        call writei(qcstauv,maxstaa,m,
     &   ' QC Status : 100000 failed comb, 10000 MD , 1000 Ob chk,',
     &   ' 100 bd ck,10 gs err, 1 Ob msng    ',atime)
        call writev(uc,maxstab,1,m,' OPTIMUM U ',atime,on,0.)
        call writev(vc,maxstab,1,m,' OPTIMUM V  ',atime,on,0.)
        call writev(ue,maxstaa,1,m,' QCd U OBS    ',atime,on,0.)
        call writev(ve,maxstaa,1,m,' QCd V OBS    ',atime,on,0.)
c     
c  write out QC status and obs
c
        call writec(stn,maxstaa,m,' Stations   ',atime)
        call writev(pmsl,maxstaa,1,m, ' RAW PMSL OBS',atime,on,1000.)
        call writev(pmslb,maxstab,1,m, ' OBS GUESS PMSL',atime,
     &          on,1000.)
        call writev(pmsla,maxstab,1,m, ' MOD GUESS PMSL',atime,
     &          on,1000.)
        call writei(qcstapm,maxstaa,m,
     &   ' QC Status : 100000 failed comb, 10000 MD , 1000 Ob chk,',
     &   ' 100 bd ck,10 gs err, 1 Ob msng    ',atime)
        call writev(pmslc,maxstab,1,m,' OPTIMUM PMSL ',atime,
     &          on,1000.)
        call writev(pmsle,maxstaa,1,m,' QCd PMSL OBS    ',atime,
     &          on,1000.)
c
c  write out QC status and obs
c
        call writec(stn,maxstaa,m,' Stations   ',atime)
        call writev(alt,maxstaa,1,m, ' RAW ALT OBS  ',atime,on,1000.)
        call writev(altb,maxstab,1,m, 'OBS ALT GUESS',atime,on,1000.)
        call writev(alta,maxstab,1,m, 'MOD ALT GUESS',atime,on,1000.)
        call writei(qcstal,maxstaa,m,
     &   ' QC Status : 100000 failed comb, 10000 MD , 1000 Ob chk,',
     &   ' 10 bd ck,10 gs err, 1 Ob msng    ',atime)
        call writev(altc,maxstab,1,m,' OPTIMUM ALT ',atime,on,1000.)
        call writev(alte,maxstaa,1,m,' QCd ALT OBS    ',atime,
     &          on,1000.)
        write (*,3004) atime,it
 3004   format (1x,a24,3x, 'it= ',i4)
        write(*,3006) maxstaa
 3006   format(1x,i3)
c     
c.....  Write out the NetCDF QC output file.  Convert units to
c.....  the LAPS standard output units first.
c

c************convert units
c
c
c
        do i=1,maxstaa
           stations(i)(1:20) = '                    '
           stations(i)(1:5) = stn(i)(1:5)
        enddo !i
c        
        ext_out = 'lsq'
        call get_directory('lsq', dir_out, len)
c        

        call write_qc_cdf(i4time, dir_out, ext_out, m, maxstaa, 
     1       stn, provider, reptype, lat, lon, elev,
     1       qcstat,  t,    tb,    ta,    tc,    te,    tf,  
     1       qcstatd, td,   tdb,   tda,   tdc,   tde,   tdf,  
     1       qcstauv, u,    ub,    ua,    uc,    ue,    uf,  
     1                v,    vb,    va,    vc,    ve,    vf,  
     1       qcstapm, pmsl, pmslb, pmsla, pmslc, pmsle, pmslf, 
     1       qcstal,  alt,  altb,  alta,  altc, alte,   altf, 
     1       status)
c
c.....  Now convert the units back for other storage files.
c


c*********************


c
c  put these in ta,tda for next cycle
c
        call replace(te,ta,maxstaa,1,m,1)
        call replace(tde,tda,maxstaa,1,m,1)
        call replace(ue,ua,maxstaa,1,m,1)
        call replace(ve,va,maxstaa,1,m,1)
        call replace(pmsle,pmsla,maxstaa,1,m,1)
        call replace(alte,alta,maxstaa,1,m,1)
c
c.....  Write out the updated monster file.
c

        call s_len(monfile, len)
        print*,'Writing file ',monfile(1:len)
        open(16,file=monfile(1:len),
     &         form='unformatted',status='unknown')
        write(16) atime_cur, i4time
        write(16) vit,vitd,vitu,vitv,vitpm,vital
        write(16) wit,witd,witu,witv,witpm,wital
        write(16) pt,ptd,pu,pv,ppm,pal
        write(16) monster
        write(16) ta,tda,ua,va,pmsla,alta
        close(16)
        print*,' complete'
c
c.....  Write out a "clean" LSO file.
c

        go to 9999   !skip until can check on the arrays.


        do i=1,maxstaa
           rtime = float( time(i) )
           store_1(i,1) = lat(i)          ! station latitude
           store_1(i,2) = lon(i)          ! station longitude
           store_1(i,3) = elev(i)         ! station elevation
           store_1(i,4) = rtime           ! observation time
c
           store_2(i,1) = te(i)           ! temperature (deg f)
           store_2(i,2) = tde(i)          ! dew point (deg f)
           store_2(i,3) = badflag         ! Relative Humidity
c
           if(u(i).eq.badflag .or. abs(u(i)).gt.200. .or. 
     &        v(i).eq.badflag .or. abs(v(i)).gt.200.) then
              spd = badflag
              dir = badflag
c
           elseif(u(i).eq.0.0 .and. v(i).eq.0.0) then
              spd = 0.0
              dir = 0.0                      !Undefined
c
           else
              spd = sqrt(u(i)*u(i) + v(i)*v(i) )   !speed
              dir = 57.2957795 * (atan2(u(i),v(i))) + 180.   !dir
           endif

           store_3(i,1) = dir             ! wind dir (deg)
           store_3(i,2) = spd             ! wind speed (kt)
           store_3(i,3) = badflag         ! wind gust dir (deg)
           store_3(i,4) = badflag         ! wind gust speed (kt)
c
           store_4(i,1) = alte(i)         ! altimeter setting (mb)
           store_4(i,2) = pstn(i)         ! station pressure (mb)
           store_4(i,3) = pmsle(i)        ! MSL pressure (mb)
           store_4(i,4) = float(delpch(i))! 3-h press change character
           store_4(i,5) = delp(i)         ! 3-h press change (mb)
c
           store_5(i,1) = vis(i)          ! visibility (miles)
           store_5(i,2) = solar(i)        ! solar radiation 
           store_5(i,3) = sfct(i)         ! soil/water temperature
           store_5(i,4) = sfcm(i)         ! soil moisture
c     
           store_6(i,1) = pcp1(i)         ! 1-h precipitation
           store_6(i,2) = pcp3(i)         ! 3-h precipitation
           store_6(i,3) = pcp6(i)         ! 6-h precipitation
           store_6(i,4) = pcp24(i)        ! 24-h precipitation
           store_6(i,5) = snow(i)         ! snow cover
c
           store_7(i,1) = float(kkk_s(i)) ! number of cloud layers
           store_7(i,2) = max24t(i)       ! 24-h max temperature
           store_7(i,3) = min24t(i)       ! 24-h min temperature
c
           store_2ea(i,1) = t_ea(i)
           store_2ea(i,2) = td_ea(i)
           store_2ea(i,3) = rh_ea(i)
           store_3ea(i,1) = dd_ea(i)
           store_3ea(i,2) = ff_ea(i)
           store_4ea(i,1) = p_ea(i)
           store_4ea(i,2) = alt_ea(i)
           store_5ea(i,1) = vis_ea(i)
           store_5ea(i,2) = solar_ea(i)
           store_5ea(i,3) = sfct_ea(i)
           store_5ea(i,4) = sfcm_ea(i)
           store_6ea(i,1) = pcp_ea(i)
           store_6ea(i,2) = snow_ea(i)
c
        enddo !nn
c     
        call write_surface_obs(atime_cur,outfile,n_obs_g,
     &    n_obs_b,wmoid,stations,provider,wx,reptype,autostntype,
     &    store_1,store_2,store_3,store_4,store_5,store_6,store_7,
     &    store_2ea,store_3ea,store_4ea,store_5ea,store_6ea,
     &    store_cldamt,store_cldht,maxsta,jstatus)
c
c.....  That's it.
c     
 9999   continue
        return
        end
      

