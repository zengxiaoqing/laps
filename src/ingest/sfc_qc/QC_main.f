
       Program QC_main
c
       include 'lapsparms.cmn'
       include 'surface_obs.inc'
       character laps_domain*9
       integer istatus
       logical lrunqc
c
       laps_domain = 'nest7grid'
       call get_laps_config(laps_domain,istatus)
       if(istatus .ne. 1) then
          write(6,*) 'QC_main: Error getting domain dimensions'
          stop
       endif
c
       call get_sfcqc_nl(lrunqc,istatus)
       if (istatus .ne. 1) then
         write(6,*) 'QC_main: Error getting sfcqc_nl file'
         stop
       endif
       if (lrunqc) then
         call QC_main_sub(NX_L_CMN,NY_L_CMN,nk_laps,maxstns_cmn,
     &                  laps_cycle_time_cmn,badflag)
       else
         write(6,*) 'QC_main: lrunqc false...exiting.'
       endif
 
       end
c
c
       subroutine QC_main_sub(ni,nj,nk,m,laps_cycle_time,badflag)
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
c       08 Dec 1999  John McGinley and Peter Stamus, NOAA/FSL
c          New version.  This manifestation is based on an improved 
c          model with ob, mo, and buddy estimates with the Kalman 
c          using obs as orginally done.  Also some housekeeping 
c          changes.
c
c     Notes:
c
c*********************************************************************
c
c.....  Dimensions for monster.dat file
c
        parameter( nvar = 10 )  
        parameter( ncycles = 50 )    !number of cycles to carry
c
c.....  set up tabled exponential and freq used constants
c
        common tab(10000),pi,re,rdpdg,reorpd
c
c  arrays for reading and converting observations
c
c ....lso reader arrays
c
        real lat(m), lon(m), elev(m)
        real t(m), t_ea(m), max24t(m), min24t(m)
        real td(m), td_ea(m), rh(m), rh_ea(m)
        real dd(m), ddg(m), dd_ea(m)
        real ff(m), ffg(m), ff_ea(m)
        real alt(m), pmsl(m), pstn(m)
        real alt_ea(m), delp(m), p_ea(m)
        real vis(m), vis_ea(m)
        real solar(m), solar_ea(m)
        real sfct(m), sfct_ea(m)
        real sfcm(m), sfcm_ea(m)
        real pcp1(m), pcp3(m), pcp6(m), pcp24(m)
        real snow(m), snow_ea(m), pcp_ea(m)
        real store_cldht(m,5),pctmin
c
        integer i4time, wmoid(m), jstatus
        integer time(m), delpch(m), kkk_s(m)
        integer index(m), indexa(m), indexb(m)
        integer minobs,ihr
c
        character store_cldamt(m,5)*4
        character stations(m)*20, provider(m)*11
        character stations_out(m)*20
        character reptype(m)*6, autostntype(m)*6
c
c.....  LSO write arrays
c
        real  store_1(m,4), 
     &        store_2(m,3), store_2ea(m,3),
     &        store_3(m,4), store_3ea(m,2),
     &        store_4(m,5), store_4ea(m,2),
     &        store_5(m,4), store_5ea(m,4),
     &        store_6(m,5), store_6ea(m,2),
     &        store_7(m,3)
        real  store_chtout(m,5)
        integer   wmoid_out(m)
        character outfile*256
        character store_camtout(m,5)*4, wx_out(m)*25
c
c.....  Need two sets of past obs  and a raw array -- a and b, and r
c
        real   ta(m),tda(m),dda(m),ffa(m),ua(m),va(m),fff(m)
        real   tb(m),tdb(m),ddb(m),ffb(m),ub(m),vb(m)
        real   tr(m),tdr(m),ddr(m),ffr(m),ur(m),vw(m)
        real   lata(m),lona(m),eleva(m)
        real   latb(m),lonb(m),elevb(m)
        real   pstna(m),pmsla(m),alta(m)
        real   pstnb(m),pmslb(m),altb(m)
        real   pstnr(m),pmslr(m),altr(m)
c
c.....  Utility arrays
c
        real   u(m), v(m)
        real   tc(m),tdc(m),uc(m),vc(m)
        real   te(m),tde(m),ue(m),ve(m)
        real   tf(m),tdf(m),uf(m),vf(m)
        real   pmslc(m),altc(m)
        real   pmsle(m),alte(m)
        real   pmslf(m),altf(m)
        real   ncmt(m),ncmtd(m),ncmu(m),ncmv(m),ncmpm(m),ncmalt(m)
c
c.....  Arrays for Kalman
c
c  error covariances for each variable
        real   w(m,m),vv(m,m),pt(m,m),ptd(m,m),pu(m,m),pv(m,m)
        real   ppm(m,m),pal(m,m)
c
c  input arrays-model and analysis weighting arrays
        real   dt(m),dta(m),dtd(m),dtda(m)
        real   du(m),dua(m),dv(m),dva(m),dpm(m),dpma(m),dalt(m),dalta(m)
        real   F(m,m),mwt(m,m),dwt(m,m)
        real   wmt(m),wot(m),wbt(m),wmtd(m),wotd(m),wbtd(m)
        real   wmu(m),wou(m),wbu(m),wmv(m),wov(m),wbv(m)
        real   wmp(m),wop(m),wbp(m),wma(m),woa(m),wba(m)
        real   wkt(m),wktd(m),wku(m),wkv(m),wkp(m),wka(m)
        real   tcr(nvar),tmr(nvar),tor(nvar),tcb(nvar),tmb(nvar)
        real   tob(nvar),totr(nvar),totb(nvar)
c
c  output arrays
        real   xt(m),xtt(m)
c
c  taylor series, buddy, and observation estimates 
        real   monster(m,ncycles,nvar),yta(m),ytda(m),yua(m)
        real   yva(m),ypma(m),yalta(m)
        real   byta(m),bytda(m),byua(m),byva(m),bypma(m),byalta(m)
        real   bmonster(m,24,nvar)
c
c  arrays for error processing for each variable by iteration
        real   ar(m),wr(m),vr(m),wit(m,m),witd(m,m),witu(m,m),witv(m,m)
        real   witpm(m,m),wital(m,m),vit(m,m),vitd(m,m),vitu(m,m)
        real   vitv(m,m),vitpm(m,m),vital(m,m)
c
c  theta arrays for times a and b
        real   theta(m),thetaa(m),thetab(m),thetaf(m)
c
c  utility arrays
        real   b(m),c(m,m),zot(m),msar(m)
c
c  error variables
        real   oberr,moderr
c
c arrays for data input
        integer  obstime(m),kloud(m),idp3(m)
c
c  quality control status for each variable
        integer qcstat(m),qcstatd(m),qcstapm(m),qcstal(m),qcstau(m),
     &          qcstav(m),qmonster(m,ncycles,nvar),ht(m),htd(m)
     &          ,hu(m),hv(m),hp(m),ha(m)
c
c switches 
        integer on,off
c
        character providera(m)*11, providerb(m)*11
        character reptypea(m)*6, reptypeb(m)*6
c
        logical exists, flagstart
        data exists/.false./
        data flagstart/.false./
c        
c character arrys for file names, station names, time, wx symbols
c
        real        lat_grid(ni,nj), lon_grid(ni,nj), grid_spacing
        real        fon, nof
        parameter(  fon = 5./9., nof = 9./5. )
	integer     i4prev1, i4prev2, laps_cycle_time, cycle
        integer     bgmodel
	character*9 filename, fname1, fname2
        Character   atime*24, atime_cur*24, stn(m)*5, wx(m)*25
        character   stna(m)*5, stnb(m)*5
        character   monfile*256, atime_mon*24, monjrfile*256
        character   dir_mon*256, dir_out*256, ext_out*31,dir_qcr*256
        character   dir_s*256,ext_s*31,units*10,comment*125,var_s*3
        character   nanvar*10
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
        cycler=float(cycle)/3600.
	i4prev1 = i4time - cycle
	i4prev2 = i4prev1 - cycle
        call make_fnam_lp(i4prev1, fname1, istatus)
        call make_fnam_lp(i4prev2, fname2, istatus)
c  set minimum number of obs allowed and minimum pct of previous list
        minobs=10
        pctmin=.3 
c number of standard deviations meeting bad criteria in grosserr
        badthr=5.5
c
c  observation and model error assumptions
c *** expected obs accuracy for each variable at each station
c      is calculated in the LSO file, and is passed in by
c      the read_surface_data routine.  Can use these someday.***
c
c
c observation and model error assumptions
        oberr=1.0*9./5.!C error converted to F
        oberrt=0.5*9./5.
        oberrtd=1.0*9./5.
        oberruv=0.515/.515!m/sec error conv to knts
        oberrpm=0.50
        oberral=0.50
        moderr=3.0
        moderrt=3.0
        moderrtd=3.0
        moderru=3.0
        moderrpm=1.5
        moderral=1.5
        grosst=30.
        grosstd=35.
        grossuv=40.
        grosspm=5.
        grossal=5.
c  offset for model operator (>> than ob-fg)
c
        offset=10000.
        on=1
        off=0
c
c..... Set constants and zero some arrays
c
        re=6371122.
        pi = 4.0 * atan(1.0)
        rdpdg=pi/180.
        reorpd=re*rdpdg  
c
        do i=1,m
           index(i)  = 0
           indexa(i) = 0
           indexb(i) = 0
        enddo !i
c
c  averaging period for error for w and vv
c
        ia=24
c
c averaging period for model weights in cycles; days for bias error
c
        length=24
c
c  set up exponential lookup table
c
        do i=1,10000
           tab(i)= exp (-float((i-1)*(i-1))/1000.)
        enddo !i
c
c.....  Get the grid lat/lons from the static file.
c
        call get_directory('static', dir_s, len)
        ext_s = 'nest7grid'
        var_s = 'LAT'
        call rd_laps_static(dir_s,ext_s,ni,nj,1,var_s,units,comment,
     &                      lat_grid ,grid_spacing,istatus)
        var_s = 'LON'
        call rd_laps_static(dir_s,ext_s,ni,nj,1,var_s,units,comment,
     &                      lon_grid ,grid_spacing,istatus)
c
c.....  Check for a 'monster' file in the 'LSQ' output
c.....  directory.  If there's one there, read it.
c.....  If not, fill the error arrays here.
c
        call get_directory('lsq', dir_mon, len)
        monfile = dir_mon(1:len) // 'monster.dat'
        monjrfile = dir_mon(1:len) // 'monster_jr.dat'
        call s_len(monfile, len)
        INQUIRE(FILE=monfile(1:len), EXIST=exists)
        if(.not. exists) then
           flagstart = .true.  !no monster file
        else 
           flagstart = .false. !yes monster file
        endif
c
 500    continue

        IF(flagstart) THEN
c
c.....  No monster file, initilize error arrays
c
           call zero(w,m,m)
           call zero(pt,m,m)
           call zero(ptd,m,m)
           call zero(pu,m,m)
           call zero(pv,m,m)
           call zero(ppm,m,m)
           call zero(pal,m,m)
           call zero(vv,m,m)

           do i=1,m
              zot(i)=0. 
              msar(i)=badflag
c
c initial kalmod errors rms
              wot(i)=3.
              wmt(i)=3.
              wbt(i)=3.
              wotd(i)=3.
              wmtd(i)=3.
              wbtd(i)=3.
              wou(i)=3.
              wmu(i)=3.
              wbu(i)=3.
              wov(i)=3.
              wmv(i)=3.
              wbv(i)=3.
              wop(i)=1.
              wmp(i)=1.
              wbp(i)=1.
              woa(i)=1.0
              wma(i)=1.0
              wba(i)=1.0 
c         
              do k=1,nvar
              do l=1,24
               bmonster(i,l,k)=0.
c work in first guess temp bias based on past performance
               bmonster(i,l,1)=3.5*cos(6.28318*(float(l)-2.5)/24.)
     &                    +     0.5*cos(6.28318*(float(l)-8.5)/24.)
              enddo!l
              do j=1,ncycles
                 monster(i,j,k) = badflag
              enddo !j
              enddo !k
c
                 iiiii=-i4time+i
              do j=1,2
                 if(iiiii.ge.0)iiiii=-1-iiiii
c
c     random number intial seed
                 vit(i,j)=oberrt*(-1)**int(10.*ffz(iiiii,1))
                 wit(i,j)=moderrt*(-1)**int(10.*ffz(iiiii,1))
                 witd(i,j)=moderrtd*(-1)**int(10.*ffz(iiiii,1))
                 vitd(i,j)=oberrtd*(-1)**int(10.*ffz(iiiii,1))
                 witu(i,j)=moderru*(-1)**int(10.*ffz(iiiii,1))
                 witv(i,j)=moderru*(-1)**int(10.*ffz(iiiii,1))
                 vitv(i,j)=oberruv*(-1)**int(10.*ffz(iiiii,1))
                 vitu(i,j)=oberruv*(-1)**int(10.*ffz(iiiii,1))
                 witpm(i,j)=moderrpm*(-1)**int(10.*ffz(iiiii,1))
                 vitpm(i,j)=oberrpm*(-1)**int(10.*ffz(iiiii,1))
                 wital(i,j)=moderral*(-1)**int(10.*ffz(iiiii,1))
                 vital(i,j)=oberral*(-1)**int(10.*ffz(iiiii,1))
              enddo !on j
c
              Pt(i,i)=moderrt*moderrt
              Ptd(i,i)=moderrtd*moderrtd
              Pu(i,i)=moderru*moderru
              Pv(i,i)=moderru*moderru
              Ppm(i,i)=moderrpm*moderrpm
              Pal(i,i)=moderral*moderral
           enddo !i
c
           ite=0
c
c.....  "it" is the index that records how many cycles have been 
c        written to monster
c.....  start it off with 0 for start up
c
           it=0
c
c.....  Get the LSO data file for 2 cycles ago.
c
           call read_surface_data(i4prev2,atime,n_obs_g,n_obs_b,time,
     &     wmoid,stations,providerb,wx,reptypeb,autostntype,latb,lonb,
     &     elevb,tb,tdb,rh,ddb,ffb,ddg,ffg,altb,pstnb,pmslb,delpch,delp,
     &     vis,solar,sfct,sfcm,pcp1,pcp3,pcp6,pcp24,snow,kkk_s,max24t,
     &     min24t,t_ea,td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,
     &     solar_ea,sfct_ea,sfcm_ea,pcp_ea,snow_ea,store_cldamt,
     &     store_cldht,m,jstatus)
c
           iflag=0
           if(jstatus .ne. 1) then
              print *,
     &             ' No LSO data for (Hr-2). '
              iflag=1
           else
              print *,' Found LSO data (Hr-2) at ', atime
           endif
c
           if (iflag.eq.1) go to  234
           maxstab=n_obs_b
           do k=1,maxstab
              stnb(k)(1:5)=stations(k)(1:5)
              indexb(k) = k
           enddo !k
c
           print*, maxstab, 'obs read for cycle '    
c
c
c.....  Now get LSO data for 1 cycle ago.
c
  234      continue
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
             if(iflag.eq.1) then ! bset for hr-2 was missing
              write(6,*) 'Data for Hr-2 was missing so using Hr-1 only'
              maxstab=maxstaa
              call replace(ta,tb,m,1,m,1)
              call replace(tda,tdb,m,1,m,1)
              call replace(ffa,ffb,m,1,m,1)
              call replace(dda,ddb,m,1,m,1)
              call replace(pmsla,pmslb,m,1,m,1)
              call replace(alta,altb,m,1,m,1)
              call replace(lata,latb,m,1,m,1)
              call replace(lona,lonb,m,1,m,1)
              call replace(eleva,elevb,m,1,m,1)
             endif
           pct=float(maxstaa)/float(maxstab)
           if(pct.lt.pctmin.or.maxstaa.lt.minobs) then
            print*, 'Frac of obs A/B is ',pct,' min frac is ',pctmin
            print*, 'Insufficient count ',maxstaa,' for min of ',minobs
            print*,' No Kalman QC is possible...try next cycle'
c           stop
           endif
           do k=1,maxstaa
              stna(k)(1:5)=stations(k)(1:5)
              indexa(k) = k
           enddo !k
           if(iflag.eq.1) then ! b set was missing so put stna -> stnb
            do k=1,maxstaa
              stnb(k)(1:5)=stna(k)(1:5)
              indexb(k) = indexa(k)
            enddo !k
           endif
           print*, maxstaa, 'obs read for cycle '
c trim data lists for both a and b data for those stations with all 
c missing data
           call trimlist(tb,tdb,ddb,ffb,latb,lonb,elevb,pstnb,pmslb,
     &       altb,stnb,providerb,reptypeb,indexb,maxstab,m,badflag)
c do gross error checks
           call grosserr(tb,tb,stnb,maxstab,m,badflag,badthr,grosst
     &                  ,qcstat,off)
           call grosserr(tdb,tdb,stnb,maxstab,m,badflag,badthr,grosstd
     &                  ,qcstatd,off)
           call grosserr(ffb,ffb,stnb,maxstab,m,badflag,badthr,grossuv
     &                  ,qcstau,off)
           call replacei(qcstau,qcstav,maxstab,m)
           call grosserr(pmslb,pmslb,stnb,maxstab,m,badflag,badthr,
     &                 grosspm ,qcstapm,off)
           call grosserr(altb,altb,stnb,maxstab,m,badflag,badthr,
     &                 grossal,qcstal,off)
           call trimlist(ta,tda,dda,ffa,lata,lona,eleva,pstna,pmsla,
     &       alta,stna,providera,reptypea,indexa,maxstaa,m,badflag)
c do gross error checks
           call grosserr(ta,ta,stna,maxstaa,m,badflag,badthr,grosst
     &                  ,qcstat,off)
           call grosserr(tda,tda,stna,maxstaa,m,badflag,badthr,
     &                 grosstd ,qcstatd,off)
           call grosserr(ffa,ffa,stna,maxstaa,m,badflag,badthr,grossuv
     &                  ,qcstau,off)
           call replacei(qcstau,qcstav,maxstaa,m)
           call grosserr(pmsla,pmsla,stna,maxstaa,m,badflag,badthr,
     &                 grosspm ,qcstapm,off)
           call grosserr(alta,alta,stna,maxstaa,m,badflag,badthr,
     &                  grossal,qcstal,off)
c reorder a and b
           call reorder(ta,tda,dda,ffa,lata,lona,eleva,pstna,pmsla,alta,
     &       stna,providera,reptypea,indexa,
     &       tb,tdb,ddb,ffb,latb,lonb,elevb,pstnb,pmslb,altb,stnb,
     &       providerb,reptypeb,indexb,
     &       maxstaa,maxstab,m,badflag)
           call convuv(ddb,ffb,ub,vb,maxstab,m,badflag)
           call convuv(dda,ffa,ua,va,maxstaa,m,badflag)
c perform grosserror checks
c
c.....  Convert to theta from raw temps
c
           call thet(tb,thetab,maxstaa,eleva,m)
c
c.....  Replace missing values of theta from neighbors in b set
c
           print*,'First fillqc on b array: potl temp'
           call fillqc(stnb,thetab,latb,lonb,elevb,thetab,sl,
     &       maxstaa,m,8./1000.,-1./1000.,200000.,30.,grosst )         
           call thet2T(tb,thetab,maxstaa,eleva,m)
c
c.....  Fill the b arrays so that no data is missing
c
           print*,'dewpoint'
           call fillqc(stnb,tdb,latb,lonb,elevb,thetab,vl,
     &          maxstab,m,15./1000.,-35./1000.,200000.,30.,grosstd)        
           print*,'u       '
           call fillqc(stnb,ub,latb,lonb,elevb,thetab,vl,
     &          maxstab,m,10./1000.,-10./1000.,200000.,30.,grossuv)        
           print*,'v       '
           call fillqc(stnb,vb,latb,lonb,elevb,thetab,vl,
     &       maxstab,m,10./1000.,-10./1000.,200000.,30.,grossuv)         
           print*,'mslp    '
           call fillqc(stnb,pmslb,latb,lonb,elevb,thetab,vl,
     &       maxstab,m,2./1000.,-2./1000.,200000.,1000.,grosspm)        
           print*,'alstg   '
           call fillqc(stnb,altb,latb,lonb,elevb,thetab,vl,
     &       maxstab,m,2./1000.,-2./1000.,200000.,1000.,grossal)         
           call thet(ta,thetaa,maxstaa,eleva,m)
           print*,'First fill on a array: potl temp'
           call fillqc(stna,thetaa,lata,lona,eleva,thetaa,sl,
     &             maxstaa,m,8./1000.,-1./1000.,200000.,30.,grosst)       
c
c.....  Now for first guess ta, use thetas and interp to msng stn locations
c
           call thet2T(ta,thetaa,maxstaa,eleva,m)
c
c.....  Same thing for dewpoint and winds in a set
c
           print*,'dewpoint'
           call fillqc(stna,tda,lata,lona,eleva,thetaa,vl,
     &          maxstaa,m,15./1000.,-35./1000.,200000.,30.,grosstd)               
           print*,'u       '
           call fillqc(stna,ua,lata,lona,eleva,thetaa,vl,
     &          maxstaa,m,10./1000.,-10./1000.,200000.,30.,grossuv)        
           print*,'v       '
           call fillqc(stna,va,lata,lona,eleva,thetaa,vl,
     &       maxstaa,m,10./1000.,-10./1000.,200000.,30.,grossuv)     
           print*,'pmsl    '
           call fillqc(stna,pmsla,lata,lona,eleva,thetaa,vl,
     &       maxstaa,m,2./1000.,-2./1000.,200000.,1000.,grosspm)    
           print*,'alstg   '
           call fillqc(stna,alta,lata,lona,eleva,thetaa,vl,
     &       maxstaa,m,2./1000.,-2./1000.,200000.,1000.,grossal)
           icnt=1
c if the b set has much fewer obs than the a set just write a into b
           if(pct.gt.(1./pctmin)) then
            print*,'insufficient obs in the b set: replacing with a'
            call replace(ta,tb,maxstaa,1,m,1)  
            call replace(tda,tdb,maxstaa,1,m,1)  
            call replace(ua,ub,maxstaa,1,m,1)  
            call replace(va,vb,maxstaa,1,m,1)  
            call replace(pmsla,pmslb,maxstaa,1,m,1)  
            call replace(alta,altb,maxstaa,1,m,1)  
           endif
           it=1
           call writemon(tb,tdb,ub,vb,pmslb,altb,nvar,maxstaa,m,
     &          ncycles,monster,it) 
           call writeqc(qcstat,qcstatd,qcstau,qcstav,qcstapm,
     &              qcstal,nvar,maxstaa,m,
     &          ncycles,qmonster,it) 
c
        ELSE
c
           call s_len(monfile, len)
           open(15,file=monfile(1:len),
     &         form='unformatted',status='old')
           read(15) atime_mon, i4time_mon
           atime=atime_mon
c
c.....  Check the time in the monster file.  If its too old,
c.....  bag it and start over.
c
           print *,' Found monster file valid for: ', atime_mon
           if( abs(i4time-i4time_mon) .gt. laps_cycle_time) then
              print *,
     &      '    But this monster file too old. Starting over...'
              flagstart = .true.
              close(15)
              go to 500
           endif
c
           read(15) vit,vitd,vitu,vitv,vitpm,vital
           read(15) wit,witd,witu,witv,witpm,wital
           read(15) pt,ptd,pu,pv,ppm,pal
           read(15) monster,bmonster
           read(15) stna,lata,lona,eleva
           read(15) providera, reptypea
           read(15) ta,tda,ua,va,pmsla,alta
           read(15) ncmt,ncmtd,ncmu,ncmv,ncmpm,ncmalt
           read(15) qcstat,qcstatd,qcstau,qcstav,qcstapm,qcstal
           read(15) qmonster
           read(15) it,maxstaa
           close(15)
c
c.....  Read the other storage file.
c
           call s_len(monjrfile, len)
           open(18,file=monjrfile(1:len),
     &         form='unformatted',status='old')
c
           read(18) wot,wotd,wou,wov,wop,woa
           read(18) wbt,wbtd,wbu,wbv,wbp,wba
           read(18) wmt,wmtd,wmu,wmv,wmp,wma
           read(18) tor,tcr,tmr,tob,tcb,tmb
           read(18) totr,totb
           close(18)

              
        ENDIF
c read time of either monster or lso and convert hour to a integer
c need to add one since this atime is from an hour ago
        read(atime(12:14),1098) ihr
        ihr=ihr+1
        if(ihr.eq.24)ihr=0
        print*,'GMT hour is ',ihr 
 1098   format(i3)
        if(ihr.eq.0) ihr=24
c
c.....  This is the start of routine
c
        maxstab=maxstaa
c
c.....  Put the obs in the monster stack
c
        it=it+1
           call writemon(ta,tda,ua,va,pmsla,alta,nvar,maxstaa,m,
     &          ncycles,monster,it) 
        call writeqc(qcstat,qcstatd,qcstau,qcstav,qcstapm,
     &              qcstal,nvar,maxstaa,m,
     &          ncycles,qmonster,it) 
c 
c.....  Compute distance/wind/theta defined weighting functions for 
c.....  error cov
c     
        call weights(mwt,dwt,ua,va,thetaa,lata,lona,cycle,maxstaa,m,
     &               300000.,30.)
        write(6,*) 'Process temp********************************'
        call project(2,yta,monster,1,nvar,
     &             maxstaa,m,ncycles,nn,atime,it,icyc,oberrt,badthr,
     &             bmonster,ihr)
c
c.....  Input NWP model change
c.....  Here, the model routine gets the current and previous values
c.....  of each variable at each station location in 3-d.  Then it
c.....  calculates a trend (dta, etc) for each variable, which are
c.....  used below.
c
c.....  First, set the background model flag.  A '0' forces the 
c.....  'compute_sfc_bgfields' routine (in 'model') to use the 3-d
c.....  specific humidity to calculate a station elevation sp hum,
c.....  which is used to compute a station pressure and dew point.
c
        bgmodel = 0
c
c.....  Now call the model routine.
c
        call model(bgmodel,lat_grid,lon_grid,lata,lona,eleva,
     &             i4time,laps_cycle_time,maxstaa,ni,nj,nk,
     &             m,dta,dtda,dua,dva,dpma,dalta,badflag)
c
c..... Use the trend to process this (and the other) variables.
c
        call perturb(yta,ta,yta,maxstaa,m,0.,on)
c
c make the linear model for this iteration
        call kalmod(F,yta,byta,dta,ta,wmt,wot,wbt,offset,
     1    maxstaa,mwt,m)
        call perturb(ta,zot,ta,maxstaa,m,offset,on)      
c
c make tf all missing flag values
        call replace(msar,tf,m,1,m,1)
        call mvmult(F,ta,tf,maxstaa,maxstaa,1,m)
        call perturb(ta,zot,ta,maxstaa,m,offset,off)      
        call replace(ta,tc,maxstaa,1,m,1)  
        call perturb(tf,zot,tf,maxstaa,m,offset,off)      
c validate the kalman estimate and usekalman buddy if it is not consistent
c with last cycles values
        write(6,*) 'Process dewpoint***************************'
        call project(2,ytda,monster,2,nvar,
     &               maxstaa,m,ncycles,nn,atime,it,icyc,oberrtd,badthr,
     &             bmonster,ihr)
        call perturb(ytda,tda,ytda,maxstaa,m,0.,on)
c
c make the linear model for this iteration
        call kalmod(F,ytda,bytda,dtda,tda,wmtd,wotd,wbtd,offset,
     1       maxstaa,mwt,m)
        call perturb(tda,zot,tda,maxstaa,m,offset,on)      
c
c make tdf all missing flag values
        call replace(msar,tdf,m,1,m,1)
        call mvmult(F,tda,tdf,maxstaa,maxstaa,1,m)
        call perturb(tda,zot,tda,maxstaa,m,offset,off)      
        call replace(tda,tdc,maxstaa,1,m,1)  
        call perturb(tdf,zot,tdf,maxstaa,m,offset,off)      
c validate the kalman estimate and usekalman buddy if it is not consistent
c with last cycles values
        write(6,*) 'Process u wind*******************************'
        call project(2,yua,monster,3,nvar,
     &                maxstaa,m,ncycles,nn,atime,it,icyc,oberruv,badthr,
     &             bmonster,ihr)
        call perturb(yua,ua,yua,maxstaa,m,0.,on)
c
c make the linear model for this iteration
        call kalmod(F,yua,byua,dua,ua,wmu,wou,wbu,offset,
     1       maxstaa,mwt,m)
        call perturb(ua,zot,ua,maxstaa,m,offset,on)      
c
c make uf all missing flag values
        call replace(msar,uf,m,1,m,1)
        call mvmult(F,ua,uf,maxstaa,maxstaa,1,m)
        call perturb(ua,zot,ua,maxstaa,m,offset,off)      
        call replace(ua,uc,maxstaa,1,m,1)  
        call perturb(uf,zot,uf,maxstaa,m,offset,off)      
c validate the kalman estimate and usekalman buddy if it is not consistent
c with last cycles values
        write(6,*) 'Process v wind*******************************'
        call project(2,yva,monster,4,nvar,
     &                maxstaa,m,ncycles,nn,atime,it,icyc,oberruv,badthr,
     &             bmonster,ihr)
        call perturb(yva,va,yva,maxstaa,m,0.,on)
c
c make the linear model for this iteration
        call kalmod(F,yva,byva,dva,va,wmv,wov,wbv,offset,
     1       maxstaa,mwt,m)
        call perturb(va,zot,va,maxstaa,m,offset,on)      
c
c make vf all missing flag values
        call replace(msar,vf,m,1,m,1)
        call mvmult(F,va,vf,maxstaa,maxstaa,1,m)
        call perturb(va,zot,va,maxstaa,m,offset,off)      
        call replace(va,vc,maxstaa,1,m,1)  
        call perturb(vf,zot,vf,maxstaa,m,offset,off)      
c validate the kalman estimate and usekalman buddy if it is not consistent
c with last cycles values
        write(6,*) 'Process msl P****************************** '
        call project(2,ypma,monster,5,nvar,
     &                maxstaa,m,ncycles,nn,atime,it,icyc,oberrpm,badthr,
     &             bmonster,ihr)
        call perturb(ypma,pmsla,ypma,maxstaa,m,0.,on)
c
c make the linear model for this iteration
        call kalmod(F,ypma,bypma,dpma,pmsla,wmp,wop,wbp,offset,
     1       maxstaa,mwt,m)
        call perturb(pmsla,zot,pmsla,maxstaa,m,offset,on)      
c
c make pmslf all missing flag values
        call replace(msar,pmslf,m,1,m,1)
        call mvmult(F,pmsla,pmslf,maxstaa,maxstaa,1,m)
        call perturb(pmsla,zot,pmsla,maxstaa,m,offset,off)      
        call replace(pmsla,pmslc,maxstaa,1,m,1)  
        call perturb(pmslf,zot,pmslf,maxstaa,m,offset,off)      
c validate the kalman estimate and usekalman buddy if it is not consistent
c with last cycles values
        write(6,*) 'Process Alstg*********************************'
        call project(2,yalta,monster,6,nvar,
     &       maxstaa,m,ncycles,nn,atime,it,icyc,oberral,badthr,
     &             bmonster,ihr)
        call perturb(yalta,alta,yalta,maxstaa,m,0.,on)
c
c make the linear model for this iteration
        call kalmod(F,yalta,byalta,dalta,alta,wma,woa,wba,offset,
     1       maxstaa,mwt,m)
        call perturb(alta,zot,alta,maxstaa,m,offset,on)      
c
c make altf all missing flag values
        call replace(msar,altf,m,1,m,1)
        call mvmult(F,alta,altf,maxstaa,maxstaa,1,m)
        call perturb(alta,zot,alta,maxstaa,m,offset,off)      
        call replace(alta,altc,maxstaa,1,m,1)  
        call perturb(altf,zot,altf,maxstaa,m,offset,off)      
c validate the kalman estimate and usekalman buddy if it is not consistent
c with last cycles values
c
c.....  Read in new (current) obs.....maxstaa will go up
c      
c
        call read_surface_data(i4time,atime_cur,n_obs_g,n_obs_b,time,
     &     wmoid,stations,provider,wx,reptype,autostntype,lat,lon,elev,
     &     t,td,rh,dd,ff,ddg,ffg,alt,pstn,pmsl,delpch,delp,vis,solar,
     &     sfct,sfcm,pcp1,pcp3,pcp6,pcp24,snow,kkk_s,max24t,min24t,t_ea,
     &     td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,sfct_ea,
     &     sfcm_ea,pcp_ea,snow_ea,store_cldamt,store_cldht,
     &     m,jstatus)
        iflag=0
        if(jstatus .ne. 1) then
           print *,' No current LSO data for ', filename,'.Setting   ',
     &             'data to Kalman estimates.'
     &      
           call replace(tf,t,m,1,m,1)
           call replace(tdf,td,m,1,m,1)
           call replace(uf,u,m,1,m,1)
           call replace(vf,v,m,1,m,1)
           call replace(pmslf,pmsl,m,1,m,1)
           call replace(altf,alt,m,1,m,1)
           call replace(msar,dd,m,1,m,1)
           call replace(msar,ff ,m,1,m,1)
           call replace(lata,lat,m,1,m,1)
           call replace(lona,lon,m,1,m,1)
           call replace(eleva,elev,m,1,m,1)
           do k=1,maxstaa
            stations(k)(1:5)=stna(k)(1:5)
            index(k) = k
           enddo !k
           maxsta=maxstaa
           iflag=1
        else
           print *,' Found LSO data (current) at ', atime_cur
        endif
        atime = atime_cur
c
        maxsta=n_obs_b
        pct=float(maxsta)/float(maxstaa)
        print*, 'Ob count for new cycle is ',pct,' of total list'
        do k=1,maxsta
           stn(k)(1:5)=stations(k)(1:5)
           index(k) = k
        enddo !k
        if(iflag.eq.0)print*, maxsta,'obs read for cycle ',it+1
        write(*,1066) atime
 1066   format(1x,a24)
c eliminate observations with all missing data 
        if(iflag.eq.0) then
        call trimlist(t,td,dd,ff,lat,lon,elev,pstn,pmsl,alt,
     &       stn,provider,reptype,index,maxsta,m,badflag)
c perform gross error checks
        call reorder(t,td,dd,ff,lat,lon,elev,pstn,pmsl,alt,
     &          stn,provider,reptype,index,
     &          ta,tda,dda,ffa,lata,lona,eleva,pstna,pmsla,alta,
     &          stna,providera,reptypea,indexa,
     &          maxsta,maxstaa,m,badflag)
        endif

c put data into "r"arrays to save raw data
        if(iflag.eq.0) then
           call replace(t,tr,maxsta,1,m,1)
           call replace(td,tdr,maxsta,1,m,1)
           call replace(dd,ddr,maxsta,1,m,1)
           call replace(ff,ffr,maxsta,1,m,1)
           call replace(pmsl,pmslr,maxsta,1,m,1)
           call replace(alt,altr,maxsta,1,m,1)
           call replace(t,te,maxsta,1,m,1)
           call replace(td,tde,maxsta,1,m,1)
           call replace(ff,ue,maxsta,1,m,1)
           call replace(pmsl,pmsle,maxsta,1,m,1)
           call replace(alt,alte,maxsta,1,m,1)
c these grosserror checks are against the kalman estimate
c ensure kalman values are missing for new stations
           do i=maxstab+1,maxstaa
            tf(i)=badflag
            tdf(i)=badflag
            uf(i)=badflag
            vf(i)=badflag
            pmslf(i)=badflag
            altf(i)=badflag
           enddo 
           Write(6,*) 'New T gross error check'
c   te is now the buddy value at the stn without using stn value
           call grosserr(t,tf,stn,maxsta,m,badflag,badthr,grosst
     &           ,qcstat,on)
           call fillone(stn,te,lat,lon,elev,theta,sl,maxsta,m
     &                  ,8./1000.,-1./1000.,200000.,30.)
           call recover(stn,wx,index,t,tr,te,ncmt,qcstat,
     &          tr,tf,tdr,tdf,ddr,ffr,
     &          uf,vf,pmslr,pmslf,badflag,grosst/2.,maxsta,m)

           Write(6,*) 'New TD gross error check'
           call grosserr(td,tdf,stn,maxsta,m,badflag,badthr,grosstd
     &           ,qcstatd,on)
           call fillone(stn,tde,lat,lon,elev,theta,sl,maxsta,m
     &                  ,8./1000.,-1./1000.,200000.,30.)
           call recover(stn,wx,index,td,tdr,tde,ncmtd,qcstatd,
     &          tr,tf,tdr,tdf,ddr,ffr,
     &          uf,vf,pmslr,pmslf,badflag,grosstd/2.,maxsta,m)
           Write(6,*) 'New WIND SPEED gross error check'
           do i=1,maxsta
            if(uf(i).ne.badflag) then
             fff(i)=sqrt(uf(i)*uf(i)+vf(i)*vf(i))
             if(fff(i).ge.99.) fff(i)=badflag
            else
             fff(i)=badflag
            endif
           enddo
           call grosserr(ff,fff,stn,maxsta,m,badflag,badthr,grossuv
     &           ,qcstau,on)
           call fillone(stn,ue,lat,lon,elev,theta,sl,maxsta,m
     &                  ,10./1000.,-10/1000.,200000.,30.)
           call recover(stn,wx,index,ff,ffr,ue,ncmu,qcstau,
     &          tr,tf,tdr,tdf,ddr,ffr,
     &          uf,vf,pmslr,pmslf,badflag,grossuv/2.,maxsta,m)
           call replacei(qcstau,qcstav,maxsta,m)
           Write(6,*) 'New PMSL gross error check'
           call grosserr(pmsl,pmslf,stn,maxsta,m,badflag,badthr,
     &           grosspm,qcstapm,on)
           call fillone(stn,pmsle,lat,lon,elev,theta,sl,maxsta,m
     &                  ,2./1000.,-2/1000.,200000.,1000.)
           call recover(stn,wx,index,pmsl,pmslr,pmsle,ncmpm,qcstapm,
     &          tr,tf,tdr,tdf,ddr,ffr,
     &          uf,vf,pmslr,pmslf,badflag,grosspm/2.,maxsta,m)
           Write(6,*) 'New ALT gross error check'
           call grosserr(alt,altf,stn,maxsta,m,badflag,badthr,
     &           grossal,qcstal,on)
           call fillone(stn,alte,lat,lon,elev,theta,sl,maxsta,m
     &                  ,2./1000.,-2/1000.,200000.,1000.)
           call recover(stn,wx,index,alt,altr,alte,ncmalt,qcstal,
     &          tr,tf,tdr,tdf,ddr,ffr,
     &          uf,vf,pmslr,pmslf,badflag,grossal/2.,maxsta,m)
c reorder new obs
c
c.....  Convert new winds to u,v
c
        call convuv(dd,ff,u,v,maxstaa,m,badflag)
        call convuv(ddr,ffr,ur,vr,maxstaa,m,badflag)
       else !all obs are missing
           call replace(msar,tr,maxsta,1,m,1)
           call replace(msar,tdr,maxsta,1,m,1)
           call replace(msar,ur,maxsta,1,m,1)
           call replace(msar,vr,maxsta,1,m,1)
           call replace(msar,pmslr,maxsta,1,m,1)
           call replace(msar,altr,maxsta,1,m,1)
       endif
c.....  Move last cycle to tc; new obs to ta; ready the winds 
c
        
        call replace(ta,tc,maxstaa,1,m,1)
        call replace(tf,tb,maxstaa,1,m,1)
        call replace(t,te,maxstaa,1,m,1)
        call replace(t,ta,maxstaa,1,m,1)
        call replace(tda,tdc,maxstaa,1,m,1)
        call replace(tdf,tdb,maxstaa,1,m,1)
        call replace(td,tde,maxstaa,1,m,1)
        call replace(td,tda,maxstaa,1,m,1)
        call replace(dd,dda,maxstaa,1,m,1)
        call replace(ff,ffa,maxstaa,1,m,1)
        call replace(ua,uc,maxstaa,1,m,1)
        call replace(uf,ub,maxstaa,1,m,1)
        call replace(u,ue,maxstaa,1,m,1)
        call replace(u,ua,maxstaa,1,m,1)
        call replace(va,vc,maxstaa,1,m,1)
        call replace(vf,vb,maxstaa,1,m,1)
        call replace(v,ve,maxstaa,1,m,1)
        call replace(v,va,maxstaa,1,m,1)
        call replace(pmsla,pmslc,maxstaa,1,m,1)
        call replace(pmslf,pmslb,maxstaa,1,m,1)
        call replace(pmsl,pmsle,maxstaa,1,m,1)
        call replace(pmsl,pmsla,maxstaa,1,m,1)
        call replace(alta,altc,maxstaa,1,m,1)
        call replace(altf,altb,maxstaa,1,m,1)
        call replace(alt,alte,maxstaa,1,m,1)
        call replace(alt,alta,maxstaa,1,m,1)
c
c.....  Replace missing data with kalman estimate for this cycle
c
        write(6,*) 
     &    'Fill raw msng obs w/ buddy est in e array: Potl Temp '
        call thet(te,thetaa,maxstaa,eleva,m)
        call fill(stna,thetaa,lata,lona,eleva,thetaa,sl,
     &             maxstaa,m,8./1000.,-1./1000.,200000.,30.)          
        call thet2T(te,thetaa,maxstaa,eleva,m)
        call thet(tb,thetab,maxstaa,eleva,m)
c fill kalman values (with bud est) the new obs where we don't have a kalman
        call fill(stna,thetab,lata,lona,eleva,thetaa,sl,
     &             maxstaa,m,8./1000.,-1./1000.,200000.,30.)          
        call thet2T(tb,thetab,maxstaa,eleva,m)
        write(6,*) 'Dewpoints '
        call fill(stna,tde,lata,lona,eleva,thetaa,sl,
     &             maxstaa,m,8./1000.,-1./1000.,200000.,30.)          
        call fill(stna,tdb,lata,lona,eleva,thetaa,sl,
     &             maxstaa,m,8./1000.,-1./1000.,200000.,30.)          
        write(6,*) 'U-winds   '
        call fill(stna,ue,lata,lona,eleva,thetaa,sl,
     &       maxstaa,m,10./1000.,-10./1000.,200000.,30.)          
        call fill(stna,ub,lata,lona,eleva,thetaa,sl,
     &             maxstaa,m,10./1000.,-10./1000.,200000.,30.)          
        write(6,*) 'V-winds   '
        call fill(stna,ve,lata,lona,eleva,thetaa,sl,
     &             maxstaa,m,10./1000.,-10./1000.,200000.,30.)          
        call fill(stna,vb,lata,lona,eleva,thetaa,sl,
     &             maxstaa,m,10./1000.,-10./1000.,200000.,30.)          
        write(6,*) 'MSL Press '
        call fill(stna,pmsle,lata,lona,eleva,thetaa,sl,
     &             maxstaa,m,2./1000.,-2./1000.,200000.,1000.)          
        call fill(stna,pmslb,lata,lona,eleva,thetaa,sl,
     &             maxstaa,m,2./1000.,-2./1000.,200000.,1000.)          
        write(6,*) 'ALSTG     '
        call fill(stna,alte,lata,lona,eleva,thetaa,sl,
     &             maxstaa,m,2./1000.,-2./1000.,200000.,1000.)          
        call fill(stna,altb,lata,lona,eleva,thetaa,sl,
     &             maxstaa,m,2./1000.,-2./1000.,200000.,1000.)
c
c.....  There may be times when all variables are missing so fill won't work
c.....  thus put the fully filled kalman estimate "b" into the "e" variable 
c
        call replacemsng(tb,te,maxstaa,1,m,1,badflag)
        call replacemsng(tdb,tde,maxstaa,1,m,1,badflag)
        call replacemsng(ub,ue,maxstaa,1,m,1,badflag)
        call replacemsng(vb,ve,maxstaa,1,m,1,badflag)
        call replacemsng(pmslb,pmsle,maxstaa,1,m,1,badflag)
        call replacemsng(altb,alte,maxstaa,1,m,1,badflag)
c
c.....  For all new stations appearing use the above buddy estimates for 
c.....  the kalman xt and previous cycle value tc; assume all ob, buddy, and 
c.....  nwp trends are zero
c
        do i=maxstab+1,maxstaa
           tf(i)=te(i)
           tc(i)=te(i)
           yta(i)=0.
           dta(i)=0.
           tdf(i)=tde(i)
           tdc(i)=tde(i)
           ytda(i)=0.
           dtda(i)=0.
           uf(i)=ue(i)
           uc(i)=ue(i)
           yua(i)=0.
           dua(i)=0.
           vf(i)=ve(i)
           vc(i)=ve(i)
           yva(i)=0.
           dva(i)=0.
           pmslf(i)=pmsle(i)
           pmslc(i)=pmsle(i)
           ypma(i)=0.
           dpma(i)=0.
           altf(i)=alte(i)
           altc(i)=alte(i)
           yalta(i)=0.
           dalta(i)=0.
        enddo !i
c
c.....  Run the gross error check with the kalman xt estimate; there should be 
c.....  no missings in any "f" variable
c
c save all ncms
        do i=1,maxstaa
         ht(i)=ncmt(i)
         htd(i)=ncmtd(i)
         hu(i)=ncmu(i)
         hv(i)=ncmv(i)
         hp(i)=ncmpm(i)
         ha(i)=ncmalt(i)
        enddo
        
        print*,'Gross error checks: T'
        call qcset(t,tf,ncmt,qcstat,m,maxstaa,badflag,grosst)
        print*,' dewpoint '
        call qcset(td,tdf,ncmtd,qcstatd,m,maxstaa,badflag,grosstd)
        print*, 'U wind'
        call qcset(u,uf,ncmu,qcstau,m,maxstaa,badflag,grossuv)
        print*, 'V wind'
        call qcset(v,vf,ncmv,qcstav,m,maxstaa,badflag,grossuv)
        print*, 'PMSL'
        call qcset(pmsl,pmslf,ncmpm,qcstapm,m,maxstaa,badflag,grosspm)
        print*, 'Alt'
        call qcset(alt,altf,ncmalt,qcstal,m,maxstaa,badflag,grossal)
        if(iflag.eq.1) then ! data is missing so set all q**** to 1
         do i=1,maxstaa
          qcstat(i)=1
          qcstatd(i)=1
          qcstau(i)=1
          qcstav(i)=1
          qcstapm(i)=1
          qcstal(i)=1
          ncmt(i)=ht(i)+1
          ncmtd(i)=htd(i)+1
          ncmu(i)=hu(i)+1
          ncmv(i)=hv(i)+1
          ncmpm(i)=hp(i)+1
          ncmalt(i)=ha(i)+1 
         enddo
        endif
c
c.....  For missing observations use a weighted combination of the kalman (f)
c.....  estimate and the buddy anal(e) assuming half weight for 'half' cycles
c.....  consecutive missings (ncm)
c
        half=4.
c this ensures that the a arrays are filled for all maxstaa
        call subob(ta,te,tf,ncmt,half,maxstaa,m,badflag)
        call subob(tda,tde,tdf,ncmtd,half,maxstaa,m,badflag)
        call subob(ua,ue,uf,ncmu,half,maxstaa,m,badflag)
        call subob(va,ve,vf,ncmv,half,maxstaa,m,badflag)
        call subob(pmsla,pmsle,pmslf,ncmpm,half,maxstaa,m,badflag)
        call subob(alta,alte,altf,ncmalt,half,maxstaa,m,badflag)
c
c.....  Compute new weights based on new obs and obs appearing for the first time 
c
        call weights(mwt,dwt,ua,va,thetaa,lata,lona,cycle,maxstaa,m,
     &                300000.,30.)
        write(6,*) '...................................................'
        write(6,*) '...................................................'
        write(6,*) '............ITERATION ',it,' SUMMARY...............'
c
c.....  Set qc status and gross error check
c
        write(6,*) 'Temperature Stage 2**************'
c
c.....  Redo the model
c
        call kalmod(F,yta,byta,dta,ta,wmt,wot,wbt,offset,
     1       maxstaa,mwt,m)
        print*,'temp model before kalman call'
        do i=1,maxstaa
         write(6,3333) 'stn,i,F,o,b,m,bc ',stna(i),i,F(i,i),yta(i),
     &                  byta(i),dta(i),bmonster(i,ihr,1)
 3333 format(1x,a17,a5,i5,f10.7,3f7.2,f9.5)
        enddo
c
c.....  Prepare the fields for the kalman filter
c
        call perturb(ta,zot,ta,maxstaa,m,offset,on)      
        call perturb(tc,zot,tc,maxstaa,m,offset,on)      
c
c.....  Create the v and w matrices using the error history
c
        Call avgdiagerr(wit,B,c,w,maxstaa,ia,m,it-1)
        Call avgdiagerr(Vit,B,c,vv,maxstaa,ia,m,it-1)
        write(6,*) 'Kalman for TEMP'
        call kalman(F,tc,ta,pt,it, w,vv,xt,xtt,maxstaa,m,atime,stna)
c
c.....  Find ob difference
c
        call perturb(tc,zot,tc,maxstaa,m,offset,off)
        call perturb(t,tc,dt,maxstaa,m,0.,on)      
c
c.....  Fill the observed temp difference msngs from surrounding stns
c.....  assume that vertical lapse rate of dt is very small
c
        print*,'fill in observed temp change'
        call fill(stna,dt,lata,lona,eleva,thetaa,sl,
     &             maxstaa,m,0.5/1000.,-0.5/1000.,200000.,100.)          
        call replacemsng(zot,dt,maxstaa,1,m,1,badflag)
        write(6,*) 'ERROR PROCESSING FOR T'
        call perturb(xt,zot,xt,maxstaa,m,offset,off)      
c
c.....  Put kalman optimum in ta
c
        call replace(xt,ta,maxstaa,1,m,1)
        call perturb(xtt,zot,xtt,maxstaa,m,offset,off)      
        call replace(xtt,tf,maxstaa,1,m,1)
        call perturb(xt,tc,xt,maxstaa,m,0.,on)      
        call perturb(xtt,tc,xtt,maxstaa,m,0.,on)      
        call errorproc(dt,yta,byta,dta,xtt,xt,wr,vr,ar,maxstaa,m,qcstat
     & ,oberrt,wkt,wmt,wot,wbt,length,badthr,atime,icnt,nvar,1,i4time  
     &  ,tor,tcr,tmr,tob,tcb,tmb,totr,totb,bmonster,ihr)
        call perturb(tb,yta,tc,maxstaa,m,0.,off)
        call perturb(te,dta,tc,maxstaa,m,0.,off)
        call writemon(wr,wr,wr,wr,wr,wr,1,maxstaa,m,ncycles,wit,it)
        call writemon(vr,vr,vr,vr,vr,vr,1,maxstaa,m,ncycles,vit,it)
c
c.....  Write out observations
c.....  variable key: no suffix are raw obs, b=ob trend,
c.....  c= buddy,e=NWP, f=kalman fcst, a=Kalman optimum
c.....  Temperature.  .  .  .  .  .  .  ..  .  .  .  .  .  .  .  .  .  .
c
        write(6,*) 'OUTPUT for T '
        if(iflag.eq.1) 
     &   write(6,*) 'WARNING: ALL OBS MISSING!!!!OBS ASSUMED AS KALMAN'
        call writec(stn,maxstaa,m,' Stations   ',atime)
        call writev(tc,maxstaa,1,m, 'TA LAST CYC ',atime,on,0.)
        call writev(tf,maxstaa,1,m, ' XT from KAL  ',atime,on,0.)
        call writev(tr,maxstaa,1,m, ' RAW T OBS  ',atime,on,0.)
        call perturb(tf,t,ffg ,maxstaa,m,0.,on)
        call writev(ffg,maxstaa,1,m, 'XT - RAW T  ',atime,on,0.)
        call writev(tb,maxstaa,1,m, ' OB TREND T',atime,on,0.)
        call perturb(tc,byta,tc,maxstaa,m,0.,off)
        call writev(tc,maxstaa,1,m, ' BUD GUESS T',atime,on,0.)
        call writev(te,maxstaa,1,m, ' NWP GUESS T',atime,on,0.)
        call writev(ta,maxstaa,1,m,' KALMAN  T ',atime,on,0.)
        call writei(qcstat,maxstaa,m,
     &   'QC FAILURES: 100000 FX=XT,10000NWP,1000BTND,100TND,10GRS,1MSG'
     &   ,atime)
        call writev(ncmt,maxstaa,1,m,' NumConMsg ',atime,on,0.)
c
c.....  Redo the model
c
        write(6,*) 'Dewpoint    Stage 2**************'
        call kalmod(F,ytda,bytda,dtda,tda,wmtd,wotd,wbtd,offset,
     1              maxstaa,mwt,m)
        print*,'dewpoint model before kalman call'
        do i=1,maxstaa
        write(6,3333) 'stn,i,F,o,b,m,bc ',stna(i),i,F(i,i),ytda(i),
     &              bytda(i),dtda(i),bmonster(i,ihr,2)
        enddo
c
c.....  Prepare the fields for the kalman filter
c
        call perturb(tda,zot,tda,maxstaa,m,offset,on)      
        call perturb(tdc,zot,tdc,maxstaa,m,offset,on)      
c
c.....  Create the v and w matrices using the error history
c
        Call avgdiagerr(witd,B,c,w,maxstaa,ia,m,it-1)
        Call avgdiagerr(Vitd,B,c,vv,maxstaa,ia,m,it-1)
        write(6,*) 'Kalman for Dewpoint'
        call kalman(F,tdc,tda,ptd,it, w,vv,xt,xtt,maxstaa,m,
     &              atime,stna)
c
c.....  Find ob difference
c
        call perturb(tdc,zot,tdc,maxstaa,m,offset,off)
        call perturb(td,tdc,dtd,maxstaa,m,0.,on)      
c
c.....  Fill the observed temp difference msngs from surrounding stns
c.....  assume that vertical lapse rate of dt is very small
c
        print*,'fill in observed dewpoint change'
        call fill(stna,dtd,lata,lona,eleva,thetaa,sl,
     &             maxstaa,m,0.5/1000.,-0.5/1000.,200000.,100.)          
        call replacemsng(zot,dtd,maxstaa,1,m,1,badflag)
        write(6,*) 'ERROR PROCESSING FOR TD'
        call perturb(xt,zot,xt,maxstaa,m,offset,off)      
c
c.....  Put kalman optimum in ta
c
        call replace(xt,tda,maxstaa,1,m,1)
        call perturb(xtt,zot,xtt,maxstaa,m,offset,off)      
        call replace(xtt,tdf,maxstaa,1,m,1)
        call perturb(xt,tdc,xt,maxstaa,m,0.,on)      
        call perturb(xtt,tdc,xtt,maxstaa,m,0.,on)      
        call errorproc(dtd,ytda,bytda,dtda,xtt,xt,wr,vr,ar,maxstaa,m,
     & qcstatd,oberrtd,wktd,wmtd,wotd,wbtd,length,badthr,atime,icnt,
     &  nvar,2,i4time  
     &  ,tor,tcr,tmr,tob,tcb,tmb,totr,totb,bmonster,ihr)
        call perturb(tdb,ytda,tdc,maxstaa,m,0.,off)
        call perturb(tde,dtda,tdc,maxstaa,m,0.,off)
        call writemon(wr,wr,wr,wr,wr,wr,1,maxstaa,m,ncycles,witd,it)
        call writemon(vr,vr,vr,vr,vr,vr,1,maxstaa,m,ncycles,vitd,it)
c
c.....  Write out observations
c.....  variable key: no suffix are raw obs, b=ob trend,c prev cycle
c.....  d= buddy,e=NWP, a=Kalman optimum
c.....  Dewpoint....
c
        write(6,*) 'OUTPUT for TD'
        if(iflag.eq.1) 
     &   write(6,*) 'WARNING: ALL OBS MISSING!!!!OBS ASSUMED AS KALMAN'
        call writec(stn,maxstaa,m,' Stations   ',atime)
        call writev(tdc,maxstaa,1,m, 'TDA LAST CYC ',atime,on,0.)
        call writev(tdf,maxstaa,1,m, ' XT from KAL  ',atime,on,0.)
        call writev(tdr,maxstaa,1,m, ' RAW TD OBS ',atime,on,0.)
        call perturb(tdf,td,ffg ,maxstaa,m,0.,on)
        call writev(ffg,maxstaa,1,m, 'XTD- RAW TD ',atime,on,0.)
        call writev(tdb,maxstaa,1,m, ' OB TREND TD',atime,on,0.)
        call perturb(tdc,bytda,tdc,maxstaa,m,0.,off)
        call writev(tdc,maxstaa,1,m, ' BUD GUESS TD',atime,on,0.)
        call writev(tde,maxstaa,1,m, ' NWP GUESS TD',atime,on,0.)
        call writev(tda,maxstaa,1,m,' KALMAN  TD',atime,on,0.)
        call writei(qcstatd,maxstaa,m,
     &   'QC FAILURES: 100000 FX=XT,10000NWP,100OBTND,100BUD,10GRS,1MSG' 
     &   ,atime)
        call writev(ncmtd,maxstaa,1,m,' NumConMsg ',atime,on,0.)
c
c.....  Redo the model
c
        write(6,*) 'U-Wind      Stage 2**************'
        call kalmod(F,yua,byua,dua,ua,wmu,wou,wbu,offset,
     1              maxstaa,mwt,m)
        print*,'u-wind model before kalman call '
        do i=1,maxstaa
         write(6,3333) 'stn,i,F,o,b,m,bc ',stna(i),i,F(i,i),yua(i),
     &                  byua(i),dua(i),bmonster(i,ihr,3)
        enddo
c
c.....  Prepare the fields for the kalman filter
c
        call perturb(ua,zot,ua,maxstaa,m,offset,on)      
        call perturb(uc,zot,uc,maxstaa,m,offset,on)      
c
c.....  Create the v and w matrices using the error history
c
        Call avgdiagerr(witu,B,c,W,maxstaa,ia,m,it-1)
        Call avgdiagerr(Vitu,B,c,vV,maxstaa,ia,m,it-1)
        write(6,*) 'Kalman for U'
        call kalman(F,uc,ua,pu,it, w,vv,xt,xtt,maxstaa,m,atime,stna)
c
c.....  Find ob difference
c
        call perturb(uc,zot,uc,maxstaa,m,offset,off)
        call perturb(u,uc,du,maxstaa,m,0.,on)      
c
c.....  Fill the observed temp difference msngs from surrounding stns
c.....  assume that vertical lapse rate of du is very small
c
        print*,'fill in observed U change'
        call fill(stna,du,lata,lona,eleva,thetaa,sl,
     &             maxstaa,m,0.5/1000.,-0.5/1000.,200000.,50.)          
        call replacemsng(zot,du,maxstaa,1,m,1,badflag)
        write(6,*) 'ERROR PROCESSING FOR U'
        call perturb(xt,zot,xt,maxstaa,m,offset,off)      
c
c.....  Put kalman optimum in ta
c
        call replace(xt,ua,maxstaa,1,m,1)
        call perturb(xtt,zot,xtt,maxstaa,m,offset,off)      
        call replace(xtt,uf,maxstaa,1,m,1)
        call perturb(xt,uc,xt,maxstaa,m,0.,on)      
        call perturb(xtt,uc,xtt,maxstaa,m,0.,on)      
        call errorproc(du,yua,byua,dua,xtt,xt,wr,vr,ar,maxstaa,m,qcstau
     & ,oberruv,wku,wmu,wou,wbu,length,badthr,atime,icnt,nvar,3,i4time  
     &  ,tor,tcr,tmr,tob,tcb,tmb,totr,totb,bmonster,ihr)
        call perturb(ub,yua,uc,maxstaa,m,0.,off)
        call perturb(ue,dua,uc,maxstaa,m,0.,off)
        call writemon(wr,wr,wr,wr,wr,wr,1,maxstaa,m,ncycles,witu,it)
        call writemon(vr,vr,vr,vr,vr,vr,1,maxstaa,m,ncycles,vitu,it)
c
c.....  Write out observations
c.....  variable key: no suffix are raw obs, b=ob trend,c prev cycle
c.....  d= buddy,e=NWP, a=Kalman optimum
c.....  U-Wind 
c
        write(6,*) 'OUTPUT for U '
        if(iflag.eq.1) 
     &   write(6,*) 'WARNING: ALL OBS MISSING!!!!OBS ASSUMED AS KALMAN'
        call writec(stn,maxstaa,m,' Stations   ',atime)
        call writev(uc,maxstaa,1,m, 'UA LAST CYC ',atime,on,0.)
        call writev(uf,maxstaa,1,m, ' XT from KAL  ',atime,on,0.)
        call writev(ur,maxstaa,1,m, ' RAW U OBS  ',atime,on,0.)
        call perturb(uf,u,ffg ,maxstaa,m,0.,on)
        call writev(ffg,maxstaa,1,m, 'XU - RAW U  ',atime,on,0.)
        call writev(ub,maxstaa,1,m, ' OB TREND U',atime,on,0.)
        call perturb(uc,byua,uc,maxstaa,m,0.,off)
        call writev(uc,maxstaa,1,m, ' BUD GUESS U',atime,on,0.)
        call writev(ue,maxstaa,1,m, ' NWP GUESS U',atime,on,0.)
        call writev(ua,maxstaa,1,m,' KALMAN  U ',atime,on,0.)
        call writei(qcstau,maxstaa,m,
     &   'QC FAILURES: 100000 FX=XT,10000NWP,100OBTND,100BUD,10GRS,1MSG' 
     &   ,atime)
        call writev(ncmu,maxstaa,1,m,' NumConMsg ',atime,on,0.)
        write(6,*) 'V- Wind     Stage 2**************'
c
c.....  Redo the model
c
        call kalmod(F,yva,byva,dva,va,wmv,wov,wbv,offset,
     1              maxstaa,mwt,m)
        print*,'v-wind model'
        do i=1,maxstaa
         write(6,3333) 'stn,i,F,o,b,m,bc ',stna(i),i,F(i,i),yva(i),
     &                  byva(i),dva(i),bmonster(i,ihr,4)
        enddo
c
c.....  Prepare the fields for the kalman filter
c
        call perturb(va,zot,va,maxstaa,m,offset,on)      
        call perturb(vc,zot,vc,maxstaa,m,offset,on)      
c
c.....  Create the v and w matrices using the error history
c
        Call avgdiagerr(witv,B,c,W,maxstaa,ia,m,it-1)
        Call avgdiagerr(Vitv,B,c,vV,maxstaa,ia,m,it-1)
        write(6,*) 'Kalman for V   '
        call kalman(F,vc,va,pv,it, w,vv,xt,xtt,maxstaa,m,atime,stna)
c
c.....  Find ob difference
c
        call perturb(vc,zot,vc,maxstaa,m,offset,off)
        call perturb(v,vc,dv,maxstaa,m,0.,on)      
c
c.....  Fill the observed temp difference msngs from surrounding stns
c.....  assume that vertical lapse rate of dt is very small
c
        print*,'fill in observed V change'
        call fill(stna,dv,lata,lona,eleva,thetaa,sl,
     &             maxstaa,m,0.5/1000.,-0.5/1000.,200000.,50.)          
        call replacemsng(zot,dv,maxstaa,1,m,1,badflag)
        write(6,*) 'ERROR PROCESSING FOR V'
        call perturb(xt,zot,xt,maxstaa,m,offset,off)      
c
c.....  Put kalman optimum in ta
c
        call replace(xt,va,maxstaa,1,m,1)
        call perturb(xtt,zot,xtt,maxstaa,m,offset,off)      
        call replace(xtt,vf,maxstaa,1,m,1)
        call perturb(xt,vc,xt,maxstaa,m,0.,on)      
        call perturb(xtt,vc,xtt,maxstaa,m,0.,on)      
        call errorproc(dv,yva,byva,dva,xtt,xt,wr,vr,ar,maxstaa,m,qcstav
     & ,oberruv,wkv,wmv,wov,wbv,length,badthr,atime,icnt,nvar,4,i4time  
     &  ,tor,tcr,tmr,tob,tcb,tmb,totr,totb,bmonster,ihr)
        call perturb(vb,yva,vc,maxstaa,m,0.,off)
        call perturb(ve,dva,vc,maxstaa,m,0.,off)
        call writemon(wr,wr,wr,wr,wr,wr,1,maxstaa,m,ncycles,witv,it)
        call writemon(vr,vr,vr,vr,vr,vr,1,maxstaa,m,ncycles,vitv,it)
c
c.....  Write out observations
c.....  variable key: no suffix are raw obs, b=ob trend,c prev cycle
c.....  d= buddy,e=NWP, a=Kalman optimum
c.....  V wind
c
        write(6,*) 'OUTPUT for V '
        if(iflag.eq.1) 
     &   write(6,*) 'WARNING: ALL OBS MISSING!!!!OBS ASSUMED AS KALMAN'
        call writec(stn,maxstaa,m,' Stations   ',atime)
        call writev(vc,maxstaa,1,m, 'VA LAST CYC ',atime,on,0.)
        call writev(vf,maxstaa,1,m, ' XT from KAL  ',atime,on,0.)
        call writev(vr,maxstaa,1,m, ' RAW V OBS  ',atime,on,0.)
        call perturb(vf,v,ffg ,maxstaa,m,0.,on)
        call writev(ffg,maxstaa,1,m, 'XV - RAW V  ',atime,on,0.)
        call writev(vb,maxstaa,1,m, ' OB TREND V',atime,on,0.)
        call perturb(vc,byva,vc,maxstaa,m,0.,off)
        call writev(vc,maxstaa,1,m, ' BUD GUESS V',atime,on,0.)
        call writev(ve,maxstaa,1,m, ' NWP GUESS V',atime,on,0.)
        call writev(va,maxstaa,1,m,' KALMAN  V ',atime,on,0.)
        call writei(qcstav,maxstaa,m,
     &   'QC FAILURES: 100000 FX=XT,10000NWP,100OBTND,100BUD,10GRS,1MSG' 
     &   ,atime)
        call writev(ncmv,maxstaa,1,m,' NumConMsg ',atime,on,0.)
        write(6,*) 'MSL PressureStage 2**************'
c
c.....  Redo the model
c
        call kalmod(F,ypma,bypma,dpma,pmsla,wmp,wop,wbp,offset,
     1              maxstaa,mwt,m)
        print*,'pmsl   model before kalman call'
        do i=1,maxstaa
        write(6,3333) 'stn,i,F,o,b,m,bc ',stna(i),i,F(i,i),ypma(i),
     &             bypma(i),dpma(i),bmonster(i,ihr,5)
        enddo
c
c.....  Prepare the fields for the kalman filter
c
        call perturb(pmsla,zot,pmsla,maxstaa,m,offset,on)      
        call perturb(pmslc,zot,pmslc,maxstaa,m,offset,on)      
c
c.....  Create the v and w matrices using the error history
c
        call avgdiagerr(witpm,B,c,W,maxstaa,ia,m,it-1)
        Call avgdiagerr(Vitpm,B,c,vV,maxstaa,ia,m,it-1)
        write(6,*) 'Kalman for PMSL'
        call kalman(F,pmslc,pmsla,ppm,it, w,vv,xt,xtt,maxstaa,m,
     &            atime,stna)
c
c.....  Find ob difference
c
        call perturb(pmslc,zot,pmslc,maxstaa,m,offset,off)
        call perturb(pmsl,pmslc,dpm,maxstaa,m,0.,on)      
c
c.....  Fill the observed temp difference msngs from surrounding stns
c.....  assume that vertical lapse rate of dt is very small
c
        print*,'fill in observed PMSL change'
        call fill(stna,dpm,lata,lona,eleva,thetaa,sl,
     &             maxstaa,m,0.5/1000.,-0.5/1000.,200000.,1000.)          
        call replacemsng(zot,dpm,maxstaa,1,m,1,badflag)
        write(6,*) 'ERROR PROCESSING FOR PMSL'
        call perturb(xt,zot,xt,maxstaa,m,offset,off)      
c
c.....  Put kalman optimum in ta
c
        call replace(xt,pmsla,maxstaa,1,m,1)
        call perturb(xtt,zot,xtt,maxstaa,m,offset,off)      
        call replace(xtt,pmslf,maxstaa,1,m,1)
        call perturb(xt,pmslc,xt,maxstaa,m,0.,on)      
        call perturb(xtt,pmslc,xtt,maxstaa,m,0.,on)      
        call errorproc(dpm,ypma,bypma,dpma,xtt,xt,wr,vr,ar,maxstaa,
     &  m,qcstapm,oberrpm,wkp,wmp,wop,wbp,length,badthr,atime,icnt,nvar
     &  ,5,i4time  
     &  ,tor,tcr,tmr,tob,tcb,tmb,totr,totb,bmonster,ihr)
        call perturb(pmslb,ypma,pmslc,maxstaa,m,0.,off)
        call perturb(pmsle,dpma,pmslc,maxstaa,m,0.,off)
        call perturb(pmslc,bypma,pmslc,maxstaa,m,0.,off)
        call writemon(wr,wr,wr,wr,wr,wr,1,maxstaa,m,ncycles,witpm,it)
        call writemon(vr,vr,vr,vr,vr,vr,1,maxstaa,m,ncycles,vitpm,it)
c
c.....  Write out observations
c.....  variable key: no suffix are raw obs, b=ob trend,c prev cycle
c.....  d= buddy,e=NWP, a=Kalman optimum
c.....  Temperature
c
        write(6,*) 'OUTPUT for PMSL '
        if(iflag.eq.1) 
     &   write(6,*) 'WARNING: ALL OBS MISSING!!!!OBS ASSUMED AS KALMAN'
        call writec(stn,maxstaa,m,' Stations   ',atime)
        call writev(pmslc,maxstaa,1,m, 'PMSLA LAST CYC ',atime,on,1000.)
        call writev(pmslf,maxstaa,1,m, ' XT from KAL  ',atime,on,1000.)
        call writev(pmslr,maxstaa,1,m, ' RAW PMSL OBS  ',atime,on,1000.)
        call perturb(pmslf,pmsl,ffg ,maxstaa,m,0.,on)
        call writev(ffg,maxstaa,1,m, 'XPM- RAW PM ',atime,on,0.)
        call writev(pmslb,maxstaa,1,m, ' OB TREND PMSL',atime,on,1000.)
        call perturb(pmslc,bypma,pmslc,maxstaa,m,0.,off)
        call writev(pmslc,maxstaa,1,m, ' BUD GUESS PMSL',atime,on,1000.)
        call writev(pmsle,maxstaa,1,m, ' NWP GUESS PMSL',atime,on,1000.)
        call writev(pmsla,maxstaa,1,m,' KALMAN  PMSL ',atime,on,1000.)
        call writei(qcstapm,maxstaa,m,
     &   'QC FAILURES: 100000 FX=XT,10000NWP,100OBTND,100BUD,10GRS,1MSG' 
     &   ,atime)
        call writev(ncmpm,maxstaa,1,m,' NumConMsg ',atime,on,0.)
        write(6,*) 'Altimeter   Stage 2**************'
c
c.....  Redo the model
c
        call kalmod(F,yalta,byalta,dalta,alta,wma,woa,wba,offset,
     1              maxstaa,mwt,m)
        print*,'alstg  model before kalman call'
        do i=1,maxstaa
         write(6,3333) 'stn,i,F,o,b,m,bc ',stna(i),i,F(i,i),yalta(i),
     &            byalta(i),dalta(i),bmonster(i,ihr,6)
        enddo
c
c.....  Prepare the fields for the kalman filter
c
        call perturb(alta,zot,alta,maxstaa,m,offset,on)      
        call perturb(altc,zot,altc,maxstaa,m,offset,on)      
c
c.....  Create the v and w matrices using the error history
c
        Call avgdiagerr(wital,B,c,W,maxstaa,ia,m,it-1)
        Call avgdiagerr(Vital,B,c,vV,maxstaa,ia,m,it-1)
        write(6,*) 'Kalman for ALSTG'
        call kalman(F,altc,alta,pal,it, w,vv,xt,xtt,maxstaa,m,
     &                 atime,stna)
c
c.....  Find ob difference
c
        call perturb(altc,zot,altc,maxstaa,m,offset,off)
        call perturb(alt,altc,dalt,maxstaa,m,0.,on)      
c
c.....  Fill the observed temp difference msngs from surrounding stns
c.....  assume that vertical lapse rate of dt is very small
c
        print*,'fill in observed ALSG change'
        call fill(stna,dalt,lata,lona,eleva,thetaa,sl,
     &             maxstaa,m,0.5/1000.,-0.5/1000.,200000.,1000.)          
        call replacemsng(zot,dalt,maxstaa,1,m,1,badflag)
        write(6,*) 'ERROR PROCESSING FOR ALSTG'
        call perturb(xt,zot,xt,maxstaa,m,offset,off)      
c
c.....  Put kalman optimum in ta
c
        call replace(xt,alta,maxstaa,1,m,1)
        call perturb(xtt,zot,xtt,maxstaa,m,offset,off)      
        call replace(xtt,altf,maxstaa,1,m,1)
        call perturb(xt,altc,xt,maxstaa,m,0.,on)      
        call perturb(xtt,altc,xtt,maxstaa,m,0.,on)      
        call errorproc(dalt,yalta,byalta,dalta,xtt,xt,wr,vr,ar,maxstaa,
     & m,qcstal,oberral,wka,wma,woa,wba,length,badthr,atime,icnt,nvar,6
     &     ,i4time  
     &  ,tor,tcr,tmr,tob,tcb,tmb,totr,totb,bmonster,ihr)
        call perturb(altb,yalta,altc,maxstaa,m,0.,off)
        call perturb(alte,dalta,altc,maxstaa,m,0.,off)
        call writemon(wr,wr,wr,wr,wr,wr,1,maxstaa,m,ncycles,wital,it)
        call writemon(vr,vr,vr,vr,vr,vr,1,maxstaa,m,ncycles,vital,it)
c
c.....  Write out observations
c.....  variable key: no suffix are raw obs, b=ob trend,c prev cycle
c.....  d= buddy,e=NWP, a=Kalman optimum
c.....  Temperature
c
        write(6,*) 'OUTPUT for ALSTG '
        if(iflag.eq.1) 
     &   write(6,*) 'WARNING: ALL OBS MISSING!!!!OBS ASSUMED AS KALMAN'
        call writec(stn,maxstaa,m,' Stations   ',atime)
        call writev(altc,maxstaa,1,m, 'ALTA LAST CYC  ',atime,on,1000.)
        call writev(altf,maxstaa,1,m, ' XT from KAL  ',atime,on,1000.)
        call writev(altr,maxstaa,1,m, ' RAW ALSTG OBS  ',atime,on,1000.)
        call perturb(altf,alt,ffg ,maxstaa,m,0.,on)
        call writev(ffg,maxstaa,1,m, 'XAL- RAW AL ',atime,on,0.)
        call writev(altb,maxstaa,1,m, ' OB TREND ALSTG',atime,on,1000.)
        call perturb(altc,byalta,altc,maxstaa,m,0.,off)
        call writev(altc,maxstaa,1,m, ' BUD GUESS ALSTG',atime,on,1000.)
        call writev(alte,maxstaa,1,m, ' NWP GUESS ALSTG',atime,on,1000.)
        call writev(alta,maxstaa,1,m,' KALMAN ALSTG  ',atime,on,1000.)
        call writei(qcstal ,maxstaa,m,
     &   'QC FAILURES: 100000 FX=XT,10000NWP,100OBTND,100BUD,10GRS,1MSG' 
     &   ,atime)
        call writev(ncmalt,maxstaa,1,m,' NumConMsg ',atime,on,0.)
        print *,' '
        write(6,3004) atime,it
        write(6,3006) maxstaa
        print *,' '


 3004   format (1x,a24,3x, 'it= ',i12)
 3005   format(1x,a4,3x,2f7.2,f7.0,3x,6(i7,6f6.1))
 3006   format(1x,i3)

        write(6,*) '...................................................'
        write(6,*) '...................................................'
        do i=1,maxstaa
           if( t(i) .ne. badflag)  t(i) = ((  t(i) - 32.) * fon) +273.15!conv F to K
           if(ta(i) .ne. badflag) ta(i) = (( ta(i) - 32.) * fon) +273.15  ! conv F to K
           if(tb(i) .ne. badflag) tb(i) = (( tb(i) - 32.) * fon) +273.15  ! conv F to K
           if(tc(i) .ne. badflag) tc(i) = (( tc(i) - 32.) * fon) +273.15  ! conv F to K
           if(te(i) .ne. badflag) te(i) = (( te(i) - 32.) * fon) +273.15  ! conv F to K
           if(tf(i) .ne. badflag) tf(i) = (( tf(i) - 32.) * fon) +273.15  ! conv F to K
c
           if( td(i) .ne. badflag)  td(i) = ((td(i)-32.)  * fon) +273.15  ! conv F to K
           if(tda(i) .ne. badflag) tda(i) = ((tda(i)-32.) * fon) +273.15  ! conv F to K
           if(tdb(i) .ne. badflag) tdb(i) = ((tdb(i)-32.) * fon) +273.15  ! conv F to K
           if(tdc(i) .ne. badflag) tdc(i) = ((tdc(i)-32.) * fon) +273.15  ! conv F to K
           if(tde(i) .ne. badflag) tde(i) = ((tde(i)-32.) * fon) +273.15  ! conv F to K
           if(tdf(i) .ne. badflag) tdf(i) = ((tdf(i)-32.) * fon) +273.15  ! conv F to K

           if( u(i) .ne. badflag)  u(i) =  u(i) * 0.514791          ! conv kt to m/s
           if(ua(i) .ne. badflag) ua(i) = ua(i) * 0.514791          ! conv kt to m/s
           if(ub(i) .ne. badflag) ub(i) = ub(i) * 0.514791          ! conv kt to m/s
           if(uc(i) .ne. badflag) uc(i) = uc(i) * 0.514791          ! conv kt to m/s
           if(ue(i) .ne. badflag) ue(i) = ue(i) * 0.514791          ! conv kt to m/s
           if(uf(i) .ne. badflag) uf(i) = uf(i) * 0.514791          ! conv kt to m/s
c
           if( v(i) .ne. badflag)  v(i) =  v(i)* 0.514791           ! conv kt to m/s
           if(va(i) .ne. badflag) va(i) = va(i)* 0.514791           ! conv kt to m/s
           if(vb(i) .ne. badflag) vb(i) = vb(i)* 0.514791           ! conv kt to m/s
           if(vc(i) .ne. badflag) vc(i) = vc(i)* 0.514791           ! conv kt to m/s
           if(ve(i) .ne. badflag) ve(i) = ve(i)* 0.514791           ! conv kt to m/s
           if(vf(i) .ne. badflag) vf(i) = vf(i)* 0.514791           ! conv kt to m/s
c
           if( pmsl(i) .ne. badflag)  pmsl(i) =  pmsl(i) * 100.     ! conv mb to Pa
           if(pmsla(i) .ne. badflag) pmsla(i) = pmsla(i) * 100.     ! conv mb to Pa  
           if(pmslb(i) .ne. badflag) pmslb(i) = pmslb(i) * 100.     ! conv mb to Pa
           if(pmslc(i) .ne. badflag) pmslc(i) = pmslc(i) * 100.     ! conv mb to Pa
           if(pmsle(i) .ne. badflag) pmsle(i) = pmsle(i) * 100.     ! conv mb to Pa
           if(pmslf(i) .ne. badflag) pmslf(i) = pmslf(i) * 100.     ! conv mb to Pa
c
           if( alt(i) .ne. badflag)  alt(i) =  alt(i) * 100.        ! conv mb to Pa
           if(alta(i) .ne. badflag) alta(i) = alta(i) * 100.        ! conv mb to Pa  
           if(altb(i) .ne. badflag) altb(i) = altb(i) * 100.        ! conv mb to Pa
           if(altc(i) .ne. badflag) altc(i) = altc(i) * 100.        ! conv mb to Pa
           if(alte(i) .ne. badflag) alte(i) = alte(i) * 100.        ! conv mb to Pa
           if(altf(i) .ne. badflag) altf(i) = altf(i) * 100.        ! conv mb to Pa
c
        enddo !i
c        
        ext_out = 'lsq'
        call get_directory('lsq', dir_out, len)
c        
        call write_qc_cdf(i4time, dir_out, ext_out, m, maxstaa, 
     1       stna, providera, reptypea, lata, lona, eleva,
     1       qcstat,  t,    tb,    ta,    tc,    te,    tf,  
     1       qcstatd, td,   tdb,   tda,   tdc,   tde,   tdf,  
     1       qcstau, u,    ub,    ua,    uc,    ue,    uf,  
     1                v,    vb,    va,    vc,    ve,    vf,  
     1       qcstapm, pmsl, pmslb, pmsla, pmslc, pmsle, pmslf, 
     1       qcstal,  alt,  altb,  alta,  altc, alte,   altf, 
     1       status)
c
c.....  Now convert the units back for other storage files.
c
        do i=1,maxstaa
c
           if( t(i) .ne. badflag)  t(i) = ((  t(i) - 273.15) * nof) +32.    ! conv F to K
           if(ta(i) .ne. badflag) ta(i) = (( ta(i) - 273.15) * nof) +32.    ! conv F to K
           if(tb(i) .ne. badflag) tb(i) = (( tb(i) - 273.15) * nof) +32.    ! conv F to K
           if(tc(i) .ne. badflag) tc(i) = (( tc(i) - 273.15) * nof) +32.    ! conv F to K
           if(te(i) .ne. badflag) te(i) = (( te(i) - 273.15) * nof) +32.    ! conv F to K
           if(tf(i) .ne. badflag) tf(i) = (( tf(i) - 273.15) * nof) +32.    ! conv F to K
c
           if( td(i) .ne. badflag)  td(i) = ((  td(i) -273.15)*nof) +32.    ! conv F to K
           if(tda(i) .ne. badflag) tda(i) = (( tda(i) -273.15)*nof) +32.    ! conv F to K
           if(tdb(i) .ne. badflag) tdb(i) = (( tdb(i) -273.15)*nof) +32.    ! conv F to K
           if(tdc(i) .ne. badflag) tdc(i) = (( tdc(i) -273.15)*nof) +32.    ! conv F to K
           if(tde(i) .ne. badflag) tde(i) = (( tde(i) -273.15)*nof) +32.    ! conv F to K
           if(tdf(i) .ne. badflag) tdf(i) = (( tdf(i) -273.15)*nof) +32.    ! conv F to K

           if( u(i) .ne. badflag)  u(i) =  u(i) * 1.94254         ! conv m/s to kt
           if(ua(i) .ne. badflag) ua(i) = ua(i) * 1.94254         ! conv m/s to kt
           if(ub(i) .ne. badflag) ub(i) = ub(i) * 1.94254         ! conv m/s to kt
           if(uc(i) .ne. badflag) uc(i) = uc(i) * 1.94254         ! conv m/s to kt
           if(ue(i) .ne. badflag) ue(i) = ue(i) * 1.94254         ! conv m/s to kt
           if(uf(i) .ne. badflag) uf(i) = uf(i) * 1.94254         ! conv m/s to kt
c
           if( v(i) .ne. badflag)  v(i) =  v(i) * 1.94254         ! conv m/s to kt
           if(va(i) .ne. badflag) va(i) = va(i) * 1.94254         ! conv m/s to kt
           if(vb(i) .ne. badflag) vb(i) = vb(i) * 1.94254         ! conv m/s to kt
           if(vc(i) .ne. badflag) vc(i) = vc(i) * 1.94254         ! conv m/s to kt
           if(ve(i) .ne. badflag) ve(i) = ve(i) * 1.94254         ! conv m/s to kt
           if(vf(i) .ne. badflag) vf(i) = vf(i) * 1.94254         ! conv m/s to kt
c
           if( pmsl(i) .ne. badflag)  pmsl(i) =  pmsl(i) * 0.01   ! conv Pa to mb
           if(pmsla(i) .ne. badflag) pmsla(i) = pmsla(i) * 0.01   ! conv Pa to mb
           if(pmslb(i) .ne. badflag) pmslb(i) = pmslb(i) * 0.01   ! conv Pa to mb
           if(pmslc(i) .ne. badflag) pmslc(i) = pmslc(i) * 0.01   ! conv Pa to mb
           if(pmsle(i) .ne. badflag) pmsle(i) = pmsle(i) * 0.01   ! conv Pa to mb
           if(pmslf(i) .ne. badflag) pmslf(i) = pmslf(i) * 0.01   ! conv Pa to mb
c
           if( alt(i) .ne. badflag)  alt(i) =  alt(i) * 0.01      ! conv Pa to mb
           if(alta(i) .ne. badflag) alta(i) = alta(i) * 0.01      ! conv Pa to mb
           if(altb(i) .ne. badflag) altb(i) = altb(i) * 0.01      ! conv Pa to mb
           if(altc(i) .ne. badflag) altc(i) = altc(i) * 0.01      ! conv Pa to mb
           if(alte(i) .ne. badflag) alte(i) = alte(i) * 0.01      ! conv Pa to mb
           if(altf(i) .ne. badflag) altf(i) = altf(i) * 0.01      ! conv Pa to mb
c
        enddo !i
c
c.....  Write out the updated monster file, and the other storage file.
c.....  First, check the arrays for bad values.
c
        ifound_nan = 0
        nanvar = '          '
        print *,' Checking monster.dat arrays for NaNs...'
c
        call checknan_2d(vit,m,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'vit'
           go to 800
        endif
        call checknan_2d(vitd,m,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'vitd'
           go to 800
        endif
        call checknan_2d(vitu,m,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'vitu'
           go to 800
        endif
        call checknan_2d(vitv,m,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'vitv'
           go to 800
        endif
        call checknan_2d(vitpm,m,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'vitpm'
           go to 800
        endif
        call checknan_2d(vital,m,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'vital'
           go to 800
        endif
c
        call checknan_2d(wit,m,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'wit'
           go to 800
        endif
        call checknan_2d(witd,m,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'witd'
           go to 800
        endif
        call checknan_2d(witu,m,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'witu'
           go to 800
        endif
        call checknan_2d(witv,m,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'witv'
           go to 800
        endif
        call checknan_2d(witpm,m,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'witpm'
           go to 800
        endif
        call checknan_2d(wital,m,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'wital'
           go to 800
        endif
c
        call checknan_2d(pt,m,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'pt'
           go to 800
        endif
        call checknan_2d(ptd,m,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'ptd'
           go to 800
        endif
        call checknan_2d(pu,m,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'pu'
           go to 800
        endif
        call checknan_2d(pv,m,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'pv'
           go to 800
        endif
        call checknan_2d(ppm,m,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'ppm'
           go to 800
        endif
        call checknan_2d(pal,m,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'pal'
           go to 800
        endif
c
        do i=1,m
        do j=1,ncycles
        do k=1,nvar
           call checknan_pt(monster(i,j,k), nan_flag)
           if(nan_flag .ne. 1) then
              ifound_nan = 1
              nanvar = 'monster'
              go to 800
           endif
        enddo !k
        enddo !j
        enddo !i    
c
        call checknan_1d(lata,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'lata'
           go to 800
        endif
        call checknan_1d(lona,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'lona'
           go to 800
        endif
        call checknan_1d(eleva,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'eleva'
           go to 800
        endif
c
        call checknan_1d(ta,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'ta'
           go to 800
        endif
        call checknan_1d(tda,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'tda'
           go to 800
        endif
        call checknan_1d(ua,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'ua'
           go to 800
        endif
        call checknan_1d(va,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'va'
           go to 800
        endif
        call checknan_1d(pmsla,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'pmslaa'
           go to 800
        endif
        call checknan_1d(alta,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'alta'
           go to 800
        endif
c
        call checknan_1d(wot,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'wot'
           go to 800
        endif
        call checknan_1d(wotd,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'wotd'
           go to 800
        endif
        call checknan_1d(wou,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'wou'
           go to 800
        endif
        call checknan_1d(wov,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'wov'
           go to 800
        endif
        call checknan_1d(wop,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'wop'
           go to 800
        endif
        call checknan_1d(woa,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'woa'
           go to 800
        endif
c
        call checknan_1d(wbt,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'wbt'
           go to 800
        endif
        call checknan_1d(wbtd,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'wbtd'
           go to 800
        endif
        call checknan_1d(wbu,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'wbu'
           go to 800
        endif
        call checknan_1d(wbv,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'wbv'
           go to 800
        endif
        call checknan_1d(wbp,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'wbp'
           go to 800
        endif
        call checknan_1d(wba,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'wba'
           go to 800
        endif
c
        call checknan_1d(wmt,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'wmt'
           go to 800
        endif
        call checknan_1d(wmtd,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'wmtd'
           go to 800
        endif
        call checknan_1d(wmu,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'wmu'
           go to 800
        endif
        call checknan_1d(wmv,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'wmv'
           go to 800
        endif
        call checknan_1d(wmp,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'wmp'
           go to 800
        endif
        call checknan_1d(wma,m,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'wma'
           go to 800
        endif
c
        call checknan_1d(tor,nvar,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'tor'
           go to 800
        endif
        call checknan_1d(tcr,nvar,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'tcr'
           go to 800
        endif
        call checknan_1d(tmr,nvar,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'tmr'
           go to 800
        endif
        call checknan_1d(tob,nvar,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'tob'
           go to 800
        endif
        call checknan_1d(tcb,nvar,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'tcb'
           go to 800
        endif
        call checknan_1d(tmb,nvar,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'tmb'
           go to 800
        endif
c
        call checknan_1d(totr,nvar,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'totr'
           go to 800
        endif
        call checknan_1d(totb,nvar,nan_flag)
        if(nan_flag .ne. 1) then
           ifound_nan = 1
           nanvar = 'totb'
           go to 800
        endif
c
c.....  If no NaN's, write out the monster.dat and monster_jr.dat
c.....  files.  If there are, don't.
c
 800    continue
        if(ifound_nan .eq. 0) then
           call s_len(monfile, len)
	   print *,' No NaNs found...'
           print*,'Writing file ',monfile(1:len)
           open(16,file=monfile(1:len),
     &          form='unformatted',status='unknown')
           write(16) atime_cur, i4time
           write(16) vit,vitd,vitu,vitv,vitpm,vital
           write(16) wit,witd,witu,witv,witpm,wital
           write(16) pt,ptd,pu,pv,ppm,pal
           write(16) monster,bmonster
           write(16) stna,lata,lona,eleva
           write(16) providera, reptypea
           write(16) ta,tda,ua,va,pmsla,alta
           write(16) ncmt,ncmtd,ncmu,ncmv,ncmpm,ncmalt
           write(16) qcstat,qcstatd,qcstau,qcstav,qcstapm,qcstal
           write(16) qmonster
           write(16) it,maxstaa
           close(16)
           print *,' Done.'
c
           call s_len(monjrfile, len)
           print*,'Writing file ',monjrfile(1:len)
           open(19,file=monjrfile(1:len),
     &          form='unformatted',status='unknown')
c
           write(19) wot,wotd,wou,wov,wop,woa
           write(19) wbt,wbtd,wbu,wbv,wbp,wba
           write(19) wmt,wmtd,wmu,wmv,wmp,wma
           write(19) tor,tcr,tmr,tob,tcb,tmb
           write(19) totr,totb
           close(19)
           print *,' Done.'
c
        else
c
           print *,' **WARNING. Found a NaN in monster/monster_jr.dat ',
     &             'variable: ', nanvar
           print *,'            Skipping monster/monster_jr.dat writes.'
           print *,'            QC will automatically restart          '
c
        endif
c
        print*,' Kalman complete'
c
c.....  Write out a LSO file using the Kalman estimates
c.....  ...call it LSO_QC
c
        do i=1,maxstaa
           ii = index(i)
           rtime = float( time(ii) )
           stations_out(i) = '                    '
           stations_out(i)(1:5) = stna(i)(1:5)
           store_1(i,1) = lat(i)          ! station latitude
           store_1(i,2) = lon(i)          ! station longitude
           store_1(i,3) = elev(i)         ! station elevation
           store_1(i,4) = rtime           ! observation time
c
           store_2(i,1) = ta(i)           ! temperature (deg f)
           store_2(i,2) = tda(i)          ! dew point (deg f)
           store_2(i,3) = badflag         ! Relative Humidity
c
           if(ua(i).eq.badflag .or. abs(ua(i)).gt.200. .or. 
     &        va(i).eq.badflag .or. abs(va(i)).gt.200.) then
              spd = badflag
              dir = badflag
c 
           elseif(ua(i).eq.0.0 .and. va(i).eq.0.0) then
              spd = 0.0
              dir = 0.0                      !Undefined
c
           else
              spd = sqrt(ua(i)*ua(i) + va(i)*va(i) )   !speed
              dir = 57.2957795 * (atan2(ua(i),va(i))) + 180.   !dir

           endif
c
           store_3(i,1) = dir               ! wind dir (deg)
           store_3(i,2) = spd               ! wind speed (kt)
           store_3(i,3) = badflag           ! wind gust dir (deg)
           store_3(i,4) = badflag           ! wind gust speed (kt)
c
           store_4(i,1) = alta(i)           ! altimeter setting (mb)
           store_4(i,2) = badflag           ! station pressure (mb)
           store_4(i,3) = pmsla(i)          ! MSL pressure (mb)
           store_4(i,4) = float(delpch(ii)) ! 3-h press change character
           store_4(i,5) = delp(ii)          ! 3-h press change (mb)
c
           store_5(i,1) = vis(ii)           ! visibility (miles)
           store_5(i,2) = solar(ii)         ! solar radiation 
           store_5(i,3) = sfct(ii)          ! soil/water temperature
           store_5(i,4) = sfcm(ii)          ! soil moisture
c     
           store_6(i,1) = pcp1(ii)          ! 1-h precipitation
           store_6(i,2) = pcp3(ii)          ! 3-h precipitation
           store_6(i,3) = pcp6(ii)          ! 6-h precipitation
           store_6(i,4) = pcp24(ii)         ! 24-h precipitation
           store_6(i,5) = snow(ii)          ! snow cover
c
           kkk = kkk_s(ii)
           store_7(i,1) = float(kkk)        ! number of cloud layers
           store_7(i,2) = max24t(ii)        ! 24-h max temperature
           store_7(i,3) = min24t(ii)        ! 24-h min temperature
c
           store_2ea(i,1) = t_ea(ii)
           store_2ea(i,2) = td_ea(ii)
           store_2ea(i,3) = rh_ea(ii)
           store_3ea(i,1) = dd_ea(ii)
           store_3ea(i,2) = ff_ea(ii)
           store_4ea(i,1) = p_ea(ii)
           store_4ea(i,2) = alt_ea(ii)
           store_5ea(i,1) = vis_ea(ii)
           store_5ea(i,2) = solar_ea(ii)
           store_5ea(i,3) = sfct_ea(ii)
           store_5ea(i,4) = sfcm_ea(ii)
           store_6ea(i,1) = pcp_ea(ii)
           store_6ea(i,2) = snow_ea(ii)
c
           wmoid_out(i) = wmoid(ii)
           wx_out(i) = wx(ii)
           if(kkk .gt. 0) then
              do k=1,kkk
                 store_chtout(i,k) = store_cldht(ii,k)
                 store_camtout(i,k) = store_cldamt(ii,k)
              enddo !k
           endif
c
        enddo !i
c
        n_obs_g = 0
        n_obs_b = maxstaa
c     
        call get_directory('lso',outfile,len)
        outfile = outfile(1:len) // filename // '.lso_qc'
        call s_len(outfile, len)
        call write_surface_obs(atime_cur,outfile(1:len),n_obs_g,
     &    n_obs_b,wmoid_out,stations_out,provider,wx_out,reptype,
     &    autostntype,
     &    store_1,store_2,store_3,store_4,store_5,store_6,store_7,
     &    store_2ea,store_3ea,store_4ea,store_5ea,store_6ea,
     &    store_camtout,store_chtout,m,jstatus)
c Issue LAPS QC SUMMARY

           call writeqc(qcstat,qcstatd,qcstau,qcstav,qcstapm,
     &              qcstal,nvar,maxstaa,m,
     &          ncycles,qmonster,it+1) 
        call get_directory('log', dir_qcr, len)
        outfile = dir_qcr(1:len) //'qc_rpt.log'
        open(17,file=outfile(1:len+10),
     &         form='formatted',status='unknown')
       ncy=ncycles
       if(it+1.lt.ncycles) ncy=it+1
       write(17,3001) atime,ncy,it+1,maxstaa
 3001 format(1x,'LAPS QUALITY CONTROL REPORT FOR ',a24,
     &   'FOR THE LAST ', i3,' DATA CYCLES'/1x,
     &   ' QC ROUTINE HAS RUN FOR ',i10, 'CYCLES SINCE LAST', 
     &     ' RESTART'/1x,' THERE ARE ',i5,' STATIONS IN THE LIST')
       write(17,3055) 'KEY:N=number,E=Error,F=Failures,G=gross,'
       write(17,3055) 'S=self trend check, P=NWP Model check,  ' 
       write(17,3055) 'B=Buddy check, K=Kalman check, E=error  '
       write(17,3055) 'A=Average, %=percent, M=missing, O=Obs  ' 
 3055 format(1x,a40)
      do i=1,maxstaa
       write(17,3002) stna(i),providera(i),reptypea(i),lat(i),lon(i),
     &              elev(i)
 3002 format(1x,a5,1x,a11,1x,a7,1x,f5.2,1x,f7.2,f7.0)
       do n=1,nvar   
         nkf=0
         nnf=0
         nbf=0
         nsf=0
         ngf=0
         nms=0
         if(n.eq.1) then
          fname1='TEMPRTURE'
          b(1)=wot(i)
          b(2)=wbt(i)
          b(3)=wmt(i)
          b(4)=wkt(i)
         endif
         if(n.eq.2) then
          fname1='DEWPOINT '
          b(1)=wotd(i)
          b(2)=wbtd(i)
          b(3)=wmtd(i)
          b(4)=wktd(i)
         endif
         if(n.eq.3) then
          fname1='U-WIND   '
          b(1)=wou(i)
          b(2)=wbu(i)
          b(3)=wmu(i)
          b(4)=wku(i)
         endif
         if(n.eq.4) then
          fname1='V-WIND   '
          b(1)=wov(i)
          b(2)=wbv(i)
          b(3)=wmv(i)
          b(4)=wkv(i)
         endif
         if(n.eq.5) then
          fname1='MSL PRESS'
          b(1)=wop(i)
          b(2)=wbp(i)
          b(3)=wmp(i)
          b(4)=wkp(i)
         endif
         if(n.eq.6) then
          fname1='ALSTG    '
          b(1)=woa(i)
          b(2)=wba(i)
          b(3)=wma(i)
          b(4)=wka(i)
         endif
         if(n.ge.7) go to 8000
         do l=1,ncy 
c         break out qc information
           if(qmonster(i,l,n)-100000.ge.0) then
            nkf=nkf+1
            qmonster(i,l,n)=qmonster(i,l,n)-100000
           endif
           if(qmonster(i,l,n)-10000.ge.0) then
            nnf=nnf+1
            qmonster(i,l,n)=qmonster(i,l,n)-10000
           endif
           if(qmonster(i,l,n)-1000.ge.0) then
            nbf=nbf+1
            qmonster(i,l,n)=qmonster(i,l,n)-1000
           endif
           if(qmonster(i,l,n)-100.ge.0) then
            nsf=nsf+1
            qmonster(i,l,n)=qmonster(i,l,n)-100
           endif
           if(qmonster(i,l,n)-10.ge.0) then
            ngf=ngf+1
            qmonster(i,l,n)=qmonster(i,l,n)-10
           endif
           if(qmonster(i,l,n)-1.eq.0) then
            nms=nms+1
           endif
         enddo! on l
         write(17,3003) fname1,ncy,nms,int(float(nms)/float(ncy)*100),
     &                      ngf,int(float(ngf)/float(ncy)*100),
     &  nsf,int(float(nsf)/float(ncy)*100),b(1), 
     &  nbf,int(float(nbf)/float(ncy)*100),b(2),
     &  nnf,int(float(nnf)/float(ncy)*100),b(3),
     &  nkf,int(float(nkf)/float(ncy)*100),b(4) 
 3003  format(3x,a9,' NO ',i3,' NM ',i3,' %M ',i3,' GF ',i3,' %GF ',i3,
     & ' SF ',i3,' %SF ',i3,' ASE ',f5.1/3x,' BF ',i3,' %BF ',i3,
     & ' ABE ',f5.1,' PF ',i3,' %PF ',i3,' APE ',f5.1,' KF ',i3,
     & ' %KF ',i3,' AKE ',f5.1)                             
        enddo! on n
 8000  enddo!on i
       print*,'QC Report written to qc_rpt.log'
       close(17)
            
c
c...  That's it.
c     
        return
        end
      

