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
        subroutine read_surface_obs(infile,maxsta,atime,n_meso_g,
     &   n_meso_pos,n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,
     &   n_obs_pos_g,n_obs_b,n_obs_pos_b,stn,obstype,lat,lon,elev,wx,
     &   t,td,dd,ff,ddg,ffg,pstn,pmsl,alt,kloud,ceil,lowcld,cover,rad,
     &   idp3,store_emv,store_amt,store_hgt,vis,obstime,istatus)
c
c*******************************************************************************
c
c       Routine to read mesonet and SAO data for LAPS that has been written
c       into the .LSO file by the 'get_surface_obs' routine.
c
c       Changes:
c               P. Stamus  12-30-92  Original version.
c                          01-07-93  Add/change obs counters.
c                          01-08-93  Add read_header entry.
c                          12-21-95  Change format in prep of C*5 stn arrays.
c
c                      *** 09-03-98  MAJOR CHANGE...Now reading new format LSO
c                                        and passing arrays back with things as
c                                        close to old version as possible.
c
c       Input/Output:
c
c        Variable        Var type   I/O   Description
c       ----------      ---------- ----- -------------
c        infile            A*70      I    Directory where LSO file is.
c        maxsta             I        I    Max Number of stations allowed.
c        atime             A*24      O    Data time in dd-mmm-yyyy hh:mm
c        n_meso_g           I        O    Number of FSL mesonet stations
c        n_meso_pos         I        O    Total number mesonet stations psbl
c        n_sao_g            I        O    Number of SAOs in the laps grid
c        n_sao_pos_g        I        O    Total num. of SAOs psbl in laps grid
c        n_sao_b            I        O    Number of SAOs in the box
c        n_sao_pos_b        I        O    Total num of SAOs psbl in the box
c        n_obs_g            I        O    Number of obs in the laps grid
c        n_obs_pos_g        I        O    Total num of obs psbl in the laps grid
c        n_obs_b            I        O    Number of obs in the box
c        n_obs_pos_b        I        O    Total num of obs possible in the box
c        stn               A*3 A     O    Station names (array)
c        lat                RA       O    Station latitude (deg)
c        lon                RA       O    Station longitude (deg)
c        elev               RA       O    Station elevation (m)
c        obstype           A*8 A     O    Observation type (SA, SP, ASOS, etc)
c        obstime            IA       O    Time of observation (hhmm)
c        wx                A*8 A     O    Observed weather
c        t                  RA       O    Temperature (F)
c        td                 RA       O    Dewpoint (F)
c        dd                 RA       O    Wind direction (deg)
c        ff                 RA       O    Wind speed (kt)
c        ddg                RA       O    Gust wind direction (deg)
c        ffg                RA       O    Gust wind speed (kt)
c        pstn               RA       O    Station pressure (mb)
c        pmsl               RA       O    MSL pressure (mb)
c        alt                RA       O    Altimeter setting (mb)
c        kloud              IA       O    Number of cloud layers...max of 5.
c        ceil               RA       O    Ceiling height (m)
c        lowcld             RA       O    Height lowest cloud (m)
c        cover              RA       O    Cloud cover (tenths)
c        vis                RA       O    Visibility (miles)
c        rad                RA       O    Solar radiation.
c        idp3               IA       O    3-h coded pressure change (e.g.,608)
c        store_emv         A*1 A     O    Cloud descriptors: ea. layer, ea. stn
c        store_amt         A*4 A     O    Cloud layer coverage.
c        store_hgt          RA       O    Height of each cloud layer.
c        istatus            I        O    Status flag: 1 = normal
c                                                     -1 = file not found
c                                                     -2 = Arrays too small
c
c       User Notes:
c
c       1.  Arrays should be dimensioned 'maxsta' in the calling program,
c           with maxsta *at least* 120 (for CO domain).
c
c       2.  Pressures are stored as reported, except that altimeters are
c           converted to millibars.
c
c       3.  The 'kloud' variable tells whether there are clouds and how
c           many layers if there are:
c               a) kloud = 0    means   No cloud DATA (but NOT "no clouds").
c               b) kloud = 1    means   CLR or 1 cloud layer.  A height is
c                                       given for CLR which is the maximum valid
c                                       height of the observation (automatic
c                                       stations have limited valid heights).
c               c) kloud = 2-5  means   Two to five cloud layers.
c
c       4.  Thin obscured (-X) is a cloud layer and is given a 'badflag'
c           height, since it is not supposed to have a height (you're supposed
c           to be able to see other clouds and/or sky).
c
c*******************************************************************************
c
        real*4          badflag
        parameter       (badflag = -99.9)
c
c.....  Input arrays (for new format LSO)
c
	real*4 lat(maxsta), lon(maxsta), elev(maxsta)
	real*4 t(maxsta), t_ea(maxsta), max24t(maxsta), min24t(maxsta)
	real*4 td(maxsta), td_ea(maxsta), rh(maxsta), rh_ea(maxsta)
	real*4 dd(maxsta), ddg(maxsta), dd_ea(maxsta)
	real*4 ff(maxsta), ffg(maxsta), ff_ea(maxsta)
	real*4 alt(maxsta), alt_ea(maxsta), delp(maxsta)
	real*4 pstn(maxsta), pmsl(maxsta), p_ea(maxsta)
	real*4 vis(maxsta), vis_ea(maxsta)
	real*4 rad(maxsta), solar_ea(maxsta)
	real*4 sfct(maxsta), sfct_ea(maxsta)
	real*4 sfcm(maxsta), sfcm_ea(maxsta)
	real*4 pcp1(maxsta), pcp3(maxsta), pcp6(maxsta), pcp24(maxsta)
	real*4 snow(maxsta), snow_ea(maxsta), pcp_ea(maxsta)
	real*4 store_hgt(maxsta,5)
c
	integer*4 i4time, wmoid(maxsta), jstatus
	integer*4 time(maxsta), delpch(maxsta), kloud(maxsta)
c
	character filetime*9, infile*256, btime*24
	character stations(maxsta)*20, provider(maxsta)*11
        character reptype(maxsta)*6, atype(maxsta)*6
	character wx_in(maxsta)*25, store_amt(maxsta,5)*4
c
c.....  Output arrays (as old format LSO if different)
c
        real*4   ceil(maxsta),lowcld(maxsta),cover(maxsta)
c
        Integer*4   obstime(maxsta),idp3(maxsta)
c
        Character   atime*24,stn(maxsta)*3,obstype(maxsta)*8
        character   store_emv(maxsta,5)*1, wx(maxsta)*8
c
c
c.....  Start here.  Set the status to nothing, zero out the cloud storage
c.....  and character arrays.
c
        istatus = 0
        ibadflag = int(badflag)
c
        do j=1,5
        do i=1,maxsta
          store_emv(i,j) = ' '
          store_amt(i,j) = '    '
          store_hgt(i,j) = badflag
        enddo !i
        enddo !j
c
        do i=1,maxsta
           stn(i)(1:3) = '   '
           obstype(i)(1:8) = '        '
           wx(i)(1:8) = '        '
        enddo !i
c
c.....  Figure out the i4time and call the read routine.
c
        call s_len(infile, len)
        filetime(1:9) = infile(len-12:len-4)
        call i4time_fname_lp(filetime,i4time,status)
c
	call read_surface_data(i4time,atime,n_obs_g,n_obs_b,time,
     &    wmoid,stations,provider,wx_in,reptype,atype,lat,lon,
     &    elev,t,td,rh,dd,ff,ddg,ffg,alt,pstn,pmsl,delpch,delp,vis,rad,
     &    sfct,sfcm,pcp1,pcp3,pcp6,pcp24,snow,kloud,max24t,min24t,t_ea,
     &    td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,sfct_ea,
     &    sfcm_ea,pcp_ea,snow_ea,store_amt,store_hgt,
     &    maxsta,jstatus)
c
        if(jstatus .ne. 1) then
           print *,' ERROR.  No LSO file found for ', filename
           istatus = -1
           return
        endif
c
c.....  Shuffle data for the differences between old and new formats.
c.....  First, the header.
c
        n_meso_g = 0           ! # of mesonet stations
        n_meso_pos = 0         ! total # mesonet stations possible
        n_sao_g = 0            ! # of saos in the laps grid
        n_sao_pos_g = 0        ! total # of saos possible in laps grid
        n_sao_b = n_obs_b      ! # of saos in the box
        n_sao_pos_b = 0        ! total # of saos possible in the box
        n_obs_pos_g = 0        ! total # of obs psbl in the laps grid
        n_obs_pos_b = 0        ! total # of obs possible in the box
c
c.....  Now the station data.
c
        do i=1,n_obs_b
           idp3(i) = ibadflag
           cover(i) = badflag
           lowcld(i) = badflag
           ceil(i) = badflag
c
           wx(i)(1:8) = wx_in(i)(1:8)
c
           if(reptype(i)(1:4) .eq. 'LDAD') then
              if(provider(i)(1:4) .eq. 'CDOT') then
                 stn(i)(1:1) = stations(i)(1:1)
                 stn(i)(2:2) = stations(i)(3:3) ! skip the '-'
                 if(ichar(stations(i)(4:4)) .lt. 32) then
                    stn(i)(3:3) = ' '
                 else
                    stn(i)(3:3) = stations(i)(4:4)
                 endif
              else
                 stn(i)(1:3) = stations(i)(1:3)
              endif
           else
              stn(i)(1:3) = stations(i)(2:4)
           endif
c
           if(reptype(i)(1:4) .eq. 'LDAD') then
              obstype(i)(1:8) = provider(i)(1:8)
              do k=8,1,-1
                 if(ichar(obstype(i)(k:k)) .eq. 0) then
                    obstype(i)(k:k) = ' '
                 else
                    go to 150
                 endif
              enddo !k
 150          continue
           else
              obstype(i)(1:5) = reptype(i)(1:5)
              if(atype(i)(1:1) .eq. 'A') obstype(i)(8:8) = 'A'
              if(atype(i)(3:3) .eq. '1') obstype(i)(7:7) = '1'
              if(atype(i)(3:3) .eq. '2') obstype(i)(7:7) = '2'
           endif
c
           do j=1,5
              if(store_amt(i,j)(2:3) .eq. 'VV') 
     &                               store_amt(i,j)(1:4) = ' X  '
           enddo !j
c
        enddo !i
c
c..... End of data gathering. Let's go home...
c
        istatus = 1             ! everything's ok...
        print *, ' Normal completion of new READ_SURFACE_OBS'
c
        return
        end
c
c
	subroutine read_surface_data(i4time,btime,n_obs_g,n_obs_b,time,
     &    wmoid,stations,provider,wx,reptype,autostntype,lat,lon,elev,
     &    t,td,rh,dd,ff,ddg,ffg,alt,stnp,mslp,delpch,delp,vis,solar,
     &    sfct,sfcm,pcp1,pcp3,pcp6,pcp24,snow,kkk_s,max24t,min24t,t_ea,
     &    td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,sfct_ea,
     &    sfcm_ea,pcp_ea,snow_ea,store_cldamt,store_cldht,
     &    maxsta,jstatus)
c
c*****************************************************************************
c
c	Routine to read the LAPS LSO surface data file with the expanded 
c       format.   The data is passed back to the calling routine in 1-d arrays.
c
c	Changes:
c		P. Stamus  04-13-98  Original version (from write_surface_LS2).
c                          05-01-98  Added soil moisture variables.
c               J. Edwards 09-16-98  moved all format definitions to 
c                                    src/include/lso_formats.inc
c                                    changed 909 definition to allow for 
c                                    missing data
c
c*****************************************************************************
c
	real*4 lat(maxsta), lon(maxsta), elev(maxsta)
	real*4 t(maxsta), t_ea(maxsta), max24t(maxsta), min24t(maxsta)
	real*4 td(maxsta), td_ea(maxsta), rh(maxsta), rh_ea(maxsta)
	real*4 dd(maxsta), ddg(maxsta), dd_ea(maxsta)
	real*4 ff(maxsta), ffg(maxsta), ff_ea(maxsta)
	real*4 alt(maxsta), alt_ea(maxsta), delp(maxsta)
	real*4 stnp(maxsta), mslp(maxsta), p_ea(maxsta)
	real*4 vis(maxsta), vis_ea(maxsta)
	real*4 solar(maxsta), solar_ea(maxsta)
	real*4 sfct(maxsta), sfct_ea(maxsta)
	real*4 sfcm(maxsta), sfcm_ea(maxsta)
	real*4 pcp1(maxsta), pcp3(maxsta), pcp6(maxsta), pcp24(maxsta)
	real*4 snow(maxsta), snow_ea(maxsta), pcp_ea(maxsta)
	real*4 store_cldht(maxsta,5)
c
	integer*4 i4time, wmoid(maxsta), jstatus
	integer*4 time(maxsta), delpch(maxsta), kkk_s(maxsta)
c
	character filetime*9, infile*256, btime*24
	character stations(maxsta)*20, provider(maxsta)*11
        character reptype(maxsta)*6, autostntype(maxsta)*6
	character wx(maxsta)*25, store_cldamt(maxsta,5)*4
c
c
c.....  Blank out the character arrays.
c
	jstatus = 0
	do i=1,maxsta
	   stations(i) = '                    '
	   provider(i) = '           '
	   reptype(i)  = '      '
	   autostntype(i)  = '      '
	   wx(i) = '                         '
	   do j=1,5
	      store_cldamt(i,j) = '    '
	   enddo !j
	enddo !i
c
c.....	Get the file.
c
	call make_fnam_lp(i4time, filetime, istatus)
	call get_directory('lso', infile, len)
	infile = infile(1:len) // filetime(1:9) // '.lso'
c
	open(11,file=infile,status='old',err=999)
c
c.....	Read the header.
c
	read(11,900) btime,		! time
     &               n_obs_g,		! # of obs in the laps grid
     &               n_obs_b		! # of obs in the box
c
c.....	Read the station data.
c
	do k=1,n_obs_b
c
	   read(11,901)  stations(k),              !station id
     &                   wmoid(k),                 !WMO id number
     &                   provider(k),              !data provider
     &                   lat(k), lon(k), elev(k),  !lat, lon, elev
     &                   time(k)		   !obs time
c
	  read(11,903)   reptype(k),               !station report type
     &                   autostntype(k),           !station type (manual/auto)
     &                   wx(k)                     !present weather
c
	  read(11,905)   t(k), t_ea(k),            !temp, temp expected accuracy
     &                   td(k), td_ea(k),          !dew point, dew point exp. accuracy
     &                   rh(k), rh_ea(k)           !Rel hum, rh expected accuracy
c
	  read(11,907)   dd(k), ff(k),             !wind dir, wind speed
     &                   ddg(k), ffg(k),           !wind gust dir, wind gust speed
     &                   dd_ea(k), ff_ea(k)        !dir expected accuracy, spd exp accuracy
c
	  read(11,909)   alt(k),                   !altimeter
     &                   stnp(k),                  !station pressure
     &                   mslp(k),                  !MSL pressure
     &                   delpch(k),                !3-h press change character
     &                   delp(k),                  !3-h pressure change
     &                   p_ea(k), alt_ea(k)        !pressure exp accuracy, alt exp accuracy
c
	  read(11,911)   vis(k), vis_ea(k),        !visibility, vis exp accuracy
     &                   solar(k), solar_ea(k),    !solar, solar exp accuracy
     &                   sfct(k), sfct_ea(k),      !soil/water temp, soil/water temp exp accuracy
     &                   sfcm(k), sfcm_ea(k)       !soil moist, soil moist temp exp accuracy
c
	  read(11,913)   pcp1(k),                  !1-h precipitation
     &                   pcp3(k),                  !3-h precipitation
     &                   pcp6(k),                  !6-h precipitation
     &                   pcp24(k),                 !24-h precipitation
     &                   snow(k),                  !snow depth
     &                   pcp_ea(k), snow_ea(k)     !precip and snow exp accuracy
c
	  read(11,915)  kkk_s(k),                  !num cld layers 
     &                  max24t(k),                 !24-h max temperature
     &                  min24t(k)                  !24-h min temperature
c
c.....	Read the cloud data if we have any.
c
	  if(kkk_s(k) .gt. 0) then
	    do ii=1,kkk_s(k)
  	      read(11,917) store_cldamt(k,ii), store_cldht(k,ii)   !layer cloud amount and height
	    enddo !ii
	  endif
c
	enddo !k
c
	endfile(11)
	close(11)	
c
c..... End of data reading.  Let's go home...
c
	jstatus = 1
	return
999     continue
	print *,' ++ ERROR opening LSO file in READ_SURFACE_DATA ++'
        jstatus = -1
	return
        include 'lso_formats.inc'
	end
c
c
	subroutine read_sfc_metadata(i4time,n_obs_g,n_obs_b,
     &    wmoid,stations,provider,lat,lon,elev,maxsta,jstatus)
c
c*****************************************************************************
c
c	Routine to read the LAPS LSO surface data file and return the station
c       metadata.   The data is passed back to the calling routine in 1-d 
c       arrays.
c
c	Changes:
c		P. Stamus  04-13-98  Original version (from write_surface_LS2).
c
c*****************************************************************************
c
	real*4 lat(maxsta), lon(maxsta), elev(maxsta)
c
	integer*4 i4time, wmoid(maxsta), jstatus
c
	character filetime*9, infile*256, btime*24
	character stations(maxsta)*20, provider(maxsta)*11
	character dum*132
c
	jstatus = 0
c
c.....	Get the file.
c
	call make_fnam_lp(i4time, filetime, istatus)
	call get_directory('lso', infile, len)
	infile = infile(1:len) // filetime(1:9) // '.lso'
c
	open(11,file=infile,status='old',err=999)
c
c.....	Read the header.
c
	read(11,900) btime,		! time
     &               n_obs_g,		! # of obs in the laps grid
     &               n_obs_b		! # of obs in the box
c
c.....	Read the station data.
c
	do k=1,n_obs_b
c
	   read(11,901)  stations(k),              !station id
     &                   wmoid(k),                 !WMO id number
     &                   provider(k),              !data provider
     &                   lat(k), lon(k), elev(k),  !lat, lon, elev
     &                   idummy        		   !obs time
c
	  read(11,919) dum                       
c
	  read(11,919) dum
c
	  read(11,919) dum
c
	  read(11,919) dum
c
	  read(11,919) dum
c
	  read(11,919) dum
c
	  read(11,915)  kkk_s,                     !num cld layers 
     &                  dummy,                     !24-h max temperature
     &                  dummy                      !24-h min temperature
c
c.....	Read the cloud data if we have any.
c
	  if(kkk_s .gt. 0) then
	    do ii=1,kkk_s
	      read(11,919) dum
	    enddo !ii
	  endif
c
	enddo !k
c
	endfile(11)
	close(11)	
c
c..... End of data reading.  Let's go home...
c
	jstatus = 1
	return
999     continue
	print *,' ++ ERROR opening LSO file in READ_SFC_METADATA ++'
        jstatus = -1
	return
        include 'lso_formats.inc'
	end
c
c
	subroutine read_sfc_state(i4time,btime,n_obs_g,n_obs_b,
     &    stations,provider,lat,lon,elev,t,td,dd,ff,alt,stnp,mslp,
     &    maxsta,jstatus)
c
c*****************************************************************************
c
c	Routine to read the LAPS LSO surface data file and return the state
c       variables (wind, temp, dewpt, altimeter & pressure).   The data is 
c       passed back to the calling routine in 1-d arrays.
c
c	Changes:
c		P. Stamus  04-13-98  Original version (from write_surface_LS2).
c               J. Edwards 09-16-98  moved all format definitions to 
c                                    src/include/lso_formats.inc
c                                    changed 909 definition to allow for 
c                                    missing data
c
c*****************************************************************************
c
	real*4 lat(maxsta), lon(maxsta), elev(maxsta)
	real*4 t(maxsta)
	real*4 td(maxsta)
	real*4 dd(maxsta)
	real*4 ff(maxsta)
	real*4 alt(maxsta)
	real*4 stnp(maxsta), mslp(maxsta)
c
	integer*4 i4time, jstatus
c
	character filetime*9, infile*256, btime*24
	character stations(maxsta)*20, provider(maxsta)*11
	character dum*132
c
	jstatus = 0
c
c.....	Get the file.
c
	call make_fnam_lp(i4time, filetime, istatus)
	call get_directory('lso', infile, len)
	infile = infile(1:len) // filetime(1:9) // '.lso'
c
	open(11,file=infile,status='old',err=999)
c
c.....	Read the header.
c
	read(11,900) btime,		! time
     &               n_obs_g,		! # of obs in the laps grid
     &               n_obs_b		! # of obs in the box
c
c.....	Read the station data.
c
	do k=1,n_obs_b
c
	   read(11,901)  stations(k),              !station id
     &                   idummy,                   !WMO id number
     &                   provider(k),              !data provider
     &                   lat(k), lon(k), elev(k),  !lat, lon, elev
     &                   dummy                     !obs time
c
	  read(11,919) dum   
c
	  read(11,905)   t(k), dummy,              !temp, temp expected accuracy
     &                   td(k), dummy,             !dew point, dew point exp. accuracy
     &                   dummy, dummy              !Rel hum, rh expected accuracy
c
	  read(11,907)   dd(k), ff(k),             !wind dir, wind speed
     &                   dummy, dummy,             !wind gust dir, wind gust speed
     &                   dummy, dummy              !dir expected accuracy, spd exp accuracy
c
	  read(11,909)   alt(k),                   !altimeter
     &                   stnp(k),                  !station pressure
     &                   mslp(k),                  !MSL pressure
     &                   dummy,                    !3-h press change character
     &                   dummy,                    !3-h pressure change
     &                   dummy, dummy              !pressure exp accuracy, alt exp accuracy
c
	  read(11,919) dum   
c
	  read(11,919) dum   
c
	  read(11,915)  kkk_s,                     !num cld layers 
     &                  dummy,                     !24-h max temperature
     &                  dummy                      !24-h min temperature
c
c.....	Read the cloud data if we have any.
c
	  if(kkk_s .gt. 0) then
	    do ii=1,kkk_s
  	      read(11,919) dum                     !layer cloud amount and height
	    enddo !ii
	  endif
c
	enddo !k
c
	endfile(11)
	close(11)	
c
c..... End of data reading.  Let's go home...
c
	jstatus = 1
	return
999     continue
	print *,' ++ ERROR opening LSO file in READ_SFC_STATE ++'
        jstatus = -1
	return
        include 'lso_formats.inc'
	end
c
c
	subroutine read_sfc_temp(i4time,btime,n_obs_g,n_obs_b,
     &    stations,provider,lat,lon,elev,t,maxsta,jstatus)
c
c*****************************************************************************
c
c	Routine to read the LAPS LSO surface data file and return the
c       temperature.   The data is passed back to the calling routine in 
c       1-d arrays.
c
c	Changes:
c		P. Stamus  04-13-98  Original version (from write_surface_LS2).
c
c*****************************************************************************
c
	real*4 lat(maxsta), lon(maxsta), elev(maxsta)
	real*4 t(maxsta)

c
	integer*4 i4time, jstatus
c
	character filetime*9, infile*256, btime*24
	character stations(maxsta)*20, provider(maxsta)*11
	character dum*132
c
	jstatus = 0
c
c.....	Get the file.
c
	call make_fnam_lp(i4time, filetime, istatus)
	call get_directory('lso', infile, len)
	infile = infile(1:len) // filetime(1:9) // '.lso'
c
	open(11,file=infile,status='old',err=999)
c
c.....	Read the header.
c
	read(11,900) btime,		! time
     &               n_obs_g,		! # of obs in the laps grid
     &               n_obs_b		! # of obs in the box
c
c.....	Read the station data.
c
	do k=1,n_obs_b
c
	   read(11,901)  stations(k),              !station id
     &                   idummy,                   !WMO id number
     &                   provider(k),              !data provider
     &                   lat(k), lon(k), elev(k),  !lat, lon, elev
     &                   dummy                     !obs time
c
	  read(11,919) dum
c
	  read(11,905)   t(k), dummy,              !temp, temp expected accuracy
     &                   dummy, dummy,             !dew point, dew point exp. accuracy
     &                   dummy, dummy              !Rel hum, rh expected accuracy
c
	  read(11,919) dum
c
	  read(11,919) dum
c
	  read(11,919) dum
c
	  read(11,919) dum
c
	  read(11,915)  kkk_s,                     !num cld layers 
     &                  dummy,                     !24-h max temperature
     &                  dummy                      !24-h min temperature
c
c.....	Read the cloud data if we have any.
c
	  if(kkk_s .gt. 0) then
	    do ii=1,kkk_s
  	      read(11,919) dum                     !layer cloud amount and height
	    enddo !ii
	  endif
c
	enddo !k
c
	endfile(11)
	close(11)	
c
c..... End of data reading.  Let's go home...
c
	jstatus = 1
	return
999     continue
	print *,' ++ ERROR opening LSO file in READ_SFC_TEMP ++'
        jstatus = -1
	return
        include 'lso_formats.inc'
	end
c
c
	subroutine read_sfc_wind(i4time,btime,n_obs_g,n_obs_b,
     &    stations,provider,lat,lon,elev,dd,ff,maxsta,jstatus)
c
c*****************************************************************************
c
c	Routine to read the LAPS LSO surface data file and return wind data. 
c       The data is passed back to the calling routine in 1-d arrays.
c
c	Changes:
c		P. Stamus  04-13-98  Original version (from write_surface_LS2).
c
c*****************************************************************************
c
	real*4 lat(maxsta), lon(maxsta), elev(maxsta)
	real*4 dd(maxsta)
	real*4 ff(maxsta)
c
	integer*4 i4time, jstatus
c
	character filetime*9, infile*256, btime*24
	character stations(maxsta)*20, provider(maxsta)*11
	character dum*132
c
	jstatus = 0
c
c.....	Get the file.
c
	call make_fnam_lp(i4time, filetime, istatus)
	call get_directory('lso', infile, len)
	infile = infile(1:len) // filetime(1:9) // '.lso'
c
	open(11,file=infile,status='old',err=999)
c
c.....	Read the header.
c
	read(11,900) btime,		! time
     &               n_obs_g,		! # of obs in the laps grid
     &               n_obs_b		! # of obs in the box
c
c.....	Read the station data.
c
	do k=1,n_obs_b
c
	   read(11,901)  stations(k),              !station id
     &                   idummy,                   !WMO id number
     &                   provider(k),              !data provider
     &                   lat(k), lon(k), elev(k),  !lat, lon, elev
     &                   dummy                     !obs time
c
	  read(11,919) dum
c
	  read(11,919) dum
c
	  read(11,907)   dd(k), ff(k),             !wind dir, wind speed
     &                   dummy, dummy,             !wind gust dir, wind gust speed
     &                   dummy, dummy              !dir expected accuracy, spd exp accuracy
c
	  read(11,919) dum
c
	  read(11,919) dum
c
	  read(11,919) dum
c
	  read(11,915)  kkk_s,                     !num cld layers 
     &                  dummy,                     !24-h max temperature
     &                  dummy                      !24-h min temperature
c
c.....	Read the cloud data if we have any.
c
	  if(kkk_s .gt. 0) then
	    do ii=1,kkk_s
  	      read(11,919) dum
	    enddo !ii
	  endif
c
	enddo !k
c
	endfile(11)
	close(11)	
c
c..... End of data reading.  Let's go home...
c
	jstatus = 1
	return
999     continue
	print *,' ++ ERROR opening LSO file in READ_SFC_WIND ++'
        jstatus = -1
	return
        include 'lso_formats.inc'
	end
c
c
	subroutine read_sfc_press(i4time,btime,n_obs_g,n_obs_b,
     &    stations,provider,lat,lon,elev,alt,stnp,mslp,maxsta,jstatus)
c
c*****************************************************************************
c
c	Routine to read the LAPS LSO surface data file and return pressure
c       data.   The data is passed back to the calling routine in 1-d arrays.
c
c	Changes:
c		P. Stamus  04-13-98  Original version (from write_surface_LS2).
c               J. Edwards 09-16-98  moved all format definitions to 
c                                    src/include/lso_formats.inc
c                                    changed 909 definition to allow for 
c                                    missing data
c
c*****************************************************************************
c
	real*4 lat(maxsta), lon(maxsta), elev(maxsta)
	real*4 alt(maxsta)
	real*4 stnp(maxsta), mslp(maxsta)
c
	integer*4 i4time, jstatus
c
	character filetime*9, infile*256, btime*24
	character stations(maxsta)*20, provider(maxsta)*11
	character dum*132
c
	jstatus = 0
c
c.....	Get the file.
c
	call make_fnam_lp(i4time, filetime, istatus)
	call get_directory('lso', infile, len)
	infile = infile(1:len) // filetime(1:9) // '.lso'
c
	open(11,file=infile,status='old',err=999)
c
c.....	Read the header.
c
	read(11,900) btime,		! time
     &               n_obs_g,		! # of obs in the laps grid
     &               n_obs_b		! # of obs in the box
c
c.....	Read the station data.
c
	do k=1,n_obs_b
c
	   read(11,901)  stations(k),              !station id
     &                   idummy,                   !WMO id number
     &                   provider(k),              !data provider
     &                   lat(k), lon(k), elev(k),  !lat, lon, elev
     &                   dummy                     !obs time
c
	  read(11,919) dum
c
	  read(11,919) dum
c
	  read(11,919) dum
c
	  read(11,909)   alt(k),                   !altimeter
     &                   stnp(k),                  !station pressure
     &                   mslp(k),                  !MSL pressure
     &                   dummy,                    !3-h press change character
     &                   dummy,                    !3-h pressure change
     &                   dummy, dummy              !pressure exp accuracy, alt exp accuracy
c
	  read(11,919) dum
c
	  read(11,919) dum
c
	  read(11,915)  kkk_s,                     !num cld layers 
     &                  dummy,                     !24-h max temperature
     &                  dummy                      !24-h min temperature
c
c.....	Read the cloud data if we have any.
c
	  if(kkk_s .gt. 0) then
	    do ii=1,kkk_s
  	      read(11,919) dum
	    enddo !ii
	  endif
c
	enddo !k
c
	endfile(11)
	close(11)	
c
c..... End of data reading.  Let's go home...
c
	jstatus = 1
	return
999     continue
	print *,' ++ ERROR opening LSO file in READ_SFC_PRESS ++'
        jstatus = -1
	return
        include 'lso_formats.inc'
	end
c
c
	subroutine read_sfc_precip(i4time,btime,n_obs_g,n_obs_b,
     &    stations,provider,lat,lon,elev,pcp1,pcp3,pcp6,pcp24,snow,
     &    maxsta,jstatus)
c
c*****************************************************************************
c
c	Routine to read the LAPS LSO surface data file and return precipitation
c       data.   The data is passed back to the calling routine in 1-d arrays.
c
c	Changes:
c		P. Stamus  04-13-98  Original version (from write_surface_LS2).
c               J. Edwards 09-16-98  moved all format definitions to 
c                                    src/include/lso_formats.inc
c                                    changed 909 definition to allow for 
c                                    missing data
c
c*****************************************************************************
c
	real*4 lat(maxsta), lon(maxsta), elev(maxsta)
	real*4 pcp1(maxsta), pcp3(maxsta), pcp6(maxsta), pcp24(maxsta)
	real*4 snow(maxsta)
c
	integer*4 i4time, jstatus
c
	character filetime*9, infile*256, btime*24
	character stations(maxsta)*20, provider(maxsta)*11
	character dum*132
c
	jstatus = 0
c
c.....	Get the file.
c
	call make_fnam_lp(i4time, filetime, istatus)
	call get_directory('lso', infile, len)
	infile = infile(1:len) // filetime(1:9) // '.lso'
c
	open(11,file=infile,status='old',err=999)
c
c.....	Read the header.
c
	read(11,900) btime,		! time
     &               n_obs_g,		! # of obs in the laps grid
     &               n_obs_b		! # of obs in the box
c
c.....	Read the station data.
c
	do k=1,n_obs_b
c
	   read(11,901)  stations(k),              !station id
     &                   idummy,                   !WMO id number
     &                   provider(k),              !data provider
     &                   lat(k), lon(k), elev(k),  !lat, lon, elev
     &                   dummy                     !obs time
c
	  read(11,919) dum
c
	  read(11,919) dum
c
	  read(11,919) dum
c
	  read(11,919) dum
c
	  read(11,919) dum
c
	  read(11,913)   pcp1(k),                  !1-h precipitation
     &                   pcp3(k),                  !3-h precipitation
     &                   pcp6(k),                  !6-h precipitation
     &                   pcp24(k),                 !24-h precipitation
     &                   snow(k),                  !snow depth
     &                   dummy, dummy              !precip and snow exp accuracy
c
	  read(11,915)  kkk_s,                     !num cld layers 
     &                  dummy,                     !24-h max temperature
     &                  dummy                      !24-h min temperature
c
c.....	Read the cloud data if we have any.
c
	  if(kkk_s .gt. 0) then
	    do ii=1,kkk_s
  	      read(11,919) dum
	    enddo !ii
	  endif
c
	enddo !k
c
	endfile(11)
	close(11)	
c
c..... End of data reading.  Let's go home...
c
	jstatus = 1
	return
 999    continue
	print *,' ++ ERROR opening LSO file in READ_SFC_PRECIP ++'
        jstatus 	= -1
	return
        include 'lso_formats.inc'
	end
