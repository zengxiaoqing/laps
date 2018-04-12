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
cdoc    Routine to read the newer LAPS LSO surface data file with the expanded 
cdoc    variable list. The data is passed back to the calling routine in 1-d 
cdoc    arrays.
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
	real lat(maxsta), lon(maxsta), elev(maxsta)
	real t(maxsta), t_ea(maxsta), max24t(maxsta), min24t(maxsta)
	real td(maxsta), td_ea(maxsta), rh(maxsta), rh_ea(maxsta)
	real dd(maxsta), ddg(maxsta), dd_ea(maxsta)
	real ff(maxsta), ffg(maxsta), ff_ea(maxsta)
	real alt(maxsta), alt_ea(maxsta), delp(maxsta)
	real stnp(maxsta), mslp(maxsta), p_ea(maxsta)
	real vis(maxsta), vis_ea(maxsta)
	real solar(maxsta), solar_ea(maxsta)
	real sfct(maxsta), sfct_ea(maxsta)
	real sfcm(maxsta), sfcm_ea(maxsta)
	real pcp1(maxsta), pcp3(maxsta), pcp6(maxsta), pcp24(maxsta)
	real snow(maxsta), snow_ea(maxsta), pcp_ea(maxsta)
	real store_cldht(maxsta,5)
c
	integer i4time, wmoid(maxsta), jstatus
	integer time(maxsta), delpch(maxsta), kkk_s(maxsta)
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
	open(11,file=infile,status='old',form='formatted',err=999)
c
c.....	Read the header.
c
	read(11,900) btime,		! time
     &               n_obs_g,		! # of obs in the laps grid
     &               n_obs_b		! # of obs in the box
c
c.....	Read the station data.
c
        if(n_obs_b.eq.0) then
         print*, 'no obs in file ',infile,' returning'
         jstatus=-1
         return
        endif
	do k=1,n_obs_b
c
	   read(11,901)  stations(k),            !station id
     &                   wmoid(k),               !WMO id number
     &                   provider(k),            !data provider
     &                   lat(k), lon(k), elev(k),!lat/lon (deg), elev (meters MSL)
     &                   time(k)		 !obs time
c
	  read(11,903)   reptype(k),             !station report type
     &                   autostntype(k),         !station type (manual/auto)
     &                   wx(k)                   !present weather
c
	  read(11,905)   t(k), t_ea(k),          !temp, temp expected accuracy (Deg F)
     &                   td(k), td_ea(k),        !dew point, dew point exp. accuracy (Deg F)
     &                   rh(k), rh_ea(k)         !Rel hum, rh expected accuracy (%)
c
	  read(11,907)   dd(k), ff(k),           !wind dir (deg), wind speed (knots)
     &                   ddg(k), ffg(k),         !wind gust dir (deg), wind gust speed (knots)
     &                   dd_ea(k), ff_ea(k)      !dir expected accuracy (deg), spd exp accuracy (knots)
c
	  read(11,909)   alt(k),                 !altimeter (mb)
     &                   stnp(k),                !station pressure (mb)
     &                   mslp(k),                !MSL pressure (mb)
     &                   delpch(k),              !3-h press change character (FMH-1 Manual, Sec 12.7.2)
     &                   delp(k),                !3-h pressure change (mb)
     &                   p_ea(k), alt_ea(k)      !pressure exp accuracy, alt exp accuracy
c
	  read(11,911)   vis(k), vis_ea(k),      !visibility, vis exp accuracy (miles)
     &                   solar(k), solar_ea(k),  !solar, solar exp accuracy (watts/meter**2)
     &                   sfct(k), sfct_ea(k),    !soil/water temp, soil/water temp exp accuracy (Deg F)
     &                   sfcm(k), sfcm_ea(k)     !soil moist, soil moist temp exp accuracy (%)
c
	  read(11,913)   pcp1(k),                !1-h precipitation (inches)
     &                   pcp3(k),                !3-h precipitation (inches)
     &                   pcp6(k),                !6-h precipitation (inches)
     &                   pcp24(k),               !24-h precipitation 
     &                   snow(k),                !snow depth (inches)
     &                   pcp_ea(k), snow_ea(k)   !precip and snow exp accuracy (inches)
c
	  read(11,915)  kkk_s(k),                !num cld layers 
     &                  max24t(k),               !24-h max temperature (Deg F)
     &                  min24t(k)                !24-h min temperature (Deg F)
c
c.....	Read the cloud data if we have any.
c
	  if(kkk_s(k) .gt. 0) then
	    do ii=1,kkk_s(k)
  	      read(11,917) store_cldamt(k,ii), store_cldht(k,ii) ! layer cloud amount and height
                                                                 ! (meters MSL)
	    enddo !ii
	  endif
c
	enddo !k
c
!       endfile(11)
	close(11)	
c
c..... End of data reading.  Let's go home...
c
	jstatus = 1
	return
999     continue
        print *,'Could not open file: ',infile
	print *,' ++ ERROR opening LSO file in READ_SURFACE_DATA ++'
        jstatus = -1
	return
        include 'lso_formats.inc'
	end
c
c
	subroutine read_sfc_metadata(i4time,n_obs_g,n_obs_b,time
     1              ,wmoid,stations,provider,reptype,autostntype
     1              ,lat,lon,elev,maxsta,jstatus)
c
c*****************************************************************************
c
cdoc    Routine to read the LAPS LSO surface data file and return the station
cdoc    metadata.   The data is passed back to the calling routine in 1-d 
cdoc    arrays.
c
c	Changes:
c		P. Stamus  04-13-98  Original version (from write_surface_LS2).
c      
c               S. Albers  2000      Expanded for OSSE use
c
c*****************************************************************************
c
	real lat(maxsta), lon(maxsta), elev(maxsta)
c
	integer i4time, wmoid(maxsta), jstatus
        integer time(maxsta)                 ! integer representation of HHMM
c
	character filetime*9, infile*256, btime*24
	character stations(maxsta)*20, provider(maxsta)*11
        character reptype(maxsta)*6, autostntype(maxsta)*6
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
	open(11,file=infile,status='old',form='formatted',err=999)
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
     &                   time(k)                   !obs time
c
	  read(11,903)   reptype(k),               !station report type
     &                   autostntype(k)            !station type (manual/auto)
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
!       endfile(11)
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
	subroutine read_sfc_state(i4time,ext,btime,n_obs_g,n_obs_b,
     &    stations,provider,lat,lon,elev,t,td,rh,dd,ff,alt,stnp,mslp,
     &    maxsta,jstatus)
c
c*****************************************************************************
c
cdoc    Routine to read the LAPS LSO surface data file and return the state
cdoc    variables (wind, temp, dewpt, altimeter & pressure).   The data is 
cdoc    passed back to the calling routine in 1-d arrays.
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
	real lat(maxsta), lon(maxsta), elev(maxsta)
	real t(maxsta)
	real td(maxsta)
	real rh(maxsta)
	real dd(maxsta)
	real ff(maxsta)
	real alt(maxsta)
	real stnp(maxsta), mslp(maxsta)
c
	integer i4time, jstatus
c
	character filetime*9, infile*256, btime*24
	character stations(maxsta)*20, provider(maxsta)*11
	character dum*132
        character ext*(*)
c
	jstatus = 0
c
c.....	Get the file.
c
	call make_fnam_lp(i4time, filetime, istatus)
	call get_directory('lso', infile, len)
        call s_len(ext,len_ext)
	infile = infile(1:len) // filetime(1:9) // '.' // ext(1:len_ext)       
c
	open(11,file=infile,status='old',form='formatted',err=999)
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
     &                   idummy                    !obs time
c
	  read(11,919) dum   
c
	  read(11,905)   t(k), dummy,              !temp, temp expected accuracy
     &                   td(k), dummy,             !dew point, dew point exp. accuracy
     &                   rh(k), dummy              !Rel hum, rh expected accuracy
c
	  read(11,907)   dd(k), ff(k),             !wind dir, wind speed
     &                   dummy, dummy,             !wind gust dir, wind gust speed
     &                   dummy, dummy              !dir expected accuracy, spd exp accuracy
c
	  read(11,909)   alt(k),                   !altimeter
     &                   stnp(k),                  !station pressure
     &                   mslp(k),                  !MSL pressure
     &                   idummy,                   !3-h press change character
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
!       endfile(11)
	close(11)	
c
c..... End of data reading.  Let's go home...
c
	jstatus = 1
	return
999     continue
	print *,' ++ ERROR opening ',ext(1:len_ext),
     1          ' file in READ_SFC_STATE ++'
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
cdoc    Routine to read the LAPS LSO surface data file and return the
cdoc    temperature.   The data is passed back to the calling routine in 
cdoc    1-d arrays.
c
c	Changes:
c		P. Stamus  04-13-98  Original version (from write_surface_LS2).
c
c*****************************************************************************
c
	real lat(maxsta), lon(maxsta), elev(maxsta)
	real t(maxsta)

c
	integer i4time, jstatus
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
	open(11,file=infile,status='old',form='formatted',err=999)
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
     &                   idummy                    !obs time
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
!       endfile(11)
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
        subroutine read_sfc_wind(i4time,ext,n_obs_g,n_obs_b,obstime
     1                          ,stations,provider,lat,lon,elev
     1                          ,dd,ff,dd_ea,ff_ea,maxsta,jstatus)
c
c*****************************************************************************
c
cdoc    Routine to read the LAPS LSO surface data file and return wind data. 
cdoc    The data is passed back to the calling routine in 1-d arrays.
c
c	Changes:
c		P. Stamus  04-13-98  Original version (from write_surface_LS2).
c
c*****************************************************************************
c
	real lat(maxsta), lon(maxsta), elev(maxsta)
	real dd(maxsta)
	real ff(maxsta)
	real dd_ea(maxsta)
	real ff_ea(maxsta)
c
	integer i4time, jstatus
        integer obstime(maxsta)
c
	character filetime*9, infile*256, btime*24
	character stations(maxsta)*20, provider(maxsta)*11
	character dum*132
        character*(*) ext
c
	jstatus = 0
c
c.....	Get the file.
c
	call make_fnam_lp(i4time, filetime, istatus)
	call get_directory('lso', infile, len)
	infile = infile(1:len) // filetime(1:9) // '.' // ext
c
	open(11,file=infile,status='old',form='formatted',err=999)
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
     &                   obstime(k)                !obs time
c
	  read(11,919) dum
c
	  read(11,919) dum
c
	  read(11,907)   dd(k), ff(k),             !wind dir, wind speed (kt)
     &                   dummy, dummy,             !wind gust dir, wind gust speed
     &                   dd_ea(k), ff_ea(k)        !dir expected accuracy, spd exp accuracy
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
!       endfile(11)
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
cdoc    Routine to read the LAPS LSO surface data file and return pressure
cdoc    data.   The data is passed back to the calling routine in 1-d arrays.
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
	real lat(maxsta), lon(maxsta), elev(maxsta)
	real alt(maxsta)
	real stnp(maxsta), mslp(maxsta)
c
	integer i4time, jstatus
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
	open(11,file=infile,status='old',form='formatted',err=999)
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
     &                   idummy                    !obs time
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
     &                   idummy,                   !3-h press change character
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
!       endfile(11)
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
cdoc    Routine to read the LAPS LSO surface data file and return precipitation
cdoc    data.   The data is passed back to the calling routine in 1-d arrays.
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
	real lat(maxsta), lon(maxsta), elev(maxsta)
	real pcp1(maxsta), pcp3(maxsta), pcp6(maxsta), pcp24(maxsta)
	real snow(maxsta)
c
	integer i4time, jstatus
c
	character filetime*9, infile*256, btime*24
	character stations(maxsta)*20, provider(maxsta)*11
	character dum*132
c
	jstatus = 0
        n_obs_g = 0
        n_obs_b = 0
c
c.....	Get the file.
c
	call make_fnam_lp(i4time, filetime, istatus)
	call get_directory('lso', infile, len)
	infile = infile(1:len) // filetime(1:9) // '.lso'
c
	open(11,file=infile,status='old',form='formatted',err=999)
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
     &                   idummy                    !obs time
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
!       endfile(11)
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
c
c
	subroutine read_surface_dataqc(i4time,btime,n_obs_g,n_obs_b,time,
     &    wmoid,stations,provider,wx,reptype,autostntype,lat,lon,elev,
     &    t,td,rh,dd,ff,ddg,ffg,alt,stnp,mslp,delpch,delp,vis,solar,
     &    sfct,sfcm,pcp1,pcp3,pcp6,pcp24,snow,kkk_s,max24t,min24t,t_ea,
     &    td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,sfct_ea,
     &    sfcm_ea,pcp_ea,snow_ea,store_cldamt,store_cldht,
     &    maxsta,jstatus)
c
c*****************************************************************************
c
cdoc    Routine to read the LAPS LSO_QC surface data file with the expanded 
cdoc    variable list. The data is passed back to the calling routine in 1-d 
cdoc    arrays.
c
c	Changes:
c		P. Stamus  12-17-98  Original version (from read_surface_data).
c
c*****************************************************************************
c
	real lat(maxsta), lon(maxsta), elev(maxsta)
	real t(maxsta), t_ea(maxsta), max24t(maxsta), min24t(maxsta)
	real td(maxsta), td_ea(maxsta), rh(maxsta), rh_ea(maxsta)
	real dd(maxsta), ddg(maxsta), dd_ea(maxsta)
	real ff(maxsta), ffg(maxsta), ff_ea(maxsta)
	real alt(maxsta), alt_ea(maxsta), delp(maxsta)
	real stnp(maxsta), mslp(maxsta), p_ea(maxsta)
	real vis(maxsta), vis_ea(maxsta)
	real solar(maxsta), solar_ea(maxsta)
	real sfct(maxsta), sfct_ea(maxsta)
	real sfcm(maxsta), sfcm_ea(maxsta)
	real pcp1(maxsta), pcp3(maxsta), pcp6(maxsta), pcp24(maxsta)
	real snow(maxsta), snow_ea(maxsta), pcp_ea(maxsta)
	real store_cldht(maxsta,5)
c
	integer i4time, wmoid(maxsta), jstatus
	integer time(maxsta), delpch(maxsta), kkk_s(maxsta)
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
	infile = infile(1:len) // filetime(1:9) // '.lso_qc'
c
	open(11,file=infile,status='old',form='formatted',err=999)
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
  	      read(11,917) store_cldamt(k,ii), store_cldht(k,ii) ! layer cloud amount and height
                                                                 ! (meters MSL)
	    enddo !ii
	  endif
c
	enddo !k
c
!       endfile(11)
	close(11)	
c
c..... End of data reading.  Let's go home...
c
	jstatus = 1
	return
999     continue
        print *,'Could not open file: ',infile
	print *,' ++ ERROR opening LSO_QC file in READ_SURFACE_DATAQC ++'
        jstatus = -1
	return
        include 'lso_formats.inc'
	end
