cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   
cdis
c
c
        subroutine read_surface_sa(infile,maxsta,atime,n_obs_g,
     &   n_obs_b,stn,reptype,atype,lat,lon,elev,wx,
     &   t,td,dd,ff,ddg,ffg,pstn,pmsl,alt,kloud,ceil,lowcld,cover,rad,
     &   sfct,idp3,store_emv,store_amt,store_hgt,vis,obstime,istatus)       
c
cdoc    This routine calls 'read_surface_data' and is used primarily to read
cdoc    in cloud info along the lines of the arrays in the "old" LSO format.
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
c               S. Albers        99  Rename to 'read_surface_sa'
c
c       Input/Output:
c
c        Variable        Var type   I/O   Description
c       ----------      ---------- ----- -------------
c        infile            A*70      I    Directory where LSO file is.
c        maxsta             I        I    Max Number of stations allowed.
c        atime             A*24      O    Data time in dd-mmm-yyyy hh:mm
c        n_obs_g            I        O    Number of obs in the laps grid
c        n_obs_b            I        O    Number of obs in the box
c        stn               A*3 A     O    Station names (array)
c        lat                RA       O    Station latitude (deg)
c        lon                RA       O    Station longitude (deg)
c        elev               RA       O    Station elevation (m)
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
	character filetime*9, infile*256 
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
        Character   atime*24,stn(maxsta)*3     
        character   store_emv(maxsta,5)*1, wx(maxsta)*8
c
c
c.....  Start here.  Set the status to nothing, zero out the cloud storage
c.....  and character arrays.
c
        call get_sfc_badflag(badflag,istatus)

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
        if(jstatus .ne. 1 .and. jstatus .ne. -1) then
           print *,' ERROR: No valid LSO file found for ', filename       
           istatus = -1
           return
        endif
c
c.....  Shuffle data for the differences between old and new formats.
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
              call right_justify(stations(i))
              stn(i)(1:3) = stations(i)(18:20)

           endif
c
        enddo !i
c
c..... End of data gathering. Let's go home...
c
        istatus = 1             ! everything's ok...
        print *, ' Normal completion of new READ_SURFACE_SA'
c
        return
        end

