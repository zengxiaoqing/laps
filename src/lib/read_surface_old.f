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
        subroutine read_surface_old(infile,maxsta,atime,n_meso_g,
     &   n_meso_pos,n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,
     &   n_obs_pos_g,n_obs_b,n_obs_pos_b,stn,obstype,lat,lon,elev,wx,
     &   t,td,dd,ff,ddg,ffg,pstn,pmsl,alt,kloud,ceil,lowcld,cover,rad,
     &   idp3,store_emv,store_amt,store_hgt,vis,obstime,istatus)
c
c*******************************************************************************
c
cdoc    Routine to read surface data for LAPS that has been written into
cdoc    the LSO/LSOQC file by the 'write_surface_obs' routine. Return arguments
cdoc    are more in tune with an earlier lso variable list. Calls 
cdoc    'read_surface_data' and 'read_surface_dataqc'.
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
        real          badflag
        parameter       (badflag = -99.9)
c
c.....  Input arrays (for new format LSO)
c
	real lat(maxsta), lon(maxsta), elev(maxsta)
	real t(maxsta), t_ea(maxsta), max24t(maxsta), min24t(maxsta)
	real td(maxsta), td_ea(maxsta), rh(maxsta), rh_ea(maxsta)
	real dd(maxsta), ddg(maxsta), dd_ea(maxsta)
	real ff(maxsta), ffg(maxsta), ff_ea(maxsta)
	real alt(maxsta), alt_ea(maxsta), delp(maxsta)
	real pstn(maxsta), pmsl(maxsta), p_ea(maxsta)
	real vis(maxsta), vis_ea(maxsta)
	real rad(maxsta), solar_ea(maxsta)
	real sfct(maxsta), sfct_ea(maxsta)
	real sfcm(maxsta), sfcm_ea(maxsta)
	real pcp1(maxsta), pcp3(maxsta), pcp6(maxsta), pcp24(maxsta)
	real snow(maxsta), snow_ea(maxsta), pcp_ea(maxsta)
	real store_hgt(maxsta,5)
c
	integer i4time, wmoid(maxsta), jstatus
	integer time(maxsta), delpch(maxsta), kloud(maxsta)
c
	character filetime*9, infile*256 
	character stations(maxsta)*20, provider(maxsta)*11
        character reptype(maxsta)*6, atype(maxsta)*6
	character wx_in(maxsta)*25, store_amt(maxsta,5)*4
c
c.....  Output arrays (as old format LSO if different)
c
        real   ceil(maxsta),lowcld(maxsta),cover(maxsta)
c
        Integer   obstime(maxsta),idp3(maxsta)
c
        Character   atime*24,stn(maxsta)*3,obstype(maxsta)*8
        character   store_emv(maxsta,5)*1, wx(maxsta)*8
c
c
c.....  Start here.  Set the status to nothing, zero out the cloud storage
c.....  and character arrays.
c
        write(6,*)' Subroutine read_surface_old...'
 
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
c
        write(6,*)' infile = ',infile(1:len)

        if(infile(len-1:len) .eq. 'qc')then
           filetime(1:9) = infile(len-15:len-7)
           call i4time_fname_lp(filetime,i4time,status)
           write(6,*)' Calling read_surface_dataqc for ',filetime

	   call read_surface_dataqc(i4time,atime,n_obs_g,n_obs_b,time,
     &     wmoid,stations,provider,wx_in,reptype,atype,lat,lon,
     &     elev,t,td,rh,dd,ff,ddg,ffg,alt,pstn,pmsl,delpch,delp,vis,rad,
     &     sfct,sfcm,pcp1,pcp3,pcp6,pcp24,snow,kloud,max24t,min24t,t_ea,
     &     td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,sfct_ea,
     &     sfcm_ea,pcp_ea,snow_ea,store_amt,store_hgt,
     &     maxsta,jstatus)
c
           if(jstatus .ne. 1) then
              print *,' ERROR.  No LSOQC file found for ', filetime
              istatus = -1
              return
           endif

        else
           filetime(1:9) = infile(len-12:len-4)
           call i4time_fname_lp(filetime,i4time,status)
           write(6,*)' Calling read_surface_data for ',filetime

	   call read_surface_data(i4time,atime,n_obs_g,n_obs_b,time,
     &     wmoid,stations,provider,wx_in,reptype,atype,lat,lon,
     &     elev,t,td,rh,dd,ff,ddg,ffg,alt,pstn,pmsl,delpch,delp,vis,rad,
     &     sfct,sfcm,pcp1,pcp3,pcp6,pcp24,snow,kloud,max24t,min24t,t_ea,
     &     td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,sfct_ea,
     &     sfcm_ea,pcp_ea,snow_ea,store_amt,store_hgt,
     &     maxsta,jstatus)
c
           if(jstatus .ne. 1) then
              print *,' ERROR.  No LSO file found for ', filetime
              istatus = -1
              return
           endif

        endif ! lso_qc infile
c
c.....  Shuffle data for the differences between old and new formats.
c.....  First, the header.
c
        n_meso_g = 0           ! # of mesonet stations
        n_meso_pos = 0         ! total # mesonet stations possible
        n_sao_g = 0            ! # of saos in the laps grid
        n_sao_pos_g = 0        ! total # of saos possible in laps grid
        n_sao_b = n_obs_b      ! # of saos in the box
        n_sao_pos_b = n_obs_b  ! total # of saos possible in the box
        n_obs_pos_g = n_obs_g  ! total # of obs psbl in the laps grid
        n_obs_pos_b = n_obs_b  ! total # of obs possible in the box
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

           if(reptype(i)(1:4) .eq. 'LDAD') then
              wx(i) = 'UNKNOWN'
           endif
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

           else ! right justify the string
              call left_justify(stations(i))
              call s_len(stations(i),len_sta)
              len_sta = max(3,len_sta)
              stn(i)(1:3) = stations(i)(len_sta-2:len_sta)

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
        enddo !i
c
c..... End of data gathering. Let's go home...
c
        istatus = 1             ! everything's ok...
        print *, ' Normal completion of new READ_SURFACE_OLD'
c
        return
        end

        subroutine left_justify(string)

        character*(*) string

        len_string = len(string)
        call filter_string(string)

        index_first = 0
        do j = len_string,1,-1
            if(string(j:j) .ne. ' ')then
                index_first = j
            endif
        enddo
        
        if(index_first .gt. 1)then
            index_last_out = len_string-index_first+1
            string(1:index_last_out) 
     1                         = string(index_first:len_string)

            if(index_last_out .lt. len_string)then
                do j = index_last_out+1,len_string
                    string(j:j) = ' '
                enddo ! j
            endif

        endif

        return
        end

        subroutine right_justify(string)

        character*(*) string

        call left_justify(string)

        call s_len(string,len1)
        len2 = len(string)

        iarg = len2-len1+1
        string(iarg:len2) = string(1:len1)

        if(iarg-1 .ge. 1)then
            do i = 1,iarg-1
                string(i:i) = ' '
            enddo ! i
        endif

        return
        end        

