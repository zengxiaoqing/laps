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
        Real*4   lat(*),lon(*),elev(*),t(*),td(*),dd(*),ff(*),ddg(*)
        real*4   ffg(*),pstn(*),pmsl(*),alt(*),store_hgt(maxsta,5)
        real*4   ceil(*),lowcld(*),cover(*),vis(*),rad(*)
c
        Integer*4   obstime(*),kloud(*),idp3(*)
c
        Character   infile*70,atime*24,stn(*)*3,obstype(*)*8,wx(*)*8
        character   store_emv(maxsta,5)*1,store_amt(maxsta,5)*4
c
c.....  Start here.  Set the status to nothing, zero out the cloud storage.
c
        istatus = 0
        do j=1,5
        do i=1,maxsta
          store_emv(i,j) = ' '
          store_amt(i,j) = '    '
          store_hgt(i,j) = badflag
        enddo !i
        enddo !j
c
c.....  Open the file.  Check for a 'file not found' or other problem.
c
        

        open(1,iostat=ios,file=infile,status='old',access='sequential',
     &       form='formatted')
        if(ios .ne. 0) then     ! error during read
          istatus = -1
          write(6,650) infile
650       format(' +++ ERROR opening: ',a70,' +++')
          write(6,651) ios
651       format('     IOS code = ',i5)
          return
        endif
c
c.....  File open...first read the header.
c
        read(1,900) atime,              ! data time
     &               n_meso_g,          ! # of mesonet stations
     &               n_meso_pos,        ! total # mesonet stations possible
     &               n_sao_g,           ! # of saos in the laps grid
     &               n_sao_pos_g,       ! total # of saos possible in laps grid
     &               n_sao_b,           ! # of saos in the box
     &               n_sao_pos_b,       ! total # of saos possible in the box
     &               n_obs_g,           ! # of obs in the laps grid
     &               n_obs_pos_g,       ! total # of obs psbl in the laps grid
     &               n_obs_b,           ! # of obs in the box
     &               n_obs_pos_b        ! total # of obs possible in the box
900     format(1x,a24,10(1x,i4))
c
c.....  Error trapping for too many stations for array size.
c
        if(n_obs_b .gt. maxsta) then
          print 990, maxsta,n_obs_b,atime
990       format(' +++ ERROR in READ_SURFACE_OBS: maxstns = ',i8,/,
     & ' but there are ',i8,' stations in the ',a24,' obs file.',/)
          print *,'    Increase the value of "maxstns" and try again.'
          istatus = -2
          return
        endif
c
c.....  Now read the station data.
c
        do k=1,n_obs_b
          read(1,901) stn(k),lat(k),lon(k),elev(k),obstype(k),
     & obstime(k),wx(k)
901       format(1x,a3,2x,f6.2,1x,f7.2,1x,f5.0,1x,a8,1x,i4,1x,a8)
c
          read(1,903) t(k),td(k),dd(k),ff(k),ddg(k),ffg(k),pstn(k),
     &                pmsl(k),alt(k)
903       format(4x,2(f6.1,1x),4(f5.0,1x),3(f6.1,1x))
c
          read(1,905) kloud(k),ceil(k),lowcld(k),cover(k),vis(k),rad(k),
     & idp3(k)
905       format(4x,i2,2(1x,f7.1),1x,f5.1,1x,f7.3,1x,f6.1,1x,i4)
c
c.....  Read the cloud data if we have any.
c
          if(kloud(k) .gt. 0) then
            do ii=1,kloud(k)
              read(1,907) store_emv(k,ii),store_amt(k,ii),store_hgt(k,ii
     1)
907           format(5x,a1,1x,a4,1x,f7.1)
            enddo !ii
          endif
c
        enddo !k
c
        rewind(1)
        close(1)
c
c..... End of data gathering. Let's go home...
c
        istatus = 1             ! everything's ok...
        print *, ' Normal completion of READ_SURFACE_OBS'
c
        return
c
c
        entry read_surface_header(infile,atime,n_meso_g,n_meso_pos,
     &   n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,n_obs_pos_g,
     &   n_obs_b,n_obs_pos_b,istatus)
c
c.....  Entry to read and pass back just the header info from the lso file.
c
c.....  Open the file.  Check for a 'file not found' or other problem.
c
        istatus = 0
        open(1,iostat=ios,file=infile,status='old',access='sequential',
     &       form='formatted')
        if(ios .ne. 0) then     ! error during read
          istatus = -1
          write(6,650) infile
          write(6,651) ios
          return
        endif
c
c.....  File open...first read the header.
c
        read(1,900) atime,              ! data time
     &               n_meso_g,          ! # of mesonet stations
     &               n_meso_pos,        ! total # mesonet stations possible
     &               n_sao_g,           ! # of saos in the laps grid
     &               n_sao_pos_g,       ! total # of saos possible in laps grid
     &               n_sao_b,           ! # of saos in the box
     &               n_sao_pos_b,       ! total # of saos possible in the box
     &               n_obs_g,           ! # of obs in the laps grid
     &               n_obs_pos_g,       ! total # of obs psbl in the laps grid
     &               n_obs_b,           ! # of obs in the box
     &               n_obs_pos_b        ! total # of obs possible in the box
c
c.....  Rewind and close the file so we can call this again in the same program.
c
        rewind(1)
        close(1)
        istatus = 1
c
        return
        end
