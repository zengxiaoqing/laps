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
        subroutine read_obs_i(infile,maxstns,atime,num_meso,
     &   num_saos,num_sfc,stations,lat_s,lon_s,elev_s,wx_s,cover_s,
     &   hgt_ceil,hgt_low,t_s,td_s,dd_s,ff_s,ddg_s,ffg_s,pr_s,sr_s,
     &   istatus)
c
c*******************************************************************************
c
c       Routine to read SAO and Mesonet surface obs written by the LAPS
c       'lapsdata' program.
c
c       Changes:
c               P.A. Stamus     11-29-88        Original version.
c                               02-01-90        Version for interactive MDAT.
c                               02-14-91        Add solar radiation.
c
c       Input/Output:
c
c          Variable      Var type    I/O    Description
c         ----------    ----------  -----  -------------
c          infile          A*70       I     Directory where input data is.
c          maxstns          I         I     Max number of stations allowed
c          atime           A*24       O     Data time: dd-mmm-yyyy hh:mm
c          num_meso         I         O     Number of mesonet stations in file
c          num_saos         I         O     Number of SAO stations in file
c          num_sfc          I         O     Total number of surface obs.
c          stations        A*3 A      O     Array of the station names
c          lat_s            RA        O     Latitude of the stations
c          lon_s            RA        O     Longitude of the stations
c          elev_s           RA        O     Elevation of the stations (m)
c          wx_s            A*8 A      O     Array of observed weather
c          cover_s          RA        O     Cloud cover (tenths)
c          hgt_ceil         RA        O     Ceiling height (m)
c          hgt_low          RA        O     Height lowest cloud (m)
c          t_s              RA        O     Temperature (F)
c          td_s             RA        O     Dewpoint (F)
c          dd_s             RA        O     Wind direction (deg)
c          ff_s             RA        O     Wind speed (kt)
c          ddg_s            RA        O     Gust wind direction (deg)
c          ffg_s            RA        O     Gust wind speed (kt)
c          pr_s             RA        O     Pressure variable - see note #2.
c          sr_s             RA        O     Solar radiation.
c          istatus          I         O     Status flag: 1 = normal
c                                                       -1 = file not found
c
c       User Notes:
c
c       1.  Arrays should be dimensioned 'maxstns' in the calling program,
c           with maxstns >= 60 for this routine.
c
c       2.  The pressure variable for the Mesonet data is station pressure
c           in millibars.  The pressure variable for SAOs is altimeter setting
c           in millibars.  The mesonet stations are always the first 'num_meso'
c           stations in the pr_s array, and SAOs are the rest.  No corrections,
c           changes, reductions, etc., are made to the data in this routine.
c
c       3.  INFILE is now passed in, complete and ready for the open statement.
c
c*******************************************************************************
c
        real lat_s(maxstns), lon_s(maxstns), elev_s(maxstns)
        real cover_s(maxstns), hgt_ceil(maxstns), hgt_low(maxstns)
        real t_s(maxstns), td_s(maxstns), pr_s(maxstns), sr_s(maxstns)
        real dd_s(maxstns), ff_s(maxstns), ddg_s(maxstns), ffg_s(maxst
     1ns)
c
        character stations(maxstns)*3, wx_s(maxstns)*8
        character atime*24, infile*256
c
        istatus = 0
c
c.....  Open the OBS data file.  Check for a
c.....  'file not found', and notify the user if necessary.
c
        open(1,iostat=ios,file=infile,status='old')
        if(ios .eq. 29) then            !file not found
          print *,
     &    ' +++++ OBS file (SAO & Mesonet data) not available. +++++'
          istatus = -1
          return
        endif
c
c.....  Now read the time and number of stations in the file, then read
c.....  the data.
c
        read(1,901) atime,num_meso,num_sfc
901     format(1x,a24,2i6)
c
        do k=1,num_sfc
          read(1,902)stations(k),lat_s(k),lon_s(k),elev_s(k),wx_s(k),
     &           cover_s(k),hgt_ceil(k),hgt_low(k),t_s(k),td_s(k),
     &           dd_s(k),ff_s(k),ddg_s(k),ffg_s(k),pr_s(k),sr_s(k)
        enddo !k
902     format(1x,A3,2f7.2,1x,f5.0,1x,a8,1x,f5.1,2(1x,f7.1),1x,2f6.1,
     &         1X,4(1X,f5.0),1X,F6.1,1x,f6.1)
c
        num_saos = num_sfc - num_meso
        istatus = 1                     ! normal return
c
        return
        end
