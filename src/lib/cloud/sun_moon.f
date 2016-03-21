
      subroutine sun_eclipse_parms(i4time,rlat_r4,rlon_r4,ht,idebug
     1                            ,alt_r4,azi_r4,dist_m_r4
     1                            ,earth_radius,elgms_r4
     1                            ,r4_mag,r4_obsc,obsc_limbc)

      IMPLICIT REAL*8(A,B,C,D,E,F,G,H,O,P,Q,R,S,T,U,V,W,X,Y,Z)
      include '../../include/astparms.for'
      include 'wa.inc'

      ANGDIF(X,Y)=DMOD(X-Y+9.4247779607694D0,6.2831853071796D0)
     .-3.1415926535897932D0

      ATAN3(X,Y)=DMOD((DATAN2(X,Y)+6.2831853071796D0),6.2831853071796D0)

      CHARACTER BLANK,TWI,MOON,RISE,SET,signm,signs
      CHARACTER*3 MNTH(12)
      REAL*8 M,MAG,LON,LHMSH,
     1  MXG,MYG,MZG,MXR,MYR,MZR,LAT
      REAL*8 NX,NY,NZ,EX,EY,EZ,LST
      INTEGER RAH,DECD,FRAME,ELGMC,ALTDK1
      Character*5 CTT,CTR,CTS,BTT,BMT,C5_BLANK
      character c8_appm*8,c8_apps*8
      DATA MODE/1/,IFL/0/,TIME/0.0D0/,TIMEOL/0.D0/,IPRINT/1/
      DATA FRAME/2/
      DATA RISE/'R'/,SET/'S'/,BLANK/' '/,C5_BLANK/'     '/
      DATA MNTH/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP'
     .,'OCT','NOV','DEC'/

      integer i4time,istatus
      integer i4time_last /0/
      real rlat_r4,rlon_r4,alt_r4,azi_r4,dist_m_r4
      real alm_r4,azm_r4,elgms_r4,r4_mag,ht,earth_radius

      real r4_ratio,r4_obsc,obsc_limbc(nc),solar_eclipse_magnitude

      save SXG,SYG,SZG,MXG,MYG,MZG,i4time_last,ET,UT,LST
      save NX,NY,NZ,EX,EY,EZ,ZNX,ZNY,ZNZ,TX,TY,TZ,TMAG
      save topo_flag

      n_loc = 1

      if(idebug .ge. 1)then
        write(6,*)
        write(6,*)' subroutine sun_eclipse_parms:',i4time,i4time_last
      endif

      rlon = rlon_r4

      if(i4time .ne. i4time_last)then
        write(6,*)
        write(6,*)' initialize sun_eclipse_parms:',i4time,i4time_last

        i4time_last = i4time
 
        call i4time_to_jd(i4time,tb,istatus)
        tf = tb

        ti=ti/1440.d0 ! input in minutes
C
C ENTER LOOPS
        L = 1
        RSN = 1.0
        ITER=3

C
C ENTER TIME LOOPS
        t = tb

        DELTAT=delta_t(T)
        UT = T
        ET = UT + DELTAT

!                         day  deg rad
        call sidereal_time(UT,rlon,LST)

        do 1000 i = 1,iter
C
C CALCULATE POSITION OF EARTH (1950 coordinates - antedated)
840       continue
!         CALL POSINT(ET-RSN/C,1,SXG,SYG,SZG)
          CALL POSIN(ET-RSN/C,1,SXG,SYG,SZG)
          call xyz_to_polar_r(-SXG,-SYG,-SZG,DECS,RAS,RSN)
C         write(13,843)R,SI
843       FORMAT(6F10.6)

1000    CONTINUE
C
C CALCULATE COORDINATES OF SUN (coordinates of date)
        CALL PRECES(T1950,ET,SXG,SYG,SZG,1)
C
C CALCULATE POSITION OF MOON (topocentric coordinates of date)
        CALL MOON_BRWN(ET,MXG,MYG,MZG)

        write(6,*)
        write(6,*)' Sun coords  = ',SXG,SYG,SZG
        write(6,*)' Moon coords = ',MXG,MYG,MZG

!       Calculate Eclipse Conditions at each grid point
        phi = rlat_r4*rpd
!       topo_flag = ht/earth_radius   ! based on mean radius
        topo_flag = ht/(r_e_km*1000.) ! based on equatorial radius
c
c calculate orientation of observer's zenith (coordinates of date)
        CALL TOPO_ff1(phi,rlon,UT,TXZ,TYZ,TZZ)
        TZMAG = SQRT(TXZ**2 + TYZ**2 + TZZ**2) ! based on equ radius
c
c calculate geocentric topo position including flattening factor
        CALL TOPO(phi,rlon,UT,TX,TY,TZ)
        TMAG = SQRT(TX**2 + TY**2 + TZ**2)

!       Observer location relative to Earth center
!       Note that 'ht' is measured along surface normal
        TX = TX + topo_flag * TXZ
        TY = TY + topo_flag * TYZ
        TZ = TZ + topo_flag * TZZ

!       Get direction cosines based on alt/azi (north along the horizon)
!       RA is 180 degrees from RA of meridian
!       DEC is 90 - latitude
        RAN = 180d0 + LST/rpd
        DECN = 90d0 - rlat_r4       
        NX = COSD(DECN) * COSD(RAN)
        NY = COSD(DECN) * SIND(RAN)
        NZ = SIND(DECN)

!       Direction cosines of zenith
        ZNX = TXZ / TZMAG
        ZNY = TYZ / TZMAG
        ZNZ = TZZ / TZMAG

!       East horizon is N horizon cross product with zenith unit vector
        call crossproduct(NX,NY,NZ,ZNX,ZNY,ZNZ,EX,EY,EZ)

        write(6,*)' LST,RLON-LST ',LST/rpd,RLON-LST/rpd ! deg
        write(6,*)' DECN,RAN ',DECN,RAN                 ! deg
        write(6,*)' Zenith unit vector    ',ZNX,ZNY,ZNZ
        write(6,*)' N horizon unit vector ',NX,NY,NZ
        write(6,*)' E horizon unit vector ',EX,EY,EZ

      endif ! new time

      SDIST_M = dist_m_r4 
!     SDIST_AU = (SDIST_M / earth_radius) * TMAG
      SDIST_AU = (SDIST_M / 1000d0) / KM_PER_AU 

      SINALT = sind(alt_r4)
      COSALT = cosd(alt_r4)
      SINAZI = sind(azi_r4)
      COSAZI = cosd(azi_r4)

!     Component of ray along northern horizon
      RX = TX + SDIST_AU * COSAZI * COSALT * NX
      RY = TY + SDIST_AU * COSAZI * COSALT * NY
      RZ = TZ + SDIST_AU * COSAZI * COSALT * NZ

      if(idebug .ge. 2)then
         write(6,*)' Topo (observer) coords = ',TX,TY,TZ
         write(6,*)' Ray coords (N) =         ',RX,RY,RZ
      endif

!     Component of ray along eastern horizon
      RX = RX + SDIST_AU * SINAZI * COSALT * EX
      RY = RY + SDIST_AU * SINAZI * COSALT * EY
      RZ = RZ + SDIST_AU * SINAZI * COSALT * EZ

      if(idebug .ge. 2)then
         write(6,*)' Ray coords (E) = ',RX,RY,RZ
      endif

!     Component of ray along zenith
      RX = RX + SDIST_AU * SINALT * ZNX
      RY = RY + SDIST_AU * SINALT * ZNY
      RZ = RZ + SDIST_AU * SINALT * ZNZ

      if(idebug .ge. 1)then

!        Ray location distance from Earth center (AU)
         RAYMAG = SQRT(RX**2 + RY**2 + RZ**2)

         RAYLAT = ASIND(RZ/RAYMAG) ! geocentric
         RAYLST = ATAN2D(RY,RX)
         RAYLON = RAYLST + (RLON-LST/rpd) ! deg
         if(RAYLON .gt. +180d0)RAYLON = RAYLON - 360d0
         if(RAYLON .lt. -180d0)RAYLON = RAYLON + 360d0
         RAYMSL = ((RAYMAG-TMAG)*KM_PER_AU) * 1000d0

         write(6,*)' rlat,rlon,UT,ht,topo_flag = '
     1              ,rlat_r4,rlon,UT,ht,topo_flag
         write(6,*)' SDIST_M,SDIST_AU = ',SDIST_M,SDIST_AU
         write(6,*)' Ray  coords = ',RX,RY,RZ
         write(6,*)' Ray  lat/lst/lon = ',RAYLAT,RAYLST,RAYLON
         write(6,*)' Ray  mag/tmag/msl = ',RAYMAG,TMAG,RAYMSL
      endif

      SXR = SXG - RX
      SYR = SYG - RY
      SZR = SZG - RZ

      call xyz_to_polar_r(-SXR,-SYR,-SZR,DECS,RAS,RSN)
      HAS=angdif(RAMR,RAS)
!     write(13,*)SZG,rsn,DECS/rpd,RAS/rpd,has/rpd

C
C CALCULATE ALT AND AZ of SUN
      call equ_to_altaz_r(DECS,HAS,PHI,ALS,AZS)
      als = als/rpd
      call refract(als,apps,pres)

      if(als .lt. -1.0)then
          c8_apps = '        '
      else
          write(c8_apps,1002)apps
1002      format(f8.2)
      endif

      azs = azs/rpd

      MXR=MXG-RX
      MYR=MYG-RY
      MZR=MZG-RZ
      call xyz_to_polar_r(MXR,MYR,MZR,DECM,RAM,RMN)

      HAM=angdif(RAMR,RAM)

      CALL anglevectors(-MXR,-MYR,-MZR,SXR,SYR,SZR,ELGARG)
      ELGMST = ELGARG/RPD

      CALL anglevectors(-MXG,-MYG,-MZG,SXG,SYG,SZG,ELGARG)
      ELGMSG = ELGARG/RPD

      if(idebug .ge. 2)then
         write(6,*)' DECS / RAS / HAS = '
     1              ,DECS/rpd,RAS/rpd,HAS/rpd
         write(6,*)' DECM / RAM / HAM / ELSMST = '
     1              ,DECM/rpd,RAM/rpd,HAM/rpd,ELGMST
      endif

!     Calculate Apparent Positions
      call equ_to_altaz_r(DECM,HAM,PHI,ALM,AZM)

      alm = alm/rpd
      call refract(alm,appm,pres)

      if(alm .lt. -1.0)then
            c8_appm = '       '
      else
          write(c8_appm,1001)appm
1001      format(f8.2)
      endif

      azm = azm/rpd

      CALL CJYMD(UT,IYEAR,MONTH,DATE)

      idt = int(date)
      frac = (date - int(date)) * 2.d0 * pi
      CALL CLOCK(frac,CTT)

      call coords(decm,signm,dec_dg_m,dec_mn_m,dec_sc_m,
     1            ram,       ra_hr_m, ra_mn_m, ra_sc_m)

      call coords(decs,signs,dec_dg_s,dec_mn_s,dec_sc_s,
     1            ras,       ra_hr_s, ra_mn_s, ra_sc_s)


      alm_r4 = alm
      azm_r4 = azm
!     elgms_r4 = 10.
      elgms_r4 = elgmst

      phase_angle_deg = 180. - elgmst ! approximate

!     call phase_func_moon(phase_angle_deg,mode,area_rel,area_es
!    1                    ,phase_corr_moon)

      if(.false.)then
!         Calculate Solar Eclipse Magnitude
!         call magnitude(0,0,SXT,SYT,SZT,0.,0.,0.,amag
!    1                                          ,diam_sun)
!         call magnitude(1,1,SXT,SYT,SZT,
!    1      SXG+MXG,SYG+MYG,SZG+MZG,amag,diam_moon)

      endif ! .true. 

      diam_sun = 1800.
      diam_moon = 1852.2
      elgsec = elgmst * 3600.

      overlap_sec = -(elgsec - 0.5 * (diam_sun + diam_moon))
      if(overlap_sec .gt. 0.)then
          solar_eclipse_magnitude = overlap_sec / diam_sun
!         Larger term will increase exponent
!         fracmag = solar_eclipse_magnitude**30.
!         exponent = 1.35*(1.-fracmag) + 0.5*fracmag
!         r4_obsc = solar_eclipse_magnitude**exponent
!         r4_obsc = min(r4_obsc,1.0)
          r4_ratio = diam_moon/diam_sun
          call get_obscuration(solar_eclipse_magnitude,r4_ratio
     1                        ,r4_obsc,obsc_limbc)
          if(r4_obsc/r4_obsc .ne. 1.0)then
              write(6,*)' ERROR in sun_eclipse_parms r4_obsc ='
     1                     ,r4_obsc,solar_eclipse_magnitude,r4_ratio
              stop
          endif
      else
          solar_eclipse_magnitude = 0.
          r4_obsc = 0.
      endif

      r4_mag = solar_eclipse_magnitude
  
      if(idebug .ge. 2)then
          write(6,*)' diam_sun,diam_moon,overlap_sec,r4_mag ',
     1                diam_sun,diam_moon,overlap_sec,r4_mag
          write(6,*)' return from sun_eclipse_parms'
      endif
C
      RETURN
      END

      subroutine sun_moon(i4time,lat_2d,lon_2d,ni,nj,is,js,alm_r4,azm_r4
     1                   ,idebug                           ! I
     1                   ,ht,earth_radius                  ! I
     1                   ,elgms_r4,r4_mag,r4_rmn           ! O
     1                   ,solar_eclipse_magnitude,r4_obsc,obsc_limb) ! O   

      IMPLICIT REAL*8(A,B,C,D,E,F,G,H,O,P,Q,R,S,T,U,V,W,X,Y,Z)
!     include '../util/utilparms.for'
      include '../../include/astparms.for'
      include 'wa.inc'

      ANGDIF(X,Y)=DMOD(X-Y+9.4247779607694D0,6.2831853071796D0)
     .-3.1415926535897932D0

      ATAN3(X,Y)=DMOD((DATAN2(X,Y)+6.2831853071796D0),6.2831853071796D0)

      CHARACTER BLANK,TWI,MOON,RISE,SET,signm,signs
      CHARACTER*3 MNTH(12)
      REAL*8 I,M,MAG,LON,LHMSH,
     1  MXG,MYG,MZG,MXT,MYT,MZT,LAT
      INTEGER RAH,DECD,FRAME,ELGMC,ALTDK1
      Character*5 CTT,CTR,CTS,BTT,BMT,C5_BLANK
      character c8_appm*8,c8_apps*8
      character*20 c20_site
      DIMENSION LAT(9),LON(9)
      DATA MODE/1/,IFL/0/,TIME/0.0D0/,TIMEOL/0.D0/,IPRINT/0/
      DATA FRAME/2/
      DATA RISE/'R'/,SET/'S'/,BLANK/' '/,C5_BLANK/'     '/
      DATA MNTH/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP'
     .,'OCT','NOV','DEC'/

      integer i4time,istatus
      real lat_2d(ni,nj)
      real lon_2d(ni,nj)
      real alm_r4,azm_r4,elgms_r4,r4_mag,r4_rmn
      real solar_eclipse_magnitude,ht,earth_radius
      real r4_ratio,r4_obsc,obsc_limb,obsc_limbc(nc)
C
      n_loc = 1
      lat = lat_2d(is,js)
      lon = lon_2d(is,js)

      call i4time_to_jd(i4time,tb,istatus)
      tf = tb
      c20_site = 'las'

!     Read parameters
!     open(11,file='sun_moon.parms',status='old')
!     read(11,*)n_loc                           ! Coordinates
!     do idum = 1,n_loc
!       read(11,101)c20_site
!101     format(a20)
!       read(11,*)LAT(idum),LON(idum),pres
!     enddo
!     read(11,*)TB                              ! Time interval
!     read(11,*)TF
!     read(11,*)TI
!     read(11,*)topo_flag

      ti=ti/1440.d0 ! input in minutes
C
C ENTER LOOPS
      L = 1
      PHI=LAT(L)*RPD
      RSN = 1.0
      SINPHI=DSIN(PHI)
      COSPHI=DCOS(PHI)
      ITER=3

C
C ENTER TIME LOOPS
      t = tb
C
C WRITE HEADINGS
      IF(IPRINT.EQ.1.AND.IFL.EQ.0)write(6,11)c20_site,LAT(L),LON(L)
11    FORMAT(1H1,30X,a20,'  LAT=',F7.2,'  LON=',F7.2/)
10    IF(IPRINT.EQ.1.AND.IFL.EQ.0)write(6,1)
1     FORMAT(
     1 ' YEAR MON DT   UT                   MOON         ',
     1 '     Coordinates of Date        SUN',
     1           26x,'        ELONG  ALTDIF'/
     1            21x,'Alt    App     Az        Dec          RA',
     1      10x,'Alt    App     Az        Dec          RA',7x,'(mag)')


      DELTAT=delta_t(T)
      UT = T
      ET = UT + DELTAT
c
c calculate orientation of observer's zenith (coordinates of date)
      CALL TOPO_ff1(PHI,LON(L),UT,TXZ,TYZ,TZZ)
      RAMR=ATAN3(TYZ,TXZ)

      IF(IPRINT.EQ.1)write(6,*)' lat,lon = ',LAT(L),LON(L)
      IF(IPRINT.EQ.1)write(6,*)' TXZ,TYZ,TZZ,RAMR',TXZ,TYZ,TZZ,RAMR/rpd

      topo_flag = 1.0 + ht/earth_radius

c
c calculate geocentric topo position including flattening factor
      CALL TOPO(PHI,LON(L),UT,TX,TY,TZ)

      TX = TX * topo_flag
      TY = TY * topo_flag
      TZ = TZ * topo_flag

      do 1000 i = 1,iter
C
C CALCULATE POSITION OF EARTH (1950 coordinates - antedated)
840   continue
!     CALL POSINT(ET-RSN/C,1,SXG,SYG,SZG)
      CALL POSIN(ET-RSN/C,1,SXG,SYG,SZG)
      SXT = SXG - TX
      SYT = SYG - TY
      SZT = SZG - TZ
C     write(13,843)R,SI
843   FORMAT(6F10.6)

1000  CONTINUE

C
C CALCULATE COORDINATES OF SUN (coordinates of date)
      CALL PRECES(T1950,ET,SXG,SYG,SZG,1)
      CALL PRECES(T1950,ET,SXT,SYT,SZT,1)

      call xyz_to_polar_r(-SXT,-SYT,-SZT,DECS,RAS,RSN)
      HAS=angdif(RAMR,RAS)
!     write(13,*)SZG,rsn,DECS/rpd,RAS/rpd,has/rpd

C
C CALCULATE ALT AND AZ of SUN
      call equ_to_altaz_r(DECS,HAS,PHI,ALS,AZS)
      als = als/rpd
      call refract(als,apps,pres)

      if(als .lt. -1.0)then
          c8_apps = '        '
      else
          IF(IPRINT.EQ.1)write(c8_apps,1002)apps
1002      format(f8.2)
      endif


      azs = azs/rpd

C
C CALCULATE POSITION OF MOON (topocentric coordinates of date)
      CALL MOON_BRWN(ET,MXG,MYG,MZG)
      MXT=MXG-TX
      MYT=MYG-TY
      MZT=MZG-TZ
      call xyz_to_polar_r(MXT,MYT,MZT,DECM,RAM,RMN)
      r4_rmn = RMN

!     Insert star in lunar position
!     DECM = ?
!     RAM = ?
!     CALL PRECES(2433282.423357D0,ET,DECM,RAM,RDUM,2)

      HAM=angdif(RAMR,RAM)

      CALL anglevectors(-MXT,-MYT,-MZT,SXT,SYT,SZT,ELGARG)
      ELGMST = ELGARG/RPD

      CALL anglevectors(-MXG,-MYG,-MZG,SXG,SYG,SZG,ELGARG)
      ELGMSG = ELGARG/RPD

      IF(IPRINT.EQ.1)write(6,*)' DECM / RAM / HAM = '
     1                          ,DECM/rpd,RAM/rpd,HAM/rpd

!     Calculate Apparent Positions
      call equ_to_altaz_r(DECM,HAM,PHI,ALM,AZM)

      alm = alm/rpd
      call refract(alm,appm,pres)

      if(alm .lt. -1.0)then
          c8_appm = '       '
      else
          IF(IPRINT.EQ.1)write(c8_appm,1001)appm
1001      format(f8.2)
      endif

      azm = azm/rpd

      CALL CJYMD(UT,IYEAR,MONTH,DATE)

      idt = int(date)
      frac = (date - int(date)) * 2.d0 * pi
      CALL CLOCK(frac,CTT)

      call coords(decm,signm,dec_dg_m,dec_mn_m,dec_sc_m,
     1          ram,       ra_hr_m, ra_mn_m, ra_sc_m)

      call coords(decs,signs,dec_dg_s,dec_mn_s,dec_sc_s,
     1          ras,       ra_hr_s, ra_mn_s, ra_sc_s)


      alm_r4 = alm
      azm_r4 = azm
      elgms_r4 = elgmst

      phase_angle_deg = 180. - elgmst ! approximate

      call phase_func_moon(phase_angle_deg,mode,area_rel,area_es
     1                    ,phase_corr_moon)

      r4_mag = -12.74 + phase_corr_moon

      if(.true.)then

!         Calculate Solar Eclipse Magnitude
!         call magnitude(0,0,SXT,SYT,SZT,0.,0.,0.,amag
!    1                                          ,diam_sun)
!         call magnitude(1,1,SXT,SYT,SZT,
!    1      SXG+MXG,SYG+MYG,SZG+MZG,amag,diam_moon)

          diam_sun = 1800.
          diam_moon = 1852.2
          elgsec = elgmst * 3600.

          overlap_sec = -(elgsec - 0.5 * (diam_sun + diam_moon))
          solar_eclipse_magnitude = overlap_sec / diam_sun

          if(solar_eclipse_magnitude .gt. 0.)then
!             Larger term will increase exponent
              fracmag = solar_eclipse_magnitude**30.
              exponent = 1.35*(1.-fracmag) + 0.5*fracmag
!             r4_obsc = solar_eclipse_magnitude**exponent
!             r4_obsc = min(r4_obsc,1.0)
              r4_ratio = diam_moon/diam_sun
              call get_obscuration(solar_eclipse_magnitude,r4_ratio
     1                            ,r4_obsc,obsc_limbc)
              obsc_limb = obsc_limbc(2)
          else
              r4_obsc = 0.
          endif

          if(idebug .ge. 2)then
             write(6,*)
             write(6,*)' sun_moon debug:'
             write(6,*)' i4time,tb,lat,lon',i4time,tb,lat(1),lon(1)
             write(6,*)' deltat/sec',DELTAT,DELTAT*86400.
             write(6,*)' elgmst/sec',elgmst,elgmst*3600.
             write(6,*)' diam_sun,diam_moon,overlap_sec,mag,r4_obsc ',
     1diam_sun,diam_moon,overlap_sec,solar_eclipse_magnitude,r4_obsc
             write(6,*)' DECS / RAS / HAS = '
     1                  ,DECS/rpd,RAS/rpd,HAS/rpd
             write(6,*)' DECM / RAM / HAM / ELSMST = '
     1                  ,DECM/rpd,RAM/rpd,HAM/rpd,ELGMST
             write(6,*)' ET / MXT / MYT / MZT = '
     1                  ,ET,MXT,MYT,MZT
             write(6,*)
          endif

      endif ! .true.

!     Correct moon's magnitude for lunar eclipse
      if(elgmsg .gt. 177.)then

!         Calculate Lunar Eclipse Magnitude
          delta_moon_au = mag(MXG,MYG,MZG)
          delta_sun_au  = mag(SXG,SYG,SZG)
          delta_moon_km = delta_moon_au * km_per_au
          delta_sun_km  = delta_sun_au  * km_per_au

          diam_moon_km = r_m_km * 2d0

          elgsupl = 180d0 - elgmsg

          x_moon =     delta_moon_km * cosd(elgsupl)
          y_moon = abs(delta_moon_km * sind(elgsupl))

          umbral_length_km = delta_sun_km / (R_S_km / R_E_km - 1d0)

          angrad = asin(r_e_km/umbral_length_km)
          umbral_height_km = (umbral_length_km - x_moon)
     1                              * tan(angrad) * 1.02

          dist_center_km = (y_moon - umbral_height_km) * cos(angrad)
          umbral_magnitude = -(dist_center_km / diam_moon_km) + 0.5

          penumbral_length_km = delta_sun_km / (R_S_km / R_E_km + 1d0)

          angrad = asin(r_e_km/penumbral_length_km)
          penumbral_height_km  = (penumbral_length_km + x_moon)
     1                              * tan(angrad)

          dist_center_km = (y_moon - penumbral_height_km) * cos(angrad)
          penumbral_magnitude = -(dist_center_km / diam_moon_km) + 0.5

          penumbral_width = penumbral_magnitude - umbral_magnitude

          bri_2nd = .0002

          if(umbral_magnitude .gt. 1.)then          ! Total    
              IF(IPRINT.EQ.1)write(6,*)' Total'           
              IF(IPRINT.EQ.1)write(6,*)' umbral_magnitude = '
     1                                  ,umbral_magnitude
              frac_bri = bri_2nd / umbral_magnitude**4
          elseif(umbral_magnitude .gt. 0.)then      ! Partial
              IF(IPRINT.EQ.1)then
                write(6,*)' Partial'         
                write(6,*)' umbral_magnitude = ',umbral_magnitude
                write(6,*)' penumbral_magnitude = ',penumbral_magnitude
                write(6,*)' penumbral_width = ',penumbral_width
              endif
              frac_bri = ((1.0 - umbral_magnitude) * 0.5) 
     1                 / penumbral_width
!             frac_bri = max(frac_bri,bri_2nd) 
              frac_bri = frac_bri + umbral_magnitude * bri_2nd
          elseif(penumbral_magnitude .gt. 0.)then   ! Penumbral
              IF(IPRINT.EQ.1)then
                write(6,*)' Penumbral'       
                write(6,*)' penumbral_magnitude = ',penumbral_magnitude
                write(6,*)' penumbral_width = ',penumbral_width
              endif
              frac_bri = 1.0 - (penumbral_magnitude * 0.5)
     1                       / penumbral_width
          else
              frac_bri = 1.0
          endif

          rmag_corr = -log10(frac_bri) * 2.5

          r4_mag = r4_mag + rmag_corr

          IF(IPRINT.EQ.1)write(6,451)frac_bri,rmag_corr,r4_mag
451       format(' Lunar Eclipse: frac_bri/rmag_corr/r4_mag'
     1          ,f11.6,2f9.3)

      endif ! Lunar eclipse magnitude correction

C
C WRITE OUT DATA
500   continue
      IF(IPRINT.EQ.1)write(6,2)IYEAR,MNTH(MONTH),idt,ctt,
     1  alm,c8_appm(2:8),azm,signm,int(abs(dec_dg_m)),int(dec_mn_m),dec_
     1sc_m,
     1                        int(ra_hr_m),      int(ra_mn_m), ra_sc_m,
     1  als,c8_apps(2:8),azs,signs,int(abs(dec_dg_s)),int(dec_mn_s),dec_
     1sc_s,
     1                        int(ra_hr_s),      int(ra_mn_s), ra_sc_s,
     1                          elgmst,alm-als
2     format(1x,i4,1x,a3,1x,i2,1x,a5,
     1  f7.2,a7,f8.2,2x,a1,i2,1x,i2,f5.1,2x,i2,1x,i2,f6.2,
     1          2x,
     1  f7.2,a7,f8.2,2x,a1,i2,1x,i2,f5.1,2x,i2,1x,i2,f6.2,
     1          1x,f7.2,f7.2)

C
9999  CONTINUE

      RETURN
      END


        subroutine coords(dec_rd,sign,dec_dg,dec_mn,dec_sc,
     1                  ra_rd,       ra_hr, ra_mn, ra_sc)

        IMPLICIT REAL*8(A,B,C,D,E,F,G,H,O,P,Q,R,S,T,U,V,W,X,Y,Z)
        include '../../include/astparms.for'

        character*1 sign,plus/'+'/,minus/'-'/

        if(dec_rd .ge. 0)then
            sign = plus
        else
            sign = minus
        endif

        dec_dg = abs(dec_rd / rpd)
        dec_mn = (dec_dg - int(dec_dg)) * 60d0
        dec_sc = (dec_mn - int(dec_mn)) * 60d0

        ra_hr = ra_rd / rpd / 15d0
        ra_mn = (ra_hr - int(ra_hr)) * 60d0
        ra_sc = (ra_mn - int(ra_mn)) * 60d0

        return
        end

        subroutine get_obscuration(rmag,ratio,obscuration,obsc_limbc)

!       Eclipse obscuration where ratio is moon/sun angular diameter

!       Limb darkening reference (Hestroffer and Magnan, 1997)
!       www.physics.hmc.edu/faculty/esin/a101/limbdarkening.pdf

!       Optional switch to first Pierce Method in Astrophysical Quantities 
        include 'wa.inc'

        parameter (method = 3) ! first Pierce Method

        real obsc_limbc(nc)
        real a_a(nc)  / 0.78, 0.74, 0.54/
        real b_a(nc)  / 0.39, 0.43, 0.60/
        real c_a(nc)  /-0.57,-0.56,-0.44/

!       Mean brightness divided by central brightness
        real cb(nc)   /.855 ,.850, .820/ ! larger increases obsc

        parameter (pi = 3.14159265)
        seg_area(theta,r) = 0.5 * r**2 * (theta - sin(theta))

        real mu
        mu(r) = sqrt(1.-r**2)
        rlimbdk1(r,u,alphal) = 1. - (u * (1. - (mu(r))**alphal))
        rlimbdk3(r,a,b,c) = a + b*mu(r) 
     1                        + c * (1.0-mu(r) * log(1.0 + 1./mu(r)))

        u = 1.00

        if(rmag .ge. 0.99999)then
            obscuration = 1.0
            obsc_limbc(:) = 1.0         
            return
        elseif(rmag .eq. 0.)then 
            obscuration = 0.0
            obsc_limbc(:) = 0.0            
        endif

        rs = 1.0
        rm = rs*ratio

        oa = -(rmag * (2.*rs) - (rs + rm))
        a = (rs**2 - rm**2 + oa**2) / (2. * oa)

        alpha = 2. * acos(a/rs)
        beta  = 2. * acos((oa-a)/rm)

        seg1 = seg_area(alpha,rs)
        seg2 = seg_area(beta,rm)  

        obscuration = (seg1 + seg2) / pi

        if(rmag .le. 0.5)then
            obsc_limbc(:) = obscuration
        else
!           calculate mean radius of "quadratic" segment
            rmean = 1. - ((1.-rmag)*(2./3.)) 
            do ic = 1,nc
!               alphal = 0.8 
                alphal = -0.023 + .292 / wa(ic)
                if(method .eq. 1)then
                    rmean_ill = rlimbdk1(rmean,u,alphal)
                else
                    rmean_ill = rlimbdk3(rmean,a_a(ic),b_a(ic),c_a(ic))
                    rmean_ill = rmean_ill / cb(ic)
                endif
                obsc_limbc(ic) = 1.0 - ((1.0 - obscuration) * rmean_ill)
            enddo ! ic
        endif

!       write(6,*)'rmag/obscuration/obsc_limbc = '
!    1            ,rmag,obscuration,obsc_limbc
!       write(6,*)'rmean/rmean_ill/mu = ',rmean,rmean_ill,mu(rmean)
!       write(6,*)'mu(r)/mu(r)**alphal',mu(rmean),mu(rmean)**alphal
!       write(6,*)'rmean_ill derived',(1.-obsc_limbc)/(1.-obscuration)
!       stop

!       write(6,1)rmag,ratio,oa,a
!1      format('rmag/ratio/oa/a',4f10.5)
!       write(6,2)alpha,beta,seg1,seg2,obscuration
!2      format('alpha/beta/seg1/seg2/obsc',5f10.5)
!       stop

        return
        end
