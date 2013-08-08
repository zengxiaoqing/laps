
      subroutine sun_moon(i4time,lat_2d,lon_2d,ni,nj,is,js,alm_r4,azm_r4
     1                   ,elgms_r4,r4_mag)                                     

      IMPLICIT REAL*8(A,B,C,D,E,F,G,H,O,P,Q,R,S,T,U,V,W,X,Y,Z)
!     include '../util/utilparms.for'
      include '../../include/astparms.for'

      ANGDIF(X,Y)=DMOD(X-Y+9.4247779607694D0,6.2831853071796D0)
     .-3.1415926535897932D0

      ATAN3(X,Y)=DMOD((DATAN2(X,Y)+6.2831853071796D0),6.2831853071796D0)

      CHARACTER BLANK,TWI,MOON,RISE,SET,signm,signs
      DIMENSION MNTH(12)
      REAL*8 I,M,MAG,LON,LHMSH,
     1  MXG,MYG,MZG,MXT,MYT,MZT,LAT
      INTEGER RAH,DECD,FRAME,ELGMC,ALTDK1
      Character*5 CTT,CTR,CTS,BTT,BMT,C5_BLANK
      character c8_appm*8,c8_apps*8
      character*20 c20_site
      DIMENSION LAT(9),LON(9)
      DATA MODE/1/,IFL/0/,TIME/0.0D0/,TIMEOL/0.D0/,IPRINT/1/
      DATA FRAME/2/
      DATA RISE/'R'/,SET/'S'/,BLANK/' '/,C5_BLANK/'     '/
      DATA MNTH/3HJAN,3HFEB,3HMAR,3HAPR,3HMAY,3HJUN,3HJUL,3HAUG,3HSEP
     .,3HOCT,3HNOV,3HDEC/

      integer i4time,istatus
      real lat_2d(ni,nj)
      real lon_2d(ni,nj)
      real alm_r4,azm_r4,elgms_r4,r4_mag
C
      n_loc = 1
      lat = lat_2d(is,js)
      lon = lon_2d(is,js)
      topo_flag = 1.0

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
c calculate position of observer's zenith (coordinates of date)
      CALL TOPO_ff1(PHI,LON(L),UT,TXZ,TYZ,TZZ)
      RAMR=ATAN3(TYZ,TXZ)

      write(6,*)' lat,lon = ',LAT(L),LON(L)
      write(6,*)' TXZ,TYZ,TZZ,RAMR',TXZ,TYZ,TZZ,RAMR/rpd

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
          write(c8_apps,1002)apps
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

!     Insert star in lunar position
!     DECM = ?
!     RAM = ?
!     CALL PRECES(2433282.423357D0,ET,DECM,RAM,RDUM,2)


      HAM=angdif(RAMR,RAM)

      CALL anglevectors(-MXT,-MYT,-MZT,SXT,SYT,SZT,ELGARG)
      ELGMST = ELGARG/RPD

      CALL anglevectors(-MXG,-MYG,-MZG,SXG,SYG,SZG,ELGARG)
      ELGMSG = ELGARG/RPD

      write(6,*)' DEC / RA / HA = ',DECM/rpd,RAM/rpd,HAM/rpd

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

      if(.false.)then

!     Calculate Solar Eclipse Magnitude
      call magnitude(0,0,SXT,SYT,SZT,0.,0.,0.,amag
     1                                          ,diam_sun)
      call magnitude(1,1,SXT,SYT,SZT,
     1  SXG+MXG,SYG+MYG,SZG+MZG,amag,diam_moon)

      overlap_sec = -(elgmst * 3600. - 0.5 * (diam_sun + diam_moon))
      solar_eclipse_magnitude = overlap_sec / diam_sun

      if(solar_eclipse_magnitude .gt. 0.)then
          elgmst = solar_eclipse_magnitude * 100.
      endif

      if(elgmsg .lt. 90.)goto500

!     Calculate Lunar Eclipse Magnitude
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
     1                          * tan(angrad) * 1.02

      dist_center_km = (y_moon - umbral_height_km) * cos(angrad)
      umbral_magnitude = -(dist_center_km / diam_moon_km) + 0.5


      penumbral_length_km = delta_sun_km / (R_S_km / R_E_km + 1d0)

      angrad = asin(r_e_km/penumbral_length_km)
      penumbral_height_km  = (penumbral_length_km + x_moon)
     1                          * tan(angrad)

      dist_center_km = (y_moon - penumbral_height_km) * cos(angrad)
      penumbral_magnitude = -(dist_center_km / diam_moon_km) + 0.5




      if(penumbral_magnitude .gt. 0.)then
          elgmst = penumbral_magnitude * 100.
!         elgmst = umbral_magnitude * 100.
      endif

      if(umbral_magnitude .gt. 0.)then
          elgmst = umbral_magnitude * 100.
      endif

      endif ! .false.

C
C WRITE OUT DATA
500   write(6,2)IYEAR,MNTH(MONTH),idt,ctt,
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
