
      subroutine sun_planet(i4time,n_planet,rlat,rlon
     1                     ,dec_r4,ra_r4,alm_r4,azm_r4,elgms_r4,r4_mag) 

      IMPLICIT REAL*8(A,B,C,D,E,F,G,H,O,P,Q,R,S,T,U,V,W,X,Y,Z)

      include '../../include/astparms.for'
!     include '../planets/posint_inc.for'

      ANGDIF(X,Y)=DMOD(X-Y+9.4247779607694D0,6.2831853071796D0)
     .-3.1415926535897932D0

      ATAN3(X,Y)=DMOD((DATAN2(X,Y)+6.2831853071796D0),6.2831853071796D0)

      CHARACTER BLANK,TWI,MOON,RISE,SET,SIGN,PLUS,MINUS
      DIMENSION R(3),RR(3),RHO(3),RRI(3),P(3),Q(3),W(3),PP(3),QQ(3)
     +,WW(3),D1(3),D2(3),I1(3),I2(3)
      CHARACTER*3 MNTH(12)
      REAL*8 I,M,MAG,LNOD,MU,NLC,LON,MDEGTL,LHMSH,MX,MY,MZ,LAT
      REAL*8 MX_1950, MY_1950, MZ_1950
      INTEGER RAH,DECD,FRAME,ELGMC,ALTDK1
      Character*5 CTT,CTR,CTS,BTT,BMT,C5_BLANK
      character c8_appm*8,c8_apps*8
      character*1 c_observe
      character*20 c20_site
      character*8  names(13)
!     character*8  names(13)/'Earth','Mercury','Venus','Mars','Jupiter','Saturn',' ',' ',' ',' ',' ',' ',' '/
      DIMENSION LAT(9),LON(9)

      real*8 maglimd_r8,maglimt_r8,maglimn_r8,maglim_r8
      real alm_r4,azm_r4,r4_mag,elgms_r4,rlat,rlon,dec_r4,ra_r4

!     CUBERT(X)=DEXP(DLOG(DABS(X))/3.)*X/DABS(X)

      DATA MODE/1/,IFL/0/,TIME/0.0D0/,TIMEOL/0.D0/,IPRINT/1/
      DATA FRAME/2/
      DATA RISE/'R'/,SET/'S'/,BLANK/' '/,C5_BLANK/'     '/
      DATA PLUS/'+'/,MINUS/'-'/
      DATA MNTH/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP'
     .,'OCT','NOV','DEC'/
      Data names/'Earth','Mercury','Venus','Mars','Jupiter','Saturn',
     1         'Uranus','Neptune','Pluto',4*' '/

!     Saturn Ring Direction Cosines
!     data zx,zy,zz/-.0912836831D0, -.0724471759D0, -.9931861335D0/

!     Read parameters
!     open(11,file='sun_planet.parms',status='old')
!     read(11,*)n_loc,n_planet                          ! Coordinates
!     do idum = 1,n_loc
!       read(11,101)c20_site
101     format(a20)
!       read(11,*)LAT(idum),LON(idum),pres
!     enddo
!     read(11,*)TB                              ! Time interval
!     read(11,*)TF
!     read(11,*)TI

      lat(1) = rlat
      lon(1) = rlon

      call i4time_to_jd(i4time,tb,istatus)
      tf = tb
      c20_site = 'las'

      ti=ti/1440.d0 ! input in minutes
      RPH=PI/12.
      ARGP=ARGP*RPD
      LNOD=LNOD*RPD
      I=I*RPD
      O=OBLIQ1950*RPD
C
C ENTER LOOPS
      L = 1
      PHI=LAT(L)*RPD
      RSN = 1.0
      SINPHI=DSIN(PHI)
      COSPHI=DCOS(PHI)
      ITER=3

C ENTER TIME LOOPS
      t = tb
C
C WRITE HEADINGS
      IF(IPRINT.EQ.1.AND.IFL.EQ.0)write(13,11)c20_site,LAT(L),LON(L)
11    FORMAT(1H1,30X,a20,'  LAT=',F7.2,'  LON=',F7.2/)
10    IF(IPRINT.EQ.1.AND.IFL.EQ.0)write(13,1)names(n_planet)
1     FORMAT(
     1 '  YEAR MON DT   UT         ',a8,'        ',
     1 '          Coordinates of Date       SUN',
     1            21x,'   ELONG   ALTDIF ALTDIF2   MAG    Ill    Diam'
     1               ,'   MGLM    VIS'/
     1 '                       Alt     App     Az       Dec    RA  ',
     1          '        Alt     App     Az       Dec    RA')


400   if((t - tf) * ti .gt. 0.)goto9999

      DELTAT=DELTA_T(T)
      UT = T
      ET = UT + DELTAT
c
c calculate position of observer's zenith (coordinates of date)
      CALL TOPO_ff1(PHI,LON(L),UT,TXZ,TYZ,TZZ)
      RAMR=ATAN3(TYZ,TXZ)

      do iter = 1,3 ! Light time iteration

C CALCULATE POSITION OF EARTH (coordinates of date)
          TC = ET - Delta_planet/C
          CALL POSINT(TC,1,RRI(1),RRI(2),RRI(3))
!         CALL POSIN(TC,1,RRI(1),RRI(2),RRI(3))
          CALL PRECES(T1950,ET,RRI(1),RRI(2),RRI(3),1)
C         write(13,843)R,RRI
843       FORMAT(6F10.6)


C CALCULATE POSITION OF PLANET (coordinates of date)
          TC = ET - Delta_planet/C
          CALL POSINT(TC,n_planet,R(1),R(2),R(3))
!         CALL POSIN(TC,n_planet,R(1),R(2),R(3))
          CALL PRECES(T1950,ET,R(1),R(2),R(3),1)
          delta_planet = sqrt((r(1)-rri(1))**2 + (r(2)-rri(2))**2
     1                                         + (r(3)-rri(3))**2 )
C
C CALCULATE COORDINATES OF SUN (coordinates of date)
          TC = ET - RSN/C
          CALL POSINT(TC,1,RRI1,RRI2,RRI3)
!         CALL POSIN(TC,1,RRI1,RRI2,RRI3)
          CALL PRECES(T1950,ET,RRI1,RRI2,RRI3,1)
          call xyz_to_polar_r(-RRI1,-RRI2,-RRI3,DECS,RAS,RSN)

      enddo ! light time iter

      HAS=angdif(RAMR,RAS)
!     write(13,*)rri(3),rsn,DECS/rpd,RAS/rpd,has/rpd

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

!     Relative position is planet vector minus earth vector (coords of date)
      CALL TOPO(PHI,LON(L),UT,TX,TY,TZ)
      MX=R(1)-(RRI(1)+TX)
      MY=R(2)-(RRI(2)+TY)
      MZ=R(3)-(RRI(3)+TZ)

      MX_1950 = MX
      MY_1950 = MY
      MZ_1950 = MZ

!     Precess planet minus Earth vector to coordinates of date
!     CALL PRECES(T1950,ET,MX,MY,MZ,1)

      call xyz_to_polar_r(MX,MY,MZ,DECM,RAM,RMN)

!     Insert star in lunar position
!     DECM = ?
!     RAM = ?
!     CALL PRECES(2433282.423357D0,ET,DECM,RAM,RDUM,2)

      HAM=angdif(RAMR,RAM)
      CALL anglevectors(-MX,-MY,-MZ,RRI(1),RRI(2),RRI(3),ELGARG)
      ELGMS = ELGARG/RPD

      call equ_to_altaz_r(DECM,HAM,PHI,ALM,AZM)

      alm = alm/rpd
      call refract(alm,appm,pres)

      if(alm .lt. -1.0)then
          c8_appm = '        '
      else
          write(c8_appm,1001)appm
1001      format(f8.2)
      endif

      azm = azm/rpd

      CALL CJYMD(UT,IYEAR,MONTH,DATE)

      idt = int(date)
      frac = (date - int(date)) * 2.d0 * pi
      CALL CLOCK(frac,CTT)

      call MAGNitude(N_planet,0,RRI(1),RRI(2),RRI(3),r(1),r(2),r(3)
     1                                  ,amag,diam_sec)
      call phase(RRI(1),RRI(2),RRI(3),r(1),r(2),r(3),phase_angle,r_ill)

      if(n_planet .eq. 6)then
!         Calculate Saturn's Ring Opening
          call anglevectors(ZX,ZY,ZZ,R(1),R(2),R(3),BPRIME_R)
          call anglevectors(ZX,ZY,ZZ,MX_1950,MY_1950,MZ_1950,B_R)
          BPRIME_D = BPRIME_R / rpd - 90.
          B_D      = B_R      / rpd - 90.
!         write(6,*)B_D,BPRIME_D
          if(B_D * BPRIME_D .lt. 0.)then ! Dark side visible
              product = -B_D * BPRIME_D
              write(6,3)IYEAR,MNTH(MONTH),idt,ctt,B_D,BPRIME_D,product
 3            format(1x,i5,1x,a3,1x,i2,1x,a5,f9.3,f9.3,f9.4)

              write(14,3)IYEAR,MNTH(MONTH),idt,ctt,B_D,BPRIME_D,product

              write(13,2)IYEAR,MNTH(MONTH),idt,ctt,
     1          alm,c8_appm,azm,decm/rpd,ram/rpd,
     1          als,c8_apps,azs,decs/rpd,ras/rpd,
     1                  elgms,alm-als,amag,r_ill

          endif
      endif

!     Convert to real
      alm_r4 = alm
      azm_r4 = azm
      r4_mag = amag
      elgms_r4 = elgms
      dec_r4 = decm/rpd
      ra_r4 = ram/rpd

      return

C     CALL VISIBILITY ROUTINE
      altdf = alm - als
      call vi(amag,elgms,altdf,al1_r8,al2_r8,0.D0,180.D0
     1       ,maglimd_r8,maglimt_r8,maglimn_r8,maglim_r8
     1       ,vis_r8,c_observe)

C
C WRITE OUT DATA
      if(n_planet .ne. 6)then

          altdf2 = (alm-als) - (amag*2.0)

          write(13,2)IYEAR,MNTH(MONTH),idt,ctt,
     1          alm,c8_appm,azm,decm/rpd,ram/rpd,
     1          als,c8_apps,azs,decs/rpd,ras/rpd,
     1                  elgms,alm-als,altdf2,amag,r_ill,
     1                  diam_sec,
     1                  maglim_r8,vis_r8,c_observe
2         format(1x,i5,1x,a3,1x,i2,1x,a5,
     1          f8.2,a8,f8.2,f8.3,f8.3,
     1          4x,
     1          f8.2,a8,f8.2,2x,f6.2,f8.2,
     1          1x,f7.1,2f8.3,f8.2,f7.3,f7.2,f7.1,f7.1,a1)

      endif

!     t = t + ti
!     goto400

C
9999  return
      END
