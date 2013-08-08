

      SUBROUTINE POSINT(T,PLANET,R1,R2,R3)
      IMPLICIT REAL*8 (A,B,C,D,E,F,G,H,O,P,Q,R,S,T,U,V,W,X,Y,Z)
      include '../../include/astparms.for'
      integer*4 aw,pl,dg
      parameter (aw = 59)
      parameter (pl = 13)
      parameter (dg = 13)
      INTEGER S,ARYCEN,PLANET,APOINT,ARYHW,DEGREE,DEGM1,TOPOF,APNTP1
     .,APNTP2,APNTP3,APNTM1,APNTM2
      REAL*8 LMA,KSQ,LAT,LON,LST,MX,MY,MZ
      DIMENSION A1(dg),AP1(dg),B1(dg),G1(dg),BLC(dg),ER(40),EV(39),S(14)
      COMMON X(pl,aw),Y(pl,aw),Z(pl,aw),XDD(pl,aw),YDD(pl,aw),ZDD(pl,aw)
     .,XD(pl),YD(pl),ZD(pl),KSQ,A(dg),AP(dg),B(dg),G(dg),TCENT,H
     .,ARG3                ,IEGREE,NPINMS
     .,INIT,MAXPNT,MINPNT,ARYCEN,ARYHW,NUMBER,NRECEN,NUMB1,IARG1,IARG3

      Character*14 C_FILENAME,c_fileref
      Character*7  C_file1 
      DATA C_FILENAME /'posint.2440800'/
      DATA C_file1    /'posout.'/

      DATA INIT1/0/,NPLAN/13/,NMASS/12/,MOON/1/,TOPOF/0/,NREL/9/
      DATA EMRATP1/82.3007D0/
      DATA AP1/4621155471343.D0,-13232841914856.D0,47013743726958.D0
     .,-114321700672600.D0,202271967611865.D0,-266609549584656.D0
     .,264429021895332.D0,-197106808276656.D0,108982933333425.D0
     .,-43427592828040.D0,11807143978638.D0,-1962777574776.D0
     .,150653570023.D0/
      DATA B1/9072652009253.D0,-39726106418680.D0,140544566352762.D0
     .,-344579280210120.D0,613137294629235.D0,-811345852376496.D0
     .,807012356281740.D0,-602852367932304.D0,333888089374395.D0
     .,-133228219027160.D0,36262456774618.D0,-6033724094760.D0
     .,463483373517.D0/
      DATA A1/150653570023.D0,2662659061044.D0,-1481863453062.D0
     .,3926822700380.D0,-6604398106155.D0,8380822992264.D0
     .,-8088023425188.D0,5907495735864.D0,-3215663657055.D0
     .,1265630766980.D0,-340671801462.D0,56165516844.D0,-4281164477.D0/
      DATA G1/8153167962181.D0,-30282470168220.D0,101786970400854.D0
     .,-244754720525140.D0,430901143447815.D0,-566360031923352.D0
     .,560707147204260.D0,-417423287061288.D0,230582338000515.D0
     .,-91816356051340.D0,24948866958246.D0,-4145485036740.D0
     .,318065528209.D0/
      DATA BLC/13.,-78.,286.,-715.,1287.,-1716.,1716.,-1287.,715.,-286.,
     .78.,-13.,1./
      DATA ER/2440800.5D0
     .,0.6388226331871791D+0,-.7234663691390341D+0,-.3137200168890611D+0
     .,-.3472826956196578D+0,-.2538198621886492D+0,-.1002526349881503D+0
     .,-.2729258630496380D+0,-.6191874090265185D+0,-.2617889156852868D+0
     .,-.1043274450714903D+1,0.1149609567068845D+1,0.5556917256839074D+0
     .,-.4241139166447770D+1,-.3149106786784197D+1,-.1247304720843395D+1
     .,0.6459610758313189D+1,0.6101130124894796D+1,0.2243382247325508D+1
     .,-.1814311868095884D+2,-.2460531337350525D+1,-.8217265400205982D+0
     .,-.1530833404391927D+2,-.2435427381300139D+2,-.9591188954641517D+1
     .,-.3031628535964815D+2,-.1792426075360332D+1,0.8620253359830489D+1
     .,0.2789435760038078D+1,0.8315178079235009D+0,-.1770597331965483D+0
     .,0.2838518810188359D+1,-.1731274670350053D+1,0.1349237040557463D+0
     .,-.2229853716668692D+1,-.3544645600066997D+0,0.1505957250879900D+0
     .,0.2055647541706483D+1,0.4646808549219674D+0,0.571194552728441D-2/
      DATA EV/
     . 0.1308702020817010D-1,0.9880582391240620D-2,0.4284464891607220D-2
     .,0.1164351062650545D-1,-.1796891116169349D-1,-.1081911400871376D-1
     .,0.1860072818427010D-1,-.6590509148918250D-2,-.4144660354448460D-2
     .,-.1028575666404713D-1,-.7076375256663700D-2,-.2972120784045060D-2
     .,0.4626311466814210D-2,-.5059626431851830D-2,-.2283653988263900D-2
     .,-.4252234340822560D-2,0.3559265916616560D-2,0.1655565829382040D-2
     .,0.5261924464533600D-3,-.3735782236086040D-2,-.1644317466285200D-2
     .,0.2693302696796940D-2,-.1429477504289330D-2,-.6534280574085100D-3
     .,0.3986115842506400D-3,-.3146654475399760D-2,-.1114235221665060D-2
     .,-.2763866219486554D-2,0.8274414179803879D-2,0.4463057993782906D-2
     .,0.3424650232882695D-2,0.7498263250144949D-2,-.1733648725951915D-2
     .,0.2299023705728994D-2,-.1052177497695365D-1,-.4500463020367574D-2
     .,-.4717858754284870D-2,0.1190825966047914D-1,0.242366396205390D-2/
      DATA CLTY,CLTZ/-.3978811766D0,.9174369566D0/
!     DATA K/.01720209895D0/,DIFMAX/0.D0/,LAT/40.72D0/,LON/-73.87D0/
      DATA DIFMAX/0.D0/,LAT/40.72D0/,LON/-73.87D0/

      IF(INIT.NE.0)GO TO 1000

!     Read in initial conditions
      read(11,11)c_filename
11    format(a14)
      if(c_filename .ne. c_fileref)then
          write(6,*)' Opening planet file for initial conditions: '
     1                ,c_filename
          open(12,file='/home/fab/albers/ast/planets/'
     1                 //c_filename,status='old')
          READ(12,12)ER,EV
          write(6,12)ER,EV
12        FORMAT(//1X,D25.10/26(/1X,3D26.19))
          close(12)
      else
          write(6,*)' Using default initial conditions'
      endif ! We are reading from an external file (not using default data)
      read(11,*)TOPOF,LAT,LON
      write(6,*)' TOPOF,lat,lon',topof,lat,lon

      H=1.D0
      ARYCEN=aw/2 + 1
      DEGREE=13
      DEGM1=DEGREE-1
      HSQ=H*H
      ARYHW=ARYCEN-1
      IARG2=DEGREE/2
      NRECEN=ARYCEN+IARG2
      IARG1=IARG2+1
      IARG3=ARYCEN-IARG1
      ARG3=FLOAT(IARG3)*H
      NUMBER=ARYCEN
      NUMB2=NUMBER+1
      NUMB3=NUMBER+DEGM1
      NUMB1=NUMBER-1
      KSQ=K*K
      NPINMS=MIN0(NPLAN-1,NMASS)
      PHI=LAT*RPD

      DO I=1,13
          A(I)=A1(I)/2615348736000.D0*HSQ
          AP(I)=AP1(I)/2615348736000.D0*HSQ
          B(I)=B1(I)/1743565824000.D0*HSQ
          G(I)=G1(I)/5230697472000.D0*HSQ
      enddo

      IDIR=0
C
C INITIALIZE THE INTEGRATION
      TCENT=ER(1)
      ISCR=1

      DO J=1,NPLAN
          X(J,NUMBER)=ER(ISCR+1)
          Y(J,NUMBER)=ER(ISCR+2)
          Z(J,NUMBER)=ER(ISCR+3)

          XD(J)=EV(ISCR)
          YD(J)=EV(ISCR+1)
          ZD(J)=EV(ISCR+2)

          ISCR=ISCR+3

      enddo ! J

      CALL ACCEL(NUMBER,NPLAN,NREL,NMASS,NPINMS,TCENT)
      GO TO 100
C
C EXTRAPOLATE THE ACCELERATIONS BACKWARDS
25    IDIR=1
      DIFMX1=DIFMAX
      DIFMAX=0.
      DO NP=1,DEGM1
          N=NUMBER-NP
          N1=N+1

          DO J=1,NPLAN
              XDDOLD=XDD(J,N)
              YDDOLD=YDD(J,N)
              ZDDOLD=ZDD(J,N)
              XDD(J,N)=BLC(1)*XDD(J,N1)
              YDD(J,N)=BLC(1)*YDD(J,N1)
              ZDD(J,N)=BLC(1)*ZDD(J,N1)

              DO I=2,DEGREE
                  XDD(J,N)=XDD(J,N)+BLC(I)*XDD(J,N+I)
                  YDD(J,N)=YDD(J,N)+BLC(I)*YDD(J,N+I)
                  ZDD(J,N)=ZDD(J,N)+BLC(I)*ZDD(J,N+I)

              enddo ! I

              IF(INIT1.GE.5)THEN
                  DIFMAX=
     1     DMAX1(DIFMAX,ABS(XDD(J,N)-XDDOLD),ABS(YDD(J,N)-YDDOLD),
     .             ABS(ZDD(J,N)-ZDDOLD))
              endif

          enddo ! J

      enddo ! NP

      WRITE(6,9101)DIFMAX,DIFMX1
      IF(INIT1.GT.25)then
          write(6,*)' Maximum iterations exceeded'
          STOP
      ENDIF
      IF(DIFMAX.GT.DIFMX1.AND.DIFMAX.LT.1D-9.AND.INIT1.GT.5)GO TO 990
100   INIT1=INIT1+1
C
C CALCULATE X(-1), Y(-1), Z(-1)
      DO J=1,NPLAN
          X(J,NUMB1)=X(J,NUMBER)-H*EV(J*3-2)
          Y(J,NUMB1)=Y(J,NUMBER)-H*EV(J*3-1)
          Z(J,NUMB1)=Z(J,NUMBER)-H*EV(J*3)

          DO N=1,DEGREE
              IARG=NUMBER+(N-1)*IDIR
              X(J,NUMB1)=X(J,NUMB1)+G(N)*XDD(J,IARG)
              Y(J,NUMB1)=Y(J,NUMB1)+G(N)*YDD(J,IARG)
              Z(J,NUMB1)=Z(J,NUMB1)+G(N)*ZDD(J,IARG)

          enddo ! N
      enddo ! J
C
C INTEGRATE FORWARD 12 STEPS
      DO N=NUMB2,NUMB3
          IF(INIT1.NE.1)THEN

              DO M=1,DEGREE
                  S(M)=N-M+1
              enddo ! M

          else

              DO M=1,DEGREE
                  S(M)=NUMBER
              enddo ! M

          endif ! init1 .ne. 1


          NM1=N-1
          NM2=N-2

          DO J=1,NPLAN
              X(J,N)=X(J,NM1)+X(J,NM1)-X(J,NM2)
              Y(J,N)=Y(J,NM1)+Y(J,NM1)-Y(J,NM2)
              Z(J,N)=Z(J,NM1)+Z(J,NM1)-Z(J,NM2)

              DO M=1,DEGREE
                  X(J,N)=X(J,N)+A(M)*XDD(J,S(M))
                  Y(J,N)=Y(J,N)+A(M)*YDD(J,S(M))
                  Z(J,N)=Z(J,N)+A(M)*ZDD(J,S(M))

              enddo ! M
          enddo ! J
      enddo ! N

      DO N=NUMB2,NUMB3

          IF(NREL.NE.0)THEN
              DO J=1,NREL
                  XD(J)=X(J,N-1)-X(J,N-2)
                  YD(J)=Y(J,N-1)-Y(J,N-2)
                  ZD(J)=Z(J,N-1)-Z(J,N-2)

                  DO M=1,DEGREE
                      NMM=N-M
                      IF(INIT1.EQ.1)NMM=NUMBER
                      XD(J)=XD(J)+B(M)*XDD(J,NMM)
                      YD(J)=YD(J)+B(M)*YDD(J,NMM)
                      ZD(J)=ZD(J)+B(M)*ZDD(J,NMM)

                  enddo ! M

                  XD(J)=XD(J)/H
                  YD(J)=YD(J)/H
                  ZD(J)=ZD(J)/H

              enddo ! J
          endif ! NREL .NE. 0

          T1=ER(1)+FLOAT(N-ARYCEN)*H
610       CALL ACCEL(N,NPLAN,NREL,NMASS,NPINMS,T1)

      enddo ! N

      GO TO 25

990   INIT=1
      MAXPNT=NUMB3-ARYCEN
      MINPNT=NUMBER-ARYCEN
C
C DETERMINE WHETHER TO INTEGRATE OR INTERPOLATE
1000  RSTEPS=(T-TCENT)/H
      ISTEPS=INT(RSTEPS)

      IF(ISTEPS.eq.0)then
          IF(RSTEPS.ne.0.D0)then
              IDPL=INT(RSTEPS/ABS(RSTEPS))

          else
              IDPL=-1

          endif
      else
          IDPL=ISTEPS/IABS(ISTEPS)

      endif ! ISTEPS .eq. 0

      P=ABS(RSTEPS-FLOAT(ISTEPS))

      IF(PLANET.EQ.0)then
          MAXINT=ISTEPS
          IF(IDPL.EQ.-1)MAXINT=MAXINT+DEGREE
          MININT=MAXINT-DEGREE
          IDR1=IDPL

      else
          MAXINT=ISTEPS+(5+IDPL)/2
          MININT=ISTEPS-(5-IDPL)/2

      endif ! PLANET .EQ. 0


1050  IF(MAXINT.gt.MAXPNT)then
          MAXPNT=MAXPNT+1
          IPOINT=MAXPNT
          IDIR=1

      else
          IF(MININT.lt.MINPNT)then
              MINPNT=MINPNT-1
              IPOINT=MINPNT
              IDIR=-1
          else
              goto4000
          endif

      endif

      IF(IABS(IPOINT).GT.ARYHW)then
C
C REARRANGE ARRAYS
          N1=ARYCEN*(1-IDIR)
          N2=ARYCEN-IARG1*IDIR

          DO N=1,NRECEN
              N1=N1+IDIR
              N2=N2+IDIR
              DO J=1,NPLAN
                  X(J,N1)=X(J,N2)
                  Y(J,N1)=Y(J,N2)
                  Z(J,N1)=Z(J,N2)
                  XDD(J,N1)=XDD(J,N2)
                  YDD(J,N1)=YDD(J,N2)
                  ZDD(J,N1)=ZDD(J,N2)
              enddo
          enddo

          IF(IDIR.LE.0)THEN
              IPOINT=-IARG1
              TCENT=TCENT-ARG3
              ISTEPS=ISTEPS+IARG3
              MAXINT=MAXINT+IARG3
              MININT=MININT+IARG3
              MAXPNT=ARYHW
              MINPNT=IPOINT

          else
              IPOINT=IARG1
              TCENT=TCENT+ARG3
              ISTEPS=ISTEPS-IARG3
              MAXINT=MAXINT-IARG3
              MININT=MININT-IARG3
              MAXPNT=IPOINT
              MINPNT=-ARYHW

          endif ! idir .le. 0

      endif ! Rearrange arrays

C
C INTEGRATE FORWARD OR BACKWARD ONE STEP
      S(1)=IPOINT+ARYCEN
      S(2)=S(1)-IDIR
      S(3)=S(2)-IDIR
      S(4)=S(3)-IDIR
      S(5)=S(4)-IDIR
      S(6)=S(5)-IDIR
      S(7)=S(6)-IDIR
      S(8)=S(7)-IDIR
      S(9)=S(8)-IDIR
      S(10)=S(9)-IDIR
      S(11)=S(10)-IDIR
      S(12)=S(11)-IDIR
      S(13)=S(12)-IDIR
      S(14)=S(13)-IDIR
      DO 2500 I=1,NPLAN
      X(I,S(1))=X(I,S(2))+X(I,S(2))-X(I,S(3))+AP(1)*XDD(I,S(2))
     .+AP(2)*XDD(I,S(3))+AP(3)*XDD(I,S(4))+AP(4)*XDD(I,S(5))
     .+AP(5)*XDD(I,S(6))+AP(6)*XDD(I,S(7))+AP(7)*XDD(I,S(8))
     .+AP(8)*XDD(I,S(9))+AP(9)*XDD(I,S(10))+AP(10)*XDD(I,S(11))
     .+AP(11)*XDD(I,S(12))+AP(12)*XDD(I,S(13))+AP(13)*XDD(I,S(14))
      Y(I,S(1))=Y(I,S(2))+Y(I,S(2))-Y(I,S(3))+AP(1)*YDD(I,S(2))
     .+AP(2)*YDD(I,S(3))+AP(3)*YDD(I,S(4))+AP(4)*YDD(I,S(5))
     .+AP(5)*YDD(I,S(6))+AP(6)*YDD(I,S(7))+AP(7)*YDD(I,S(8))
     .+AP(8)*YDD(I,S(9))+AP(9)*YDD(I,S(10))+AP(10)*YDD(I,S(11))
     .+AP(11)*YDD(I,S(12))+AP(12)*YDD(I,S(13))+AP(13)*YDD(I,S(14))
      Z(I,S(1))=Z(I,S(2))+Z(I,S(2))-Z(I,S(3))+AP(1)*ZDD(I,S(2))
     .+AP(2)*ZDD(I,S(3))+AP(3)*ZDD(I,S(4))+AP(4)*ZDD(I,S(5))
     .+AP(5)*ZDD(I,S(6))+AP(6)*ZDD(I,S(7))+AP(7)*ZDD(I,S(8))
     .+AP(8)*ZDD(I,S(9))+AP(9)*ZDD(I,S(10))+AP(10)*ZDD(I,S(11))
     .+AP(11)*ZDD(I,S(12))+AP(12)*ZDD(I,S(13))+AP(13)*ZDD(I,S(14))
      IF(NREL.LT.I)GO TO 2500
      XD(I)=(X(I,S(2))-X(I,S(3))+B(1)*XDD(I,S(2))+B(2)*XDD(I,S(3))
     .+B(3)*XDD(I,S(4))+B(4)*XDD(I,S(5))+B(5)*XDD(I,S(6))
     .+B(6)*XDD(I,S(7))+B(7)*XDD(I,S(8))+B(8)*XDD(I,S(9))
     .+B(9)*XDD(I,S(10))+B(10)*XDD(I,S(11))+B(11)*XDD(I,S(12))
     .+B(12)*XDD(I,S(13))+B(13)*XDD(I,S(14)))/H
      YD(I)=(Y(I,S(2))-Y(I,S(3))+B(1)*YDD(I,S(2))+B(2)*YDD(I,S(3))
     .+B(3)*YDD(I,S(4))+B(4)*YDD(I,S(5))+B(5)*YDD(I,S(6))
     .+B(6)*YDD(I,S(7))+B(7)*YDD(I,S(8))+B(8)*YDD(I,S(9))
     .+B(9)*YDD(I,S(10))+B(10)*YDD(I,S(11))+B(11)*YDD(I,S(12))
     .+B(12)*YDD(I,S(13))+B(13)*YDD(I,S(14)))/H
      ZD(I)=(Z(I,S(2))-Z(I,S(3))+B(1)*ZDD(I,S(2))+B(2)*ZDD(I,S(3))
     .+B(3)*ZDD(I,S(4))+B(4)*ZDD(I,S(5))+B(5)*ZDD(I,S(6))
     .+B(6)*ZDD(I,S(7))+B(7)*ZDD(I,S(8))+B(8)*ZDD(I,S(9))
     .+B(9)*ZDD(I,S(10))+B(10)*ZDD(I,S(11))+B(11)*ZDD(I,S(12))
     .+B(12)*ZDD(I,S(13))+B(13)*ZDD(I,S(14)))/H
C     PRINT*,I,S(1),X(I,S(1))
2500  CONTINUE

      T1=TCENT+H*FLOAT(IPOINT)
C     PRINT*,NPLAN,S(1),X(1,S(1)),X(13,S(1))
      CALL ACCEL(S(1),NPLAN,NREL,NMASS,NPINMS,T1)
      GO TO 1050

4000  APOINT=ISTEPS+ARYCEN
      IF(PLANET.EQ.0)GO TO 9000
      PSQ=P*P
      PCUBD=PSQ*P
      APNTM1=APOINT-IDPL
      APNTP1=APOINT+IDPL
      APNTP2=APNTP1+IDPL
C
C INTERPOLATE TO FIND POSITION OF PLANET AT TIME T
      B2=.25*(PSQ-P)
      B3=(.5*P-1.5*PSQ+PCUBD)/6.D0
      DX01=X(PLANET,APNTP1)-X(PLANET,APOINT)
      DY01=Y(PLANET,APNTP1)-Y(PLANET,APOINT)
      DZ01=Z(PLANET,APNTP1)-Z(PLANET,APOINT)
      DX02=X(PLANET,APNTP1)-2.*X(PLANET,APOINT)+X(PLANET,APNTM1)
      DY02=Y(PLANET,APNTP1)-2.*Y(PLANET,APOINT)+Y(PLANET,APNTM1)
      DZ02=Z(PLANET,APNTP1)-2.*Z(PLANET,APOINT)+Z(PLANET,APNTM1)
      DX12=X(PLANET,APNTP2)-2.*X(PLANET,APNTP1)+X(PLANET,APOINT)
      DY12=Y(PLANET,APNTP2)-2.*Y(PLANET,APNTP1)+Y(PLANET,APOINT)
      DZ12=Z(PLANET,APNTP2)-2.*Z(PLANET,APNTP1)+Z(PLANET,APOINT)
      DX03=DX12-DX02
      DY03=DY12-DY02
      DZ03=DZ12-DZ02
      R1=X(PLANET,APOINT)+P*DX01+B2*(DX02+DX12)+B3*DX03
      R2=Y(PLANET,APOINT)+P*DY01+B2*(DY02+DY12)+B3*DY03
      R3=Z(PLANET,APOINT)+P*DZ01+B2*(DZ02+DZ12)+B3*DZ03
      IF(PLANET.NE.2)GO TO 4400
      APNTP3=APNTP2+IDPL
      APNTM2=APNTM1-IDPL
      B4=P*(PCUBD-PSQ-PSQ-P+2.D0)/48.D0
      DDX4=X(PLANET,APNTP3)-2.*X(PLANET,APNTP2)+X(PLANET,APNTP1)-DX12
     .    +X(PLANET,APOINT)-2.*X(PLANET,APNTM1)+X(PLANET,APNTM2)-DX02
      DDY4=Y(PLANET,APNTP3)-2.*Y(PLANET,APNTP2)+Y(PLANET,APNTP1)-DY12
     .    +Y(PLANET,APOINT)-2.*Y(PLANET,APNTM1)+Y(PLANET,APNTM2)-DY02
      DDZ4=Z(PLANET,APNTP3)-2.*Z(PLANET,APNTP2)+Z(PLANET,APNTP1)-DZ12
     .    +Z(PLANET,APOINT)-2.*Z(PLANET,APNTM1)+Z(PLANET,APNTM2)-DZ02
      R1=R1+B4*DDX4
      R2=R2+B4*DDY4
      R3=R3+B4*DDZ4
4400  CONTINUE


      IF(PLANET.EQ.1.AND.MOON.EQ.1)THEN

!         ADD LUNAR PERTURBATION TO EARTH'S POSITION

!         Get Earth Moon vector for equinox of date
          CALL MOON_BRWN(T,MX,MY,MZ)

!         Convert vector to 1950 coordinates
          CALL PRECES(T,T1950,MX,MY,MZ,1)

          BX2 = -MX/EMRATP1
          BY2 = -MY/EMRATP1
          BZ2 = -MZ/EMRATP1

          R1=R1+BX2
          R2=R2+BY2
          R3=R3+BZ2


          IF(TOPOF.EQ.1)then
!             CONVERT TO TOPOCENTRIC POSITION
              ut1 = t - 48./86400.
              call TOPO(PHI,LON,ut1,TX,TY,TZ)
!             Convert Topocentric Position to 1950 Coordinates
              CALL PRECES(t,T1950,TX,TY,TZ,1)
              R1=R1+TX
              R2=R2+TY
              R3=R3+TZ
          endif

      endif ! Planet .eq. 1 .and. moon .eq. 1


!     call rotate(R1,R2,R3,-23.45)

      RETURN


C
C WRITE OUT POSITIONS AND VELOCITIES FOR A TIME NEAR TIME T
9000  TT=TCENT+FLOAT(ISTEPS)*H
      write(c_filename,9001)c_file1,int(TT - 0.5)
9001  format(a7,i7)

      open(2,file='/home/fab/albers/ast/planets/'
     1            //c_filename//'_new',status='new')

      WRITE(6,9100)TT
      WRITE(2,9100)TT
9100  FORMAT(//1X,F25.10/)
9101  FORMAT(1X,3D26.19)

      DO J=1,NPLAN
          WRITE(2,9101)x(j,apoint),y(j,apoint),z(j,apoint)
          WRITE(6,9101)x(j,apoint),y(j,apoint),z(j,apoint)

      enddo ! J

      DO J=1,NPLAN
          XVEL=X(J,APOINT-IDR1)-X(J,APOINT-IDR1-IDR1)
          YVEL=Y(J,APOINT-IDR1)-Y(J,APOINT-IDR1-IDR1)
          ZVEL=Z(J,APOINT-IDR1)-Z(J,APOINT-IDR1-IDR1)

          DO I=1,DEGREE
              XVEL=XVEL+B(I)*XDD(J,APOINT-I*IDR1)
              YVEL=YVEL+B(I)*YDD(J,APOINT-I*IDR1)
              ZVEL=ZVEL+B(I)*ZDD(J,APOINT-I*IDR1)
          enddo ! I


          XVEL=XVEL/H*FLOAT(IDR1)
          YVEL=YVEL/H*FLOAT(IDR1)
          ZVEL=ZVEL/H*FLOAT(IDR1)

          WRITE(2,9101)XVEL,YVEL,ZVEL
          WRITE(6,9101)XVEL,YVEL,ZVEL
      enddo ! J

      DO I=2,NPLAN
          CALL MAGNitude(I,0,X(1,APOINT),Y(1,APOINT),Z(1,APOINT)
     .    ,X(I,APOINT),Y(I,APOINT),Z(I,APOINT),value,dum)
          WRITE(6,9600)I,VALUE
9600      FORMAT(1X,I4,F15.8)
      enddo ! I

      close(2)


9999  RETURN
      END

