

      SUBROUTINE ACCEL(ISCRPT,NPLAN,NREL,NMASS,NPINMS,T)
      IMPLICIT REAL*8 (A,B,C,D,E,F,G,H,O,P,Q,R,S,T,U,V,W,X,Y,Z)
      DOUBLE PRECISION MASS,KSQ,KSDM15,LMA
      DIMENSION MASS(13),X(13),Y(13),Z(13),XDD(13),YDD(13),ZDD(13)
     .,UX(13),UY(13),UZ(13)
      REAL*8 MX,MY,MZ
      COMMON XP(13,59),YP(13,59),ZP(13,59),XDDP(13,59),YDDP(13,59)
     .,ZDDP(13,59),XD(13),YD(13),ZD(13),KSQ
      DATA MASS/3.040432924D-6,1.660136795D-7,2.447839598D-6
     .         ,3.227149362D-7,9.547861040D-4,2.858367872D-4
     .         ,4.372731645D-5,5.177591384D-5,3.333333333D-7
     .         ,5.9       D-10,1.3       D-10,1.2       D-10,0.D0/
      DATA FOURM/3.948251502D-8/,CM2/3.335661215D-5/
      DATA EMRAT/81.3007D0/,F1/.9878494351D0/,F2/.0121505649/
      DATA CLTY,CLTZ/-.3978811766D0,.9174369566D0/
      DATA EMRATP1/82.3007D0/,T1950/2433282.423D0/

      DO I=1,NPLAN
          X(I)=XP(I,ISCRPT)
          Y(I)=YP(I,ISCRPT)
          Z(I)=ZP(I,ISCRPT)
      enddo
C
C CALCULATE MUTUAL ACCELERATIONS BETWEEN SUN AND EACH PLANET
      DO I=1,NPLAN
          RLCM1=1.D0/SQRT(X(I)*X(I)+Y(I)*Y(I)+Z(I)*Z(I))
          RLCM3=RLCM1*RLCM1*RLCM1

          IF(I.eq.1)then

!             CALCULATE NEWTONIAN SOLAR ACCELERATION OF EARTH MOON SYSTEM

!             Get Earth Moon vector for equinox of date
              CALL MOON_BRWN(T,MX,MY,MZ)

!             Convert vector to 1950 coordinates
              CALL PRECES(T,T1950,MX,MY,MZ,1)

              PX = -MX/EMRATP1
              PY = -MY/EMRATP1
              PZ = -MZ/EMRATP1

              XE=X(1)+PX
              YE=Y(1)+PY
              ZE=Z(1)+PZ

              XM=X(1)-PX*EMRAT
              YM=Y(1)-PY*EMRAT
              ZM=Z(1)-PZ*EMRAT

              ARGE=-KSQ/SQRT(XE*XE+YE*YE+ZE*ZE)**3
              ARGM=-KSQ/SQRT(XM*XM+YM*YM+ZM*ZM)**3

              XDD(I)=ARGE*XE*F1+ARGM*XM*F2
              YDD(I)=ARGE*YE*F1+ARGM*YM*F2
              ZDD(I)=ARGE*ZE*F1+ARGM*ZM*F2

          else
              ARG=-KSQ*RLCM3
              XDD(I)=ARG*X(I)
              YDD(I)=ARG*Y(I)
              ZDD(I)=ARG*Z(I)

          endif ! I .eq. 1

30        IF(NREL.GE.I)THEN

!             ADD PERTURBATIVE ACCELERATION DUE TO GENERAL RELATIVITY
              RDRD=X(I)*XD(I)+Y(I)*YD(I)+Z(I)*ZD(I)
              VSQ=XD(I)*XD(I)+YD(I)*YD(I)+ZD(I)*ZD(I)
              A=FOURM*RLCM1-VSQ*CM2
              B=FOURM*RDRD*RLCM3
              XDD(I)=XDD(I)*(1.D0-A)+B*XD(I)
              YDD(I)=YDD(I)*(1.D0-A)+B*YD(I)
              ZDD(I)=ZDD(I)*(1.D0-A)+B*ZD(I)

          endif

          UX(I)=XDD(I)*MASS(I)
          UY(I)=YDD(I)*MASS(I)
          UZ(I)=ZDD(I)*MASS(I)

      enddo ! N = 1,nplan

C
C ADD IN MUTUAL ACCELERATIONS BETWEEN PAIRS OF PLANETS
      DO I=1,NPINMS
          I1=I+1
          DO J=I1,NPLAN
              XDELT=X(J)-X(I)
              YDELT=Y(J)-Y(I)
              ZDELT=Z(J)-Z(I)
              KSDM15=KSQ/SQRT(XDELT*XDELT+YDELT*YDELT+ZDELT*ZDELT)**3
              ARG=KSDM15*MASS(J)
              XDD(I)=XDD(I)+ARG*XDELT
              YDD(I)=YDD(I)+ARG*YDELT
              ZDD(I)=ZDD(I)+ARG*ZDELT
              ARG=-KSDM15*MASS(I)
              XDD(J)=XDD(J)+ARG*XDELT
              YDD(J)=YDD(J)+ARG*YDELT
              ZDD(J)=ZDD(J)+ARG*ZDELT

          enddo

      enddo
C
C CALCULATE TOTAL ACCELERATION OF SUN
      SUMUX=UX(1)
      SUMUY=UY(1)
      SUMUZ=UZ(1)

      DO I=2,NMASS
          SUMUX=SUMUX+UX(I)
          SUMUY=SUMUY+UY(I)
          SUMUZ=SUMUZ+UZ(I)

      enddo
C
C COMPUTE ACC. RELATIVE TO SUN = TOTAL ACC. - ACC. OF SUN
      DO I=1,NPLAN
          XDD(I)=XDD(I)+SUMUX
          YDD(I)=YDD(I)+SUMUY
          ZDD(I)=ZDD(I)+SUMUZ
      enddo

C
C PLACE ACCELERATIONS IN TWO DIMENSIONAL ARRAY
      DO I=1,NPLAN
          XDDP(I,ISCRPT)=XDD(I)
          YDDP(I,ISCRPT)=YDD(I)
          ZDDP(I,ISCRPT)=ZDD(I)
      enddo

      RETURN
      END
