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

        subroutine get_tw_approx_2d(t_k,td_k,p_pa,ni,nj,tw_k)

!       Steve Albers 1991

!       This routine is fast but only accurate near 0 degrees C (273K)

        real*4 t_k(ni,nj)     ! Input
        real*4 td_k(ni,nj)    ! Input
        real*4 p_pa(ni,nj)    ! Input
        real*4 tw_k(ni,nj)    ! Output

        const = alog(0.5)

        if(t_k(1,1) .eq. 0. .or. td_k(1,1) .eq. 0.)then
            write(6,*)' Bad input data to get_tw_approx_2d'
            return
        endif

        do j = 1,nj
        do i = 1,ni
            t_f =  t_k(i,j)  * 1.8 - 459.67
            td_f = td_k(i,j) * 1.8 - 459.67
            tmid_f = 0.5 * (t_f + td_f)
            depress50 = 13.4 + tmid_f * 0.1
            rh = 0.5 ** ((t_f - td_f)/depress50)
            start = 240. - (t_k(i,j) - 240.) / 6.
            ratio = (t_k(i,j) - start)/100. * (1.0 + 0.9 * (1.0 - sqrt(r
     1h)))
            dtw_dtd = ratio  ! (t = td)
            drh_dtd = const / depress50 ! (t = td)
            dtw_drh = dtw_dtd / drh_dtd
            tw_f = t_f - (rh - 1.0) * dtw_drh
            tw_k(i,j) = (tw_f + 459.67) / 1.8
        enddo ! i
        enddo ! j

        return
        end


        subroutine get_tw_2d(t_k,td_k,p_pa,ni,nj,tw_k)

!       Steve Albers 1991
!       WARNING: This routine may not work because it calls tw_fast

        real*4 t_k(ni,nj)     ! Input
        real*4 td_k(ni,nj)    ! Input
        real*4 p_pa(ni,nj)    ! Input
        real*4 tw_k(ni,nj)    ! Output

        do j = 1,nj
        do i = 1,ni
            tw_k(i,j) = tw_fast(t_k(i,j)-273.15,td_k(i,j)-273.15,p_pa(i,
     1j)*.01)
     1                             + 273.15
        enddo ! i
        enddo ! j

        return
        end

        subroutine get_tw_2d_orig(t_k,td_k,p_pa,ni,nj,tw_k)

!       Steve Albers 1991

        real*4 t_k(ni,nj)     ! Input
        real*4 td_k(ni,nj)    ! Input
        real*4 p_pa(ni,nj)    ! Input
        real*4 tw_k(ni,nj)    ! Output

        do j = 1,nj
        do i = 1,ni
            tw_k(i,j) = tw(t_k(i,j)-273.15,td_k(i,j)-273.15,p_pa(i,j)*.0
     11)
     1                             + 273.15
        enddo ! i
        enddo ! j

        return
        end

        FUNCTION TW_fast(T,TD,P)
C
!       WARNING: This routine may not work because it calls TSA_fast
C
C   THIS FUNCTION RETURNS THE WET-BULB TEMPERATURE TW (CELSIUS)
C   GIVEN THE TEMPERATURE T (CELSIUS), DEW POINT TD (CELSIUS)
C   AND PRESSURE P (MB).  SEE P.13 IN STIPANUK (1973), REFERENCED
C   ABOVE, FOR A DESCRIPTION OF THE TECHNIQUE.
C
C       BAKER,SCHLATTER 17-MAY-1982     Original version
C
C   DETERMINE THE MIXING RATIO LINE THRU TD AND P.
        AW = W_fast(TD,P)
C
C   DETERMINE THE DRY ADIABAT THRU T AND P.
        AO = O(T,P)
        PI = P
C
C   ITERATE TO LOCATE PRESSURE PI AT THE INTERSECTION OF THE TWO
C   CURVES .  PI HAS BEEN SET TO P FOR THE INITIAL GUESS.
        DO 4 I= 1,10
           X= .02*(TMR(AW,PI)-TDA(AO,PI))
           IF (ABS(X).LT.0.01) GO TO 5
 4         PI= PI*(2.**(X))
C   FIND THE TEMPERATURE ON THE DRY ADIABAT AO AT PRESSURE PI.
 5      TI= TDA(AO,PI)
C
C   THE INTERSECTION HAS BEEN LOCATED...NOW, FIND A SATURATION
C   ADIABAT THRU THIS POINT. FUNCTION OS RETURNS THE EQUIVALENT
C   POTENTIAL TEMPERATURE (K) OF A PARCEL SATURATED AT TEMPERATURE
C   TI AND PRESSURE PI.
        AOS= OS_fast(TI+273.15,PI)-273.15
C   FUNCTION TSA RETURNS THE WET-BULB TEMPERATURE (C) OF A PARCEL AT
C   PRESSURE P WHOSE EQUIVALENT POTENTIAL TEMPERATURE IS AOS.
        TW_fast = TSA_fast(AOS,P)
        RETURN
        END

        FUNCTION W_fast(T,P) ! Saturation mixing ratio wrt water
C
C   THIS FUNCTION RETURNS THE MIXING RATIO (GRAMS OF WATER VAPOR PER
C   KILOGRAM OF DRY AIR) GIVEN THE TEMPERATURE T (CELSIUS) AND PRESSURE
C   (MILLIBARS). THE FORMULA IS QUOTED IN MOST METEOROLOGICAL TEXTS.
C
C       BAKER,SCHLATTER 17-MAY-1982     Original version
C       Albers                 1992     modified for laps
C
        X= ESLO(T)
        W_fast= 622.*X/(P-X)
        RETURN
        END

        FUNCTION Wice_fast(T,P) ! Saturation mixing ratio wrt ice
C
C   THIS FUNCTION RETURNS THE MIXING RATIO (GRAMS OF WATER VAPOR PER
C   KILOGRAM OF DRY AIR) GIVEN THE TEMPERATURE T (CELSIUS) AND PRESSURE
C   (MILLIBARS). THE FORMULA IS QUOTED IN MOST METEOROLOGICAL TEXTS.
C
C       BAKER,SCHLATTER 17-MAY-1982     Original version
C       Albers                 1993     modified for laps
C
        X= ESILO(T)
        Wice_fast= 622.*X/(P-X)
        RETURN
        END

        FUNCTION OS_fast(TK,P)
C
C   THIS FUNCTION RETURNS THE EQUIVALENT POTENTIAL TEMPERATURE OS
C   (K) FOR A PARCEL OF AIR SATURATED AT TEMPERATURE T (K)
C   AND PRESSURE P (MILLIBARS).
C
C       BAKER,SCHLATTER 17-MAY-1982     Original version
C
        DATA B/2.6518986/
C   B IS AN EMPIRICAL CONSTANT APPROXIMATELY EQUAL TO THE LATENT HEAT
C   OF VAPORIZATION FOR WATER DIVIDED BY THE SPECIFIC HEAT AT CONSTANT
C   PRESSURE FOR DRY AIR.
        TC = TK - 273.15

!       From W routine
        X= ESLO(TC)
        W= 622.*X/(P-X)

        OS_fast= TK*((1000./P)**.286)*(EXP(B*W/TK))

        RETURN
        END



        subroutine laps_be(ni,nj,nk
     1  ,t_sfc_k,td_sfc_k,p_sfc_pa,t_3d_k,ht_3d_m,p_1d_pa,topo
     1                  ,pbe_2d,nbe_2d)

!       1991    Steve Albers
!       Returns PBE and NBE in Joules, Parcel is lifted from lowest level
!                                                            i.e. sfc

        include 'lapsparms.inc' ! nothing

        real*4 t_sfc_k(ni,nj)
        real*4 td_sfc_k(ni,nj)
        real*4 p_sfc_pa(ni,nj)
        real*4 topo(ni,nj)
        real*4 t_3d_k(ni,nj,nk)
        real*4 ht_3d_m(ni,nj,nk)
        real*4 p_1d_pa(nk)
        real*4 p_1d_mb(nk)                   ! Local
        real*4 pbe_2d(ni,nj)
        real*4 nbe_2d(ni,nj)

        COMMON/INDX/ P(70),T(70),TD(70),HT(70),PBECR(20,4),TDFCR(20,2)
     1  ,VEL(20)
     1  ,temdif(70),partem(70),pbe(70)
     #  ,DD85,FF85,DD50,FF50
        REAL LCL,LI,K_INDEX

        EXTERNAL LIB$INIT_TIMER,
     1           LIB$SHOW_TIMER,
     1           my_show_timer

        ISTAT = LIB$INIT_TIMER()
C
        do k = 1,nk
            p_1d_mb(k) = p_1d_pa(k) / 100.
        enddo ! k

        do i = 1,ni
c       write(6,*)' i = ',i
        do j = 1,nj

            do k = 1,nk
                if(p_1d_pa(k) .lt. p_sfc_pa(i,j))then ! First level above sfc
                    n_first_level = k
                    goto100
                endif
            enddo

100         p(1) = p_sfc_pa(i,j) / 100.       ! pa to mb
            t(1) = t_sfc_k(i,j) - 273.15      ! K to C
            td(1) = td_sfc_k(i,j) - 273.15    ! K to C
            ht(1) = topo(i,j)

            n = 1

            do k = n_first_level,nk
                n = n + 1
                P(N) = p_1d_mb(k)
                T(N) = t_3d_k(i,j,k) - 273.15 ! K to C
                TD(N)= t(n)
                HT(N)= ht_3d_m(i,j,k)
            enddo ! k

            nlevel = n

            CALL SINDX(NLEVEL,LI,SI,BLI,TT,SWEAT,HWB0,PLCL,LCL,CCL,TCONV
     1,IO
     #      ,ICP,ICT,K_INDEX,TMAX,PBENEG,PBEPOS,T500,PBLI,VELNEG
     #      ,WATER,IHOUR)

            pbe_2d(i,j) = PBEPOS
            nbe_2d(i,j) = PBENEG

            IF(i .eq. ni/2 .and. j .eq. nj/2)then
!           IF(i .eq. i/8*8 .and. j .eq. 12)then

                write(6,*)
     1  ' n    p         t         td         ht       tdif       be'
                do n = 1,nlevel
                    write(6,301)n,p(n),t(n),td(n),ht(n),temdif(n),pbe(n)
301                 format(1x,i2,f8.1,f10.2,f10.2,f10.0,f10.2,f10.1)
                enddo

                SWEAT=0.0
!               IF(IFLAG1*IFLAG2.EQ.0)SWEAT=0.0

                 WRITE(6,430)DD85,FF85,DD50,FF50,T500
 430            FORMAT(' 850MB WIND',2F5.0,'         500MB WIND',2F5.0
     #         ,'    500MB TEMP= ',F6.1)

!               IF(IO2.EQ.1)GOTO1000

                WRITE(6,60)LI,SI,BLI
 60             FORMAT(' LI= ',F5.1,20X,'SI= ',F5.1,15X,'BLI= ',F5.1)

                WRITE(6,61)TT,SWEAT,HWB0
 61             FORMAT(' TOTAL TOTALS=',F6.1,10X,'SWEAT= ',F7.1,10X,'W.
     1BULB ZERO='
     #           ,F5.1,' KFT AGL')

                ITCONV=INT(TCONV+0.5)

                WRITE(6,62)LCL*.003281,CCL,ITCONV
 62             FORMAT(' LCL= ',F5.1,' KFT AGL',11X,'CCL=',F5.1,' KFT AG
     1L',7X
     #          ,'CONVECTIVE TEMP=',I4,' DEG F')

                KK=NINT(K_INDEX)
                ITMAX=NINT(TMAX)

                WRITE(6,63)KK,ITMAX,WATER
 63             FORMAT(' K INDEX =',I4,16X,'TMAX =',I4,' F',12X,'PRECIP.
     1 WATER='
     #         ,F5.2,' IN.')
                WRITE(6,*)
C
                IF(ICT+ICP.GT.0.OR.IO.GE.2)WRITE(6,71)P(1),PLCL
 71             FORMAT(' ENERGY ANALYSIS - LIFTING A PARCEL FROM',F6.0,'
     1 MB'
     1                  ,'   LCL=',f7.1,' MB')
C
                DO 600 II=1,ICP
 600            WRITE(6,64)PBECR(II,1),PBECR(II,2),VEL(II),PBECR(II,4)
 64             FORMAT(' ENERGY AT',F6.0,'MB=',F10.2,' JOULES/KG       V
     1ELOCITY='
     #        ,F6.2,'M/S   HT=',F7.0,' M')
C
                WRITE(6,*)
                DO 700 II=1,ICT
 700            WRITE(6,5)TDFCR(II,1),TDFCR(II,2)
 5              FORMAT(5X,'ENVIRONMENTAL MINUS PARCEL TEMPERATURE AT',F6
     1.0
     1           ,'MB  =',F6.1)
C
                IF(BLI.LE.0.0)WRITE(6,65)PBENEG,VELNEG
                IF(BLI.LE.0.0)WRITE(6,66)PBEPOS
 65             FORMAT('     CAP STRENGTH =',F10.1,' JOULES/KG.'
     +            ,' VELOCITY NEEDED',F5.1,'M/S')
 66             FORMAT('     POSITIVE AREA=',F10.1,' JOULES/KG.')
                GOTO2000
C
C  ALTERNATIVE OUTPUTTING
 1000           CONTINUE
                WRITE(6,1001)STA,T500,VELNEG,PBEPOS
 1001           FORMAT(4X,A3/F5.1,F7.1,F7.1)

2000        endif

        enddo ! j
        enddo ! i
C
 9999   CONTINUE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                                 ISTAT = LIB$SHOW_TIMER(my_show_timer)

        RETURN
        END

!
        SUBROUTINE SINDX(NLEVEL,LI,SI,BLI,TT,SWEAT,HWB0,PLCL_PBE,LCL,CCL  
     1   ,TCONV,IO,ICP,ICT,K,TMAX,PBENEG,PBEPOS,TMAN50,PBLI,VELNEG
     #   ,WATER,IHOUR)

!       1991    Steve Albers
!       Only the PBENEG and PBEPOS active outputs in this version

        DIMENSION Q(70),W(70)! ,WB(70)
        COMMON/INDX/ P(70),T(70),TD(70),HT(70),PBECR(20,4),TDFCR(20,2)
     1              ,VEL(20),temdif(70),partem(70),pbe(70)
     #              ,DD85,FF85,DD50,FF50
        REAL LI,K,LCL_AGL,LCL_PBE
        ES(X)=6.1078+X*(.443652+X*(.014289+X*(2.65065E-4+X*
     1 (3.03124E-6+X*(2.034081E-8+X*(6.13682E-11))))))
!       TDEW(E)=237.7/((7.5/ALOG10(E/6.11))-1.)
        DATA EPSILN/.6220/,G/9.80665/
        DATA BLTHCK/50.0/
        DATA RPD/.0174532925063/
 1      format('    MIXED PARCEL IS:',F10.1,'MB',F15.3,'C       '
     #          ,'MIXING RATIO=',F9.7)

c       WRITE(6,15)
 15     FORMAT('          PRESSURE     TEMP.    DEW PT.     '
     1          ,'Q      WET BULB')

        IOUT=IO

 3      FORMAT(' LVL(',I2,')',3F10.1,F11.6,F10.2)

        THETAE=OE_FAST(T(1),TD(1),P(1)) + 273.15

        CALL LCL_fast(P(1),T(1),TD(1),LCL_AGL,TLCL_PBE,PLCL_PBE)

        LCL_PBE = LCL_AGL + HT(1)

        CALL POTBE(Q,NLEVEL,P(1),T(1),W(1),PLCL_PBE
     #   ,TLCL_PBE,LCL_PBE,0.,THETAE,ICP,ICT,IO,PBENEG,PBEPOS
     #   ,VELNEG)
C
        IF(IO.GE.2)WRITE(6,3)

        DO 600 I=1,ICP
            VEL(I)=SQRT(ABS(2.*PBECR(I,2)))
            IF(PBECR(I,2).LT.0.)VEL(I)=-VEL(I)
c           WRITE(6,4)PBECR(I,1),PBECR(I,2),VEL(I)
 4          format(' ENERGY AT',F6.0,'MB=',F15.2,'JOULES       VELOCITY=
     1'
     #         ,F6.1,'MS')
 600    CONTINUE

c       DO 700 I=1,ICT
c           WRITE(6,5)TDFCR(I,1),TDFCR(I,2)
c700    CONTINUE

 5      format(' ENVIRONMENTAL - PARCEL TEMPERATURE AT',F7.0
     1          ,'MB=',F8.2,' DEG C')

        return

        END
C
C
C
c        block data
c        COMMON/INDX/ P(70),T(70),TD(70),HT(70),PBECR(20,4),TDFCR(20,2)
c     1              ,VEL(20),temdif(70),partem(70),pbe(70)
c     #              ,DD85,FF85,DD50,FF50
c        DATA PBE/70*0./
c        end
C
C
C
C
        SUBROUTINE POTBE(Q,NLEVEL,PMEAN,TMEAN,WMEAN,PLCL
     #   ,TLCL,LCL,HBLAYR,BLTHTE,ICP,ICT,IO,PBENEG,pos_area_max
     #   ,VELNEG)

!       Steve Albers 1991

        COMMON/INDX/ P(70),T(70),TD(70),HT(70),PBECR(20,4),TDFCR(20,2),V
     1EL(20)
     1 ,temdif(70),partem(70),pbe(70)
     #   ,DD85,FF85,DD50,FF50
        REAL LCL,nbe_min
        DIMENSION Q(70),DELTAH(70)
        real*4 GAMMA
        parameter (GAMMA = .009760) ! Dry Adiabatic Lapse Rate Deg/m
        DATA G/9.80665/

        include 'lapsparms.inc' ! r_missing_data

        nbe_min=0.
        pos_area_max=0.

        partem(1) = t(1)
        TEMDIF(1) = 0.
        PBE(1) = 0.

c       WRITE(6,606)PMEAN,TMEAN,WMEAN,PLCL,TLCL,LCL
c 606   format(' PMEAN,TMEAN,WMEAN,PLCL,TLCL,LCL',2F10.2,F10.5,2F10.2
c     +   ,F5.1)

        DO 800 N=2,NLEVEL
            DELTAP=P(N-1)-P(N)
            n_low=N-1

            deltah(n) = HT(n) - HT(n-1)

            IF(PLCL .lt. P(N-1))then ! Lower level is below LCL
                IF(PLCL .lt. P(N))then ! Upper level is below LCL
c                   WRITE(6,*)' DRY CASE'
                    PARTEM(N)=PARTEM(N-1)-GAMMA*DELTAH(N)
                    TEMDIF(N)=PARTEM(N)-T(N)
                    PBE(N)=PBE(N-1)+G*(TEMDIF(N)+TEMDIF(N-1))/(T(N)
     #                          +T(N-1)+546.30)*DELTAH(N)

                else ! Upper level is above LCL
C                   BRACKETING CASE AROUND LCL - DRY ADIABATIC PART
c                   WRITE(6,*)' DRY ADIABATIC PART'

                    delta_ht_dry = lcl - HT(n-1)

c                   WRITE(6,307)TLCL
 307                format(' PARCEL TEMP AT LCL= ',F10.3)

                    CALL ITPLV(P,T,NLEVEL,PLCL,SNTLCL,IO)

                    t_dif_lcl=TLCL-SNTLCL

                    pbe_dry=G*(t_dif_lcl+TEMDIF(N-1))/(SNTLCL+T(N-1)+546
     1.30)
     +                  *delta_ht_dry

 951                format(' pos_area_max,PSI,NEG,NNG,PBECR',5F9.1,I3)

!d                  WRITE(6,777)N,P(N),TLCL,SNTLCL,t_dif_lcl
!d    #                    ,TEMDIF(N-1),delta_ht_dry,HT(N),PBE(N-1)+pbe_dry

C                   MOIST ADIABATIC PART
c                   WRITE(6,*)' MOIST ADIABATIC PART'
                    delta_ht_wet=DELTAH(N)-delta_ht_dry

                    partem(n) = tmlaps_fast(blthte,P(n))
                    TEMDIF(N)=PARTEM(N)-T(N)

                    pbe_wet = G*(TEMDIF(N)+t_dif_lcl)/(T(N)+SNTLCL+546.3
     10)
     #                           *delta_ht_wet

                    PBE(N)=PBE(N-1) + pbe_dry + pbe_wet

                endif ! Upper level below LCL (Dry or bracket)

            else ! Lower Level is above LCL
c               WRITE(6,*)' GETTING PARCEL TEMPERATURE FOR MOIST CASE'
                partem(n) = tmlaps_fast(blthte,P(n))

          !     Add Layer
                TEMDIF(N)=PARTEM(N)-T(N)
                PBE(N)=PBE(N-1)+G*(TEMDIF(N)+TEMDIF(N-1))/(T(N)
     #                          +T(N-1)+546.30)*DELTAH(N)

            endif


c           WRITE(6,777)N,P(N),PARTEM(N),T(N)
c    #                ,TEMDIF(N),TEMDIF(N-1),DELTAH(N),HT(N),PBE(N)
 777        format(' PBE DATA',I3,F6.1,6F8.2,F12.3)

 800    CONTINUE
C
C  FLAG SIGNIFICANT LEVELS FOR ENERGY
        ICP=0
        ICT=0
        N1 = 2
c       WRITE(6,382)N1,NLEVEL
 382    format(' N1,NLEVEL',2I3)
C
C  DETERMINE ENERGY EXTREMA - NEUTRAL BUOYANCY
        DO 1100 N=N1,NLEVEL
c           WRITE(6,940)N
 940        format(' LOOKING FOR NEUTRAL BUOYANCY - ENERGY EXTREMUM, LEV
     1EL',I3)

            IF((TEMDIF(N)*TEMDIF(N-1)).lt.0.)then
                ICT=ICT+1
                SLOPE=(TEMDIF(N)-TEMDIF(N-1))/ALOG((P(N)/P(N-1)))
                FRAC=TEMDIF(N-1)/(TEMDIF(N-1)-TEMDIF(N))
                TDFCR(ICT,1)=EXP(ALOG(P(N-1))-TEMDIF(N-1)/SLOPE)
                TDFCR(ICT,2)=0.
c               WRITE(6,948)ICT,N,TEMDIF(N-1),TEMDIF(N),SLOPE,FRAC
c    #           ,TDFCR(ICT,1)
 948            format(' NEUT BUOY',2I3,5F12.5)
                ICP=ICP+1
                TMID=(1.-FRAC)*T(N-1)+FRAC*T(N)
                PBECR(ICP,1)=TDFCR(ICT,1)
                PBECR(ICP,2)=PBE(N-1)+G*TEMDIF(N-1)/(T(N-1)+TMID+546.3)
     &             *DELTAH(N)*FRAC
C
                if(p(n) .ge. 500.)then
                    nbe_min=AMIN1(PBECR(ICP,2),nbe_min)
                endif

                pos_area=PBECR(ICP,2)-nbe_min
                pos_area_max=max(pos_area,pos_area_max)

!d              WRITE(6,951)pos_area_max,pos_area,PBENEG,nbe_min,PBECR(ICP,2),ICP

                PBECR(ICP,3)=PBECR(ICP,2)-PBECR(max(1,ICP-1),2)
                PBECR(ICP,4)=HT(N-1)+FRAC*DELTAH(N)
C
c               WRITE(6,949)ICP,PBECR(ICP,1),PBECR(ICP,2),PBECR(ICP,3)
c    &          ,PBECR(ICP,4),pos_area_max,nbe_min
 949            format(' PBECR',I3,4F10.3,2F10.2)

            endif

            goto1100
C
C  FIND BUOYANCY MAXIMA
c           IF(N.lt.NLEVEL)then
c               IF(TEMDIF(N) .gt. TEMDIF(N-1)
c       1                               .and. TEMDIF(N) .gt. TEMDIF(N+1))then
c                   IF(ABS(TEMDIF(N)).ge.ABS(TEMDIF(N-1)).or.
c     &                 TEMDIF(N)*TEMDIF(N-1).le.0)then
c                       ICT=ICT+1
c                       TDFCR(ICT,1)=P(N)
c                       TDFCR(ICT,2)=-TEMDIF(N)
c                       PBECR(1,3)=PBECR(1,2)
c                       WRITE(6,999)ICT,TEMDIF(N)
c 999                   format(' BUOYANCY MAX, ICT TEMDIF(N)',I3,F5.1)
c                    endif
c                endif
c            endif
C
C  LOOK FOR AN LFC (ZERO ENERGY)
c 1000      CONTINUE

C  FIND TEMPERATURE DIFFERENCE AS A FUNCTION OF HEIGHT THEN
C  FIND PBE AS A FUNCTION OF HEIGHT (ABOVE SIG LEVEL)
C  USE QUADRATIC FORMULA TO FIND ZERO PBE
c           TAVE=0.5*(T(N)+T(N-1)+546.3)
c           X=TEMDIF(N-1)
c           Y=(TEMDIF(N)-TEMDIF(N-1))/DELTAH(N)
c           TWOA=G*Y/TAVE
c           B=G*X/TAVE
c           C=PBE(N-1)
c           ARG=-B/TWOA
c           ARGDSC=B*B-2.*TWOA*C
C
c           IF(ARGDSC.ge.0)then
c               DISC=SQRT(ARGDSC)/TWOA
c               WRITE(6,352)TWOA,B,C,X,Y,ARG,DISC,ARGDSC
c 352           FORMAT(8F15.5)

c               DO ISIGN=-1,1,2
c                   HH=ARG+ABS(DISC)*ISIGN
c                   FRAC=HH/DELTAH(N)

c                   IF(ABS(FRAC-0.5).le.0.5)then ! Zero crossing detected
c                       ICP=ICP+1
c                       PBECR(ICP,1)=
c       1                   EXP((ALOG(P(N))-ALOG(P(N-1)))*FRAC+ALOG(P(N-1)))
c                       PBECR(ICP,2)=0.
c                       PBECR(ICP,3)=0.
c                       PBECR(ICP,4)=HT(N-1)+HH

c                       WRITE(6,1089)N,PBE(N),PBE(N-1),HH,PBECR(ICP,1)
c     &                          ,PBECR(ICP,4),HT(N),HH,DELTAH(N),FRAC
c 1089                  FORMAT(I3,9F11.3)
c                   endif

c               enddo

c           endif

1100    CONTINUE

c       WRITE(6,464)ICP,ICT,N1,NLEVEL
 464    format(' ICP,ICT,N1,NLEVEL',4I5)
C

!       Case when equlibrium level is above top of domain
        pos_area_max = max(pos_area_max,PBE(NLEVEL) - nbe_min)


!       At least one region of positive area in sounding
        IF(pos_area_max .gt. 0.0)then
            PBENEG=nbe_min
            VELNEG=SQRT(2.0*ABS(PBENEG))

        else ! Case when no positive area exists anywhere in sounding
            pos_area_max=-1.0
            PBENEG = r_missing_data ! -1e6
            VELNEG = 0.

        endif


c       WRITE(6,485)pos_area_max,PBENEG,VELNEG
 485    format(' pos_area_max',F10.1,' PBENEG',F10.1,' VELNEG',F10.1)
        RETURN
        END
C
!
C       1991    Steve Albers
C
        FUNCTION WTBLB(P,TC,W,IOUT)
        THETAE=THAEK(P,TC,W)
        CALL MSAD5(WTBLB_arg,P,THETAE,TC,20.,SLOPE,I1,I2,IA,IOUT)
        WTBLB = WTBLB_arg
        RETURN
        END
C
C
C
        SUBROUTINE BLAYR(P,T,Q,PMEAN,TMEAN,WMEAN,THKNES,NLEVEL,HH,
     +                    LOWEST,IO)

!       Steve Albers 1991

        REAL INTLOG
        DIMENSION P(70),T(70),Q(70)
        TVIRT(TT,QQ)=TT/(1.-QQ*.37803)
        THICK(P1,P2,TC,QQ)=ALOG(P1/P2)*TVIRT((TC+273.15),QQ)*.09604
        IF(NLEVEL.LT.2)GOTO9000
        HH=0.
        SUMWT=0.
        SUMT=0.
        SUMQ=0.
        IF(IO.GE.2)WRITE(6,1)LOWEST,P(LOWEST),T(LOWEST),Q(LOWEST),SUMWT
     +                        ,SUMWT,SUMT,SUMQ
C
        NLOW=LOWEST+1

        DO 100 I=NLOW,NLEVEL
        IF((P(LOWEST)-P(I))-THKNES)50,150,200
 50     WT=P(I)-P(I-1)
        SUMWT=SUMWT+WT
        AT=T(I-1)
        AQ=Q(I-1)
        ARGI=1./ALOG(P(I)/P(I-1))
        BT=(T(I)-T(I-1))*ARGI
        BQ=(Q(I)-Q(I-1))*ARGI
        INTLOG=P(I)*(ALOG(P(I))-1.)-P(I-1)*(ALOG(P(I-1))-1.)
        PDELT=P(I)-P(I-1)
        ARGT=(AT-(BT*ALOG(P(I-1))))*PDELT+BT*INTLOG
        ARGQ=(AQ-(BQ*ALOG(P(I-1))))*PDELT+BQ*INTLOG
        SUMT=SUMT+ARGT
        SUMQ=SUMQ+ARGQ
        ARGTP=ARGT/PDELT
        ARGQP=ARGQ/PDELT
        HH=HH+THICK(P(I-1),P(I),ARGTP,ARGQP)
        IF(IO.GE.2)WRITE(6,1)I,P(I),T(I),Q(I),WT,SUMWT,SUMT,SUMQ,BT,BQ
     &  ,ZERO,HH
 100    CONTINUE
        GOTO9000
C
 150    WT=P(I)-P(I-1)
        GOTO250
C
 200    CONTINUE
        WT=-(THKNES-P(LOWEST)+P(I-1))
 250    CONTINUE
        PBOUND=P(LOWEST)-THKNES
        SUMWT=SUMWT+WT
        AT=T(I-1)
        AQ=Q(I-1)
        ARGI=1./ALOG(P(I)/P(I-1))
        BT=(T(I)-T(I-1))*ARGI
        BQ=(Q(I)-Q(I-1))*ARGI
        INTLOG=PBOUND*(ALOG(PBOUND)-1.)-P(I-1)*(ALOG(P(I-1))-1.)
        PDELT=PBOUND-P(I-1)
        ARGT=(AT-(BT*ALOG(P(I-1))))*PDELT+BT*INTLOG
        ARGQ=(AQ-(BQ*ALOG(P(I-1))))*PDELT+BQ*INTLOG
        SUMT=SUMT+ARGT
        SUMQ=SUMQ+ARGQ
        HH=HH+THICK(P(I-1),PBOUND,ARGT/PDELT,ARGQ/PDELT)
        IF(IO.GE.2)WRITE(6,1)I,P(I),T(I),Q(I),WT,SUMWT,SUMT,SUMQ,BT,BQ
     +              ,INTLOG,HH
 1      FORMAT(1x,I2,F6.0,F6.1,F6.4,2F7.1,2F11.4,F11.6,F9.6,F9.3,F6.2)
C
        PMEAN=P(LOWEST)-.5*THKNES
        SUMWTI=1./SUMWT
        TMEAN=SUMT*SUMWTI
        QMEAN=SUMQ*SUMWTI
        WMEAN=QMEAN/(1.-QMEAN)
        GOTO9999
C
 9000   IF(IO.GE.2)WRITE(6,9001)
 9001   FORMAT(' NOT ENOUGH LEVELS TOO OBTAIN A BOUNDARY LAYER')
 9999   RETURN
        END
C
C
C
C
C
        SUBROUTINE LLCL(LCL,TLCL,PLCL,P,TC,W,IO) ! In Meters

!       Steve Albers 1991

        REAL LCL,KAPPA
        ES(X)=6.1078+X*(.443652+X*(.014289+X*(2.65065E-4+X*
     1 (3.03124E-6+X*(2.034081E-8+X*(6.13682E-11))))))
        DATA EPSILN/.62197/,GAMMAI/102.4596/ ! Lapse Rate m/deg
        TK=273.15+TC
        TLCL=TK
        KAPPA=.28613105*(1.-.23*W)
        SLOPE=5.
        E=(W*P)/(EPSILN+W)
        ITER=1
        GOTO140
 100    CONTINUE
        ITER=ITER+1
        IF(IO.GE.1)WRITE(6,25)TLCL,RES,SLOPE,DELTA
 25     FORMAT(' TLCL=',5F12.5)
        IF(ABS(RES)-.05)150,150,130
 130    TLCL=TLCL+DELTA
 140    RES=(ES(TLCL-273.15)/E)**KAPPA*TK-TLCL
        IF(ITER.NE.1)SLOPE=(RES-RESOLD)/DELTA
        DELTA=-RES/SLOPE
        RESOLD=RES
        GOTO100
 150    TLCL=TLCL-273.15
        PLCL=P*ES(TLCL)/E
        LCL=(TC-TLCL)*GAMMAI
        RETURN
        END
C
C
        SUBROUTINE LCL_fast(P,TC,TD,HLCL,TLCL,PLCL)

!       Steve Albers 1991

        REAL LCL,KAPPA
        ES(X)=6.1078+X*(.443652+X*(.014289+X*(2.65065E-4+X*
     1 (3.03124E-6+X*(2.034081E-8+X*(6.13682E-11))))))
        DATA EPSILN/.62197/,GAMMAI/102.4596/        ! Lapse Rate m/deg
        KAPPA=.28613105*(1.-.23*W)

        TLCL = TD - (0.212 + 0.001571 * TD - 0.000436 * TC) * (TC - TD)
        PLCL=P*ES(TLCL)/ES(TD)
        HLCL=(TC-TLCL)*GAMMAI

        RETURN
        END
C
C
C
C
C
        SUBROUTINE ITPLV(P,PARAM,NLEVEL,PINT,PARMAN,IO)

!       Steve Albers 1991

        DIMENSION P(70),PARAM(70)
        DO 100 N=1,NLEVEL
        IF(P(N)-PINT)300,200,100
 100    CONTINUE
C
 200    PARMAN=PARAM(N)
        GOTO999
C
 300    FRAC=ALOG(PINT/P(N-1))/ALOG(P(N)/P(N-1))
        PARMAN=PARAM(N-1)+FRAC*(PARAM(N)-PARAM(N-1))
        IF(IO.GE.2)WRITE(6,666)N,PINT,FRAC,PARMAN,P(N),P(N-1),PARAM(N)
     #  ,PARAM(N-1)
 666    FORMAT(' INTERPOLATING TO LEVEL',I3,F10.3,F8.5,F8.3,2F7.1,2F6.1)
C
 999    RETURN
        END
C
        SUBROUTINE NEWTN(X,XOLD,Y,YOLD,SLOPE,ITER,IO,FENCE)

!       Steve Albers 1991

        IF(ITER.GT.0)SLOPE=(Y-YOLD)/(X-XOLD)
        DELTA=-Y/SLOPE
        DELTA=AMAX1(-FENCE,AMIN1(DELTA,FENCE))
        XOLD=X
        X=X+DELTA
        ITER=ITER+1
        IF(IO.GE.2)WRITE(6,1)ITER,X,XOLD,Y,YOLD,SLOPE,DELTA
 1      FORMAT(' NEWTON ITER',I4,2F11.4,2F11.7,2E11.3)
        YOLD=Y
        RETURN
        END
C
C
        SUBROUTINE MSAD5
     ^  (TEMNEW,PRESNW,THETAE,TGUESS,SLOPEG,SLOPE,I1,I2,IA,IO)

!       Steve Albers 1991

        DIMENSION TEMPNW(4)
        I1=0
        I2=0
        IA=0
C
C  ESTIMATE PARCEL TEMP AT 500MB
        IF(PRESNW.NE.500.)GOTO725
        TEMPNW(1)=THETAE-(307.260+72.122*EXP((THETAE-382.635)*.0141672))
        TEMPNW(1)=TEMPNW(1)+0.65*EXP(-(.077*(TEMPNW(1)+27.))**2)
        TEMPNW(1)=TEMPNW(1)-0.50*EXP(-(.200*(TEMPNW(1)+4.0))**2)
        IF(THETAE.LT.254.90.OR.THETAE.GT.361.53)GOTO725
        TEMNEW=TEMPNW(1)
        GOTO900
C
 725    TEMPNW(1)=TGUESS
        SLOPE=SLOPEG
C  ENTER ITERATIVE LOOP IF THETAE IS OUT OF RANGE OF APPROXIMATE FORMULA
C  DETERMINE LATENT HEAT RELEASED ABOVE 500MB
 750    ITER=1
 760    EPSILN=THAE(TEMPNW(ITER),TEMPNW(ITER),PRESNW)-THETAE
        I1=I1+1
C  TEST FOR CONVERGENCE
        IF(ABS(EPSILN).LT..01)GOTO890
        IF(IO.GE.3.AND.ITER.EQ.1)
     ^  WRITE(6,666)ITER,TEMPNW(ITER),SLOPE,EPSILN,DELTA,EPSOLD,DELOLD
     ^  ,SLOPEG
        IF(ITER.GT.1)SLOPE=(EPSILN-EPSOLD)/DELOLD
        DELTA=-EPSILN/SLOPE
        ITER=ITER+1
        IF(IO.GE.3)
     ^  WRITE(6,666)ITER,TEMPNW(ITER),SLOPE,EPSILN,DELTA,EPSOLD,DELOLD
 666    FORMAT(' MSAD5',I2,2X,7F11.5)
        TEMPNW(ITER)=TEMPNW(ITER-1)+DELTA
        I2=I2+1
        IF(ITER.EQ.4)GOTO850
        EPSOLD=EPSILN
        DELOLD=DELTA
        GOTO760
C
C  USE AITKEN'S FORMULA TO ACCELERATE CONVERGENCE
 850    RATIO=DELTA/DELOLD
        RATABS=ABS(RATIO)
        RATIO=AMIN1(RATABS,.5)*RATIO/RATABS
        IF(IO.GE.3)WRITE(6,666)
     ^  ITER,TEMPNW(ITER),DELTA,DELOLD,RATIO
        TEMPNW(1)=TEMPNW(3)+DELTA/(1.-RATIO)
        IA=IA+1
        GOTO750
C
 890    TEMNEW=TEMPNW(ITER)
        IOUT=0
        IF(IOUT.EQ.0.OR.IO.LT.3)GOTO900
        WRITE(6,666)ITER,TEMNEW,SLOPE,EPSILN
        WRITE(6,203)
 203    FORMAT('  ')
 900    CONTINUE
C       WRITE(6,901)THETAE,PRESNW,TEMNEW,I1,I2,IA
C901    FORMAT(' TRAVELLED DOWN MOIST ADIABAT, THETAE=',F7.2
C     #  ,' NEW PRESSURE=',F7.2,' NEW TEMP=',F7.2,3I3)
        RETURN
        END
C
!
C      1991     Steve Albers
C
       FUNCTION THAE(TC,TD,P)
C
C   COMPUTES THE EQUIVALENT POTENTIAL TEMPURATURE (K).
C    (USING THE ROSSBY DEFINITION)
        ES(X)=6.1078+X*(.443652+X*(.014289+X*(2.65065E-4+X*
     1 (3.03124E-6+X*(2.034081E-8+X*(6.13682E-11))))))
        T=TC+273.15
        E=ES(TD)
        PMEI=1./(P-E)
        W=(E*.62197)*PMEI
        T=T*(1.+W)/(1.+W*.62197)
       CP=0.2396
       THAE=THD(P,T,W,PMEI) * EXP(RL(TC)*W/(CP*T))
       RETURN
       END
C
       FUNCTION THD(P,T,W,PMEI)
C   COMPUTES THE DRY AIR POTENTIAL TEMPERATURE
       real*4 AK
       parameter (AK=.28613105)

       AKS=AK * (1.0+1.608*W)/(1.0 + 1.941569*W)
       THD=T * ((1000.0*PMEI)**AKS)
       RETURN
       END

       FUNCTION RL(TM2)
C   LATENT HEAT OF EVAPORATION
C      TM2=T-273.15
C      RL=597.31-0.589533*TM2+0.001005333*(TM2*TM2)
       RL=597.31-((0.589533+0.001005333*TM2)*TM2)
       RETURN
       END
C
C
       SUBROUTINE DRYAD(P1,P2,T1,T2,W,IO)

!       Steve Albers 1991

C   COMPUTES THE DRY AIR POTENTIAL TEMPERATURE
       AK=.28613105
C       AKS=AK * (1.0+1.608*W)/(1.0 + 1.941569*W)
C       E=W*P1/(0.62197 +W)
C       PMEI=1./(P1-E)
C       T2=T1 * ((P2*PMEI)**AKS)
        T2=T1*(P2/P1)**AK
        IF(IO.GE.2)WRITE(6,1)P1,P2,T1,T2,W
 1      FORMAT(' P1,P2,T1,T2,W',4F10.3,F10.5)
       RETURN
       END
C
       FUNCTION THETE(T,TD,ALT,HGT)
C  T TEMP (C), TD DEW PT (C), ALT (ALTIMETER SETTING IN.)
C  HGT HEIGHT ASL (M).
C  CONVERT PA FROM INCHES TO MB
       PA=ALT*33.8639
C  SFC PRESSURE
       PA=PA*(1.0/(1.0+HGT*2.2222E-05))**5.25530
        THETE=THAE(T,TD,PA)
        RETURN
        END
C
C
C
       FUNCTION XMXRAT(PRES,DEWP)
C   COMPUTE MIXING RATIO (GM/GM) GIVEN DEW POINT TEMP
C   AND THE PRESSURE (MB)
       RATMIX=EXP(21.16-5415.0/DEWP)
       RATMIX=RATMIX/PRES
       IF(RATMIX.LT.(5.0E-05)) RATMIX=5.0E-05
       XMXRAT=RATMIX
       RETURN
       END
C
C
C
       FUNCTION THAEK(P,TC,W)
C   COMPUTES THE EQUIVALENT POTENTIAL TEMPURATURE (K).
C    (USING THE ROSSBY DEFINITION)
        Q=W/(1.0+W)
        E=(P*Q)/.62197
        PMEI=1./(P-E)
        T=TC+273.15
        T=T*(1.+W)/(1.+W*.62197)
       CP=0.2396
       THAEK=THD(P,T,W,PMEI) * EXP(RL(TC)*W/(CP*T))
       RETURN
       END
C

        function oe_fast(t,td,pres)

!       Steve Albers 1991

!       t  in deg c  (Input) Valid input temp range is -60C to +60C
!       td in deg c  (Input) Max dewpoint depression allowed is 30
!       pres in mb   (Input) Min allowed is 500mb
!       oe_fast in C (Output)

!       Quick way to get Theta(e) using lookup table

!       1991    Steve Albers

        character*31 ext
        character*150 directory
        integer*4 len_dir

        integer*4 pres2_low,pres2_high,pres2_interval,n_pres2
        parameter (pres2_low = 500)
        parameter (pres2_high = 1100)
        parameter (pres2_interval = 10)

        parameter (n_pres2 = (pres2_high - pres2_low) / pres2_interval)

        integer*4 t_low,t_high,t_interval,n_t
        parameter (t_low  = -60)
        parameter (t_high = +60)
        parameter (t_interval = 2)

        parameter (n_t = (t_high - t_low) / t_interval)

        integer*4 tdprs_low,tdprs_high,tdprs_interval,n_tdprs
        parameter (tdprs_low  = 0)
        parameter (tdprs_high = +30)
        parameter (tdprs_interval = 1)

        parameter (n_tdprs = (tdprs_high - tdprs_low) / tdprs_interval)

        real*4 thetae_lut(0:n_t,0:n_tdprs,0:n_pres2)

        save init,thetae_lut
        data init/0/

        if(init .eq. 0)then
            ext = 'dat'
            call get_directory(ext,directory,len_dir)

            write(6,*)' Reading thetae_lut.dat from file'
            open(11,file=directory(1:len_dir)//'thetae_lut.dat'
     1                ,form='unformatted',status='old',err=10)
            read(11,err=10)thetae_lut
            close(11)
            init = 1
            goto20
10          write(6,*)' Generating thetae_lut.dat - no valid file exists
     1'

            open(12,file=directory(1:len_dir)//'thetae_lut.dat'
     1                          ,form='unformatted',status='new')
            i = 0
            do t = t_low,t_high,t_interval

                j = 0
                do tdprs = tdprs_low,tdprs_high,tdprs_interval

                    k = 0
                    do pres2 = pres2_low,pres2_high,pres2_interval

                        td = t-tdprs
                        thetae_lut(i,j,k) = OE(t,td,pres2)

                        if(i/10*10 .eq. i .and. j/10*10 .eq. j .and.
     1                                        k/10*10 .eq. k)then
                            write(6,201)t,td,pres2,thetae_lut(i,j,k)
201                         format(1x,2f10.0,f10.0,f10.2)
                        endif

                        k = k + 1

                        enddo
                    j = j + 1

                enddo
                i = i + 1

            enddo

            write(12)thetae_lut
            close(12)

20      endif

        tdprs = t - td

        if(t     .ge. float(t_high))    t     = float(t_high)     - 1e-5
        if(tdprs .ge. float(tdprs_high))tdprs = float(tdprs_high) - 1e-5
        if(pres  .ge. float(pres2_high))pres  = float(pres2_high) - 1e-2

        ri = (t -     float(t_low)    ) / float(t_interval)
        rj = (tdprs - float(tdprs_low)) / float(tdprs_interval)
        rk = (pres -  float(pres2_low)) / float(pres2_interval)

        i = int(ri)
        j = int(rj)
        k = int(rk)

        fraci = ri - i
        fracj = rj - j
        frack = rk - k

        Z1=thetae_lut(i  , j  ,k)
        Z2=thetae_lut(i+1, j  ,k)
        Z3=thetae_lut(i+1, j+1,k)
        Z4=thetae_lut(i  , j+1,k)

        thetae_low =  Z1+(Z2-Z1)*fraci+(Z4-Z1)*fracj
     1                    - (Z2+Z4-Z3-Z1)*fraci*fracj

        Z1=thetae_lut(i  , j  ,k+1)
        Z2=thetae_lut(i+1, j  ,k+1)
        Z3=thetae_lut(i+1, j+1,k+1)
        Z4=thetae_lut(i  , j+1,k+1)

        thetae_high =  Z1+(Z2-Z1)*fraci+(Z4-Z1)*fracj
     1                    - (Z2+Z4-Z3-Z1)*fraci*fracj

        oe_fast = thetae_low * (1. - frack) + thetae_high * frack

c       residual = abs(oe_fast - oe(t,t-tdprs,pres))
c       write(6,101)t,t-tdprs,pres,oe(t,t-tdprs,pres),residual
c101    format(1x,3f10.2,f10.4)

        return
        end

        function tmlaps_fast(thetae,pres)

!       Steve Albers 1991

!       Quick way to get temperature along a moist adiabat given theta(e)
!       This uses a lookup table

!       thetae in K  (Input)
!       pres in mb   (Input)
!       t moist in K (Output)

!       1991    Steve Albers

        character*31 ext
        character*150 directory
        integer*4 len_dir

        integer*4 thetae_low,thetae_high,thetae_interval,n_thetae
        parameter (thetae_low  = 210)
        parameter (thetae_high = 410)
        parameter (thetae_interval = 1)

        parameter (n_thetae = (thetae_high - thetae_low) / thetae_interv
     1al)

        integer*4 pres_low,pres_high,pres_interval,n_pres
        parameter (pres_low = 50)
        parameter (pres_high = 1100)
        parameter (pres_interval = 5)

        parameter (n_pres = (pres_high - pres_low) / pres_interval)

        real*4 t_moist(0:n_thetae,0:n_pres)

        save init,t_moist
        data init/0/

        if(init .eq. 0)then
            ext = 'dat'
            call get_directory(ext,directory,len_dir)

            write(6,*)' Reading tmlaps_lut.dat from file'
            open(11,file=directory(1:len_dir)//'tmlaps_lut.dat'
     1          ,form='unformatted',status='old',err=10)
            read(11,err=10)t_moist
            close(11)
            init = 1
            goto20
10          write(6,*)' Generating tmlaps_lut.dat - no valid file exists
     1'

            open(11,file=directory(1:len_dir)//'tmlaps_lut.dat'
     1                             ,form='unformatted',status='new')

            i = 0
            do thetae = thetae_low,thetae_high,thetae_interval
                j = 0

                do pres = pres_low,pres_high,pres_interval

                    t_moist(i,j) = tmlaps(thetae-273.15,pres)

                    if(i/5*5 .eq. i .and. j/5*5 .eq. j)then
                        write(6,201)pres,thetae,t_moist(i,j)
201                     format(1x,f10.0,f10.0,f10.2)
                    endif
                    j = j + 1

                enddo
                i = i + 1

            enddo

            write(11)t_moist
            close(11)

20      endif

        ri = (thetae - float(thetae_low)) / float(thetae_interval)
        rj = (pres -   float(pres_low)  ) / float(pres_interval)

        if(ri .ge. float(n_thetae))ri = float(n_thetae) - 1e-2
        if(rj .ge. float(n_pres))  rj = float(n_pres)   - 1e-2

        i = int(ri)
        j = int(rj)

        fraci = ri - i
        fracj = rj - j

        Z1=t_moist(i  , j  )
        Z2=t_moist(i+1, j  )
        Z3=t_moist(i+1, j+1)
        Z4=t_moist(i  , j+1)

        tmlaps_fast =  Z1+(Z2-Z1)*fraci+(Z4-Z1)*fracj
     1                    - (Z2+Z4-Z3-Z1)*fraci*fracj

c       residual = abs(tmlaps_fast - tmlaps(thetae-273.15,pres))
c       write(6,101)tmlaps_fast,tmlaps(thetae-273.15,pres),residual
c101    format(1x,2f10.2,f10.4)

        return
        end

        FUNCTION DWPT_laps(T,RH)
C
C   THIS FUNCTION RETURNS THE DEW POINT (CELSIUS) GIVEN THE TEMPERATURE
C   (CELSIUS) AND RELATIVE HUMIDITY (%).
C
C       BAKER,SCHLATTER 17-MAY-1982     Original version
C
C       THE FORMULA IS USED IN THE
C   PROCESSING OF U.S. RAWINSONDE DATA AND IS REFERENCED IN PARRY, H.
C   DEAN, 1969: "THE SEMIAUTOMATIC COMPUTATION OF RAWINSONDES,"
C   TECHNICAL MEMORANDUM WBTM EDL 10, U.S. DEPARTMENT OF COMMERCE,
C   ENVIRONMENTAL SCIENCE SERVICES ADMINISTRATION, WEATHER BUREAU,
C   OFFICE OF SYSTEMS DEVELOPMENT, EQUIPMENT DEVELOPMENT LABORATORY,
C   SILVER SPRING, MD (OCTOBER), PAGE 9 AND PAGE II-4, LINE 460.
        X = 1.-0.01*RH
C   COMPUTE DEW POINT DEPRESSION.
        DPD =(14.55+0.114*T)*X+((2.5+0.007*T)*X)**3+(15.9+0.117*T)*X**14
        DWPT_laps = T-DPD
        RETURN
        END


        function twet_fast(t_c,td_c,pres_mb)

!       Steve Albers 1991
!       This is a fast routine using lookup tables
!       WARNING: This routine is only active below 500mb due to size restriction
!       of the lookup table. Max dewpoint depression allowed is 30
!       Further approximation used (twet = t_c) when t_c is outside valid range

!       twet_fast in C (Output)

        if(pres_mb .ge. 500. .and. t_c .ge. -60. 
     1                       .and. t_c .le. +60.)then ! Valid ranges to input
            thetae_k = oe_fast(t_c,td_c,pres_mb) + 273.15
            twet_fast = tmlaps_fast(thetae_k,pres_mb)

        else
            twet_fast = t_c

        endif

        return
        end

       subroutine li_laps(tc,td,pr,t500laps,i4time,imax,jmax,li,flag,ist
     1atus)

!      1991     Steve Albers

       include 'lapsparms.inc' ! for r_missing_data

       real*4 tc(imax,jmax) ! Input T  in deg F
       real*4 td(imax,jmax) ! Input Td in deg F
       real*4 pr(imax,jmax) ! Input Pr in MB
       real*4 li(imax,jmax) ! Output Li in Deg K/C

       real*4 t500laps(imax,jmax) ! Used Locally Only

       character*13 filename13

!      Get 500 mb temp field
       if(flag .eq. 0.0)then ! Use LT1 (or equivalent) file
           write(6,*)' Reading 500 temp from LT1 (or equivalent) file'
           k_level = 500
           call get_temp_2d(i4time,7200,i4time_nearest
     1                          ,k_level,imax,jmax,t500laps,istatus)

           if(istatus .ne. 1)then
               write(6,*)' Aborting Li calculation, no 500 temps'
               goto999
           endif

       elseif(abs(flag) .ge. 1e6)then ! Use Sounding file

           i4time_round = ((i4time+21600)/43200) * 43200
           open(52,file='disk$laps_s90x:'
     1          //filename13(i4time_round,'dmn'),status='old')
           read(52,*)
           read(52,*)
           read(52,*)
           read(52,*)
           read(52,*)

10         read(52,*,end=20)pres,temp

           if(pres .eq. 500.)then
               write(6,*)' Using sounding 500mb temp',temp

               DO i=1,imax
               DO j=1,jmax
                   t500laps(i,j) = temp + 273.15
               enddo
               enddo

               goto50

           else
               goto10

           endif

20         write(6,*)' Could not find 500mb temp'
           istatus = 0
           close(52)
           return

50         close(52)

       else ! Use input temp
           write(6,*)' USING INPUT 500MB TEMP'
           DO i=1,imax
           DO j=1,jmax
                t500laps(i,j) = flag
           enddo
           enddo

       endif

!      Calculate Fields
       do j=1,jmax
          do i=1,imax

              if(abs(tc(i,j)) .le. 1e5
     1   .and. abs(td(i,j)) .le. 1e5
     1   .and. abs(pr(i,j)) .le. 1e5
     1   .and. abs(t500laps(i,j)) .le. 1e5)then
                  tcij_c = (tc(i,j)-32.)/1.8
                  tdij_c = (td(i,j)-32.)/1.8
                  t500laps_c = t500laps(i,j) - 273.15
                  li(i,j) = func_li(TCij_c,TDij_c,PR(i,j),t500laps_c
     1                             ,r_missing_data)
              else
                  li(i,j) = r_missing_data
              endif

          end do ! j
      end do ! i

999   return
      end

        function func_li(t_c,td_c,psta_mb,t500_c,r_missing_data)

        real*4 li

        td_in = min(td_c,t_c)

        thetae = THAE(t_c,td_in,psta_mb)

!       t500_parcel = t500(thetae)
        call thetae_to_t500(thetae,t500_parcel,istatus)

        if(istatus .eq. 1)then
            func_li = t500_c - t500_parcel

        else
            write(6,1,err=2)psta_mb,t_c,td_c
1           format(' WARNING: setting li to missing - p,t,td= ',3f10.0)       
2           func_li = r_missing_data

        endif

        return
        end


!       function t500(thetae)
        subroutine thetae_to_t500(thetae,t500,istatus)

        real diff(10)
C
C CALCULATE LI
        ITER=1
        ITERC1=1
        ITERC2=1

C ESTIMATE PARCEL TEMP AT 500MB
        T500NW=THETAE-(307.260+72.122*EXP((THETAE-382.635)*.0141672))
        T500NW=T500NW+0.65*EXPM(-(.077*(T500NW+27.))**2)
        T500NW=T500NW-0.50*EXPM(-(.200*(T500NW+4.0))**2)
!       IF(OUTPUT.EQ.1)WRITE(31,*)THETAE,' T500 ',T500NW
        IF(THETAE.LT.254.90.OR.THETAE.GT.361.53)THEN
            T500OL=T500NW
c           WRITE(6,*)' USING ITERATIVE METHOD TO GET 500MB TEMP'
c           WRITE(31,*)' USING ITERATIVE METHOD TO GET 500MB TEMP'
C
C ENTER ITERATIVE LOOP IF THETAE IS OUT OF RANGE OF APPROXIMATE FORMULA
C DETERMINE LATENT HEAT RELEASED ABOVE 500MB
750         FUDGE=THAE(T500OL,T500OL,500.)-1.21937*(T500OL+273.15)
            THP5=THETAE-FUDGE
            T500NW=(THP5/1.21937)-273.15
!           WRITE(31,666)
!       1   ITER,iok,jok,XOK,YOK,IUS,JUS,XUS,YUS,T500NW,T500OL,FUDGE,RATIO
!666        FORMAT(I1,2I3,F6.2,F4.1,2I3,2F7.3,2X,4F10.5)
            ITERC1=ITERC1+1

            if(ITERC1 .gt. 1000)then
                write(6,*)
     1               ' WARNING in thetae_to_t500, too many iterations, '       
     1              ,' thetae = ',thetae
                istatus = 0
                return
            endif

C TEST FOR CONVERGENCE
            DIFF(ITER)=T500NW-T500OL
            IF(ABS(DIFF(ITER)).LT..10)GOTO900
            ITER=ITER+1
            ITERC2=ITERC2+1
C
C USE AITKEN'S FORMULA TO ACCELERATE CONVERGENCE
            IF(ITER.EQ.3)THEN
                RATIO=DIFF(2)/DIFF(1)
!               WRITE(31,666)
!       1       ITER,iok,jok,XOK,YOK,IUS,JUS,XUS,YUS,T500NW,DIFF(2)
!       1           ,DIFF(1),RATIO
                T500NW=T500OL+DIFF(2)/(1.-RATIO)
                T500OL=T500NW
                ITER=1
                ITERAI=ITERAI+1
                GOTO750
            endif

800         T500OL=T500NW
            GOTO750

        ENDIF

900     t500 = t500nw

        istatus = 1
        return
        end


        function expm(x)

        if(x .ge. -70.)then
            expm = exp(x)
        else
            expm = 0.
        endif

        end


        function devirt_td(t_k,td_k,p_pa)

!       This function yields the approximate temperature given the virtual temp
!       tv from the mthermo library is called. Please suggest improvements
!       to this routine to Steve Albers at FSL.

!       t_k       Input temp in K
!       td_k      Input dew point temp in K
!       p_pa      Input pressure as pascals
!       devirt_td Output devirtualized temp as K

        t_c  = t_k  - 273.15
        td_c = td_k - 273.15
        p_mb = p_pa / 100.

        tv_c = tv(t_c,td_c,p_mb)

        tv_k = tv_c + 273.15

        devirt_td = t_k - (tv_k-t_k)

        return
        end

        function devirt_sh(t_k,sh,p_pa)

!       This function yields the approximate temperature given the virtual temp
!       tv from the mthermo library is called. Please suggest improvements
!       to this routine to Steve Albers at FSL.

!       t_k       Input temp in K
!       sh        Input specific humidity (dimensionless)
!       p_pa      Input pressure as pascals
!       devirt_td Output devirtualized temp as K

        t_c  = t_k  - 273.15
        p_mb = p_pa / 100.

        tv_c = tv_sh(t_c,sh,p_mb)

        tv_k = tv_c + 273.15

        devirt_sh = t_k - (tv_k-t_k)

        return
        end

        function devirt_rh(t_k,rh,p_pa)

!       This function yields the approximate temperature given the virtual temp
!       devirt_td from the laps library is called. Please suggest improvements
!       to this routine to Steve Albers at FSL.

!       t_k       Input temp in K
!       rh        Input rh as fraction
!       p_pa      Input pressure as pascals
!       devirt_rh Output devirtualized temp as K

        t_c  = t_k  - 273.15

        rh_pct = rh * 100.

        td_c = dwpt(t_c,rh_pct)
        td_k = td_c + 273.15

        devirt_k = devirt_td(t_k,td_k,p_pa)

        devirt_rh = devirt_k

        return
        end

        FUNCTION TV_SH(T,SH,P)
C
C   THIS FUNCTION RETURNS THE VIRTUAL TEMPERATURE TV (CELSIUS) OF
C   A PARCEL OF AIR AT TEMPERATURE T (CELSIUS), DEW POINT TD
C   (CELSIUS), AND PRESSURE P (MILLIBARS). THE EQUATION APPEARS
C   IN MOST STANDARD METEOROLOGICAL TEXTS.
C
C       BAKER,SCHLATTER 17-MAY-1982     Original version
C       ALBERS                 1994     Modified for SH Input
C
        DATA CTA,EPS/273.16,0.62197/
C   CTA = DIFFERENCE BETWEEN KELVIN AND CELSIUS TEMPERATURES.
C   EPS = RATIO OF THE MEAN MOLECULAR WEIGHT OF WATER (18.016 G/MOLE)
C         TO THAT OF DRY AIR (28.966 G/MOLE)
        TK = T+CTA
C   CALCULATE THE DIMENSIONLESS MIXING RATIO.
        W = SH / (1. - SH)
        TV_SH = TK*(1.+W/EPS)/(1.+W)-CTA
        RETURN
        END

      FUNCTION TSA_fast(OS,P)
C
C   THIS FUNCTION RETURNS THE TEMPERATURE TSA (CELSIUS) ON A SATURATION
C   ADIABAT AT PRESSURE P (MILLIBARS). OS IS THE EQUIVALENT POTENTIAL
C   TEMPERATURE OF THE PARCEL (CELSIUS). SIGN(A,B) REPLACES THE
C   ALGEBRAIC SIGN OF A WITH THAT OF B.
C
C       BAKER,SCHLATTER 17-MAY-1982     Original version
C       Modification for better convergence, Keith Brewster, Feb 1994.
C
C   B IS AN EMPIRICAL CONSTANT APPROXIMATELY EQUAL TO THE LATENT HEAT
C   OF VAPORIZATION FOR WATER DIVIDED BY THE SPECIFIC HEAT AT CONSTANT
C   PRESSURE FOR DRY AIR.
C
      REAL B
      PARAMETER (B=2.6518986)
      A= OS+273.15
C
C   Above 200 mb figure all the moisture is wrung-out, so
C   the temperature is that which has potential temp of theta-e.
C   Otherwise iterate to find combo of moisture and temp corresponding
C   to thetae.
C
      IF( p.lt.200.) THEN
        TQ=A*((P/1000.)**.286)
      ELSE
C   D IS AN INITIAL VALUE USED IN THE ITERATION BELOW.
        D= 120.
C   TQ IS THE FIRST GUESS FOR TSA.
        TQ= 253.15
        x = 0.
C
C   ITERATE TO OBTAIN SUFFICIENT ACCURACY....SEE TABLE 1, P.8
C   OF STIPANUK (1973) FOR EQUATION USED IN ITERATION.
        DO 1 I= 1,25
           TQC= TQ-273.15
           D= 0.5*D

           x_last = x

           X= A*EXP(-B*W_fast(TQC,P)/TQ)-TQ*((1000./P)**.286)

c          write(6,*)' tsa_fast: t,err= ',i,tqc,x
           IF (ABS(X).LT.1E-3) GO TO 2
c
           IF (x_last * x .lt. 0.) THEN
               slope = (x-x_last) / (tq - tq_last)
               delta = - x / slope
               ad = amin1(abs(delta),d)
               tq_last = tq
               TQ = TQ + SIGN(ad,delta)
           ELSE
               tq_last = tq
               TQ= TQ+SIGN(D,X)
           END IF

 1      CONTINUE
        END IF
 2      TSA_fast = TQ-273.15
c       write(6,*)' tsa_fast: i,t,err= ',i,tsa_fast,x
        RETURN
        END
c
