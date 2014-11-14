
        SUBROUTINE MOON_BRWN(T,MX,MY,MZ)

c       T is the Julian Date (TDT)
c       Mx, My, Mz are the Geocentric Coordinates of the Moon using the
c               Mean Equinox and Ecliptic of Date (AU)

        IMPLICIT REAL*8 (A-Z)

        include '.././include/astparms.for'

c       INTEGER DECG,RAH,MINPX,LATG,LONGG

        DATA    LL0  /   270.434164 D0/,
     1          LL1  /481267.8831   D0/,
     1  LL2  /     -.001133 D0/,
     1  LL3  /      .0000019D0/

        DATA    W0  /   334.3295556D0/,
     1  W1  /  4069.0340333D0/,
     1  W2  /     -.010325 D0/,
     1  W3  /     -.0000125D0/

        DATA    N0  /   259.183275 D0/,
     1  N1  / -1934.1420083D0/,
     1  N2  /      .0020777D0/,
     1  N3  /      .0000022D0/

        DATA    LLP0 /   279.6966777D0/,
     1  LLP1 / 36000.768925 D0/,
     1  LLP2 /      .0003027D0/

        DATA    WP0 /   281.2208444D0/,
     1  WP1 /     1.719175 D0/,
     1  WP2 /      .000452 D0/

        O = Obliquity(T)

c       write(6,*)' O = ',O

        T1=(T-2415020.0D0)/36525.D0
        T2=T1**2
        T3=T1**3

        LL =  LL0 +  LL1*T1 +  LL2*T2 +  LL3*T3
        W  =   W0 +   W1*T1 +   W2*T2 +   W3*T3
        N  =   N0 +   N1*T1 +   N2*T2 +   N3*T3
        LLP= LLP0 + LLP1*T1 + LLP2*T2
        WP =  WP0 +  WP1*T1 +  WP2*T2

c       L      =  LL - W
c       LP     =  LLP - WP
c       F      =  LL - N
c       D      =  LL - LLP

c       L_DEG  =  DMOD(L+360.D4,360.D0)
c       LP_DEG =  DMOD(LP+360.D4,360.D0)
c       F_DEG =   DMOD(F+360.D4,360.D0)
c       D_DEG  =  DMOD(D+360.D4,360.D0)
c       LL_DEG = DMOD(LL,360.D0)
c       W_DEG  = DMOD( W,360.D0)
c       N_DEG  = DMOD( N+360.D4,360.D0)
c       WP_DEG = DMOD(WP,360.D0)


c       WRITE(6,1)LL_DEG,W_DEG,N_DEG,LP_DEG,WP_DEG,D_DEG,L_DEG,F_DEG
c1      FORMAT(1X,'Mean Values: LL,W,N,LP,WP,D,L,F'/8F10.4/)

!       Additive terms
        VENUS_TERM = 346.560 + 132.870 * T1 - .0091731 * T2

        NADD =  95.96 * SIND(N)
     1      + 15.58 * SIND(N -2.3 * (T1-18.5)+276.2)

        N = N + NADD/3600.

        LLADD = + 14.27  * SIND(VENUS_TERM)
     1  + 10.71  * SIND(140.0   * (T1-18.5)+170.7)
     1  +  7.261 * SIND(N)

        LL = LL + LLADD/3600.

        WP = WP
     1  +  2.076 * SIND(N) / 3600.

        L      =  LL - W
        LP     =  LLP - WP
        F      =  LL - N
        D      =  LL - LLP

c       L_DEG  =  DMOD(L+360.D4,360.D0)
c       LP_DEG =  DMOD(LP+360.D4,360.D0)
c       F_DEG =   DMOD(F+360.D4,360.D0)
c       D_DEG  =  DMOD(D+360.D4,360.D0)
c       LL_DEG = DMOD(LL,360.D0)
c       W_DEG  = DMOD( W,360.D0)
c       N_DEG  = DMOD( N+360.D4,360.D0)
c       WP_DEG = DMOD(WP,360.D0)


c       WRITE(6,2)LL_DEG,W_DEG,N_DEG,LP_DEG,WP_DEG,D_DEG,L_DEG,F_DEG
c2      FORMAT(1X,'After additive terms: LL,W,N,LP,WP,D,L,F'/8F10.4/)



        L  = L  * RPD
        W  = W  * RPD
        N  = N  * RPD
        LP = LP * RPD
        WP = WP * RPD
        D  = D  * RPD
        F  = F  * RPD

        LONGTERMS = 0.
     1  +          13.902 *  SIN(4.*D )
     1  +        2369.902 *  SIN(2.*D )

        LONGTERMS = LONGTERMS
     1  +         191.953 *  SIN(L + 2.*D )
     1  +       22639.500 *  SIN(L        )
     1  -        4586.426 *  SIN(L - 2.*D )
     1  -          38.428 *  SIN(L - 4.*D )

        LONGTERMS = LONGTERMS
     1  -          24.420 *  SIN(   LP + 2.*D)
     1  -         668.111 *  SIN(   LP)
     1  -         165.145 *  SIN(   LP -2.*D )
     1  -         125.154 *  SIN(    D)

        LONGTERMS = LONGTERMS
     1  +          14.387 *  SIN(2.*L + 2.*D )
     1  +         769.016 *  SIN(2.*L )
     1  -         211.656 *  SIN(2.*L - 2.*D )
     1  -          30.773 *  SIN(2.*L - 4.*D )
     1  -         109.667 *  SIN(L + LP)
     1  -         205.962 *  SIN(L+LP - 2.*D )

        LONGTERMS = LONGTERMS
     1  +          14.577 *  SIN(L - LP + 2.*D)
     1  +         147.693 *  SIN(L - LP)
     1  +          28.475 *  SIN(L - LP - 2.*D)

        LONGTERMS = LONGTERMS
     1  -           7.486 *  SIN(2.*LP)
     1  -           8.096 *  SIN(2.*LP - 2.*D)

        LONGTERMS = LONGTERMS
     1  -           5.471 *  SIN(2.*F + 2.*D )
     1  -         411.608 *  SIN(2.*F )
     1  -          55.173 *  SIN(2.*F - 2.*D )

        LONGTERMS = LONGTERMS
     1  +          18.023 *  SIN(LP +    D)

        LONGTERMS = LONGTERMS
     1  -           8.466 *  SIN(L +    D)
     1  +          18.609 *  SIN(L -    D)

        LONGTERMS = LONGTERMS
     1  +          36.124 *  SIN(3.*L)
     1  -          13.193 *  SIN(3.*L - 2.*D )

        LONGTERMS = LONGTERMS
     1  -           7.649 *  SIN(2.*L + LP)
     1  -           8.627 *  SIN(2.*L + LP -2.*D)
     1  +           9.703 *  SIN(2.*L - LP)
     1  -           7.412 *  SIN(L + 2.*LP - 2.*D)

        LONGTERMS = LONGTERMS
     1  -          45.099 *  SIN(L + 2.*F)
     1  -           6.382 *  SIN(L - 2.*F + 2.*D)
     1  +          39.532 *  SIN(L - 2.*F)
     1  +           9.366 *  SIN(L - 2.*F - 2.*D)

        LONGTERMS = LONGTERMS/3600.D0
        LONG = LL + LONGTERMS
        LONG = LONG -.0072 ! correct longitude
        LONG = DMOD(LONG+360.D0,360.D0)

C       CALCULATE LATITUDE

        STERMS = 0.
     1  -         112.79 *   SIN(   D)
     1  +        2373.36 *   SIN(2.*D)
     1  +          14.06 *   SIN(4.*D)
     1  +         192.72 *   SIN(L + 2.*D)
     1  +       22609.07 *   SIN(L       )
     1  -        4578.13 *   SIN(L - 2.*D)
     1  -          38.64 *   SIN(L - 4.*D)

        STERMS = STERMS
     1  +          14.78 *   SIN(2.*L + 2.*D)
     1  +         767.96 *   SIN(2.*L       )
     1  -         152.53 *   SIN(2.*L - 2.*D)
     1  -          34.07 *   SIN(2.*L - 4.*D)

        STERMS = STERMS
     1  +          50.64 *   SIN(3.*L       )
     1  -          16.40 *   SIN(3.*L - 2.*D)

        STERMS = STERMS
     1  -          25.10 *   SIN(LP + 2.*D)
     1  +          17.93 *   SIN(LP +    D)
     1  -         126.98 *   SIN(LP       )
     1  -         165.06 *   SIN(LP - 2.*D)

        STERMS = STERMS
     1  -          16.35 *   SIN(2.*LP - 2.*D)

        STERMS = STERMS
     1  -          11.75 *   SIN(L+LP + 2.*D)
     1  -         115.18 *   SIN(L+LP       )
     1  -         182.36 *   SIN(L+LP - 2.*D)
     1  -           9.66 *   SIN(L+LP - 4.*D)

        STERMS = STERMS
     1  -          23.59 *   SIN(LP-L + 2.*D)
     1  -         138.76 *   SIN(LP-L)
     1  -          31.70 *   SIN(LP-L - 2.*D)

        STERMS = STERMS
     1  -          10.56 *   SIN(2.*L +LP)
     1  -           7.59 *   SIN(2.*L +LP - 2.*D)

        STERMS = STERMS
     1  +          11.67 *   SIN(2.*L -LP)

        STERMS = STERMS
     1  -           6.12 *   SIN(L + 2.*LP -2.*D)

        STERMS = STERMS
     1  -          52.14 *   SIN(2.*F -2.*D)

        STERMS = STERMS
     1  -           9.52 *   SIN(L + 2.*F -2.*D)

        NTERMS = 0.
     1  -         526.069 *  SIN(F - 2.*D)
     1  -           3.352 *  SIN(F - 4.*D)
     1  +          44.297 *  SIN(F+L-2.*D)
     1  -          30.598 *  SIN(F-L-2.*D)
     1  -          24.649 *  SIN(F-2.*L)
     1  -          22.571 *  SIN(F+LP-2.*D)
     1  +          20.599 *  SIN(F-L)
     1  -           6.000 *  SIN(F+L-4.*D)
     1  -           2.000 *  SIN(F-2.*L-2.*D)
     1  +          10.985 *  SIN(F-LP-2.*D)

        GAMMA1_C = -1.540 * COS(-L +LP -2.*D)
        OMEGA1 = .0004664 * COS(N)
        OMEGA2 = .0000754 * COS(N+ 275.05 - 2.30*T1)

        S = F + STERMS/3600.d0*RPD

        LAT = (18518.511+1.189+GAMMA1_C) * SIN(S)
     1          - 6.241 * SIN(3.*S) + NTERMS
        LAT = LAT * (1.-OMEGA1-OMEGA2)/3600.
        LAT = LAT + .0026 ! correct latitude (.0025)

C       CALCULATE DISTANCE

        SINEPX = 0.
     1  +             .2607 * COS(4.*D)
     1  +           28.2333 * COS(2.*D)
     1  +         3422.7000

        SINEPX = SINEPX
     1  +            3.0861 * COS(L + 2.*D)
     1  +          186.5398 * COS(L       )
     1  +           34.3117 * COS(L - 2.*D)
     1  +             .6008 * COS(L - 4.*D)

        SINEPX = SINEPX
     1  -             .3000 * COS(LP+ 2.*D)
     1  -             .3997 * COS(LP      )
     1  +            1.9178 * COS(LP- 2.*D)

        SINEPX = SINEPX
     1  -             .9781 * COS(D)

        SINEPX = SINEPX
     1  +             .2833 * COS(2.*L + 2.*D)
     1  +           10.1657 * COS(2.*L       )
     1  -             .3039 * COS(2.*L - 2.*D)
     1  +             .3722 * COS(2.*L - 4.*D)

        SINEPX = SINEPX
     1  -             .9490 * COS(L+LP)
     1  +            1.4437 * COS(L+LP-2.*D)
     1  +             .2302 * COS(L-LP+2.*D)
     1  -             .2257 * COS(L-LP-2.*D)

        SINEPX = SINEPX
     1  -             .1052 * COS(2.*F - 2.*D)

        SINEPX = SINEPX
     1  -             .1093 * COS(L+D)

        SINEPX = SINEPX
     1  +             .1494 * COS(LP + D)

        SINEPX = SINEPX
     1  +            1.1528 * COS(L-LP)

        SINEPX = SINEPX
     1  +             .6215 * COS(3.*L)
     1  -             .1187 * COS(3.*L - 2.*D)

        SINEPX = SINEPX
     1  -             .1038 * COS(2.*L + LP)
     1  +             .1268 * COS(2.*L - LP)
     1  -             .0833 * COS(L + 2.*F - 2.*D)

        SINEPX = SINEPX
     1  -             .7136 * COS(L-2.*F)

c       MINPX = IDINT(SINEPX/60.D0)
c       SECPX = SINEPX - MINPX*60.D0
c       WRITE(6,20)LONG,LAT,MINPX,SECPX
c20     FORMAT(1X,'LONG,LAT,SINEPX ',2F12.6,I4,F8.3)

c       LATD = DABS(LAT)
c       LATG = IDINT(LATD)
c       LATM = (LATD - LATG) * 60.D0
c       LATG = LATG * LAT/DABS(LAT)
c       LATS = (LATM - int(LATM)) * 60.D0

c       LONGD = DABS(LONG)
c       LONGG = IDINT(LONGD)
c       LONGM = (LONGD - LONGG) * 60.D0
c       LONGG = LONGG * LONG/DABS(LONG)
c       LONGS = (LONGM - int(LONGM)) * 60.D0

c       WRITE(6,29)LATG,int(LATM),LATS,LONGG,int(LONGM),LONGS
c29     FORMAT(/1X,'Coords of Date:   LAT',
c       1       I4,i3,F6.2,'   LONG',I4,i3,F6.2/)

        LONG = LONG * RPD
        LAT  = LAT  * RPD
        SINEPX=(SINEPX/3600.)*RPD
        R     =   1./SINEPX * 4.2635D-5

        DEC  = ASIN (COS(O)*SIN(LAT)+SIN(O)*COS(LAT)*SIN(LONG))
        RA   =-ASIN((SIN(O)*SIN(LAT)-COS(O)*COS(LAT)*SIN(LONG))/COS(DEC)
     1)
        IF(COS(LONG).LT.0.)RA = PI-RA
        RA = DMOD(RA+2.d0*PI,2.d0*PI)

c       DECD = DABS(DEC) / RPD
c       DECG = IDINT(DECD)
c       DECM = (DECD - DECG) * 60.D0
c       DECG = DECG * DEC/DABS(DEC)
c       DECS = (DECM - int(DECM)) * 60.D0

c       RAD  = (RA / RPD)/15.D0
c       RAH  = IDINT(RAD)
c       RAM  = (RAD - RAH)*60.D0
c       RAS  = (RAM - int(RAM)) * 60.D0

c       WRITE(6,30)DECG,int(DECM),DECS,RAH,int(RAM),RAS,R
c30     FORMAT(/1X,'Coords of Date:   DEC',
c       1       I4,i3,F6.2,'   RA',I4,i3,F6.2,'  R',f15.6/)

        MX = R * COS(DEC) * COS(RA)
        MY = R * COS(DEC) * SIN(RA)
        MZ = R * SIN(DEC)

        RETURN
        END

