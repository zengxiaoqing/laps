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

        subroutine put_stability(
     1           i4time_needed                   ! Input
     1          ,NX_L,NY_L,NZ_L                  ! Input
     1          ,heights_3d                      ! Input
     1          ,topo                            ! Input
     1          ,laps_cycle_time                 ! Input
     1          ,temp_3d                         ! Input
     1          ,rh_3d_pct                       ! Input
     1          ,temp_sfc_k                      ! Input
     1          ,pres_sfc_pa                     ! Input
     1          ,twet_snow                       ! Input
     1          ,td_3d_k                         ! Output
     1          ,istatus)                        ! Output

cdoc    Calculate and write out set of 2-D stability grids

!       Arrays passed in
        real temp_3d(NX_L,NY_L,NZ_L)
        real rh_3d_pct(NX_L,NY_L,NZ_L)
        real heights_3d(NX_L,NY_L,NZ_L)
        real temp_sfc_k(NX_L,NY_L)
        real pres_sfc_pa(NX_L,NY_L)
        real topo(NX_L,NY_L)

!       Output
        real td_3d_k(NX_L,NY_L,NZ_L)

!       Local declarations for stability 
        real t_sfc_f(NX_L,NY_L)
        real td_sfc_k(NX_L,NY_L)
        real td_sfc_f(NX_L,NY_L)

        real pbe_2d(NX_L,NY_L)
        real nbe_2d(NX_L,NY_L)
        real si_2d(NX_L,NY_L)
        real tt_2d(NX_L,NY_L)
        real k_2d(NX_L,NY_L)
        real lcl_2d(NX_L,NY_L)
        real wb0_2d(NX_L,NY_L)
        real wb1_2d(NX_L,NY_L)

        real li(NX_L,NY_L)

        real pres_sfc_mb(NX_L,NY_L)
        real pres_3d(NX_L,NY_L,NZ_L)

        integer nfields
        parameter (nfields=8)

        character*10 units_2d_a(nfields)
        character*125 comment_2d_a(nfields)
        character*3 var_2d_a(nfields)
        real out_multi_2d(NX_L,NY_L,nfields)

        character*31 EXT

        character*3 var_2d
        character*10  units_2d
        character*125 comment_2d

        real k_to_f

!       Read in surface dewpoint data
        var_2d = 'TD'
        ext = 'lsx'
        call get_laps_2dgrid(i4time_needed,laps_cycle_time/2
     1                      ,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L       
     1                      ,td_sfc_k,0,istatus)
        if(istatus .ne. 1)then
            write(6,*)' LAPS Sfc Dewpoint not available'
            write(6,*)' Abort put_stability routine'
            return
        endif

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)return

!       call get_pres_3d(i4time_needed,NX_L,NY_L,NZ_L,pres_3d,istatus)       
!       if(istatus .ne. 1)return

!       Convert RH to Td
        do k = 1,NZ_L
        do j = 1,NY_L
        do i = 1,NX_L
            if(rh_3d_pct(i,j,k) .ge. 0. .and.
     1         rh_3d_pct(i,j,k) .le. 100.       )then    ! RH in valid range

                t_c = temp_3d(i,j,k)-273.15
                td_c = dwpt(t_c,rh_3d_pct(i,j,k))
                td_3d_k(i,j,k) = td_c + 273.15

!               ew = pres_3d(i,j,k)/100. * sh_3d(i,j,k)
!               td_c = dpt(ew)                        ! Check valid input range

!               td_c = make_td(pres_3d(i,j,k)/100.
!    1                        ,temp_3d(i,j,k)-273.15
!    1                        ,sh_3d(i,j,k)*1000.
!    1                        ,-100.)

            elseif(rh_3d_pct(i,j,k) .gt. 100. .and.
     1             rh_3d_pct(i,j,k) .le. 101.       )then ! RH slightly high

                td_3d_k(i,j,k) = temp_3d(i,j,k)
                write(6,*)' WARNING: RH out of bounds',rh_3d_pct(i,j,k)       
     1                   ,' at ',i,j,k,' setting td = t'

            else ! invalid value of rh to pass into dwpt
                td_3d_k(i,j,k) = r_missing_data
                write(6,*)' ERROR: RH out of bounds',rh_3d_pct(i,j,k)       
     1                   ,' at ',i,j,k
                istatus = 0
                return

            endif

        enddo ! i
        enddo ! j
        enddo ! k            

        call laps_be(NX_L,NY_L,NZ_L,twet_snow
     1              ,temp_sfc_k,td_sfc_k,pres_sfc_pa
     1              ,temp_3d,td_3d_k,heights_3d,topo,blayr_thk_pa
     1              ,pbe_2d,nbe_2d,si_2d,tt_2d,k_2d,lcl_2d,wb0_2d
     1              ,wb1_2d,r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' WARNING, Bad istatus returned from laps_be'
            return
        endif

!       Fill pres_sfc_mb
        call move(pres_sfc_pa,pres_sfc_mb,NX_L,NY_L)
        call multcon(pres_sfc_mb,.01,NX_L,NY_L)

!       Convert T, Td to F
        do i = 1,NX_L
        do j = 1,NY_L
            t_sfc_f(i,j)  = k_to_f(temp_sfc_k(i,j))
            td_sfc_f(i,j) = k_to_f(td_sfc_k(i,j))
        enddo ! j
        enddo ! i

        flag = 0.0
        call li_laps(t_sfc_f,td_sfc_f,pres_sfc_mb
     1              ,i4time_needed,NX_L,NY_L,li,flag,istatus)

!       call move
        call move(pbe_2d,out_multi_2d(1,1,1),NX_L,NY_L)
        call move(nbe_2d,out_multi_2d(1,1,2),NX_L,NY_L)
        call move(    li,out_multi_2d(1,1,3),NX_L,NY_L)
        call move( si_2d,out_multi_2d(1,1,4),NX_L,NY_L)
        call move( tt_2d,out_multi_2d(1,1,5),NX_L,NY_L)
        call move(  k_2d,out_multi_2d(1,1,6),NX_L,NY_L)
        call move(lcl_2d,out_multi_2d(1,1,7),NX_L,NY_L)
        call move(wb0_2d,out_multi_2d(1,1,8),NX_L,NY_L)

!       add var arrays
        ext = 'lst'

        var_2d_a(1) = 'PBE'
        var_2d_a(2) = 'NBE'
        var_2d_a(3) = 'LI'
        var_2d_a(4) = 'SI'
        var_2d_a(5) = 'TT'
        var_2d_a(6) = 'K'
        var_2d_a(7) = 'LCL'
        var_2d_a(8) = 'WB0'

        units_2d_a(1) = 'J/KG'
        units_2d_a(2) = 'J/KG'
        units_2d_a(3) = 'K'
        units_2d_a(4) = 'K'
        units_2d_a(5) = 'K'
        units_2d_a(6) = 'K'
        units_2d_a(7) = 'M'
        units_2d_a(8) = 'M'

        comment_2d_a(1) = 'CAPE'
        comment_2d_a(2) = 'CIN'
        comment_2d_a(3) = 'Lifted_Index'
        comment_2d_a(4) = 'Showalter_Index'
        comment_2d_a(5) = 'Total_Totals'
        comment_2d_a(6) = 'K_Index'
        comment_2d_a(7) = 'LCL'
        comment_2d_a(8) = 'Wet_Bulb_Zero'

        call put_laps_multi_2d(i4time_needed,ext,var_2d_a,units_2d_a
     1                        ,comment_2d_a,out_multi_2d,NX_L,NY_L,8    
     1                        ,istatus)
        if(istatus .ne. 1)then
            write(6,*)' LST output error'
        else
            write(6,*)' Successfully wrote LST'
        endif

        return

        end


        subroutine laps_be(ni,nj,nk,twet_snow
     1        ,t_sfc_k,td_sfc_k,p_sfc_pa,t_3d_k,td_3d_k,ht_3d_m,topo       
     1        ,blayr_thk_pa,pbe_2d,nbe_2d,si_2d,tt_2d,k_2d,lcl_2d,wb0_2d
     1        ,wb1_2d,r_missing_data,istatus)

!       1991    Steve Albers
cdoc    Returns 2-D PBE and NBE in Joules, Parcel is lifted from lowest level
!                                                                    i.e. sfc

        real t_sfc_k(ni,nj)
        real td_sfc_k(ni,nj)
        real p_sfc_pa(ni,nj)
        real topo(ni,nj)
        real t_3d_k(ni,nj,nk)
        real td_3d_k(ni,nj,nk)
        real ht_3d_m(ni,nj,nk)
        real p_1d_pa(nk)
        real p_1d_mb(nk)                   ! Local
        real pres_3d(ni,nj,nk)             ! Local
        real pbe_2d(ni,nj)
        real nbe_2d(ni,nj)

        real si_2d(ni,nj)
        real tt_2d(ni,nj)
        real k_2d(ni,nj)
        real lcl_2d(ni,nj)
        real wb0_2d(ni,nj)
        real wb1_2d(ni,nj)
        
        include 'lapsparms.for'
        integer MXL
        parameter (MXL=MAX_LVLS+1) ! number of 3D levels plus the sfc level

        COMMON/INDX/ P(MXL),T(MXL),TD(MXL),HT(MXL),PBECR(20,4)
     1  ,TDFCR(20,2),VEL(20)
     1  ,temdif(MXL),partem(MXL),pbe(MXL)
     #  ,DD85,FF85,DD50,FF50
        REAL LCL,LI,K_INDEX

!       Initialize pbe array
        do i = 1,MXL
            pbe(i) = 0.
        enddo ! i

        call get_systime_i4(i4time,istatus)
        if(istatus .ne. 1)stop

        call get_pres_3d(i4time,ni,nj,nk,pres_3d,istatus)
        if(istatus .ne. 1)stop

        if(nj .gt. 600)then
            jint = 40
        else
            jint = 20
        endif

        do j = 1,nj
        do i = 1,ni
c       write(6,*)' i = ',i

            IF(j .eq. (j/jint)*jint .AND. i .eq. ni/2)then
                idebug = 2 ! 1
                write(6,*)
                write(6,*)
     1          ' --- stability index debugging at gridpoint ---',i,j
            else 
                idebug = 0
            endif
        
            do k = 1,nk
                p_1d_pa(k) = pres_3d(i,j,k)
                p_1d_mb(k) = p_1d_pa(k) / 100.
            enddo ! k

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
                T(N) = t_3d_k(i,j,k)  - 273.15 ! K to C
                TD(N)= td_3d_k(i,j,k) - 273.15 ! K to C
!               TD(N)= t(n)
                HT(N)= ht_3d_m(i,j,k)
            enddo ! k

            nlevel = n

            IO = 0

            CALL SINDX(NLEVEL,LI,SI,BLI,TT,SWEAT
     1                ,twet_snow                                        ! I
     1                ,HWB0,HWB_snow                                    ! O
     1                ,PLCL,LCL,CCL
     1                ,TCONV,IO
     1                ,ICP,ICT,K_INDEX,TMAX,PBENEG,PBEPOS,T500,PBLI
     1                ,VELNEG,WATER,IHOUR,idebug,istatus)
            if(istatus .ne. 1)then
                write(6,*)' WARNING: Bad istatus returned from SINDX'
                return
            endif

            pbe_2d(i,j) = PBEPOS
            nbe_2d(i,j) = PBENEG
          
            si_2d(i,j) = SI
            tt_2d(i,j) = TT
            k_2d(i,j)  = K_INDEX

            if(LCL .ne. r_missing_data)then
                lcl_2d(i,j) = LCL                        ! M MSL
            else
                lcl_2d(i,j) = r_missing_data
            endif

            if(HWB0 .ne. r_missing_data)then
                wb0_2d(i,j) = HWB0*304.8006 + topo(i,j)      ! KFT AGL to M MSL
            else
                wb0_2d(i,j) = r_missing_data
            endif

            if(HWB_snow .ne. r_missing_data)then
                wb1_2d(i,j) = HWB_snow*304.8006 + topo(i,j)  ! KFT AGL to M MSL
            else
                wb1_2d(i,j) = r_missing_data
            endif

            iwarn = 0
            if(nanf(wb0_2d(i,j)) .eq. 1)then
                write(6,*)' Warning: HWB0 Nan ',i,j
                iwarn = 1
            endif

            IF(idebug .ge. 1 .OR. iwarn .eq. 1)then

                write(6,*)
                write(6,*)' Indices at:',i,j

                write(6,*)' n    p         t         td'
     1                   ,'         ht       tdif       be'
                do n = 1,nlevel
                    write(6,301)n,p(n),t(n),td(n),ht(n),temdif(n),pbe(n)
301                 format(1x,i2,f8.1,f10.2,f10.2,f10.0,f10.2,f10.1)
                enddo

                SWEAT=0.0
!               IF(IFLAG1*IFLAG2.EQ.0)SWEAT=0.0

                WRITE(6,420)PBENEG,PBEPOS
 420            FORMAT(' PBENEG',F8.1,'               PBEPOS',F8.1)

                WRITE(6,430)DD85,FF85,DD50,FF50,T500
 430            FORMAT(' 850MB WIND',2F5.0,'         500MB WIND',2F5.0
     #         ,'    500MB TEMP= ',F6.1)

!               IF(IO2.EQ.1)GOTO1000

                WRITE(6,60)LI,SI,BLI
 60             FORMAT(' LI= ',F5.1,20X,'SI= ',F5.1,15X,'BLI= ',F5.1)

                WRITE(6,61,err=161)TT,SWEAT,HWB0
 61             FORMAT(' TOTAL TOTALS=',F6.1,10X,'SWEAT= ',F7.1
     1                ,10X,'W. BULB ZERO=',F5.1,' KFT AGL')

 161            ITCONV=INT(TCONV+0.5)

                WRITE(6,62)LCL*.003281,CCL,ITCONV
 62             FORMAT(' LCL= ',F5.1,' KFT MSL',11X,'CCL=',F5.1
     1                ,' KFT AGL',7X,'CONVECTIVE TEMP=',I4,' DEG F')

                ITMAX=NINT(TMAX)

                WRITE(6,63,err=163)K_INDEX,ITMAX,WATER
 63             FORMAT(' K INDEX =',f7.1,13X,'TMAX =',I4,' F',12X
     1                ,'PRECIP. WATER=',F5.2,' IN.')
 163            WRITE(6,*)
C
                IF(ICT+ICP.GT.0.OR.IO.GE.2)WRITE(6,71)P(1),PLCL
 71             FORMAT(' ENERGY ANALYSIS - LIFTING A PARCEL FROM'
     1                ,F6.0,'MB','   LCL=',f7.1,' MB')
C
                DO 600 II=1,ICP
 600            WRITE(6,64)PBECR(II,1),PBECR(II,2),VEL(II),PBECR(II,4)
 64             FORMAT(' ENERGY AT',F6.0,'MB=',F10.2
     1                ,' JOULES/KG       VELOCITY=',F6.2,'M/S   HT='
     1                ,F7.0,' M')
C
                WRITE(6,*)
                DO 700 II=1,ICT
 700            WRITE(6,5)TDFCR(II,1),TDFCR(II,2)
 5              FORMAT(5X,'ENVIRONMENTAL MINUS PARCEL TEMPERATURE AT'
     1                ,F6.0,'MB  =',F6.1)
C
                IF(BLI.LE.0.0)WRITE(6,65)PBENEG,VELNEG
                IF(BLI.LE.0.0)WRITE(6,66)PBEPOS
 65             FORMAT('     CAP STRENGTH =',F10.1,' JOULES/KG.'
     +                ,' VELOCITY NEEDED',F5.1,'M/S')
 66             FORMAT('     POSITIVE AREA=',F10.1,' JOULES/KG.')
                GOTO2000
C
C  ALTERNATIVE OUTPUTTING
 1000           CONTINUE
                WRITE(6,1001)STA,T500,VELNEG,PBEPOS
 1001           FORMAT(4X,A3/F5.1,F7.1,F7.1)

2000        endif

        enddo ! i
        enddo ! j
C
 9999   CONTINUE

        RETURN
        END

!
        SUBROUTINE SINDX(NLEVEL,LI,SI,BLI,TT,SWEAT
     1   ,twet_snow,HWB0,HWB_snow
     1   ,PLCL_PBE,LCL_PBE_MSL,CCL  
     1   ,TCONV,IO,ICP,ICT,K,TMAX,PBENEG,PBEPOS,TMAN50
     1   ,PBLI,VELNEG,WATER,IHOUR,idebug,istatus)

cdoc    Calculate a variety of stability indices from an input sounding

!       1991    Steve Albers
!       1999    Steve Albers    Adding more indices to active output

!       Note that a surface parcel is currently used for the indices

        include 'lapsparms.for'
        integer MXL
        parameter (MXL=MAX_LVLS+1) ! number of 3D levels plus the sfc level

        DIMENSION Q(MXL),W(MXL),WB(MXL)
        COMMON/INDX/ P(MXL),T(MXL),TD(MXL),HT(MXL),PBECR(20,4)
     1              ,TDFCR(20,2),VEL(20),temdif(MXL),partem(MXL)
     1              ,pbe(MXL),DD85,FF85,DD50,FF50
        REAL LI,K,LCL_AGL,LCL_PBE_MSL
        ES(X)=6.1078+X*(.443652+X*(.014289+X*(2.65065E-4+X*
     1 (3.03124E-6+X*(2.034081E-8+X*(6.13682E-11))))))
!       TDEW(E)=237.7/((7.5/ALOG10(E/6.11))-1.)

        logical l_large_domain

        DATA EPSILN/.6220/,G/9.80665/
        DATA BLTHCK/50.0/
        DATA RPD/.0174532925063/
 1      format('    MIXED PARCEL IS:',F10.1,'MB',F15.3,'C       '
     1          ,'MIXING RATIO=',F10.7)

c       WRITE(6,15)
 15     FORMAT('          PRESSURE     TEMP.    DEW PT.     '
     1          ,'Q      WET BULB')

        IOUT=IO

        l_large_domain = .false.

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            return
        endif

!       Fill WB array (imported from old code)
 	DO 100 N=1,NLEVEL                                            
            IF(TD(N).EQ.-99.0)THEN                                         
                Q(N)=-99.                                                      
                WB(N)=-99.                                                     
  	        IF(IO.GE.2)WRITE(6,3)N,P(N),T(N),TD(N),Q(N),WB(N)    
                IFL=1                                                
                TD(N)=T(N)-15.0                                      
                IF(T(N).LT.-40.)TD(N)=T(N)                           
            ENDIF                                                           
                                                                          
 	    Q(N)=(ES(TD(N))*EPSILN)/P(N)                             
 	    Q(N)=max(Q(N),0.0001)                                    
 	    W(N)=Q(N)/(1.-Q(N))                                      
 	    IOUT=IO                                                  
 	    IF(N.LE.1)IOUT=IO                                        

 	    WB(N)=WTBLB(P(N),T(N),W(N),0,istatus)                            
            if(istatus .ne. 1)then
                write(6,*)' Bad istatus return from WTBLB: P,T,W = '
     1                   ,P(N),T(N),W(N)
                return
            endif

 	    IF(WB(N).LT.TD(N))WB(N)=TD(N)                            
 	    IF(WB(N).GT.T(N)) WB(N)=T(N)                             
 	    IF(IO.GE.2.AND.IFL.EQ.0)
     1         WRITE(6,3)N,P(N),T(N),TD(N),Q(N),WB(N)     
            IFL=0                                                            
 3	    FORMAT(' LVL(',I2,')',3F10.1,F11.6,F10.2)                      
 100	CONTINUE                                                             

!       Here is the new section imported (from old code) for TT,SI,K
C                                                                         
C  CALCULATE LIFTED INDEX                                                 
 	CALL ITPLV(P,T,NLEVEL,500.,TMAN50,IO,istatus)                         
!WNI Modified to check istatus of ITPLV.  When running LAPS
!WNI for East Asia, some of the mountains extended above the 500mb
!WNI level, which cause a complete abort of put stability.  Now,
!WNI we just fill indices dependent on 500mb to missing if this happens.
!WNI Brent Shaw, WNI, Dec 2006
!	WREQ50=ES(TP500)/500.                                                   
!	IF(WREQ50.LT.WMEAN)GOTO40                                               

C                                                                         
C  CALCULATE SHOWALTER INDEX                                              
        SI=r_missing_data
        TT=r_missing_data
        SWEAT=r_missing_data
        K=r_missing_data
        IF (istatus .EQ. 1) THEN  !WNI .. 500mb below ground
          IF(P(1).GE.850.0)THEN
            CALL ITPLV(P,T ,NLEVEL,850.,TMAN85,IO,istatus)
 	    CALL ITPLV(P,TD,NLEVEL,850.,TDMN85,IO,istatus)
 	    THETAE=THAE(TMAN85,TDMN85,850.)
 	    CALL MSAD5(TP500,500.,THETAE,25.,20.,SLOPE,I1,I2,IA,0
     1                ,istatus)    
            if(istatus .ne. 1)then
                write(6,*)' Error, skipping Showlater in SINDX'
            else
!	        IF(WREQ50.LT.WMEAN)GOTO50                                 
!	        RH=WMEAN/WREQ50                                           
 50	        SI=TMAN50-TP500                                           
            endif ! valid results from MSAD5

C                                                                         
C  CALCULATE TOTAL TOTALS AND SWEAT AND K INDICIES                        
 	    CALL ITPLV(P,T,NLEVEL,700.,TMAN70,IO,istatus)        
 	    CALL ITPLV(P,TD,NLEVEL,700.,TDMN70,IO,istatus)       
 	    TT=TMAN85+TDMN85-2.*TMAN50                   
 	    A=max(TDMN85,0.)                             
 	    B=max(TT-49.,0.)                             
 	    C=0.                                         
 	    IF(FF85.LT.15.OR.FF50.LT.15.)   GOTO500      
 	    IF(DD85.GT.DD50)                GOTO500      
 	    IF(DD85.LT.130..OR.DD85.GT.250.)GOTO500      
 	    IF(DD50.LT.210..OR.DD50.GT.310.)GOTO500      
 	    C=SIN((DD50-DD85)*RPD)+.2                    
 500	    SWEAT=12.*A+20.*B+2.*FF85+FF50+125.*C        
 	    K=TMAN85-TMAN50+TDMN85-TMAN70+TDMN70         

          ENDIF
        ELSE  !WNI 
          PRINT *,"No 500mb level, some indices set to missing" !WNI
          istatus = 1 ! WNI
        ENDIF ! WNI

C                                                                         
C  CALCULATE WET BULB ZERO LEVEL                                          
        twet0 = 0.       ! Deg C
        call wtblb_lvl(twet0,P,T,Q,WB,MXL,NLEVEL,PWB0,HWB0)

        call wtblb_lvl(twet_snow,P,T,Q,WB,MXL,NLEVEL,PWB_snow,HWB_snow)       

!       Calculate theta(e) based on sfc parcel
        if(l_large_domain)then
            THETAE=OE_FAST(T(1),TD(1),P(1)) + 273.15
        else
            THETAE=OE(T(1),TD(1),P(1)) + 273.15
        endif

        if(idebug .ge. 2)then
            if(l_large_domain)then
                THETAE2=OE(T(1),TD(1),P(1)) + 273.15
                write(6,510)P(1),T(1),TD(1),THETAE,THETAE2
 510            format(' P/T/TD/THETAE',5f8.2)
                if(abs(THETAE-THETAE2) .gt. 1.0)then
                    write(6,*)' ERROR: large difference in THETA values'
                    stop
                endif
            else
                write(6,510)P(1),T(1),TD(1),THETAE
            endif
        endif

!       Calculate "fast" LCL based on sfc parcel
        CALL LCL_fast(P(1),T(1),TD(1),LCL_AGL,TLCL_PBE,PLCL_PBE)

!       This LCL is for the surface parcel passed in
        LCL_PBE_MSL = LCL_AGL + HT(1)

!       Calculate CAPE/CIN based on sfc parcel
        CALL POTBE(Q,NLEVEL,P(1),T(1),W(1),PLCL_PBE
     1   ,TLCL_PBE,LCL_PBE_MSL,THETAE,ICP,ICT,IO,PBENEG,PBEPOS
     1   ,VELNEG,idebug)
C
        DO 600 I=1,ICP
            VEL(I)=SQRT(ABS(2.*PBECR(I,2)))
            IF(PBECR(I,2).LT.0.)VEL(I)=-VEL(I)
c           WRITE(6,4)PBECR(I,1),PBECR(I,2),VEL(I)
 4          format(' ENERGY AT',F6.0,'MB='
     1            ,F15.2,'JOULES       VELOCITY=',F6.1,'MS')
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
c        COMMON/INDX/ P(MXL),T(MXL),TD(MXL),HT(MXL),PBECR(20,4),TDFCR(20,2)
c     1              ,VEL(20),temdif(MXL),partem(MXL),pbe(MXL)
c     #              ,DD85,FF85,DD50,FF50
c        DATA PBE/MXL*0./
c        end
C
C
C
C
C
C
C
        SUBROUTINE POTBE(Q,NLEVEL,PMEAN,TMEAN,WMEAN,PLCL
     #   ,TLCL,LCL,BLTHTE,ICP,ICT,IO,PBENEG,pos_area_max
     #   ,VELNEG,idebug)

cdoc    Calculate a PBE/LCL related indices from an input sounding

!       Steve Albers 1991

        include 'lapsparms.for'
        integer MXL
        parameter (MXL=MAX_LVLS+1) ! number of 3D levels plus the sfc level

        COMMON/INDX/ P(MXL),T(MXL),TD(MXL),HT(MXL),PBECR(20,4)
     1              ,TDFCR(20,2),VEL(20)
     1              ,temdif(MXL),partem(MXL),pbe(MXL)
     #              ,DD85,FF85,DD50,FF50
        REAL LCL,nbe_min
        DIMENSION Q(MXL),DELTAH(MXL)
        real GAMMA
        parameter (GAMMA = .009760) ! Dry Adiabatic Lapse Rate Deg/m
        DATA G/9.80665/

        if(idebug .ge. 2)then
            WRITE(6,*)
            WRITE(6,*)' SUBROUTINE POTBE'
        endif

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading r_missing_data: STOP'
            stop
        endif

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
                    if(idebug .ge. 2)WRITE(6,*)' DRY CASE'
                    PARTEM(N)=PARTEM(N-1)-GAMMA*DELTAH(N)
                    TEMDIF(N)=PARTEM(N)-T(N)
                    PBE(N)=PBE(N-1)+G*(TEMDIF(N)+TEMDIF(N-1))/(T(N)
     #                          +T(N-1)+546.30)*DELTAH(N)

                else ! Upper level is above LCL
C                   BRACKETING CASE AROUND LCL - DRY ADIABATIC PART
                    if(idebug .ge. 2)WRITE(6,*)' DRY ADIABATIC PART'

                    delta_ht_dry = lcl - HT(n-1)

                    if(idebug .ge. 2)WRITE(6,307)TLCL
 307                format(' PARCEL TEMP AT LCL= ',F10.3)

                    CALL ITPLV(P,T,NLEVEL,PLCL,SNTLCL,IO,istatus)

                    t_dif_lcl=TLCL-SNTLCL

                    pbe_dry=G*(t_dif_lcl+TEMDIF(N-1))/
     1                        (SNTLCL+T(N-1)+546.30) * delta_ht_dry

 951                format(' pos_area_max,PSI,NEG,NNG,PBECR',5F9.1,I3)

!d                  WRITE(6,777)N,P(N),TLCL,SNTLCL,t_dif_lcl
!d    #                    ,TEMDIF(N-1),delta_ht_dry,HT(N),PBE(N-1)+pbe_dry

C                   MOIST ADIABATIC PART
                    if(idebug .ge. 2)WRITE(6,*)' MOIST ADIABATIC PART'
                    delta_ht_wet=DELTAH(N)-delta_ht_dry

                    partem(n) = tmlaps_fast(blthte,P(n))
                    TEMDIF(N)=PARTEM(N)-T(N)

                    pbe_wet = G*(TEMDIF(N)+t_dif_lcl)/
     1                          (T(N)+SNTLCL+546.30) * delta_ht_wet

                    PBE(N)=PBE(N-1) + pbe_dry + pbe_wet

                endif ! Upper level below LCL (Dry or bracket)

            else ! Lower Level is above LCL
                if(idebug .ge. 2)then
                  WRITE(6,*)' GETTING PARCEL TEMPERATURE FOR MOIST CASE'
     1                     ,blthte,P(n)
                endif
                partem(n) = tmlaps_fast(blthte,P(n))

          !     Add Layer
                TEMDIF(N)=PARTEM(N)-T(N)
                PBE(N)=PBE(N-1)+G*(TEMDIF(N)+TEMDIF(N-1))/(T(N)
     #                          +T(N-1)+546.30)*DELTAH(N)

            endif


            if(idebug .ge. 2)then
                WRITE(6,777)N,P(N),PARTEM(N),T(N)
     #                ,TEMDIF(N),TEMDIF(N-1),DELTAH(N),HT(N),PBE(N)
            endif
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
            if(idebug .ge. 2)WRITE(6,940)N
 940        format(' LOOKING FOR NEUTRAL BUOYANCY - ENERGY EXTREMUM, LEV
     1EL',I3)

            IF((TEMDIF(N)*TEMDIF(N-1)).lt.0.)then
                ICT=ICT+1
                SLOPE=(TEMDIF(N)-TEMDIF(N-1))/ALOG((P(N)/P(N-1)))
                FRAC=TEMDIF(N-1)/(TEMDIF(N-1)-TEMDIF(N))
                TDFCR(ICT,1)=EXP(ALOG(P(N-1))-TEMDIF(N-1)/SLOPE)
                TDFCR(ICT,2)=0.
                if(idebug .ge. 2)then
                    WRITE(6,948)ICT,N,TEMDIF(N-1),TEMDIF(N),SLOPE,FRAC
     #                         ,TDFCR(ICT,1)
                endif
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

                if(idebug .ge. 2)WRITE(6,951)pos_area_max
     1                       ,pos_area,PBENEG,nbe_min,PBECR(ICP,2),ICP

                PBECR(ICP,3)=PBECR(ICP,2)-PBECR(max(1,ICP-1),2)
                PBECR(ICP,4)=HT(N-1)+FRAC*DELTAH(N)
C
                if(idebug .ge. 2)then
                  WRITE(6,949)ICP,PBECR(ICP,1),PBECR(ICP,2),PBECR(ICP,3)
     &          ,PBECR(ICP,4),pos_area_max,nbe_min
 949              format(' PBECR',I3,4F10.3,2F10.2)
                endif

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

        if(idebug .ge. 2)then
          WRITE(6,485)pos_area_max,PBENEG,VELNEG
 485      format(' pos_area_max',F10.1,' PBENEG',F10.1,' VELNEG',F10.1)
        endif

        RETURN
        END
C
!
C       1991    Steve Albers
C
        FUNCTION WTBLB(P,TC,W,IOUT,istatus)

cdoc    Calculate Wet Bulb, given P,T,W

        THETAE=THAEK(P,TC,W)
        CALL MSAD5(WTBLB_arg,P,THETAE,TC,20.,SLOPE,I1,I2,IA,IOUT
     1            ,istatus)      

        if(istatus .ne. 1)then
            write(6,*)' Error in WTBLB',P,TC,W,IOUT
!           stop
        endif

        WTBLB = WTBLB_arg
        RETURN
        END
C
C
C
        SUBROUTINE BLAYR(P,T,Q,PMEAN,TMEAN,WMEAN,THKNES,NLEVEL,HH,
     +                    LOWEST,IO)

cdoc    Calculate boundary layer mean values from an input sounding

!       Steve Albers 1991

        include 'lapsparms.for'
        integer MXL
        parameter (MXL=MAX_LVLS+1) ! number of 3D levels plus the sfc level

        REAL INTLOG
        DIMENSION P(MXL),T(MXL),Q(MXL)
        TVIRT(TT,QQ)=TT/(1.-QQ*.37803)
        THICK(P1,P2,TC,QQ)=ALOG(P1/P2)*TVIRT((TC+273.15),QQ)*.09604

        IF(NLEVEL.LT.2)GOTO9000
        ZERO=0.
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

        if(pdelt .ne. 0.)then
            HH=HH+THICK(P(I-1),PBOUND,ARGT/PDELT,ARGQ/PDELT)
        endif

        IF(IO.GE.2)WRITE(6,1)I,P(I),T(I),Q(I),WT,SUMWT,SUMT,SUMQ,BT,BQ
     +              ,INTLOG,HH
! Hongli Jiang: W>=D+3 from F6.4 to F7.4 10/14/2013
 1      FORMAT(1x,I2,F6.0,F6.1,F7.4,2F7.1,2F11.4,F11.6,F9.6,F9.3,F6.2)
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

cdoc    Calculate LCL properties from an input parcel

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

cdoc    Calculate LCL properties from an input sounding (efficiently)

!       Steve Albers 1991

        ES(X)=6.1078+X*(.443652+X*(.014289+X*(2.65065E-4+X*
     1 (3.03124E-6+X*(2.034081E-8+X*(6.13682E-11))))))
        DATA EPSILN/.62197/,GAMMAI/102.4596/        ! Lapse Rate m/deg

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
        SUBROUTINE ITPLV(P,PARAM,NLEVEL,PINT,PARMAN,IO,istatus)

cdoc    Interpolate any parameter from a pressure sounding to a specific pres

!       Steve Albers 1991

        include 'lapsparms.for'
        integer MXL
        parameter (MXL=MAX_LVLS+1) ! number of 3D levels plus the sfc level

        DIMENSION P(MXL),PARAM(MXL)

        if(p(1) .lt. pint)then
            write(6,*)' Error in ITPLV: p(1) < pint',p(1),pint
            istatus = 0
            print *, p(:)  ! WNIDB
            return
        endif

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
 999    istatus = 1
        RETURN
        END
C
        SUBROUTINE NEWTN(X,XOLD,Y,YOLD,SLOPE,ITER,IO,FENCE,istatus)

cdoc    Newton iteration

!       Steve Albers 1991

        IF(ITER.GT.0)SLOPE=(Y-YOLD)/(X-XOLD)
        DELTA=-Y/SLOPE
        DELTA=AMAX1(-FENCE,AMIN1(DELTA,FENCE))
        XOLD=X
        X=X+DELTA
        ITER=ITER+1
        if(iter .gt. 100)then
            write(6,*)' Error in NEWTN: too many iterations',iter
            istatus = 0
            return
        endif

        IF(IO.GE.2)WRITE(6,1)ITER,X,XOLD,Y,YOLD,SLOPE,DELTA
 1      FORMAT(' NEWTON ITER',I4,2F11.4,2F11.7,2E11.3)
        YOLD=Y

        istatus = 1
        RETURN
        END
C
C
        SUBROUTINE MSAD5
     ^  (TEMNEW,PRESNW,THETAE,TGUESS,SLOPEG,SLOPE,I1,I2,IA,IO,istatus)       

cdoc    Calculate along a moist adiabat. Solve for T, given ThetaE and P

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

        if(I1 .gt. 200000)then
            write(6,*)' Error, too many I1 loops in MSAD5'
            istatus = 0
            return
        endif

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
        istatus = 1
        RETURN
        END
C
!
C      1991     Steve Albers
C
       FUNCTION THAE(TC,TD,P)

cdoc   Calculate Theta(e), given T, Td, P
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

cdoc   COMPUTES THE DRY AIR POTENTIAL TEMPERATURE

       real AK
       parameter (AK=.28613105)

       AKS=AK * (1.0+1.608*W)/(1.0 + 1.941569*W)
       THD=T * ((1000.0*PMEI)**AKS)
       RETURN
       END

       FUNCTION RL(TM2)
cdoc   LATENT HEAT OF EVAPORATION
C      TM2=T-273.15
C      RL=597.31-0.589533*TM2+0.001005333*(TM2*TM2)
       RL=597.31-((0.589533+0.001005333*TM2)*TM2)
       RETURN
       END
C
C
       SUBROUTINE DRYAD(P1,P2,T1,T2,W,IO)

!       Steve Albers 1991

cdoc   COMPUTES THE DRY AIR POTENTIAL TEMPERATURE (theta), GIVEN P and T
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

cdoc   Compute theta(e), given T, TD, Altimeter setting, and elevation (HGT)

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
cdoc   COMPUTE MIXING RATIO (GM/GM) GIVEN DEW POINT TEMP AND THE PRESSURE (MB)
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
cdoc   COMPUTES THE EQUIVALENT POTENTIAL TEMPURATURE (K).
cdoc   (USING THE ROSSBY DEFINITION)
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

cdoc    Quick way to get Theta(e) using lookup table

!       1991    Steve Albers

        character*31 ext
        character*150 directory
        integer len_dir

        integer pres2_low,pres2_high,pres2_interval,n_pres2
        parameter (pres2_low = 250)
        parameter (pres2_high = 1100)
        parameter (pres2_interval = 10)

        parameter (n_pres2 = (pres2_high - pres2_low) / pres2_interval)

        integer t_low,t_high,t_interval,n_t
        parameter (t_low  = -60)
        parameter (t_high = +60)
        parameter (t_interval = 2)

        parameter (n_t = (t_high - t_low) / t_interval)

        integer tdprs_low,tdprs_high,tdprs_interval,n_tdprs
        parameter (tdprs_low  = 0)
        parameter (tdprs_high = +30)
        parameter (tdprs_interval = 1)

        parameter (n_tdprs = (tdprs_high - tdprs_low) / tdprs_interval)

        real thetae_lut(0:n_t,0:n_tdprs,0:n_pres2)

        logical l_write_lut

        save init,thetae_lut
        data init/0/

        l_write_lut = .false.

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
            close(11)

            if(l_write_lut)then
                open(12,file=directory(1:len_dir)//'thetae_lut.dat'
     1                         ,form='unformatted',status='unknown')
            endif
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

            if(l_write_lut)then
                write(12)thetae_lut
                close(12)
            endif

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

cdoc    Quick way to get temperature along a moist adiabat given theta(e)
cdoc    This uses a lookup table

!       thetae in K  (Input)
!       pres in mb   (Input)
!       t moist in K (Output)

!       1991    Steve Albers

        character*31 ext
        character*150 directory
        integer len_dir

        integer thetae_low,thetae_high,thetae_interval,n_thetae
        parameter (thetae_low  = 210)
        parameter (thetae_high = 410)
        parameter (thetae_interval = 1)

        parameter (n_thetae = (thetae_high - thetae_low) / thetae_interv
     1al)

        integer pres_low,pres_high,pres_interval,n_pres
        parameter (pres_low = 50)
        parameter (pres_high = 1100)
        parameter (pres_interval = 5)

        parameter (n_pres = (pres_high - pres_low) / pres_interval)

        real t_moist(0:n_thetae,0:n_pres)

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
cdoc THIS FUNCTION RETURNS THE DEW POINT (CELSIUS) GIVEN THE TEMPERATURE
cdoc (CELSIUS) AND RELATIVE HUMIDITY (%).
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

C   This has been updated to test for out of bounds values (2004 Steve Albers)

        if(rh .lt. 0. .or. rh .gt. 100.)then
            call get_r_missing_data(r_missing_data,istatus)
            DWPT_laps = r_missing_data
            return
        endif

        if(t .lt. -100. .or. t .gt. +100.)then
            call get_r_missing_data(r_missing_data,istatus)
            DWPT_laps = r_missing_data
            return
        endif

        X = 1.-0.01*RH
C   COMPUTE DEW POINT DEPRESSION.
        DPD =(14.55+0.114*T)*X+((2.5+0.007*T)*X)**3+(15.9+0.117*T)*X**14
        DWPT_laps = T-DPD
        RETURN
        END


        function twet_fast(t_c,td_c,pres_mb)

!       Steve Albers 1991
cdoc    This is a fast approximate Wet Bulb routine using lookup tables
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

       subroutine li_laps(tc,td,pr,i4time,imax,jmax,li,flag
     !                   ,istatus)

cdoc   Compute 2-D grid of LI, given a grid of parcels to launch

!      1991     Steve Albers

       real tc(imax,jmax) ! Input T  in deg F
       real td(imax,jmax) ! Input Td in deg F
       real pr(imax,jmax) ! Input Pr in MB
       real li(imax,jmax) ! Output Li in Deg K/C

       real t500laps(imax,jmax) ! Used Locally Only

       character*13 filename13

       call get_r_missing_data(r_missing_data,istatus)
       if(istatus .ne. 1)then
           write(6,*)' Error reading r_missing_data'
           return
       endif

       call constant(li,r_missing_data,imax,jmax)
       call get_laps_cycle_time(laps_cycle_time,istatus)
!      Get 500 mb temp field
       if(flag .eq. 0.0)then ! Use LT1 (or equivalent) file
           write(6,*)' Reading 500 temp from LT1 (or equivalent) file'
           k_level = 500
           call get_temp_2d(i4time,laps_cycle_time,i4time_nearest
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

cdoc    Calculate LI given an input parcel

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


        subroutine thetae_to_t500(thetae,t500,istatus)

cdoc    Given theta(e), what is T-500mb?

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

cdoc    Calculate exp function with check to avoid underflow with large inputs.

        if(x .ge. -70.)then
            expm = exp(x)
        else
            expm = 0.
        endif

        end


        function devirt_td(t_k,td_k,p_pa)

cdoc    This function yields the approximate temperature given the virtual temp
cdoc    tv from the mthermo library is called. 

!       Please suggest improvements to this routine to Steve Albers at FSL.

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

cdoc    This function yields the approximate temperature given the virtual temp
cdoc    tv from the mthermo library is called. 

!       Please suggest improvements to this routine to Steve Albers at FSL.

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

cdoc    This function yields the approximate temperature given the virtual temp
cdoc    devirt_td from the laps library is called. 

!       Please suggest improvements to this routine to Steve Albers at FSL.

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
cdoc    THIS FUNCTION RETURNS THE VIRTUAL TEMPERATURE TV (CELSIUS) OF
cdoc    A PARCEL OF AIR AT TEMPERATURE T (CELSIUS), DEW POINT TD
cdoc    (CELSIUS), AND PRESSURE P (MILLIBARS). THE EQUATION APPEARS
cdoc    IN MOST STANDARD METEOROLOGICAL TEXTS.
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
cdoc  THIS FUNCTION RETURNS THE TEMPERATURE TSA (CELSIUS) ON A SATURATION
cdoc  ADIABAT AT PRESSURE P (MILLIBARS). OS IS THE EQUIVALENT POTENTIAL
cdoc  TEMPERATURE OF THE PARCEL (CELSIUS). SIGN(A,B) REPLACES THE
cdoc  ALGEBRAIC SIGN OF A WITH THAT OF B.
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

        subroutine get_tw_approx_2d(t_k,td_k,p_pa,ni,nj,tw_k)

cdoc    Calculate Wet Bulb, using a fast approximate method

!       Steve Albers 1991

!       This routine is fast but only accurate near 0 degrees C (273K)

        real t_k(ni,nj)     ! Input
        real td_k(ni,nj)    ! Input
        real p_pa(ni,nj)    ! Input
        real tw_k(ni,nj)    ! Output

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
            ratio = (t_k(i,j) - start)/100. * 
     1                                   (1.0 + 0.9 * (1.0 - sqrt(rh)))        
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

cdoc    Calculate 2-D grid of Wet Bulb, using 'tw_fast'

!       Steve Albers 1991
!       WARNING: This routine may not work because it calls tw_fast

        real t_k(ni,nj)     ! Input
        real td_k(ni,nj)    ! Input
        real p_pa(ni,nj)    ! Input
        real tw_k(ni,nj)    ! Output

        do j = 1,nj
        do i = 1,ni
            tw_k(i,j) = 
     1      tw_fast(t_k(i,j)-273.15,td_k(i,j)-273.15,p_pa(i,j)*.01)
     1                                                     + 273.15
        enddo ! i
        enddo ! j

        return
        end

        subroutine get_tw_2d_orig(t_k,td_k,p_pa,ni,nj,tw_k)

cdoc    Calculate 2-D grid of Wet Bulb, using 'tw'

!       Steve Albers 1991

        real t_k(ni,nj)     ! Input
        real td_k(ni,nj)    ! Input
        real p_pa(ni,nj)    ! Input
        real tw_k(ni,nj)    ! Output

        do j = 1,nj
        do i = 1,ni
            tw_k(i,j) = 
     1      tw(t_k(i,j)-273.15,td_k(i,j)-273.15,p_pa(i,j)*.01)
     1                                                     + 273.15
        enddo ! i
        enddo ! j

        return
        end

        FUNCTION TW_fast(T,TD,P)
C
C
cdoc    THIS FUNCTION RETURNS THE WET-BULB TEMPERATURE TW (CELSIUS)
cdoc    GIVEN THE TEMPERATURE T (CELSIUS), DEW POINT TD (CELSIUS)
cdoc    AND PRESSURE P (MB).  SEE P.13 IN STIPANUK (1973), REFERENCED
cdoc    ABOVE, FOR A DESCRIPTION OF THE TECHNIQUE.

cdoc    WARNING: This routine may not work because it calls TSA_fast
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
cdoc    THIS FUNCTION RETURNS THE MIXING RATIO (GRAMS OF WATER VAPOR PER
cdoc    KILOGRAM OF DRY AIR) GIVEN THE TEMPERATURE T (CELSIUS) AND PRESSURE
cdoc    (MILLIBARS). THE FORMULA IS QUOTED IN MOST METEOROLOGICAL TEXTS.
cdoc    Note this is a faster version done by Steve Albers.
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
cdoc    THIS FUNCTION RETURNS THE MIXING RATIO (GRAMS OF WATER VAPOR PER
cdoc    KILOGRAM OF DRY AIR) GIVEN THE TEMPERATURE T (CELSIUS) AND PRESSURE
cdoc    (MILLIBARS). THE FORMULA IS QUOTED IN MOST METEOROLOGICAL TEXTS.
cdoc    Note this is a faster version done by Steve Albers
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
cdoc    THIS FUNCTION RETURNS THE EQUIVALENT POTENTIAL TEMPERATURE OS
cdoc    (K) FOR A PARCEL OF AIR SATURATED AT TEMPERATURE T (K)
cdoc    AND PRESSURE P (MILLIBARS).
cdoc    Note this is a faster version done by Steve Albers

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


        subroutine wtblb_lvl(twet_c,P,T,Q,WB,MXL,NLEVEL,PWB0,HWB0)

        real twet_c         ! Input:  wet bulb temperature in degrees C 
        real P(MXL)         ! Input:  pressure sounding (mb)
        real T(MXL)         ! Input:  temperature sounding (C)
        real Q(MXL)         ! Input:  specific humidity sounding
        real WB(MXL)        ! Input:  wet bulb temperature sounding
        integer MXL         ! Input:  size of P,T,Q arrays
        integer nlevel      ! Input:  number of levels in vertical arrays
        real PWB0           ! Output: pressure of the wet bulb level
        real HWB0           ! Output: height of the wet bulb level (kft agl)

        IO = 0
 	IOUT=MIN(IO,1)                                                          
 	IF(WB(1).GE.twet_c)GOTO150                                                  
!	IF(IO.GE.2)WRITE(6,87)                                                  
 87	format(' SURFACE WETBULB TEMP IS BELOW twet_c')                          
 	HWB0=0.                                                                 
 	PWB0=0.                                                                 
 	GOTO390                                                                 
C                                                                         
!       Test for bracketing of the Wet Bulb Zero
 150	DO 200 N=2,NLEVEL                                                    
! 	    IF(WB(N)*WB(N-1))250,250,200                                            
            if( (wb(n) - twet_c) * (wb(n-1) - twet_c) )250,250,200 
 200    CONTINUE                                                           

        write(6,*)' Warning, no wet bulb detected in loop'
        goto390

 250    frac = (twet_c - WB(N-1)) / (WB(N) - WB(N-1))
        PWB0 = P(N-1) * (1. - frac) + P(N) * frac

!       Integrate hydrostatically to get height from pressure
 	CALL BLAYR(P,T,Q,DUM1,DUM2,DUM3,P(1)-PWB0,NLEVEL,HWB0,1,IO) 
!	IF(IO.GE.1)WRITE(6,351)PWB0,HWB0                            
 351	format(' WETBULB ZERO IS AT',F6.1,'MB     OR',F6.2,'KFT  AGL')        

        goto400                      ! Normal condition

 390    HWB0 = r_missing_data        ! Indeterminate condition
        PWB0 = r_missing_data

 400    CONTINUE

        return
        end
