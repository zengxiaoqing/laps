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

        subroutine get_heights_hydrostatic
     1              (temp_3d_k          ! Input (3d temp K)
     1              ,pres_sfc_pa        ! Input (surface pressure pa)
     1              ,pres_3d            ! Input (3d pressure pa)
     1              ,sh_3d              ! Input (3d specific humidity)
     1              ,topo               ! Input (terrain m)
     1              ,ni,nj,nk           ! Input (Dimensions)
     1              ,heights_3d         ! Output (m)
     1              ,istatus)           ! Output

!           1991        Steve Albers
!       Nov 1991        Steve Albers    Change to Sfc P Input + Terrain
!           1996        Steve Albers    Expanded 'esat_lut' array from -100C
!                                       to the colder value of -120
!
!       For each horizontal gridpoint, the 'pres_sfc_pa' at the lowest level 
!       must be greater than or equal to the reference (surface) pressure.

        real esat_lut(-120:+100)

        real temp_3d_k(ni,nj,nk)
        real pres_sfc_pa(ni,nj)
        real topo(ni,nj)
        real heights_3d(ni,nj,nk)
        real sh_3d(ni,nj,nk)

        real heights_below(ni,nj)
        real a_below(ni,nj)
        real pres_sfc_mb(ni,nj)
        real z_correction(ni,nj)

        real pres_3d(ni,nj,nk)                    

        real pres_3d_mb(ni,nj,nk)               ! Local
        real p_1d_mb(nk)                        ! Local
        real alog_array(ni,nj,nk)               ! Local

        real make_td, k_to_c

        include 'constants.inc'
        real C1
        PARAMETER (C1 = EP_1) 

        real C2
        PARAMETER (C2 = r_d / (2. * grav) )

cdoc    Generate heights of grid points using laps data and hydrostatic equation
        write(6,*)' Generating 3D height grid'
        write(6,*)' Integrating hydrostatic equation over the LAPS grid'
c       write(6,*)' Initialize and calculate first level'

        do it = -120,+100
            esat_lut(it) = esat(float(it))
        enddo

        call get_r_missing_data(r_missing_data,istatus)
        z_correction = r_missing_data

        do k = 1,nk
            do j = 1,nj
            do i = 1,ni
!               p_1d_mb(k) = zcoord_of_level(k)/100.
                pres_3d_mb(i,j,k) = pres_3d(i,j,k) / 100.
                if(k .gt. 1)then
                    alog_array(i,j,k-1) = 
     1                  alog(pres_3d_mb(i,j,k-1)/pres_3d_mb(i,j,k))
                endif
            enddo ! i
            enddo ! j
        enddo


        do j = 1,nj
        do i = 1,ni
            pres_sfc_mb(i,j) = pres_sfc_pa(i,j) / 100.
            heights_below(i,j) = 0.
            t_2d_c_below = k_to_c(temp_3d_k(i,j,1))
            sh_below = sh_3d(i,j,1)
            w_below = sh_below / (1. - sh_below)
            A_below(i,j)= temp_3d_k(i,j,1) * (1. + C1 * w_below)
        enddo ! i
        enddo ! j

        do k = 2,nk ! Integrate upward one level at a time.

cCCCCCCCCCCCCCCCCCC                         ISTAT = LIB$SHOW_TIMER(,,,)

            do j = 1,nj
            do i = 1,ni
                sh_value = sh_3d(i,j,k)
                w_value = sh_value / (1. - sh_value)

                A1= temp_3d_k(i,j,k) * (1. + C1 * w_value)

                Z_add = C2 * (A1+A_below(i,j)) * (alog_array(i,j,k-1))       

                heights_3d(i,j,k) = heights_below(i,j) + z_add

!               Test whether this layer contains the surface
                if(pres_3d_mb(i,j,k  ) .le. pres_sfc_mb(i,j) .and.
     1             pres_3d_mb(i,j,k-1) .ge. pres_sfc_mb(i,j)      )then       

                    alog_factor = alog( pres_3d_mb(i,j,k-1)
     1                                 /pres_sfc_mb(i,j)   )
                    z_to_sfc = z_add * (  alog_factor 
     1                                  / alog_array(i,j,k-1) )

                    z_below = topo(i,j) - z_to_sfc

!                   z_below = topo(i,j) -
!       1           z_thk(pres_sfc_mb(i,j),pres_3d_mb(i,j,k-1),t_1d_c,td_1d_c
!       1                                               ,alog_array,esat_lut,2)

                    z_correction(i,j) = z_below - heights_below(i,j)
                endif

                heights_below(i,j) = heights_3d(i,j,k)
                a_below(i,j) = a1

            enddo ! i
            enddo ! j

c           write(6,*)' Completed level',k

        enddo ! k

cCCCCCCCCCCCCCCCCCC                         ISTAT = LIB$SHOW_TIMER(,,,)

        ndiffp = 0
        diffp_max = -abs(r_missing_data)
        diffp_min = +abs(r_missing_data)
        iwrite = 0

      ! Apply Constant of Integration to the heights
        do j = 1,nj
        do i = 1,ni
            z_correction_orig = z_correction(i,j)
            if(z_correction(i,j) .eq. r_missing_data)then
                z_correction(i,j) = 0. ! Level 1 becomes the reference height
                ndiffp = ndiffp + 1
            endif

            heights_3d(i,j,1) = z_correction(i,j)

!           Obtain diffp information
            diffp = pres_sfc_mb(i,j) - pres_3d_mb(i,j,1) 
            diffp_max = max(diffp,diffp_max)
            diffp_min = min(diffp,diffp_min)
            iwrite = iwrite + 1
            if(iwrite .le. 10)then
                write(6,*)i,j,pres_sfc_mb(i,j),pres_3d_mb(i,j,1)
     1               ,diffp,diffp_max,diffp_min,ndiffp,z_correction_orig       
            endif
        enddo ! i
        enddo ! j

        do k = 2,nk
            do j = 1,nj
            do i = 1,ni
                heights_3d(i,j,k) = heights_3d(i,j,k) + z_correction(i,j
     1)
            enddo ! j
            enddo ! i
        enddo ! k

c       write(6,*)' Added in constant of integration to the heights'
cCCCCCCCCCCCCCCCCCC                         ISTAT = LIB$SHOW_TIMER(,,,)

        write(6,*)' Range of diffp ( Psfc - P[1] ) = ',diffp_min
     1                                                ,diffp_max,'mb'

!       Return with diagnostic information
        ngrid = ni*nj
        if(ndiffp .eq. ngrid)then
            write(6,*)' WARNING: surface extends outside '
     1                ,'3D grid at all grid points'    
            write(6,*)' Output heights are relative to level 1'
            istatus = -1
            return
        elseif(ndiffp .gt. 0)then
            write(6,*)' ERROR: surface extends outside '
     1                ,'3D grid at some grid points',ndiffp    
            istatus = 0
            return
        else
            write(6,*)' Success in get_heights_hydrostatic'
            istatus = 1
            return
        endif

        end


        FUNCTION Z_thk(PT,P,T,TD,alog_array,esat_lut,N)
C
cdoc    THIS FUNCTION RETURNS THE THICKNESS OF A LAYER BOUNDED BY PRESSURE
cdoc    P(1) AT THE BOTTOM AND PRESSURE PT AT THE TOP.
C
C       BAKER,SCHLATTER 17-MAY-1982     Original version
C       Albers                 1990     Restructured
C       Albers                 1996     Expanded 'esat_lut' array from
C                                       -100C to the colder value of -120
C
C   ON INPUT:
C       P = PRESSURE (MB).  NOTE THAT P(I).GT.P(I+1).
C       T = TEMPERATURE (CELSIUS)
C       TD = DEW POINT (CELSIUS)
C       N = NUMBER OF LEVELS IN THE SOUNDING AND THE DIMENSION OF
C           P, T AND TD
C   ON OUTPUT:
C       Z = GEOMETRIC THICKNESS OF THE LAYER (M)
C
C   THE ALGORITHM INVOLVES NUMERICAL INTEGRATION OF THE HYDROSTATIC
C   EQUATION FROM P(1) TO PT. IT IS DESCRIBED ON P.15 OF STIPANUK
C   (1973).
C
        DIMENSION T(N),P(N),TD(N),TK(100),alog_array(N)
C       C1 = .001*(1./EPS-1.) WHERE EPS = .62197 IS THE RATIO OF THE
C                             MOLECULAR WEIGHT OF WATER TO THAT OF
C                             DRY AIR. THE FACTOR 1000. CONVERTS THE
C                             MIXING RATIO W FROM G/KG TO A DIMENSION-
C                             LESS RATIO.
C       C2 = R/(2.*G) WHERE R IS THE GAS CONSTANT FOR DRY AIR
C                     (287 KG/JOULE/DEG K) AND G IS THE ACCELERATION
C                     DUE TO THE EARTH'S GRAVITY (9.8 M/S**2). THE
C                     FACTOR OF 2 IS USED IN AVERAGING TWO VIRTUAL
C                     TEMPERATURES.
        real esat_lut(-120:+100)

        include 'constants.inc' 
        real C1,C2
        PARAMETER (C1 = .001 * EP_1)
        PARAMETER (C2 = r_d / (2. * grav) )

        DO 5 I= 1,N
           TK(I)= c_to_k(T(I))
    5   CONTINUE

        do i = 1,N-1
            j = I + 1

            IF (P(J) .gt. PT)then             ! We are still in complete layers
                A1= TK(J)*(C2 + W_laps(TD(J),P(J),esat_lut))
                A2= TK(I)*(C2 + W_laps(TD(I),P(I),esat_lut))
                Z_thk = Z_thk + (A1+A2)*(alog_array(i))

            elseif(P(J) .eq. PT)then          ! Finish up on Complete Layer
                A1= TK(J)*(C2 + W_laps(TD(J),P(J),esat_lut))
                A2= TK(I)*(C2 + W_laps(TD(I),P(I),esat_lut))
                Z_thk = Z_thk + (A1+A2)*(alog_array(i))
                return

            else ! P(J) .lt. PT               ! Finish up on partial layer
                A1= TK(J)*(C2 + W_laps(TD(J),P(J),esat_lut))
                A2= TK(I)*(C2 + W_laps(TD(I),P(I),esat_lut))
                Z_thk = Z_thk + (A1+A2)*(ALOG(P(I)/PT))
                RETURN

            endif

        enddo ! i


        write(6,*)' Error, PT out of bounds in Z_THK'
        Z_thk = -1.0

        RETURN
        END


        FUNCTION W_laps(T,P,esat_lut)

cdoc    Convert T(C) and P to W. This is a fast approximate routine.
cdoc    Output W_laps is dimensionless mixing ratio

!       This function really only works when T > -50C but the efficiency will
!       outweigh the error when T < -50C in this application

        real esat_lut(-120:+100)

        include 'constants.inc'
        real C1,C2,const
        PARAMETER (C1 = .001 * EP_1)
        PARAMETER (C2 = r_d / (2. * grav) )

        parameter (const = 622. * c1 * c2)

C       ES(X)=6.1078+X*(.443652+X*(.014289+X*(2.65065E-4+X*
C    1 (3.03124E-6+X*(2.034081E-8+X*(6.13682E-11))))))
C
C
        X= esat_lut(nint(T))
        W_laps= const*X/(P-X)
        RETURN
        END


