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

        SUBROUTINE FFLXC(NI,NJ,M,SCALE,uh,vh,field,flxcnv,lat,lon 
     1  ,flu,flv,sigma,r_missing_data)

cdoc    Calculates convergence given a wind field rotated to the map projection

!       Original Version S. Albers                                      1986
!       Made more general                                               1989
!       Made even more general (Adjustable Dims)                        1991
!       Added new map projections for sigma calc                        1997
!       Input M is changed, hopefully balances changes in get_sigma     1998

        include 'trigd.inc'
        REAL M
        DIMENSION FIELD(ni,nj),FLU(ni,nj),FLV(ni,nj),SIGMA(ni,nj)
        real UH(ni,nj),VH(ni,nj),FLXCNV(ni,nj),lat(ni,nj),lon(ni,nj)
        DATA  RPD/.0174532925/

c       WRITE(6,18)
c18     FORMAT('        CALCULATING FLUX CONVERGENCE FIELD')

C
C  CALCULATE FIELD FLUX AND INITIALIZE ARRAYS
C  UNITS FOR FIELD FLUX CONVERGENCE ARE 1/SCALE * S**-1
C      COEFF=-1/
C      (SCALE/(2*(DELTA*.0254MPIN/M) METERS/GRIDPOINT) * .5147818MSPKT)
!       COEFF=-1./(SCALE/(2.*(DELTA*.0254/M))*.5147818)
        COEFF=-1./(SCALE/(2./M))
c       write(6,*)' coeff = ',coeff
        NIM1=NI-1
        NIM2=NI-2
        NJM1=NJ-1
        NJM2=NJ-2
C
C       CALCULATE SIGMA
C       THIS VALUE REPRESENTS  (M*RHO(EARTH)*(1.+SIN(PHI0))/MPIN)**2

        DO 100 I=1,NI
        DO 100 J=1,NJ
!           SIGMA(I,J) = (1.+sind(PHI0)) / (1.+ sind(lat(i,j))) 
!           sigma(i,j)=1.0

            call get_sigma(lat(i,j),lon(i,j),sigma(i,j),istatus)
            if(istatus .ne. 1)stop
C
            FL=FIELD(I,J)/(COEFF*SIGMA(I,J))
            FLU(I,J)=FL* UH(I,J)
            FLV(I,J)=FL* VH(I,J)
            FLXCNV(I,J)=0.
100     continue

C
C
        DFLUDX=-3.*FLU(   1,   1)+4.*FLU(   2,   1)-FLU(   3,   1)
        DFLVDY=-3.*FLV(   1,   1)+4.*FLV(   1,   2)-FLV(   1,   3)
 1000   FLXCNV(1,1)=(DFLUDX+DFLVDY)*SIGMA(1,1)**2
C
        DO 2000 J=2,NJM1
        DFLUDX=-3.*FLU(   1,   J)+4.*FLU(   2,   J)-FLU(   3,   J)
        DFLVDY=    FLV(   1, J+1)  - FLV(   1, J-1)
 2000   FLXCNV(1,J)=(DFLUDX+DFLVDY)*SIGMA(1,J)**2
C
        DFLUDX=-3.*FLU(   1,  NJ)+4.*FLU(   2,  NJ)-FLU(   3,  NJ)
        DFLVDY= 3.*FLV(   1,  NJ)-4.*FLV(   1,NJM1)+FLV(   1,NJM2)
 3000   FLXCNV(1,NJ)=(DFLUDX+DFLVDY)*SIGMA(1,NJ)**2
C
        DO 4000 I=2,NIM1
        DFLUDX=    FLU( I+1,   1)  - FLU( I-1,   1)
        DFLVDY=-3.*FLV(   I,   1)+4.*FLV(   I,   2)-FLV(   I,   3)
 4000   FLXCNV(I,1)=(DFLUDX+DFLVDY)*SIGMA(I,1)**2
C
        DO 5000 I=2,NIM1
        IP1=I+1
        IM1=I-1

        DO 5000 J=2,NJM1
        DFLUDX=    FLU( IP1,   J)  - FLU( IM1,   J)
        DFLVDY=    FLV(   I, J+1)  - FLV(   I, J-1)
        FLXCNV(I,J)=(DFLUDX+DFLVDY)*SIGMA(I,J)**2
C       CALL OUTPTF(I,J,FLU,FLV,DFLUDX,DFLVDY,FLXCNV,UH,VH,SIGMA
C      .      ,NIM2,NJM2,FIELD,SCALE,COEFF)
 5000   CONTINUE
C
        DO 6000 I=2,NIM1
        DFLUDX=    FLU( I+1,  NJ)  - FLU( I-1,  NJ)
        DFLVDY= 3.*FLV(   I,  NJ)-4.*FLV(   I,NJM1)+FLV(   I,NJM2)
 6000   FLXCNV(I,NJ)=(DFLUDX+DFLVDY)*SIGMA(I,NJ)**2
C
        DFLUDX= 3.*FLU(  NI,   1)-4.*FLU(NIM1,   1)+FLU(NIM2,   1)
        DFLVDY=-3.*FLV(  NI,   1)+4.*FLV(  NI,   2)-FLV(  NI,   3)
 7000   FLXCNV(NI,1)=(DFLUDX+DFLVDY)*SIGMA(NI,1)**2
C
        DO 8000 J=2,NJM1
        DFLUDX= 3.*FLU(  NI,   J)-4.*FLU(NIM1,   J)+FLU(NIM2,   J)
        DFLVDY=    FLV(  NI, J+1)  - FLV(  NI, J-1)
 8000   FLXCNV(NI,J)=(DFLUDX+DFLVDY)*SIGMA(NI,J)**2
C
        DFLUDX= 3.*FLU(  NI,  NJ)-4.*FLU(NIM1,  NJ)+FLU(NIM2,  NJ)
        DFLVDY= 3.*FLV(  NI,  NJ)-4.*FLV(  NI,NJM1)+FLV(  NI,NJM2)
 9000   FLXCNV(NI,NJ)=(DFLUDX+DFLVDY)*SIGMA(NI,NJ)**2
        RETURN
        END

        subroutine get_sigma(rlat,rlon,sigma,istatus)

!       Steve Albers

cdoc    Sigma (map factor) is defined to be one when we are located in the 
cdoc    intersection of the projection plane and the earth's surface. It varies
cdoc    from unity elsewhere.

!       Equations from Principles of Meteorological Analysis, Walter Saucier
!       Pages 32,33

        include 'trigd.inc'
        real n
        character*6 c6_maproj

        integer init
        save init
        data init/0/

        save slat1,slat2,slon,c6_maproj

        if(init .eq. 0)then
            call get_standard_latitudes(slat1,slat2,istatus)
            if(istatus .ne. 1)then
                write(6,*)'get_sigma: bad istatus'
                return
            endif

            call get_standard_longitude(slon,istatus)
            if(istatus .ne. 1)then
                write(6,*)'get_sigma: bad istatus'
                return
            endif

            call get_c6_maproj(c6_maproj,istatus)
            if(istatus .ne. 1)then
                write(6,*)'get_sigma: bad istatus'
                return
            endif

            init = 1

        endif

        colat0 = 90. - slat1
        colat1 = 90. - slat1
        colat2 = 90. - slat2
        colat  = 90. - rlat

        if(c6_maproj .eq. 'plrstr')then ! polar stereo

!           Rotate to arbitrary pole
            polat = slat2
            polon = slon
            rlat_ge = rlat
            rlon_ge = rlon
            call GETOPS(rlat_ps,rlon_ps,rlat_ge,rlon_ge,polat,polon)       
!                          O       O        I      I      I     I

            call get_grid_spacing(grid_spacing_m,istatus)
            if(istatus .ne. 1)then
                write(6,*)'get_sigma: bad istatus'
                return
            endif

            call get_ps_parms(slat1,slat2,grid_spacing_m,phi0
     1                       ,grid_spacing_proj_m)

!           phi0 = slat1
            phi  = rlat_ps

!           Check to see if you are at the opposite pole
            if(phi .eq. -90.)then      
                sigma = 0.
                istatus = 0
                write(6,*)'get_sigma: sigma undefined at opposite pole'       
                return
            endif

            sigma = (1. + sind(phi0)) / (1. + sind(phi))               ! eq. 13

        elseif(c6_maproj .eq. 'lambrt')then ! lambert

!           Check to see if you are at a pole
            if(abs(rlat) .eq. 90.)then      
                sigma = 0.
                istatus = 0
                write(6,*)'get_sigma: sigma undefined at the pole'
                return
            endif

            if(slat1 .eq. slat2)then
                n = cosd(colat0)
                arg   =  tand(colat/2.) / tand(colat0/2.)
                sigma = (sind(colat0)   / sind(colat)   ) * arg**n     ! eq. 3

            else

                n = alog(sind(colat1)   /sind(colat2)   ) /            ! eq. 9
     1              alog(tand(colat1/2.)/tand(colat2/2.))

                arg   =  tand(colat/2.) / tand(colat1/2.)
                sigma = (sind(colat1)   / sind(colat)   ) * arg**n     ! eq. 10

            endif

        elseif(c6_maproj .eq. 'merctr')then ! mercator

!           Check to see if you are at a pole
            if(abs(rlat) .eq. 90.)then      
                sigma = 0.
                istatus = 0
                write(6,*)'get_sigma: sigma undefined at the pole'
                return
            endif

            sigma = sind(colat1) / sind(colat)                         ! eq. 11

        else
            write(6,*)' Invalid map projection in get_sigma: ',c6_maproj
            istatus = 0
            return

        endif

        istatus = 1
        return
        end
