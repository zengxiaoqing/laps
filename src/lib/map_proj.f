



C
C********************** MAP PROJECTION ROUTINES ****************************
C
C
C
        SUBROUTINE XYTOPS(X,Y,PLA,PLO,ERAD)
C
C     This convert x,y-polar stereographic coordinates to
C     lat/lon values.
C     longitude:   0 - 180  ; Greenwich and east, 0 - -180 Greenwich and west
C     latitude : -90 -  90  ; positive for northern hemisphere
C     It is assumed that the x-axis point towards the east why the
C     longitude is rotated relative to the 'standard pol.ste.' location
C            
      PI180=3.14159265/180.
      VDIST = ERAD*2.0
C
C     calculate distance from (0,0) and define (0,0) as 90,0 (90,-90 in
C     the rotated system)
C
      DIST=SQRT(X**2+Y**2)
      IF(DIST.EQ.0) THEN
         PLA= 90.0
         PLO=-90.0
      ELSE
C
C     calculate the latitude by means of atan
C
         PLA=ATAN(DIST/VDIST)/PI180
         PLA=90.0-2.0*PLA
C
C     calculate the longitude taking the directions into account
C
         IF(X.EQ.0.0) THEN
            IF(Y.GT.0.0) THEN
               PLO= 90.0
            ELSE
               PLO=-90.0
            END IF
         ELSE
            IF(X.GT.0.0) THEN
               PLO=ATAN(Y/X)/PI180
            ELSE
               PLO=ATAN(Y/X)/PI180+180.0
            END IF
         END IF
      END IF
C
C     rotate the longitude
C
      PLO=AMOD(PLO+450.0,360.0)
      RETURN
      END
C
C     *******************************************************************
C
      SUBROUTINE PSTOGE(PLA,PLO,GLAT,GLON,RLAT,WLON1)
C
C     Convert polar stereographic coordinates to geographical lat/lon
C     ditto with the pol.ste. pole at rlat,wlon1
C     
C     longitude:   0 ; 360 positive east (on input)
C               -180 ; 180 positive east (on output)
C     latitude : -90 ;  90 posive on northern hemisphere
C     It is assumed that the polar stereographic coordinates have been
C     rotated to the standard format with 0 degrees longitude along wlon1
C
C     TSP 21 JUNE 89
C
C     set flag for n/s hemisphere
C
      PI180=3.14159265/180.
C      
      IF(RLAT.GE.0.0) THEN
         HSIGN= 1.0
      ELSE
         HSIGN=-1.0
      END IF
C
C     test for a n/s pole case
C
      IF(RLAT.EQ.90.0) THEN
	 GLAT=PLA
         GLON=MOD(PLO+WLON1,360.0)
         GO TO 2000
      END IF
      IF(RLAT.EQ.-90.0) THEN
         GLAT=-PLA
         GLON=MOD(PLO+WLON1,360.0)
         GO TO 2000
      END IF
C
C     test for longitude on 'greenwich or date line'
C
      IF(PLO.EQ.0) THEN
         GLAT=RLAT-90.0+PLA
         IF(GLAT.LT.-90.0) THEN
            GLAT=-180.0-GLAT
            GLON=MOD(WLON1+180.0,360.0)
         ELSE
            GLON=WLON1
         END IF
         GO TO 2000
      END IF      
      IF(PLO.EQ.180.0) THEN
         GLAT=RLAT+90.0-PLA
         IF(GLAT.GT.90.0) THEN
            GLAT=180.0-GLAT
            GLON=MOD(WLON1+180.0,360.0)
         ELSE
            GLON=WLON1
         END IF
         GO TO 2000         
      END IF
C
C     Determine longitude distance relative to wlon1 so it belongs to
C     the absolute interval 0 - 180
C
      ARGU1=PLO
      IF(PLO.GT.180.0) ARGU1 = PLO-360.0
C
C     Get the latitude, the help circle BB and the longitude by first
C     calculating the argument and legalize it - then take the inverse fct.
C
      IF(HSIGN.GT.0.0) THEN
         ARG2A = SIN(PLA*PI180)*SIN(HSIGN*RLAT*PI180)+
     +        COS(PLA*PI180)*COS(RLAT*PI180)*COS((180.0-ARGU1)*PI180)
      ELSE
         ARG2A = SIN(PLA*PI180)*SIN(HSIGN*RLAT*PI180)+
     +        COS(PLA*PI180)*COS(RLAT*PI180)*COS(ARGU1*PI180)
      END IF
      ARG2A = MIN(ARG2A, 1.0)
      ARG2A = MAX(ARG2A,-1.0)
      GLAT  = HSIGN*ASIN(ARG2A)
C
      IF(HSIGN.GT.0.0) THEN
         ARG2A = COS(RLAT*PI180)*SIN(PLA*PI180)+
     +        SIN(RLAT*PI180)*COS(PLA*PI180)*COS(ARGU1*PI180)
      ELSE
         ARG2A = COS(RLAT*PI180)*SIN(PLA*PI180)+
     +       SIN(-RLAT*PI180)*COS(PLA*PI180)*COS((180.0-ARGU1)*PI180)
      END IF
      ARG2A = MIN(ARG2A, 1.0)
      ARG2A = MAX(ARG2A,-1.0)      
      BB    = ACOS(ARG2A)
C
!     ARG2A = COS(GLAT)*COS(BB)/(1.0-SIN(GLAT)**2)
      ARG2A = COS(BB)/COS(GLAT)              ! 1997 Steve Albers
      ARG2A = MIN(ARG2A, 1.0)
      ARG2A = MAX(ARG2A,-1.0)      
      GLON  = ACOS(ARG2A)
C     
C     convert the radians to degrees 
C
        GLAT = GLAT/PI180
        GLON = GLON/PI180
C
C       the solution is symmetric so the direction must be if'ed
C
        IF(ARGU1.LT.0.0) THEN
           GLON = 360.0-GLON
        END IF
        GLON=AMOD(GLON+WLON1,360.0)
C
 2000 CONTINUE
C
C     the resultant longitude must be in the interval from -180, 180
C      
      IF(GLON.GT.180.0) GLON=GLON-360.0
      RETURN
      END




C     ****************************************************************
C
      SUBROUTINE GETOPS(PLA_r4,PLO_r4,GLAT_r4,GLON_r4,RLAT_r4,WLON1_r4)       
C                        Out     Out    In      In      In       In
C
C     Convert geographical lat/lon coordinates to polar stereographic
C     ditto with the pol.ste. pole at RLAT,WLON1
C     
C     longitude:-180 ; 180 positive east (on input)
C              :   0 ; 360 positive east (on output)
C     latitude : -90 ;  90 posive on northern hemisphere
C     The result is rotated 270 degrees relative to 'standard pol.ste.'
C     WLON1 is defined in the same way as the input
C     approach so as to get the x-axis to point towards the east, and the
C     y-axis towards the north along 0 degrees (at NP south along 180)
C
C     TSP 20/06-89
C
C     constants
C
      implicit real*8 (a-h,o-z)

      real pla_r4,plo_r4,glat_r4,glon_r4,rlat_r4,wlon1_r4

      PLA = pla_r4
      PLO = plo_r4
      GLAT = glat_r4
      GLON = glon_r4
      RLAT = rlat_r4
      WLON1 = wlon1_r4

      PI180 = 3.1415926535897932/180.0
C
C     Set flag for N/S hemisphere and convert longitude to <0 ; 360> interval
C
      IF(RLAT.GE.0.0) THEN
         HSIGN= 1.0
      ELSE
         HSIGN=-1.0
      END IF
      GLOR=GLON
      IF(GLOR.LT.0.0) GLOR=360.0+GLOR
      RWLON1=WLON1
      IF(RWLON1.LT.0.0) RWLON1=360.0+WLON1
C
C     Test for a N/S pole case
C
      IF(RLAT.EQ.90.0) THEN
         PLA=GLAT
         PLO=DMOD(GLOR+270.0D0,360.0D0)
         GO TO 2000
      END IF
      IF(RLAT.EQ.-90.0) THEN
         PLA=-GLAT
         PLO=DMOD(GLOR+270.0D0,360.0D0)
         GO TO 2000
      END IF
C
C     Test for GE coordinates at PS pole (Steve Albers - 1998)
C
      if(GLAT .eq. RLAT .and. GLON .eq. WLON1)then
         PLA = 90.0
         PLO = 0.
         GOTO2000
      endif
C
C     Test for longitude on 'Greenwich or date line'
C
      IF(GLOR.EQ.RWLON1) THEN
         IF(GLAT.GT.RLAT) THEN
            PLA=90.0-GLAT+RLAT
            PLO=90.0
         ELSE
            PLA=90.0-RLAT+GLAT
            PLO=270.0
         END IF
         GO TO 2000
      END IF      
      IF(DMOD(GLOR+180.0D0,360.0D0).EQ.RWLON1) THEN
         PLA=RLAT-90.0+GLAT
         IF(PLA.LT.-90.0) THEN
            PLA=-180.0-PLA
            PLO=270.0
         ELSE
            PLO= 90.0
         END IF
         GO TO 2000         
      END IF
C
C     Determine longitude distance relative to RWLON1 so it belongs to
C     the absolute interval 0 - 180
C
      ARGU1 = GLOR-RWLON1
      IF(ARGU1.GT. 180.0) ARGU1 = ARGU1-360.0
      IF(ARGU1.LT.-180.0) ARGU1 = ARGU1+360.0
C
C     1. Get the help circle BB and angle ALPHA (legalize arguments)
C
      ARG2A = COS(GLAT*PI180)*COS(ARGU1*PI180)
      ARG2A = MAX(ARG2A,-1.0D0)
      ARG2A = MIN(ARG2A, 1.0D0)         
      BB    = ACOS(ARG2A)
C
      ARG2A = SIN(HSIGN*GLAT*PI180)/SIN(BB)
      ARG2A = MAX(ARG2A,-1.0D0)
      ARG2A = MIN(ARG2A, 1.0D0)
      ALPHA = ASIN(ARG2A)
C
C     2. Get PLA and PLO (still legalizing arguments)
C
      ARG2A = COS(RLAT*PI180)*COS(BB)+
     +        SIN(HSIGN*RLAT*PI180)*SIN(HSIGN*GLAT*PI180)
      ARG2A = MAX(ARG2A,-1.0D0)
      ARG2A = MIN(ARG2A, 1.0D0)         
      PLA   = ASIN(ARG2A)
C
      ARG2A = SIN(BB)*COS(ALPHA)/COS(PLA)
      ARG2A = MAX(ARG2A,-1.0D0)
      ARG2A = MIN(ARG2A, 1.0D0)
      PLO   = ASIN(ARG2A)
C
C    Test for passage of the 90 degree longitude (duallity in PLO)
C         Get PLA for which PLO=90 when GLAT is the latitude
C
      ARG2A = SIN(HSIGN*GLAT*PI180)/SIN(HSIGN*RLAT*PI180)
      ARG2A = MAX(ARG2A,-1.0D0)
      ARG2A = MIN(ARG2A, 1.0D0)         
      PLA90 = ASIN(ARG2A)
C
C         Get help arc BB and angle ALPHA
C
      ARG2A = COS(RLAT*PI180)*SIN(PLA90)
      ARG2A = MAX(ARG2A,-1.0D0)
      ARG2A = MIN(ARG2A, 1.0D0)
      BB    = ACOS(ARG2A)

      ARG2A = SIN(HSIGN*GLAT*PI180)/SIN(BB)
      ARG2A = MAX(ARG2A,-1.0D0)
      ARG2A = MIN(ARG2A, 1.0D0)        
      ALPHA = ASIN(ARG2A)
C
C         Get GLOLIM - it is nesc. to test for the existence of solution
C
      ARGU2  = COS(GLAT*PI180)*COS(BB)/
     +            (1.-SIN(HSIGN*GLAT*PI180)*SIN(BB)*SIN(ALPHA))
      IF(ABS(ARGU2).GT.1.0) THEN
      GLOLIM = 999.0
      ELSE
        GLOLIM = ACOS(ARGU2)/PI180
      END IF
C
C     Modify (if nesc.) the PLO solution
C
      IF((ABS(ARGU1).GT.GLOLIM.AND.GLAT.LE.RLAT).OR.
     +   GLAT.GT.RLAT) THEN
            PLO = PI180*180.0 - PLO
      END IF
C
C     The solution is symmetric so the direction must be if'ed
C
      IF(ARGU1.LT.0.0) THEN
         PLO = -PLO
      END IF
C
C     Convert the radians to degrees
C
      PLA = PLA/PI180        
      PLO = PLO/PI180
C
C     To obtain a rotated value (ie so x-axis in pol.ste. points east)
C     add 270 to longitude
C
      PLO=DMOD(PLO+270.0D0,360.0D0)
C
 2000 CONTINUE      

      pla_r4 = PLA
      plo_r4 = PLO
!     glat_r4 = GLAT
!     glon_r4 = GLON
      rlat_r4 = RLAT
      wlon1_r4 = WLON1

      RETURN
      END                                  
C
C     ******************************************************************
C
      SUBROUTINE PSTOXY(X,Y,PLA,PLO,ERAD)
C
C     This program convert polar stereographic coordinates to x,y ditto
C     longitude:   0 - 360  ; positive to the east
C     latitude : -90 -  90  ; positive for northern hemisphere
C     it is assumed that the x-axis point towards the east and
C     corresponds to longitude = 0
C
C     TSP 20/06-89
C
C     constants and functions
C            
      FAC(PLA) = ERAD*2.0/(1.0+SIN(PLA*PI180))*COS(PLA*PI180)
      XC(PLA,PLO) = FAC(PLA)*COS(PLO*PI180)
      YC(PLA,PLO) = FAC(PLA)*SIN(PLO*PI180)      
      PI180=3.14159265/180.0
C
C     Calculate the coordinates
C
      X = XC(PLA,PLO)
      Y = YC(PLA,PLO)
C
      RETURN
      END

      subroutine latlon_to_xy_old(glat,glon,rlat,wlon1,erad,x,y)

!     combines getops + pstoxy

!                  o   o   i    i    i     i
      call GETOPS(PLA,PLO,GLAT,GLON,RLAT,WLON1)

      plo = plo - wlon1

!                 o o  i   i    i
      call PSTOXY(X,Y,PLA,PLO,ERAD)

      return
      end

      subroutine latlon_to_xy(glat,glon,erad,x,y)

!     combines getops + pstoxy

      call latlon_to_uv(glat,glon,u,v,istatus)

      call uv_to_xy(u,v,erad,x,y)

      return
      end

      subroutine xy_to_latlon_old(x,y,erad,rlat,wlon1,glat,glon)

!     combines xytops + pstoge

!                 i i  o   o    i
      call XYTOPS(X,Y,PLA,PLO,ERAD)

!                  i   i    o    o    i    i
      call PSTOGE(PLA,PLO,GLAT,GLON,RLAT,WLON1)

      return
      end


      subroutine xy_to_latlon(x,y,erad,glat,glon)

!     combines xytops + pstoge

      call xy_to_uv(x,y,erad,u,v)

      call uv_to_latlon(u,v,glat,glon,istatus)

      return
      end


      subroutine xy_to_uv(x,y,erad,u,v)

      use mem_namelist, ONLY: c6_maproj,slat1=>standard_latitude
     1                                 ,slat2=>standard_latitude2
     1                                 ,grid_spacing_m

      include 'trigd.inc'

      real n 

      if(c6_maproj .eq. 'plrstr')then     ! Haltiner & Williams 1-21
          call get_ps_parms(slat1,slat2,grid_spacing_m,phi0
     1                                 ,grid_spacing_proj_m)
          factor = 1. + sind(phi0)

          u = x / (factor * erad)
          v = y / (factor * erad)

      elseif(c6_maproj .eq. 'lambrt')then 
          call lambert_parms(slat1,slat2,n,s,rconst)

          u = x / (erad * rconst)
          v = y / (erad * rconst)

      elseif(c6_maproj .eq. 'merctr')then ! Haltiner & Williams 1-8-2
          u = x / erad
          v = y / erad

      elseif(c6_maproj .eq. 'latlon')then 
          u = x / erad
          v = y / erad

      else
          write(6,*)'xy_to_uv - Error: invalid map projection '
     1             ,c6_maproj             
          stop

      endif

      return
      end


      subroutine uv_to_xy(u,v,erad,x,y)

      use mem_namelist, ONLY: c6_maproj,slat1=>standard_latitude
     1                                 ,slat2=>standard_latitude2
     1                                 ,grid_spacing_m
      include 'trigd.inc'

      real n 

      if(c6_maproj .eq. 'plrstr')then     ! Haltiner & Williams 1-21
          call get_ps_parms(slat1,slat2,grid_spacing_m,phi0
     1                                 ,grid_spacing_proj_m)
          factor = 1. + sind(phi0)

          x = u * factor * erad
          y = v * factor * erad

      elseif(c6_maproj .eq. 'lambrt')then 
          call lambert_parms(slat1,slat2,n,s,rconst)

          x = u * (erad * rconst)
          y = v * (erad * rconst)

      elseif(c6_maproj .eq. 'merctr')then ! Haltiner & Williams 1-8-2
          x = u * erad
          y = v * erad

      elseif(c6_maproj .eq. 'latlon')then 
          x = u * erad
          y = v * erad

      else
          write(6,*)'uv_to_xy - Error: invalid map projection '
     1             ,c6_maproj       
          stop

      endif

      return
      end
