      subroutine polar_stereographic
c
c *** Routines to convert from geographical lat, lon to 
c        polar stereographic grid i, j and vice-versa.
c     Equations and code mostly obtained from RAMS.
c     snook (12/20/95)
c
ccc      implicit none
c
      integer np,n
c
      real*4 glat(np),glon(np),    !Earth lat (deg N), lon (deg +E)
     .       pslat,pslon,          !Pol ste. lat, lon (deg N, deg +E)
     .       psi(np),psj(np),      !Pol ste. i, j
     .       xmin,ymin,dx,dy
c
      common /pscorner/xmin,ymin,dx,dy
c
c===============================================================================
c
      entry latlon_2_psij(np,glat,glon,psi,psj)
c_______________________________________________________________________________
c
      call ps_param(xmin,ymin,dx,dy)
      do n=1,np
         call geoll_2_psll(glat(n),glon(n),pslat,pslon)
         call psll_2_psij(pslat,pslon,psi(n),psj(n))
      enddo
c
      return
c
c===============================================================================
c
      entry psij_2_latlon(np,psi,psj,glat,glon)
c_______________________________________________________________________________
c
c
      do n=1,np
         call psij_2_psll(psi(n),psj(n),pslat,pslon)
         call psll_2_geoll(pslat,pslon,glat(n),glon(n))
      enddo
c
      return
c
      end
c
c===============================================================================
c
      subroutine geoll_2_psll(glat,glon,PLA,PLO)
C
C     Convert geographical lat/lon coordinates to polar stereographic
C     ditto with the pol.ste. pole at RLAT,WLON1 (these names are
C     used to be compatible to the RAMSIN-parameters)
C     longitude:-180 ; 180 positive east (on input)
C              :   0 ; 360 positive east (on output)
C     latitude : -90 ;  90 posive on northern hemisphere
C     The result is rotated 270 degrees relative to 'standard pol.ste.'
C     WLON1 is defined in the same way as the input
C     approach so as to get the x-axis to point towards the east, and the
C     y-axis towards the north along 0 degrees (at NP south along 180)
C
C     TSP 20/06-89
      double precision pi180,c1,c2,c3,c4,c5,c6,arg2a,bb,pla1,alpha
     +   ,plo1,pla90,argu2
c
      integer nx,ny,nz           !No. of PS domain grid points
      real*4 RLAT,WLON1,rota,       !Pol ste. std lat, lon and rotation
     .       sw(2),ne(2)           !SW lat, lon, NE lat, lon
      common /psgrid/nx,ny,nz,RLAT,WLON1,rota,sw,ne
c_______________________________________________________________________________
C
C     constants
C
      c1=1.
      PI180 = dasin(c1)/90.
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
         PLO=AMOD(GLOR+270.0-WLON1,360.0)
         GO TO 2000
      END IF
      IF(RLAT.EQ.-90.0) THEN
         PLA=-GLAT
         PLO=AMOD(GLOR+270.0,360.0)
         GO TO 2000
      END IF
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
      IF(AMOD(GLOR+180.0,360.0).EQ.RWLON1) THEN
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
      c2=glat*pi180
      c3=argu1*pi180
      ARG2A = dCOS(c2)*dCOS(c3)
      ARG2A = dMAX1(ARG2A,-c1)
      ARG2A = dMIN1(ARG2A, c1)         
      BB    = dACOS(ARG2A)
C
      c4=hsign*glat*pi180
      ARG2A = dSIN(c4)/dSIN(BB)
      ARG2A = dMAX1(ARG2A,-c1)
      ARG2A = dMIN1(ARG2A, c1)
      ALPHA = dASIN(ARG2A)
C
C     2. Get PLA and PLO (still legalizing arguments)
C
      c5=rlat*pi180
      c6=hsign*rlat*pi180
      ARG2A = dCOS(c5)*dCOS(BB)+
     +        dSIN(c6)*dSIN(c4)
      ARG2A = dMAX1(ARG2A,-c1)
      ARG2A = dMIN1(ARG2A, c1)         
      PLA1   = dASIN(ARG2A)
C
      ARG2A = dSIN(BB)*dCOS(ALPHA)/dCOS(PLA1)
      ARG2A = dMAX1(ARG2A,-c1)
      ARG2A = dMIN1(ARG2A, c1)
      PLO1   = dASIN(ARG2A)
C
C    Test for passage of the 90 degree longitude (duallity in PLO)
C         Get PLA for which PLO=90 when GLAT is the latitude
C
      ARG2A = dSIN(c4)/dSIN(c6)
      ARG2A = dMAX1(ARG2A,-c1)
      ARG2A = dMIN1(ARG2A, c1)         
      PLA90 = dASIN(ARG2A)
C
C         Get help arc BB and angle ALPHA
C
      ARG2A = dCOS(c5)*dSIN(PLA90)
      ARG2A = dMAX1(ARG2A,-c1)
      ARG2A = dMIN1(ARG2A, c1)
      BB    = dACOS(ARG2A)

      ARG2A = dSIN(c4)/dSIN(BB)
      ARG2A = dMAX1(ARG2A,-c1)
      ARG2A = dMIN1(ARG2A, c1)        
      ALPHA = dASIN(ARG2A)
C
C         Get GLOLIM - it is nesc. to test for the existence of solution
C
      ARGU2  = dCOS(c2)*dCOS(BB)/
     +            (1.-dSIN(c4)*dSIN(BB)*dSIN(ALPHA))
      IF(dABS(ARGU2).GT.c1) THEN
      GLOLIM = 999.0
      ELSE
        GLOLIM = dACOS(ARGU2)/PI180
      END IF
C
C     Modify (if nesc.) the PLO solution
C
      IF((ABS(ARGU1).GT.GLOLIM.AND.GLAT.LE.RLAT).OR.
     +   GLAT.GT.RLAT) THEN
            PLO1 = PI180*180.0 - PLO1
      END IF
C
C     The solution is symmetric so the direction must be if'ed
C
      IF(ARGU1.LT.0.0) THEN
         PLO1 = -PLO1
      END IF
C
C     Convert the radians to degrees
C
      PLA = PLA1/PI180        
      PLO = PLO1/PI180
C
C     To obtain a rotated value (ie so x-axis in pol.ste. points east)
C     add 270 to longitude
C
      PLO=AMOD(PLO+270.0,360.0)
C
 2000 CONTINUE      
      RETURN
      END                                  
c
c===============================================================================
c
      subroutine psll_2_psij(pslat,pslon,psi,psj)
c
ccc      implicit none
c
      real*4 pslat,pslon,      !Pol ste. lat, lon (deg N, deg +E)
     .       psi,psj,          !Pol ste. i,j
     .       x,y,xmin,ymin,dx,dy,
     .       mag
c
      common /pscorner/xmin,ymin,dx,dy
c_______________________________________________________________________________
c
      mag=2./(1.+sind(pslat))*cosd(pslat)
      x=mag*cosd(pslon)
      y=mag*sind(pslon)
      psi=(x-xmin)/dx+1.
      psj=(y-ymin)/dy+1.
C
      return
      end
c
c===============================================================================
c
      subroutine psij_2_psll(psi,psj,pslat,pslon)
c
ccc      implicit none
c
      real*4 psi,psj,          !Pol ste. i,j
     .       pslat,pslon,      !Pol ste. lat, lon (deg N, deg +E)
     .       x,y,dist,
     .       xmin,ymin,dx,dy 
c
      common /pscorner/xmin,ymin,dx,dy
c_______________________________________________________________________________
c
      x=(psi-1.)*dx+xmin
      y=(psj-1.)*dy+ymin
      dist=sqrt(x**2+y**2)
      if (dist .eq. 0) then
         pslat=90.
         pslon=-90.
      else
         pslat=atand(dist/2.)
         pslat=90.-2.*pslat
c
         if (x .eq. 0.) then
            pslon=90.
         else
            if (x .gt. 0.) then
               pslon=atand(y/x)
            else
               pslon=atand(y/x)+180.
            endif
         endif
      endif
c
      pslon=amod(pslon+450.,360.)
c
      return
      end
c
c===============================================================================
c
      subroutine psll_2_geoll(pla,plo,glat,glon)
C
C     Convert polar stereographic coordinates to geographical lat/lon
C     ditto with the pol.ste. pole at rlat,wlon1 (these names are
C     used to be compatible to the ramsin-parameters)
C     longitude:   0 ; 360 positive east (on input)
C               -180 ; 180 positive east (on output)
C     latitude : -90 ;  90 posive on northern hemisphere
C     It is assumed that the polar stereographic coordinates have been
C     rotated to the standard format with 0 degrees longitude along wlon1
C
C     TSP 21 JUNE 89
c
      integer nx,ny,nz           !No. of PS domain grid points
      real*4 RLAT,WLON1,rota,      !Pol ste. std lat, lon and rotation
     .       sw(2),ne(2)           !SW lat, lon, NE lat, lon
      common /psgrid/nx,ny,nz,RLAT,WLON1,rota,sw,ne
c_______________________________________________________________________________
C
C     set flag for n/s hemisphere
C
      c1 = 1.
      PI180 = asin(c1)/90.
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
      ARG2A = COS(GLAT)*COS(BB)/(1.0-SIN(GLAT)**2)
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
c
c===============================================================================
c
      subroutine ps_param(xmin,ymin,dx,dy)
c
ccc      implicit none
c
      real*4 pslat,pslon,
     .       xmin,xmax,ymin,ymax,
     .       dx,dy,mag
c
      integer nx,ny,nz           !No. of PS domain grid points
      real*4 lat0,lon0,rota,       !Pol ste. std lat, lon and rotation
     .       sw(2),ne(2)           !SW lat, lon, NE lat, lon
      common /psgrid/nx,ny,nz,lat0,lon0,rota,sw,ne
c
c_______________________________________________________________________________
c
      call geoll_2_psll(sw(1),sw(2),pslat,pslon)
      mag=2./(1.+sind(pslat))*cosd(pslat)
      xmin=mag*cosd(pslon)
      ymin=mag*sind(pslon)
      call geoll_2_psll(ne(1),ne(2),pslat,pslon)
      mag=2./(1.+sind(pslat))*cosd(pslat)
      xmax=mag*cosd(pslon)
      ymax=mag*sind(pslon)
      dx=(xmax-xmin)/float(nx-1)
      dy=(ymax-ymin)/float(ny-1)
c
      return
      end
c
c===============================================================================
c===============================================================================
c
      subroutine conical_equidistant
c
c *** Routines to convert from geographical lat, lon to 
c        conical-equidistant grid i, j vice-versa.
c     Equations obtained from Adrian Marroquin and Tom Black (NCEP).
c     snook (12/20/95)
c
ccc      implicit none
c
      integer np,n
c
      real*4 glat(np),glon(np),      !Earth lat, lon (deg N, deg +E)
     .       celat,celon,            !Con eq. lat, lon (deg N, deg +E)
     .       cei(np),cej(np)         !Con eq. i,j
c
c===============================================================================
c
      entry latlon_2_coneqij(np,glat,glon,cei,cej)
c_______________________________________________________________________________
c
      do n=1,np
         call geoll_2_coneqll(glat(n),glon(n),celat,celon)
         call coneqll_2_coneqij(celat,celon,cei(n),cej(n))
      enddo
      return
c
c===============================================================================
c
      entry coneqij_2_latlon(np,cei,cej,glat,glon)
c_______________________________________________________________________________
c
      do n=1,np
         call coneqij_2_coneqll(cei(n),cej(n),celat,celon)
         call coneqll_2_geoll(celat,celon,glat(n),glon(n))
      enddo
      return
c
      end
c
c===============================================================================
c
      subroutine geoll_2_coneqll(glat,glon,celat,celon)
c
ccc      implicit none
c
      real*4 glat,glon,      !Earth lat, lon (deg N, deg +E)
     .       celat,celon,    !Con eq. lat, lon (deg N, deg +E)
     .       x,y,z
c
      integer nx,ny,nz
      real*4 lat0,lon0,dphi,dlam
      common /coneqgrid/nx,ny,nz,lat0,lon0,dphi,dlam
c_______________________________________________________________________________
c
      x=cosd(lat0)*cosd(glat)*cosd(lon0-glon)+sind(lat0)*sind(glat)
      y=-cosd(glat)*sind(lon0-glon)
      z=-sind(lat0)*cosd(glat)*cosd(lon0-glon)+cosd(lat0)*sind(glat)
c
      celat=atand(z/(x**2+y**2)**0.5)
      celon=atand(y/x)
c
      return
c
      end
c
c===============================================================================
c
      subroutine coneqll_2_coneqij(celat,celon,cei,cej)
c
ccc      implicit none
c
      real*4 celat,celon,      !Con eq. lat, lon (deg N, deg +E)
     .       cei,cej           !Con eq. i,j
c
      integer nxt
c
      integer nx,ny,nz
      real*4 lat0,lon0,dphi,dlam
      common /coneqgrid/nx,ny,nz,lat0,lon0,dphi,dlam
c_______________________________________________________________________________
c
      nxt=2*nx-1
      cei=float((nxt-1)/2+1)+celon/dlam
      cej=float((ny -1)/2+1)+celat/dphi
c
      return
      end
c
c===============================================================================
c
      subroutine coneqij_2_coneqll(cei,cej,celat,celon)
c
ccc      implicit none
c
      real*4 cei,cej,          !Con eq. i,j
     .       celat,celon       !Con eq. lat, lon (deg N, deg +E)
c
      integer nyt
c
      integer nx,ny,nz
      real*4 lat0,lon0,dphi,dlam
      common /coneqgrid/nx,ny,nz,lat0,lon0,dphi,dlam
c_______________________________________________________________________________
c
      nyt=ny/2+1
      celon=(cei-float(nx ))*dlam
      celat=(cej-float(nyt))*dphi
c
      return
      end
c
c===============================================================================
c
      subroutine coneqll_2_geoll(celat,celon,glat,glon)
c
ccc      implicit none
c
      real*4 celat,celon,    !Con eq. lat, lon (deg N, deg +E)
     .       glat,glon       !Earth lat, lon (deg N, deg +E)
c
      integer nx,ny,nz
      real*4 lat0,lon0,dphi,dlam
      common /coneqgrid/nx,ny,nz,lat0,lon0,dphi,dlam
c_______________________________________________________________________________
c
      glat=asind(sind(celat)*cosd(lat0)+
     .           cosd(celat)*sind(lat0)*cosd(celon))
      if (celon .lt. 0.) then
         glon=lon0-acosd(cosd(celat)*cosd(celon)/cosd(glat)/cosd(lat0)-
     .                   tand(glat)*tand(lat0))
      else
         glon=lon0+acosd(cosd(celat)*cosd(celon)/cosd(glat)/cosd(lat0)-
     .                   tand(glat)*tand(lat0))
      endif
c
      return
      end
c
c===============================================================================
c===============================================================================
c
      subroutine lambert_conformal
c
c *** Routines to convert from geographical lat, lon to 
c        Lambert-conformal grid i, j and vice-versa.
c     Equations obtained from NCAR graphics documentation.
c     snook (12/20/95)
c
ccc      implicit none
c
      integer np,n
c
      real*4 glat(np),glon(np),    !Earth lat (deg N), lon (deg +E)
     .       lci(np),lcj(np),      !Lambert-confomal i, j
     .       s,cone,r,
     .       xmin,ymin,dx,dy,x,y
c
      real*4 lat1,lat2,lon0,       !Lambert-conformal std lat1, lat2, lon
     .       sw(2),ne(2)           !SW lat, lon, NE lat, lon
      integer nx,ny,nz           !No. of LC domain grid points
      common /lcgrid/nx,ny,nz,lat1,lat2,lon0,sw,ne
c
c===============================================================================
c
      entry latlon_2_lcij(np,glat,glon,lci,lcj)
c_______________________________________________________________________________
c
      call lc_param(s,cone,xmin,ymin,dx,dy)
c
      do n=1,np
         r=(tand(45.-s*glat(n)/2.))**cone
         x=r*sind(cone*(glon(n)-lon0))
         y=-s*r*cosd(cone*(glon(n)-lon0))
         lci(n)=(x-xmin)/dx+1
         lcj(n)=(y-ymin)/dy+1
      enddo
c
      return
c
c===============================================================================
c
      entry lcij_2_latlon(np,lci,lcj,glat,glon)
c_______________________________________________________________________________
c
      call lc_param(s,cone,xmin,ymin,dx,dy)
c
      do n=1,np
         x=(lci(n)-1)*dx+xmin
         y=(lcj(n)-1)*dy+ymin
         glon(n)=lon0+atand(-s*x/y)/cone
         glat(n)=(90.-
     .            2.*atand((x/sind(cone*(glon(n)-lon0)))**(1./cone)))/s
      enddo
c
      return
c
      end
c
c===============================================================================
c
      subroutine lc_param(s,cone,xmin,ymin,dx,dy)
c
ccc      implicit none
c
      real*4 s,cone,r,
     .       xmin,xmax,ymin,ymax,
     .       dx,dy
c
      real*4 lat1,lat2,lon0,       !Lambert-conformal std lat1, lat2, lon
     .       sw(2),ne(2)           !SW lat, lon, NE lat, lon
      integer nx,ny,nz           !No. of LC domain grid points
      common /lcgrid/nx,ny,nz,lat1,lat2,lon0,sw,ne
c_______________________________________________________________________________
c
      if (lat1 .ge. 0.) then
         s=1.
      else
         s=-1.
      endif
      if (lat1 .ne. lat2) then
         cone=alog(cosd(lat1)/cosd(lat2))/
     .        alog(tand(45.-s*lat1/2.)/tand(45.-s*lat2/2.))
      else
         cone=cosd(90.-s*lat1)
      endif
c
      r=(tand(45.-s*sw(1)/2.))**cone
      xmin=r*sind(cone*(sw(2)-lon0))
      ymin=-s*r*cosd(cone*(sw(2)-lon0))
      r=(tand(45.-s*ne(1)/2.))**cone
      xmax=r*sind(cone*(ne(2)-lon0))
      ymax=-s*r*cosd(cone*(ne(2)-lon0))
      dx=(xmax-xmin)/float(nx-1)
      dy=(ymax-ymin)/float(ny-1)
c
      return
      end
c
c===============================================================================
c===============================================================================
c
      subroutine lat_lon
c
c *** Routines to convert from geographical lat, lon to 
c        lat-lon grid i, j and vice-versa.
c     snook (11/5/96)
c
ccc      implicit none
c
      integer np,n
c
      real*4 glat(np),glon(np),      !Earth lat, lon (deg N, deg +E)
     .       lli(np),llj(np),        !Lat-lon grid i,j
     .       diff
c
      integer nx,ny,nz             !No. of LL domain grid points
      real*4 lat0,lon0,dlat,dlon     !SW corner lat, lon, lat, lon spacing
      common /llgrid/nx,ny,nz,lat0,lon0,dlat,dlon
c
c===============================================================================
c
      entry latlon_2_llij(np,glat,glon,lli,llj)
c_______________________________________________________________________________
c
      do n=1,np
         diff=glon(n)-lon0
         if (diff .lt. 0.) diff=diff+360.
         if (diff .ge. 360.) diff=diff-360.
         lli(n)=diff/dlon+1.
         llj(n)=(glat(n)-lat0)/dlat+1.
      enddo
      return
c
c===============================================================================
c
      entry llij_2_latlon(np,lli,llj,glat,glon)
c_______________________________________________________________________________
c
      do n=1,np
         glon(n)=(lli(n)-1.)*dlon+lon0
         glat(n)=(llj(n)-1.)*dlat+lat0
      enddo
      return
c
      end
