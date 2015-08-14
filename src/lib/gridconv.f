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
      real glat(np),glon(np),    !Earth lat (deg N), lon (deg +E)
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
      include 'trigd.inc'
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
      real RLAT,WLON1,rota,       !Pol ste. std lat, lon and rotation
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
      include 'trigd.inc'
      real pslat,pslon,      !Pol ste. lat, lon (deg N, deg +E)
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
      include 'trigd.inc'
      real psi,psj,          !Pol ste. i,j
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
      include 'trigd.inc'
C     rotated to the standard format with 0 degrees longitude along wlon1
C
C     TSP 21 JUNE 89
c
      integer nx,ny,nz           !No. of PS domain grid points
      real RLAT,WLON1,rota,      !Pol ste. std lat, lon and rotation
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
      include 'trigd.inc'
      real pslat,pslon,
     .       xmin,xmax,ymin,ymax,
     .       dx,dy,mag
c
      integer nx,ny,nz           !No. of PS domain grid points
      real lat0,lon0,rota,       !Pol ste. std lat, lon and rotation
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
      real glat(np),glon(np),      !Earth lat, lon (deg N, deg +E)
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
      include 'trigd.inc'
      real glat,glon,      !Earth lat, lon (deg N, deg +E)
     .       celat,celon,    !Con eq. lat, lon (deg N, deg +E)
     .       x,y,z
c
      integer nx,ny,nz
      real lat0,lon0,dphi,dlam
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
      real celat,celon,      !Con eq. lat, lon (deg N, deg +E)
     .       cei,cej           !Con eq. i,j
c
      integer nxt
c
      integer nx,ny,nz
      real lat0,lon0,dphi,dlam
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
      real cei,cej,          !Con eq. i,j
     .       celat,celon       !Con eq. lat, lon (deg N, deg +E)
c
      integer nyt
c
      integer nx,ny,nz
      real lat0,lon0,dphi,dlam
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
      include 'trigd.inc'
      real celat,celon,    !Con eq. lat, lon (deg N, deg +E)
     .       glat,glon       !Earth lat, lon (deg N, deg +E)
c
      integer nx,ny,nz
      real lat0,lon0,dphi,dlam
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
      include 'trigd.inc'
      integer np,n
c
      real glat(np),glon(np),    !Earth lat (deg N), lon (deg +E)
     .       lci(np),lcj(np),      !Lambert-confomal i, j
     .       s,cone,r,
     .       xmin,ymin,dx,dy,x,y
c
      real lat1,lat2,lon0,       !Lambert-conformal std lat1, lat2, lon
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
      include 'trigd.inc'
      real s,cone,r,
     .       xmin,xmax,ymin,ymax,
     .       dx,dy
c
      real lat1,lat2,lon0,       !Lambert-conformal std lat1, lat2, lon
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
c     smart/pw li (04/1/03)
c
      integer np,n
      integer nx,ny,nz             !No. of LL domain grid points
 
      real   glat(np),glon(np)     !Earth lat, lon (deg N, deg +E)
      real   lli(np),llj(np)       !Lat-lon grid i,j
      real   diff
c     real   sw(2),ne(2)
c     common /llgrid/nx,ny,sw,ne,cgrddef

      real   lat0,lon0,dlat,dlon     !SW corner lat, lon, lat, lon spacing
      character*1 cgrddef
      common /llgrid/nx,ny,nz,lat0,lon0,dlat,dlon,cgrddef
c
c===============================================================================
c
      entry latlon_2_llij(np,glat,glon,lli,llj)
c_______________________________________________________________________________
c

c pw li code for integration later on.
c     dlat=(ne(2)-sw(2))/(nx-1)
c     dlon=(ne(1)-sw(1))/(ny-1)
c     do i=1,n
c        lli(i)=(glon(i)-sw(2))/dlon + 1.
c        llj(i)=(glat(i)-sw(1))/dlat + 1.
c     enddo
 
      dlond=dlon
      dlatd=dlat
      if(cgrddef.eq.'S')then
         do n=1,np
            diff=glon(n)-lon0
            if (diff .lt. 0.) diff=diff+360.
            if (diff .ge. 360.) diff=diff-360.
            lli(n)=diff/dlond+1.
            llj(n)=(glat(n)-lat0)/dlatd+1.
         enddo
      elseif(cgrddef.eq.'N')then
         do n=1,np
            diff=glon(n)-lon0
            if (diff .lt. 0.) diff=diff+360.
            if (diff .ge. 360.) diff=diff-360.
            lli(n)=diff/dlon+1.
            llj(n)=(lat0-glat(n))/dlat+1.
         enddo

      else
         print*,'you must specify whether the standard
     .           lat is Southern or Northern boundary'
      endif

      return
c
c===============================================================================
c
      entry llij_2_latlon(np,lli,llj,glat,glon)
c_______________________________________________________________________________
c
c Note: if lat0=northern boundary of ll grid set dlat=-dlat
      do n=1,np
         glon(n)=(lli(n)-1.)*dlon+lon0
         glat(n)=(llj(n)-1.)*dlat+lat0
      enddo
      return
c
      end
c
c===============================================================================
c===============================================================================
c

      subroutine cylindrical_equidistant
c
c routine to convert from cylindrical equidistant to grid ri/rj
c used for mapping the WSI radar data to laps.
c
c     J. Smart 11-19-98   Original Working Version
c                         Basic equation set from
c                         Map Projections Used by the U.S. Geological Survey
c                         (Snyder, J.P. 1983 Geological Survey Bulletin - 1532)
c
      include 'trigd.inc'
      implicit none
      integer np,n
 
      real glat(np),glon(np),      !Earth lat, lon (deg N, deg +E)
     .       lli(np),llj(np),        !Lat-lon grid i,j
     .       nw(2),se(2)             !NW grid lat/lon; SE grid lat/lon

      double precision     diff,x,y
     .                    ,xmin,ymin
     .                    ,xmax,ymax
     .                    ,pi,dg2rd

      double precision     coslatc
     .                    ,dx,dy
     .                    ,drlatc,drlonc
     .                    ,dglat, dglon

!     double precision     dcosd

      real  r
 
      integer nx,ny,nz             !No. of LL domain grid points
      real  rlatc,rlonc          !Grid center lat, lon

c     real coslatc

      data r/6.3712e6/             !this earth radius is = to that in lapsgrid.f
      common /cegrid/nx,ny,nz,nw,se,rlatc,rlonc
c
c===============================================================================
c
      entry latlon_2_ceij(np,glat,glon,lli,llj)
c_______________________________________________________________________________
c
c Note: WSI grid (1,1) is NW corner. Adjustment made for rj using ny since
c       equation set assumes y-axis is on equator.
c
      pi=acos(-1.0)
      dg2rd=pi/180.
      drlatc=rlatc
      drlonc=rlonc

      xmin=r*((nw(2)-drlonc)*dg2rd)*dcosd(drlatc)
      ymin=r*se(1)*dg2rd
      xmax=r*((se(2)-drlonc)*dg2rd)*dcosd(drlatc)
      ymax=r*nw(1)*dg2rd

      diff=xmax-xmin
      if(diff.lt.0.000001)then
         xmax=dabs(xmin)
      endif

      dx=(xmax-xmin)/(nx-1)
      dy=(ymax-ymin)/(ny-1)

      coslatc=dcosd(drlatc)
      do n=1,np
         diff=glon(n)-drlonc
         if (diff .lt.-180.) diff=diff+360.
         if (diff .ge. 180.) diff=diff-360.
         diff=diff*dg2rd
         x=r*diff*coslatc
         y=r*(glat(n)*dg2rd)
         lli(n)=(x-xmin)/dx + 1.
         llj(n)=float(ny)-((y-ymin)/dy) + 1.
      enddo

      return
c
c===============================================================================
c
      entry ceij_2_latlon(np,lli,llj,glat,glon)
c
c equation set: rlat (phi) = y/r
c               rlon (lambda) = rlonc + x/(R cos(rlatc))
c        where:
c               r= earth radius
c               rlonc= central meridian
c               rlatc= central parallel

      pi=acos(-1.0)
      dg2rd=pi/180.
      drlatc=rlatc
      drlonc=rlonc

      xmin=r*((nw(2)-drlonc)*dg2rd)*dcosd(drlatc)
      ymin=r*se(1)*dg2rd
      xmax=r*((se(2)-drlonc)*dg2rd)*dcosd(drlatc)
      ymax=r*nw(1)*dg2rd

      diff=xmax-xmin
      if(diff.lt.0.000001)then
         xmax=dabs(xmin)
      endif

      dx=(xmax-xmin)/(nx-1)
      dy=(ymax-ymin)/(ny-1)

      coslatc=dcosd(drlatc)

      do n=1,np

         dglon=rlonc*dg2rd +
     &         (xmin+(lli(n)+1.0)*dx)/(r*coslatc)
         dglat=(ymin+(ny-llj(n)+1.0)*dy)/r

         glon(n)=dglon/dg2rd
         glat(n)=dglat/dg2rd

      enddo

      return
      end
c
c ===========================================================================
c
      subroutine latlon_2_mcij(n,rlat,rlon,ri,rj)

cc     implicit none
      include 'trigd.inc'

      integer n,i
      integer nx,ny
      real rlat(n),rlon(n),ri(n),rj(n)
      real rlonc,dlon,rlatc,dlat
      real sw(2),ne(2)
      real x,y
      real xmax,ymax
      real xmin,ymin
      real dx,dy
      real R,PI
      real deg2rad
      common /mcgrid/rlonc,rlatc,nx,ny,sw,ne

      call get_earth_radius(R,istatus)

      PI=acos(-1.)
      deg2rad=PI/180.

      ymax=R*(log(tan(PI/4.+0.5*(ne(1)-rlatc)*deg2rad)))
      ymin=R*(log(tan(PI/4.+0.5*(sw(1)-rlatc)*deg2rad)))

      xmax=R*(deg2rad*(ne(2)-rlonc))
      xmin=R*(deg2rad*(sw(2)-rlonc))

      dx=(xmax-xmin)/float(nx-1)
      dy=(ymax-ymin)/float(ny-1)

      do i=1,n
         dlon=rlon(i)-rlonc
         dlat=rlat(i)-rlatc
         if(dlon.gt.180.)dlon=dlon-360.
         if(dlon.lt.-180.)dlon=dlon+360.
         x=R*(deg2rad*dlon)
         y=R*(log(tan(PI/4.+0.5*(dlat*deg2rad))))
         ri(i)=(x-xmin)/dx + 1.
         rj(i)=(y-ymin)/dy + 1.
      enddo

      return
      end
c
c===============================================================================
c
      subroutine latlon_2_npij(n,rlat,rlon,ri,rj)
c
c this is conversion of lat/lon to ri/rj in a
c lat-lon grid. Very similar to latlon_2_llij.
c
      integer n,i
      integer nx,ny
      real rlat(n),rlon(n),ri(n),rj(n)
      real sw(2),ne(2)
      real dx,dy
      common /npgrid/nx,ny,sw,ne

c     print *, ' Inside latlon_2_npij'
c     print *, nx,ny,dx,dy
c     print *, sw(1), sw(2),ne(1),ne(2)

      dx=(ne(2)-sw(2))/(nx-1)
      dy=(ne(1)-sw(1))/(ny-1)
 
      do i=1,n
        ri(i)=(rlon(i)-sw(2))/dx + 1.
        rj(i)=(rlat(i)-sw(1))/dy + 1.
      enddo
 
      return
      end

c -------------------------------------------------------

      subroutine init_gridconv_cmn(gproj,nxbg,nybg,nzbg
     &,dlat,dlon,cenlat,cenlon,Lat0,Lat1,Lon0
     &,sw1,sw2,ne1,ne2,cgrddef,istatus)
c
c JS 4-01
c
      implicit none

      character*(*)  gproj
      character*(*)  cgrddef
      integer        istatus
      integer        nxbg,nybg,nzbg
      real           dlat,dlon
      real           Lat0,Lat1
      real           Lon0,Lon1,Lon2
      real           sw1,ne1
      real           sw2,ne2
      real           cenlat,cenlon

c
c *** Common block variables for lat-lon grid.
c
      integer   nx_ll,ny_ll,nz_ll
      real    lat0_ll,lon0_ll,d_lat,d_lon
      character*1 cgrddef_ll
      common /llgrid/nx_ll,ny_ll,nz_ll,lat0_ll,lon0_ll
     &,d_lat,d_lon,cgrddef_ll
c
c *** Common block variables for lambert-conformal grid.
c
      integer   nx_lc,ny_lc,nz_lc
      real    lat1_lc,lat2_lc,lon0_lc,sw_lc(2),ne_lc(2)
      common /lcgrid/nx_lc,ny_lc,nz_lc,lat1_lc,lat2_lc
     &,lon0_lc,sw_lc,ne_lc
c
c *** Common block variables for cyclindrical equidistant grid.
c
      integer   nx,ny,nz
      real    rlatc,rlonc,nw(2),se(2),dx,dy
      common /cegrid/nx,ny,nz,nw,se,rlatc,rlonc
c
c *** Common block variables for polar stereographic grid.
c
      integer nx_ps,ny_ps,nz_ps    !No. of PS domain grid points
      real lat0_ps,lon0_ps,rota  !Pol ste. std lat, lon and rotation
     .      ,sw_ps(2),ne_ps(2)     !SW lat, lon, NE lat, lon
      common /psgrid/nx_ps,ny_ps,nz_ps,lat0_ps,lon0_ps
     .              ,rota,sw_ps,ne_ps

      integer nx_np,ny_np
      real    sw_np(2),ne_np(2)
      common /npgrid/nx_np,ny_np,sw_np,ne_np
c
c *** Common block variables for mercator grid.
c
      integer nx_mc,ny_mc          !No. of domain grid points
      real    sw_mc(2),ne_mc(2)    !SW lat, lon, NE lat, lon
      real    rlatc_mc,rlonc_mc    !Center lat/lon of domain
      common /mcgrid/rlonc_mc,rlatc_mc,nx_mc,ny_mc,sw_mc,ne_mc

      if(gproj.eq.'LC')then
         nx_lc=nxbg
         ny_lc=nybg
         nz_lc=nzbg
         lat1_lc=Lat0
         lat2_lc=Lat1
         lon0_lc=Lon0
         sw_lc(1)=sw1
         sw_lc(2)=sw2
         ne_lc(1)=ne1
         ne_lc(2)=ne2
      endif

      if(gproj.eq.'LL')then
         nx_ll=nxbg
         ny_ll=nybg
         nz_ll=nzbg
         lat0_ll=Lat0
         lon0_ll=Lon0
         d_lat=dlat
         d_lon=dlon
         cgrddef_ll=cgrddef
      endif

      if(gproj.eq.'LE')then
         nx=nxbg
         ny=nybg
         nz=nzbg
         rlatc=cenlat
         rlonc=cenlon
         se(1)=sw1
         se(2)=sw2
         nw(1)=ne1
         nw(2)=ne2
      endif

      if(gproj.eq.'PS')then
         nx_ps=nxbg
         ny_ps=nybg
         nz_ps=nzbg
         lat0_ps=Lat0
         lon0_ps=Lon0
         sw_ps(1)=sw1
         sw_ps(2)=sw2
         ne_ps(1)=ne1
         ne_ps(2)=ne2
      endif

      if(gproj.eq.'MC')then
         nx_mc=nxbg
         ny_mc=nybg
         rlatc_mc=cenlat
         rlonc_mc=cenlon
         sw_mc(1)=sw1
         sw_mc(2)=sw2
         ne_mc(1)=ne1
         ne_mc(2)=ne2
      endif

      return
      end
c
c===============================================================================
c
      subroutine init_hinterp(nx_bg,ny_bg,nx_laps,ny_laps,gproj,
     .     lat,lon,grx,gry,bgmodel,cmodel,wrapped)

c
      implicit none
c
      integer nx_bg,ny_bg,nx_laps,ny_laps,bgmodel
c
      real lat(nx_laps,ny_laps),lon(nx_laps,ny_laps),
     .       grx(nx_laps,ny_laps),gry(nx_laps,ny_laps) 

      logical wrapped
c
      integer i,j,k
      integer istatus
c
      character*132 cmodel
      character*2   gproj

      integer lenc

      integer nxc,nyc,nzc
      integer nx,ny
      logical lprintmessage
      real sw(2),ne(2),rota,lat0,lon0
      real nw(2),se(2),rlatc,rlonc
      real tolx,toly
      real grxdifsum1,grxdifsum2
      real grydifsum1,grydifsum2
      real r_missing_data
c     parameter (tol=0.10)
      common /psgrid/nxc,nyc,nzc,lat0,lon0,rota,sw,ne
c________________________________________________________________________________
c
      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne. 1)then
         print*,'Error getting r_missing_data - init_hinterp'
         return
      endif
c
c *** Determine location of LAPS grid point in background data i,j space.
c
      if (gproj .eq. 'PS') then

c      print*,nxc,nyc,nzc,lat0,lon0,rota,sw,ne

         call latlon_2_psij(nx_laps*ny_laps,lat,lon,grx,gry)
      elseif (gproj .eq. 'LC') then
         call latlon_2_lcij(nx_laps*ny_laps,lat,lon,grx,gry)
      elseif (gproj .eq. 'CE') then
         call latlon_2_coneqij(nx_laps*ny_laps,lat,lon,grx,gry)
      elseif (gproj .eq. 'LL') then
         call latlon_2_llij(nx_laps*ny_laps,lat,lon,grx,gry)
      elseif (gproj .eq. 'LE') then
         call latlon_2_ceij(nx_laps*ny_laps,lat,lon,grx,gry)
      elseif (gproj .eq. 'NP') then
         call latlon_2_npij(nx_laps*ny_laps,lat,lon,grx,gry)
      elseif (gproj .eq. 'MC') then
         call latlon_2_mcij(nx_laps*ny_laps,lat,lon,grx,gry)
      else
         print*,"Error: Unknown gproj spec in gridconv ",gproj
      endif
c
c *** Check that all LAPS grid points are within the background data coverage.
c

c set tolerance based upon the grid spacing as a function of grx/gry
      grxdifsum1=0.0
      grxdifsum2=0.0
      grydifsum1=0.0
      grydifsum2=0.0
      do j=1,ny_laps
         grxdifsum1=grxdifsum1+(grx(2,j)-grx(1,j))
         grxdifsum2=grxdifsum2+(grx(nx_laps-1,j)-grx(nx_laps,j))
      enddo
      do i=1,nx_laps
         grydifsum1=grydifsum1+(gry(i,2)-gry(i,1))
         grydifsum2=grydifsum2+(gry(i,ny_laps-1)-gry(i,ny_laps))
      enddo

      tolx=(abs(grxdifsum1)/ny_laps+abs(grxdifsum2)/ny_laps)*0.5
      toly=(abs(grydifsum1)/nx_laps+abs(grydifsum2)/nx_laps)*0.5

      print*,'horiz mapping tolerance x/y: ',tolx,toly
      lprintmessage=.true.
c
c *** First, check for wrapping if a global data set.
c
cWNI WNI COMMENT:  THIS SECTION NOT NEEDED BECAUSE
cWNI  THE LATLON_2_LLIJ ALREADY HANDLES THIS
      call s_len(cmodel,lenc)

      if ( bgmodel .eq. 6 .or. 
     .     bgmodel .eq. 8) then
CWNI     do j=1,ny_laps
CWNI        do i=1,nx_laps
CWNI           if (grx(i,j) .lt. 1.) grx(i,j)=grx(i,j)+float(nx_bg)
CWNI           if (grx(i,j) .gt. nx_bg) grx(i,j)=grx(i,j)-float(nx_bg)
CWNI           if (gry(i,j) .lt. 1.) then
CWNI              gry(i,j)=2.-gry(i,j)
CWNI              grx(i,j)=grx(i,j)-float(nx_bg/2)
CWNI              if (grx(i,j) .lt. 1.) grx(i,j)=grx(i,j)+float(nx_bg)
CWNI              if (grx(i,j).gt.nx_bg) grx(i,j)=grx(i,j)-float(nx_bg)
CWNI           endif
CWNI           if (gry(i,j) .gt. ny_bg) then
CWNI              gry(i,j)=float(2*ny_bg)-gry(i,j)
CWNI              grx(i,j)=grx(i,j)-float(nx_bg/2)
CWNI              if (grx(i,j) .lt. 1.) grx(i,j)=grx(i,j)+float(nx_bg)
CWNI              if (grx(i,j).gt.nx_bg) grx(i,j)=grx(i,j)-float(nx_bg)
CWNI           endif
CWNI           if (grx(i,j) .lt. 1) then
CWNI              grx(i,j)=grx(i,j)+1.
CWNI           endif
CWNI        enddo
CWNI     enddo
CWNI  elseif(bgmodel.eq.4.and.cmodel(1:lenc).eq.'AVN_SBN_CYLEQ')then
      elseif( (bgmodel.eq.4.and.cmodel(1:lenc).eq.'AVN_SBN_CYLEQ') .OR. !WNI
     .        (bgmodel .eq. 10 .and. cmodel .eq. 'GFS_ISO')) then      !WNI

CWNI     do j=1,ny_laps
CWNI        do i=1,nx_laps
CWNI           if (grx(i,j) .lt. 1.) grx(i,j)=grx(i,j)+float(nx_bg)
CWNI           if (grx(i,j) .gt. nx_bg) grx(i,j)=grx(i,j)-float(nx_bg)
CWNI           if (gry(i,j) .lt. 1.) then
CWNI              gry(i,j)=2.-gry(i,j)
CWNI              grx(i,j)=grx(i,j)-float(nx_bg/2)
CWNI              if (grx(i,j) .lt. 1.) grx(i,j)=grx(i,j)+float(nx_bg)
CWNI              if (grx(i,j).gt.nx_bg) grx(i,j)=grx(i,j)-float(nx_bg)
CWNI           endif
CWNI           if (gry(i,j) .gt. ny_bg) then
CWNI              gry(i,j)=float(2*ny_bg)-gry(i,j)
CWNI              grx(i,j)=grx(i,j)-float(nx_bg/2)
CWNI              if (grx(i,j) .lt. 1.) grx(i,j)=grx(i,j)+float(nx_bg)
CWNI              if (grx(i,j).gt.nx_bg) grx(i,j)=grx(i,j)-float(nx_bg)
CWNI           endif
CWNI        enddo
CWNI     enddo

      elseif( wrapped ) then ! SCA
              


c
c ****** If not a global data set, then check that LAPS domain is fully
c           within background domain.
c
      else
         do j=1,ny_laps
            do i=1,nx_laps
c
c LAPS must fit into model grid which must also fit into LAPS grid thus we
c introduce a small fudge factor on the grid boundaries.
c               

               if(grx(i,j).gt.1.-tolx) grx(i,j) = max(1.,grx(i,j))
               if(gry(i,j).gt.1.-toly) gry(i,j) = max(1.,gry(i,j))

               if(grx(i,j).lt.float(nx_bg)+tolx) 
     +              grx(i,j) = min(float(nx_bg)-tolx,grx(i,j))
               if(gry(i,j).lt.float(ny_bg)+toly) 
     +              gry(i,j) = min(float(ny_bg)-toly,gry(i,j))

               if (grx(i,j) .lt. 1 .or. grx(i,j) .gt. nx_bg .or.
     .             gry(i,j) .lt. 1 .or. gry(i,j) .gt. ny_bg) then

                  grx(i,j) = r_missing_data
                  gry(i,j) = r_missing_data

           if(lprintmessage)then
              print*,'Domain gridpt outside of bkgd data coverage.'
              print*,'   data i,j,lat,lon - ',i,j,lat(i,j),lon(i,j)
              print*,'   grx, gry:',grx(i,j),gry(i,j)
              lprintmessage=.false.
c                 stop 'init_hinterp'
           endif

               endif
            enddo
         enddo
      endif
c
      return
      end

      subroutine hinterp_field(nx_bg,ny_bg,nx_laps,ny_laps,nz,
     .     grx,gry,fvi,flaps,wrapped)
c
      implicit none
c
      integer nx_bg,ny_bg,nx_laps,ny_laps,nz        

      logical wrapped
c
c *** Input vertically interpolated field - 2D should be slightly faster than
c     'hinterp_field_3d'
c *** Output Laps field
c
      real fvi(nx_bg,ny_bg,nz),
     .       flaps(nx_laps,ny_laps,nz),
     .       grx(nx_laps,ny_laps),gry(nx_laps,ny_laps)
c
      integer i,j,k
      integer istatus
      real    r_missing_data
c

      write(6,*)' Subroutine hinterp_field'

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)then
         print*,'Error getting r_missing_data - hinterp_field'
         return
      endif
c
c *** Horizontally interpolate variable.
c
      do k=1,nz
         do j=1,ny_laps
         do i=1,nx_laps
            if(grx(i,j).lt.r_missing_data.and. 
     .         gry(i,j).lt.r_missing_data)then
               call gdtost_i(fvi(1,1,k),nx_bg,ny_bg,
     .              grx(i,j),gry(i,j),flaps(i,j,k),wrapped)
            else
               flaps(i,j,k)=r_missing_data
            endif
         enddo
         enddo
      enddo
c
      return
      end

      subroutine hinterp_field_3d(nx_bg,ny_bg,nx_laps,ny_laps,nz,
     .     grx,gry,fvi,flaps,wrapped)
c
!     implicit none
c
      include 'bgdata.inc'

      integer nx_bg,ny_bg,nx_laps,ny_laps,nz        

      dimension r(4,nz),scr(4,nz),staval(nz)

      logical wrapped
c
c *** Input vertically interpolated field - optimized for 3D.
c *** Output Laps field
c
      real fvi(nx_bg,ny_bg,nz),
     .       flaps(nx_laps,ny_laps,nz),
     .       grx(nx_laps,ny_laps),gry(nx_laps,ny_laps)
c
      integer i,j,k
      integer istatus
      real    r_missing_data
c
      write(6,*)' Subroutine hinterp_field_3d'

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)then
         print*,'Error getting r_missing_data - hinterp_field'
         return
      endif
c
c *** Horizontally interpolate variable.
c
      ix = nx_bg
      iy = ny_bg

      do jl=1,ny_laps
      do il=1,nx_laps
        if(grx(il,jl).lt.r_missing_data.and. 
     .     gry(il,jl).lt.r_missing_data)then
!         call gdtost_i(fvi(1,1,k),nx_bg,ny_bg,
!    .         grx(i,j),gry(i,j),flaps(i,j,k),wrapped)
!         subroutine gdtost_i(ab,ix,iy,stax,stay,staval,wrapped)

          stax = grx(il,jl)
          stay = gry(il,jl)

            
!         dimension ab(ix,iy),r(4),scr(4)
!         logical wrapped ! WNI added
!         include 'bgdata.inc'
c_______________________________________________________________________________
c
          iy1=int(stay)-1
          if(stay.lt.1.0)iy1=1.0
          iy2=iy1+3
          ix1=int(stax)-1
          if(stax.lt.1.0) THEN 
            IF (wrapped) THEN   ! WNI
              ix1 = ix1 + ix
            ELSE ! WNI
              ix1=1.0
            ENDIF ! WNI
          endif ! WNI
          ix2=ix1+3
          fiym2=float(iy1)-1
          fixm2=float(ix1)-1

            ii=0

            do iii=ix1,ix2
               i=iii
c
c ****** Account for wrapping around effect of global data at Greenwich.
c
               if (wrapped) then    ! WNI
                  if (i .lt. 1) i=i+ix  
                  if (i .gt. ix) i=i-ix
               endif

               ii=ii+1
               if (i .ge. 1 .and. i .le. ix) then 
                  jj=0
                  do jjj=iy1,iy2
                     j=jjj
c
c ************ Account for N-S wrapping effect of global data.
c
                     if (wrapped) THEN  !WNI
                        if (j .lt. 1) then
                           j=2-j
                           i=i-ix/2
                           if (i .lt. 1) i=i+ix
                           if (i .gt. ix) i=i-ix
                        endif
                        if (j .gt. iy) then
                           j=2*iy-j
                           i=i-ix/2
                           if (i .lt. 1) i=i+ix
                           if (i .gt. ix) i=i-ix
                        endif
                     endif

                     jj=jj+1

                     if (j .ge. 1 .and. j .le. iy) then
                        r(jj,:)=fvi(i,j,:) ! ab(i,j)
                     else
                        r(jj,:)=missingflag
                     endif
                  enddo
                  yy=stay-fiym2

                  if (yy .eq. 2.0) then
                    scr(ii,:)=r(2,:)
                  else
      
                    staval(:)=missingflag
                    xxx=yy
                    wt1=(xxx-3.)/(-1.)
                    wt2=1.0-wt1
                    yz22_temp=wt1*(xxx-4.)/(-2.)
                    yz23_temp=wt2*(xxx-4.)/(-1.)
                    yz24_temp=(xxx-2.)*(xxx-3.)/2.

                    yz11_temp=(xxx-2.)*(xxx-3.)/(2.)
                    yz12_temp=wt1*(xxx-1.)/(1.)
                    yz13_temp=wt2*(xxx-1.)/(2.)

                    do k=1,nz

!                    call binom(1.,2.,3.,4.,r(1),r(2),r(3),r(4),yy,scr(ii)) (now inlined)
!                    subroutine binom(x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy)
                     y1 = r(1,k)
                     y2 = r(2,k)
                     y3 = r(3,k)
                     y4 = r(4,k)
c
                     yyy=missingflag
                     if ( .not. (y2 .gt. 1.e19 .or. y3 .gt. 1.e19) )then
c
c
                       if (y4 .lt. 1.e19) then
                         yz22=yz22_temp
                         yz23=yz23_temp
                         yz24=yz24_temp
                       else
                         yz22=wt1
                         yz23=wt2
                         yz24=0.0
                       endif
c
                       if (y1 .lt. 1.e19) then
                         yz11=yz11_temp
                         yz12=yz12_temp
                         yz13=yz13_temp
                       else
                         yz11=0.0
                         yz12=wt1
                         yz13=wt2
                       endif
c
                       if (yz11 .eq. 0. .and. yz24 .eq. 0.) then
                         yyy=wt1*y2+wt2*y3
                       else
                         yyy=wt1*(yz11*y1+yz12*y2+yz13*y3)
     1                      +wt2*(yz22*y2+yz23*y3+yz24*y4)
                       endif
                     endif

                     scr(ii,k) = yyy
!                    return
!                    end
                    enddo ! k
                  endif
               else 
                  scr(ii,:)=missingflag
               endif ! i is valid
            enddo ! iii

            xx=stax-fixm2
            if (xx .eq. 2.0) then
              staval(:)=scr(2,:)
            else

!             call binom(1.,2.,3.,4.,scr(1),scr(2),scr(3),scr(4),xx,staval) (now inlined)
!             subroutine binom(x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy)

              xxx=xx
              wt1=(xxx-3.)/(-1.)
              wt2=1.0-wt1
              yz22_temp=wt1*(xxx-4.)/(-2.)
              yz23_temp=wt2*(xxx-4.)/(-1.)
              yz24_temp=(xxx-2.)*(xxx-3.)/2.

              yz11_temp=(xxx-2.)*(xxx-3.)/(2.)
              yz12_temp=wt1*(xxx-1.)/(1.)
              yz13_temp=wt2*(xxx-1.)/(2.)

              do k=1,nz

               y1 = scr(1,k)
               y2 = scr(2,k)
               y3 = scr(3,k)
               y4 = scr(4,k)
c
               yyy=missingflag
               if ( .not. (y2 .gt. 1.e19 .or. y3 .gt. 1.e19) )then
c
                 wt1=(xxx-3.)/(-1.)
                 wt2=1.0-wt1
c
                 if (y4 .lt. 1.e19) then
                   yz22=yz22_temp
                   yz23=yz23_temp
                   yz24=yz24_temp
                 else
                   yz22=wt1
                   yz23=wt2
                   yz24=0.0
                 endif
c
                 if (y1 .lt. 1.e19) then
                   yz11=yz11_temp
                   yz12=yz12_temp
                   yz13=yz13_temp
                 else
                   yz11=0.0
                   yz12=wt1
                   yz13=wt2
                 endif
c
                 if (yz11 .eq. 0. .and. yz24 .eq. 0.) then
                   yyy=wt1*y2+wt2*y3
                 else
                   yyy=wt1*(yz11*y1+yz12*y2+yz13*y3)
     1                +wt2*(yz22*y2+yz23*y3+yz24*y4)
                 endif
               endif

               staval(k) = yyy
!              return
!              end

             enddo ! k

            endif ! xx = 2
c
            if(staval(1).eq.missingflag)then
               print*,'hinterp_field_3d: staval = missingflag ',il,jl
            endif

            flaps(il,jl,:) = staval(:)
!           return
!           end

        else
          flaps(il,jl,:)=r_missing_data
        endif

      enddo ! il
      enddo ! jl
c
      return
      end


      subroutine hinterp_field_2d(nx_bg,ny_bg,nx_laps,ny_laps,nz,
     .     grx,gry,fvi,flaps,wrapped)
c
!     implicit none
c
      include 'bgdata.inc'

      integer nx_bg,ny_bg,nx_laps,ny_laps,nz        

      dimension r(4),scr(4)

      logical wrapped
c
c *** Input vertically interpolated field.
c *** Output Laps field
c
      real fvi(nx_bg,ny_bg,nz),
     .       flaps(nx_laps,ny_laps,nz),
     .       grx(nx_laps,ny_laps),gry(nx_laps,ny_laps)
c
      integer i,j,k
      integer istatus
      real    r_missing_data
c
      write(6,*)' Subroutine hinterp_field_2d'

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)then
         print*,'Error getting r_missing_data - hinterp_field'
         return
      endif
c
c *** Horizontally interpolate variable.
c
      ix = nx_bg
      iy = ny_bg

      do k=1,nz
        do jl=1,ny_laps
        do il=1,nx_laps
          if(grx(il,jl).lt.r_missing_data.and. 
     .       gry(il,jl).lt.r_missing_data)then
!           call gdtost_i(fvi(1,1,k),nx_bg,ny_bg,
!    .           grx(i,j),gry(i,j),flaps(i,j,k),wrapped)
!           subroutine gdtost_i(ab,ix,iy,stax,stay,staval,wrapped)

            stax = grx(il,jl)
            stay = gry(il,jl)
            
!           dimension ab(ix,iy),r(4),scr(4)
!           logical wrapped ! WNI added
!           include 'bgdata.inc'
c_______________________________________________________________________________
c
            iy1=int(stay)-1
            if(stay.lt.1.0)iy1=1.0
            iy2=iy1+3
            ix1=int(stax)-1
            if(stax.lt.1.0) THEN 
              IF (wrapped) THEN   ! WNI
                ix1 = ix1 + ix
              ELSE ! WNI
                ix1=1.0
              ENDIF ! WNI
            endif ! WNI
            ix2=ix1+3
            staval=missingflag
            fiym2=float(iy1)-1
            fixm2=float(ix1)-1
            ii=0
            do iii=ix1,ix2
               i=iii
c
c ****** Account for wrapping around effect of global data at Greenwich.
c
               if (wrapped) then    ! WNI
                  if (i .lt. 1) i=i+ix  
                  if (i .gt. ix) i=i-ix
               endif

               ii=ii+1
               if (i .ge. 1 .and. i .le. ix) then 
                  jj=0
                  do jjj=iy1,iy2
                     j=jjj
c
c ************ Account for N-S wrapping effect of global data.
c
                     if (wrapped) THEN  !WNI
                        if (j .lt. 1) then
                           j=2-j
                           i=i-ix/2
                           if (i .lt. 1) i=i+ix
                           if (i .gt. ix) i=i-ix
                        endif
                        if (j .gt. iy) then
                           j=2*iy-j
                           i=i-ix/2
                           if (i .lt. 1) i=i+ix
                           if (i .gt. ix) i=i-ix
                        endif
                     endif

                     jj=jj+1
                     if (j .ge. 1 .and. j .le. iy) then
                        r(jj)=fvi(i,j,k) ! ab(i,j)
                     else
                        r(jj)=missingflag
                     endif
                  enddo
                  yy=stay-fiym2
                  if (yy .eq. 2.0) then
                     scr(ii)=r(2)
                  else
      
!                    call binom(1.,2.,3.,4.,r(1),r(2),r(3),r(4),yy,scr(ii)) (now inlined)
!                    subroutine binom(x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy)
                     y1 = r(1)
                     y2 = r(2)
                     y3 = r(3)
                     y4 = r(4)
                     xxx=yy
c
                     yyy=missingflag
                     if ( .not. (y2 .gt. 1.e19 .or. y3 .gt. 1.e19) )then
c
                       wt1=(xxx-3.)/(-1.)
                       wt2=1.0-wt1
c
                       if (y4 .lt. 1.e19) then
                         yz22=wt1*(xxx-4.)/(-2.)
                         yz23=wt2*(xxx-4.)/(-1.)
                         yz24=(xxx-2.)*(xxx-3.)/2.
                       else
                         yz22=wt1
                         yz23=wt2
                         yz24=0.0
                       endif
c
                       if (y1 .lt. 1.e19) then
                         yz11=(xxx-2.)*(xxx-3.)/(2.)
                         yz12=wt1*(xxx-1.)/(1.)
                         yz13=wt2*(xxx-1.)/(2.)
                       else
                         yz11=0.0
                         yz12=wt1
                         yz13=wt2
                       endif
c
                       if (yz11 .eq. 0. .and. yz24 .eq. 0.) then
                         yyy=wt1*y2+wt2*y3
                       else
                         yyy=wt1*(yz11*y1+yz12*y2+yz13*y3)
     1                      +wt2*(yz22*y2+yz23*y3+yz24*y4)
                       endif
                     endif

                     scr(ii) = yyy
!                    return
!                    end

                  endif
               else 
                  scr(ii)=missingflag
               endif
            enddo
            xx=stax-fixm2
            if (xx .eq. 2.0) then
               staval=scr(2)
            else
!              call binom(1.,2.,3.,4.,scr(1),scr(2),scr(3),scr(4),xx,staval) (now inlined)
!              subroutine binom(x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy)
               y1 = scr(1)
               y2 = scr(2)
               y3 = scr(3)
               y4 = scr(4)
               xxx=xx
c
               yyy=missingflag
               if ( .not. (y2 .gt. 1.e19 .or. y3 .gt. 1.e19) )then
c
                 wt1=(xxx-3.)/(-1.)
                 wt2=1.0-wt1
c
                 if (y4 .lt. 1.e19) then
                   yz22=wt1*(xxx-4.)/(-2.)
                   yz23=wt2*(xxx-4.)/(-1.)
                   yz24=(xxx-2.)*(xxx-3.)/(+2.)
                 else
                   yz22=wt1
                   yz23=wt2
                   yz24=0.0
                 endif
c
                 if (y1 .lt. 1.e19) then
                   yz11=(xxx-2.)*(xxx-3.)/(+2.)
                   yz12=wt1*(xxx-1.)/(1.)
                   yz13=wt2*(xxx-1.)/(2.)
                 else
                   yz11=0.0
                   yz12=wt1
                   yz13=wt2
                 endif
c
                 if (yz11 .eq. 0. .and. yz24 .eq. 0.) then
                   yyy=wt1*y2+wt2*y3
                 else
                   yyy=wt1*(yz11*y1+yz12*y2+yz13*y3)
     1                +wt2*(yz22*y2+yz23*y3+yz24*y4)
                 endif
               endif

               staval = yyy
!              return
!              end

            endif
c
            if(staval.eq.missingflag)then
               print*,'hinterp_field_2d: staval = missingflag',il,jl
            endif

            flaps(il,jl,k) = staval
!           return
!           end

          else
            flaps(il,jl,k)=r_missing_data
          endif
        enddo
        enddo
      enddo
c
      return
      end
c

c
c===============================================================================
c
      subroutine gdtost(ab,ix,iy,stax,stay,staval,wrapped)
c
c *** Subroutine to return stations back-interpolated values(staval)
c        from uniform grid points using overlapping-quadratics.
c        gridded values of input array a dimensioned ab(ix,iy), where
c        ix = grid points in x, iy = grid points in y.  Station
c        location given in terms of grid relative station x (stax)
c        and station column.
c *** Values greater than 1.0e30 indicate missing data.
c
      dimension ab(ix,iy),r(4),scr(4)
      logical wrapped ! WNI added
      include 'bgdata.inc'
c_______________________________________________________________________________
c
      iy1=int(stay)-1
      if(stay.lt.1.0)iy1=1.0
      iy2=iy1+3
      ix1=int(stax)-1
      if(stax.lt.1.0) THEN 
        IF (wrapped) THEN   ! WNI
          ix1 = ix1 + ix
        ELSE ! WNI
          ix1=1.0
        ENDIF ! WNI
      endif ! WNI
      ix2=ix1+3
      staval=missingflag
      fiym2=float(iy1)-1
      fixm2=float(ix1)-1
      ii=0
      do iii=ix1,ix2
         i=iii
c
c ****** Account for wrapping around effect of global data at Greenwich.
c
         if (wrapped) then    ! WNI
            if (i .lt. 1) i=i+ix  
            if (i .gt. ix) i=i-ix
         endif

         ii=ii+1
         if (i .ge. 1 .and. i .le. ix) then 
            jj=0
            do jjj=iy1,iy2
               j=jjj
c
c ************ Account for N-S wrapping effect of global data.
c
               if (wrapped) THEN  !WNI
                  if (j .lt. 1) then
                     j=2-j
                     i=i-ix/2
                     if (i .lt. 1) i=i+ix
                     if (i .gt. ix) i=i-ix
                  endif
                  if (j .gt. iy) then
                     j=2*iy-j
                     i=i-ix/2
                     if (i .lt. 1) i=i+ix
                     if (i .gt. ix) i=i-ix
                  endif
               endif

               jj=jj+1
               if (j .ge. 1 .and. j .le. iy) then
                  r(jj)=ab(i,j)
               else
                  r(jj)=missingflag
               endif
            enddo
            yy=stay-fiym2
            if (yy .eq. 2.0) then
               scr(ii)=r(2)
            else
               call binom(1.,2.,3.,4.,r(1),r(2),r(3),r(4),yy,scr(ii))
            endif
         else 
            scr(ii)=missingflag
         endif
      enddo
      xx=stax-fixm2
      if (xx .eq. 2.0) then
         staval=scr(2)
      else
         call binom(1.,2.,3.,4.,scr(1),scr(2),scr(3),scr(4),xx,staval)
      endif
c
      if(staval.eq.missingflag)then
         print*,'gdtost: ',stax,stay,' staval = missingflag' 
      endif
      return
      end

c
c===============================================================================
c
      subroutine gdtost_i(ab,ix,iy,stax,stay,staval,wrapped)
c
c *** Subroutine to return stations back-interpolated values(staval)
c        from uniform grid points using overlapping-quadratics.
c        gridded values of input array a dimensioned ab(ix,iy), where
c        ix = grid points in x, iy = grid points in y.  Station
c        location given in terms of grid relative station x (stax)
c        and station column.
c *** Values greater than 1.0e30 indicate missing data.
c
      dimension ab(ix,iy),r(4),scr(4)
      logical wrapped ! WNI added
      include 'bgdata.inc'
c_______________________________________________________________________________
c
      iy1=int(stay)-1
      if(stay.lt.1.0)iy1=1.0
      iy2=iy1+3
      ix1=int(stax)-1
      if(stax.lt.1.0) THEN 
        IF (wrapped) THEN   ! WNI
          ix1 = ix1 + ix
        ELSE ! WNI
          ix1=1.0
        ENDIF ! WNI
      endif ! WNI
      ix2=ix1+3
      staval=missingflag
      fiym2=float(iy1)-1
      fixm2=float(ix1)-1
      ii=0
      do iii=ix1,ix2
         i=iii
c
c ****** Account for wrapping around effect of global data at Greenwich.
c
         if (wrapped) then    ! WNI
            if (i .lt. 1) i=i+ix  
            if (i .gt. ix) i=i-ix
         endif

         ii=ii+1
         if (i .ge. 1 .and. i .le. ix) then 
            jj=0
            do jjj=iy1,iy2
               j=jjj
c
c ************ Account for N-S wrapping effect of global data.
c
               if (wrapped) THEN  !WNI
                  if (j .lt. 1) then
                     j=2-j
                     i=i-ix/2
                     if (i .lt. 1) i=i+ix
                     if (i .gt. ix) i=i-ix
                  endif
                  if (j .gt. iy) then
                     j=2*iy-j
                     i=i-ix/2
                     if (i .lt. 1) i=i+ix
                     if (i .gt. ix) i=i-ix
                  endif
               endif

               jj=jj+1
               if (j .ge. 1 .and. j .le. iy) then
                  r(jj)=ab(i,j)
               else
                  r(jj)=missingflag
               endif
            enddo
            yy=stay-fiym2
            if (yy .eq. 2.0) then
               scr(ii)=r(2)
            else

!              call binom(1.,2.,3.,4.,r(1),r(2),r(3),r(4),yy,scr(ii)) (now inlined)
!              subroutine binom(x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy)
               y1 = r(1)
               y2 = r(2)
               y3 = r(3)
               y4 = r(4)
               xxx=yy
c
               yyy=missingflag
               if ( .not. (y2 .gt. 1.e19 .or. y3 .gt. 1.e19) )then
c
                 wt1=(xxx-3.)/(-1.)
                 wt2=1.0-wt1
c
                 if (y4 .lt. 1.e19) then
                   yz22=wt1*(xxx-4.)/(-2.)
                   yz23=wt2*(xxx-4.)/(-1.)
                   yz24=(xxx-2.)*(xxx-3.)/2.
                 else
                   yz22=wt1
                   yz23=wt2
                   yz24=0.0
                 endif
c
                 if (y1 .lt. 1.e19) then
                   yz11=(xxx-2.)*(xxx-3.)/(2.)
                   yz12=wt1*(xxx-1.)/(1.)
                   yz13=wt2*(xxx-1.)/(2.)
                 else
                   yz11=0.0
                   yz12=wt1
                   yz13=wt2
                 endif
c
                 if (yz11 .eq. 0. .and. yz24 .eq. 0.) then
                   yyy=wt1*y2+wt2*y3
                 else
                   yyy=wt1*(yz11*y1+yz12*y2+yz13*y3)
     1                +wt2*(yz22*y2+yz23*y3+yz24*y4)
                 endif
               endif

               scr(ii) = yyy
!              return
!              end

            endif
         else 
            scr(ii)=missingflag
         endif
      enddo
      xx=stax-fixm2
      if (xx .eq. 2.0) then
         staval=scr(2)
      else
!        call binom(1.,2.,3.,4.,scr(1),scr(2),scr(3),scr(4),xx,staval) (now inlined)
!        subroutine binom(x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy)
         y1 = scr(1)
         y2 = scr(2)
         y3 = scr(3)
         y4 = scr(4)
         xxx=xx
c
         yyy=missingflag
         if ( .not. (y2 .gt. 1.e19 .or. y3 .gt. 1.e19) )then
c
           wt1=(xxx-3.)/(-1.)
           wt2=1.0-wt1
c
           if (y4 .lt. 1.e19) then
             yz22=wt1*(xxx-4.)/(-2.)
             yz23=wt2*(xxx-4.)/(-1.)
             yz24=(xxx-2.)*(xxx-3.)/(+2.)
           else
             yz22=wt1
             yz23=wt2
             yz24=0.0
           endif
c
           if (y1 .lt. 1.e19) then
             yz11=(xxx-2.)*(xxx-3.)/(+2.)
             yz12=wt1*(xxx-1.)/(1.)
             yz13=wt2*(xxx-1.)/(2.)
           else
             yz11=0.0
             yz12=wt1
             yz13=wt2
           endif
c
           if (yz11 .eq. 0. .and. yz24 .eq. 0.) then
             yyy=wt1*y2+wt2*y3
           else
             yyy=wt1*(yz11*y1+yz12*y2+yz13*y3)
     1          +wt2*(yz22*y2+yz23*y3+yz24*y4)
           endif
         endif

         staval = yyy
!        return
!        end

      endif
c
C      if(staval.eq.missingflag)then
C         print*,'gdtost_i: staval = missingflag'
C      endif
      return
      end

c
c===============================================================================
c
      SUBROUTINE GDTOST2(A,IX,IY,STAX,STAY,STAVAL)
*  SUBROUTINE TO RETURN STATIONS BACK-INTERPOLATED VALUES(STAVAL)
*  FROM UNIFORM GRID POINTS USING OVERLAPPING-QUADRATICS.
*  GRIDDED VALUES OF INPUT ARRAY A DIMENSIONED A(IX,IY),WHERE
*  IX=GRID POINTS IN X, IY = GRID POINTS IN Y .  STATION
*  LOCATION GIVEN IN TERMS OF GRID RELATIVE STATION X (STAX)
*  AND STATION COLUMN.
*  VALUES GREATER THAN 1.0E30 INDICATE MISSING DATA.
*
      real A(IX,IY),R(4),SCR(4),stax,stay,staval
     +  ,fixm2,fiym2,yy,xx
      IY1=INT(STAY)-1
      IY2=IY1+3
      IX1=INT(STAX)-1
      IX2=IX1+3
      STAVAL=1E30
      FIYM2=FLOAT(IY1)-1
      FIXM2=FLOAT(IX1)-1
      II=0
      DO 100 I=IX1,IX2
      II=II+1
      IF(I.GE.1.AND.I.LE.IX) GO TO 101
      SCR(II)=1E30
      GO TO 100
101   JJ=0
      DO 111 J=IY1,IY2
      JJ=JJ+1
      IF(J.GE.1.AND.J.LE.IY) GO TO 112
      R(JJ)=1E30
      GO TO 111
112   R(JJ)=A(I,J)
111   CONTINUE
      YY=STAY-FIYM2
      CALL BINOM(1.,2.,3.,4.,R(1),R(2),R(3),R(4),YY,SCR(II))
100   CONTINUE
      XX=STAX-FIXM2
      CALL BINOM(1.,2.,3.,4.,SCR(1),SCR(2),SCR(3),SCR(4),XX,STAVAL)
      RETURN
      END
c
cc ------------------------------------------------------------------
c

      subroutine binom(x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy)
c
      include 'bgdata.inc'
      yyy=missingflag
      if ( .not. (y2 .gt. 1.e19 .or. y3 .gt. 1.e19) )then
c
         wt1=(xxx-x3)/(x2-x3)
         wt2=1.0-wt1
c
         if (y4 .lt. 1.e19) then
c           yz22=(xxx-x3)*(xxx-x4)/((x2-x3)*(x2-x4))
            yz22=wt1*(xxx-x4)/(x2-x4)
c           yz23=(xxx-x2)*(xxx-x4)/((x3-x2)*(x3-x4))
            yz23=wt2*(xxx-x4)/(x3-x4)
            yz24=(xxx-x2)*(xxx-x3)/((x4-x2)*(x4-x3))
         else
            yz22=wt1
            yz23=wt2
            yz24=0.0
         endif
c
         if (y1 .lt. 1.e19) then
            yz11=(xxx-x2)*(xxx-x3)/((x1-x2)*(x1-x3))
c           yz12=(xxx-x1)*(xxx-x3)/((x2-x1)*(x2-x3))
            yz12=wt1*(xxx-x1)/(x2-x1)
c           yz13=(xxx-x1)*(xxx-x2)/((x3-x1)*(x3-x2))
            yz13=wt2*(xxx-x1)/(x3-x1)
         else
            yz11=0.0
            yz12=wt1
            yz13=wt2
         endif
c
         if (yz11 .eq. 0. .and. yz24 .eq. 0.) then
            yyy=wt1*y2+wt2*y3
         else
            yyy=wt1*(yz11*y1+yz12*y2+yz13*y3)
     1         +wt2*(yz22*y2+yz23*y3+yz24*y4)
         endif
c
      endif

      return
      end
