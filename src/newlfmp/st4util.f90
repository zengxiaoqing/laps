module st4util

implicit none

save

integer :: ncid

end module

!===============================================================================

subroutine get_st4_dims(fname,nx,ny,nz)

use st4util

implicit none

include 'netcdf.inc'

integer :: nx,ny,nz,nid,icode
character(len=*) :: fname
logical :: there

! Open stage IV file, and leave open for future use.

inquire(file=trim(fname),exist=there)
if (there) then
   icode=nf_open(trim(fname),nf_nowrite,ncid)
   if (ncid <= 0) then
      print*,'Could not open st4 file: ',trim(fname)
      stop
   endif
else
   print*,'Could not find nmm file: ',trim(fname)
   stop
endif

! Read stage IV grid dimensions.

icode=nf_inq_dimid(ncid,'x',nid)
icode=nf_inq_dimlen(ncid,nid,nx)
icode=nf_inq_dimid(ncid,'y',nid)
icode=nf_inq_dimlen(ncid,nid,ny)
icode=nf_inq_dimid(ncid,'record',nid)
icode=nf_inq_dimlen(ncid,nid,nz)

return
end

!===============================================================================

subroutine fill_st4_grid

use lfmgrid
use st4util
use constants

implicit none

include 'netcdf.inc'

integer :: nid,icode,i,j
real, parameter :: lat1=23.117,lon1=240.977

! Stage IV data uses the NCEP 255 grid, which is a user-defined grid.
! Grid parameter settings as of 11 May 2004 are given at the following URLs.
! http://www.emc.ncep.noaa.gov/mmb/ylin/pcpanl/QandA/gridchange/gridchange.htm
! http://wwwt.emc.ncep.noaa.gov/mmb/ylin/pcpanl/QandA/#GRIDINFO
! Polar Stereographic true at 60N
! Y-axis is parallel to 105W.
! Grid point spacing is 4.7625 at 60N
! lat1=23.117, lon1=240.977
! lat1, lon1 are the latitude and longitude (degrees east) of gridpoint (1,1)

nprojection='POLAR STEREOGRAPHIC'
ntruelat1 = 60.0
ntruelat2 = ntruelat1
nstdlon = -105.0
ngrid_spacingx = 4.7625
ngrid_spacingy = ngrid_spacingx
!icode=nf_get_att_real(ncid,NF_GLOBAL,'Dx',ngrid_spacingx)
!icode=nf_get_att_real(ncid,NF_GLOBAL,'Dy',ngrid_spacingy)

! Caculate latitude/longitude locations of Stage IV grid points
do j=1,ny
do i=1,nx
 call w3fb07(float(i),float(j),lat1,lon1,ngrid_spacingx*1000.,&
             nstdlon+360.,nlat(i,j),nlon(i,j))
 nlon(i,j) = nlon(i,j)-360.
enddo
enddo

! Read accumulated precip.
! And convert from mm to m.

icode=nf_inq_varid(ncid,'APCP',nid)
if (icode .eq. 0) then
   icode=nf_get_var_real(ncid,nid,npcp_tot)
endif
npcp_tot = npcp_tot/1000.

return
end

!===============================================================================
       SUBROUTINE W3FB07(XI,XJ,ALAT1,ALON1,DX,ALONV,ALAT,ALON)
!$$$   SUBPROGRAM  DOCUMENTATION  BLOCK
!
! SUBPROGRAM:  W3FB07        GRID COORDS TO LAT/LON FOR GRIB
!   PRGMMR: STACKPOLE        ORG: NMC42       DATE:88-04-05
!
! ABSTRACT: CONVERTS THE COORDINATES OF A LOCATION ON EARTH GIVEN IN A
!   GRID COORDINATE SYSTEM OVERLAID ON A POLAR STEREOGRAPHIC MAP PRO-
!   JECTION TRUE AT 60 DEGREES N OR S LATITUDE TO THE
!   NATURAL COORDINATE SYSTEM OF LATITUDE/LONGITUDE
!   W3FB07 IS THE REVERSE OF W3FB06.
!   USES GRIB SPECIFICATION OF THE LOCATION OF THE GRID
!
! PROGRAM HISTORY LOG:
!   88-01-01  ORIGINAL AUTHOR:  STACKPOLE, W/NMC42
!   90-04-12  R.E.JONES   CONVERT TO CRAY CFT77 FORTRAN
!
! USAGE:  CALL W3FB07(XI,XJ,ALAT1,ALON1,DX,ALONV,ALAT,ALON)
!   INPUT ARGUMENT LIST:
!     XI       - I COORDINATE OF THE POINT  REAL*4
!     XJ       - J COORDINATE OF THE POINT  REAL*4
!     ALAT1    - LATITUDE  OF LOWER LEFT POINT OF GRID (POINT 1,1)
!                LATITUDE <0 FOR SOUTHERN HEMISPHERE; REAL*4
!     ALON1    - LONGITUDE OF LOWER LEFT POINT OF GRID (POINT 1,1)
!                  EAST LONGITUDE USED THROUGHOUT; REAL*4
!     DX       - MESH LENGTH OF GRID IN METERS AT 60 DEG LAT
!                 MUST BE SET NEGATIVE IF USING
!                 SOUTHERN HEMISPHERE PROJECTION; REAL*4
!                   190500.0 LFM GRID,
!                   381000.0 NH PE GRID, -381000.0 SH PE GRID, ETC.
!     ALONV    - THE ORIENTATION OF THE GRID.  I.E.,
!                THE EAST LONGITUDE VALUE OF THE VERTICAL MERIDIAN
!                WHICH IS PARALLEL TO THE Y-AXIS (OR COLUMNS OF
!                THE GRID) ALONG WHICH LATITUDE INCREASES AS
!                THE Y-COORDINATE INCREASES.  REAL*4
!                   FOR EXAMPLE:
!                   255.0 FOR LFM GRID,
!                   280.0 NH PE GRID, 100.0 SH PE GRID, ETC.
!
!   OUTPUT ARGUMENT LIST:
!     ALAT     - LATITUDE IN DEGREES (NEGATIVE IN SOUTHERN HEMI.)
!     ALON     - EAST LONGITUDE IN DEGREES, REAL*4
!
!   REMARKS: FORMULAE AND NOTATION LOOSELY BASED ON HOKE, HAYES,
!     AND RENNINGER'S "MAP PROJECTIONS AND GRID SYSTEMS...", MARCH 1981
!     AFGWC/TN-79/003
!
! ATTRIBUTES:
!   LANGUAGE: CRAY CFT77 FORTRAN
!   MACHINE:  CRAY Y-MP8/832
!
!$$$
!
         DATA  RERTH /6.3712E+6/,PI/3.1416/
         DATA  SS60  /1.86603/
!
!        PRELIMINARY VARIABLES AND REDIFINITIONS
!
!        H = 1 FOR NORTHERN HEMISPHERE; = -1 FOR SOUTHERN
!
!        REFLON IS LONGITUDE UPON WHICH THE POSITIVE X-COORDINATE
!        DRAWN THROUGH THE POLE AND TO THE RIGHT LIES
!        ROTATED AROUND FROM ORIENTATION (Y-COORDINATE) LONGITUDE
!        DIFFERENTLY IN EACH HEMISPHERE
!
         IF (DX.LT.0) THEN
           H      = -1.0
           DXL    = -DX
           REFLON = ALONV - 90.0
         ELSE
           H      = 1.0
           DXL    = DX
           REFLON = ALONV - 270.0
         ENDIF
!
         RADPD  = PI    / 180.0
         DEGPRD = 180.0 / PI
         REBYDX = RERTH / DXL
!
!        RADIUS TO LOWER LEFT HAND (LL) CORNER
!
         ALA1 =  ALAT1 * RADPD
         RMLL = REBYDX * COS(ALA1) * SS60/(1. + H * SIN(ALA1))
!
!        USE LL POINT INFO TO LOCATE POLE POINT
!
         ALO1 = (ALON1 - REFLON) * RADPD
         POLEI = 1. - RMLL * COS(ALO1)
         POLEJ = 1. - H * RMLL * SIN(ALO1)
!
!        RADIUS TO THE I,J POINT (IN GRID UNITS)
!
         XX =  XI - POLEI
         YY = (XJ - POLEJ) * H
         R2 =  XX**2 + YY**2
!
!        NOW THE MAGIC FORMULAE
!
         IF (R2.EQ.0) THEN
           ALAT = H * 90.
           ALON = REFLON
         ELSE
           GI2    = (REBYDX * SS60)**2
           ALAT   = DEGPRD * H * ASIN((GI2 - R2)/(GI2 + R2))
           ARCCOS = ACOS(XX/SQRT(R2))
           IF (YY.GT.0) THEN
             ALON = REFLON + DEGPRD * ARCCOS
           ELSE
             ALON = REFLON - DEGPRD * ARCCOS
           ENDIF
         ENDIF
         IF (ALON.LT.0) ALON = ALON + 360.
!
      RETURN
      END
