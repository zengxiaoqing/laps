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
        Subroutine RAzm_Lat_Lon_GM(r4_SLat,r4_SLon,r4_Range,
     1          r4_Azimuth,r4_TLat,r4_TLon,Status)
        include 'trigd.inc'
C***Given a range and azimuth from some site, calculate the latitude and
C   longitude. Range is considered to be a great circle arc length.

C       J. Wakefield    25 Apr 86       Original version.
C       S. Albers       8  May 90       Calculations done in Double Precision
C                                       This prevents 'singularity' at 0 AZ.
C       S. Albers      24  Aug 92       Change Status from .true. to 1

C Argument      I/O     Type                    Description
C --------      ---     ----    -----------------------------------------------
C SLat,SLon      I      R*4     Site location (deg).
C Range          I      R*4     Target range in km from the site.
C Azimuth        I      R*4     Target azimuth in degrees.
C TLat,TLon      O      R*4     Location of target (deg).
C Status         O      I*4     Standard system status.


!       Assumed earth radius is 6370000m
        real*8          km_per_Deg
        Parameter      (km_per_Deg=111.17747d0) ! km per degree of latitude

C***Parameter list variables
        Real          r4_SLat,r4_SLon,r4_Range,r4_Azimuth,r4_TLat,r4_T
     1Lon
        Real*8          SLat,SLon,Range,Azimuth,TLat,TLon
        Integer       Status

C***Local variables
        Real*8          Dist,CosDLon,DLon

C***Library symbols
!       EXTERNAL SS$_NORMAL, EDF__IVARG
!       Logical*4       LTest_Diag_GG

C***Begin RAzm_Lat_Lon_GM ------------------------------------------------------

C***Check input arguments
        If(Abs(r4_SLat).gt.90. .or. Abs(r4_SLon).gt.180. .or.
     1   r4_Range.lt..0   .or.
     1   r4_Azimuth.lt..0 .or. r4_Azimuth.gt.360.)Then
         Status=0
!        If(LTest_Diag_GG())Call Output_Diag_GG('RAzm_Lat_Lon_GM',Status)
         Return
        Else
         Status=1
        EndIf

        azimuth = r4_azimuth
        range = r4_range
        slat = r4_slat
        slon = r4_slon

C***Do it.  Note that calculations are done in degrees.
        Dist=Range/km_per_Deg

        TLat=dASinD(dCosD(Azimuth)*dSinD(Dist)*dCosD(SLat) 
     1           + dSinD(SLat)*dCosD(Dist))

        CosDLon=(dCosD(Dist) - dSinD(SLat)*dSinD(TLat)) 
     1            / (dCosD(SLat)*dCosD(TLat))
        If(Abs(CosDLon).gt.1.)CosDLon=Sign(1.d0,CosDLon)
        DLon=dACosD(CosDLon)
        If(Azimuth.ge..0.and.Azimuth.le.180.)Then
         TLon=SLon+DLon
        Else
         TLon=SLon-DLon
        EndIf

        r4_tlat = tlat
        r4_tlon = tlon

        Return
        End
