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
        subroutine windcnvrt(uwind,vwind,direction,speed)
c
c-----  Given wind components, calculate the corresponding speed and direction.
c-----  Hacked up from the windcnvrt_gm program.

c
C Argument      I/O     Type                    Description
C --------      ---     ----    -----------------------------------------------
C UWind          I      R*4     U-component of wind
C VWind          I      R*4     V-component of wind
C Direction      O      R*4     Wind direction (meteorological degrees)
C Speed          O      R*4     Wind speed (same units as input arguments)
c
C-----  If magnitude of UWind or VWind > 1E18, Speed and Direction set to -99.
c
        real*4          Flag
        Parameter      (Flag=1.e37)
c
        Real*4          UWind,VWind,Direction,Speed
c
        If(Abs(UWind).gt.1E18.or.Abs(VWind).gt.1E18)Then
         Speed=Flag
         Direction=Flag
        ElseIF(Uwind.eq.0.0.and.VWind.eq.0.0)Then
         Speed=0.0
         Direction=0.0                                  !Undefined
        Else
         Speed=SqRt(UWind*UWind+VWind*VWind)            !Wind speed
         Direction=57.2957795*(ATan2(UWind,VWind))+180. !Wind direction (deg)
        EndIf
c
        Return
        End
