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
       Subroutine Calc_Evap(imax,jmax,
     &                       Laps_u,
     &                       Laps_v,
     &                       Laps_T,
     &                       Laps_TD,
     &                       Laps_Evap,
     &                       IStatus)

C       Subroutine to calculate the pan evaporation rates for the 
C       Laps Grid. 10 km every Hour
C       Uses the Assumption of Non Radiation Limiting as In pg 166
C       Hydrology for Engineers  by Linsley, Kohler and Paulhus
C       2/5/93
C


      integer*4 imax,jmax
      Include 'soilm.inc'

      Real    Laps_u(Imax,Jmax)
      Real    Laps_v(Imax,Jmax)
      Real    Laps_T(Imax,Jmax)
      Real    Laps_TD(Imax,Jmax)
      Real    Laps_Evap(Imax,Jmax)
      Real    Val1, Val2, Val3
      Real    WindSpeed, TempDegF, DewPoint_DegF
      Integer*4 IStatus

      Do J = 1, Jmax
        Do I = 1, Imax
           WindSpeed = (Laps_u(I,J)**2 + Laps_v(I,J)**2)**0.5 !m/s
           WindSpeed = WindSpeed * 53.6979     ! miles per day
           TempDegF = 1.8 * (Laps_T(I,J) - 273.15) + 32.0  !Deg F
           DewPoint_DegF = 1.8 * (Laps_TD(I,J) - 273.15) + 32.0
           Val1 = ( 0.0041 * TempDegF + 0.676) ** 8
           Val2 = ( 0.0041 * DewPoint_DegF + 0.676) ** 8
           Val3 = (0.37 + 0.0041 * WindSpeed)

           Laps_Evap(I,J) = (Val1 -Val2)**0.88 * Val3
        Enddo
      Enddo
      IStatus = 1
c     Write(6,*) 'Completed Evaporation Calc'

      Return
      End
