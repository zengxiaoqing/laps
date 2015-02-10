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

      Real Function ZtoPsa(Z)

C*  This routine converts a height in meters into a pressure in a standard
C*  atmosphere in millibars.

      Implicit None

      include 'constants.inc' ! psamslpa, gamma

      Real T0,p0,p11,z11,c1,c2,z,Flag,Flg

      Data Flag,Flg/1e-30,199998./
      Data T0/288./
      Data c1,c2/5.256,14600./
      Data z11,p11/11000.,226.0971/

      p0 = psamslpa / 100.

      If (Z.gt.Flg) Then
          ZtoPsa=Flag
      Else If (Z.lt.z11) Then
          ZtoPsa=p0*((T0-gamma*Z)/T0)**c1
      Else
          ZtoPsa=p11*10.**((z11-Z)/c2)
      End If

      Return
      End
