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

      Real Function PsaToZ(psain)

C*  This routine converts a pressure in a standard atmosphere in millibars
C*  into a height in meters

      psa = psain
      scale_high = 8600.
      scale_low  = 6000.
      scale_height = scale_low + (scale_high - scale_low) * psa/1000.

!     Initial Guess
      scale_2 = scale_high * 0.35 + scale_height * 0.65
      scale_heights = log(1012./psa)
      psatoz = scale_heights * scale_2

      icount = 0

10    icount = icount + 1

      error = ztopsa(psatoz)/psa - 1.0

      correction = error * scale_height

c     write(6,1)icount,psa,psatoz,error,correction
1     format(i4,f8.1,f10.2,e15.4,f10.2)

      psatoz = psatoz + correction

      if(icount .ge. 20)then
          write(6,*)' Too many iterations in psatoz, input= ',psain
          write(6,1)icount,psa,psatoz,error,correction
          stop
      endif

      if(abs(error) .gt. 1e-6)goto10

c     write(6,1)icount,psa,psatoz,error,correction

      return
      end

