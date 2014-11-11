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
        Subroutine Cloud_bogus_w (dx, cloud_type, height, nk       ! I
     1                           ,vv_to_height_ratio_Cu            ! I
     1                           ,vv_to_height_ratio_Sc            ! I
     1                           ,vv_for_St                        ! I
     1                           ,l_deep_vv                        ! I
     1                           ,w)                               ! O

!Original version October 1990.

!Modified May 1991 when we realized that a grid box is a lot bigger than
!any updraft.  We reduced the maximum vv in the parabolic profiles in cumulus
!clouds by a fairly large amount (30 m/s in 10 km Cu to 5 m/s), and the vv max
!for stratocumulus by a smaller amount (50 cm/s in 4-km Sc to 20 cm/s).

!  Modified June 2002 - Once again reduced the maximum cloud vv magnitude
!                       and made it dependent on grid spacing.  Also
!                       changed parabolic vv profile for cumuliform clouds
!                       to only go down to the cloud base, rather than 
!                       1/3rd of the cloud depth below base to try and
!                       and improve elevated convection cases.

!Can be used with either regular LAPS analysis grid or the cloud analysis grid.
        Implicit none
        Integer nk, cloud_type(nk)
        Real dx, height(nk), w(nk) ! dx is M and w is M/S

!The following specifies the maximum vv in two cloud types as functions
!of cloud depth.  Make parabolic vv profile, except for stratiform clouds,
!which get a constant value.  The values are tuned to give values that
!an NWP model would typically produce.  
        Real vv_to_height_ratio_Cu
        Real vv_to_height_ratio_Sc
        Real vv_for_St

        Real ratio, vv, Parabolic_vv_profile

        Integer k, k1, kbase, ktop
        Real zbase, ztop

        Logical l_deep_vv

!   Cloud Type      /'  ','St','Sc','Cu','Ns','Ac','As','Cs','Ci','Cc','Cb'/
!   Integer Value     0     1    2    3    4    5    6    7    8    9   10

!Zero out return vector.
        Do k = 1, nk
         w(k) = 0.
        End do

!Put in the vv's for cumuliform clouds (Cu or Cb) first.
        ratio = vv_to_height_ratio_Cu / dx
        Do k = 1, nk
         If (cloud_type(k) .eq. 3  .OR.  cloud_type(k) .eq. 10) then
          kbase = k
          Go to 10
         End if
        End do
        Go to 100

10      Do k = kbase, nk
         If(l_deep_vv)then ! We are using Adan's change
          If (cloud_type(k) .ne. 0) then ! change to the cloudtop by Adan
           ktop = k
          Else
           Go to 20
          End if
         else ! Older strategy with shallower parabolic profiles
          If (cloud_type(k) .eq. 3  .OR.  cloud_type(k) .eq. 10) then
           ktop = k
          Else
           Go to 20
          End if
         endif ! l_deep_vv
        End do

20      k1 = k          ! save our place in the column
        zbase = height(kbase)
        ztop  = height(ktop)
        Do k = 1, nk
         vv = Parabolic_vv_profile (zbase, ztop, ratio, height(k))
         If (vv .gt. 0.) then
          w(k) = vv
         Else
          w(k) = 0. ! could wipe out VV from another layer
         End if
        End do
        k1 = k1 + 1
        If (k1 .ge. nk) go to 100

!Try for another level of Cu.
        Do k = k1, nk
         If (cloud_type(k) .eq. 3  .OR.  cloud_type(k) .eq. 10) then
          kbase = k
          Go to 10
         End if
        End do

!Now do the stratocumulus or similar clouds (Sc, Ac, Cc, Ns).
100     ratio = vv_to_height_ratio_Sc/dx
        Do k = 1, nk
         If (cloud_type(k) .eq. 2  .OR.  cloud_type(k) .eq. 4  .OR.
     1     cloud_type(k) .eq. 5  .OR.  cloud_type(k) .eq. 9) then
          kbase = k
          Go to 110
         End if
        End do
        Go to 200

110     Do k = kbase, nk
         If (cloud_type(k) .eq. 2  .OR.  cloud_type(k) .eq. 4  .OR.
     1     cloud_type(k) .eq. 5  .OR.  cloud_type(k) .eq. 9) then
          ktop = k
         Else
          Go to 120
         End if
        End do

120     k1 = k          ! save our place in the column
        zbase = height(kbase)
        ztop  = height(ktop)
        Do k = 1, nk
         vv = Parabolic_vv_profile (zbase, ztop, ratio, height(k))
         If (vv .gt. w(k)) w(k) = vv
        End do
        k1 = k1 + 1
        If (k1 .ge. nk) go to 200       ! try for stratiform clouds

!Try for another level of Sc.
        Do k = k1, nk
         If (cloud_type(k) .eq. 2  .OR.  cloud_type(k) .eq. 4  .OR.
     1     cloud_type(k) .eq. 5  .OR.  cloud_type(k) .eq. 9) then
          kbase = k
          Go to 110
         End if
        End do

!Make sure there is non-zero vv wherever there are clouds of any kind.
!Also, return missing-data value for any non-bogussed vv value.
200     Do k = 1, nk
         If (cloud_type(k).ne.0 .AND. w(k).lt.vv_for_St) w(k) = vv_for_S
     1t
         If (w(k) .eq. 0.) w(k) = 1E37
        End do

        Return
        End

!-------------------------------------------------------------------
        Real Function Parabolic_vv_profile (zbase, ztop, ratio, z)
!The vertical velocity is zero at cloud top, peaks one third of the way up
!from the base, and extends below the base by one third of the cloud depth.

!  JUNE 2002 - No longer extending profile to below cloud base.

        Implicit none
        Real zbase, ztop, ratio, z
        Real depth, vvmax, vvspan, halfspan, height_vvmax, x

        depth = ztop - zbase
        If (depth .le. 0.) then
         Parabolic_vv_profile = 0.
         Return
        End if

        vvmax = ratio * depth
        vvspan = depth * 1.1    
        halfspan = vvspan / 2.
        height_vvmax = ztop - halfspan
        x = -vvmax/(halfspan*halfspan)

        Parabolic_vv_profile = x * (z-height_vvmax)**2 + vvmax

        Return
        End


! The variant below is called from 'get_radar_deriv.f/radar_bogus_w'

!-------------------------------------------------------------------
        Real Function Parabolic_vv_profile1 (zbase, ztop, ratio, z)
!The vertical velocity is zero at cloud top, peaks one third of the way up
!from the base, and extends below the base by one third of the cloud depth.

!  JUNE 2002 - No longer extending profile to below cloud base.

        Implicit none
        Real zbase, ztop, ratio, z
        Real depth, vvmax, vvspan, halfspan, height_vvmax, x

        depth = ztop - zbase
        If (depth .le. 0.) then
         Parabolic_vv_profile1 = 0.
         Return
        End if

        vvmax = ratio * depth
        vvspan = depth
        halfspan = vvspan / 2.
        height_vvmax = ztop - halfspan
        x = -vvmax/(halfspan*halfspan)

        Parabolic_vv_profile1 = x * (z-height_vvmax)**2 + vvmax

        Return
        End
