cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
               Subroutine OPtrans(
     &		Nlev,
     &		P,
     &		T,
     &		Q,
     &		Pw,
     &		O3,
     &		Angle,
     &		Nchannels,
     &		Channels,
     &		Ptrans
     &                     )

*	Version 4.0  25 Oct 96
*	
*	Name:	Subroutine Optrans
*
*
*	Purpose: Computes the transmittance using the optical path
*		 regression algorithm.
*	
*	Input:  
*		I*4 Nlev      number of levels
*		R*4 P(Nlev)   Pressure levels
*		R*4 T(Nlev)   Temperature (K)
*		R*4 Q(Nlev)   Mixing ratio 
*		R*4 Pw(Nlev)  Precipitable water (cm) 
*		R*4 O3(Nlev)  Integrated Ozone Profile
*		R*4 Angle   Zenith Angle to use (Degrees)
*		I*4 Nchannels Number of channels to use
*		I*4 Channels(Nchannels) Array of channels for which to 
*		    compute transmittance.
*
*		from for042 wet coefficients
*		from for043 dry coefficients
*		from for044 ozo coefficients
*		From for061 wet control file
*		From for062 dry control file
*		From for063 ozo control file
*
*
*	Output: R*4 Ptrans(0:Nlev,Nchannels) Transmittances at Nlev levels
*
*	Author:	Thomas J. Kleespies
*		Physics Branch
*		Satellite Research Laboratory
*		Office of Research and Applications
*		NOAA/NESDIS
*		301-763-8136
*		301-763-8108 FAX
*
*	Mailing Address:
*		810 NSC E/RA-14
*		NOAA/NESDIS
*		Washington, D.C. 20233
*
*	Email:
*		Kleespies@NZMS.WWB.NOAA.GOV
*
*

	Implicit None
        save

	Include 'constants_optran.inc'

*	Input
	Integer Nlev  	! number of levels
	Real P(Nlev)  	! Pressure (mb)
	Real T(Nlev)	! temperature (K)
	Real Q(Nlev)	! mixing ratio
	Real Pw(Nlev) 	! precipitable water profile (cm)
	Real O3(Nlev)   ! integrated ozone profile
	Real Angle	! zenith angle (degrees)
	Integer Nchannels ! number of channels to process
	Integer Channels(Nchannels) ! channel list

*	Output
	! estimated transmittance in pressure space
	Real Ptrans(0:Nlevel,Nchannels)

*	Local
		
	Real*8 Bw(Nwet+1, Nw ,Nchan)  ! Vector of wet regression coefficients
	Real*8 Bd(Ndry+1, Nw ,Nchan)  ! Vector of dry regression coefficients
	Real*8 Bo(NOzo+1, Nw ,Nchan)  ! Vector of dry regression coefficients

	Real WW(0:Nw)		! wet absorber amount space
	Real AA(0:Nw)		! dry absorber amount space
	Real OO(0:Nw)		! Ozo absorber amount space
	Real PwMax 
        data PwMax /13.87981/   ! max wet absorber at WW(Nw)
	Real PMax  
        data PMax /2100./       ! max dry absorber at AA(Nw)
	Real OzMax  
        data OzMax /1.1403912/  ! max ozo absorber at OO(Nw)
c        Real PwMax /12.33761/   ! maximum value for wet absorber space
c	Real PMax  /2000./	! max dry absorber at AA(Nw)
c        Real OzMax /153.4007/   ! maximum value for Ozo absorber space
	Real SlantPw(Nlevel)	! wet absorber on a generalized slant path
	Real SlantP(Nlevel)	! dry absorber on a generalized slant path
	Real SlantOz(Nlevel)	! Ozo absorber on a generalized slant path

	Real ODWet(Nlevel,Nchan)    ! wet optical depth on pressure levels
	Real ODDry(Nlevel,Nchan)    ! dry optical depth on pressure levels
	Real ODOzo(Nlevel,Nchan)    ! ozo optical depth on pressure levels

	Real Dtr 
        data Dtr /1.74532952e-2/ ! degrees to radians conversion

	Real Secant		! secant of angle


	Integer Level

	Integer i ! local index
	Integer Ichan 	! channel index for read
	Integer Ipred	! predictor index for read

	Integer Pred_Index_Wet(0:Nwet,Nchan) ! wet predictors
	Integer Pred_Index_Dry(0:NDry,Nchan) ! dry predictors
	Integer Pred_Index_Ozo(0:NOzo,Nchan) ! dry predictors

c	Integer IOs ! io return from open


*  This may need a SAVE statement on some systems
	Logical First 
        data First /.true./     ! flag to initialize coefficients and such
				! first time routine is called.
c     LAPS APPS
        character*200 fname
        integer len

	If(First) Then ! first time called, initialize variables and parameters


           call get_directory ('static',fname,len)

	
	Open(42,file=fname(1:len)//'optranlib/wet_coeff.dat',
     &          form='unformatted')
	Open(43,file=fname(1:len)//'optranlib/dry_coeff.dat',
     &          form='unformatted')
	Open(44,file=fname(1:len)//'optranlib/ozo_coeff.dat',
     &          form='unformatted')
	

	 Read(42) Bw	! read wet coefficients
	 Close(42)
	 Read(43) Bd	! read dry coefficients
	 Close(43)
	 Read(44) Bo	! read ozo coefficients
	 Close(44)

*	Read the control files

 3100     format(14i3)

         Open(61,
     &		file=fname(1:len)//'optranlib/Wet_Control_File.dat',
     &          form='Formatted',
     &          status='old',
     &          err=310)

         Do Ichan = 1 , Nchan
          Read(61,3100)  (Pred_Index_Wet(Ipred,Ichan), Ipred=0,Nwet)
         EndDo

         Close(61)

         Open(62,
     &		file=fname(1:len)//'optranlib/Dry_Control_File.dat',
     &          form='Formatted',
     &          status='old',
     &          err=310)

         Do Ichan = 1 , Nchan
          Read(62,3100)  (Pred_Index_Dry(Ipred,Ichan), Ipred=0,NDry)
         EndDo

         Close(62)

         Open(63,
     &		file=fname(1:len)//'optranlib/Ozo_Control_File.dat',
     &          form='Formatted',
     &          status='old',
     &          err=310)

         Do Ichan = 1 , Nchan
          Read(63,3100)  (Pred_Index_Ozo(Ipred,Ichan), Ipred=0,NOzo)
         EndDo

         Close(63)

*	Construct the absorber amount spaces for wet and dry cases

	 Call Make_Wlevels(1, PwMax, 1.0e-7, WW)
         Call Make_Wlevels(1, Pmax ,   0.01, AA)
         Call Make_Wlevels(1, Ozmax ,2.0e-5, OO)

	 First = .False.  ! assure that we don't come here again

	EndIf ! first


*	This is where the actual work starts

*	It may be desirable to put in a check to assure that Nlev <= Nlevel

	Secant = 1.0 / Cos(Angle*DTR)
	Do level= 1 , Nlev	! compute slant path absorber amounts
         SlantP(level) = P(level)*Secant - 0.005
	 SlantPw(level) = Pw(level)*Secant
	 SlantOz(level) = O3(level)*Secant
	EndDo

*	Get Wet Optical Depth profile
           Call Optrans_Species     (
     &		Nlev,
     &		Nwet+1,
     &		Nwet,
     &          Bw,
     &          WW,
     &          P,
     &          T,
     &		Q,
     &          SlantPw,
     &          Nchannels,
     &          Channels,
     &		Pred_Index_Wet,
     &          Odwet
     &                          )
*	Get Dry Optical Depth profile
           Call Optrans_Species     (
     &		Nlev,
     &		Ndry+1,
     &		Ndry,
     &          Bd,
     &          AA,
     &          P,
     &          T,
     &		Q,
     &          SlantP,
     &          Nchannels,
     &          Channels,
     &		Pred_Index_Dry,
     &          Oddry
     &                          )
*	Get Ozone Optical Depth profile
           Call Optrans_Species     (
     &		Nlev,
     &		NOzo+1,
     &		NOzo,
     &          Bo,
     &          OO,
     &          P,
     &          T,
     &		Q,
     &          SlantOz,
     &          Nchannels,
     &          Channels,
     &		Pred_Index_Ozo,
     &          OdOzo
     &                          )

*	Compute total transmittance profile

	Do i = 1 ,  Nchannels ! now process only requested channels
	 Ichan = Channels(i)
	 Ptrans(0,i) = 1.0
	 Do level = 1 , Nlev
	  Ptrans(level,i) = Ptrans(level-1,i)
     &		    *exp(-(
     &			   OdWet(level,i)
     &			  +OdDry(level,i)
     &			  +OdOzo(level,i)
     &			) )
	 EndDo ! Nlev

	EndDo  ! channel

*	That's all folks!

	Return

  310   Continue
	Write(6,*) 'Error on control file open '
c	Write(6,*) IOs
	Stop
	end
	Subroutine Optrans_Species	(
     &		Nlev,
     &		NCoeff,
     &		Npred,
     &		B,
     &		WW,
     &		P,
     &		T,
     &		Q,
     &		Absorber,
     &		NChannels,
     &		Channels,
     &		Pred_Index,
     &		OD
     &					)

*
*	Version 4: 25 October 1995
*	
*	Name: OPtrans_Species
*
*
*	Purpose:
*	This is the routine that determines the optical depth profile
*	given the pressure profile, temperature profile, and absorber
*	profile.  This is generalized such that it will work with any
*	species.
*	
*	Input:
*	Integer Nlev		! number of levels in input atmosphere
*	Integer Ncoeff		! number of coefficients
*	Integer Npred		! number of predictors
*	Real B(Ncoeff,Nw,Nchan)	! regression coefficients for this species
*	Real WW(0:Nw)		! Absorber space for this species
*	Real P(Nlev)		! atmospheric pressure profile
*	Real T(Nlev)		! atmospheric temperature profile
*	Real Q(Nlev)		! atmospheric mixing ratio profile
*	Real Absorber(Nlev)	! slant path absorber profile
*	Integer Nchannels	! Number of channels 
*	Integer Channels(Nchannels) ! vector of channels to use
*       Integer Pred_Index(0:Npred,Nchan)  ! index set of predictors to use
*					! for each channel.
*
*	Output variables
*	Real OD(Nlev,Nchannels) ! Optical depth in pressure space
*
*	Author:	Thomas J. Kleespies
*		Physics Branch
*		Satellite Research Laboratory
*		Office of Research and Applications
*		NOAA/NESDIS
*		301-763-8136
*		301-763-8108 FAX
*
*	Mailing Address:
*		810 NSC E/RA-14
*		NOAA/NESDIS
*		Washington, D.C. 20233
*
*	Email:
*		Kleespies@NZMS.WWB.NOAA.GOV
*
	Implicit None

	Include 'constants_optran.inc'

*	Input variables
	Integer Nlev		! number of levels in input atmosphere
	Integer Ncoeff		! number of coefficients
	Integer Npred		! number of predictors
	Real*8 B(Ncoeff,Nw,Nchan)! regression coefficients
	Real WW(0:Nw)		! absorber space
	Real P(Nlev)		! atmospheric pressure profile
	Real T(Nlev)		! atmospheric temperature profile
	Real Q(Nlev)		! atmospheric mixing ratio profile
	Real Absorber(Nlev)	! atmospheric absorber profile
	Integer Nchannels	! Number of channels 
	Integer Channels(Nchannels) ! vector of channels to use
        Integer Pred_Index(0:NPred,Nchan) ! index set of predictors to use
     &                                    ! for each channel.

*	Output variables
	Real OD(Nlevel,Nchannels)	! Optical Depth

*	Local variables

	Real*8 Sum		! summation 

	Real Average_Absorber	! average absorber amount

        Real Kw(Nw)             ! k = predictand in absorber space
	
	Real Kp(Nlevel)		! absorption coeff in pressure space

	Real DP(Nlevel)		! absorber within layer
	
	Real Factor(Nlevel)	! linear interpolation factor

	Integer i,level,kt	! utility variables
	Integer Ichan,Jchan
	Integer K2,K1
	Integer Ipred

	Integer MaxWlevels	! Actual number of levels used in W space

	Integer Wlevels(Nw)	! index array pointing to the W space
				! levels actually used

        Integer M		! index for interpolation

        Logical Search		! .true. for first call to Search_Plevel_Linear

        Real*8 xx(Nw,15) ! Presently hard coded for 15 potential predictors 

*	Find values in W space that bracket the absorber levels
	Search = .True.

	Do level = 1 , Nlev
	 Kt = 2*(level-1) + 1
	 
	 If(level .eq. 1) Then
	  Average_Absorber = Absorber(level) / 2.
	 Else
	  Average_Absorber = (Absorber(level)+Absorber(level-1))/2.
	 EndIf

	 Call Search_Plevel_Linear(Search,WW,Average_Absorber,Nw,M)
	 Wlevels(kt) = M
	 Wlevels(kt+1) = M+1
	 Factor(level) = (WW(M+1)-Average_Absorber)/(WW(M+1)-WW(M))
	 If(level .eq.1) Then
	  Dp(Level) = Absorber(level)
	 Else
	  Dp(Level) = (Absorber(level)-Absorber(level-1))
	 EndIf
	EndDo

	MaxWlevels = kt + 1

*	Compute the predictors for this atmosphere

        Call Get_Predictors_All(
     &          P,
     &          Absorber,
     &          T,
     &		Q,
     &		Nlev,
     &          WW,
     &          MaxWlevels,
     &		Wlevels,
     &		xx
     &                            )

	Do  300 Jchan = 1 , Nchannels

	 Ichan = Channels(Jchan)

	 If(Pred_Index(0,Ichan) .lt. 1) Then
	  Do level = 1 , Nlev ! no absorption for this channel
	   OD(level,Jchan) = 0.0
	  EndDo
	  Go To 300
	 EndIf

*	Zero Kw array

	 Do i = 1 , Nw
	  Kw(i) = 0.0
	 EndDo

*	Compute k in W space

	 Do 200 i = 1 , MaxWlevels 

	  K1 = Wlevels(i)
		
	  If(kw(K1) .ne. 0.0) Go to 200 ! already calculated

           Sum = B(1,K1,Ichan)

           Do Ipred = 2 , NCoeff
            Sum = Sum + B(Ipred,K1,Ichan)
     &                      * xx(K1,Pred_Index(Ipred-1,Ichan))
           EndDo

	   Kw(K1) = Sum

  200    Continue

*       Interpolate back to P space

	 Do level = 1 , Nlev
	  Kt = 2*(level-1) + 1
	  K1 = Wlevels(kt)
	  K2 = wlevels(kt+1)
	  Kp(level) = Kw(K2) - (Kw(k2)-Kw(k1))*Factor(level)
	  If(Kp(level) .lt. 0.0) Kp(level) = 0.0
	 EndDo

*	Compute optical Depth

         OD(1,Jchan) = 0.0

         Do level = 1 , Nlev

	  OD(level,Jchan) = Dp(Level)*kp(level) 

	 EndDo

  300   Continue ! Jchan

	Return
	End
	Subroutine Get_Predictors_All(
     &		Press,
     &		Absorb,
     &		Temp,
     &		MixR,
     &		Nlev,
     &		WW,
     &		MaxWlevels,
     &		Wlevels,
     &		x )
*	
*	Name: Subroutine Get_Predictors_All
*
*	Version 4: 25 October 95
*
*	Purpose:
*	Gets all predictors out of set of 14
*	
*	Input:
*	Press(Nlev)	Atmospheric pressure levels
*	Absorb(Nlev)	Atmospheric absorber profile on slant path
*	Temp(Nlev)	Atmospheric temperature profile
*	MixR(Nlev)	Atmospheric mixing ratio profile
*	Nlev		Number of atmospheric levels
*	WW(0:Nw)	Standard absorber space
*	MaxWLevels	Number of levels to be filled in absorber space
*	Wlevels(MaxWlevels)index of absorber space levels to be interpolated to 
*
*	Output:
*
*	X - various predictors
*
*	Note in the following description 'x' is multiplication,
*					  '^' is exponention
*					  '*' is a "star'd" quantity
*	X(1)	T	
*	X(2)	P	
*	X(3)	T^2
*	X(4)	P^2
*	X(5)    TxP
*	X(6)	T^2xP
*	X(7)	TxP^2
*	X(8)	T^2xP^2
*	X(9)	T*
*	X(10)	P*
*	X(11)	T**
*	X(12)	P**
*	X(13)	T***
*	X(14)	P***	
*	X(15)   q
*
*	Author:	Thomas J. Kleespies
*		Physics Branch
*		Satellite Research Laboratory
*		Office of Research and Applications
*		NOAA/NESDIS
*               301-763-8136
*               301-763-8108 FAX
*
*	Mailing Address:
*		810 NSC E/RA-14
*		NOAA/NESDIS
*		Washington, D.C. 20233
*
*	Email:
*		Kleespies@NZMS.WWB.NOAA.GOV
*
*

	Implicit None

	Include 'constants_optran.inc'

*	Input
	Integer Nlev		! number of atmospheric levels
	Real Press(Nlev)	! Atmospheric pressure levels
	Real Absorb(Nlev)		! atmospheric absorber profile
	Real Temp(Nlev)		! atmospheric temperature profile
	Real MixR(Nlev)		! atmospheric mixing ratio profile
	Real WW(0:Nw)		! standard Absorber profile

	Integer MaxWLevels	! number of levels to be filled in absorber space
	Integer Wlevels(MaxWlevels) ! index of absorber space levels to be used

*	Output

	Real*8 x(Nw,15)

	Real*8 delA

	Real*8 s1,s2,s3,s4,s5,s6

*	Local

	Integer i,j,k,kt,level

	Logical Search_Absorber
	Integer M1

	Real*8 Linear_Interpolation_D ! function to do DP linear interpolation

	Real*8 xx1(Nlevel)	! these are the star'd quantities
	Real*8 xx2(Nlevel)	! in pressure space
	Real*8 xx3(Nlevel)
	Real*8 xx4(Nlevel)
	Real*8 xx5(Nlevel)
	Real*8 xx6(Nlevel)

	Integer MaxLevel
	Parameter (MaxLevel=100)
	Real*8 P(MaxLevel)
	Real*8 T(MaxLevel)
	Real*8 Absorber(MaxLevel)
	Real*8 Q(MaxLevel)

	Search_Absorber = .True.

*	Zero arrays

	Do i = 1 , Nw
	 Do j = 1 , 14
	   X(i,j) = 0.0D0
	 EndDo
	EndDo	 

*	Make atmospheric variables double precision

	If(Nlev .gt. MaxLevel) Then
	 Write(6,*) 'Too many levels in input profile'
	 Write(6,*) MaxLevel , ' is maximum allowed'
	 Write(6,*) 'Change parameter MaxWlevels in Get_Predictors_All'
	 Return
	EndIf

	Do level = 1 , Nlev
	 T(level) = Temp(level)
	 P(level) = Press(level)
	 Absorber(level) = Absorb(level)
	 Q(level) = MixR(level)
	EndDo

*	First compute star'd terms for the specified atmosphere
*	i.e. on pressure levels

	If(Absorber(1) .gt. 0.0) Then
	 xx1(1)  = T(1) / 2. ! if we treat T(0) = 0, Absorber(0) = 0
	 xx2(1)  = P(1) / 2. ! if we treat P(0) = 0, Absorber(0) = 0
	 xx3(1) = T(1) / 2. ! if we treat T(0) = 0, Absorber(0) = 0
	 xx4(1) = P(1) / 2. ! if we treat P(0) = 0, Absorber(0) = 0
	 xx5(1) = T(1) / 2. ! if we treat T(0) = 0, Absorber(0) = 0
	 xx6(1) = P(1) / 2. ! if we treat P(0) = 0, Absorber(0) = 0
 	Else
	 xx1(1) = 0.0
	 xx2(1) = 0.0
	 xx3(1) = 0.0
	 xx4(1) = 0.0
	 xx5(1) = 0.0
	 xx6(1) = 0.0
	EndIf

	s1 = 0.0
	s2 = 0.0
	s3 = 0.0
	s4 = 0.0
	s5 = 0.0
	s6 = 0.0

	Do k = 2 , Nlev
	 delA = Absorber(k) - Absorber(k-1)
	 s1 = s1 + (T(k)+T(k-1))*delA	! T*
	 s2 = s2 + (P(k)+P(k-1))*delA	! P*
	 s3 = s3 + (Absorber(k)*T(k) + Absorber(k-1)*T(k-1))*delA !T**
	 s4 = s4 + (Absorber(k)*P(k) + Absorber(k-1)*P(k-1))*delA !P**
	 s5 = s5 + (Absorber(k)*Absorber(k)*T(k) + 
     &		    Absorber(k-1)*Absorber(k-1)*T(k-1))*delA !T***
	 s6 = s6 + (Absorber(k)*Absorber(k)*P(k) + 
     &		    Absorber(k-1)*Absorber(k-1)*P(k-1))*delA !P***

	 If(s1 .eq. 0.0) Then 
	   xx1(k) = 0.0
	 Else
	   xx1(k) = s1/(2.0*Absorber(k))
	 EndIf
	 If(s2 .eq. 0.0) Then 
	   xx2(k) = 0.0
	 Else
	   xx2(k) = s2/(2.0*Absorber(k))
	 EndIf
	 If(s3 .eq. 0.0) Then 
	   xx3(k) = 0.0
	 Else
	   xx3(k) = s3/(Absorber(k)*Absorber(k))
	 EndIf
	 If(s4 .eq. 0.0) Then 
	   xx4(k) = 0.0
	 Else
	   xx4(k) = s4/(Absorber(k)*Absorber(k))
	 EndIf
	 If(s5 .eq. 0.0) Then 
	   xx5(k) = 0.0
	 Else
	   xx5(k) = 1.5*s5/(Absorber(k)*Absorber(k)*Absorber(k))
	 EndIf
	 If(s6 .eq. 0.0) Then 
	   xx6(k) = 0.0
	 Else
	   xx6(k) = 1.5*s6/(Absorber(k)*Absorber(k)*Absorber(k))
	 EndIf

	EndDo

*	Next compute for the standard absorber levels
	 
	Do 200 k = 1 , MaxWlevels	

	  kt = Wlevels(k)
	  If(x(kt,1) .ne. 0) Go to 200 ! already did this level


          Call Search_Wlevel_Linear(
     &		Search_Absorber,Absorb,Nlev,WW(kt),M1)

          x(kt,1)    = Linear_Interpolation_D(Absorber(M1),T(M1), 	! T
     &                                  Absorber(M1+1),T(M1+1),
     &                                  DBLE(WW(kt)))

          x(kt,2)    = Linear_Interpolation_D(Absorber(M1),P(M1), 	! P
     &                                  Absorber(M1+1),P(M1+1),
     &                                  DBLE(WW(kt)))

	  x(kt,3)    =  x(kt,1)*x(kt,1) 			! T**2
	  x(kt,4)    =  x(kt,2)*x(kt,2) 			! P**2
	  x(kt,5)    =  x(kt,1)*x(kt,2) 			! T*P
	  x(kt,6)    =  x(kt,3)*x(kt,2) 			! T**2 P
	  x(kt,7)    =  x(kt,1)*x(kt,4) 			! T P**2
	  x(kt,8)    =  x(kt,3)*x(kt,4) 			! T**2 P**2

          x(kt,9)= 	Linear_Interpolation_D(Absorber(M1),xx1(M1),
     &                                  Absorber(M1+1),xx1(M1+1),
     &                                  DBLE(WW(kt)))

          x(kt,10)= 	Linear_Interpolation_D(Absorber(M1),xx2(M1),
     &                                  Absorber(M1+1),xx2(M1+1),
     &                                  DBLE(WW(kt)))

          x(kt,11)= 	Linear_Interpolation_D(Absorber(M1),xx3(M1),
     &                                  Absorber(M1+1),xx3(M1+1),
     &                                  DBLE(WW(kt)))

          x(kt,12)= 	Linear_Interpolation_D(Absorber(M1),xx4(M1),
     &                                  Absorber(M1+1),xx4(M1+1),
     &                                  DBLE(WW(kt)))

          x(kt,13)= 	Linear_Interpolation_D(Absorber(M1),xx5(M1),
     &                                  Absorber(M1+1),xx5(M1+1),
     &                                  DBLE(WW(kt)))

          x(kt,14)= 	Linear_Interpolation_D(Absorber(M1),xx6(M1),
     &                                  Absorber(M1+1),xx6(M1+1),
     &                                  DBLE(WW(kt)))

          x(kt,15)= 	Linear_Interpolation_D(Absorber(M1),Q(M1),
     &                                  Absorber(M1+1),Q(M1+1),
     &                                  DBLE(WW(kt)))

  200   Continue	! MaxWlevels loop

	Return
	End
	Real*8 Function Linear_Interpolation_D(x1,y1,x2,y2,x)
*
*	Performs Simple Linear Interpolation
*	
******* WARNING: Don't Forget to declare REAL in calling routine **********
*
*	Input-
*
*	x1,y1,x2,y2 - Pairs bounding the interval
*	x - independent variable - we want to find the y value associated to x
*
*	Returned value: dependent variable associated with x
*
*       3 Dec 90
*
*
*	Author:	Thomas J. Kleespies
*		Physics Branch
*		Satellite Research Laboratory
*		Office of Research and Applications
*		NOAA/NESDIS
*               301-763-8136
*               301-763-8108 FAX
*
*	Mailing Address:
*		810 NSC E/RA-14
*		NOAA/NESDIS
*		Washington, D.C. 20233
*
*	Email:
*		Kleespies@NZMS.WWB.NOAA.GOV

	Implicit NONE

	Real*8 x1,y1,x2,y2,x
	Real*8 Denominator

	Denominator = x2-x1

	If(Denominator .ne. 0.0) Then
	   Linear_Interpolation_D = y2 - (y2-y1)*(x2-x)/Denominator
	Else
	   Linear_Interpolation_D = y2
	EndIf

	Return
	End
	Subroutine Search_Wlevel_Linear(First,Absorber,Nlev,WW,M)
*	
*	Name:	Subroutine Search_Wlevel_Linear
*
*	Version 4: 25 October 1995
*
*	Purpose:
*	Searches the atmospheric absorber profile for
*	the 'best' bracketing values for subsequent linear interpolation.
*	
*	Input:
*	First		set .true. for first call for a given profile
*	Absorber(Nlev)  Absorber amount profile
*	Nlev		number of levels in Absorber
*	WW		value of absorber to bracket
*
*	Output:
*	M	 	returned value of M which represents
*			the first point in the interpolation array.
*
*	Author:	Thomas J. Kleespies
*		Physics Branch
*		Satellite Research Laboratory
*		Office of Research and Applications
*		NOAA/NESDIS
*		301-763-8136
*		301-763-8108 FAX
*
*	Mailing Address:
*		810 NSC E/RA-14
*		NOAA/NESDIS
*		Washington, D.C. 20233
*
*	Email:
*		Kleespies@NZMS.WWB.NOAA.GOV
*

	Implicit None

	Logical First
	Integer Nlev		! number of levels
	Real Absorber(Nlev)	! Absorber amount profile
	Real WW			! value of absorber to bracket
	Integer M		! returned value of M which represents
				! the first point in the interpolation array.

	Integer i,j		! utility variables
	Integer Start_Level	! starting level for search

	If(First) Then
	  Start_Level = 1
	  First = .false.
	Else
	  Start_Level = M  ! M should be retained by calling routine
	EndIf

	Do i = Start_Level , Nlev
	 If(WW .lt. Absorber(i)) Then
	  j = i
	  Go to 110
	 EndIf
	EndDo

	j = Nlev

  110   Continue

	If(j .le. 1) Then 	! extrapolate above atmosphere
	  M = 1 
	Else If (j .eq. Nlev) Then ! extrapolate below atmosphere
	  M = Nlev-1
	Else
	  M = j-1		! bracket within atmosphere
	EndIf

	Return
	End
	Subroutine Search_Plevel_Linear(First,WW,Absorber,Wlevels,M)
*	
*	Name: 	Subroutine Search_Plevel_Linear
*
*	Version 4: 25 October 1995
*		
*
*	Purpose:
*	Searches the standard absorber profile for
*	the 'best' bracketing values for subsequent linear interpolation.
*	
*	
*	Input:
*	First		set .true. for first call for a given profile
*	WW(0:Nw)	Standard absorber profile
*	Absorber	Absorber amount to bracket
*	Wlevels		Actual number of levels used in WW space
*
*	Output:
*	M		returned value of M which represents
*			the first point in the interpolation array.
*
*	Author:	Thomas J. Kleespies
*		Physics Branch
*		Satellite Research Laboratory
*		Office of Research and Applications
*		NOAA/NESDIS
*		301-763-8136
*		301-763-8108 FAX
*
*	Mailing Address:
*		810 NSC E/RA-14
*		NOAA/NESDIS
*		Washington, D.C. 20233
*
*	Email:
*		Kleespies@NZMS.WWB.NOAA.GOV
*
*

	Implicit None

	Include 'constants_optran.inc'

	Logical First 		! set .true. for initial invocation
	Real WW(0:Nw)		! standard absorber profile
	Real Absorber		! absorber amount to bracket
	Integer Wlevels		! Actual number of levels used in WW space
	Integer M		! returned value of M which represents
				! the first point in the interpolation array.

	Integer i,j		! utility variables

	Integer Start_Level	! starting level for search

	If(First) Then
	  Start_Level = 0
	  First = .False.
	Else
	  Start_Level = M  ! M should be retained by calling routine
	EndIf

	Do i = Start_Level , Wlevels
	 If(Absorber .lt. WW(i)) Then
	  j = i
	  Go to 110
	 EndIf
	EndDo

	j = Wlevels

  110   Continue

	If (j .gt. 1) Then
	 M = j-1 
	Else
	 M = 1
	EndIf

	Return
	End
	Subroutine Make_Wlevels(Nprofiles,MaxAbsorber,MinAbsorber,W)
*	
*	Name: Subroutine Make_Wlevels
*
*	Version 4: 25 October 1995
*
*	Purpose:
*	Makes the absorber space of Nw levels,
*	identified as W(0:Nw).  The most important value for
*	defining the parameter 'alpha' is W(200), which must be
*	at least the maximum possible absorber amount. The code 
*	is structured to produce different standard profiles for each channel.
*	
*	Input:
*	Nprofiles	number of profiles to use
*			this is the dimension of MaxAbsorber
*			and the second dimension of W.
*	MaxAbsorber(Nprofiles)	! absorber amount at W(200)
*	MinAbsorber	! absorber amount at W(1), W(0) == 0.0
*
*	Output:
*	W(0:Nw,Nprofiles) ! standard absorber profiles
*
*	Author:	Thomas J. Kleespies
*		Physics Branch
*		Satellite Research Laboratory
*		Office of Research and Applications
*		NOAA/NESDIS
*		301-763-8136
*		301-763-8108 FAX
*
*	Mailing Address:
*		810 NSC E/RA-14
*		NOAA/NESDIS
*		Washington, D.C. 20233
*
*	Email:
*		Kleespies@NZMS.WWB.NOAA.GOV
*
*

	Implicit None
        save

	Include 'constants_optran.inc'

	Integer Nprofiles	! number of profiles to use
				! this is the dimension of MaxAbsorber
				! and the second dimension of W.
	Real MaxAbsorber(Nprofiles)	! absorber amount at W(200)
	Real MinAbsorber	! absorber amount at W(1), W(0) == 0.0

	Real W(0:Nw,Nprofiles) ! standard absorber profiles

	Real*8 alpha    ! exponential parameter

	Integer k     ! utility variables

	Real*8 Tolerance 
        data Tolerance /1.0e-10/ ! iteration-to-iteration convergence criterion

	Integer MaxIter 
        data MaxIter /100/      ! arbitrary value to test non-convergence

	Real*8 x1,x2,fx,fp ! Newton's method variables
			   ! fx is the function
			   ! fp is the function primed (as in it's derivative)
			   ! x1 and x2 are previous and new approximations

	Integer Ichan,Iter ! local variables

	Do Ichan = 1 , NProfiles  ! profile loop

	 W(0,Ichan) = 0		  ! top of atmosphere always has zero absorber
	 W(1,Ichan) = MinAbsorber
	 W(Nw,Ichan) = MaxAbsorber(Ichan)
	
*	solve for alpha by Newton's method

* these are somewhat arbitrary initial guesses
* the lower one is important for the wet atmospheres, and the
* upper one is important for the dry atmospheres.

	 x2 = 0.3	
	 x1 = 3000.
	 iter = 0
	 Do While (abs(x2-x1) .ge. tolerance)

	  iter = iter + 1
	  If(Iter .gt. MaxIter) Then
	   Write(6,*) 'Make_Levels failed to converge '
	   stop
	  EndIf

	  x1 = x2
! 	the function
	  fx = (exp(Nw*x1)-1.) - (W(Nw,Ichan)/W(1,Ichan))*(exp(x1)-1.)
! 	function primed 
	  fp = Nw*exp(Nw*x1) - (W(Nw,Ichan)/W(1,Ichan))* exp(x1)   
	  x2 = x1 - fx/fp

	 EndDo

	 alpha = x2

*	Now that we have alpha, construct the profile 

	 Do k = 2 , Nw
	  W(k,Ichan) = W(1,Ichan) * (exp(k*alpha)-1.)/(exp(alpha)-1.)
	 EndDo

	EndDo ! profile loop

	return
	end





	Real Function CompBright(Vnu,T,Tau,Tskin,N)
c	Computes brightness temperature for a temperature and
c	transmittance profile.

*	Input
*
*	R*4 Vnu		Wavenumber
*	R*4 T(N)	Temperature profile
*	R*4 Tau(N)	Transmittance profile
*	R*4 Tskin	Skin temperature
*	I*4 N		Number of levels

	Implicit None
	
	Real Vnu
	Integer N
	Real T(N)
	Real Tau(N)
	Real Tskin

	Real Sum,b1,b2,bs,tau1,tau2
	Integer i

	Real C1 
        data C1 /1.1905e-5/
	Real C2 
        data C2 /1.4385/

	Real Planck,Bright,V,Temp,Radiance
	Planck(V,Temp) = (C1*v*v*v)/(Exp(C2*V/Temp) - 1)
        Bright(V,Radiance) = C2*V / Log(C1*V*V*V/Radiance + 1.0)

	Sum=0.

	b1=Planck(Vnu,T(1))
	tau1=tau(1)

	do i=2,N
	 b2=Planck(Vnu,T(i))
	 TAU2 = TAU(I)
	 Sum=Sum+.5*(b1+b2)*(tau1-tau2)
	 b1=b2
	 tau1=tau2
	EndDo

	bs=Planck(Vnu,Tskin)
	Sum=Sum+bs*tau(N)

	CompBright = 0.0
	if(Sum.gt.0.) CompBright=Bright(Vnu,Sum)

	return
	end

