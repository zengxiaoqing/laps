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
	Subroutine Precip_Water_Profile(Q,P,ZA,N,Pw)
C
C************************************************************************
C*									*
C*	Module Name:	Precip_Water_Profile				*
C*									*
C*	Language:	Fortran-77	   Library:			*
C*	Version.Rev:	1.0 30 April 96	Programmer: Kleespies		*
C*									*
C*	Calling Seq:	Call Precip_Water_Profile(Q,P,Za,N,Pw)		*
C*									*
C*	Description:	Computes precipitable water profilefrom a mixing*
C*			ratio profile and a pressure profile, along	*
C*			a slant path.					*
C*									*
C*			This code largely lifted from RADCOM.		*
C*									*
C*	Input Args:	R*4 	Q(N)	Mixing ratio profile		*
C*			R*4	P(N)	Pressure profile		*
C*			R*4	ZA	Zenith Angle			*
C*			R*4	N	Number of points in profile	*
C*									*
C*	Output Args:	R*4	Precipitable Water profile		*
C*									*
C*	Common Blks:	None						*
C*									*
C*	Include:	NOne						*
C*									*
C*	Externals:	None						*
C*									*
C*	Data Files:	None						*
C*									*
C*	Restrictions:	Q,P,Pw dimensioned in calling routine.		*
C*									*
C*	Error Codes:	None						*
C*									*
C*	Error Messages:	None						*
C*									*
C************************************************************************
C
	Implicit NONE

*	Input
	Integer N
	Real Q(N),P(N),ZA  ! modified by Dan Birkenheuer 8/10/98

*	Output
	Real Pw(N)
	Real Path, CGrav 
	Data CGrav / 5.098581e-4 / ! .5/980.665
	Real Delp

	Integer i	! local variables

	Path = CGrav / Cos(ZA*acos(-1.)/180.)

	Pw(1) = 0.0
	
	Do i = 2 , N
	 Delp = Abs(P(i) - P(i-1))
	 Pw(i) = Pw(i-1) + Path*(Q(i) + Q(i-1))*Delp
	EndDo

	Return
	End
