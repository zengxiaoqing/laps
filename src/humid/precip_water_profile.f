cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps
cdis
cdis    this software and its documentation are in the public domain and
cdis    are furnished "as is."  the united states government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  they assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  if significant modifications or enhancements
cdis    are made to this software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
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
