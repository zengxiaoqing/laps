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
cdis 	
        implicit none
	include 'constants_optran.inc'

	Real*8 Bw(Nwet+1,Nw,Nchan)
	Real*8 Bd(Ndry+1,Nw,Nchan)
	Real*8 Bo(Nozo+1,Nw,Nchan)

	Integer i,j,k,Ichan

	Open(12,file='coef.dat',
     &		form='formatted')
	Open(42,file='wet_coeff.dat',
     &          form='unformatted')
	Open(43,file='dry_coeff.dat',
     &          form='unformatted')
	Open(44,file='ozo_coeff.dat',
     &          form='unformatted')

	Do Ichan = 1 , Nchan
1200     Format(i4,6e20.12)

	 do i = 1 , Nw
	  Read(12,1200) k,(Bd(j,i,Ichan),j=1,Ndry+1)
	 EndDo
	 Do i = 1 , Nw
	  Read(12,1200) k,(Bw(j,i,Ichan),j=1,Nwet+1)
  	 EndDo
	 Do i = 1 , Nw
	  Read(12,1200) k,(Bo(j,i,Ichan),j=1,Nozo+1)
  	 EndDo

	EndDo

	Write(42) Bw
	Write(43) Bd
	Write(44) Bo

	Close(42)
	Close(43)
	Close(44)


	Stop
	end
