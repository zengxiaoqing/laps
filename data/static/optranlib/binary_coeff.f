	implicit none
c	include 'constants.inc'
********      Constants for program processing     ********************

	Integer Natm
	Parameter (Natm = 32)	! number of dependent atmospheres
	Integer Nlevel
	Parameter (Nlevel = 42) ! number of levels
	Integer Nw
	Parameter (Nw = 300)	! number of standard water levels (-1)
	Integer Nangle
	Parameter (Nangle = 6)	! number of look angles
	Integer Nchan
	Parameter (Nchan = 24)	! number of channels 

	Integer Nwet
	Parameter (Nwet = 5)   ! = number of predictors
	Integer Ndry
	Parameter (Ndry = 5)   ! = number of predictors
	Integer NOzo
	Parameter (NOzo = 5)   ! = number of predictors

	Integer MaxLevel
	Parameter (MaxLevel = 42) ! max # levels in pressure space

********              End Constants                *********************

	Real Bw(Nwet+1,Nw,Nchan)
	Real Bd(Ndry+1,Nw,Nchan)
	Real Bo(Nozo+1,Nw,Nchan)

	Integer i,j,k,Ichan

	Open(12,file='TOVS_NOAA14_COEFF.DAT;1',
     &		form='formatted',status='old')
	Open(42,file='Noaa14_wet_coeff.dat',
     &          form='unformatted',status='unknown')
	Open(43,file='Noaa14_dry_coeff.dat',
     &          form='unformatted',status='unknown')
	Open(44,file='Noaa14_ozo_coeff.dat',
     &          form='unformatted',status='unknown')

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
