	implicit none
	include 'Constants.inc'

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
