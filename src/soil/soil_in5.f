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
C   Subroutine to read in the soil type infomation
C   Created 5/2/92
C   Chandran Subramaniam
C
C
	Subroutine Soil_In5(SoilType,IStatus) 

      include 'lapsparms.for'
      parameter(imax=nx_l,jmax=ny_l)
      integer*4 imax,jmax
        include 'soilm.inc'
        integer SoilType(Imax,Jmax)
	Open(Unit = 2, File = 'Soils.dat', Status = 'Old',
     1  Access = 'Sequential', Iostat = IERR, ERR = 664)
	Do J = 1 , Jmax
            Read(2,*) (SoilType(I,J), I = 1, Imax)
        Enddo
        Close(2)
        Write(*,*) 'Read Soils Data'
	IStatus = -1
        Return
664	WRite(*,*)'Using Default Soil Types'
	Do J = 1 , Jmax
            Do I = 1, Imax
               SoilType(I,J) = 5
            Enddo
        Enddo
	IStatus = -1

	Return
	End
