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
      Subroutine Soil_In5(imax,jmax,SoilType,IStatus) 

      include 'soilm.inc'
      integer*4 imax,jmax
      integer*4 nf
      integer SoilType(imax,jmax)
      character*150 c_dir
      character*256 filename

      istatus = -1

      call get_directory('static',c_dir,lend)
      filename=c_dir(1:lend)//'soil/Soils.dat'

      open(Unit = 2, File = filename, Status = 'Old',
     1  Access = 'Sequential', Iostat = IERR, ERR = 664)

      do J = 1 , Jmax
         Read(2,*) (SoilType(I,J), I = 1, Imax)
      enddo
      close(2)

      nf = index(filename,' ')-1
      write(6,*) 'Got Soils Data from ',filename(1:nf)
      istatus = 0
      return

664   write(6,*)'Using Default Soil Types'
      do J = 1 , Jmax
         do I = 1, Imax
            soiltype(I,J) = 5
         enddo
      enddo
      istatus = 0

      Return
      End
