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
      integer imax,jmax
      integer nf
      integer SoilType(imax,jmax)
      real    r_missing_data
      real,allocatable::static_stl(:,:)
      character*150 c_dir
      character*256 filename
      character*3   var
      character*150 directory
      character*31  ext
      character*10  units
      character*125 comment
c
c Current categories for top layer soil texture types
c in static file
c ------------------------------
c  1          SAND
c  2          LOAMY SAND
c  3          SANDY LOAM
c  4          SILT LOAM
c  5          SILT
c  6          LOAM
c  7          SANDY CLAY LOAM
c  8          SILTY CLAY LOAM
c  9          CLAY LOAM
c 10          SANDY CLAY
c 11          SILTY CLAY
c 12          CLAY
c 13          ORGANIC MATERIALS
c 14          WATER
c 15          BEDROCK
c 16          OTHER (land-ice)
c
c Expected texture types for this lsm
C 1. loamy sand:
C 2. sandy loam:
C 3. loam: 
C 4. sandy clay loam:
C 5. silty clay loam:
C 6. silty clay
c
c Mapping between 16 category and 6 category.
c
c if type = 0 then type = 5 ! By JRS. Cannot have type = 0;
c                             default to original setting
c if type16 = 1 then type6 = 1
c if type16 = 2 then type6 = 1
c if type16 = 3 then type6 = 2
c if type16 = 4 then type6 = 5
c if type16 = 5 then type6 = 4
c if type16 = 6 then type6 = 3
c if type16 = 7 then type6 = 4
c if type16 = 8 then type6 = 5
c if type16 = 9 then type6 = 5
c if type16 =10 then type6 = 6
c if type16 =11 then type6 = 6
c if type16 =12 then type6 = 6
c if type16 =13 then type6 = 1
c if type16 =14 then type6 = 5 !absurd but was = 5 prior to this
c if type16 =15 then type6 = 5 !absurd but was = 5 prior to this
c if type16 =16 then type6 = 5 !absurd but was = 5 prior to this
c

      istatus = -1

      allocate(static_stl(imax,jmax))
      ext='static'
      call get_directory(ext,c_dir,lend)
      call get_r_missing_data(r_missing_data,istatus)

c     filename=c_dir(1:lend)//'soil/Soils.dat'
c     open(Unit = 2, File = filename, Status = 'Old',
c    1  Access = 'Sequential', Iostat = IERR, ERR = 664)
c     do J = 1 , Jmax
c        Read(2,*) (SoilType(I,J), I = 1, Imax)
c     enddo
c     close(2)
c     nf = index(filename,' ')-1
c     write(6,*) 'Got Soils Data from ',filename(1:nf)
c     istatus = 0
c     return

      var='STL'
      ext='nest7grid'
      call rd_laps_static(c_dir,ext,imax,jmax,1,var,units
     .,comment,static_stl,gridspace,istatus)
      if(istatus.ne.1)then
         print*,'Error reading static file for soil type'
         return
      endif

c664   write(6,*)'Using Default Soil Types'
      print*,' Using static STL soil texture '
      do J = 1 , Jmax
         do I = 1, Imax
            if(static_stl(i,j).eq.0..or.
     +         static_stl(i,j).eq.r_missing_data)then
             soiltype(i,j) = 5. !default to 5 if no type for this grid point
            elseif(static_stl(i,j).eq.1.)then
             soiltype(i,j) = 1.
            elseif(static_stl(i,j).eq.2.)then
             soiltype(i,j) = 1.
            elseif(static_stl(i,j).eq.3.)then
             soiltype(i,j) = 2.
            elseif(static_stl(i,j).eq.4.)then
             soiltype(i,j) = 5.
            elseif(static_stl(i,j).eq.5.)then
             soiltype(i,j) = 4.
            elseif(static_stl(i,j).eq.6.)then
             soiltype(i,j) = 3.
            elseif(static_stl(i,j).eq.7.)then
             soiltype(i,j) = 4.
            elseif(static_stl(i,j).eq.8.)then
             soiltype(i,j) = 5.
            elseif(static_stl(i,j).eq.9.)then
             soiltype(i,j) = 5.
            elseif(static_stl(i,j).eq.10.)then
             soiltype(i,j) = 6.
            elseif(static_stl(i,j).eq.11.)then
             soiltype(i,j) = 6.
            elseif(static_stl(i,j).eq.12.)then
             soiltype(i,j) = 6.
            elseif(static_stl(i,j).eq.13.)then
             soiltype(i,j) = 1.
            elseif(static_stl(i,j).eq.14.)then
             soiltype(i,j) = 5.
            elseif(static_stl(i,j).eq.15.)then
             soiltype(i,j) = 5.
            elseif(static_stl(i,j).eq.16.)then
             soiltype(i,j) = 5.
            endif
         enddo
      enddo

      do j=1,jmax
      do i=1,imax
         if(soiltype(i,j).eq.0)then
           print*,'i/j= ',i,j,'soil type = 0 at i/j'
           print*,'Soil type value = 0 in soil_in5.f'
           print*,'at i/j = ',i,j
           print*,'!! Terminating !!'
           return
         endif
      enddo
      enddo

      deallocate (static_stl)

      istatus = 0

      Return
      End
