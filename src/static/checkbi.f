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
c
c
      program checkbi
c
c     Program to read and print the binary static files.
c
c     P. Stamus     02-20-95
c
c
      parameter(ni = 61, nj = 61)
      real*4 a(ni,nj)
      character infile*200
c
      write(6,200)
 200  format(1x,'Enter filename: ',$)
      read(5,300) infile
 300  format(a)
      open(51,file=infile,
     &         form='unformatted',status='old')
      read(51) a
      close(51)
c
      amax = -1.e30
      amin =  1.e30
      do j=1,nj
      do i=1,ni
         write(6,900) i,j,a(i,j)
         if(a(i,j) .gt. amax) then
            amax = a(i,j)
            i_mx = i
            j_mx = j
         endif
         if(a(i,j) .lt. amin) then
            amin = a(i,j)
            i_mn = i
            j_mn = j
         endif
      enddo !i
      enddo !j
 900  format(2i4,2x,f25.5)
c
      print *,' '
      write(6,910) amax, i_mx, j_mx
 910  format(1x,' Max of ',f25.5,' at ',2i4)
      write(6,920) amin, i_mn, j_mn
 920  format(1x,' Min of ',f25.5,' at ',2i4)
c
      stop
      end
