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
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
      subroutine analz_gps (lat,lon,wt,nn,glat,glon,
     1                       data_out,data_weights,ii,jj,istatus)

      implicit none

      integer ii,jj,nn,istatus
      real lat(nn),lon(nn),wt(nn)
      real data_out(ii,jj)
      real data_weights(ii,jj)
      real glat(ii,jj)
      real glon(ii,jj)

c     volitile array points
      real points (nn,3)
      integer mask (ii,jj)
      real r50
      real ri, rj


      integer i,j,n,ncount
      real sum

c
      istatus =  0

c     initialaize data_out on first call  (-1 used as missing)
c     since zero water is valid number

      do i = 1,ii
         do j = 1,jj
            data_out(i,j) = -1.
            mask (i,j) = 0
            data_weights(i,j) = 0.0
         enddo
      enddo



c     foreach n element of wt, determine its location in ii,jj space

      ncount = 0

      do n = 1,nn

         call  latlon_to_rlapsgrid(lat(n),lon(n),glat,glon,ii,jj,
     1        ri,rj, istatus)

         if (istatus.eq.1 .and. wt(n) .gt. 0.0) then
            ncount = ncount + 1
            i = nint (ri)
            j = nint (rj)
            data_out (i,j) = wt(n)
            mask(i,j) = 1
            data_weights (i,j) = 1.
            points(ncount,1) = wt(n)
            points(ncount,2) = i
            points(ncount,3) = j

            write(6,*)lat(n),lon(n), glat(i,j), glon(i,j), wt(n)
         else
            continue
         endif

      enddo

c     now that data_out is as full as it is going to get with the water
c     process the remainder.

c     compute the fraction of data_out that is empty

      sum = 0.

      do i = 1,ii
         do j = 1,jj
            if (data_out (i,j) .eq. -1.) then
               sum = sum + 1
               data_out(i,j) = points(1,1) ! helps converge
            endif
         enddo
      enddo

      if(ncount.eq.0) then      !abort gvap here
         write(6,*) 'No GVAP/GPS data avail to process... abort'
         istatus = 0
         return
      endif

      write (6,*) ncount, ' out of ', n, ' total avial data used'
      write (6,*) 'gvap/gps field fraction empty = ', sum/ii/jj

c      if (sum/ii/jj .gt. 0.75) return ! return if not enough data

      
c     now have fairly full data array.  now analyze

      call prep_grid (ii,jj,data_out,nn,points,ncount,istatus)
      if (istatus == 0) then
         write (6,*) 'error in prep_grid'
         return
      endif
      call slv_laplc (data_out,mask,ii,jj)

c     prep the weighting array for the above analyzed sheet

      r50 = 15.e+3

      call weight_field (data_weights, mask,  ii,jj,r50 , istatus)

      if (istatus .ne. 1) then! test weight_field
         write (6,*) 'Failure in weight_field from analz_gvap'
         return
      endif

c      call slv_laplc (data_weights, mask, ii,jj)

c     test NaN values coming out of data_out and data_weights

      call check_nan2 (data_out,ii,jj,istatus)
      if (istatus.ne.1) then 
         write (6,*) 'NaN detected in data_out in analz_gvap.f'
         return
      endif
      call check_nan2 (data_weights,ii,jj,istatus)
      if (istatus.ne.1) then 
         write (6,*) 'NaN detected in data_weights in analz_gvap.f'
         return
      endif

      istatus = 1

      return

      end


