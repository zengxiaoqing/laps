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
      subroutine analz_gvap (lat,lon,wt,nn,glat,glon,
     1                       data_out,data_weights,ii,jj,istatus)

      integer ii,jj,nn,istatus
      real lat(nn),lon(nn),wt(nn)
      real data_out(ii,jj)
      real data_weights(ii,jj)
      real glat(ii,jj)
      real glon(ii,jj)

c     volitile array points
      real points (nn,3)
      integer mask (ii,jj)


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
      call slv_laplc (data_out,mask,ii,jj)

c     prep the weighting array for the above analyzed sheet

      call weight_field (data_weights, mask,  ii,jj, 15.e+3, istatus)

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
      call check_nan2 (data_out,ii,jj,istatus)
      if (istatus.ne.1) then 
         write (6,*) 'NaN detected in data_weights in analz_gvap.f'
         return
      endif

      istatus = 1

      return

      end


